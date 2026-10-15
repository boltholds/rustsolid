//! G4b1: source-backed cross-face radial pairing without fabricated shells.
use rustsolid::{
    BrepBody,BrepModel,BooleanOperation,FaceId,FacePairGraph,FacePairId,
    FacePatchAssemblyPlan,FacePatchAssemblyResult,GeometryTolerance,
    GlobalEdgeStitchPlan,GlobalEdgeStitchResult,GlobalStitchHold,
    GraphOperand,HoleParameters,IntersectionGraphBuilder,MultiTraceUvPlan,
    MultiTraceUvResult,Point3,Solid,StitchPatchInput,StitchedEdgeOrigin,
    WholeSolidRegionPlan,WholeSolidRegionResult,Surface3,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn as_box(o:Point3,dx:f64,dy:f64,dz:f64)->BrepBody{
    BrepBody::from(Solid::block(o,dx,dy,dz).unwrap())
}
fn annulus()->BrepBody{
    let base=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(rustsolid::ThroughHoleBrep::from_box(
        &base,HoleParameters{center_x:0.,center_z:0.,radius:1.}
    ).unwrap())
}
fn cutter()->BrepBody{as_box(p(0.,4.,-4.),1.7,2.,8.)}
fn face_in_plane(model:&BrepModel,axis:usize,level:f64)->FaceId{
    model.faces.iter().find_map(|face|{
        let Surface3::Plane(surface)=model.face_surface(face.id)? else{
            return None;
        };
        let pos=surface.frame.origin.array();
        let normal=surface.frame.normal.array();
        ((pos[axis]-level).abs()<1e-9 && normal[axis].abs()>0.9)
            .then_some(face.id)
    }).expect("required real B-rep source planar face")
}
fn source_pairs(
    left:&BrepModel,right:&BrepModel,operand:GraphOperand,face:FaceId,
)->Vec<FacePairGraph>{
    let builder=IntersectionGraphBuilder::new(tol()).unwrap();
    match operand{
        GraphOperand::First=>right.faces.iter().enumerate().map(|(index,f)|{
            FacePairGraph {
                id:FacePairId(index as u32),first_face:face,second_face:f.id,
                graph:builder.intersect_faces(left,face,right,f.id).unwrap(),
            }
        }).collect(),
        GraphOperand::Second=>left.faces.iter().enumerate().map(|(index,f)|{
            FacePairGraph {
                id:FacePairId(index as u32),first_face:f.id,second_face:face,
                graph:builder.intersect_faces(left,f.id,right,face).unwrap(),
            }
        }).collect(),
    }
}
fn prepared(
    left:&BrepBody,right:&BrepBody,operand:GraphOperand,face:FaceId,
    op:BooleanOperation,
)->(Vec<FacePairGraph>,FacePatchAssemblyPlan){
    let a=left.shared().unwrap();
    let b=right.shared().unwrap();
    let pairs=source_pairs(&a,&b,operand,face);
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,operand,face,tol()).unwrap();
    let MultiTraceUvResult::Ready(uv)=result else {
        panic!("G2 supported bounded pair required, got: {result:?}");
    };
    let result=WholeSolidRegionPlan::build(
        left,right,&pairs,&uv,op,tol()).unwrap();
    let WholeSolidRegionResult::Ready(g3)=result else {
        panic!("G3 supported full-solid region required, got: {result:?}");
    };
    let result=FacePatchAssemblyPlan::build(
        left,right,&pairs,&uv,&g3,tol()).unwrap();
    let FacePatchAssemblyResult::Ready(patch)=result else {
        panic!("G4a validated patch required, got: {result:?}");
    };
    (pairs,patch)
}
fn verify_plan(
    left:&BrepBody,right:&BrepBody,
    source_patches:&[StitchPatchInput<'_>],
)->GlobalEdgeStitchPlan {
    let response=GlobalEdgeStitchPlan::build(
        left,right,source_patches,tol()).unwrap();
    let GlobalEdgeStitchResult::Preview(plan)=response else {
        panic!("Expected G4b1 provenance ledger: {response:?}");
    };
    plan
}

#[test]
fn single_face_cannot_be_misrepresented_as_a_closed_boolean_shell(){
    let a=annulus();let b=cutter();
    let (pairs,patch)=prepared(
        &a,&b,GraphOperand::First,FaceId(1),BooleanOperation::Difference);
    let input=[StitchPatchInput{plan:&patch,pairs:&pairs}];
    let response=GlobalEdgeStitchPlan::build(&a,&b,&input,tol()).unwrap();
    match response {
        GlobalEdgeStitchResult::Preview(plan)=>{
            assert!(!plan.can_apply && !plan.closed_shell_proven);
            assert!(!plan.radially_paired);
            assert!(!plan.open_uses.is_empty());
            assert!(!plan.missing_source_faces.is_empty());
            assert_eq!(plan.selected_coedge_uses,patch.uses.iter()
                .filter(|u|u.kept).count());
        },
        // A G2 cut which ends inside an opposite face is not a source-
        // certified physical *cross-face* endpoint. Held is also correct.
        GlobalEdgeStitchResult::Held(
            GlobalStitchHold::UncertifiedCutEndpoint|
            GlobalStitchHold::UncertifiedIntermediateCutSplit
        )=>{},
        other=>panic!("G4b1 must not publish an incomplete shell: {other:?}"),
    }
}

#[test]
fn duplicate_source_face_is_a_nonmanifold_error(){
    let a=annulus();let b=cutter();
    let (pairs,patch)=prepared(
        &a,&b,GraphOperand::First,FaceId(1),BooleanOperation::Difference);
    let input=[
        StitchPatchInput{plan:&patch,pairs:&pairs},
        StitchPatchInput{plan:&patch,pairs:&pairs},
    ];
    assert!(matches!(
        GlobalEdgeStitchPlan::build(&a,&b,&input,tol()).unwrap(),
        GlobalEdgeStitchResult::Held(GlobalStitchHold::DuplicateSourceFace)
    ));
}

#[test]
fn tampered_local_coedge_coordinates_rejected_before_stitching(){
    let a=annulus();let b=cutter();
    let (pairs,mut patch)=prepared(
        &a,&b,GraphOperand::First,FaceId(1),BooleanOperation::Difference);
    patch.uses[0].start.x+=0.2;
    let input=[StitchPatchInput{plan:&patch,pairs:&pairs}];
    assert!(GlobalEdgeStitchPlan::build(&a,&b,&input,tol()).is_err());
}

#[test]
fn two_full_crossing_planar_source_faces_can_pair_only_their_true_cut(){
    // The two source rectangles intersect transversely along the Z axis.
    // The shared source-face intersection is bounded at z=+-2 on both
    // supports, with original analytic graph boundary witnesses.
    let a=as_box(p(-2.,0.,-2.),4.,4.,4.);
    let b=as_box(p(0.,-2.,-2.),4.,4.,4.);
    let am=a.shared().unwrap();
    let bm=b.shared().unwrap();
    let af=face_in_plane(&am,1,0.);
    let bf=face_in_plane(&bm,0,0.);
    let (apairs,apatch)=prepared(
        &a,&b,GraphOperand::First,af,BooleanOperation::Difference);
    let (bpairs,bpatch)=prepared(
        &a,&b,GraphOperand::Second,bf,BooleanOperation::Difference);
    let inputs=[
        StitchPatchInput{plan:&apatch,pairs:&apairs},
        StitchPatchInput{plan:&bpatch,pairs:&bpairs},
    ];
    let plan=verify_plan(&a,&b,&inputs);
    assert!(plan.topology_preflight.is_none(),
        "partial source-face coverage must not gain a radial shell certificate");
    assert!(plan.edges.iter().all(|edge|
        matches!(edge.source,StitchedEdgeOrigin::Cut{..})),
        "coincident edges from different operands are not original-edge twins");
    let cross:Vec<_>=plan.edges.iter().filter(|edge|
        matches!(edge.source,StitchedEdgeOrigin::Cut{..})).collect();
    assert_eq!(cross.len(),1,"one real bounded transverse cross-face cut");
    let edge=cross[0];
    assert_ne!(edge.source_faces[0].operand,edge.source_faces[1].operand);
    assert_eq!(edge.source_faces[0].face,af);
    assert_eq!(edge.source_faces[1].face,bf);
    assert_ne!(edge.radial_uses[0].patch,edge.radial_uses[1].patch);
    assert!(!plan.can_apply && !plan.closed_shell_proven);
    assert!(!plan.missing_source_faces.is_empty());
    plan.validate_against(&a,&b,&inputs,tol()).unwrap();
}

#[test]
fn empty_patch_set_refuses_to_declare_a_closed_shell(){
    let a=annulus();let b=cutter();
    assert_eq!(
        GlobalEdgeStitchPlan::build(&a,&b,&[],tol()).unwrap(),
        GlobalEdgeStitchResult::Held(GlobalStitchHold::BudgetExceeded)
    );
}

#[test]
fn forged_original_boundary_witness_cannot_be_silently_discarded_by_g2(){
    let a=as_box(p(-2.,0.,-2.),4.,4.,4.);
    let b=as_box(p(0.,-2.,-2.),4.,4.,4.);
    let am=a.shared().unwrap();
    let bm=b.shared().unwrap();
    let face=face_in_plane(&am,1,0.);
    let mut pairs=source_pairs(&am,&bm,GraphOperand::First,face);
    let witness=pairs.iter_mut().find(|pair|pair.graph.arcs.iter().any(|arc|
        arc.contact==rustsolid::TrimContact::BoundaryOverlap)).unwrap();
    let node=witness.graph.nodes.iter_mut()
        .find(|node|!node.boundary_sources.is_empty())
        .expect("real physical overlap must have source boundary witnesses");
    node.boundary_sources.clear();
    assert!(MultiTraceUvPlan::build_selected(
        &am,&bm,&pairs,GraphOperand::First,face,tol()).is_err(),
        "forged boundary overlap cannot pass the source graph check");
}
