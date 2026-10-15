//! G4a: real source-bound local B-rep patch/coedge previews.
//! No physical EdgeId allocation, closed global shell or persisted edit.
use rustsolid::{
    BrepBody, BrepModel, BooleanOperation, FaceId, FacePairGraph, FacePairId,
    FacePatchAssemblyHold, FacePatchAssemblyPlan, FacePatchAssemblyResult,
    FaceSense, GeometryTolerance, GraphOperand, HoleParameters,
    IntersectionGraphBuilder, MultiTraceUvPlan, MultiTraceUvResult,
    MultiUvBoundaryProvenance, Point3, Solid, ThroughHoleBrep,
    WholeSolidRegionPlan, WholeSolidRegionResult,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn dist2(a:Point3,b:Point3)->f64 {
    (a.x-b.x).powi(2)+(a.y-b.y).powi(2)+(a.z-b.z).powi(2)
}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn source()->BrepBody{
    let base=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(
        &base,HoleParameters{center_x:0.,center_z:0.,radius:1.},
    ).unwrap())
}
fn cutter()->BrepBody {
    BrepBody::from(Solid::block(p(0.,4.,-4.),1.7,2.,8.).unwrap())
}
fn pairs(
    left:&BrepModel,right:&BrepModel,
    operand:GraphOperand,selected:FaceId,
)->Vec<FacePairGraph> {
    let builder=IntersectionGraphBuilder::new(tol()).unwrap();
    match operand {
        GraphOperand::First=>right.faces.iter().enumerate().map(|(i,face)|{
            FacePairGraph {
                id:FacePairId(i as u32),first_face:selected,
                second_face:face.id,
                graph:builder.intersect_faces(left,selected,right,face.id).unwrap(),
            }
        }).collect(),
        GraphOperand::Second=>left.faces.iter().enumerate().map(|(i,face)|{
            FacePairGraph {
                id:FacePairId(i as u32),first_face:face.id,
                second_face:selected,
                graph:builder.intersect_faces(left,face.id,right,selected).unwrap(),
            }
        }).collect(),
    }
}
fn plans(
    left:&BrepBody,right:&BrepBody,
    operand:GraphOperand,operation:BooleanOperation,
)->(Vec<FacePairGraph>, MultiTraceUvPlan, WholeSolidRegionPlan) {
    let a=left.shared().unwrap();
    let b=right.shared().unwrap();
    let face=FaceId(1); // Analytical annular planar face, source on chosen side.
    let selected_pairs=pairs(&a,&b,operand,face);
    let uv=MultiTraceUvPlan::build_selected(
        &a,&b,&selected_pairs,operand,face,tol()).unwrap();
    let MultiTraceUvResult::Ready(uv)=uv else {
        panic!("G2 did not construct the bounded patch: {uv:?}");
    };
    let result=WholeSolidRegionPlan::build(
        left,right,&selected_pairs,&uv,operation,tol()).unwrap();
    let WholeSolidRegionResult::Ready(g3)=result else {
        panic!("G3 did not classify bounded patch: {result:?}");
    };
    (selected_pairs,uv,g3)
}
fn prepared(
    left:&BrepBody,right:&BrepBody,
    pairs:&[FacePairGraph],uv:&MultiTraceUvPlan,g3:&WholeSolidRegionPlan,
)->FacePatchAssemblyPlan {
    let result=FacePatchAssemblyPlan::build(
        left,right,pairs,uv,g3,tol()).unwrap();
    let FacePatchAssemblyResult::Ready(plan)=result else {
        panic!("G4a did not prepare source-bound coedges: {result:?}");
    };
    plan
}

#[test]
fn first_operand_difference_keeps_only_material_exterior_patches_and_real_cut_uses(){
    let a=source();
    let b=cutter();
    let left_before=format!("{a:?}");
    let right_before=format!("{b:?}");
    let (pairs,uv,g3)=plans(&a,&b,GraphOperand::First,BooleanOperation::Difference);
    assert_eq!(uv.cells.len(),3);
    assert_eq!(g3.kept_cells,2);
    let plan=prepared(&a,&b,&pairs,&uv,&g3);
    assert_eq!(plan.source_operand,GraphOperand::First);
    assert_eq!(plan.source_face,FaceId(1));
    assert_eq!(plan.cut_pairs.len(),uv.cuts.len());
    assert_eq!(plan.cut_pairs.len(),3);
    assert_eq!(plan.open_cut_uses,3);
    assert!(plan.open_source_uses>0);
    assert!(!plan.can_apply);
    assert!(plan.kept_area_uv>0.0 && plan.kept_area_uv<uv.original_face_area);
    assert!(plan.retained_loops.len()>=2);
    for cut in &plan.cut_pairs {
        assert_eq!(cut.counterpart_operand,GraphOperand::Second);
        assert_ne!(cut.selected_use,cut.discarded_use);
        let selected=&plan.uses[cut.selected_use];
        let discarded=&plan.uses[cut.discarded_use];
        assert!(selected.kept);
        assert!(!discarded.kept);
        assert!(matches!(selected.provenance,MultiUvBoundaryProvenance::Cut{..}));
        assert_eq!(selected.provenance,discarded.provenance);
        assert!(dist2(selected.start,discarded.end)<1e-12,
            "opposite G2 uses must meet on the same analytic 3D cut");
    }
    let selected_source=plan.uses.iter().filter(|u|u.kept &&
        matches!(u.provenance,MultiUvBoundaryProvenance::Source{..})).count();
    assert_eq!(selected_source,plan.open_source_uses);
    plan.validate_against(&a,&b,&pairs,tol()).unwrap();
    assert_eq!(format!("{a:?}"),left_before);
    assert_eq!(format!("{b:?}"),right_before);
}

#[test]
fn second_operand_difference_reverses_selected_patch_normal_and_loop_order(){
    let left=cutter();
    let right=source();
    let (pairs,uv,g3)=plans(&left,&right,GraphOperand::Second,
        BooleanOperation::Difference);
    assert_eq!(g3.kept_cells,1,
        "subtractor's inside wall patches must be kept and reversed");
    let plan=prepared(&left,&right,&pairs,&uv,&g3);
    let model=right.shared().unwrap();
    let source_reversed=model.face_sense(FaceId(1)).unwrap()==FaceSense::Reversed;
    assert_eq!(plan.source_operand,GraphOperand::Second);
    assert_eq!(plan.open_cut_uses,3);
    assert_eq!(plan.cut_pairs.len(),3);
    assert!(plan.retained_loops.iter().all(|cycle|
        cycle.reversed_for_shell==!source_reversed));
    assert!(plan.uses.iter().filter(|u|u.kept).all(|u|
        u.reversed_for_shell==!source_reversed));
    for cycle in &plan.retained_loops {
        let first=&plan.uses[*cycle.uses.first().unwrap()];
        let last=&plan.uses[*cycle.uses.last().unwrap()];
        assert!(dist2(first.start,last.end)<1e-12);
        for win in cycle.uses.windows(2) {
            let a=&plan.uses[win[0]];
            let b=&plan.uses[win[1]];
            assert!(dist2(a.end,b.start)<1e-12,
                "reversed shell loops must remain directed and closed");
        }
    }
    plan.validate_against(&left,&right,&pairs,tol()).unwrap();
}

#[test]
fn tampering_g3_selection_or_g2_span_cannot_publish_candidate_coedges(){
    let a=source();
    let b=cutter();
    let (pairs,uv,g3)=plans(&a,&b,GraphOperand::First,BooleanOperation::Difference);
    let plan=prepared(&a,&b,&pairs,&uv,&g3);

    let mut swapped=g3.clone();
    swapped.witnesses[0].keep=!swapped.witnesses[0].keep;
    assert!(FacePatchAssemblyPlan::build(
        &a,&b,&pairs,&uv,&swapped,tol()).is_err());

    let mut tampered=uv.clone();
    tampered.cells[0].cycles[0].spans[0].start.u+=0.05;
    assert!(FacePatchAssemblyPlan::build(
        &a,&b,&pairs,&tampered,&g3,tol()).is_err());

    let mut fake=plan.clone();
    fake.uses[0].start.x+=1.0;
    assert!(fake.validate_against(&a,&b,&pairs,tol()).is_err());
    assert!(!plan.can_apply);
}

#[test]
fn switching_operation_or_source_body_rejects_old_patch_preview(){
    let a=source();let b=cutter();
    let (pairs,uv,g3)=plans(&a,&b,GraphOperand::First,BooleanOperation::Difference);
    let plan=prepared(&a,&b,&pairs,&uv,&g3);
    let mut wrong=plan.clone();
    // Union and Difference select the *same* first-operand outside patches.
    // Intersection changes this face's retained region and must invalidate it.
    wrong.operation=BooleanOperation::Intersection;
    assert!(wrong.validate_against(&a,&b,&pairs,tol()).is_err());

    let shifted=BrepBody::from(
        Solid::block(p(0.1,4.,-4.),1.7,2.,8.).unwrap());
    assert!(plan.validate_against(&a,&shifted,&pairs,tol()).is_err());
}

#[test]
fn malformed_cut_index_never_becomes_topology(){
    let a=source();let b=cutter();
    let (pairs,mut uv,g3)=plans(&a,&b,GraphOperand::First,
        BooleanOperation::Difference);
    uv.cuts[0].segment=9999;
    assert!(FacePatchAssemblyPlan::build(&a,&b,&pairs,&uv,&g3,tol()).is_err());
}

#[test]
fn local_cut_incidence_policy_fails_closed_without_global_brep_commit(){
    let a=source();let b=cutter();
    let (pairs,uv,g3)=plans(&a,&b,GraphOperand::First,
        BooleanOperation::Intersection);
    let plan=prepared(&a,&b,&pairs,&uv,&g3);
    assert_eq!(g3.kept_cells,1);
    assert_eq!(plan.cut_pairs.len(),3);
    assert_eq!(plan.open_cut_uses,plan.cut_pairs.len());
    assert!(plan.uses.iter().any(|use_|!use_.kept));
    assert!(!plan.can_apply,
        "G4a must never treat two source-local cut uses as a closed shell");
    assert!(plan.open_source_uses>0);
    assert_ne!(FacePatchAssemblyHold::CutIncidenceMismatch,
        FacePatchAssemblyHold::NonSeparatingCut);
}
