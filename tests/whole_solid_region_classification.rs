//! RustSolid G3: classify complete planar G2 cells against the whole
//! opposite solid. No attempt is made to publish resulting B-rep topology.
use rustsolid::{
    BooleanOperation, BooleanRegionClassifier, BrepBody, BrepModel,
    CertifiedCellLocation, CylindricalBrep, FaceId, FacePairGraph,
    FacePairId, GeometryTolerance, GraphOperand, GraphState, HoleParameters,
    IntersectionGraphBuilder, MultiTraceUvPlan, MultiTraceUvResult,
    Point3, Solid, Surface3, ThroughHoleBrep, VolumePosition,
    WholeSolidRegionHold, WholeSolidRegionPlan, WholeSolidRegionResult,
};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tolerance()->GeometryTolerance{GeometryTolerance::default()}
fn box_body(min:Point3,dx:f64,dy:f64,dz:f64)->BrepBody{
    BrepBody::from(Solid::block(min,dx,dy,dz).unwrap())
}
fn box_hole()->BrepBody{
    let source=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(
        &source,HoleParameters{center_x:0.,center_z:0.,radius:1.},
    ).unwrap())
}
fn cutter()->BrepBody{box_body(p(0.,4.,-4.),1.7,2.,8.)}
fn face_x(model:&BrepModel,x:f64)->FaceId{
    model.faces.iter().find_map(|face|{
        let Surface3::Plane(plane)=model.face_surface(face.id)? else {return None;};
        ((plane.frame.origin.x-x).abs()<1e-10 &&
            plane.frame.normal.x.abs()>0.9).then_some(face.id)
    }).expect("vertical face at requested X")
}
/// Enumerate every opposite face, even if the intersection graph is Empty.
/// A local preview restricted to crossing faces is not a G3 certificate.
fn all_pairs_for_first_face(
    source:&BrepModel,face:FaceId,other:&BrepModel,
)->Vec<FacePairGraph>{
    let builder=IntersectionGraphBuilder::new(tolerance()).unwrap();
    other.faces.iter().enumerate().map(|(index,other_face)|{
        let graph=builder.intersect_faces(source,face,other,other_face.id)
            .unwrap();
        FacePairGraph {
            id:FacePairId(index as u32),
            first_face:face,second_face:other_face.id,graph,
        }
    }).collect()
}
fn uv_ready(
    left:&BrepModel,right:&BrepModel,
    pairs:&[FacePairGraph],face:FaceId,
)->MultiTraceUvPlan{
    let result=MultiTraceUvPlan::build_selected(
        left,right,pairs,GraphOperand::First,face,tolerance()).unwrap();
    let MultiTraceUvResult::Ready(plan)=result else {
        panic!("expected G2 source geometry: {result:?}");
    };
    plan
}
fn classify(
    left:&BrepBody,right:&BrepBody,
    pairs:&[FacePairGraph],plan:&MultiTraceUvPlan,
    op:BooleanOperation,
)->WholeSolidRegionResult{
    WholeSolidRegionPlan::build(
        left,right,pairs,plan,op,tolerance()).unwrap()
}
fn ready(result:WholeSolidRegionResult)->WholeSolidRegionPlan {
    let WholeSolidRegionResult::Ready(plan)=result else {
        panic!("expected G3 classified cells: {result:?}");
    };
    plan
}

#[test]
fn annular_source_and_box_cutter_have_inside_middle_and_outside_sides(){
    let source=box_hole();
    let other=cutter();
    let a=source.shared().unwrap();
    let b=other.shared().unwrap();
    let pairs=all_pairs_for_first_face(&a,FaceId(1),&b);
    assert_eq!(pairs.len(),b.faces.len());
    assert_eq!(pairs.iter().filter(|pair|pair.graph.state==GraphState::Ready)
        .count(),2);
    let uv=uv_ready(&a,&b,&pairs,FaceId(1));
    assert_eq!(uv.cells.len(),3);
    let classified=classify(&source,&other,&pairs,&uv,
        BooleanOperation::Difference);
    let plan=ready(classified);
    assert_eq!(plan.witnesses.len(),3);
    assert_eq!(plan.kept_cells,2);
    assert_eq!(plan.witnesses.iter().filter(|cell|
        cell.location==CertifiedCellLocation::Inside).count(),1);
    assert_eq!(plan.witnesses.iter().filter(|cell|
        cell.location==CertifiedCellLocation::Outside).count(),2);
    assert!(!plan.can_apply);
    assert!(plan.witnesses.iter().all(|cell|!cell.reverse_orientation));
    for evidence in &plan.witnesses {
        let observed=BooleanRegionClassifier::classify_point(
            &other,evidence.point).unwrap();
        assert_eq!(observed,match evidence.location {
            CertifiedCellLocation::Inside=>VolumePosition::Inside,
            CertifiedCellLocation::Outside=>VolumePosition::Outside,
        });
    }
    plan.validate_against(&source,&other,&pairs,tolerance()).unwrap();
    let mut modified=plan.clone();
    modified.witnesses[0].point.x+=0.25;
    assert!(modified.validate_against(&source,&other,&pairs,tolerance()).is_err());
}

#[test]
fn selection_policy_is_consistent_for_union_intersection_and_difference(){
    let source=box_hole();
    let other=cutter();
    let a=source.shared().unwrap();
    let b=other.shared().unwrap();
    let pairs=all_pairs_for_first_face(&a,FaceId(1),&b);
    let uv=uv_ready(&a,&b,&pairs,FaceId(1));
    let difference=ready(classify(&source,&other,&pairs,&uv,
        BooleanOperation::Difference));
    let union=ready(classify(&source,&other,&pairs,&uv,
        BooleanOperation::Union));
    let intersection=ready(classify(&source,&other,&pairs,&uv,
        BooleanOperation::Intersection));
    assert_eq!(difference.kept_cells,2);
    assert_eq!(union.kept_cells,2);
    assert_eq!(intersection.kept_cells,1);
    assert_eq!(difference.witnesses.iter().map(|w|w.location)
        .collect::<Vec<_>>(),intersection.witnesses.iter().map(|w|w.location)
        .collect::<Vec<_>>());
}

#[test]
fn partial_opposite_face_selection_is_held_even_with_valid_g2_cells(){
    let source=box_hole();
    let other=cutter();
    let a=source.shared().unwrap();
    let b=other.shared().unwrap();
    let pairs=all_pairs_for_first_face(&a,FaceId(1),&b)
        .into_iter().filter(|pair|pair.graph.state==GraphState::Ready)
        .collect::<Vec<_>>();
    assert_eq!(pairs.len(),2);
    let uv=uv_ready(&a,&b,&pairs,FaceId(1));
    let result=classify(&source,&other,&pairs,&uv,
        BooleanOperation::Difference);
    assert_eq!(result,WholeSolidRegionResult::Held(
        WholeSolidRegionHold::IncompleteOppositeFaceCoverage));
}

#[test]
fn upright_cylinder_classifies_three_vertical_face_strips_by_radial_bounds(){
    // The source face x=0 intersects an upright cylinder in exactly two
    // vertical generator lines z=-2 and z=2. Cylinder caps lie outside this
    // face's y=[1,4] trim. The G2 traces delimit 3 complete cells.
    let source=box_body(p(0.,1.,-3.),2.,3.,6.);
    let other=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),2.,5.,tolerance()).unwrap());
    let a=source.shared().unwrap();
    let b=other.shared().unwrap();
    let face=face_x(&a,0.);
    let pairs=all_pairs_for_first_face(&a,face,&b);
    let uv=uv_ready(&a,&b,&pairs,face);
    assert_eq!(uv.cells.len(),3);
    let plan=ready(classify(&source,&other,&pairs,&uv,
        BooleanOperation::Intersection));
    assert_eq!(plan.kept_cells,1);
    assert_eq!(plan.witnesses.iter().filter(|cell|
        cell.location==CertifiedCellLocation::Inside).count(),1);
    assert_eq!(plan.witnesses.iter().filter(|cell|
        cell.location==CertifiedCellLocation::Outside).count(),2);
    plan.validate_against(&source,&other,&pairs,tolerance()).unwrap();
}

#[test]
fn source_bodies_remain_identical_after_classification_and_validation(){
    let source=box_hole();
    let other=cutter();
    let original_source=format!("{source:?}");
    let original_other=format!("{other:?}");
    let a=source.shared().unwrap();let b=other.shared().unwrap();
    let pairs=all_pairs_for_first_face(&a,FaceId(1),&b);
    let uv=uv_ready(&a,&b,&pairs,FaceId(1));
    let plan=ready(classify(&source,&other,&pairs,&uv,
        BooleanOperation::Difference));
    plan.validate_against(&source,&other,&pairs,tolerance()).unwrap();
    assert_eq!(format!("{source:?}"),original_source);
    assert_eq!(format!("{other:?}"),original_other);
    assert!(!plan.can_apply);
}
