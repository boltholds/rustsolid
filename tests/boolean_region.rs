//! First-party solid-angle / analytic inequality and graph-cut region witnesses.
use rustsolid::{
    BooleanRegionClassifier as Classifier, BooleanOperation, BrepBody,
    CylindricalBrep, FaceId, Frame3, GeometryTolerance, GraphRegionState,
    GraphOperand, IntersectionGraphBuilder, Point2, Point3, Solid,
    VolumePosition,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance {GeometryTolerance::default()}
fn cube()->BrepBody {BrepBody::from(Solid::block(p(0.,0.,0.),2.,3.,4.).unwrap())}
fn cylinder()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}
fn classify(body:&BrepBody,point:Point3)->VolumePosition{
    Classifier::classify_point(body,point).unwrap()
}
#[test]
fn exact_cuboid_polyhedral_inside_outside_and_boundary_witnesses(){
    let a=cube();
    for point in [p(1.,1.,2.),p(0.1,0.1,0.1),p(1.99,2.99,3.99)] {
        assert_eq!(classify(&a,point),VolumePosition::Inside);
    }
    for point in [p(0.,1.,2.),p(1.,0.,2.),p(1.,3.,2.),p(2.,1.,2.),
        p(1.,2.,4.),p(0.,0.,0.)]{
        assert_eq!(classify(&a,point),VolumePosition::Boundary);
    }
    for point in [p(-0.1,1.,2.),p(1.,3.5,2.),p(1.,1.,4.1)]{
        assert_eq!(classify(&a,point),VolumePosition::Outside);
    }
}
#[test]
fn concave_polyhedral_body_uses_real_triangles_not_bbox(){
    let l=Solid::extrude_xz(&[
        Point2{x:0.,z:0.},Point2{x:3.,z:0.},Point2{x:3.,z:1.},
        Point2{x:1.,z:1.},Point2{x:1.,z:3.},Point2{x:0.,z:3.}
    ],2.).unwrap();
    let body=BrepBody::from(l);
    assert_eq!(classify(&body,p(0.5,1.,2.5)),VolumePosition::Inside);
    assert_eq!(classify(&body,p(2.,1.,0.5)),VolumePosition::Inside);
    assert_eq!(classify(&body,p(2.,1.,2.)),VolumePosition::Outside,
        "concave notch is inside bbox but outside actual solid");
    assert_eq!(classify(&body,p(1.,1.,2.)),VolumePosition::Boundary);
}
#[test]
fn cylinder_volume_is_analytic_and_caps_are_not_silent_inside(){
    let a=cylinder();
    for point in [p(0.,1.,0.),p(1.,4.,0.),p(0.,2.,1.5)]{
        assert_eq!(classify(&a,point),VolumePosition::Inside);
    }
    for point in [p(2.,2.,0.),p(0.,0.,0.),p(0.,5.,0.),
        p(2.,0.,0.),p(0.,5.,2.)]{
        assert_eq!(classify(&a,point),VolumePosition::Boundary);
    }
    for point in [p(2.1,2.,0.),p(0.,-0.01,0.),p(0.,5.01,0.)] {
        assert_eq!(classify(&a,point),VolumePosition::Outside);
    }
}
#[test]
fn rotated_cylinder_membership_uses_local_frame(){
    let frame=Frame3::from_axes(p(1.,2.,3.),p(0.,1.,0.),p(0.,0.,-1.),tol()).unwrap();
    let body=BrepBody::from(CylindricalBrep::new(frame,1.,4.,tol()).unwrap());
    assert_eq!(classify(&body,p(-1.,2.,3.)),VolumePosition::Inside);
    assert_eq!(classify(&body,p(-1.,3.,3.)),VolumePosition::Boundary);
    assert_eq!(classify(&body,p(1.,2.,3.)),VolumePosition::Boundary);
    assert_eq!(classify(&body,p(-1.,3.5,3.)),VolumePosition::Outside);
}
fn circle_graph()->(BrepBody,BrepBody,rustsolid::IntersectionGraph){
    let block=BrepBody::from(Solid::block(p(-3.,2.,-3.),6.,1.,6.).unwrap());
    let cylinder=cylinder();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&block.shared().unwrap(),FaceId(0),
            &cylinder.shared().unwrap(),FaceId(2)).unwrap();
    assert_eq!(graph.arcs.len(),1);
    assert!(graph.arcs[0].closed);
    (block,cylinder,graph)
}
#[test]
fn periodic_plane_cylinder_crossing_classifies_both_face_regions(){
    let (a,b,graph)=circle_graph();
    let plan=Classifier::classify_graph(&graph,&a,FaceId(0),&b,FaceId(2),
        BooleanOperation::Intersection).unwrap();
    assert_eq!(plan.state,GraphRegionState::ClassifiedPreview);
    assert!(plan.hold.is_none());
    assert_eq!(plan.samples.len(),6);
    assert!(!plan.can_apply());
    for patch in &plan.samples {
        assert_ne!(patch.positive_other,patch.negative_other);
        assert_ne!(patch.keep_positive,patch.keep_negative);
        assert!(!patch.reverse_orientation);
        if patch.operand==GraphOperand::First {assert_eq!(patch.face,FaceId(0));}
        else {assert_eq!(patch.face,FaceId(2));}
    }
    plan.validate_against(&graph,&a,FaceId(0),&b,FaceId(2)).unwrap();
}
#[test]
fn intersection_union_and_difference_have_explicit_per_patch_truth_tables(){
    let (a,b,graph)=circle_graph();
    let inter=Classifier::classify_graph(&graph,&a,FaceId(0),&b,FaceId(2),
        BooleanOperation::Intersection).unwrap();
    let union=Classifier::classify_graph(&graph,&a,FaceId(0),&b,FaceId(2),
        BooleanOperation::Union).unwrap();
    let difference=Classifier::classify_graph(&graph,&a,FaceId(0),&b,FaceId(2),
        BooleanOperation::Difference).unwrap();
    for (i,patch) in inter.samples.iter().enumerate() {
        let u=&union.samples[i];let d=&difference.samples[i];
        assert_eq!(patch.positive_other,u.positive_other);
        assert_eq!(patch.negative_other,u.negative_other);
        assert_eq!(patch.keep_positive,!u.keep_positive);
        assert_eq!(patch.keep_negative,!u.keep_negative);
        if patch.operand==GraphOperand::First {
            assert_eq!(u.keep_positive,d.keep_positive);
            assert_eq!(u.keep_negative,d.keep_negative);
            assert!(!d.reverse_orientation);
        } else {
            assert_eq!(patch.keep_positive,d.keep_positive);
            assert_eq!(patch.keep_negative,d.keep_negative);
            assert!(d.reverse_orientation);
        }
    }
}
#[test]
fn graph_region_witnesses_reject_tampering_without_mutation(){
    let (a,b,graph)=circle_graph();
    let mut plan=Classifier::classify_graph(&graph,&a,FaceId(0),&b,FaceId(2),
        BooleanOperation::Difference).unwrap();
    plan.samples[0].keep_positive=!plan.samples[0].keep_positive;
    assert!(plan.validate_against(&graph,&a,FaceId(0),&b,FaceId(2)).is_err());
    assert!(!plan.can_apply());
    a.validate().unwrap();b.validate().unwrap();
}
#[test]
fn nonintersecting_selected_faces_do_not_invent_cut_regions(){
    let upper=BrepBody::from(Solid::block(p(-3.,6.,-3.),6.,1.,6.).unwrap());
    let cyl=cylinder();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&upper.shared().unwrap(),FaceId(0),
            &cyl.shared().unwrap(),FaceId(2)).unwrap();
    assert!(graph.arcs.is_empty());
    let plan=Classifier::classify_graph(&graph,&upper,FaceId(0),&cyl,FaceId(2),
        BooleanOperation::Difference).unwrap();
    assert_eq!(plan.state,GraphRegionState::NoCut);
    assert!(plan.samples.is_empty());
}
#[test]
fn classification_is_stable_under_large_common_translation(){
    let shift=1e8;
    let block=BrepBody::from(Solid::block(p(shift-3.,shift+2.,shift-3.),6.,1.,6.).unwrap());
    let cylinder=BrepBody::from(CylindricalBrep::upright(
        p(shift,shift,shift),2.,5.,tol()).unwrap());
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&block.shared().unwrap(),FaceId(0),
            &cylinder.shared().unwrap(),FaceId(2)).unwrap();
    let plan=Classifier::classify_graph(&graph,&block,FaceId(0),&cylinder,FaceId(2),
        BooleanOperation::Intersection).unwrap();
    assert_eq!(plan.state,GraphRegionState::ClassifiedPreview);
    assert_eq!(plan.samples.len(),6);
}
