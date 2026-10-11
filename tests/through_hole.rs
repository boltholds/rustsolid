//! Exact analytic through-hole Boolean, genus-one topology, UV carrier
//! validation, model-space mass, atomic parametric edit and arena history.
use rustsolid::{
    BooleanKernel,BooleanOperation,BooleanOutcome,BooleanRegionClassifier,
    BrepArena,BrepBody,Curve3,FaceId,GeometryTolerance,HoleParameters,
    LoopRole,Point3,Solid,Surface3,ThroughHoleBrep,TopologyEntity,VolumePosition,
};
use std::f64::consts::PI;

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn cutter()->rustsolid::CylindricalBrep{
    rustsolid::CylindricalBrep::upright(p(0.,-1.,0.),1.,4.,
        GeometryTolerance::default()).unwrap()
}
fn hole()->ThroughHoleBrep{
    ThroughHoleBrep::from_box(&block(),HoleParameters{
        center_x:0.,center_z:0.,radius:1.0}).unwrap()
}
fn close(a:f64,b:f64){assert!((a-b).abs()<1e-9,"expected {b}; got {a}");}
fn topology(h:&ThroughHoleBrep){
    h.validate().unwrap();
    let m=h.model();m.validate().unwrap();
    let sum=m.summary();
    assert_eq!(sum.vertices,10);
    assert_eq!(sum.edges,15);
    assert_eq!(sum.faces,7);
    assert_eq!(sum.loops,9);
    assert_eq!(sum.coedges,30);
    assert_eq!(sum.shells,1);
    assert_eq!(sum.euler_characteristic,0);
    assert_eq!(sum.genus,1);
    assert_eq!(sum.cylindrical_faces,1);
    assert_eq!(sum.circular_edges,2);
    assert_eq!(sum.closed_edges,2);
    assert_eq!(sum.seam_edges,1);
    assert!(matches!(m.face_surface(FaceId(6)),Some(Surface3::Cylinder(_))));
    assert!(matches!(m.edge_curve(rustsolid::EdgeId(12)),Some(Curve3::Circle(_))));
    assert!(matches!(m.edge_curve(rustsolid::EdgeId(13)),Some(Curve3::Circle(_))));
    assert_eq!(m.faces[0].loops.len(),2);
    assert_eq!(m.faces[1].loops.len(),2);
    assert_eq!(m.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),2);
    assert_eq!(m.coedges[28].twin,rustsolid::CoedgeId(29));
    assert_eq!(m.coedges[29].twin,rustsolid::CoedgeId(28));
    assert_eq!(m.faces[6].loops.len(),1);
}
#[test]
fn direct_hole_is_true_genus_one_analytic_brep_not_facet_boolean(){
    let h=hole();topology(&h);
    let m=h.model();
    close(m.mass.volume,144.-4.*PI);
    close(m.mass.surface_area,168.+6.*PI);
    close(m.mass.centroid.x,0.);
    close(m.mass.centroid.y,1.);
    close(m.mass.centroid.z,0.);
}
#[test]
fn boolean_box_minus_axis_cylinder_produces_analytic_outcome(){
    let a=BrepBody::from(block());
    let b=BrepBody::from(cutter());
    let result=BooleanKernel::apply(&a,&b,BooleanOperation::Difference).unwrap();
    assert_eq!(result.face_sources.len(),7);
    assert_eq!(result.face_sources[6].result_face,FaceId(6));
    assert!(result.face_sources[6].reversed);
    match result.outcome{
        BooleanOutcome::Analytic(h)=>topology(&h),
        _=>panic!("expected exact analytic genus-one result"),
    }
    assert!(matches!(BooleanKernel::apply(&a,&b,BooleanOperation::Union),
        Err(rustsolid::BooleanError::UnsupportedOperand(_))));
}
#[test]
fn invalid_hole_placement_and_non_through_cutter_are_refused(){
    for q in [
        HoleParameters{center_x:2.1,center_z:0.,radius:1.},
        HoleParameters{center_x:0.,center_z:0.,radius:3.},
        HoleParameters{center_x:0.,center_z:0.,radius:0.},
        HoleParameters{center_x:f64::NAN,center_z:0.,radius:1.},
    ]{assert!(ThroughHoleBrep::from_box(&block(),q).is_err());}
    let short=rustsolid::CylindricalBrep::upright(p(0.,-1.,0.),1.,3.,
        GeometryTolerance::default()).unwrap();
    assert!(BooleanKernel::apply(&BrepBody::from(block()),
        &BrepBody::from(short),BooleanOperation::Difference).is_err());
}
#[test]
fn radius_and_center_edit_commit_atomically_with_undo_redo(){
    let mut h=hole();
    assert_eq!(h.revision(),0);
    let original=h.model().mass.volume;
    let mut tx=h.begin_edit();
    tx.set_radius(1.5).unwrap();
    tx.set_center(0.25,-0.5).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    close(h.model().mass.volume,144.-9.*PI);
    topology(&h);
    assert_eq!(h.undo_depth(),1);
    assert_eq!(h.redo_depth(),0);
    assert!(h.undo().unwrap());
    close(h.model().mass.volume,original);
    assert_eq!(h.revision(),2);
    assert!(h.redo().unwrap());
    close(h.model().mass.volume,144.-9.*PI);
    assert_eq!(h.revision(),3);
    topology(&h);
}
#[test]
fn invalid_stage_or_late_commit_failure_preserves_geometry_and_history(){
    let mut h=hole();
    let snapshot=h.model().mass;
    {
        let mut tx=h.begin_edit();
        tx.set_radius(1.3).unwrap();
        assert!(tx.set_radius(f64::NAN).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(h.revision(),0);
    assert_eq!(h.undo_depth(),0);
    close(h.model().mass.volume,snapshot.volume);
    {
        let mut tx=h.begin_edit();
        tx.set_radius(2.99).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(h.revision(),0);
    assert_eq!(h.undo_depth(),0);
    close(h.model().mass.volume,snapshot.volume);
    topology(&h);
}
#[test]
fn arena_hole_edit_preserves_handles_but_invalidates_snapshot_revision(){
    let mut arena=BrepArena::from_through_hole(hole()).unwrap();
    let original=arena.snapshot().unwrap();
    let inner_edge=arena.handle(TopologyEntity::Edge(rustsolid::EdgeId(12))).unwrap();
    let hole_face=arena.handle(TopologyEntity::Face(FaceId(6))).unwrap();
    assert_eq!(arena.update_hole(HoleParameters{
        center_x:0.1,center_z:0.2,radius:1.3}).unwrap(),1);
    assert!(arena.check_snapshot(&original).is_err());
    assert_eq!(arena.resolve(inner_edge).unwrap().entity,
        TopologyEntity::Edge(rustsolid::EdgeId(12)));
    assert_eq!(arena.resolve(hole_face).unwrap().entity,
        TopologyEntity::Face(FaceId(6)));
    assert_eq!(arena.undo_depth(),1);
    assert!(arena.undo().unwrap().is_some());
    assert_eq!(arena.revision(),2);
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(6))),Some(hole_face));
    assert!(arena.redo().unwrap().is_some());
    assert_eq!(arena.revision(),3);
    assert_eq!(arena.handle(TopologyEntity::Edge(rustsolid::EdgeId(12))),Some(inner_edge));
    arena.validate().unwrap();
}
#[test]
fn hole_membership_distinguishes_removed_material_and_curved_boundary(){
    let h=BrepBody::from(hole());
    for (point,state) in [
        (p(0.,1.,0.),VolumePosition::Outside),
        (p(1.,1.,0.),VolumePosition::Boundary),
        (p(2.,1.,0.),VolumePosition::Inside),
        (p(3.,1.,0.),VolumePosition::Boundary),
        (p(4.,1.,0.),VolumePosition::Outside),
    ] {
        assert_eq!(BooleanRegionClassifier::classify_point(&h,point).unwrap(),state);
    }
}
#[test]
fn facet_count_cannot_affect_analytic_hole_geometry(){
    let h=hole();
    let reference=h.model().mass;
    for segments in [12,16,64,128]{
        let mesh=h.tessellate(segments).unwrap();
        assert!(mesh.vertices.len()>30);
        assert_eq!(mesh.triangles.len(),mesh.triangle_faces.len());
        assert_eq!(mesh.vertices.len()%4,0);
        assert!(mesh.triangle_faces.iter().any(|id|id.0==6));
        assert!((h.model().mass.volume-reference.volume).abs()<1e-12);
    }
}
