//! Canonical analytic intersection topology: exact source carrier, physical
//! circular and vertical edges, UV pcurves, coedge links and atomic edits.
//! No Parasolid implementation or vendor outputs form any test oracle.
use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters,BrepBody,BrepSource,Curve2,Curve3,FaceId,FaceSense,
    HoleParameters,LoopRole,MixedCrossingBrep,MixedCrossingSurfacePlanner,
    MixedThroughBlindBrep,MixedThroughBlindParameters,Point3,Solid,Surface3,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius,depth}
}
fn sample()->MixedCrossingBrep{
    MixedCrossingBrep::from_box(&block(),t(0.,0.,0.5),b(0.5,0.,0.5,1.5)).unwrap()
}
fn close(actual:f64,expected:f64,eps:f64){
    assert!((actual-expected).abs()<eps,
        "{actual:.14} differs from {expected:.14} by {:.3e}",(actual-expected).abs());
}

#[test]
fn canonical_shared_shell_is_closed_genus_one_with_two_cap_inner_loops(){
    let body=sample();
    body.validate().unwrap();
    let m=body.model();
    m.validate().unwrap();
    let s=m.summary();
    assert_eq!((s.vertices,s.edges,s.coedges,s.faces,s.loops,s.shells),
        (14,23,46,11,13,1));
    assert_eq!((s.euler_characteristic,s.genus),(0,1));
    assert_eq!((s.circular_edges,s.cylindrical_faces),(7,4));
    assert_eq!((m.faces[0].loops.len(),m.faces[1].loops.len()),(2,2));
    assert_eq!(m.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),2);
    assert!(m.shells[0].closed);
    assert_eq!(m.shells[0].faces.len(),11);
    assert_eq!(m.faces[10].loops.len(),1);
    assert_eq!(m.face_sense(FaceId(10)),Some(FaceSense::Reversed));
    for fid in 6..10 {assert_eq!(m.face_sense(FaceId(fid)),Some(FaceSense::Reversed));}
    assert!(matches!(m.face_surface(FaceId(10)),Some(Surface3::Plane(_))));
    for fid in 6..10 {
        assert!(matches!(m.face_surface(FaceId(fid)),Some(Surface3::Cylinder(_))));
    }
}

#[test]
fn every_new_edge_has_opposite_twins_and_all_trim_carriers_are_analytic(){
    let body=sample();let m=body.model();
    for e in &m.edges[12..] {
        let [a,b]=e.coedges;
        let ca=&m.coedges[a.0 as usize];
        let cb=&m.coedges[b.0 as usize];
        assert_ne!(ca.reversed,cb.reversed);
        assert_eq!(ca.twin,b);assert_eq!(cb.twin,a);
        assert_eq!(ca.edge,e.id);assert_eq!(cb.edge,e.id);
        assert_ne!(ca.face,cb.face,"no degenerate same-face seam");
    }
    for c in &m.coedges[24..] {
        assert_eq!(m.next_in_loop(c.prev).unwrap(),c.id);
        let tr=m.coedge_pcurve(c.id).unwrap();
        assert!(matches!(tr,Curve2::Line(_)|Curve2::OrientedCircle(_)));
    }
    for e in &m.edges[12..19] {
        assert!(matches!(m.edge_curve(e.id),Some(Curve3::Circle(_))));
    }
    for e in &m.edges[19..23] {
        assert!(matches!(m.edge_curve(e.id),Some(Curve3::Line(_))));
    }
}

#[test]
fn independently_computed_area_volume_and_floor_material_sense_agree(){
    let body=sample();let m=body.model();
    let overlap=2.0*0.25*(PI/3.0)-0.125*(3.0_f64).sqrt();
    close(m.mass.volume,144.0-PI*0.25*4.0-(PI*0.25-overlap)*1.5,1e-9);
    close(m.mass.volume,140.14094873829234,1e-10);
    let planned=MixedCrossingSurfacePlanner::plan(
        &block(),body.through(),body.blind()).unwrap();
    close(m.mass.surface_area,planned.predicted_surface_area,1e-9);
    assert!(m.mass.centroid.is_finite());
    // The floor is a support plane with geometric normal -Y and
    // reversed material sense, hence the actual outward normal is +Y.
    let n=m.oriented_face_normal(FaceId(10),0.0,0.0).unwrap();
    close(n.y,1.0,1e-12);
    let out:BrepBody=body.clone().into();
    out.validate().unwrap();
    assert_eq!(out.shared().unwrap().summary(),m.summary());
    assert!(out.tessellate(32).is_err(),
        "never silently pretend an unverified mesh is watertight");
}

#[test]
fn source_updates_are_strong_atomic_with_undo_redo(){
    let mut body=sample();
    let initial=body.model().mass;
    assert_eq!(body.revision(),0);
    let mut tx=body.begin_edit();
    tx.replace_through(t(0.,0.,0.55)).unwrap();
    tx.replace_blind(b(0.55,0.,0.5,1.4)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    body.validate().unwrap();
    assert_ne!(body.model().mass.volume,initial.volume);
    assert!(body.undo().unwrap());
    assert_eq!(body.revision(),2);
    close(body.model().mass.volume,initial.volume,1e-12);
    assert!(body.redo().unwrap());
    assert_eq!(body.revision(),3);
    body.validate().unwrap();
}

#[test]
fn invalid_late_commit_is_never_partially_published(){
    let mut body=sample();
    let old=body.model().mass;
    {
        let mut tx=body.begin_edit();
        tx.replace_through(t(0.,0.,0.55)).unwrap();
        tx.replace_blind(b(1.,0.,0.5,1.5)).unwrap();
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=body.begin_edit();
        tx.replace_through(t(0.,0.,0.6)).unwrap();
        assert!(tx.replace_blind(b(f64::NAN,0.,0.5,1.5)).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(body.revision(),0);
    assert_eq!(body.model().mass,old);
    assert_eq!(body.undo_depth(),0);
    assert_eq!(body.redo_depth(),0);
    body.validate().unwrap();
}

#[test]
fn unsupported_arrangements_stay_denied_and_legacy_multihole_is_unchanged(){
    let s=block();
    for bad in [b(1.,0.,0.5,1.5),b(0.,0.,0.5,1.5),
        b(0.0,0.0,0.3,1.5),b(0.5,0.,0.5,4.0)] {
        assert!(MixedCrossingBrep::from_box(&s,t(0.,0.,0.5),bad).is_err());
    }
    assert!(MixedThroughBlindBrep::from_box(
        &s,MixedThroughBlindParameters{
            through:vec![t(0.,0.,0.5)],blind:vec![b(0.5,0.,0.5,1.5)]
        }).is_err(),"only a dedicated certified source may construct overlap");
}

#[test]
fn mixed_crossing_rejects_tampered_source_brep_model(){
    let body=sample();let mut m=body.model().clone();
    m.coedges[24].next=m.coedges[24].id;
    assert!(m.validate().is_err());
    let mut m=body.model().clone();
    m.geometry.coedges[24].face=FaceId(6);
    assert!(m.validate().is_err());
    let mut m=body.model().clone();
    m.faces[10].sense=FaceSense::Forward;
    assert!(m.validate().is_err());
}
