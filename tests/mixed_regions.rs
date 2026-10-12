//! First-party mathematical, topological and transactional regressions
//! for one B-rep containing actual *through* tunnels and independent
//! finite-depth blind regions, including one nonconvex crossing pair.
//! Triangulations are never used as sources of analytical mass.
use rustsolid::{
    BlindHoleParameters,BlindRegionKind,BrepBody,BrepArena,
    BooleanRegionClassifier,FaceId,FaceSense,HoleParameters,
    MixedThroughBlindBrep,MixedThroughBlindParameters,MixedRegionComposer,
    Point3,Solid,Surface3,TopologyEntity,VolumePosition,
};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn params(through:Vec<HoleParameters>,blind:Vec<BlindHoleParameters>)
    ->MixedThroughBlindParameters
{MixedThroughBlindParameters{through,blind}}
fn sample()->MixedThroughBlindBrep{
    MixedThroughBlindBrep::from_box(&block(),params(
        vec![t(-1.7,0.,0.5)],
        vec![b(1.5,0.,0.75,2.0)])).unwrap()
}
fn approx(a:f64,c:f64,eps:f64){
    assert!((a-c).abs()<=eps,"got {a}, expected {c}, residual {}",(a-c).abs());
}
fn check(body:&MixedThroughBlindBrep,through:usize,singles:usize,crossed:usize){
    body.validate().unwrap();
    let m=body.model();m.validate().unwrap();
    let sum=m.summary();
    assert_eq!(sum.vertices,8+2*through+2*singles+4*crossed);
    assert_eq!(sum.edges,12+3*through+3*singles+6*crossed);
    assert_eq!(sum.coedges,2*sum.edges);
    assert_eq!(sum.faces,6+through+2*singles+3*crossed);
    assert_eq!(sum.loops,6+3*through+3*singles+4*crossed);
    assert_eq!(sum.shells,1);
    assert_eq!(sum.genus,through);
    assert_eq!(sum.euler_characteristic,2-2*through as isize);
    assert_eq!(m.faces[0].loops.len(),1+through);
    assert_eq!(m.faces[1].loops.len(),1+through+singles+crossed);
    assert_eq!(m.loops.iter().filter(|lp|
        lp.role==rustsolid::LoopRole::Inner).count(),
        2*through+singles+crossed);
    assert_eq!(sum.circular_edges,2*through+2*singles+4*crossed);
    assert_eq!(sum.cylindrical_faces,through+singles+2*crossed);
    for face in &m.faces {
        assert_eq!(face.sense,if face.id.0<6{FaceSense::Forward}
            else{FaceSense::Reversed});
    }
    for c in body.components(){
        assert!(matches!(m.face_surface(c.floor_face),Some(Surface3::Plane(_))));
        for &face in &c.cylinder_faces{
            assert!(matches!(m.face_surface(face),Some(Surface3::Cylinder(_))));
        }
    }
}
#[test]
fn one_full_tunnel_plus_one_blind_floor_has_true_genus_one_shell(){
    let body=sample();check(&body,1,1,0);
    let m=body.model();
    assert_eq!(m.faces[0].loops.len(),2);
    assert_eq!(m.faces[1].loops.len(),3);
    assert_eq!(m.summary().faces,9);
    assert_eq!(body.components()[0].floor_face,FaceId(7));
    let remove_through=PI*0.5*0.5*4.0;
    let remove_blind=PI*0.75*0.75*2.0;
    approx(m.mass.volume,144.0-remove_through-remove_blind,1e-9);
    approx(m.mass.surface_area,
        168.0-2.0*PI*0.5*0.5+2.0*PI*0.5*4.0+
        2.0*PI*0.75*2.0,1e-9);
    approx(m.mass.centroid.y,1.0-(remove_blind*(2.0-1.0))/
        m.mass.volume,1e-9);
    let removed_x=remove_through*(-1.7)+remove_blind*1.5;
    approx(m.mass.centroid.x,-removed_x/m.mass.volume,1e-9);
}
#[test]
fn two_full_tunnels_plus_crossing_pair_and_single_floor_yield_genus_two(){
    let through=vec![t(-1.8,0.0,0.4),t(-0.4,-2.0,0.35)];
    let blind=vec![
        b(1.1,0.0,0.5,1.8),
        b(2.0,0.0,0.5,1.8),
        b(0.0,1.8,0.3,2.6),
    ];
    let body=MixedThroughBlindBrep::from_box(&block(),
        params(through,blind)).unwrap();
    check(&body,2,1,1);
    assert!(body.components().iter().any(|c|
        c.kind==BlindRegionKind::TwoCrossingCircles));
    assert!(body.components().iter().any(|c|
        c.kind==BlindRegionKind::SingleCircle));
    assert_eq!(body.model().summary().genus,2);
    assert_eq!(body.model().summary().euler_characteristic,-2);
}
#[test]
fn blind_pair_repartition_is_one_atomic_edit_without_changing_through_genus(){
    let blind=vec![
        b(0.7,0.0,0.5,1.5),
        b(1.9,0.0,0.5,1.5),
    ];
    let mut body=MixedThroughBlindBrep::from_box(&block(),
        params(vec![t(-1.5,0.0,0.45)],blind)).unwrap();
    check(&body,1,2,0);
    let original=body.model().mass;
    let mut tx=body.begin_edit();
    tx.replace_blind(1,b(1.55,0.0,0.5,1.5)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    check(&body,1,0,1);
    assert_ne!(body.model().summary().faces,11);
    assert_eq!(body.undo_depth(),1);
    assert!(body.undo().unwrap());
    check(&body,1,2,0);
    approx(body.model().mass.volume,original.volume,1e-9);
    assert_eq!(body.revision(),2);
    assert!(body.redo().unwrap());
    check(&body,1,0,1);
    assert_eq!(body.revision(),3);
}
#[test]
fn mixed_family_overlap_tangent_and_near_touch_are_rejected(){
    for bx in [-1.7,-0.5,-0.5+1e-11]{
        assert!(MixedThroughBlindBrep::from_box(&block(),params(
            vec![t(-1.7,0.0,0.6)],
            vec![b(bx,0.0,0.6,1.5)])).is_err());
    }
}
#[test]
fn independently_invalid_blind_region_cannot_be_hidden_by_through_cut(){
    for blind in [
        vec![b(2.8,0.,0.5,1.0)],
        vec![b(1.0,0.,0.7,0.0)],
        vec![b(1.0,0.,0.7,4.0)],
        vec![b(1.1,0.,0.6,1.0),b(1.6,0.,0.6,2.0)],
    ] {
        let result=MixedThroughBlindBrep::from_box(&block(),params(
            vec![t(-1.8,0.,0.35)],blind));
        assert!(result.is_err(),"invalid cuts unexpectedly accepted");
    }
}
#[test]
fn source_permutation_keeps_exact_analytical_volume_and_center(){
    let through=vec![t(-1.6,0.,0.4),t(-1.6,1.4,0.4)];
    let blind=vec![
        b(1.0,0.0,0.55,1.8),
        b(1.95,0.0,0.55,1.8),
        b(0.5,1.9,0.3,2.0),
    ];
    let mut th2=through.clone();th2.reverse();
    let mut bl2=blind.clone();bl2.reverse();
    let a=MixedThroughBlindBrep::from_box(&block(),params(through,blind)).unwrap();
    let b=MixedThroughBlindBrep::from_box(&block(),params(th2,bl2)).unwrap();
    check(&a,2,1,1);check(&b,2,1,1);
    approx(a.model().mass.volume,b.model().mass.volume,1e-9);
    approx(a.model().mass.surface_area,b.model().mass.surface_area,1e-9);
    approx(a.model().mass.centroid.x,b.model().mass.centroid.x,1e-9);
    approx(a.model().mass.centroid.y,b.model().mass.centroid.y,1e-9);
    approx(a.model().mass.centroid.z,b.model().mass.centroid.z,1e-9);
}
#[test]
fn mixed_undo_redo_across_both_families_never_changes_partial_source(){
    let mut body=sample();
    let old=body.model().mass;
    let mut tx=body.begin_edit();
    tx.replace_through(0,t(-1.8,0.15,0.45)).unwrap();
    tx.replace_blind(0,b(1.5,-0.1,0.7,2.25)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    check(&body,1,1,0);
    assert_ne!(body.model().mass.volume,old.volume);
    assert!(body.undo().unwrap());
    approx(body.model().mass.volume,old.volume,1e-9);
    assert_eq!(body.revision(),2);
    assert!(body.redo().unwrap());
    assert_eq!(body.revision(),3);
    check(&body,1,1,0);
}
#[test]
fn poisoned_and_late_invalid_mixed_transactions_preserve_mass_revision_history(){
    let mut body=sample();let original=body.model().mass;
    {
        let mut tx=body.begin_edit();
        tx.replace_through(0,t(-1.9,0.,0.4)).unwrap();
        assert!(tx.replace_blind(99,b(0.,0.,1.,2.)).is_err());
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=body.begin_edit();
        tx.replace_through(0,t(1.5,0.,0.65)).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(body.revision(),0);
    assert_eq!(body.undo_depth(),0);
    assert_eq!(body.model().mass,original);
    check(&body,1,1,0);
}
#[test]
fn mixed_brep_is_queryable_with_actual_analytic_topology_and_material_sense(){
    let body=sample();
    let m=body.model();
    let query=rustsolid::geometry_query(m,
        "query {body {kind topology {faceCount genus eulerCharacteristic circularEdges} faces(id:7) {id sense surface {kind}} mass {volume}}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(query["body"]["kind"],"analytic_mixed_through_blind_regions");
    assert_eq!(query["body"]["topology"]["genus"],1);
    assert_eq!(query["body"]["faces"][0]["sense"],"reversed");
    assert_eq!(query["body"]["faces"][0]["surface"]["kind"],"plane");
    let x=BrepBody::from(body.clone());x.validate().unwrap();
    assert_eq!(x.shared().unwrap().summary(),m.summary());
}
#[test]
fn common_volume_membership_distinguishes_through_void_blind_void_and_floor(){
    let body=BrepBody::from(sample());
    for (point,status) in [
        (p(-1.7,0.,0.),VolumePosition::Outside),
        (p(-1.7,-1.,0.),VolumePosition::Outside),
        (p(1.5,2.0,0.),VolumePosition::Outside),
        (p(1.5,1.0,0.),VolumePosition::Boundary),
        (p(1.5,0.0,0.),VolumePosition::Inside),
        (p(2.25,2.0,0.),VolumePosition::Boundary),
        (p(0.,0.,0.),VolumePosition::Inside),
        (p(3.5,0.,0.),VolumePosition::Outside),
    ]{
        assert_eq!(BooleanRegionClassifier::classify_point(&body,point).unwrap(),
            status,"wrong material at {point:?}");
    }
}
#[test]
fn read_only_arena_snapshot_does_not_claim_durable_semantic_identity_for_mixed_edits(){
    let body=sample();
    let arena=BrepArena::from_mixed_regions(body).unwrap();
    let snapshot=arena.snapshot().unwrap();
    assert_eq!(snapshot.model.summary().genus,1);
    assert_eq!(arena.revision(),0);
    assert_eq!(arena.undo_depth(),0);
    assert!(arena.handle(TopologyEntity::Face(FaceId(1))).is_some());
    arena.validate().unwrap();
}
#[test]
fn genuine_carrier_editor_accepts_mixed_analytic_curves_without_facet_standins(){
    let mut editor=rustsolid::BrepCarrierEditor::from_body(
        &BrepBody::from(sample())).unwrap();
    let h=editor.handle(rustsolid::CarrierOwner::Edge(rustsolid::EdgeId(12))).unwrap();
    let initial=editor.model().mass;
    let mut tx=editor.begin_transaction();
    tx.reparameterize_edge(h,1.0,0.3).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();
    assert_eq!(editor.model().mass,initial);
    assert_eq!(editor.model().summary().genus,1);
}
