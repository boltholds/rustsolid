//! Compound analytical region tests. No faceted Boolean approximation
//! is accepted for mathematical topology, mass or source curve ownership.
use rustsolid::{
    BlindHoleParameters,BlindRegionKind,BlindRegionComposer,BrepBody,
    CompoundBlindBrep,FaceId,FaceSense,Point3,Solid,Surface3,
};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn h(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn mixed()->Vec<BlindHoleParameters>{
    vec![h(-1.5,-0.5,0.75,1.8),h(-0.4,-0.5,0.7,1.8),
        h(1.7,1.2,0.55,2.6)]
}
fn approx(a:f64,b:f64,eps:f64){
    assert!((a-b).abs()<eps,"expected {b}, got {a}, error {}",(a-b).abs());
}
fn certify(body:&CompoundBlindBrep,groups:usize){
    body.validate().unwrap();
    let m=body.model();
    m.validate().unwrap();
    let s=m.summary();
    assert_eq!(s.genus,0);
    assert_eq!(s.euler_characteristic,2);
    assert_eq!(s.shells,1);
    assert_eq!(m.faces[1].loops.len(),1+groups);
    assert_eq!(m.faces[0].loops.len(),1);
    assert_eq!(m.loops.iter().filter(|l|l.role==rustsolid::LoopRole::Inner).count(),groups);
    for c in body.components(){
        assert_eq!(m.face_sense(c.floor_face),Some(FaceSense::Reversed));
        assert!(matches!(m.face_surface(c.floor_face),Some(Surface3::Plane(_))));
        for &id in &c.cylinder_faces{
            assert_eq!(m.face_sense(id),Some(FaceSense::Reversed));
            assert!(matches!(m.face_surface(id),Some(Surface3::Cylinder(_))));
        }
    }
}
#[test]
fn crossing_pair_plus_disjoint_single_has_two_true_top_loops_and_floor_faces(){
    let b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    certify(&b,2);
    let s=b.model().summary();
    assert_eq!((s.vertices,s.edges,s.coedges),(14,21,42));
    assert_eq!((s.faces,s.loops),(11,13));
    assert_eq!(s.circular_edges,6);
    assert_eq!(s.cylindrical_faces,3);
    assert_eq!(b.components().len(),2);
    assert_eq!(b.components()[0].kind,BlindRegionKind::TwoCrossingCircles);
    assert_eq!(b.components()[0].source_indices,vec![0,1]);
    assert_eq!(b.components()[1].kind,BlindRegionKind::SingleCircle);
    assert_eq!(b.components()[1].source_indices,vec![2]);
    assert_eq!(b.components()[0].floor_face,FaceId(8));
    assert_eq!(b.components()[1].floor_face,FaceId(10));
    let [a,c,d]=[mixed()[0],mixed()[1],mixed()[2]];
    let pair=rustsolid::CrossedDiskRegion::from_box(&block(),[a,c]).unwrap();
    let expected_volume=144.0-
        pair.disk_union_area*pair.cut_depth-PI*d.radius*d.radius*d.depth;
    let expected_area=168.0+
        pair.exposed_perimeter*pair.cut_depth+2.0*PI*d.radius*d.depth;
    approx(b.model().mass.volume,expected_volume,1e-9);
    approx(b.model().mass.surface_area,expected_area,1e-9);
}
#[test]
fn two_disjoint_crossing_pairs_have_two_nonconvex_analytic_floors(){
    let params=[
        h(-1.75,-1.2,0.6,1.3),h(-0.7,-1.2,0.6,1.3),
        h(0.7,1.2,0.6,2.6),h(1.75,1.2,0.6,2.6)
    ];
    let b=CompoundBlindBrep::from_box(&block(),params.to_vec()).unwrap();
    certify(&b,2);
    let s=b.model().summary();
    assert_eq!((s.vertices,s.edges,s.coedges),(16,24,48));
    assert_eq!((s.faces,s.loops),(12,14));
    assert_eq!(s.circular_edges,8);
    assert_eq!(s.cylindrical_faces,4);
    assert_eq!(b.components()[0].kind,BlindRegionKind::TwoCrossingCircles);
    assert_eq!(b.components()[1].kind,BlindRegionKind::TwoCrossingCircles);
    assert_eq!((b.components()[0].floor_face.0,b.components()[1].floor_face.0),
        (8,11));
}
#[test]
fn all_independent_disks_keep_one_floor_per_source_and_different_depths(){
    let holes=vec![
        h(-1.8,-1.2,0.3,1.3),
        h(0.1,1.5,0.4,2.2),
        h(1.8,-1.2,0.5,2.9),
    ];
    let b=CompoundBlindBrep::from_box(&block(),holes.clone()).unwrap();
    certify(&b,3);
    assert_eq!(b.model().summary().faces,12);
    for c in b.components(){assert_eq!(c.circle_count(),1);}
    let expected=144.0-holes.iter().map(|h|
        PI*h.radius*h.radius*h.depth).sum::<f64>();
    approx(b.model().mass.volume,expected,1e-9);
}
#[test]
fn an_independent_pair_is_not_double_subtracted_by_a_disjoint_third_cut(){
    let src=mixed();
    let m=BlindRegionComposer::build(&block(),&src).unwrap();
    let pair=rustsolid::CrossedDiskRegion::from_box(&block(),[src[0],src[1]]).unwrap();
    let expected=144.0-pair.disk_union_area*src[0].depth-
        PI*src[2].radius.powi(2)*src[2].depth;
    approx(m.0.mass.volume,expected,1e-9);
    assert_eq!(m.1.len(),2);
}
#[test]
fn permuting_sources_keeps_exact_mass_and_area_despite_dense_renumbering(){
    let holes=mixed();
    let a=CompoundBlindBrep::from_box(&block(),holes.clone()).unwrap();
    let b=CompoundBlindBrep::from_box(&block(),
        vec![holes[2],holes[1],holes[0]]).unwrap();
    certify(&a,2);certify(&b,2);
    approx(a.model().mass.volume,b.model().mass.volume,1e-9);
    approx(a.model().mass.surface_area,b.model().mass.surface_area,1e-9);
    approx(a.model().mass.centroid.x,b.model().mass.centroid.x,1e-9);
    approx(a.model().mass.centroid.y,b.model().mass.centroid.y,1e-9);
    approx(a.model().mass.centroid.z,b.model().mass.centroid.z,1e-9);
}
#[test]
fn topological_repartition_from_crossed_pair_to_three_singles_is_atomic(){
    let mut b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let old_mass=b.model().mass;
    let old_f=b.model().summary().faces;
    let mut tx=b.begin_edit();
    tx.replace(1,h(-0.4,-0.5,0.25,1.8)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    certify(&b,3);
    assert_eq!(b.model().summary().faces,12);
    assert_ne!(b.model().summary().faces,old_f);
    assert_eq!(b.undo_depth(),1);
    assert!(b.undo().unwrap());
    certify(&b,2);
    assert_eq!(b.model().mass,old_mass);
    assert_eq!(b.revision(),2);
    assert!(b.redo().unwrap());
    certify(&b,3);
    assert_eq!(b.revision(),3);
}
#[test]
fn nested_or_third_connected_circle_or_tangent_is_refused(){
    for holes in [
        vec![h(-1.,0.,0.9,2.),h(0.,0.,0.9,2.),h(1.,0.,0.9,2.)],
        vec![h(0.,0.,1.1,1.4),h(0.1,0.,0.2,1.4)],
        vec![h(-1.,0.,0.7,1.4),h(0.4,0.,0.7,1.4)],
        vec![h(-1.,0.,0.7,1.4),h(0.,0.,0.7,1.5)],
    ] {
        assert!(CompoundBlindBrep::from_box(&block(),holes).is_err());
    }
}
#[test]
fn invalid_staging_and_invalid_late_geometric_commit_leave_memento_unchanged(){
    let mut b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let original=b.model().mass;
    {
        let mut tx=b.begin_edit();
        tx.replace(0,h(-1.4,-0.5,0.75,1.8)).unwrap();
        assert!(tx.replace(123,h(0.,0.,1.,2.)).is_err());
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=b.begin_edit();
        tx.replace(0,h(-1.4,-0.5,0.75,1.8)).unwrap();
        tx.replace(1,h(-0.4,-0.5,0.7,2.0)).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(b.revision(),0);
    assert_eq!(b.undo_depth(),0);
    assert_eq!(b.model().mass,original);
    certify(&b,2);
}
#[test]
fn compound_body_is_brep_queryable_and_standalone_display_is_available(){
    let shape=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let body=BrepBody::from(shape.clone());
    body.validate().unwrap();
    let model=body.shared().unwrap();
    assert_eq!(model.summary().genus,0);
    let graph=rustsolid::geometry_query(&model,
        "query {body {kind topology {faceCount loopCount genus} mass {volume} faces(id:10){id sense surface{kind}}}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(graph["body"]["kind"],"analytic_compound_blind_regions");
    assert_eq!(graph["body"]["topology"]["faceCount"],11);
    assert_eq!(graph["body"]["faces"][0]["sense"],"reversed");
    assert!(!body.tessellate(24).unwrap().triangles.is_empty());
}

#[test]
fn sixteen_independent_blind_regions_respect_bounded_finite_topology(){
    let positions=[-2.25,-0.75,0.75,2.25];
    let mut sources=Vec::new();
    for &x in &positions {
        for (index,&z) in positions.iter().enumerate(){
            sources.push(h(x,z,0.25,1.0+0.3*index as f64));
        }
    }
    assert_eq!(sources.len(),16);
    let body=CompoundBlindBrep::from_box(&block(),sources.clone()).unwrap();
    certify(&body,16);
    let t=body.model().summary();
    assert_eq!((t.vertices,t.edges,t.coedges),(40,60,120));
    assert_eq!((t.faces,t.loops),(38,54));
    assert_eq!(t.euler_characteristic,2);
    let removed=sources.iter().map(|h|PI*h.radius*h.radius*h.depth).sum::<f64>();
    approx(body.model().mass.volume,144.0-removed,1e-9);
    let display=body.tessellate(12).unwrap();
    assert_eq!(display.top_face_loop_count,17);
    assert_eq!(display.floor_regions.len(),16);
    assert_eq!(display.mesh.vertices.len(),8+2*16*12);
    assert!(display.faceting_volume_error>0.0);
    let one_more=sources.into_iter().chain(std::iter::once(
        h(0.0,0.0,0.1,2.0))).collect::<Vec<_>>();
    assert!(CompoundBlindBrep::from_box(&block(),one_more).is_err());
}
#[test]
fn common_volume_membership_tracks_mixed_voids_and_true_floor_boundaries(){
    use rustsolid::{BooleanRegionClassifier,VolumePosition};
    let body=BrepBody::from(
        CompoundBlindBrep::from_box(&block(),mixed()).unwrap());
    for (point,position) in [
        (p(-1.5,2.0,-0.5),VolumePosition::Outside),
        (p(-0.4,2.0,-0.5),VolumePosition::Outside),
        (p(1.7,2.0,1.2),VolumePosition::Outside),
        (p(1.7,0.4,1.2),VolumePosition::Boundary),
        (p(1.7,0.0,1.2),VolumePosition::Inside),
        (p(-1.5,1.2,-0.5),VolumePosition::Boundary),
        (p(-1.5,0.0,-0.5),VolumePosition::Inside),
        (p(0.7,2.0,-0.5),VolumePosition::Inside),
        (p(0.0,-1.0,0.0),VolumePosition::Boundary),
        (p(3.5,1.0,0.0),VolumePosition::Outside),
    ] {
        assert_eq!(BooleanRegionClassifier::classify_point(&body,point).unwrap(),
            position,"compound volume classification at {point:?}");
    }
}
#[test]
fn arena_can_inspect_compound_regions_but_does_not_mutate_dense_names(){
    use rustsolid::{BrepArena,TopologyEntity};
    let body=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let arena=BrepArena::from_compound_blind(body).unwrap();
    let before=arena.snapshot().unwrap();
    let floor=arena.handle(TopologyEntity::Face(FaceId(8))).unwrap();
    assert_eq!(arena.model().unwrap().summary().faces,11);
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(8))),Some(floor));
    assert_eq!(arena.undo_depth(),0);
    assert_eq!(arena.redo_depth(),0);
    arena.validate().unwrap();
    arena.check_snapshot(&before).unwrap();
    let mut detached=arena.detached_carrier_editor().unwrap();
    let face=detached.handle(rustsolid::CarrierOwner::Face(FaceId(8))).unwrap();
    let mut edit=detached.begin_transaction();
    edit.rebase_face_uv(face,0.1,0.0).unwrap();
    edit.commit().unwrap();
    detached.validate().unwrap();
    arena.check_snapshot(&before).unwrap();
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(8))),Some(floor));
}
