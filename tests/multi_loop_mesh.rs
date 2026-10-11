//! Independent display-only mesh contracts for planar regions with
//! multiple analytic circular inner loops and closed curved sidewalls.
use rustsolid::{
    BrepStructuralEditor, BrepArena, BrepBody, FaceId, HoleParameters, HoleEdgeRole,
    StructuralEdgeRef, StructuralFaceRef, Point3, Solid, Curve3, Surface3,
};
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn params(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn built(holes:&[HoleParameters])->BrepStructuralEditor{
    let mut editor=BrepStructuralEditor::new(block()).unwrap();
    if !holes.is_empty(){
        let mut tx=editor.begin_transaction();
        for &hole in holes {tx.add_through_hole(hole).unwrap();}
        tx.commit().unwrap();
    }
    editor.validate().unwrap();
    editor
}
fn assert_close(value:f64,expected:f64,epsilon:f64){
    assert!((value-expected).abs()<=epsilon,
        "got {value}, expected {expected}, error {}, limit {epsilon}",
        (value-expected).abs());
}
fn check_facets(
    mesh:&rustsolid::Mesh,
    source:&rustsolid::BrepModel,
    holes:usize,segments:usize,
) {
    let expected_cap_triangles=2+holes*(segments+2);
    assert_eq!(mesh.vertices.len(),8+2*holes*segments);
    assert_eq!(mesh.triangles.len(),
        2*expected_cap_triangles+8+2*holes*segments);
    assert_eq!(mesh.triangle_faces.len(),mesh.triangles.len());
    let mut edge_use=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut triangle_count=BTreeMap::<u32,usize>::new();
    for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        *triangle_count.entry(owner.0).or_default()+=1;
        assert!((owner.0 as usize)<source.faces.len());
        let [a,b,c]=*tri;
        for (u,v) in [(a,b),(b,c),(c,a)]{
            edge_use.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
        assert_ne!(a,b);assert_ne!(b,c);assert_ne!(a,c);
    }
    assert!(edge_use.values().all(|pair|
        pair.len()==2 && pair[0]!=pair[1]),
        "every mesh edge must have exactly two opposite triangle uses");
    assert_eq!(triangle_count[&0],expected_cap_triangles);
    assert_eq!(triangle_count[&1],expected_cap_triangles);
    for face in 2..6 {
        assert_eq!(triangle_count[&face],2);
    }
    for i in 0..holes{
        assert_eq!(triangle_count[&(i as u32+6)],2*segments);
    }
}
#[test]
fn zero_hole_face_regions_reproduce_exact_block_mesh_contract(){
    let e=built(&[]);
    let m=e.tessellate(16).unwrap();
    assert_eq!(m.bottom_region.loop_count,1);
    assert_eq!(m.bottom_region.triangles,2);
    check_facets(&m.mesh,e.model(),0,16);
    assert_close(m.faceted_volume,144.,1e-9);
    assert_close(m.faceted_surface_area,168.,1e-9);
    assert_close(m.faceting_volume_error,0.,1e-9);
}
#[test]
fn two_hole_regions_and_cylindrical_walls_form_one_watertight_mesh(){
    let holes=[params(-1.4,-0.4,0.7),params(1.4,0.5,0.9)];
    let e=built(&holes);
    let cap=e.model().summary();
    assert_eq!(cap.genus,2);
    let r=e.tessellate(48).unwrap();
    check_facets(&r.mesh,e.model(),2,48);
    assert_eq!(r.bottom_region.loop_count,3);
    assert_eq!(r.top_region.loop_count,3);
    assert_close(r.bottom_region.polygonal_area,r.top_region.polygonal_area,1e-10);
    assert_close(r.analytic_volume,144.-4.*PI*(0.7*0.7+0.9*0.9),1e-9);
    assert!(r.faceting_volume_error>0.0);
    assert_eq!(r.circle_segments,48);
    assert!(r.bottom_region.polygonal_area>0.);
}
#[test]
fn polygonal_error_converges_quadratically_without_mutating_analytic_mass(){
    let e=built(&[params(-1.4,0.,0.75),params(1.4,0.,0.75)]);
    let before=e.snapshot();
    let coarse=e.tessellate(12).unwrap();
    let normal=e.tessellate(48).unwrap();
    let fine=e.tessellate(192).unwrap();
    assert!(coarse.faceting_volume_error>normal.faceting_volume_error);
    assert!(normal.faceting_volume_error>fine.faceting_volume_error);
    assert!(fine.faceting_volume_error<normal.faceting_volume_error*0.08);
    for m in [&coarse,&normal,&fine] {
        assert_close(m.analytic_volume,before.model.mass.volume,1e-10);
        assert!(m.faceted_volume>=m.analytic_volume);
    }
    e.check_snapshot(&before).unwrap();
    assert_eq!(e.revision(),1);
}
#[test]
fn circle_sampling_uses_canonical_edge_carriers_after_compaction_and_edits(){
    let mut e=built(&[params(-1.45,0.,0.6),params(1.4,0.,0.75)]);
    let first=e.handle(0).unwrap();
    let second=e.handle(1).unwrap();
    {
        let mut tx=e.begin_transaction();
        tx.remove_through_hole(first).unwrap();
        tx.reparameterize_edge(StructuralEdgeRef::Hole{
            hole:second,role:HoleEdgeRole::BottomCircle},1.,0.6).unwrap();
        tx.rebase_face_uv(StructuralFaceRef::Base(FaceId(0)),0.15,-0.2).unwrap();
        tx.rebase_face_uv(StructuralFaceRef::HoleSide(second),0.,0.35).unwrap();
        tx.commit().unwrap();
    }
    let m=e.tessellate(64).unwrap();
    check_facets(&m.mesh,e.model(),1,64);
    assert!(matches!(e.model().edge_curve(rustsolid::EdgeId(12)),
        Some(Curve3::Circle(_))));
    assert!(matches!(e.model().face_surface(FaceId(6)),
        Some(Surface3::Cylinder(_))));
    assert_eq!(e.model().summary().genus,1);
}
#[test]
fn three_asymmetric_holes_generate_three_exact_face_ownership_sets(){
    let e=built(&[
        params(-1.75,-0.5,0.5),
        params(0.,1.3,0.4),
        params(1.7,-0.6,0.6)
    ]);
    let m=e.tessellate(36).unwrap();
    check_facets(&m.mesh,e.model(),3,36);
    assert_eq!(e.model().summary().genus,3);
    assert_eq!(m.bottom_region.loop_count,4);
}
#[test]
fn insertion_order_changes_vertex_numbering_but_not_polyhedral_mass(){
    let holes=[params(-1.5,0.,0.4),params(1.5,0.,0.8)];
    let a=built(&holes);
    let b=built(&[holes[1],holes[0]]);
    let ma=a.tessellate(64).unwrap();
    let mb=b.tessellate(64).unwrap();
    assert_close(ma.faceted_volume,mb.faceted_volume,1e-8);
    assert_close(ma.faceted_surface_area,mb.faceted_surface_area,1e-8);
    assert_eq!(ma.mesh.triangles.len(),mb.mesh.triangles.len());
    assert_eq!(ma.bottom_region.triangles,mb.bottom_region.triangles);
}
#[test]
fn derived_mesh_does_not_mutate_source_arena_or_structural_revision(){
    let source=BrepArena::from_solid(block()).unwrap();
    let snap=source.snapshot().unwrap();
    let mut e=source.detached_structural_editor().unwrap();
    let mut tx=e.begin_transaction();
    tx.add_through_hole(params(0.,0.,1.0)).unwrap();
    tx.commit().unwrap();
    let before=e.snapshot();
    e.tessellate(16).unwrap();
    e.tessellate(64).unwrap();
    e.check_snapshot(&before).unwrap();
    assert_eq!(e.revision(),1);
    source.check_snapshot(&snap).unwrap();
    assert_eq!(source.model().unwrap().mass.volume,144.0);
}
#[test]
fn invalid_resolutions_fail_without_exposing_an_unverified_mesh(){
    let e=built(&[params(0.,0.,1.)]);
    for n in [0,2,11,513,10000] {
        assert!(e.tessellate(n).is_err(),"invalid segmentation {n} must fail");
    }
    assert_eq!(e.model().summary().genus,1);
    e.validate().unwrap();
}
#[test]
fn center_radius_sweep_preserves_watertight_mesh_and_exact_mass(){
    for radius in [0.25,0.45,0.8,1.25] {
        let e=built(&[params(0.3,-0.3,radius)]);
        for segments in [12,24,64] {
            let m=e.tessellate(segments).unwrap();
            check_facets(&m.mesh,e.model(),1,segments);
            assert_close(m.analytic_volume,144.-4.*PI*radius*radius,1e-9);
        }
    }
}
#[test]
fn far_from_origin_stays_watertight_with_local_chart_coordinates(){
    let shift=1e8;
    let source=Solid::block(
        p(shift-3.,shift-1.,shift-3.),6.,4.,6.).unwrap();
    let mut e=BrepStructuralEditor::new(source).unwrap();
    {
        let mut tx=e.begin_transaction();
        tx.add_through_hole(params(shift-1.4,shift+0.2,0.7)).unwrap();
        tx.add_through_hole(params(shift+1.4,shift-0.2,0.9)).unwrap();
        tx.commit().unwrap();
    }
    let result=e.tessellate(48).unwrap();
    check_facets(&result.mesh,e.model(),2,48);
    assert_close(result.analytic_volume,
        144.-4.*PI*(0.49+0.81),1e-9);
    assert!((result.faceted_volume-result.analytic_volume).abs()<1.0);
}
