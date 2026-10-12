//! Independent display correctness for multi-component analytical blind
//! regions (single circular floors and two-arc nonconvex union floors).
use rustsolid::{BlindHoleParameters,CompoundBlindBrep,BlindRegionKind,
    Point3,Solid,FaceId,Mesh};
use std::collections::BTreeMap;
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
fn approx(a:f64,b:f64,eps:f64){assert!((a-b).abs()<eps,"{a} vs {b}");}
fn closure(mesh:&Mesh){
    assert_eq!(mesh.triangles.len(),mesh.triangle_faces.len());
    let mut uses=BTreeMap::<(u32,u32),Vec<bool>>::new();
    for t in &mesh.triangles{
        let [a,b,c]=*t;
        assert_ne!(a,b);assert_ne!(a,c);assert_ne!(b,c);
        for (x,y) in [(a,b),(b,c),(c,a)]{
            uses.entry((x.min(y),x.max(y))).or_default().push(x<y);
        }
    }
    assert!(uses.values().all(|e|e.len()==2&&e[0]!=e[1]),
        "each display edge must have opposite use by exactly two faces");
}
#[test]
fn mixed_pair_and_single_has_two_closed_display_floors_and_one_top_face(){
    let b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let d=b.tessellate(48).unwrap();
    closure(&d.mesh);
    assert_eq!(d.segments_per_circle,48);
    assert_eq!(d.top_face_loop_count,3);
    assert_eq!(d.floor_regions.len(),2);
    assert_eq!(d.floor_regions[0].loop_vertices,96);
    assert_eq!(d.floor_regions[1].loop_vertices,48);
    assert_eq!(d.floor_regions[0].floor_face,FaceId(8));
    assert_eq!(d.floor_regions[1].floor_face,FaceId(10));
    assert_eq!(d.floor_regions[0].floor_triangles,94);
    assert_eq!(d.floor_regions[1].floor_triangles,46);
    assert_eq!(d.top_face_triangles,150);
    assert_eq!(d.mesh.vertices.len(),8+6*48);
    assert_eq!(d.mesh.triangles.len(),12*48+10);
    assert_eq!(d.mesh.triangle_faces.iter().filter(|f|f.0==8).count(),94);
    assert_eq!(d.mesh.triangle_faces.iter().filter(|f|f.0==10).count(),46);
    assert_eq!(d.mesh.triangle_faces.iter().filter(|f|f.0==1).count(),150);
    assert!(d.faceting_volume_error>0.);
    approx(d.analytic_volume,b.model().mass.volume,1e-9);
}
#[test]
fn two_crossing_pairs_are_closed_with_two_independent_compound_floors(){
    let b=CompoundBlindBrep::from_box(&block(),vec![
        h(-1.75,-1.2,0.6,1.3),h(-0.7,-1.2,0.6,1.3),
        h(0.7,1.2,0.6,2.6),h(1.75,1.2,0.6,2.6)]).unwrap();
    let display=b.tessellate(48).unwrap();
    closure(&display.mesh);
    assert_eq!(display.top_face_loop_count,3);
    assert_eq!(display.top_face_triangles,198);
    assert_eq!(display.mesh.vertices.len(),8+8*48);
    assert_eq!(display.mesh.triangles.len(),16*48+10);
    assert!(display.floor_regions.iter().all(|f|f.loop_vertices==96));
    assert!(display.faceting_volume_error>0.);
}
#[test]
fn polygonal_volume_error_converges_to_unchanged_analytical_mass(){
    let b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let analytic=b.model().mass;
    let a=b.tessellate(12).unwrap();
    let medium=b.tessellate(48).unwrap();
    let fine=b.tessellate(192).unwrap();
    for d in [&a,&medium,&fine]{
        closure(&d.mesh);
        approx(d.analytic_volume,analytic.volume,1e-9);
        assert!(d.faceted_volume>=d.analytic_volume);
        assert!(d.faceted_surface_area<b.model().mass.surface_area);
    }
    assert!(a.faceting_volume_error>medium.faceting_volume_error);
    assert!(medium.faceting_volume_error>fine.faceting_volume_error);
    assert!(fine.faceting_volume_error<medium.faceting_volume_error*0.09);
    assert_eq!(b.model().mass,analytic);
}
#[test]
fn disjoint_and_crossed_component_permutation_preserves_polygonal_area_and_volume(){
    let a=mixed();
    let first=CompoundBlindBrep::from_box(&block(),a.clone()).unwrap();
    let second=CompoundBlindBrep::from_box(&block(),
        vec![a[2],a[1],a[0]]).unwrap();
    let x=first.tessellate(48).unwrap();
    let y=second.tessellate(48).unwrap();
    closure(&x.mesh);closure(&y.mesh);
    approx(x.faceted_volume,y.faceted_volume,1e-8);
    approx(x.faceted_surface_area,y.faceted_surface_area,1e-8);
    assert_eq!(x.top_face_loop_count,y.top_face_loop_count);
    assert_eq!(x.mesh.triangles.len(),y.mesh.triangles.len());
}
#[test]
fn atomic_group_repartition_rebuilds_entire_mesh_and_restores_on_undo(){
    let mut b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let original=b.tessellate(24).unwrap();
    let mut tx=b.begin_edit();
    tx.replace(1,h(-0.4,-0.5,0.25,1.8)).unwrap();
    tx.commit().unwrap();
    assert_eq!(b.components().len(),3);
    assert!(b.components().iter().all(|c|c.kind==BlindRegionKind::SingleCircle));
    let display=b.tessellate(24).unwrap();
    closure(&display.mesh);
    assert_eq!(display.top_face_loop_count,4);
    assert_eq!(display.floor_regions.len(),3);
    assert!(b.undo().unwrap());
    let restored=b.tessellate(24).unwrap();
    closure(&restored.mesh);
    approx(restored.faceted_volume,original.faceted_volume,1e-9);
    assert_eq!(restored.mesh.triangle_faces,original.mesh.triangle_faces);
    assert!(b.redo().unwrap());
    closure(&b.tessellate(24).unwrap().mesh);
}
#[test]
fn invalid_segment_count_or_sampling_budget_never_changes_brep_state(){
    let b=CompoundBlindBrep::from_box(&block(),mixed()).unwrap();
    let original=b.model().mass;
    for segments in [0,1,11,513,8192]{
        assert!(b.tessellate(segments).is_err());
    }
    assert_eq!(b.revision(),0);
    assert_eq!(b.model().mass,original);
    b.validate().unwrap();
}
#[test]
fn moderately_shifted_world_coordinates_keep_local_faceted_oracle(){
    let shift=1e4;
    let box_=Solid::block(p(shift-3.,shift-1.,shift-3.),6.,4.,6.).unwrap();
    let hs=vec![h(shift-1.5,shift-0.5,0.75,1.8),
        h(shift-0.4,shift-0.5,0.7,1.8),
        h(shift+1.7,shift+1.2,0.55,2.6)];
    let b=CompoundBlindBrep::from_box(&box_,hs).unwrap();
    let d=b.tessellate(48).unwrap();
    closure(&d.mesh);
    approx(d.analytic_volume,b.model().mass.volume,1e-8);
    assert!(d.faceting_volume_error>0.);
}
