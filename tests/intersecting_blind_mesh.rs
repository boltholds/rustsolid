//! Independent geometry.v1 display certificates for a truly intersecting
//! two-circular-arc inner UV loop with distinct cylinder side patches and
//! a single nonconvex planar union floor.
use rustsolid::{
    BlindHoleParameters, FaceId, IntersectingBlindBrep, Mesh, Point3, Solid,
};
use std::collections::BTreeMap;
use std::f64::consts::PI;

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn h(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn shape()->IntersectingBlindBrep{
    IntersectingBlindBrep::from_box(&block(),[
        h(-0.7,0.,1.,2.),h(0.7,0.,1.,2.)
    ]).unwrap()
}
fn eq(a:f64,b:f64,eps:f64){assert!((a-b).abs()<=eps,"{a} vs {b}");}
fn watertight(mesh:&Mesh,n:usize){
    assert_eq!(mesh.vertices.len(),8+4*n);
    assert_eq!(mesh.triangles.len(),8*n+12);
    assert_eq!(mesh.triangle_faces.len(),mesh.triangles.len());
    let mut uses=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut faces=BTreeMap::<u32,usize>::new();
    for (tri,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        *faces.entry(face.0).or_default()+=1;
        let [a,b,c]=*tri;
        for (i,j) in [(a,b),(b,c),(c,a)]{
            uses.entry((i.min(j),i.max(j))).or_default().push(i<j);
        }
    }
    assert!(uses.values().all(|v|v.len()==2 && v[0]!=v[1]));
    assert_eq!(faces[&0],2);
    assert_eq!(faces[&1],4+2*n);
    for id in 2..6{assert_eq!(faces[&id],2);}
    assert_eq!(faces[&6],2*n);
    assert_eq!(faces[&7],2*n);
    assert_eq!(faces[&8],2*n-2);
}
#[test]
fn two_intersecting_cylinders_have_single_watertight_union_floor(){
    let s=shape();
    let mesh=s.tessellate(48).unwrap();
    watertight(&mesh.mesh,48);
    assert_eq!(mesh.segments_per_arc,48);
    assert_eq!(mesh.opening_triangles,100);
    assert_eq!(mesh.floor_triangles,94);
    eq(mesh.analytic_volume,s.model().mass.volume,1e-9);
    assert!(mesh.faceting_volume_error>0.);
    assert!(mesh.faceted_opening_area>0.);
    eq(mesh.faceted_opening_area+mesh.faceted_floor_area,36.,1e-8);
}
#[test]
fn faceted_volume_converges_without_changing_exact_analytic_mass(){
    let s=shape();
    let original=s.model().mass;
    let coarse=s.tessellate(12).unwrap();
    let medium=s.tessellate(48).unwrap();
    let fine=s.tessellate(192).unwrap();
    watertight(&coarse.mesh,12);
    watertight(&medium.mesh,48);
    watertight(&fine.mesh,192);
    assert!(coarse.faceting_volume_error>medium.faceting_volume_error);
    assert!(medium.faceting_volume_error>fine.faceting_volume_error);
    assert!(fine.faceting_volume_error<medium.faceting_volume_error*0.09);
    assert_eq!(s.model().mass,original);
    assert_eq!(s.revision(),0);
}
#[test]
fn crossing_mesh_respects_source_arc_carriers_after_parametric_commit(){
    let mut s=shape();
    let mut tx=s.begin_edit();
    tx.replace(0,h(-0.8,-0.1,1.1,2.5)).unwrap();
    tx.replace(1,h(0.7,0.1,0.95,2.5)).unwrap();
    tx.commit().unwrap();
    let mesh=s.tessellate(32).unwrap();
    watertight(&mesh.mesh,32);
    assert!(mesh.faceted_volume>=mesh.analytic_volume);
    assert_eq!(s.model().summary().genus,0);
    assert!(s.undo().unwrap());
    watertight(&s.tessellate(32).unwrap().mesh,32);
    assert!(s.redo().unwrap());
    watertight(&s.tessellate(32).unwrap().mesh,32);
}
#[test]
fn rotated_unequal_circle_crossing_has_closed_nonconvex_floor(){
    let s=IntersectingBlindBrep::from_box(&block(),[
        h(-0.65,-0.55,1.25,2.3),h(0.8,0.5,0.95,2.3)
    ]).unwrap();
    let mesh=s.tessellate(64).unwrap();
    watertight(&mesh.mesh,64);
    eq(mesh.analytic_volume,
        block().mass.volume-s.region().disk_union_area*2.3,1e-9);
}
#[test]
fn invalid_display_segment_count_is_fail_closed_and_revision_stable(){
    let s=shape();
    for n in [0,2,8,11,513,8192]{
        assert!(s.tessellate(n).is_err());
    }
    assert_eq!(s.revision(),0);
    s.validate().unwrap();
}
#[test]
fn sampling_resolution_never_changes_true_exposed_analytical_perimeter(){
    let s=shape();
    let r=s.region().exposed_perimeter;
    let true_surface=block().mass.surface_area+r*2.;
    for segments in [12,48,192]{
        let mesh=s.tessellate(segments).unwrap();
        assert!(mesh.faceted_surface_area<true_surface);
        eq(s.model().mass.surface_area,true_surface,1e-9);
    }
}
#[test]
fn near_world_translation_remains_local_chart_stable(){
    let shift=1e4;
    let source=Solid::block(p(shift-3.,shift-1.,shift-3.),6.,4.,6.).unwrap();
    let s=IntersectingBlindBrep::from_box(&source,[
        h(shift-0.7,shift,1.,2.),h(shift+0.7,shift,1.,2.)
    ]).unwrap();
    let mesh=s.tessellate(48).unwrap();
    watertight(&mesh.mesh,48);
    eq(mesh.analytic_volume,s.model().mass.volume,1e-8);
}
