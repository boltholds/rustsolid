//! Independent watertight facet certificates for the v0.31 crossing source.
//! The source B-rep and analytic mass are never reconstructed from triangles.
use std::collections::BTreeMap;
use rustsolid::{BlindHoleParameters,FaceId,GeometryTolerance,HoleParameters,Point3,
    Solid,ThroughBlindCrossingBrep};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn source()->Solid {Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn source_model()->ThroughBlindCrossingBrep {
    ThroughBlindCrossingBrep::from_box(&source(),t(0.,0.,0.5),b(0.5,0.,0.5,1.5)).unwrap()
}
fn check_watertight(triangles:&[[u32;3]],positions:&[Point3]){
    let mut pairs=BTreeMap::<(u32,u32),Vec<bool>>::new();
    for tri in triangles{
        assert!(tri.iter().all(|v|(*v as usize)<positions.len()));
        for (u,v) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]{pairs.entry(
            (u.min(v),u.max(v))).or_default().push(u<v);
        }
    }
    assert!(pairs.values().all(|d|d.len()==2 && d[0]!=d[1]));
}
#[test]
fn crossed_through_blind_emits_real_watertight_display_with_owner_face_ids(){
    let body=source_model();
    body.validate().unwrap();
    for n in [12,24,48]{
        let display=body.tessellate(n).unwrap();
        assert_eq!(display.segments_per_arc,n);
        assert!(display.top_triangles>0&&display.bottom_triangles>0&&display.crescent_triangles>0);
        assert_eq!(display.mesh.triangles.len(),display.mesh.triangle_faces.len());
        check_watertight(&display.mesh.triangles,&display.mesh.vertices);
        assert!(display.mesh.triangle_faces.contains(&FaceId(8)));
        assert!(display.mesh.triangle_faces.contains(&FaceId(9)));
        for f in &display.mesh.triangle_faces {assert!((f.0 as usize)<body.model().faces.len());}
        assert_eq!(body.model().summary().genus,1);
        assert_eq!(body.model().summary().euler_characteristic,0);
        assert_eq!(body.revision(),0);
        assert_eq!(body.undo_depth(),0);
        assert!(display.faceted_volume>0.0);
        assert!(display.faceted_surface_area>0.0);
        assert!((display.analytic_volume-140.14094873829234).abs()<1e-10);
    }
}

#[test]
fn faceted_mass_converges_toward_analytic_crossed_volume(){
    let body=source_model();
    let coarse=body.tessellate(12).unwrap();
    let medium=body.tessellate(48).unwrap();
    let fine=body.tessellate(192).unwrap();
    assert!(fine.faceting_volume_error.abs()<medium.faceting_volume_error.abs());
    assert!(medium.faceting_volume_error.abs()<coarse.faceting_volume_error.abs());
    assert!(fine.faceting_volume_error.abs()<0.001);
    assert_eq!(body.model().mass.volume,fine.analytic_volume);
}

#[test]
fn display_rebuild_is_deterministic_and_follows_committed_source_edits(){
    let mut shape=source_model();
    let a=shape.tessellate(24).unwrap();
    let repeat=shape.tessellate(24).unwrap();
    assert_eq!(a.mesh.vertices,repeat.mesh.vertices);
    assert_eq!(a.mesh.triangles,repeat.mesh.triangles);
    assert_eq!(a.mesh.triangle_faces,repeat.mesh.triangle_faces);
    {
        let mut tx=shape.begin_edit();
        tx.replace_blind(b(0.6,0.,0.5,2.0)).unwrap();
        tx.commit().unwrap();
    }
    shape.validate().unwrap();
    let changed=shape.tessellate(24).unwrap();
    assert_ne!(a.analytic_volume,changed.analytic_volume);
    assert!(shape.undo().unwrap());
    let reverted=shape.tessellate(24).unwrap();
    assert_eq!(a.mesh.vertices,reverted.mesh.vertices);
    assert_eq!(a.mesh.triangles,reverted.mesh.triangles);
}

#[test]
fn invalid_segment_budget_and_near_tangent_keep_current_source_unchanged(){
    let mut body=source_model();
    let old=body.model().mass;
    let old_revision=body.revision();
    for seg in [0,4,11,513,4096] {assert!(body.tessellate(seg).is_err());}
    {
        let mut tx=body.begin_edit();
        tx.replace_blind(b(1.0,0.,0.5,1.5)).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(body.model().mass,old);
    assert_eq!(body.revision(),old_revision);
    body.validate().unwrap();
}

#[test]
fn translation_is_geometrically_safe_without_global_origin_cancellation(){
    let s=Solid::block_with_tolerance(p(1e6-3.,-1.,-1e6-3.),6.,4.,6.,
        GeometryTolerance::default()).unwrap();
    let body=ThroughBlindCrossingBrep::from_box(&s,
        t(1e6,-1e6,0.5),b(1e6+0.5,-1e6,0.5,1.5)).unwrap();
    body.validate().unwrap();
    let d=body.tessellate(48).unwrap();
    assert!((d.analytic_volume-140.14094873829234).abs()<1e-8);
    check_watertight(&d.mesh.triangles,&d.mesh.vertices);
}
