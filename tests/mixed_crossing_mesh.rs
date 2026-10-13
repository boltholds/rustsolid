//! Independent faceted display certificates for analytic crossed
//! through/blind B-rep. The mesh is never an editable source of CAD truth.
use rustsolid::{
    BlindHoleParameters,HoleParameters,MixedCrossingBrep,
    MixedCrossingMesher,Point3,Solid,
};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn body()->MixedCrossingBrep{
    MixedCrossingBrep::from_box(&box_(),t(0.,0.,0.5),b(0.5,0.,0.5,1.5)).unwrap()
}
fn close(a:f64,b:f64,e:f64){assert!((a-b).abs()<e,"{a:.12} != {b:.12}");}

#[test]
fn certified_crossing_mesh_is_watertight_genus_one_with_face_owners(){
    let body=body();
    let d=MixedCrossingMesher::tessellate(&body,32).unwrap();
    assert_eq!(d.genus,1);
    assert_eq!(d.mesh.triangles.len(),d.mesh.triangle_faces.len());
    assert_eq!(d.circle_segments_per_arc,32);
    assert_eq!(d.mesh.vertices.len() as isize-d.triangle_edges as isize+
        d.mesh.triangles.len() as isize,0);
    let mut owned=[0usize;11];
    for f in &d.mesh.triangle_faces{owned[f.0 as usize]+=1;}
    assert!(owned.into_iter().all(|c|c>0),
        "each of 11 source B-rep faces contributes display facets");
    assert!(d.facet_volume.is_finite()&&d.facet_volume>0.0);
    close(d.facet_volume,d.expected_polygonal_volume,2e-6);
    close(d.analytic_volume,140.14094873829234,1e-10);
    assert!(d.facet_volume_error>= -2e-6);
}

#[test]
fn analytic_volume_is_unchanged_by_facet_density_and_facet_error_converges(){
    let body=body();let m=body.model().mass;
    let a=body.tessellate(16).unwrap();
    let b=body.tessellate(32).unwrap();
    let c=body.tessellate(64).unwrap();
    assert_eq!(body.model().mass,m);
    close(a.analytic_volume,m.volume,1e-12);
    close(b.analytic_volume,m.volume,1e-12);
    close(c.analytic_volume,m.volume,1e-12);
    assert!(a.facet_volume_error>b.facet_volume_error);
    assert!(b.facet_volume_error>c.facet_volume_error);
    assert!(c.facet_volume_error>= -1e-7);
    assert!(a.mesh.triangles.len()<b.mesh.triangles.len());
    assert!(b.mesh.triangles.len()<c.mesh.triangles.len());
}

#[test]
fn numerical_mesh_agrees_with_two_independent_polygonized_cap_sections(){
    let body=body();
    let d=body.tessellate(64).unwrap();
    let box_height=4.0;
    let upper_depth=1.5;
    let predicted=(box_height-upper_depth)*d.bottom_cap_polygon_area+
        upper_depth*d.top_cap_polygon_area;
    close(predicted,d.facet_volume,2e-6);
    assert!(d.floor_polygon_area>0.);
    let analytic_floor=PI*0.25-body.plan().footprint.overlap_area.unwrap();
    assert!(d.floor_polygon_area<analytic_floor+1e-6);
    assert!(d.facet_area>0.0);
}

#[test]
fn negative_and_excessive_segment_budgets_refuse_without_mutating_source(){
    let body=body();
    let before=(body.revision(),body.model().mass,body.undo_depth());
    for segments in [0,2,7,513,10000] {
        assert!(body.tessellate(segments).is_err());
        assert_eq!((body.revision(),body.model().mass,body.undo_depth()),before);
    }
}

#[test]
fn rebuild_after_parameter_edit_and_undo_keeps_watertight_certificates(){
    let mut body=body();
    let before=body.tessellate(32).unwrap().facet_volume;
    let mut tx=body.begin_edit();
    tx.replace_through(t(-0.1,0.,0.6)).unwrap();
    tx.replace_blind(b(0.5,0.,0.5,1.75)).unwrap();
    tx.commit().unwrap();
    let changed=body.tessellate(32).unwrap();
    assert_ne!(before,changed.facet_volume);
    assert!(changed.facet_volume_error>= -2e-6);
    assert!(body.undo().unwrap());
    let restored=body.tessellate(32).unwrap();
    close(restored.facet_volume,before,1e-10);
    assert!(body.redo().unwrap());
    let redone=body.tessellate(32).unwrap();
    close(redone.facet_volume,changed.facet_volume,1e-10);
}

#[test]
fn translated_world_coordinates_do_not_change_source_analytic_volume(){
    let base=body().tessellate(32).unwrap();
    let sx=10000.0;let sz=-10000.0;
    let source=Solid::block(p(-3.+sx,-1.,-3.+sz),6.,4.,6.).unwrap();
    let shifted=MixedCrossingBrep::from_box(&source,
        t(sx,sz,0.5),b(sx+0.5,sz,0.5,1.5)).unwrap();
    let shifted_display=shifted.tessellate(32).unwrap();
    close(shifted.model().mass.volume,base.analytic_volume,1e-9);
    close(shifted_display.facet_volume,base.facet_volume,2e-6);
}
