use rustsolid::{
    CoedgeId, Curve2, Curve3, CurvedFaceKind, CylindricalBrep, FaceId, Frame3,
    GeometryTolerance, Point3, Surface3,
};
use std::f64::consts::{PI, TAU};

fn cylinder() -> CylindricalBrep {
    CylindricalBrep::upright(Point3{x:0.0,y:0.0,z:0.0},2.0,5.0,
        GeometryTolerance::default()).unwrap()
}
fn approx(a:f64,b:f64) {assert!((a-b).abs()<=1e-9*a.abs().max(b.abs()).max(1.0),"{a} vs {b}");}
fn near(a:Point3,b:Point3) {approx(a.x,b.x);approx(a.y,b.y);approx(a.z,b.z);}

#[test]
fn native_closed_cylinder_has_three_analytic_faces_and_euler_two(){
    let solid=cylinder();
    assert_eq!((solid.vertices.len(),solid.edges.len(),solid.coedges.len(),solid.faces.len(),solid.loops.len(),solid.shells.len()),(2,3,6,3,3,1));
    assert_eq!(solid.euler_characteristic(),2);
    assert_eq!(solid.faces.iter().map(|f|f.kind).collect::<Vec<_>>(),
        vec![CurvedFaceKind::BottomCap,CurvedFaceKind::TopCap,CurvedFaceKind::CylinderSide]);
    assert!(matches!(solid.face_surface(FaceId(2)),Some(Surface3::Cylinder(_))));
    assert!(matches!(solid.edge_curve(solid.seam_edge()),Some(Curve3::Line(_))));
    assert!(matches!(solid.edge_curve(rustsolid::EdgeId(0)),Some(Curve3::Circle(_))));
    approx(solid.mass.volume,PI*4.0*5.0);
    approx(solid.mass.surface_area,8.0*PI+TAU*2.0*5.0);
    near(solid.mass.centroid,Point3{x:0.0,y:2.5,z:0.0});
    near(solid.bbox.min,Point3{x:-2.0,y:0.0,z:-2.0});
    near(solid.bbox.max,Point3{x:2.0,y:5.0,z:2.0});
    solid.validate().unwrap();
}

#[test]
fn periodic_seam_has_two_distinct_pcurves_on_same_face() {
    let s=cylinder();
    let a=&s.coedges[4];let b=&s.coedges[5];
    assert_eq!(a.edge,b.edge);assert_eq!(a.face,FaceId(2));assert_eq!(b.face,FaceId(2));
    assert_eq!(a.twin,CoedgeId(5));assert_eq!(b.twin,CoedgeId(4));
    let first=s.trim_curve(a.id).unwrap();let second=s.trim_curve(b.id).unwrap();
    let Curve2::Line(a_uv)=first else {panic!("expected seam line")};
    let Curve2::Line(b_uv)=second else {panic!("expected seam line")};
    approx(a_uv.origin.u,TAU);approx(b_uv.origin.u,0.0);
    let p=s.side_surface().evaluate(0.0,2.5).unwrap().point;
    let q=s.side_surface().evaluate(TAU,2.5).unwrap().point;
    near(p,q);
}

#[test]
fn cap_loops_are_single_closed_circle_coedges() {
    let s=cylinder();
    for face_index in 0..2 {
        let lp=&s.loops[face_index];
        let c=&s.coedges[lp.first_coedge.0 as usize];
        assert_eq!(c.next,c.id);assert_eq!(c.prev,c.id);
        assert_eq!(c.face,FaceId(face_index as u32));
        assert!(matches!(s.trim_curve(c.id),Some(Curve2::Circle(_))));
    }
}

#[test]
fn facets_are_watertight_keep_source_face_and_are_resolution_independent() {
    let s=cylinder();
    let coarse=s.tessellate(8).unwrap();
    let fine=s.tessellate(128).unwrap();
    assert_eq!((coarse.vertices.len(),coarse.triangles.len()),(18,32));
    assert_eq!((fine.vertices.len(),fine.triangles.len()),(258,512));
    for f in 0..3 {
        assert!(fine.triangle_faces.iter().any(|id|id.0==f));
    }
    approx(s.mass.volume,PI*4.0*5.0);
    s.check_mesh(&coarse).unwrap();s.check_mesh(&fine).unwrap();
}

#[test]
fn arbitrary_orientation_retains_analytic_mass_and_exact_bbox() {
    let tol=GeometryTolerance::default();
    let origin=Point3{x:6.0,y:-5.0,z:7.0};
    let frame=Frame3::from_axes(origin,Point3{x:0.0,y:1.0,z:0.0},
        Point3{x:0.0,y:0.0,z:1.0},tol).unwrap();
    let s=CylindricalBrep::new(frame,2.0,5.0,tol).unwrap();
    // Axis is +X, so only X extends by height.
    near(s.mass.centroid,Point3{x:8.5,y:-5.0,z:7.0});
    near(s.bbox.min,Point3{x:6.0,y:-7.0,z:5.0});
    near(s.bbox.max,Point3{x:11.0,y:-3.0,z:9.0});
    s.tessellate(32).unwrap();
}

#[test]
fn malformed_seam_twin_and_trim_are_rejected() {
    let mut s=cylinder();
    s.coedges[4].twin=CoedgeId(0);
    assert!(s.validate().is_err());
    let mut s=cylinder();
    s.curves2[4]=s.curves2[5];
    assert!(s.validate().is_err());
}

#[test]
fn invalid_radius_height_frame_and_subresolution_geometry_fail() {
    let tol=GeometryTolerance::default();
    let origin=Point3{x:0.0,y:0.0,z:0.0};
    assert!(CylindricalBrep::upright(origin,0.0,1.0,tol).is_err());
    assert!(CylindricalBrep::upright(origin,1.0,f64::INFINITY,tol).is_err());
    assert!(CylindricalBrep::upright(origin,1e-12,10.0,tol).is_err());
    assert!(cylinder().tessellate(3).is_err());
}

#[test]
fn normal_and_derivatives_are_analytic_at_non_facet_angle() {
    let s=cylinder();
    let d=s.side_surface().evaluate(0.37,1.9).unwrap();
    let n=Point3{x:d.du.y*d.dv.z-d.du.z*d.dv.y,y:d.du.z*d.dv.x-d.du.x*d.dv.z,z:d.du.x*d.dv.y-d.du.y*d.dv.x};
    let l=(n.x*n.x+n.y*n.y+n.z*n.z).sqrt();
    let normal=Point3{x:n.x/l,y:n.y/l,z:n.z/l};
    let expected=Point3{x:0.37f64.cos(),y:0.0,z:-0.37f64.sin()};
    near(normal,expected);
}
