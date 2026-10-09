use rustsolid::{FaceRole, GeometryError, Point2, Point3, Solid};

fn square_ccw() -> Vec<Point2> {
    vec![Point2{x:0.0,z:0.0}, Point2{x:2.0,z:0.0},
         Point2{x:2.0,z:4.0}, Point2{x:0.0,z:4.0}]
}
fn approx(a: f64,b:f64) { assert!((a-b).abs() <= 1e-9, "{a} != {b}"); }

#[test]
fn box_is_closed_and_has_analytic_mass_properties() {
    let body=Solid::block(Point3{x:1.0,y:2.0,z:3.0},2.0,3.0,4.0).unwrap();
    assert_eq!((body.vertices.len(),body.edges.len(),body.faces.len()),(8,12,6));
    assert_eq!(body.euler_characteristic(),2);
    assert_eq!(body.mesh.triangles.len(),12);
    approx(body.mass.volume,24.0);
    approx(body.mass.surface_area,52.0);
    assert_eq!(body.mass.centroid,Point3{x:2.0,y:3.5,z:5.0});
    assert_eq!(body.bbox.min,Point3{x:1.0,y:2.0,z:3.0});
    assert_eq!(body.bbox.max,Point3{x:3.0,y:5.0,z:7.0});
    body.validate().unwrap();
}

#[test]
fn concave_extrusion_triangulates_caps_without_holes() {
    // Area 5; perimeter 12; height 2 -> volume 10, area 34.
    let polygon=vec![Point2{x:0.0,z:0.0},Point2{x:3.0,z:0.0},
                     Point2{x:3.0,z:1.0},Point2{x:1.0,z:1.0},
                     Point2{x:1.0,z:3.0},Point2{x:0.0,z:3.0}];
    let body=Solid::extrude_xz(&polygon,2.0).unwrap();
    approx(body.mass.volume,10.0);
    approx(body.mass.surface_area,34.0);
    assert_eq!((body.vertices.len(),body.edges.len(),body.faces.len()),(12,18,8));
    assert_eq!(body.mesh.triangles.len(),20);
    body.validate().unwrap();
}

#[test]
fn clockwise_and_explicitly_closed_profiles_have_equivalent_geometry() {
    let ccw=Solid::extrude_xz(&square_ccw(),3.0).unwrap();
    let mut reverse=square_ccw(); reverse.reverse(); reverse.push(reverse[0]);
    let cw=Solid::extrude_xz(&reverse,3.0).unwrap();
    approx(ccw.mass.volume,cw.mass.volume);
    approx(ccw.mass.surface_area,cw.mass.surface_area);
    assert_eq!(cw.euler_characteristic(),2);
}

#[test]
fn face_hit_returns_stable_topological_face() {
    let body=Solid::extrude_xz(&square_ccw(),5.0).unwrap();
    let hit=body.ray_cast(Point3{x:1.0,y:10.0,z:2.0},Point3{x:0.0,y:-1.0,z:0.0}).unwrap();
    assert_eq!(body.faces[hit.face_id.0 as usize].role,FaceRole::TopCap);
    approx(hit.distance,5.0);
    assert_eq!(hit.normal,Point3{x:0.0,y:1.0,z:0.0});
    let non_unit_ray=body.ray_cast(Point3{x:1.0,y:10.0,z:2.0},Point3{x:0.0,y:-2.0,z:0.0}).unwrap();
    approx(non_unit_ray.distance,5.0);
}

#[test]
fn rejects_self_intersecting_and_zero_area_profiles() {
    let bow_tie=vec![Point2{x:0.0,z:0.0},Point2{x:1.0,z:1.0},
                     Point2{x:0.0,z:1.0},Point2{x:1.0,z:0.0}];
    assert!(matches!(Solid::extrude_xz(&bow_tie,2.0),Err(GeometryError::InvalidProfile(_))));
    assert!(Solid::extrude_xz(&square_ccw(),0.0).is_err());
    assert!(Solid::extrude_xz(&square_ccw(),f64::INFINITY).is_err());
    assert!(Solid::block(Point3{x:0.0,y:0.0,z:0.0},2.0,-1.0,4.0).is_err());
}

#[test]
fn removing_a_collinear_control_point_preserves_solid() {
    let mut p=square_ccw();
    p.insert(1,Point2{x:1.0,z:0.0});
    let body=Solid::extrude_xz(&p,3.0).unwrap();
    assert_eq!(body.vertices.len(),8);
    approx(body.mass.volume,24.0);
}

#[test]
fn centroid_is_stable_far_from_the_origin() {
    // A 1x1 extrusion offset by 10^12 units previously suffered cancellation
    // when polygon centroid was computed from absolute-coordinate moments.
    let x = 1e12;
    let z = -1e12;
    let profile = [Point2 { x, z }, Point2 { x: x + 1.0, z },
                   Point2 { x: x + 1.0, z: z + 1.0 }, Point2 { x, z: z + 1.0 }];
    let solid = Solid::extrude_xz(&profile, 2.0).unwrap();
    approx(solid.mass.volume, 2.0);
    approx(solid.mass.centroid.x, x + 0.5);
    approx(solid.mass.centroid.z, z + 0.5);
    approx(solid.mass.centroid.y, 1.0);
}

#[test]
fn ray_miss_and_zero_direction_are_safe() {
    let body = Solid::block(Point3 { x: 0.0, y: 0.0, z: 0.0 }, 2.0, 3.0, 4.0).unwrap();
    assert!(body.ray_cast(Point3 { x: 99.0, y: 99.0, z: 99.0 },
                          Point3 { x: 0.0, y: -1.0, z: 0.0 }).is_none());
    assert!(body.ray_cast(Point3 { x: 1.0, y: 10.0, z: 2.0 },
                          Point3 { x: 0.0, y: 0.0, z: 0.0 }).is_none());
}

#[test]
fn non_finite_and_duplicate_vertices_are_rejected() {
    let profile = [Point2 { x: 0.0, z: 0.0 }, Point2 { x: f64::NAN, z: 0.0 },
                   Point2 { x: 1.0, z: 1.0 }];
    assert!(Solid::extrude_xz(&profile, 2.0).is_err());
    let profile = [Point2 { x: 0.0, z: 0.0 }, Point2 { x: 1.0, z: 0.0 },
                   Point2 { x: 1.0, z: 0.0 }, Point2 { x: 1.0, z: 1.0 },
                   Point2 { x: 0.0, z: 1.0 }];
    assert!(Solid::extrude_xz(&profile, 2.0).is_err());
}
