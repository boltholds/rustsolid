use rustsolid::{
    face_support_plane, intersect_segment_plane, orient2d, orient3d, GeometryTolerance,
    Orientation, Plane3, PlaneIntersection, Point3, SegmentPlaneIntersection, Solid, FaceId,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn plane_z(z:f64)->Plane3 {
    Plane3::through(p(0.0,0.0,z),p(1.0,0.0,z),p(0.0,1.0,z),GeometryTolerance::default()).unwrap()
}
fn approx(a:f64,b:f64){assert!((a-b).abs()<1e-8,"{a} vs {b}");}

#[test] fn exact_orientation_signs_and_input_validation(){
    assert_eq!(orient2d([0.0,0.0],[1.0,0.0],[0.0,1.0]).unwrap(),Orientation::Positive);
    assert_eq!(orient2d([0.0,0.0],[0.0,1.0],[1.0,0.0]).unwrap(),Orientation::Negative);
    assert_eq!(orient2d([1e12,1e12],[1e12+1.0,1e12],[1e12+2.0,1e12]).unwrap(),Orientation::Zero);
    let value=orient3d(p(0.0,0.0,0.0),p(1.0,0.0,0.0),p(0.0,1.0,0.0),p(0.0,0.0,1.0)).unwrap();
    assert_ne!(value,Orientation::Zero);
    assert_eq!(orient3d(p(0.0,0.0,0.0),p(1.0,0.0,0.0),p(0.0,1.0,0.0),p(0.5,0.5,0.0)).unwrap(),Orientation::Zero);
    assert!(orient2d([0.0,0.0],[f64::INFINITY,1.0],[2.0,3.0]).is_err());
}
#[test] fn plane_relations_are_explicit(){
    let tol=GeometryTolerance::default();
    assert_eq!(plane_z(0.0).intersect(plane_z(3.0),tol).unwrap(),PlaneIntersection::Parallel);
    assert_eq!(plane_z(0.0).intersect(plane_z(0.0),tol).unwrap(),PlaneIntersection::Coincident);
    let vertical=Plane3::through(p(0.0,0.0,0.0),p(0.0,1.0,0.0),p(0.0,0.0,1.0),tol).unwrap();
    match plane_z(0.0).intersect(vertical,tol).unwrap() {
        PlaneIntersection::Line{point,direction}=>{
            approx(point.x,0.0);approx(point.z,0.0);approx(direction.y.abs(),1.0);
        }
        other=>panic!("expected line {other:?}"),
    }
}
#[test] fn segment_intersection_classifies_crossing_touching_and_coplanarity(){
    let tol=GeometryTolerance::default();let plane=plane_z(0.0);
    match intersect_segment_plane(p(2.0,3.0,-2.0),p(2.0,3.0,2.0),plane,tol).unwrap(){
        SegmentPlaneIntersection::Point{point,parameter}=>{
            approx(point.z,0.0);approx(parameter,0.5);
        },x=>panic!("{x:?}"),
    }
    assert_eq!(intersect_segment_plane(p(2.0,1.0,1.0),p(3.0,1.0,2.0),plane,tol).unwrap(),SegmentPlaneIntersection::Disjoint);
    assert_eq!(intersect_segment_plane(p(2.0,1.0,0.0),p(3.0,1.0,0.0),plane,tol).unwrap(),SegmentPlaneIntersection::Coplanar);
    assert_eq!(intersect_segment_plane(p(2.0,1.0,0.0),p(3.0,1.0,2.0),plane,tol).unwrap(),SegmentPlaneIntersection::Point{point:p(2.0,1.0,0.0),parameter:0.0});
}
#[test] fn face_supports_intersect_without_claiming_trimmed_intersection(){
    let solid=Solid::block(p(0.0,0.0,0.0),2.0,3.0,4.0).unwrap();
    let a=face_support_plane(&solid,FaceId(0)).unwrap();
    let b=face_support_plane(&solid,FaceId(2)).unwrap();
    assert!(matches!(a.intersect(b,solid.tolerance).unwrap(),PlaneIntersection::Line{..}));
}
