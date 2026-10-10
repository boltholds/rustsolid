use rustsolid::{
    intersect_segment_plane, GeometryTolerance, Orientation, Plane3, PlaneIntersection,
    Point3, SegmentPlaneIntersection, orient2d,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
#[test]
fn high_coordinate_small_plane_is_local_and_translation_stable(){
    let tolerance=GeometryTolerance::default();
    let x=1.0e12;
    let a=Plane3::through(p(x,x,x),p(x+1.0,x,x),p(x,x+1.0,x),tolerance).unwrap();
    let b=Plane3::through(p(x,x,x),p(x,x+1.0,x),p(x,x,x+1.0),tolerance).unwrap();
    match a.intersect(b,tolerance).unwrap(){
        PlaneIntersection::Line{point,direction}=>{
            assert_eq!(point.x,x);
            assert_eq!(point.z,x);
            assert!((direction.y.abs()-1.0).abs()<1e-10);
        }
        relation=>panic!("{relation:?}"),
    }
    assert_eq!(orient2d([x,x],[x+1.0,x],[x,x+1.0]).unwrap(),Orientation::Positive);
}
#[test]
fn invalid_planes_are_rejected(){
    let tol=GeometryTolerance::default();
    assert!(Plane3::through(p(0.0,0.0,0.0),p(1.0,0.0,0.0),p(2.0,0.0,0.0),tol).is_err());
    assert!(Plane3::through(p(0.0,0.0,0.0),p(f64::INFINITY,0.0,0.0),p(0.0,1.0,0.0),tol).is_err());
    let plane=Plane3::through(p(0.0,0.0,0.0),p(1.0,0.0,0.0),p(0.0,1.0,0.0),tol).unwrap();
    assert!(intersect_segment_plane(p(0.0,0.0,1.0),p(0.0,0.0,1.0),plane,tol).is_err());
    assert_eq!(intersect_segment_plane(p(0.0,0.0,1.0),p(0.0,0.0,2.0),plane,tol).unwrap(),SegmentPlaneIntersection::Disjoint);
}
