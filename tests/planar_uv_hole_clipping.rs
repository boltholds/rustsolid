//! Certified planar face trim with an outer line loop and one exact circular hole.
//! The hole is never represented by faceted approximation in the intersection kernel.
use rustsolid::{
    BrepBody, FaceId, GeometryTolerance, IntersectionGraphBuilder, GraphState,
    HoleParameters, Point3, Solid, Surface3, ThroughHoleBrep,
    TrimmedIntersectionEngine, TrimmedIntersection, TrimContact,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn eps()->GeometryTolerance{GeometryTolerance::default()}
fn solid_with_hole()->BrepBody{
    let original=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(&original,
        HoleParameters{center_x:0.,center_z:0.,radius:1.}).unwrap())
}
fn vertical_plane(x:f64)->(BrepBody,FaceId){
    let body=BrepBody::from(Solid::block(p(x,4.,-4.),2.,2.,8.).unwrap());
    let b=body.shared().unwrap();
    let fid=b.faces.iter().find_map(|f|{
        let Surface3::Plane(plane)=b.face_surface(f.id)? else{return None};
        ((plane.frame.origin.x-x).abs()<1e-10 &&
            plane.frame.normal.x.abs()>0.99).then_some(f.id)
    }).expect("expected left vertical plane face");
    (body,fid)
}
fn approx(a:f64,b:f64){assert!((a-b).abs()<1e-8,"{a} != {b}");}

#[test]
fn planar_annulus_line_is_clipped_into_two_true_nonhole_arcs(){
    let ring=solid_with_hole();
    let r=ring.shared().unwrap();
    assert_eq!(r.faces[1].loops.len(),2,"true inner loop must be present");
    let (plane,face)=vertical_plane(0.);
    let p=plane.shared().unwrap();
    let trim=TrimmedIntersectionEngine::new(eps()).unwrap()
        .intersect_faces(&r,FaceId(1),&p,face).unwrap();
    let TrimmedIntersection::Geometry{branches,contacts,..}=trim else{
        panic!("expected split crossing line by native circular inner trim: {trim:?}");
    };
    assert_eq!(branches.len(),2);
    assert!(contacts.is_empty());
    assert!(branches.iter().all(|b|b.contact==TrimContact::Crossing));
    let mut sections=branches.into_iter()
        .map(|arc|{
            for f in [0.0,0.13,0.5,0.87,1.0]{
                let q=arc.evaluate(f).unwrap();
                approx(q.x,0.);
                approx(q.y,5.);
                assert!(q.z.abs()>=1.0-1e-8);
            }
            let a=arc.start.z;let b=arc.end.z;
            [a.min(b),a.max(b)]
        }).collect::<Vec<_>>();
    sections.sort_by(|a,b|a[0].total_cmp(&b[0]));
    approx(sections[0][0],-3.);
    approx(sections[0][1],-1.);
    approx(sections[1][0],1.);
    approx(sections[1][1],3.);

    let graph=IntersectionGraphBuilder::new(eps()).unwrap()
        .intersect_faces(&r,FaceId(1),&p,face).unwrap();
    assert_eq!(graph.state,GraphState::Ready,"{graph:?}");
    assert_eq!(graph.arcs.len(),2);
    assert_eq!(graph.nodes.len(),4);
    graph.validate_against(&r,FaceId(1),&p,face).unwrap();
    // Two points belong to the physical hole edge, not a fake UV seam.
    let boundary=graph.nodes.iter().filter(|node|
        (node.point.z.abs()-1.0).abs()<1e-8).count();
    assert_eq!(boundary,2);
}
#[test]
fn outside_plane_has_no_intersection_and_tangent_is_not_accepted_as_cut(){
    let ring=solid_with_hole();
    let r=ring.shared().unwrap();
    let (plane,face)=vertical_plane(4.0);
    let p=plane.shared().unwrap();
    let result=TrimmedIntersectionEngine::new(eps()).unwrap()
        .intersect_faces(&r,FaceId(1),&p,face).unwrap();
    assert!(matches!(result,TrimmedIntersection::Empty));
    let (near,near_face)=vertical_plane(1.0+1e-10);
    let n=near.shared().unwrap();
    let result=TrimmedIntersectionEngine::new(eps()).unwrap()
        .intersect_faces(&r,FaceId(1),&n,near_face).unwrap();
    assert!(matches!(result,TrimmedIntersection::Ambiguous{..}),
        "near hole tangency must be held: {result:?}");
}
