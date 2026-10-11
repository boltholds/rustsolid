//! Independent first-party contracts for actual bounded face intersections.
//! These are analytic B-rep models; no triangulation is used as a truth source.
use rustsolid::{
    BrepModel, CylindricalBrep, FaceId, GeometryTolerance, ModelUnits, Point3,
    Solid, Surface3, TrimAmbiguity, TrimContact, TrimmedIntersection,
    TrimmedIntersectionEngine,
};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::f64::consts::{PI,TAU};

const FIXTURES:&str=include_str!("fixtures/trimmed_intersections.v1.json");
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Corpus{schema:String,title:String,cases:Vec<Fixture>}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture{id:String,scenario:String,expected:String,branches:usize,points:usize}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block(origin:[f64;3],dims:[f64;3])->BrepModel{
    BrepModel::from_polyhedral(
        &Solid::block(p(origin[0],origin[1],origin[2]),dims[0],dims[1],dims[2]).unwrap(),
        ModelUnits::default()).unwrap()
}
fn cyl(origin:[f64;3])->BrepModel{
    BrepModel::from_cylinder(&CylindricalBrep::upright(
        p(origin[0],origin[1],origin[2]),2.0,5.0,GeometryTolerance::default()).unwrap()).unwrap()
}
fn at(m:&BrepModel,axis:usize,value:f64)->FaceId{
    m.faces.iter().find(|face|{
        let Some(Surface3::Plane(plane))=m.face_surface(face.id) else{return false};
        let x=[plane.frame.origin.x,plane.frame.origin.y,plane.frame.origin.z];
        let n=[plane.frame.normal.x,plane.frame.normal.y,plane.frame.normal.z];
        (x[axis]-value).abs()<1e-10 && n[axis].abs()>0.9
    }).expect("required planar face missing").id
}
fn engine()->TrimmedIntersectionEngine{TrimmedIntersectionEngine::new(GeometryTolerance::default()).unwrap()}
fn scenario(name:&str)->(BrepModel,FaceId,BrepModel,FaceId){
    let a=block([0.0,0.0,0.0],[4.0,4.0,4.0]);
    let bottom=at(&a,1,0.0);
    match name {
        "transverse_segment"=>{
            let b=block([2.0,-1.0,1.0],[4.0,4.0,4.0]);
            let side=at(&b,0,2.0);(a,bottom,b,side)
        }
        "disjoint_trims"=>{
            let b=block([7.0,-1.0,1.0],[2.0,4.0,4.0]);
            let side=at(&b,0,7.0);(a,bottom,b,side)
        }
        "corner_contact"=>{
            let b=block([4.0,-1.0,4.0],[2.0,1.0,2.0]);
            let side=at(&b,0,4.0);(a,bottom,b,side)
        }
        "boundary_line"=>{
            let b=block([4.0,-1.0,1.0],[2.0,1.0,4.0]);
            let side=at(&b,0,4.0);(a,bottom,b,side)
        }
        "coincident_support"=>{
            let b=block([1.0,0.0,1.0],[2.0,2.0,2.0]);
            let face=at(&b,1,0.0);(a,bottom,b,face)
        }
        "full_circle"=>{
            let a=block([-3.0,2.0,-3.0],[6.0,1.0,6.0]);
            let f=at(&a,1,2.0);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        "half_circle"=>{
            let a=block([0.0,2.0,-3.0],[3.0,1.0,6.0]);
            let f=at(&a,1,2.0);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        "two_generators"=>{
            let a=block([0.0,1.0,-3.0],[2.0,2.0,6.0]);
            let f=at(&a,0,0.0);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        "tangent_generator"=>{
            let a=block([2.0,1.0,-3.0],[2.0,2.0,6.0]);
            let f=at(&a,0,2.0);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        "near_tangent_uncertain"=>{
            let x=2.0+4e-10;
            let a=block([x,1.0,-3.0],[2.0,2.0,6.0]);
            let f=at(&a,0,x);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        "cap_seam_boundary"=>{
            let c=cyl([0.0,0.0,0.0]);(c.clone(),FaceId(0),c,FaceId(2))
        }
        "plane_outside_height"=>{
            let a=block([-3.0,7.0,-3.0],[6.0,1.0,6.0]);
            let f=at(&a,1,7.0);let c=cyl([0.0,0.0,0.0]);
            (a,f,c,FaceId(2))
        }
        _=>panic!("unknown pinned fixture {name}"),
    }
}
fn xyz(a:Point3)->[f64;3]{[a.x,a.y,a.z]}
fn distance(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}
fn approx(a:f64,b:f64){assert!((a-b).abs()<1e-8,"{a} vs {b}");}
fn check(f:&Fixture){
    let (a,fa,b,fb)=scenario(&f.scenario);
    let answer=engine().intersect_faces(&a,fa,&b,fb).unwrap_or_else(|e|panic!("{}: {e}",f.id));
    assert_eq!(answer.kind(),f.expected,"{}: {answer:?}",f.id);
    assert_eq!(answer.branches().len(),f.branches,"{}: {answer:?}",f.id);
    assert_eq!(answer.contacts().len(),f.points,"{}: {answer:?}",f.id);
    let before=answer.branches().to_vec();
    for branch in &before {
        assert!(branch.parameter_end>branch.parameter_start);
        approx(distance(branch.evaluate(0.0).unwrap(),branch.start),0.0);
        approx(distance(branch.evaluate(1.0).unwrap(),branch.end),0.0);
    }
    match f.scenario.as_str() {
        "transverse_segment"=>{
            let br=&before[0];
            let mut z=[br.start.z,br.end.z];z.sort_by(f64::total_cmp);
            approx(z[0],1.0);approx(z[1],4.0);
            for t in [0.0,0.25,0.5,0.75,1.0] {
                let q=br.evaluate(t).unwrap();approx(q.x,2.0);approx(q.y,0.0);
                assert!(q.z>=1.0-1e-9 && q.z<=4.0+1e-9);
            }
        }
        "corner_contact"=>{
            let pt=answer.contacts()[0].point;
            for (x,y) in xyz(pt).into_iter().zip([4.0,0.0,4.0]) {approx(x,y);}
        }
        "boundary_line"=>assert_eq!(before[0].contact,TrimContact::BoundaryOverlap),
        "full_circle"=>{
            assert!(before[0].closed);
            approx(before[0].parameter_end-before[0].parameter_start,TAU);
            for t in [0.0,0.1,0.27,0.6,1.0] {
                let q=before[0].evaluate(t).unwrap();
                approx(q.y,2.0);approx(q.x*q.x+q.z*q.z,4.0);
            }
        }
        "half_circle"=>{
            let angle: f64=before.iter().map(|b|b.parameter_end-b.parameter_start).sum();
            approx(angle,PI);
            for br in &before {
                for t in [0.0,0.3,0.5,0.9,1.0] {
                    let q=br.evaluate(t).unwrap();
                    approx(q.y,2.0);approx(q.x*q.x+q.z*q.z,4.0);
                    assert!(q.x>=-1e-8);
                }
            }
        }
        "two_generators"|"tangent_generator"=>{
            for br in &before {
                for t in [0.0,0.3,0.5,0.75,1.0] {
                    let q=br.evaluate(t).unwrap();approx(q.x,if f.scenario=="two_generators" {0.0}else{2.0});
                    approx(q.x*q.x+q.z*q.z,4.0);
                    assert!(q.y>=1.0-1e-8 && q.y<=3.0+1e-8);
                }
            }
            if f.scenario=="tangent_generator" {assert_eq!(before[0].contact,TrimContact::Tangent);}
        }
        "near_tangent_uncertain"=>assert!(matches!(answer,TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::NearTangency})),
        "cap_seam_boundary"=>{
            assert!(before[0].closed);
            assert_eq!(before[0].contact,TrimContact::BoundaryOverlap);
            for t in [0.0,0.13,0.37,0.75,1.0] {
                let q=before[0].evaluate(t).unwrap();
                approx(q.y,0.0);approx(q.x*q.x+q.z*q.z,4.0);
            }
        }
        _=>{},
    }
}
#[test]
fn pinned_brep_trim_intersections_match_independent_analytic_oracles(){
    let corpus:Corpus=serde_json::from_str(FIXTURES).unwrap();
    assert_eq!(corpus.schema,"rustsolid.trimmed_intersections.v1");
    assert!(!corpus.title.is_empty());
    assert_eq!(corpus.cases.len(),12);
    let mut ids=BTreeSet::new();
    for fixture in &corpus.cases {
        assert!(ids.insert(fixture.id.clone()),"duplicate fixture {}",fixture.id);
        check(fixture);
    }
}
#[test]
fn face_intersection_is_symmetric_as_geometric_point_sets(){
    for name in ["transverse_segment","full_circle","half_circle","two_generators","cap_seam_boundary"]{
        let (a,fa,b,fb)=scenario(name);
        let ab=engine().intersect_faces(&a,fa,&b,fb).unwrap();
        let ba=engine().intersect_faces(&b,fb,&a,fa).unwrap();
        assert_eq!(ab.kind(),ba.kind(),"{name}");
        assert_eq!(ab.branches().len(),ba.branches().len(),"{name}");
        for branch in ab.branches(){
            let reverse=ba.branches().iter().find_map(|other|{
                let direct=distance(other.start,branch.start)
                    +distance(other.end,branch.end);
                let backwards=distance(other.end,branch.start)
                    +distance(other.start,branch.end);
                if direct<1e-7 {Some((other,false))}
                else if backwards<1e-7 {Some((other,true))}
                else {None}
            }).expect("swapped operands changed analytic trimmed segment endpoints");
            for t in [0.1,0.33,0.5,0.7,0.95] {
                let p=branch.evaluate(t).unwrap();
                let q=reverse.0.evaluate(if reverse.1 {1.0-t}else{t}).unwrap();
                assert!(distance(p,q)<1e-7,"{name} point mismatch");
            }
        }
    }
}
#[test]
fn invalid_face_and_tampered_model_are_rejected_without_mutation(){
    let (a,fa,b,fb)=scenario("transverse_segment");
    let old=format!("{a:?}");
    assert!(engine().intersect_faces(&a,FaceId(999),&b,fb).is_err());
    assert_eq!(format!("{a:?}"),old);
    let mut damaged=a.clone();
    damaged.coedges[0].next=rustsolid::CoedgeId(999);
    assert!(engine().intersect_faces(&damaged,fa,&b,fb).is_err());
    let (a,fa,b,fb)=scenario("coincident_support");
    assert!(matches!(engine().intersect_faces(&a,fa,&b,fb).unwrap(),
        TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::CoincidentSupports}));
}
