//! Pinned analytic witnesses and metamorphic checks for the intersection engine.
//! No OCCT/Parasolid outputs are used as acceptance criteria.
use rustsolid::{
    CylinderSurface, Frame3, GeometryTolerance, IntersectionAmbiguity,
    IntersectionContact, IntersectionCurve, PlaneSurface, Point3,
    IntersectionPair as Operator, ALL_INTERSECTION_PAIRS,
    Surface3, SurfaceIntersection, SurfaceIntersectionEngine,
};
use serde::Deserialize;
use std::collections::{BTreeMap,BTreeSet};
use std::f64::consts::TAU;
const RAW:&str=include_str!("fixtures/intersections.v1.json");
const SCHEMA:&str="rustsolid.intersections.v1";

#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Deserialize)]
#[serde(rename_all="snake_case")]
enum Category { Analytic, Metamorphic, Boundary }
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct FrameInput { origin:[f64;3],x:[f64;3],y:[f64;3] }
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct CylinderInput { origin:[f64;3],radius:f64,x:[f64;3],y:[f64;3] }
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture {
    id:String,operator:Operator,category:Category,
    plane_a:Option<FrameInput>,plane_b:Option<FrameInput>,
    plane:Option<FrameInput>,cylinder:Option<CylinderInput>,
    expected:String,branches:usize,
}
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct Corpus {schema:String,title:String,fixtures:Vec<Fixture>}
fn point(a:[f64;3])->Point3{Point3{x:a[0],y:a[1],z:a[2]}}
fn xyz(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn len(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn sub(a:Point3,b:Point3)->Point3{xyz(a.x-b.x,a.y-b.y,a.z-b.z)}
fn add(a:Point3,b:Point3)->Point3{xyz(a.x+b.x,a.y+b.y,a.z+b.z)}
fn mul(a:Point3,f:f64)->Point3{xyz(a.x*f,a.y*f,a.z*f)}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn frame(a:&FrameInput)->Frame3{Frame3::from_axes(point(a.origin),point(a.x),point(a.y),tol()).unwrap()}
fn surfaces(f:&Fixture)->(Surface3,Surface3){
    match f.operator {
        Operator::PlanePlane=>(Surface3::Plane(PlaneSurface{frame:frame(f.plane_a.as_ref().unwrap())}),
            Surface3::Plane(PlaneSurface{frame:frame(f.plane_b.as_ref().unwrap())})),
        Operator::PlaneCylinder=>(Surface3::Plane(PlaneSurface{frame:frame(f.plane.as_ref().unwrap())}),
            Surface3::Cylinder({let c=f.cylinder.as_ref().unwrap();
                CylinderSurface::new(Frame3::from_axes(point(c.origin),point(c.x),point(c.y),tol()).unwrap(),c.radius,tol()).unwrap()})),
    }
}
fn manifest(raw:&str)->Result<Corpus,String>{
    let out:Corpus=serde_json::from_str(raw).map_err(|e|e.to_string())?;
    if out.schema!=SCHEMA || out.title.trim().is_empty(){return Err("wrong schema/title".into());}
    let mut kinds=BTreeMap::<Operator,BTreeMap<Category,usize>>::new();
    let mut ids=BTreeSet::new();
    for f in &out.fixtures {
        let root=match f.operator {Operator::PlanePlane=>"plane_plane/",Operator::PlaneCylinder=>"plane_cylinder/"};
        if !f.id.starts_with(root) || !ids.insert(f.id.clone()){
            return Err(format!("duplicate or invalid case: {}",f.id));
        }
        if !["empty","coincident","ambiguous","crossing","tangent"].contains(&f.expected.as_str()){
            return Err(format!("unknown result {}",f.expected));
        }
        if f.branches>2 || (["empty","coincident","ambiguous"].contains(&f.expected.as_str()) && f.branches!=0)
            || (["crossing","tangent"].contains(&f.expected.as_str()) && f.branches==0){
            return Err("invalid expected branch count".into());
        }
        match f.operator {
            Operator::PlanePlane if f.plane_a.is_none() || f.plane_b.is_none() || f.plane.is_some() || f.cylinder.is_some()=>
                return Err("plane-plane requires two planar operands".into()),
            Operator::PlaneCylinder if f.plane.is_none() || f.cylinder.is_none() || f.plane_a.is_some() || f.plane_b.is_some()=>
                return Err("plane-cylinder requires both typed operands".into()),
            _=>{},
        }
        *kinds.entry(f.operator).or_default().entry(f.category).or_default()+=1;
    }
    let declared=ALL_INTERSECTION_PAIRS.into_iter().collect::<BTreeSet<_>>();
    if kinds.keys().copied().collect::<BTreeSet<_>>()!=declared || out.fixtures.len()<10 {
        return Err("intersection operator registry differs from pinned scenarios".into());
    }
    for (_,categories) in kinds {
        for (kind,minimum) in [(Category::Analytic,1),(Category::Metamorphic,1),(Category::Boundary,2)]{
            if categories.get(&kind).copied().unwrap_or(0)<minimum {
                return Err(format!("missing {kind:?} witness"));
            }
        }
    }
    Ok(out)
}
fn close(a:f64,b:f64,eps:f64){assert!((a-b).abs()<=eps,"{a} vs {b}, tolerance={eps}");}
fn space_close(a:Point3,b:Point3,eps:f64){close(len(sub(a,b)),0.,eps);}
fn line_distance(p:Point3,c:IntersectionCurve)->f64{
    let IntersectionCurve::Line(l)=c else{panic!("not a line")};
    let v=sub(p,l.origin);let t=dot(v,l.direction);
    len(sub(v,mul(l.direction,t)))
}
fn analytic_oracle(id:&str,result:&SurfaceIntersection){
    match id {
        "plane_plane/orthogonal"=>{
            for t in [-4.,0.,3.]{let p=result.curves()[0].evaluate(t).unwrap();close(p.x,0.,1e-11);close(p.z,0.,1e-11);}
        }
        "plane_cylinder/perpendicular_circle"=>{
            for t in [0.1,0.9,2.4,5.8]{
                let p=result.curves()[0].evaluate(t).unwrap();
                close(p.y,3.,1e-10);close(p.x*p.x+p.z*p.z,4.,1e-10);
            }
        }
        "plane_cylinder/oblique_section"=>{
            for t in [0.1,0.9,2.4,5.8]{
                let p=result.curves()[0].evaluate(t).unwrap();
                close(p.x+2.*(p.y-1.),0.,1e-10);close(p.x*p.x+p.z*p.z,4.,1e-10);
            }
        }
        "plane_cylinder/two_lines"=>{
            let p=result.curves().iter().map(|c|c.evaluate(0.).unwrap()).collect::<Vec<_>>();
            close(p[0].x,0.,1e-10);close(p[1].x,0.,1e-10);
            close(p[0].z.abs(),2.,1e-10);close(p[1].z.abs(),2.,1e-10);
            assert!(p[0].z*p[1].z<0.);
        }
        _=>{},
    }
}
fn translated(s:Surface3, delta:Point3)->Surface3{
    match s {
        Surface3::Plane(mut p)=>{p.frame.origin=add(p.frame.origin,delta);Surface3::Plane(p)},
        Surface3::Cylinder(mut c)=>{c.frame.origin=add(c.frame.origin,delta);Surface3::Cylinder(c)},
        Surface3::Sphere(_)=>unreachable!(),
    }
}
fn run_case(f:&Fixture){
    let (a,b)=surfaces(f);let engine=SurfaceIntersectionEngine::new(tol()).unwrap();
    let result=engine.intersect(a,b).unwrap();
    assert_eq!(result.kind(),f.expected,"{}: {result:?}",f.id);
    assert_eq!(result.curves().len(),f.branches,"{}: {result:?}",f.id);
    if !result.curves().is_empty(){
        let cert=result.certify(a,b,tol()).unwrap();
        assert!(cert.samples>=5,"{}",f.id);
        assert!(cert.max_residual_a<=cert.budget && cert.max_residual_b<=cert.budget);
    }
    if f.category==Category::Analytic {analytic_oracle(&f.id,&result);}
    if f.category==Category::Metamorphic {
        let swapped=engine.intersect(b,a).unwrap();
        assert_eq!(result.kind(),swapped.kind());
        assert_eq!(result.curves().len(),swapped.curves().len());
        if !swapped.curves().is_empty(){swapped.certify(b,a,tol()).unwrap();}
        let shift=xyz(97.,-43.,152.);
        let moved_a=translated(a,shift);let moved_b=translated(b,shift);
        let moved=engine.intersect(moved_a,moved_b).unwrap();
        assert_eq!(result.kind(),moved.kind());
        assert_eq!(result.curves().len(),moved.curves().len());
        if !moved.curves().is_empty(){moved.certify(moved_a,moved_b,tol()).unwrap();}
        if f.operator==Operator::PlanePlane {
            let point=result.curves()[0].evaluate(0.).unwrap();
            assert!(line_distance(add(point,shift),moved.curves()[0])<1e-9);
            assert!(line_distance(point,swapped.curves()[0])<1e-9);
        } else {
            // Rotate the cylinder's parameter frame without changing any
            // physical point. The same curve must be produced up to u shift.
            let Surface3::Cylinder(c)=b else{unreachable!()};
            let theta=0.39_f64;
            let nx=add(mul(c.frame.x,theta.cos()),mul(c.frame.y,theta.sin()));
            let ny=sub(mul(c.frame.y,theta.cos()),mul(c.frame.x,theta.sin()));
            let rotated=Surface3::Cylinder(CylinderSurface::new(
                Frame3::from_axes(c.frame.origin,nx,ny,tol()).unwrap(),c.radius,tol()).unwrap());
            let repar=engine.intersect(a,rotated).unwrap();
            assert_eq!(repar.kind(),result.kind());
            repar.certify(a,rotated,tol()).unwrap();
            for angle in [0.1,0.9,3.3,5.5]{
                space_close(result.curves()[0].evaluate(angle).unwrap(),
                            repar.curves()[0].evaluate(angle-theta).unwrap(),1e-9);
                space_close(add(result.curves()[0].evaluate(angle).unwrap(),shift),
                            moved.curves()[0].evaluate(angle).unwrap(),1e-9);
            }
        }
    }
    if f.id=="plane_cylinder/exact_tangent" {
        assert!(matches!(result,SurfaceIntersection::Curves{contact:IntersectionContact::Tangent,..}));
    }
    if f.id.contains("near_tangent"){
        assert!(matches!(result,SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearTangency,..}));
    }
    if f.id.contains("near_parallel"){
        assert!(matches!(result,SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearParallel,..}));
    }
}
#[test]
fn each_pinned_intersection_case_is_checked_with_analytic_or_metamorphic_math(){
    let corpus=manifest(RAW).unwrap();
    assert_eq!(corpus.fixtures.len(),14);
    for fixture in &corpus.fixtures {run_case(fixture);}
}
#[test]
fn contract_gate_rejects_missing_coverage_duplicate_ids_and_corrupted_witnesses(){
    let mut doc:serde_json::Value=serde_json::from_str(RAW).unwrap();
    doc["fixtures"].as_array_mut().unwrap().retain(|c|c["category"]!="boundary");
    assert!(manifest(&doc.to_string()).is_err());
    let mut doc:serde_json::Value=serde_json::from_str(RAW).unwrap();
    let repeated=doc["fixtures"][0]["id"].clone();
    doc["fixtures"][1]["id"]=repeated;
    assert!(manifest(&doc.to_string()).is_err());
    let mut doc:serde_json::Value=serde_json::from_str(RAW).unwrap();
    doc["fixtures"][0]["expected"]=serde_json::json!("tangent");
    let altered=manifest(&doc.to_string()).unwrap();
    let failed=std::panic::catch_unwind(||run_case(&altered.fixtures[0]));
    assert!(failed.is_err(),"forged oracle unexpectedly passed");
}
#[test]
fn unsupported_pairs_are_explicit_and_no_zero_curve_is_certified(){
    let corpus=manifest(RAW).unwrap();
    let (a,b)=surfaces(&corpus.fixtures[2]);
    let r=SurfaceIntersectionEngine::new(tol()).unwrap().intersect(a,b).unwrap();
    assert!(r.certify(a,b,tol()).is_err());
    let sphere=Surface3::Sphere(rustsolid::SphereSurface{
        frame:Frame3::from_axes(xyz(0.,0.,0.),xyz(1.,0.,0.),xyz(0.,1.,0.),tol()).unwrap(),radius:1.,
    });
    assert!(SurfaceIntersectionEngine::new(tol()).unwrap().intersect(sphere,a).is_err());
}
