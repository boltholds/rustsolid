//! Persistent, independent first-party IntersectionGraph contracts.
//! All numerical expectations follow elementary cuboid/cylinder mathematics;
//! no OCCT/Parasolid result is used as ground truth.
use rustsolid::{
    BrepModel,CylindricalBrep,FaceId,GeometryTolerance,GraphNodeKind,
    GraphOperand,GraphState,IntersectionGraph,IntersectionGraphBuilder,
    ModelUnits,Point3,Solid,Surface3,TrimContact,TrimPosition,TrimAmbiguity,
};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::f64::consts::{PI,TAU};

const CORPUS:&str=include_str!("fixtures/intersection_graph.v1.json");
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Corpus {schema:String,title:String,fixtures:Vec<Fixture>}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture {id:String,scenario:String,state:String,nodes:usize,arcs:usize}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn approx(a:f64,b:f64){assert!((a-b).abs()<1e-8,"actual {a} vs independent oracle {b}");}
fn block(origin:[f64;3],dimensions:[f64;3])->BrepModel{
    BrepModel::from_polyhedral(&Solid::block(p(origin[0],origin[1],origin[2]),
        dimensions[0],dimensions[1],dimensions[2]).unwrap(),ModelUnits::default()).unwrap()
}
fn cylinder(origin:[f64;3])->BrepModel{
    BrepModel::from_cylinder(&CylindricalBrep::upright(
        p(origin[0],origin[1],origin[2]),2.0,5.0,GeometryTolerance::default()).unwrap()).unwrap()
}
fn plane_at(model:&BrepModel,axis:usize,value:f64)->FaceId{
    model.faces.iter().find_map(|face|{
        let Surface3::Plane(plane)=model.face_surface(face.id)? else{return None};
        let x=[plane.frame.origin.x,plane.frame.origin.y,plane.frame.origin.z];
        let n=[plane.frame.normal.x,plane.frame.normal.y,plane.frame.normal.z];
        if (x[axis]-value).abs()<1e-10&&n[axis].abs()>0.9 {Some(face.id)}else{None}
    }).expect("missing requested planar face")
}
fn scenario(id:&str)->(BrepModel,FaceId,BrepModel,FaceId){
    let a=block([0.,0.,0.],[4.,4.,4.]);let bottom=plane_at(&a,1,0.);
    match id {
        "transverse_segment"=>{let b=block([2.,-1.,1.],[4.,4.,4.]);let fb=plane_at(&b,0,2.);(a,bottom,b,fb)}
        "disjoint_trims"=>{let b=block([7.,-1.,1.],[2.,4.,4.]);let fb=plane_at(&b,0,7.);(a,bottom,b,fb)}
        "corner_contact"=>{let b=block([4.,-1.,4.],[2.,1.,2.]);let fb=plane_at(&b,0,4.);(a,bottom,b,fb)}
        "boundary_line"=>{let b=block([4.,-1.,1.],[2.,1.,4.]);let fb=plane_at(&b,0,4.);(a,bottom,b,fb)}
        "coincident_support"=>{let b=block([1.,0.,1.],[2.,2.,2.]);let fb=plane_at(&b,1,0.);(a,bottom,b,fb)}
        "full_circle"=>{let plane=block([-3.,2.,-3.],[6.,1.,6.]);let face=plane_at(&plane,1,2.);
            (plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        "half_circle"=>{let plane=block([0.,2.,-3.],[3.,1.,6.]);let face=plane_at(&plane,1,2.);
            (plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        "two_generators"=>{let plane=block([0.,1.,-3.],[2.,2.,6.]);let face=plane_at(&plane,0,0.);
            (plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        "tangent_generator"=>{let plane=block([2.,1.,-3.],[2.,2.,6.]);let face=plane_at(&plane,0,2.);
            (plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        "near_tangent_uncertain"=>{let x=2.+4e-10;let plane=block([x,1.,-3.],[2.,2.,6.]);
            let face=plane_at(&plane,0,x);(plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        "cap_seam_boundary"=>{let c=cylinder([0.,0.,0.]);(c.clone(),FaceId(0),c,FaceId(2))}
        "plane_outside_height"=>{let plane=block([-3.,7.,-3.],[6.,1.,6.]);let face=plane_at(&plane,1,7.);
            (plane,face,cylinder([0.,0.,0.]),FaceId(2))}
        _=>panic!("unrecognized pinned geometry {id}"),
    }
}
fn graph(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId)->IntersectionGraph {
    IntersectionGraphBuilder::new(GeometryTolerance::default()).unwrap()
        .intersect_faces(a,fa,b,fb).unwrap()
}
fn check_fixture(case:&Fixture){
    let (a,fa,b,fb)=scenario(&case.scenario);
    let g=graph(&a,fa,&b,fb);
    assert_eq!(g.kind(),case.state,"{}: {g:?}",case.id);
    assert_eq!(g.nodes.len(),case.nodes,"{}: {g:?}",case.id);
    assert_eq!(g.arcs.len(),case.arcs,"{}: {g:?}",case.id);
    g.validate_against(&a,fa,&b,fb).unwrap();
    match case.scenario.as_str(){
        "transverse_segment"=>{
            let edge=g.arcs[0];let p0=g.nodes[edge.from.0 as usize].point;
            let p1=g.nodes[edge.to.0 as usize].point;
            let mut zs=[p0.z,p1.z];zs.sort_by(f64::total_cmp);
            approx(zs[0],1.);approx(zs[1],4.);
            for v in [p0,p1,edge.evaluate(0.5).unwrap()] {approx(v.x,2.);approx(v.y,0.);}
            assert_eq!(edge.contact,TrimContact::Crossing);
            assert_eq!(g.source_faces[0].face,fa);
            assert_eq!(g.source_faces[1].face,fb);
            assert_eq!(g.source_faces[0].operand,GraphOperand::First);
            assert_eq!(g.source_faces[1].operand,GraphOperand::Second);
        }
        "corner_contact"=>{
            assert_eq!(g.nodes[0].kind,GraphNodeKind::IsolatedContact);
            assert_eq!(g.degree(g.nodes[0].id),0);
            let q=g.nodes[0].point;approx(q.x,4.);approx(q.y,0.);approx(q.z,4.);
            assert_eq!(g.nodes[0].membership,[TrimPosition::Boundary,TrimPosition::Boundary]);
        }
        "boundary_line"=>assert_eq!(g.arcs[0].contact,TrimContact::BoundaryOverlap),
        "full_circle"|"cap_seam_boundary"=>{
            assert_eq!(g.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
            assert!(g.arcs[0].closed);
            assert_eq!(g.arcs[0].from,g.arcs[0].to);
            assert_eq!(g.degree(g.nodes[0].id),2);
            assert_eq!(g.incident_arcs(g.nodes[0].id),vec![g.arcs[0].id]);
            approx(g.arcs[0].parameter_end-g.arcs[0].parameter_start,TAU);
            assert!(g.nodes[0].boundary_sources.iter().all(|w|w.edge.0!=2));
            if case.scenario=="cap_seam_boundary" {
                assert_eq!(g.arcs[0].contact,TrimContact::BoundaryOverlap);
                assert!(g.nodes[0].boundary_sources.iter().any(|w|w.operand==GraphOperand::Second));
            }
        }
        "half_circle"=>{
            // The seam is a chart boundary, not a new physical vertex.
            assert_eq!(g.nodes.len(),2);
            assert_eq!(g.arcs.len(),1);
            approx(g.arcs[0].parameter_end-g.arcs[0].parameter_start,PI);
            for t in [0.,0.2,0.5,0.8,1.] {
                let q=g.arcs[0].evaluate(t).unwrap();
                assert!(q.x>=-1e-9);
                approx(q.y,2.);approx(q.x*q.x+q.z*q.z,4.);
            }
        }
        "two_generators"|"tangent_generator"=>{
            for edge in &g.arcs {
                for t in [0.,0.25,0.5,0.75,1.] {
                    let q=edge.evaluate(t).unwrap();
                    let x=if case.scenario=="tangent_generator"{2.}else{0.};
                    approx(q.x,x);approx(q.x*q.x+q.z*q.z,4.);
                    assert!(q.y>=1.-1e-9 && q.y<=3.+1e-9);
                }
                if case.scenario=="tangent_generator"{assert_eq!(edge.contact,TrimContact::Tangent);}
            }
        }
        "near_tangent_uncertain"=>assert_eq!(g.state,GraphState::Ambiguous(TrimAmbiguity::NearTangency)),
        _=>{},
    }
    for n in &g.nodes {
        if n.membership.contains(&TrimPosition::Boundary) {
            assert!(!n.boundary_sources.is_empty(),"{} has an unexplained B-rep boundary point",case.id);
        }
    }
}
#[test]
fn pinned_graph_provenance_and_first_party_geometric_oracles(){
    let manifest:Corpus=serde_json::from_str(CORPUS).unwrap();
    assert_eq!(manifest.schema,"rustsolid.intersection_graph.v1");
    assert!(!manifest.title.trim().is_empty());
    assert_eq!(manifest.fixtures.len(),12);
    let mut ids=BTreeSet::new();
    for fixture in &manifest.fixtures {
        assert!(ids.insert(fixture.id.clone()));
        check_fixture(fixture);
    }
}
#[test]
fn graph_source_provenance_and_carriers_detect_tampering(){
    let (a,fa,b,fb)=scenario("transverse_segment");
    let base=graph(&a,fa,&b,fb);
    let mut broken=base.clone();broken.arcs[0].source_faces[1].face=FaceId(999);
    assert!(broken.validate_against(&a,fa,&b,fb).is_err());
    let mut broken=base.clone();broken.arcs[0].parameter_end+=1.5;
    assert!(broken.validate_against(&a,fa,&b,fb).is_err());
    let mut broken=base.clone();broken.nodes[0].point.x+=0.1;
    assert!(broken.validate_against(&a,fa,&b,fb).is_err());
    let mut broken=base.clone();broken.nodes[0].boundary_sources.clear();
    assert!(broken.validate_against(&a,fa,&b,fb).is_err());
    let mut broken=base.clone();broken.nodes[0].membership=[TrimPosition::Interior,TrimPosition::Interior];
    assert!(broken.validate_against(&a,fa,&b,fb).is_err());
    let mut moved=b.clone();
    let old=moved.geometry.surfaces[0];
    if let Surface3::Plane(mut pl)=old {pl.frame.origin.x+=0.1;moved.geometry.surfaces[0]=Surface3::Plane(pl);}
    assert!(base.validate_against(&a,fa,&moved,fb).is_err());
}
#[test]
fn source_order_changes_labels_but_not_geometric_set(){
    for id in ["transverse_segment","half_circle","full_circle","two_generators","cap_seam_boundary"]{
        let (a,fa,b,fb)=scenario(id);
        let ab=graph(&a,fa,&b,fb);
        let ba=graph(&b,fb,&a,fa);
        assert_eq!((ab.nodes.len(),ab.arcs.len()),(ba.nodes.len(),ba.arcs.len()),"{id}");
        for node in &ab.nodes {
            assert!(ba.nodes.iter().any(|n|distance(n.point,node.point)<1e-8),"{id}: source swap changed point set");
        }
        for edge in &ab.arcs {
            let equivalent=ba.arcs.iter().find_map(|other|{
                let a=ab.nodes[edge.from.0 as usize].point;
                let b=ab.nodes[edge.to.0 as usize].point;
                let c=ba.nodes[other.from.0 as usize].point;
                let d=ba.nodes[other.to.0 as usize].point;
                if distance(a,c)+distance(b,d)<1e-7 {Some((other,false))}
                else if distance(a,d)+distance(b,c)<1e-7 {Some((other,true))}
                else {None}
            }).expect("operand swap changed graph endpoint identities");
            for t in [0.1,0.3,0.5,0.7,0.9] {
                let q=edge.evaluate(t).unwrap();
                let r=equivalent.0.evaluate(if equivalent.1{1.-t}else{t}).unwrap();
                assert!(distance(q,r)<1e-7,"{id}: swapped operands changed analytic curve");
            }
        }
    }
}
#[test]
fn translated_graph_preserves_curves_and_local_precision(){
    let (a,fa,b,fb)=scenario("full_circle");
    let before=graph(&a,fa,&b,fb);
    let shift=[91.,-43.,152.];
    let a2=block([-3.+shift[0],2.+shift[1],-3.+shift[2]],[6.,1.,6.]);
    let f2=plane_at(&a2,1,2.+shift[1]);
    let b2=cylinder(shift);
    let after=graph(&a2,f2,&b2,FaceId(2));
    assert_eq!(before.state,after.state);
    assert_eq!((before.nodes.len(),before.arcs.len()),(after.nodes.len(),after.arcs.len()));
    approx(before.tolerance.length,after.tolerance.length);
    for t in [0.,0.11,0.3,0.67,1.] {
        let p1=before.arcs[0].evaluate(t).unwrap();
        let p2=after.arcs[0].evaluate(t).unwrap();
        approx(p2.x-p1.x,shift[0]);approx(p2.y-p1.y,shift[1]);approx(p2.z-p1.z,shift[2]);
    }
}
#[test]
fn ambiguous_cases_never_expose_fictitious_graph_edges(){
    for scenario_id in ["near_tangent_uncertain","coincident_support"]{
        let (a,fa,b,fb)=scenario(scenario_id);
        let g=graph(&a,fa,&b,fb);
        assert!(matches!(g.state,GraphState::Ambiguous(_)));
        assert!(g.nodes.is_empty()&&g.arcs.is_empty());
        assert!(g.support_certificate.is_none());
    }
}

#[test]
fn uv_reparameterization_preserves_periodic_graph_as_one_closed_arc(){
    let (plane,face,c,fc)=scenario("full_circle");
    let original=graph(&plane,face,&c,fc);
    let tol=GeometryTolerance::default();
    let theta=0.37_f64;
    let frame=rustsolid::Frame3::from_axes(p(0.,0.,0.),
        p(theta.cos(),0.,-theta.sin()),p(-theta.sin(),0.,-theta.cos()),tol).unwrap();
    let rotated=BrepModel::from_cylinder(&CylindricalBrep::new(frame,2.,5.,tol).unwrap()).unwrap();
    let changed=graph(&plane,face,&rotated,FaceId(2));
    assert_eq!(changed.nodes.len(),1);
    assert_eq!(changed.arcs.len(),1);
    assert_eq!(changed.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
    assert!(changed.arcs[0].closed);
    for theta0 in [0.0_f64,0.2,1.1,3.5,5.0]{
        let a=original.arcs[0].curve.evaluate(theta0+theta).unwrap();
        let b=changed.arcs[0].curve.evaluate(theta0).unwrap();
        assert!(distance(a,b)<1e-8);
    }
}
