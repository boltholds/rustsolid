//! A pinned, first-party oracle suite for immutable, symbolic topology plans.
//! The expected classification and point coordinates are independent cuboid /
//! circle mathematics, never another CAD backend's topology output.
use rustsolid::{
    BrepModel, CylindricalBrep, FaceId, GeometryTolerance,
    GraphArcId, GraphNodeKind, GraphOperand, GraphState,
    IntersectionGraph, IntersectionGraphBuilder, ModelUnits, Point3,
    PlanArcRole, PlanFaceEffect, PlanRequirement, PlanState,
    PlanVertexPlacement, Solid, Surface3, TopologyDeltaPlanner,
};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

const CASES:&str=include_str!("fixtures/topology_delta_plan.v1.json");
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Corpus {schema:String,title:String,cases:Vec<Fixture>}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture {scenario:String,state:String,nodes:usize,arcs:usize,
    face_changes:usize,requirements:Vec<String>}
fn point(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn block(origin:[f64;3],dims:[f64;3])->BrepModel{
    BrepModel::from_polyhedral(&Solid::block(point(origin[0],origin[1],origin[2]),
        dims[0],dims[1],dims[2]).unwrap(),ModelUnits::default()).unwrap()
}
fn cylinder(origin:[f64;3])->BrepModel{
    BrepModel::from_cylinder(&CylindricalBrep::upright(
        point(origin[0],origin[1],origin[2]),2.,5.,GeometryTolerance::default()).unwrap()).unwrap()
}
fn plane_at(model:&BrepModel,axis:usize,value:f64)->FaceId{
    model.faces.iter().find_map(|f|{
        let Surface3::Plane(plane)=model.face_surface(f.id)? else{return None};
        let a=[plane.frame.origin.x,plane.frame.origin.y,plane.frame.origin.z];
        let n=[plane.frame.normal.x,plane.frame.normal.y,plane.frame.normal.z];
        if (a[axis]-value).abs()<1e-10 && n[axis].abs()>0.9 {Some(f.id)}else{None}
    }).expect("pinned analytic plane missing")
}
fn scenario(id:&str)->(BrepModel,FaceId,BrepModel,FaceId){
    let a=block([0.,0.,0.],[4.,4.,4.]);let bottom=plane_at(&a,1,0.);
    match id {
        "transverse_segment"=>{let b=block([2.,-1.,1.],[4.,4.,4.]);let fb=plane_at(&b,0,2.);(a,bottom,b,fb)}
        "disjoint_trims"=>{let b=block([7.,-1.,1.],[2.,4.,4.]);let fb=plane_at(&b,0,7.);(a,bottom,b,fb)}
        "corner_contact"=>{let b=block([4.,-1.,4.],[2.,1.,2.]);let fb=plane_at(&b,0,4.);(a,bottom,b,fb)}
        "boundary_line"=>{let b=block([4.,-1.,1.],[2.,1.,4.]);let fb=plane_at(&b,0,4.);(a,bottom,b,fb)}
        "coincident_support"=>{let b=block([1.,0.,1.],[2.,2.,2.]);let fb=plane_at(&b,1,0.);(a,bottom,b,fb)}
        "full_circle"=>{let a=block([-3.,2.,-3.],[6.,1.,6.]);let fa=plane_at(&a,1,2.);
            (a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        "half_circle"=>{let a=block([0.,2.,-3.],[3.,1.,6.]);let fa=plane_at(&a,1,2.);
            (a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        "two_generators"=>{let a=block([0.,1.,-3.],[2.,2.,6.]);let fa=plane_at(&a,0,0.);
            (a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        "tangent_generator"=>{let a=block([2.,1.,-3.],[2.,2.,6.]);let fa=plane_at(&a,0,2.);
            (a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        "near_tangent_uncertain"=>{let x=2.+4e-10;let a=block([x,1.,-3.],[2.,2.,6.]);
            let fa=plane_at(&a,0,x);(a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        "cap_seam_boundary"=>{let c=cylinder([0.,0.,0.]);(c.clone(),FaceId(0),c,FaceId(2))}
        "plane_outside_height"=>{let a=block([-3.,7.,-3.],[6.,1.,6.]);let fa=plane_at(&a,1,7.);
            (a,fa,cylinder([0.,0.,0.]),FaceId(2))}
        _=>panic!("unregistered pinned scenario: {id}"),
    }
}
fn make_graph(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId)->IntersectionGraph{
    IntersectionGraphBuilder::new(GeometryTolerance::default()).unwrap()
        .intersect_faces(a,fa,b,fb).unwrap()
}
fn name(r:PlanRequirement)->&'static str {
    match r {
        PlanRequirement::NumericalRefinement=>"numerical_refinement",
        PlanRequirement::GeneralizedBrepEditor=>"generalized_brep_editor",
        PlanRequirement::RegionClassification=>"region_classification",
        PlanRequirement::InnerLoopConstruction=>"inner_loop_construction",
        PlanRequirement::OpenTraceCompletion=>"open_trace_completion",
        PlanRequirement::CurvedEdgeSubdivision=>"curved_edge_subdivision",
        PlanRequirement::MultipleSplitsSameEdge=>"multiple_splits_same_edge",
        PlanRequirement::TangencyClassification=>"tangency_classification",
        PlanRequirement::BoundaryReconciliation=>"boundary_reconciliation",
    }
}
fn status(r:PlanState)->&'static str {
    match r{PlanState::NoChange=>"no_change",PlanState::PreviewOnly=>"preview_only",PlanState::Held=>"held"}
}
#[test]
fn twelve_pinned_topology_intent_cases_pass_independent_expectations(){
    let corpus:Corpus=serde_json::from_str(CASES).unwrap();
    assert_eq!(corpus.schema,"rustsolid.topology_delta_plan.v1");
    assert!(!corpus.title.is_empty());
    assert_eq!(corpus.cases.len(),12);
    let mut seen=BTreeSet::new();
    for c in &corpus.cases {
        assert!(seen.insert(&c.scenario),"duplicate case {}",c.scenario);
        let (a,fa,b,fb)=scenario(&c.scenario);
        let graph=make_graph(&a,fa,&b,fb);
        let plan=TopologyDeltaPlanner::new().plan(&graph,&a,fa,&b,fb).unwrap();
        assert_eq!(status(plan.state),c.state,"{}: {plan:?}",c.scenario);
        assert_eq!(plan.vertices.len(),c.nodes*2,"{}",c.scenario);
        assert_eq!(plan.arcs.len(),c.arcs,"{}",c.scenario);
        assert_eq!(plan.face_changes.len(),c.face_changes,"{}",c.scenario);
        let actual=plan.requirements.iter().copied().map(name).collect::<BTreeSet<_>>();
        assert_eq!(actual,c.requirements.iter().map(String::as_str).collect(),"{}",c.scenario);
        assert!(!plan.can_apply,"read-only plan must never be ready to commit");
        plan.validate_against(&graph,&a,fa,&b,fb).unwrap();
        for arc in &plan.arcs {
            assert_eq!(arc.source_faces[0].operand,GraphOperand::First);
            assert_eq!(arc.source_faces[1].operand,GraphOperand::Second);
            assert!(arc.parameter_end>arc.parameter_start);
        }
        if c.scenario=="transverse_segment" {
            assert_eq!(plan.split_edge_requests(),2);
            assert!(plan.face_changes.iter().all(|f|f.effect==PlanFaceEffect::PendingOpenTrace));
            for e in &plan.arcs {
                assert_eq!(e.role,PlanArcRole::ProposeNewCarrier);
                let m=e.analytic_curve.evaluate(0.5*(e.parameter_start+e.parameter_end)).unwrap();
                assert!((m.x-2.).abs()<1e-8 && m.y.abs()<1e-8 && (m.z-2.5).abs()<1e-8);
            }
        }
        if c.scenario=="full_circle" {
            assert_eq!(plan.vertices.len(),2);
            assert!(plan.vertices.iter().all(|v|matches!(v.placement,PlanVertexPlacement::VirtualPeriodicAnchor)));
            assert!(plan.face_changes.iter().all(|f|f.effect==PlanFaceEffect::CandidateInnerLoop));
            assert!(plan.arcs[0].closed);
            assert!((plan.arcs[0].parameter_end-plan.arcs[0].parameter_start-TAU).abs()<1e-8);
        }
        if c.scenario=="half_circle" {
            assert!(plan.face_changes.iter().any(|f|f.effect==PlanFaceEffect::CandidateBoundarySplit));
            assert!(plan.face_changes.iter().any(|f|f.effect==PlanFaceEffect::PendingOpenTrace));
            assert!((plan.arcs[0].parameter_end-plan.arcs[0].parameter_start-PI).abs()<1e-8);
        }
        if c.scenario=="corner_contact" {
            assert!(plan.vertices.iter().all(|v|matches!(v.placement,PlanVertexPlacement::ContactOnly)));
        }
        if c.scenario=="boundary_line" {
            assert_eq!(plan.arcs[0].role,PlanArcRole::ReconcileBoundary);
        }
    }
}
#[test]
fn plan_is_deterministic_immutable_and_catches_rewritten_vertex_or_provenance(){
    let (a,fa,b,fb)=scenario("half_circle");
    let a_before=format!("{a:?}");let b_before=format!("{b:?}");
    let graph=make_graph(&a,fa,&b,fb);
    let p=TopologyDeltaPlanner::new().plan(&graph,&a,fa,&b,fb).unwrap();
    let twice=TopologyDeltaPlanner::new().plan(&graph,&a,fa,&b,fb).unwrap();
    assert_eq!(p,twice);
    assert_eq!(format!("{a:?}"),a_before);
    assert_eq!(format!("{b:?}"),b_before);
    let mut tampered=p.clone();tampered.vertices[0].point.x+=0.25;
    assert!(tampered.validate_against(&graph,&a,fa,&b,fb).is_err());
    let mut tampered=p.clone();tampered.arcs[0].source_faces[1].face=FaceId(777);
    assert!(tampered.validate_against(&graph,&a,fa,&b,fb).is_err());
    let mut tampered=p.clone();tampered.requirements.clear();
    assert!(tampered.validate_against(&graph,&a,fa,&b,fb).is_err());
    let mut tampered=p.clone();tampered.can_apply=true;
    assert!(tampered.validate_against(&graph,&a,fa,&b,fb).is_err());
    let mut stale=b.clone();
    let idx=stale.faces[fb.0 as usize].geometry.surface.0 as usize;
    if let Surface3::Cylinder(mut c)=stale.geometry.surfaces[idx]{
        c.frame.origin.y+=0.1;stale.geometry.surfaces[idx]=Surface3::Cylinder(c);
    }
    assert!(p.validate_against(&graph,&a,fa,&stale,fb).is_err());
}
#[test]
fn swapped_operands_do_not_alias_local_vertex_edge_or_face_identifiers(){
    let (a,fa,b,fb)=scenario("transverse_segment");
    let ab_graph=make_graph(&a,fa,&b,fb);
    let ba_graph=make_graph(&b,fb,&a,fa);
    let ab=TopologyDeltaPlanner::new().plan(&ab_graph,&a,fa,&b,fb).unwrap();
    let ba=TopologyDeltaPlanner::new().plan(&ba_graph,&b,fb,&a,fa).unwrap();
    assert_eq!(ab.state,ba.state);
    assert_eq!(ab.arcs.len(),ba.arcs.len());
    assert_eq!(ab.vertices.len(),ba.vertices.len());
    assert_eq!(ab.split_edge_requests(),ba.split_edge_requests());
    for operand in [GraphOperand::First,GraphOperand::Second]{
        assert_eq!(ab.vertices.iter().filter(|v|v.operand==operand).count(),ab_graph.nodes.len());
        assert_eq!(ba.vertices.iter().filter(|v|v.operand==operand).count(),ba_graph.nodes.len());
    }
}
#[test]
fn graph_ambiguous_tangency_and_overlap_do_not_propose_unconditional_boolean_cut(){
    for scenario_id in ["coincident_support","near_tangent_uncertain","tangent_generator","boundary_line","cap_seam_boundary"]{
        let (a,fa,b,fb)=scenario(scenario_id);
        let g=make_graph(&a,fa,&b,fb);
        let plan=TopologyDeltaPlanner::new().plan(&g,&a,fa,&b,fb).unwrap();
        assert_eq!(plan.state,PlanState::Held,"{scenario_id}");
        assert!(!plan.is_mutation_ready());
        assert!(plan.arcs.iter().all(|a|a.role!=PlanArcRole::ProposeNewCarrier));
        if matches!(g.state,GraphState::Ambiguous(_)){assert!(plan.arcs.is_empty());}
    }
}
#[test]
fn angular_periodic_anchor_is_not_scheduled_as_physical_edge_subdivision(){
    for case in ["full_circle","cap_seam_boundary"]{
        let (a,fa,b,fb)=scenario(case);
        let g=make_graph(&a,fa,&b,fb);
        assert_eq!(g.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
        let plan=TopologyDeltaPlanner::new().plan(&g,&a,fa,&b,fb).unwrap();
        assert!(plan.vertices.iter().all(|v|matches!(v.placement,PlanVertexPlacement::VirtualPeriodicAnchor)));
        assert_eq!(plan.split_edge_requests(),0);
    }
}
