//! Independent oblique plane / finite tilted cylinder preview contracts.
//! Mathematical witnesses are first-party; no external CAD engine is used.
use rustsolid::{
    BooleanOperation,BooleanRegionClassifier,BrepBody,BrepModel,
    CylindricalBrep,FaceId,Frame3,GeometryTolerance,GraphNodeKind,
    GraphOperand,GraphRegionState,GraphState,IntersectionGraphBuilder,
    PlanRequirement,PlanState,Point2,Point3,Solid,Surface3,TopologyDeltaPlanner,
    TrimAmbiguity,TrimContact,TrimPosition,TrimmedIntersectionEngine,
};
use std::f64::consts::TAU;

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64 {(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn approx(a:f64,b:f64) {assert!((a-b).abs()<2e-8,"{a} != {b}");}
fn tol()->GeometryTolerance {GeometryTolerance::default()}

/// Horizontal finite *face* (at y=2) intersecting an oblique analytic
/// cylinder axis. One endpoint is the planar z=1 edge, the other is the
/// actual cylinder top circle (axial v=2.4). The remaining connected
/// curve passes through the artificial periodic chart seam u=0/2PI.
fn fixture()->(BrepBody,FaceId,BrepBody,FaceId){
    let block=BrepBody::from(
        Solid::block(p(-2.5,2.0,-2.2),5.0,1.0,3.2).unwrap());
    let model=block.shared().unwrap();
    let face=model.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=model.face_surface(f.id)? else{return None};
        ((s.frame.origin.y-2.0).abs()<1e-10 &&
         s.frame.normal.y.abs()>0.9).then_some(f.id)
    }).expect("missing physical horizontal plane face");

    let angle=0.4_f64.atan();
    let (si,co)=angle.sin_cos();
    let frame=Frame3::from_axes(p(0.0,0.0,0.0),
        p(co,si,0.0),p(0.0,0.0,-1.0),tol()).unwrap();
    let cylinder=BrepBody::from(CylindricalBrep::new(frame,2.0,2.4,tol()).unwrap());
    (block,face,cylinder,FaceId(2))
}

fn graph(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId)->rustsolid::IntersectionGraph{
    IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(a,fa,b,fb).unwrap()
}

#[test]
fn oblique_curve_is_clipped_at_real_edges_and_crosses_only_virtual_uv_seam(){
    let (block,plane_face,cylinder,side_face)=fixture();
    let a=block.shared().unwrap();
    let b=cylinder.shared().unwrap();
    let trimmed=TrimmedIntersectionEngine::new(tol()).unwrap()
        .intersect_faces(&a,plane_face,&b,side_face).unwrap();
    // The UV seam splits the analytic parameter intervals before
    // the graph coalesces them; it is not itself a physical endpoint.
    assert_eq!(trimmed.branches().len(),2,"{trimmed:?}");
    assert!(trimmed.branches().iter().all(|br|br.contact==TrimContact::Crossing));
    let g=graph(&a,plane_face,&b,side_face);
    assert_eq!(g.state,GraphState::Ready);
    assert_eq!(g.arcs.len(),1,"{g:?}");
    assert_eq!(g.nodes.len(),2);
    let arc=g.arcs[0];
    assert!(!arc.closed);
    assert!(arc.parameter_start>5.0,"{arc:?}");
    assert!(arc.parameter_end>TAU && arc.parameter_end<TAU+2.5);
    assert_eq!(arc.contact,TrimContact::Crossing);
    assert!(g.nodes.iter().all(|n|n.kind==GraphNodeKind::Endpoint));
    let seam_fraction=(TAU-arc.parameter_start)/(arc.parameter_end-arc.parameter_start);
    assert!(seam_fraction>0.0 && seam_fraction<1.0);
    let seam=arc.evaluate(seam_fraction).unwrap();
    let (s,c)=0.4_f64.atan().sin_cos();
    let expected_v=(2.0-2.0*s)/c;
    approx(seam.y,2.0);
    approx(seam.z,0.0);
    approx(seam.x,-2.0*s/c+2.0/c);
    assert!(expected_v>0.0 && expected_v<2.4);

    // Two different *physical* boundaries; neither is the cylindrical seam.
    let mut saw_plane_edge=false;
    let mut saw_cylinder_top=false;
    for node in &g.nodes {
        let q=node.point;
        if (q.z-1.0).abs()<2e-8 {
            saw_plane_edge=true;
            assert_eq!(node.membership,[TrimPosition::Boundary,TrimPosition::Interior]);
            assert!(node.boundary_sources.iter()
                .any(|s|s.operand==GraphOperand::First));
        } else {
            saw_cylinder_top=true;
            assert_eq!(node.membership,[TrimPosition::Interior,TrimPosition::Boundary]);
            let surface=b.face_surface(side_face).unwrap();
            let Surface3::Cylinder(cyl)=surface else{panic!("missing cylinder");};
            let relative=p(q.x-cyl.frame.origin.x,q.y-cyl.frame.origin.y,
                q.z-cyl.frame.origin.z);
            let axial=relative.x*cyl.frame.normal.x+
                relative.y*cyl.frame.normal.y+relative.z*cyl.frame.normal.z;
            approx(axial,2.4);
            assert!(node.boundary_sources.iter()
                .any(|s|s.operand==GraphOperand::Second && s.edge.0!=2));
        }
    }
    assert!(saw_plane_edge && saw_cylinder_top);
    g.validate_against(&a,plane_face,&b,side_face).unwrap();
}

#[test]
fn paired_analytic_pcurves_reconstruct_unwrapped_oblique_intersection(){
    let (block,fa,cylinder,fb)=fixture();
    let a=block.shared().unwrap();let b=cylinder.shared().unwrap();
    let g=graph(&a,fa,&b,fb);
    let arc=g.arcs[0];
    let [plane_pcurve,cylinder_pcurve]=arc.pcurves(&a,&b).unwrap();
    assert!(matches!(plane_pcurve.surface,Surface3::Plane(_)));
    assert!(matches!(cylinder_pcurve.surface,Surface3::Cylinder(_)));
    let mut previous_u=0.0;
    for (i,f) in [0.0,0.05,0.2,0.4,0.6,0.8,0.95,1.0].into_iter().enumerate(){
        let xyz=arc.evaluate(f).unwrap();
        let p0=plane_pcurve.reconstruct(f).unwrap();
        let p1=cylinder_pcurve.reconstruct(f).unwrap();
        assert!(distance(xyz,p0)<1e-8);
        assert!(distance(xyz,p1)<1e-8);
        let first_uv=plane_pcurve.evaluate(f).unwrap();
        let second_uv=cylinder_pcurve.evaluate(f).unwrap();
        assert!(first_uv.u.is_finite() && first_uv.v.is_finite());
        assert!(second_uv.u.is_finite() && second_uv.v.is_finite());
        if i>0 {assert!(second_uv.u>previous_u);}
        previous_u=second_uv.u;
    }
    let seam_fraction=(TAU-arc.parameter_start)/(arc.parameter_end-arc.parameter_start);
    approx(cylinder_pcurve.evaluate(seam_fraction).unwrap().u,TAU);
    assert!(cylinder_pcurve.evaluate(1.0).unwrap().u>TAU);
    assert!(cylinder_pcurve.evaluate(1.01).is_err());
    assert!(plane_pcurve.evaluate(-0.01).is_err());

    // Swapping operand order must preserve the geometric locus and merely
    // exchange the two different supporting charts.
    let swapped=graph(&b,fb,&a,fa);
    assert_eq!(swapped.arcs.len(),1);
    let rev=swapped.arcs[0];
    for f in [0.0,0.2,0.5,0.8,1.0] {
        assert!(distance(arc.evaluate(f).unwrap(),rev.evaluate(f).unwrap())<1e-8);
    }
    let [cyl_first,plane_second]=rev.pcurves(&b,&a).unwrap();
    assert!(matches!(cyl_first.surface,Surface3::Cylinder(_)));
    assert!(matches!(plane_second.surface,Surface3::Plane(_)));
}

#[test]
fn oblique_crossing_classifies_preview_but_never_commits_topology(){
    let (block,fa,cylinder,fb)=fixture();
    let a=block.shared().unwrap();let b=cylinder.shared().unwrap();
    let g=graph(&a,fa,&b,fb);
    let regions=BooleanRegionClassifier::classify_graph(
        &g,&block,fa,&cylinder,fb,BooleanOperation::Difference).unwrap();
    assert_eq!(regions.state,GraphRegionState::ClassifiedPreview,"{regions:?}");
    assert_eq!(regions.samples.len(),6);
    assert!(regions.samples.iter().all(|s|
        s.positive_other!=s.negative_other &&
        s.keep_positive!=s.keep_negative));
    assert!(regions.samples.iter().filter(|s|s.operand==GraphOperand::Second)
        .all(|s|s.reverse_orientation));
    regions.validate_against(&g,&block,fa,&cylinder,fb).unwrap();
    let plan=TopologyDeltaPlanner::new().plan(&g,&a,fa,&b,fb).unwrap();
    assert_eq!(plan.state,PlanState::PreviewOnly);
    assert!(!plan.can_apply);
    assert!(plan.requirements.contains(&PlanRequirement::GeneralizedBrepEditor));
    assert!(plan.requirements.contains(&PlanRequirement::OpenTraceCompletion));
    assert!(plan.requirements.contains(&PlanRequirement::CurvedEdgeSubdivision));
    a.validate().unwrap();
    b.validate().unwrap();
}

fn inclined_face_fixture()->(BrepBody,FaceId,BrepBody,FaceId){
    let angle=0.4_f64.atan();
    let (si,co)=angle.sin_cos();
    let constant=2.0*co-1.2*si;
    let profile=[
        Point2{x:constant+0.2,z:-0.2},
        Point2{x:constant-0.8,z:0.8},
        Point2{x:constant-3.8,z:-2.2},
        Point2{x:constant-2.8,z:-3.2},
    ];
    let plane=BrepBody::from(Solid::extrude_xz(&profile,4.0).unwrap());
    let plane_model=plane.shared().unwrap();
    let face=plane_model.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=plane_model.face_surface(f.id)? else{return None};
        ((s.frame.origin.x+s.frame.origin.z-constant).abs()<1e-8 &&
         s.frame.normal.x.abs()>0.5 && s.frame.normal.z.abs()>0.5)
            .then_some(f.id)
    }).expect("inclined physical plane face missing");
    let cylinder_frame=Frame3::from_axes(p(0.0,0.0,0.0),
        p(co,si,0.0),p(0.0,0.0,-1.0),tol()).unwrap();
    let cylinder=BrepBody::from(
        CylindricalBrep::new(cylinder_frame,2.0,2.4,tol()).unwrap());
    (plane,face,cylinder,FaceId(2))
}

#[test]
fn physically_inclined_finite_plane_has_dual_uv_trims_and_real_endpoints(){
    let (plane,fa,cylinder,fb)=inclined_face_fixture();
    let a=plane.shared().unwrap();let b=cylinder.shared().unwrap();
    let g=graph(&a,fa,&b,fb);
    assert_eq!(g.state,GraphState::Ready,"{g:?}");
    assert_eq!(g.arcs.len(),1,"{g:?}");
    assert_eq!(g.nodes.len(),2);
    let arc=g.arcs[0];
    assert!(arc.parameter_start>5.0 && arc.parameter_start<TAU);
    assert!(arc.parameter_end>TAU && arc.parameter_end<TAU+1.0);
    let pc=arc.pcurves(&a,&b).unwrap();
    let seam=(TAU-arc.parameter_start)/(arc.parameter_end-arc.parameter_start);
    assert!(seam>0.0 && seam<1.0);
    let mut previous=-f64::INFINITY;
    for t in [0.,0.1,0.25,0.5,0.75,0.9,1.] {
        let q=arc.evaluate(t).unwrap();
        let first=pc[0].reconstruct(t).unwrap();
        let second=pc[1].reconstruct(t).unwrap();
        assert!(distance(q,first)<1e-8);
        assert!(distance(q,second)<1e-8);
        let uv=pc[1].evaluate(t).unwrap();
        assert!(uv.u>previous);
        previous=uv.u;
    }
    approx(pc[1].evaluate(seam).unwrap().u,TAU);
    let mut plane_edge=false;
    let mut cylinder_edge=false;
    for n in &g.nodes {
        assert!(!n.boundary_sources.is_empty());
        if (n.point.z+0.2).abs()<1e-8 {
            plane_edge=true;
            assert_eq!(n.membership,[TrimPosition::Boundary,TrimPosition::Interior]);
        } else {
            cylinder_edge=true;
            assert_eq!(n.membership,[TrimPosition::Interior,TrimPosition::Boundary]);
        }
    }
    assert!(plane_edge && cylinder_edge);
    let r=BooleanRegionClassifier::classify_graph(
        &g,&plane,fa,&cylinder,fb,BooleanOperation::Intersection).unwrap();
    assert_eq!(r.state,GraphRegionState::ClassifiedPreview,"{r:?}");
    assert_eq!(r.samples.len(),6);
    let plan=TopologyDeltaPlanner::new().plan(&g,&a,fa,&b,fb).unwrap();
    assert_eq!(plan.state,PlanState::PreviewOnly);
    assert!(!plan.can_apply);
}

#[test]
fn near_axial_trim_tangency_is_held_without_fictitious_cut_edges(){
    let (block,fa,_,_)=fixture();
    let a=block.shared().unwrap();
    let angle=0.4_f64.atan();
    let (si,co)=angle.sin_cos();
    let frame=Frame3::from_axes(p(0.0,0.0,0.0),
        p(co,si,0.0),p(0.0,0.0,-1.0),tol()).unwrap();
    // Analytically the section's maximum axial coordinate occurs at u=PI.
    // A sub-tolerance perturbation of the cylinder cap creates a near
    // grazing event. Even if it is remote from the selected planar face,
    // the conservative surface/trim solver must not invent a cut.
    let height=(2.0+2.0*si)/co-4e-10;
    let cylinder=BrepBody::from(
        CylindricalBrep::new(frame,2.0,height,tol()).unwrap());
    let b=cylinder.shared().unwrap();
    let graph=graph(&a,fa,&b,FaceId(2));
    assert_eq!(graph.state,GraphState::Ambiguous(TrimAmbiguity::UnresolvedBoundaryEvent));
    assert!(graph.nodes.is_empty() && graph.arcs.is_empty());
    let region=BooleanRegionClassifier::classify_graph(
        &graph,&block,fa,&cylinder,FaceId(2),BooleanOperation::Difference).unwrap();
    assert_eq!(region.state,GraphRegionState::Held);
    assert!(region.samples.is_empty() && !region.can_apply());
    let topology=TopologyDeltaPlanner::new().plan(&graph,&a,fa,&b,FaceId(2)).unwrap();
    assert_eq!(topology.state,PlanState::Held);
    assert!(!topology.can_apply && topology.arcs.is_empty());
}
