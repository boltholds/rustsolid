//! G4b2b2: full source graph across many face pairs, globally certified
//! G1 endpoint identity and two real analytic source-face pcurves per branch.
//! The supported partial cylinder cutting curves remain *preview* only.
use rustsolid::{
    AnalyticCutRegistryHold,BrepBody,CylindricalBrep,GeometryTolerance,
    GlobalAnalyticCutPlan,GlobalAnalyticCutResult,GlobalEventKind,
    GraphOperand,Point3,Solid,Curve2,Curve3,
};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn cutter()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}
fn intersected_box()->BrepBody{
    BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
}
fn run(a:&BrepBody,b:&BrepBody)->GlobalAnalyticCutResult{
    GlobalAnalyticCutPlan::build(a,b,tol()).unwrap()
}
fn ready(s:GlobalAnalyticCutResult)->GlobalAnalyticCutPlan{
    let GlobalAnalyticCutResult::Ready(plan)=s else{
        panic!("expected G1 source-proven partial cylinder cuts: {s:?}");
    };
    plan
}
fn dist(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

#[test]
fn partial_cylinder_and_box_reuse_real_g1_shared_vertex_events(){
    let a=intersected_box();let b=cutter();
    let original_a=format!("{a:?}");
    let original_b=format!("{b:?}");
    let plan=ready(run(&a,&b));
    let am=a.shared().unwrap();
    let bm=b.shared().unwrap();
    assert_eq!(plan.audited_face_pairs,am.faces.len()*bm.faces.len());
    assert!(!plan.arcs.is_empty());
    assert!(plan.physical_event_count>=2);
    assert!(plan.global_nodes.events.iter().any(|event|
        event.kind==GlobalEventKind::BoundaryEndpoint &&
        event.members.len()>1),
        "G1 must identify genuine physical endpoints across adjacent source faces");
    assert!(!plan.can_apply && !plan.closed_shell_proven);
    assert!(plan.arcs.iter().any(|arc|
        arc.physical_endpoints.is_some()&&!arc.virtual_periodic_anchor));
    for arc in &plan.arcs{
        assert!(!arc.virtual_periodic_anchor || arc.physical_endpoints.is_none());
        assert!(matches!(arc.carrier,Curve3::Line(_)|Curve3::Harmonic(_)));
        for (face,surface) in [
            (arc.first_face,am.face_surface(arc.first_face).unwrap()),
            (arc.second_face,bm.face_surface(arc.second_face).unwrap())
        ] {
            let _=face;
            for fraction in [0.0_f64,0.1,0.25,0.5,0.75,0.9,1.0] {
                let rustsolid::ParameterRange::Bounded{start,end}=arc.parameter_domain
                    else{panic!("graph arc must have bounded carrier")};
                let t=start+(end-start)*fraction;
                let c2=if surface==am.face_surface(arc.first_face).unwrap(){
                    arc.face_pcurves[0]
                }else{arc.face_pcurves[1]};
                let uv=c2.evaluate(t).unwrap();
                let from_pcurve=surface.evaluate(uv.u,uv.v).unwrap().point;
                let from_carrier=arc.carrier.evaluate(t).unwrap();
                assert!(dist(from_pcurve,from_carrier)<1e-7);
            }
        }
    }
    plan.validate_against(&a,&b,tol()).unwrap();
    assert_eq!(format!("{a:?}"),original_a);
    assert_eq!(format!("{b:?}"),original_b);
}

#[test]
fn two_chart_pcurves_are_supported_for_line_and_harmonic_cut_families(){
    let plan=ready(run(&intersected_box(),&cutter()));
    assert!(plan.arcs.iter().any(|arc|
        matches!(arc.face_pcurves[0],Curve2::Line(_)|Curve2::Harmonic(_))));
    assert!(plan.arcs.iter().any(|arc|
        matches!(arc.face_pcurves[1],
            Curve2::Line(_)|Curve2::CylinderSection(_)|Curve2::Harmonic(_))));
    assert!(!plan.can_apply);
}

#[test]
fn tampered_graph_arc_or_global_event_is_rejected_by_source_rebuild(){
    let a=intersected_box();let b=cutter();
    let plan=ready(run(&a,&b));
    let mut tampered=plan.clone();
    tampered.arcs[0].source.pair=rustsolid::FacePairId(999);
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut tampered=plan.clone();
    tampered.global_nodes.events[0].point.x+=0.2;
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
}

#[test]
fn stale_polyhedral_revision_and_cylinder_geometry_invalidate_registry(){
    let a=intersected_box();let b=cutter();
    let plan=ready(run(&a,&b));
    let mut modified=a.clone();
    if let BrepBody::Polyhedral(ref mut s)=modified{s.revision+=1;}
    assert!(plan.validate_against(&modified,&b,tol()).is_err());
    let different=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),2.1,5.,tol()).unwrap());
    assert!(plan.validate_against(&a,&different,tol()).is_err());
    let changed=GeometryTolerance {
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8,
    };
    assert!(plan.validate_against(&a,&b,changed).is_err());
}

#[test]
fn source_empty_and_uncertain_tangency_fail_closed_without_fictitious_cut_edges(){
    let disjoint=BrepBody::from(Solid::block(
        p(10.,10.,10.),2.,2.,2.).unwrap());
    assert!(matches!(run(&disjoint,&cutter()),GlobalAnalyticCutResult::Held(
        AnalyticCutRegistryHold::NoCertifiedCrossing)));
    let nearly_tangent=BrepBody::from(Solid::block(
        p(2.+4e-10,1.,-3.),2.,2.,6.).unwrap());
    assert!(matches!(run(&nearly_tangent,&cutter()),GlobalAnalyticCutResult::Held(
        AnalyticCutRegistryHold::UpstreamGraphAmbiguity)));
}
