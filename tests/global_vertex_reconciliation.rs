//! Integration coverage for v0.41 G1: pair-scoped events and topology-led
//! endpoint reconciliation. This is read-only: no candidate B-rep is emitted.
use rustsolid::{
    BrepBody, CylindricalBrep, GeometryTolerance, GlobalNodeReconciliation,
    GlobalNodeReconciliationHold, MultiFaceGraphState,
    MultiFaceIntersectionGraphBuilder, Point3, Solid,
};

fn p(x:f64,y:f64,z:f64)->Point3 { Point3{x,y,z} }
fn tol()->GeometryTolerance { GeometryTolerance::default() }
fn block(origin:Point3,dx:f64,dy:f64,dz:f64)->BrepBody {
    BrepBody::from(Solid::block(origin,dx,dy,dz).unwrap())
}

#[test]
fn disjoint_solids_produce_empty_certified_global_node_plan() {
    let a=block(p(10.,10.,10.),2.,2.,2.).shared().unwrap();
    let b=BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol())
        .unwrap()).shared().unwrap();
    let graph=MultiFaceIntersectionGraphBuilder::new(tol()).unwrap()
        .build(&a,&b).unwrap();
    assert_eq!(graph.state,MultiFaceGraphState::Empty);
    let result=GlobalNodeReconciliation::build(&graph,&a,&b).unwrap();
    let GlobalNodeReconciliation::Ready(plan)=&result else {
        panic!("disjoint solids should have no events: {result:?}");
    };
    assert!(plan.events.is_empty());
    assert_eq!(plan.source_local_nodes,0);
    assert!(!plan.can_apply);
    result.validate_against(&graph,&a,&b).unwrap();
}

#[test]
fn overlapping_box_and_cylinder_reconcile_shared_physical_endpoints() {
    let a=block(p(-1.5,1.,-1.5),3.,3.,3.).shared().unwrap();
    let b=BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol())
        .unwrap()).shared().unwrap();
    let graph=MultiFaceIntersectionGraphBuilder::new(tol()).unwrap()
        .build(&a,&b).unwrap();
    assert_eq!(graph.state,MultiFaceGraphState::Ready);
    let result=GlobalNodeReconciliation::build(&graph,&a,&b).unwrap();
    let GlobalNodeReconciliation::Ready(plan)=&result else {
        panic!("supported transverse adjacency must reconcile: {result:?}");
    };
    let shared:Vec<_>=plan.events.iter().filter(|event|event.members.len()>1).collect();
    assert!(!shared.is_empty(),"expected endpoints shared by adjacent box faces");
    assert!(shared.iter().all(|event|event.members.iter().all(|member|
        graph.node(*member).is_some() && plan.event_for(*member)==Some(*event))));
    assert!(!plan.can_apply);
    result.validate_against(&graph,&a,&b).unwrap();

    let mut changed=result.clone();
    if let GlobalNodeReconciliation::Ready(ref mut edited)=changed{
        edited.events[0].point.x+=0.1;
    }
    assert!(changed.validate_against(&graph,&a,&b).is_err());
}

#[test]
fn uncertain_cylinder_tangency_does_not_leak_global_events() {
    let a=block(p(2.+4e-10,1.,-3.),2.,2.,6.).shared().unwrap();
    let b=BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol())
        .unwrap()).shared().unwrap();
    let graph=MultiFaceIntersectionGraphBuilder::new(tol()).unwrap()
        .build(&a,&b).unwrap();
    assert_eq!(graph.state,MultiFaceGraphState::Held);
    let result=GlobalNodeReconciliation::build(&graph,&a,&b).unwrap();
    assert_eq!(result,GlobalNodeReconciliation::Held(
        GlobalNodeReconciliationHold::UpstreamAmbiguity));
    result.validate_against(&graph,&a,&b).unwrap();
}
