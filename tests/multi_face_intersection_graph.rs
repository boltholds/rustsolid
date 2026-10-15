use rustsolid::{
    BrepBody, GeometryTolerance, MultiFaceIntersectionGraphBuilder,
    MultiFaceGraphState, FacePairId, GlobalNodeRef, GraphNodeId,
    Point3, Solid,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn eps()->GeometryTolerance{GeometryTolerance::default()}
#[test]
fn disjoint_blocks_have_no_face_intersections(){
    let a=BrepBody::from(Solid::block(p(0.,0.,0.),2.,2.,2.).unwrap());
    let b=BrepBody::from(Solid::block(p(10.,0.,0.),2.,2.,2.).unwrap());
    let first=a.shared().unwrap();
    let second=b.shared().unwrap();
    let graph=MultiFaceIntersectionGraphBuilder::new(GeometryTolerance::default())
        .unwrap().build(&first,&second).unwrap();
    assert_eq!(graph.pairs.len(),first.faces.len()*second.faces.len());
    assert_eq!(graph.state,MultiFaceGraphState::Empty);
    graph.validate_against(&first,&second).unwrap();
    assert!(graph.node(GlobalNodeRef{pair:FacePairId(999),local:GraphNodeId(0)}).is_none());
}
#[test]
fn tampered_pair_provenance_is_rejected(){
    let a=BrepBody::from(Solid::block(p(0.,0.,0.),2.,2.,2.).unwrap());
    let b=BrepBody::from(Solid::block(p(10.,0.,0.),2.,2.,2.).unwrap());
    let first=a.shared().unwrap();
    let second=b.shared().unwrap();
    let mut graph=MultiFaceIntersectionGraphBuilder::new(GeometryTolerance::default())
        .unwrap().build(&first,&second).unwrap();
    graph.pairs[0].id=FacePairId(1);
    assert!(graph.validate_against(&first,&second).is_err());
}

#[test]
fn falsified_empty_pair_is_rejected_by_source_rederivation() {
    use rustsolid::{CylindricalBrep, GraphState};
    let first=BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
        .shared().unwrap();
    let second=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),2.,5.,eps()).unwrap()).shared().unwrap();
    let mut graph=MultiFaceIntersectionGraphBuilder::new(eps()).unwrap()
        .build(&first,&second).unwrap();
    let pair=graph.pairs.iter_mut().find(|pair|pair.graph.state==GraphState::Ready)
        .expect("transverse test geometry contains a ready pair");
    pair.graph.state=GraphState::Empty;
    pair.graph.nodes.clear();
    pair.graph.arcs.clear();
    pair.graph.support_certificate=None;
    // Even if the attacker avoids a state inconsistency by changing all the
    // other pair states, source rederivation must reject the missing branch.
    assert!(graph.validate_against(&first,&second).is_err());
}
