//! G4b2b3a: actual source-edge subdivision, not a Boolean result shell.
use rustsolid::{
    AnalyticSourceEdgeSplitHold,BrepBody,CylindricalBrep,
    GeometryTolerance,GlobalAnalyticEdgeSplitPlan,GlobalAnalyticEdgeSplitResult,
    GlobalEventKind,GraphOperand,ParameterRange,Point3,Solid,Curve3,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn box_body()->BrepBody{
    BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
}
fn cutter()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}
fn plan(a:&BrepBody,b:&BrepBody)->GlobalAnalyticEdgeSplitPlan{
    let result=GlobalAnalyticEdgeSplitPlan::build(a,b,tol()).unwrap();
    let GlobalAnalyticEdgeSplitResult::Ready(plan)=result else{
        panic!("expected certified source topology splits: {result:?}");
    };
    plan
}
fn distance(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

#[test]
fn partial_box_cylinder_physically_subdivides_two_closed_source_shells(){
    let (a,b)=(box_body(),cutter());
    let source_a=a.shared().unwrap();
    let source_b=b.shared().unwrap();
    let original_a=format!("{a:?}");
    let original_b=format!("{b:?}");
    let split=plan(&a,&b);
    assert!(!split.edges.is_empty());
    assert!(split.new_source_vertices.iter().sum::<usize>()>0);
    assert!(!split.can_apply && !split.boolean_shell_proven);
    for (source,rebuilt) in [source_a,source_b]
        .iter().zip(split.split_sources.iter())
    {
        rebuilt.validate().unwrap();
        assert_eq!(rebuilt.faces.len(),source.faces.len());
        assert_eq!(rebuilt.loops.len(),source.loops.len());
        assert_eq!(rebuilt.shells.len(),source.shells.len());
        assert_eq!(rebuilt.euler_characteristic(),source.euler_characteristic());
        let delta=rebuilt.vertices.len()-source.vertices.len();
        assert_eq!(rebuilt.edges.len()-source.edges.len(),delta);
        assert_eq!(rebuilt.coedges.len()-source.coedges.len(),2*delta);
        assert_eq!(rebuilt.mass,source.mass);
    }
    assert_eq!(format!("{a:?}"),original_a);
    assert_eq!(format!("{b:?}"),original_b);
    split.validate_against(&a,&b,tol()).unwrap();
}

#[test]
fn source_edge_piece_domains_preserve_two_true_opposite_uv_trim_charts(){
    let (a,b)=(box_body(),cutter());
    let split=plan(&a,&b);
    let originals=[a.shared().unwrap(),b.shared().unwrap()];
    for edge in &split.edges {
        let idx=if edge.operand==GraphOperand::First{0}else{1};
        let source=&originals[idx];
        let produced=&split.split_sources[idx];
        let original=&source.edges[edge.original_edge.0 as usize];
        assert!(edge.pieces.len()>=2);
        assert_eq!(original.coedges,edge.original_coedges);
        assert_eq!(edge.interior_events.len()+1,edge.pieces.len());
        let mut total=0.0;
        for (i,piece) in edge.pieces.iter().enumerate(){
            let actual=&produced.edges[piece.result_edge.0 as usize];
            assert_eq!(actual.coedges,piece.result_coedges);
            assert_eq!([actual.start,actual.end],piece.vertices);
            assert!(piece.source_fractions[0]<piece.source_fractions[1]);
            if i==0 {assert_eq!(piece.source_fractions[0],0.0);}
            if i+1==edge.pieces.len(){assert_eq!(piece.source_fractions[1],1.0);}
            total+=piece.source_fractions[1]-piece.source_fractions[0];
            let ParameterRange::Bounded{start,end}=piece.curve_domain else {
                panic!("source edge must be bounded")
            };
            assert!(end>start);
            let carrier=produced.edge_curve(piece.result_edge).unwrap();
            assert!(matches!(carrier,Curve3::Line(_)));
            for fraction in [0.0,0.25,0.5,0.75,1.0]{
                let xyz=carrier.evaluate(start+(end-start)*fraction).unwrap();
                for &id in &piece.result_coedges{
                    let co=&produced.coedges[id.0 as usize];
                    let trim=produced.coedge_pcurve(id).unwrap();
                    let geometry=produced.geometry.coedges[id.0 as usize];
                    let ParameterRange::Bounded{start:uvstart,end:uvend}=geometry.domain
                        else{panic!("coedge UV domain must be bounded")};
                    let used_fraction=if co.reversed{1.0-fraction}else{fraction};
                    let uv=trim.evaluate(uvstart+(uvend-uvstart)*used_fraction).unwrap();
                    let surface=produced.face_surface(co.face).unwrap();
                    let projected=surface.evaluate(uv.u,uv.v).unwrap().point;
                    assert!(distance(xyz,projected)<1e-7);
                }
            }
        }
        assert!((total-1.0).abs()<1e-12);
    }
}

#[test]
fn one_global_event_binding_is_reused_across_source_face_pair_witnesses(){
    let (a,b)=(box_body(),cutter());
    let split=plan(&a,&b);
    assert!(split.registry.global_nodes.events.iter().any(|event|
        event.kind==GlobalEventKind::BoundaryEndpoint && event.members.len()>1));
    let mut unique=std::collections::BTreeSet::new();
    for binding in &split.event_bindings{
        assert!(unique.insert((binding.operand,binding.event)));
        assert!(!binding.source_edges.is_empty());
        let model=&split.split_sources[if binding.operand==GraphOperand::First{0}else{1}];
        assert!(model.vertices.get(binding.result_vertex.0 as usize).is_some());
        for edge in &binding.source_edges {
            let parent=&[a.shared().unwrap(),b.shared().unwrap()]
                [if binding.operand==GraphOperand::First{0}else{1}].edges[edge.0 as usize];
            assert!(parent.id==*edge);
        }
    }
    let created=split.event_bindings.iter().filter(|b|b.created_vertex).count();
    assert_eq!(created,split.new_source_vertices.iter().sum::<usize>());
}

#[test]
fn tampering_with_result_topology_or_g1_identity_invalidates_preview(){
    let (a,b)=(box_body(),cutter());
    let source=plan(&a,&b);
    let mut tampered=source.clone();
    let affected=&source.edges[0];
    let idx=if affected.operand==GraphOperand::First{0}else{1};
    let coedge=affected.pieces[0].result_coedges[0];
    tampered.split_sources[idx].coedges[coedge.0 as usize].twin =
        rustsolid::CoedgeId(u32::MAX);
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut tampered=source.clone();
    tampered.event_bindings[0].result_vertex=rustsolid::VertexId(u32::MAX);
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut tampered=source.clone();
    tampered.registry.arcs[0].source.pair=rustsolid::FacePairId(u32::MAX);
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
}

#[test]
fn source_revision_geometry_and_tolerance_changes_are_rejected(){
    let (a,b)=(box_body(),cutter());
    let split=plan(&a,&b);
    let mut modified=a.clone();
    if let BrepBody::Polyhedral(ref mut s)=modified {s.revision+=1;}
    assert!(split.validate_against(&modified,&b,tol()).is_err());
    let different=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),2.1,5.,tol()).unwrap());
    assert!(split.validate_against(&a,&different,tol()).is_err());
    let changed=GeometryTolerance{
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8
    };
    assert!(split.validate_against(&a,&b,changed).is_err());
}

#[test]
fn disjoint_sources_are_held_without_modifying_any_brep(){
    let a=BrepBody::from(Solid::block(p(15.,15.,15.),2.,2.,2.).unwrap());
    let b=cutter();
    let result=GlobalAnalyticEdgeSplitPlan::build(&a,&b,tol()).unwrap();
    assert!(matches!(result,GlobalAnalyticEdgeSplitResult::Held(
        AnalyticSourceEdgeSplitHold::UpstreamRegistry(_))));
}
