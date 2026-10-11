//! An executable correctness gate; external CAD engines are never oracles.
use rustsolid::{
    ALL_MODELING_OPERATORS, ModelingOperator, CommandBatch, CommandHistory,
    EdgeId, EulerCommand, FaceId, Point3, ProbeExpectation, ProbeModel,
    RegressionCorpus, Solid, TopologyEntity, TopologyName, VertexId,
    run_builtin_regressions, BUILTIN_REGRESSIONS,
};
use serde::Deserialize;
use std::collections::{BTreeMap,BTreeSet};

const SCHEMA:&str="rustsolid.operator_contracts.v1";
const CONTRACTS:&str=include_str!("fixtures/operator_contracts.v1.json");

#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {schema:String,title:String,operators:Vec<Contract>}
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct Contract {operator:ModelingOperator,scenarios:Vec<Scenario>}
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Deserialize)]
#[serde(rename_all="snake_case")]
enum Class {Analytic,Metamorphic,Boundary}
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
struct Scenario {
    id:String,category:Class,
    fixture:Option<String>,driver:Option<String>,variant:Option<String>,
}

fn check_manifest(raw:&str)->Result<Manifest,String>{
    let m:Manifest=serde_json::from_str(raw).map_err(|e|e.to_string())?;
    if m.schema!=SCHEMA || m.title.trim().is_empty(){return Err("invalid inventory schema/title".into());}
    let declared=ALL_MODELING_OPERATORS.iter().copied().collect::<BTreeSet<_>>();
    let actual=m.operators.iter().map(|o|o.operator).collect::<BTreeSet<_>>();
    if actual!=declared || m.operators.len()!=declared.len(){return Err("operator inventory missing, duplicated or unregistered".into());}
    let mut ids=BTreeSet::new();
    for op in &m.operators {
        let mut counts=BTreeMap::new();
        for sc in &op.scenarios {
            if !sc.id.starts_with(&format!("{}/",op.operator.as_str())) || !ids.insert(sc.id.clone()) {
                return Err(format!("invalid/duplicate scenario ID {}",sc.id));
            }
            if op.operator.is_command() {
                if sc.fixture.is_some() || sc.driver.as_deref()!=Some("euler_cuboid") {
                    return Err(format!("{} needs executable Euler driver",sc.id));
                }
                if (sc.category==Class::Boundary)!=sc.variant.is_some() {
                    return Err(format!("{} has a missing or unexpected negative variant",sc.id));
                }
            } else if sc.fixture.is_none() || sc.driver.is_some() || sc.variant.is_some() {
                return Err(format!("{} must refer to a versioned regression fixture",sc.id));
            }
            *counts.entry(sc.category).or_insert(0usize)+=1;
        }
        for (cat,min) in [(Class::Analytic,1),(Class::Metamorphic,1),(Class::Boundary,2)] {
            if counts.get(&cat).copied().unwrap_or(0)<min {
                return Err(format!("{} lacks {cat:?} coverage",op.operator.as_str()));
            }
        }
    }
    Ok(m)
}

fn primitive_has_correct_shape(op:ModelingOperator,m:&ProbeModel)->bool{
    matches!((op,m),
       (ModelingOperator::Block,ProbeModel::Block{..})
      |(ModelingOperator::ExtrudeXz,ProbeModel::Extrusion{..})
      |(ModelingOperator::AnalyticCylinder,ProbeModel::Cylinder{..}))
}
fn check_primitive(op:ModelingOperator,s:&Scenario,fixtures:&RegressionCorpus,
                   reports:&BTreeMap<String,bool>)->Result<(),String>{
    let id=s.fixture.as_deref().ok_or("missing fixture")?;
    let f=fixtures.fixtures.iter().find(|f|f.id==id).ok_or_else(||format!("unknown fixture {id}"))?;
    if !primitive_has_correct_shape(op,&f.model){return Err("fixture uses wrong primitive".into());}
    if reports.get(id)!=Some(&true){return Err(format!("mathematical fixture {id} failed"));}
    let supported=match s.category {
        Class::Analytic=>f.expectation==ProbeExpectation::Accept&&f.golden.is_some(),
        Class::Metamorphic=>f.expectation==ProbeExpectation::Accept&&!f.relations.is_empty(),
        Class::Boundary=>f.expectation==ProbeExpectation::Reject,
    };
    if !supported{return Err(format!("fixture {id} is not a {0:?} witness",s.category));}
    Ok(())
}

fn name(h:&CommandHistory,e:TopologyEntity)->TopologyName{
    h.solid().topology_name(e).unwrap().clone()
}
fn cube(origin:[f64;3],dims:[f64;3])->CommandHistory{
    let p=Point3{x:origin[0],y:origin[1],z:origin[2]};
    let s=Solid::block(p,dims[0],dims[1],dims[2]).unwrap()
        .with_feature_key("operator/contracts").unwrap();
    CommandHistory::new(s).unwrap()
}
fn command(op:ModelingOperator,h:&CommandHistory)->EulerCommand{
    let edge=name(h,TopologyEntity::Edge(EdgeId(0)));
    let face=name(h,TopologyEntity::Face(FaceId(0)));
    let ring=&h.solid().faces[0].boundary;
    let a=name(h,TopologyEntity::Vertex(ring[0]));
    let b=name(h,TopologyEntity::Vertex(ring[2]));
    match op {
        ModelingOperator::SplitEdge=>EulerCommand::SplitEdge{edge,fraction:0.5},
        ModelingOperator::MakeEdgeVertex=>EulerCommand::MakeEdgeVertex{edge,fraction:0.5},
        ModelingOperator::SplitFace=>EulerCommand::SplitFace{face,start:a,end:b},
        ModelingOperator::MakeEdgeFace=>EulerCommand::MakeEdgeFace{face,start:a,end:b},
        ModelingOperator::KillEdgeVertex=>EulerCommand::KillEdgeVertex{
            vertex:name(h,TopologyEntity::Vertex(VertexId(8)))},
        ModelingOperator::KillEdgeFace=>EulerCommand::KillEdgeFace{
            edge:name(h,TopologyEntity::Edge(EdgeId(12))),
            removed_face:name(h,TopologyEntity::Face(FaceId(6)))},
        _=>panic!("primitive is not a command"),
    }
}
fn prelude(op:ModelingOperator,h:&mut CommandHistory){
    let make=match op {
        ModelingOperator::KillEdgeVertex=>Some(ModelingOperator::MakeEdgeVertex),
        ModelingOperator::KillEdgeFace=>Some(ModelingOperator::MakeEdgeFace),
        _=>None,
    };
    if let Some(kind)=make {
        let cmd=command(kind,h);
        h.execute(CommandBatch::single(cmd)).unwrap();
    }
}
fn approx(a:f64,b:f64){
    assert!((a-b).abs()<=1.0e-9*b.abs().max(1.0),"actual={a}, oracle={b}");
}
fn verify_closed_form(body:&Solid,o:[f64;3],d:[f64;3]){
    let [w,h,depth]=d;
    approx(body.mass.volume,w*h*depth);
    approx(body.mass.surface_area,2.0*(w*h+w*depth+h*depth));
    for (actual,target) in [body.mass.centroid.x,body.mass.centroid.y,body.mass.centroid.z]
        .into_iter().zip([o[0]+w/2.0,o[1]+h/2.0,o[2]+depth/2.0]) {
        approx(actual,target);
    }
    for (actual,target) in [body.bbox.min.x,body.bbox.min.y,body.bbox.min.z]
        .into_iter().zip(o){approx(actual,target);}
    for (actual,target) in [body.bbox.max.x,body.bbox.max.y,body.bbox.max.z]
        .into_iter().zip([o[0]+w,o[1]+h,o[2]+depth]){approx(actual,target);}
    assert_eq!(body.euler_characteristic(),2);
    body.validate().unwrap();
}
#[derive(Debug)]
struct Measurement {volume:f64,area:f64,counts:(usize,usize,usize)}
fn valid_euler(op:ModelingOperator,o:[f64;3],d:[f64;3])->Measurement{
    let mut h=cube(o,d);
    prelude(op,&mut h);
    let before=(h.solid().vertices.len(),h.solid().edges.len(),h.solid().faces.len());
    let rev=h.solid().revision;
    let c=command(op,&h);
    assert_eq!(c.operator_kind(),op);
    let report=h.execute(CommandBatch::single(c)).unwrap();
    assert_eq!((report.revision_before,report.revision_after),(rev,rev+1));
    let body=h.solid();
    verify_closed_form(body,o,d);
    let counts=(body.vertices.len(),body.edges.len(),body.faces.len());
    let expected=match op {
        ModelingOperator::SplitEdge|ModelingOperator::MakeEdgeVertex=>(9,13,6),
        ModelingOperator::SplitFace|ModelingOperator::MakeEdgeFace=>(8,13,7),
        ModelingOperator::KillEdgeVertex|ModelingOperator::KillEdgeFace=>(8,12,6),
        _=>unreachable!(),
    };
    assert_eq!(counts,expected);
    let measured=Measurement {volume:body.mass.volume,area:body.mass.surface_area,counts};
    h.undo().unwrap().expect("expected undo");
    assert_eq!((h.solid().vertices.len(),h.solid().edges.len(),h.solid().faces.len()),before);
    verify_closed_form(h.solid(),o,d);
    h.redo().unwrap().expect("expected redo");
    assert_eq!((h.solid().vertices.len(),h.solid().edges.len(),h.solid().faces.len()),counts);
    verify_closed_form(h.solid(),o,d);
    measured
}
fn invalid_command(op:ModelingOperator,h:&CommandHistory,variant:&str)->Result<EulerCommand,String>{
    let edge=name(h,TopologyEntity::Edge(EdgeId(0)));
    let face=name(h,TopologyEntity::Face(FaceId(0)));
    let ring=&h.solid().faces[0].boundary;
    let a=name(h,TopologyEntity::Vertex(ring[0]));
    let b=name(h,TopologyEntity::Vertex(ring[if variant=="repeated"{0}else{1}]));
    Ok(match (op,variant) {
        (ModelingOperator::SplitEdge,"fraction_zero")=>EulerCommand::SplitEdge{edge,fraction:0.0},
        (ModelingOperator::SplitEdge,"fraction_near_endpoint")=>EulerCommand::SplitEdge{edge,fraction:1e-14},
        (ModelingOperator::MakeEdgeVertex,"fraction_zero")=>EulerCommand::MakeEdgeVertex{edge,fraction:0.0},
        (ModelingOperator::MakeEdgeVertex,"fraction_near_endpoint")=>EulerCommand::MakeEdgeVertex{edge,fraction:1e-14},
        (ModelingOperator::SplitFace,"adjacent"|"repeated")=>EulerCommand::SplitFace{face,start:a,end:b},
        (ModelingOperator::MakeEdgeFace,"adjacent"|"repeated")=>EulerCommand::MakeEdgeFace{face,start:a,end:b},
        (ModelingOperator::KillEdgeVertex,"ordinary_vertex")=>EulerCommand::KillEdgeVertex{
            vertex:name(h,TopologyEntity::Vertex(VertexId(0)))},
        (ModelingOperator::KillEdgeVertex,"unknown_vertex")=>EulerCommand::KillEdgeVertex{
            vertex:TopologyName("operator/contracts/vertex/not-found".into())},
        (ModelingOperator::KillEdgeFace,"noncoplanar")=>{
            let e=&h.solid().edges[0];
            let f=h.solid().coedges[e.coedges[0].0 as usize].face;
            EulerCommand::KillEdgeFace{edge,removed_face:name(h,TopologyEntity::Face(f))}
        },
        (ModelingOperator::KillEdgeFace,"unknown_face")=>EulerCommand::KillEdgeFace{
            edge,removed_face:TopologyName("operator/contracts/face/not-found".into())},
        _=>return Err(format!("unimplemented concrete negative variant {op:?}/{variant}")),
    })
}
fn euler_scenario(op:ModelingOperator,s:&Scenario)->Result<(),String>{
    match s.category {
        Class::Boundary=>{
            let mut h=cube([2.0,-7.0,9.0],[2.0,3.0,4.0]);
            let cmd=invalid_command(op,&h,s.variant.as_deref().unwrap_or(""))?;
            if cmd.operator_kind()!=op {return Err("operator mismatch".into());}
            let state=format!("{:?}",h.solid());let rev=h.solid().revision;
            let undo=h.undo_depth();let redo=h.redo_depth();
            if h.execute(CommandBatch::single(cmd)).is_ok(){return Err("invalid topology accepted".into());}
            if state!=format!("{:?}",h.solid()) || h.solid().revision!=rev ||
                h.undo_depth()!=undo || h.redo_depth()!=redo {
                return Err("failed edit mutated geometry/history".into());
            }
            h.solid().validate().map_err(|e|e.to_string())
        }
        Class::Analytic=>{
            valid_euler(op,[0.0,0.0,0.0],[2.0,3.0,4.0]);
            Ok(())
        }
        Class::Metamorphic=>{
            let base=valid_euler(op,[0.0,0.0,0.0],[2.0,3.0,4.0]);
            let moved=valid_euler(op,[125.0,-43.0,810.0],[2.0,3.0,4.0]);
            let scaled=valid_euler(op,[125.0,-43.0,810.0],[4.0,6.0,8.0]);
            assert_eq!(moved.counts,base.counts);
            assert_eq!(scaled.counts,base.counts);
            approx(moved.volume,base.volume);approx(moved.area,base.area);
            approx(scaled.volume,base.volume*8.0);approx(scaled.area,base.area*4.0);
            Ok(())
        }
    }
}
#[test]
fn every_operator_has_real_analytic_metamorphic_and_boundary_witnesses(){
    let manifest=check_manifest(CONTRACTS).unwrap();
    assert_eq!(manifest.operators.len(),9);
    assert_eq!(manifest.operators.iter().map(|s|s.scenarios.len()).sum::<usize>(),36);
    let fixtures:RegressionCorpus=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    let regressions=run_builtin_regressions().unwrap();
    assert!(regressions.all_passed(),"pinned mathematical regressions failed");
    let status=regressions.outcomes.iter()
        .map(|o|(o.fixture.id.clone(),o.passed)).collect::<BTreeMap<_,_>>();
    for op in &manifest.operators {
        for s in &op.scenarios {
            let result=if op.operator.is_command() {euler_scenario(op.operator,s)}
                else {check_primitive(op.operator,s,&fixtures,&status)};
            if let Err(error)=result {panic!("{} failed: {}",s.id,error);}
        }
    }
}
#[test]
fn gate_fails_on_missing_category_or_missing_operator(){
    let mut doc:serde_json::Value=serde_json::from_str(CONTRACTS).unwrap();
    doc["operators"][0]["scenarios"].as_array_mut().unwrap()
        .retain(|s|s["category"]!="boundary");
    assert!(check_manifest(&doc.to_string()).is_err());
    let mut doc:serde_json::Value=serde_json::from_str(CONTRACTS).unwrap();
    doc["operators"].as_array_mut().unwrap().pop();
    assert!(check_manifest(&doc.to_string()).is_err());
}
#[test]
fn forged_fixture_and_unknown_driver_cannot_satisfy_contract(){
    let mut doc:serde_json::Value=serde_json::from_str(CONTRACTS).unwrap();
    doc["operators"][0]["scenarios"][2]["fixture"]=serde_json::json!("known/block_2_3_4");
    let m=check_manifest(&doc.to_string()).unwrap();
    let fixed:RegressionCorpus=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    let reports=run_builtin_regressions().unwrap().outcomes.iter()
        .map(|r|(r.fixture.id.clone(),r.passed)).collect::<BTreeMap<_,_>>();
    assert!(check_primitive(m.operators[0].operator,&m.operators[0].scenarios[2],&fixed,&reports).is_err());
    let mut doc:serde_json::Value=serde_json::from_str(CONTRACTS).unwrap();
    doc["operators"][3]["scenarios"][0]["driver"]=serde_json::json!("no-op");
    assert!(check_manifest(&doc.to_string()).is_err());
}
#[test]
fn every_serialized_euler_command_has_an_exhaustive_kind_mapping(){
    let edge=TopologyName("edge".into());
    let face=TopologyName("face".into());
    let vertex=TopologyName("vertex".into());
    let commands=[
        EulerCommand::SplitEdge{edge:edge.clone(),fraction:0.5},
        EulerCommand::SplitFace{face:face.clone(),start:vertex.clone(),end:vertex.clone()},
        EulerCommand::MakeEdgeVertex{edge:edge.clone(),fraction:0.5},
        EulerCommand::MakeEdgeFace{face:face.clone(),start:vertex.clone(),end:vertex.clone()},
        EulerCommand::KillEdgeVertex{vertex:vertex.clone()},
        EulerCommand::KillEdgeFace{edge,removed_face:face},
    ];
    let kinds=commands.iter().map(EulerCommand::operator_kind).collect::<BTreeSet<_>>();
    let declared=ALL_MODELING_OPERATORS.into_iter().filter(|k|k.is_command()).collect::<BTreeSet<_>>();
    assert_eq!(kinds,declared);
    for command in commands {
        let wire=serde_json::to_value(&command).unwrap();
        assert_eq!(wire["kind"],command.operator_kind().as_str());
    }
}
