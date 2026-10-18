//! G7d project archive v3: rebuild source geometry and historical proofs.
use std::collections::BTreeSet;
use std::fs;
use rustsolid::{
    BrepBody,CylindricalBrep,G7BooleanDocument,G7Operand,
    G7PrepareResult,G7ProjectArchive,G7HandleError,
    GeometryTolerance,Point3,Solid,TopologyName,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn create()->G7BooleanDocument{
    G7BooleanDocument::new(
        "edge_notch",
        BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),1.5,5.,tol()).unwrap()),
        tol(),
    ).unwrap()
}
fn commit(doc:&mut G7BooleanDocument){
    let G7PrepareResult::Ready(prepared)=doc.prepare().unwrap()
        else{panic!("G7 source geometry expected to be supported")};
    doc.commit(prepared).unwrap();
}
fn names(doc:&G7BooleanDocument)->BTreeSet<TopologyName>{
    doc.bindings().into_iter().map(|x|x.name).collect()
}
#[test]
fn archive_v3_roundtrip_recreates_single_shell_and_84_semantic_operation_bindings(){
    let mut d=create();
    commit(&mut d);
    let original_names=names(&d);
    let original_volume=d.model().unwrap().mass.volume;
    let old=d.bindings()[15].handle;
    let archive=G7ProjectArchive::from_document(&d).unwrap();
    assert_eq!(archive.schema_version,3);
    assert_eq!(archive.document_revision,1);
    assert_eq!(archive.graph.nodes.len(),3);
    assert_eq!(archive.graph.output_bindings.len(),84);
    assert_eq!(archive.undo.len(),1);
    assert_eq!(archive.redo.len(),0);
    assert_eq!(archive.current.result.as_ref().unwrap().source_lower
        <archive.current.result.as_ref().unwrap().flux_volume,true);
    let encoded=archive.to_json().unwrap();
    assert!(encoded.contains("source_a"));
    assert!(encoded.contains("difference"));
    assert!(encoded.contains("edge_notch/difference/notch"));
    assert!(!encoded.contains("document_token"));
    assert!(!encoded.contains("generation"));
    assert!(!encoded.contains("coedges"));
    let read=G7ProjectArchive::from_json(&encoded).unwrap();
    assert_eq!(read,archive);
    let reopened=read.reopen().unwrap();
    assert_eq!(reopened.revision(),1);
    assert_eq!(reopened.model().unwrap().shells.len(),1);
    assert_eq!(reopened.model().unwrap().mass.volume,original_volume);
    assert_eq!(reopened.undo_depth(),1);
    assert_eq!(names(&reopened),original_names);
    assert!(matches!(reopened.resolve_handle(old),Err(G7HandleError::ForeignDocument)));
    reopened.validate().unwrap();
}
#[test]
fn archived_redo_state_survives_reopening_without_reactivating_old_topology_handles(){
    let mut d=create();
    commit(&mut d);
    let old=d.bindings()[0].clone();
    d.undo().unwrap();
    assert!(d.model().is_none());
    assert_eq!(d.redo_depth(),1);
    let archive=G7ProjectArchive::from_document(&d).unwrap();
    assert_eq!(archive.graph.output_bindings.len(),84,
        "semantic producer links must survive an undone result");
    let mut restored=G7ProjectArchive::from_json(&archive.to_json().unwrap())
        .unwrap().reopen().unwrap();
    assert_eq!(restored.redo_depth(),1);
    assert_eq!(restored.undo_depth(),0);
    restored.redo().unwrap();
    assert_eq!(restored.model().unwrap().faces.len(),8);
    let fresh=restored.resolve_name(&old.name).unwrap();
    assert_ne!(fresh.handle,old.handle);
    assert_eq!(fresh.name,old.name);
    assert!(matches!(restored.resolve_handle(old.handle),Err(G7HandleError::ForeignDocument)));
    restored.validate().unwrap();
}
#[test]
fn source_edit_history_can_reopen_old_result_and_commit_new_source(){
    let mut d=create();
    commit(&mut d);
    let volume=d.model().unwrap().mass.volume;
    let updated=BrepBody::from(CylindricalBrep::upright(
        p(0.1,0.,0.),1.5,5.,tol()).unwrap());
    d.replace_source(G7Operand::Second,updated).unwrap();
    assert!(d.model().is_none());
    assert_eq!(d.revision(),2);
    let snapshot=G7ProjectArchive::from_document(&d).unwrap();
    assert_eq!(snapshot.undo.len(),2);
    let mut restored=snapshot.reopen().unwrap();
    restored.undo().unwrap();
    assert_eq!(restored.model().unwrap().mass.volume,volume);
    restored.redo().unwrap();
    assert!(restored.model().is_none());
    commit(&mut restored);
    assert_eq!(restored.redo_depth(),0);
    let second=G7ProjectArchive::from_document(&restored).unwrap();
    assert_eq!(second.graph.output_bindings.len(),84);
    second.reopen().unwrap().validate().unwrap();
}
#[test]
fn foreign_or_tampered_schema_graph_mass_and_unsupported_source_are_rejected(){
    let mut d=create();
    commit(&mut d);
    let good=G7ProjectArchive::from_document(&d).unwrap();
    let mut damaged=good.clone();
    damaged.schema_version=1;
    assert!(damaged.reopen().is_err());
    damaged.schema_version=2;
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.graph.nodes.pop();
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.graph.output_bindings[0].producing_node="source_a".into();
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.current.result.as_mut().unwrap().source_upper+=0.1;
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.current.result.as_mut().unwrap().true_angular_span+=0.01;
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.document_revision=0;
    assert!(damaged.reopen().is_err());
    let mut damaged=good.clone();
    damaged.current.sources.swap(0,1);
    assert!(damaged.reopen().is_err());
    let mut value:serde_json::Value=
        serde_json::from_str(&good.to_json().unwrap()).unwrap();
    value["unknown_field"]=serde_json::json!(true);
    assert!(G7ProjectArchive::from_json(&value.to_string()).is_err());
    d.validate().unwrap();
}
#[test]
fn atomic_disk_write_and_fallback_failure_preserve_previous_project(){
    let mut d=create();
    commit(&mut d);
    let archive=G7ProjectArchive::from_document(&d).unwrap();
    let root=std::env::temp_dir().join(format!(
        "rustsolid-g7d-{}-{:?}",std::process::id(),
        std::thread::current().id(),
    ));
    fs::create_dir_all(&root).unwrap();
    let file=root.join("onewall.gefest.json");
    archive.save_atomic(&file).unwrap();
    let original=fs::read(&file).unwrap();
    assert!(original.len()>1000);
    assert_eq!(G7ProjectArchive::load_file(&file).unwrap(),archive);
    let mut fresh=G7ProjectArchive::load_file(&file).unwrap().reopen().unwrap();
    fresh.undo().unwrap();
    let updated=G7ProjectArchive::from_document(&fresh).unwrap();
    updated.save_atomic(&file).unwrap();
    assert_ne!(fs::read(&file).unwrap(),original);
    assert_eq!(G7ProjectArchive::load_file(&file).unwrap(),updated);
    assert!(updated.save_atomic(root.join("missing").join("project.json")).is_err());
    assert_eq!(G7ProjectArchive::load_file(&file).unwrap(),updated);
    fs::remove_dir_all(root).unwrap();
}
#[test]
fn source_solid_revision_without_a_durable_event_log_cannot_be_persisted(){
    let mut body=Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap();
    body.revision+=1;
    let cutter=CylindricalBrep::upright(p(0.,0.,0.),1.5,5.,tol()).unwrap();
    let doc=G7BooleanDocument::new("non_durable",
        BrepBody::Polyhedral(body),BrepBody::Cylinder(cutter),tol()).unwrap();
    assert!(G7ProjectArchive::from_document(&doc).is_err());
    assert_eq!(doc.revision(),0);
}
