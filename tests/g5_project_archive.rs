//! G5b: source-rebuilt Gefest project archive, DAG, names, Undo/Redo.
use rustsolid::{
    BrepBody,CylindricalBrep,G5BooleanDocument,G5DocumentOperand,
    G5PrepareResult,G5ProjectArchive,G5TopologyHandleError,
    GeometryTolerance,Point3,Solid,TopologyName,
};
use std::collections::BTreeSet;
use std::fs;

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tol()->GeometryTolerance {GeometryTolerance::default()}
fn make()->G5BooleanDocument{
    G5BooleanDocument::new("gef_part_a",
        BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap()),
        tol()).unwrap()
}
fn commit(doc:&mut G5BooleanDocument){
    let G5PrepareResult::Ready(stage)=doc.prepare().unwrap()
        else {panic!("source G4c should have been ready")};
    doc.commit(stage).unwrap();
}
fn snapshot(doc:&G5BooleanDocument)->String{
    format!("{:?}",doc.snapshot())
}
fn graph_names(doc:&G5BooleanDocument)->BTreeSet<TopologyName>{
    doc.bindings().into_iter().map(|b|b.name).collect()
}

#[test]
fn reopening_real_project_rebuilds_brep_graph_and_semantic_names(){
    let mut doc=make();
    commit(&mut doc);
    let source_volume=doc.model().unwrap().mass.volume;
    let old=doc.bindings()[0].handle;
    let original_names=graph_names(&doc);
    assert_eq!(original_names.len(),177);
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    assert_eq!(archive.schema_version,1);
    assert_eq!(archive.document_revision,1);
    assert_eq!(archive.graph.nodes.len(),3);
    assert_eq!(archive.graph.output_bindings.len(),177);
    assert_eq!(archive.undo.len(),1);
    let encoded=archive.to_json().unwrap();
    assert!(encoded.contains("source_a"));
    assert!(encoded.contains("source_b"));
    assert!(encoded.contains("difference"));
    assert!(encoded.contains("gef_part_a/difference/corner"));
    assert!(!encoded.contains("document_token"));
    assert!(!encoded.contains("generation"));
    assert!(!encoded.contains("coedges"));
    let read=G5ProjectArchive::from_json(&encoded).unwrap();
    assert_eq!(read,archive);
    let mut reopened=read.reopen().unwrap();
    reopened.validate().unwrap();
    assert_eq!(reopened.revision(),doc.revision());
    assert_eq!(reopened.undo_depth(),1);
    assert_eq!(reopened.redo_depth(),0);
    assert_eq!(graph_names(&reopened),original_names);
    assert_eq!(reopened.model().unwrap().mass.volume,source_volume);
    assert!(matches!(reopened.resolve_handle(old),
        Err(G5TopologyHandleError::ForeignDocument)));
    // Restored history is live, no synthetic "loaded, no undo" state.
    reopened.undo().unwrap().unwrap();
    assert!(reopened.model().is_none());
    reopened.redo().unwrap().unwrap();
    assert_eq!(graph_names(&reopened),original_names);
    reopened.validate().unwrap();
}
#[test]
fn redo_stack_survives_reopen_and_old_handles_do_not_become_valid(){
    let mut doc=make();
    commit(&mut doc);
    let old=doc.bindings()[8].clone();
    doc.undo().unwrap();
    assert!(doc.model().is_none());
    assert_eq!(doc.redo_depth(),1);
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    assert_eq!(archive.graph.output_bindings.len(),177,
        "feature output names must survive while result is undone");
    let mut opened=G5ProjectArchive::from_json(&archive.to_json().unwrap())
        .unwrap().reopen().unwrap();
    assert_eq!(opened.redo_depth(),1);
    assert_eq!(opened.undo_depth(),0);
    assert!(opened.model().is_none());
    opened.redo().unwrap().unwrap();
    let current=opened.resolve_name(&old.name).unwrap();
    assert_eq!(current.name,old.name);
    assert_ne!(current.handle,old.handle);
    assert!(matches!(opened.resolve_handle(old.handle),
        Err(G5TopologyHandleError::ForeignDocument)));
    assert!(opened.model().is_some());
    opened.validate().unwrap();
}
#[test]
fn source_replacement_with_history_roundtrips_and_restores_both_versions(){
    let mut doc=make();
    commit(&mut doc);
    let old_volume=doc.model().unwrap().mass.volume;
    let source=BrepBody::from(CylindricalBrep::upright(
        p(0.03,0.,0.),2.,5.,tol()).unwrap());
    doc.replace_source(G5DocumentOperand::Second,source).unwrap();
    assert!(doc.model().is_none());
    assert_eq!(doc.revision(),2);
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    let mut opened=G5ProjectArchive::from_json(&archive.to_json().unwrap())
        .unwrap().reopen().unwrap();
    assert_eq!(opened.revision(),2);
    assert_eq!(opened.undo_depth(),2);
    assert!(opened.model().is_none());
    opened.undo().unwrap();
    assert!((opened.model().unwrap().mass.volume-old_volume).abs()<1e-12);
    opened.redo().unwrap();
    assert!(opened.model().is_none());
    commit(&mut opened);
    assert!(opened.model().unwrap().mass.volume>0.);
    assert_eq!(opened.redo_depth(),0);
    let reread=G5ProjectArchive::from_document(&opened).unwrap();
    let final_doc=reread.reopen().unwrap();
    assert_eq!(graph_names(&final_doc).len(),177);
}
#[test]
fn corrupt_graph_stamp_schema_history_and_unknown_fields_fail_closed(){
    let mut doc=make();
    commit(&mut doc);
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    let before=snapshot(&doc);
    let mut bad=archive.clone();
    bad.schema_version=45;
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.graph.nodes.pop();
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.graph.output_bindings[0].producing_node="source_a".into();
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.current.result.as_mut().unwrap().mass_volume+=0.1;
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.current.sources.swap(0,1);
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.undo[0].result=archive.current.result.clone();
    assert!(bad.reopen().is_err());
    let mut bad=archive.clone();
    bad.document_revision=0;
    assert!(bad.reopen().is_err());
    let mut value:serde_json::Value=
        serde_json::from_str(&archive.to_json().unwrap()).unwrap();
    value["malicious_unknown_property"]=serde_json::json!(true);
    assert!(G5ProjectArchive::from_json(&value.to_string()).is_err());
    assert_eq!(snapshot(&doc),before);
}
#[test]
fn atomic_file_write_and_load_rebuild_equivalent_document_and_history(){
    let mut doc=make();
    commit(&mut doc);
    doc.undo().unwrap();
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    let root=std::env::temp_dir().join(format!(
        "rustsolid-g5b-{}-{:?}",std::process::id(),
        std::thread::current().id()));
    fs::create_dir_all(&root).unwrap();
    let path=root.join("assembly.gefest-g5.json");
    archive.save_atomic(&path).unwrap();
    let first=fs::read(&path).unwrap();
    assert!(first.len()>1000);
    let loaded=G5ProjectArchive::load_file(&path).unwrap();
    assert_eq!(loaded,archive);
    let reopened=loaded.reopen().unwrap();
    assert_eq!(reopened.redo_depth(),1);
    assert!(reopened.model().is_none());
    let mut reopened=reopened;
    reopened.redo().unwrap();
    let updated=G5ProjectArchive::from_document(&reopened).unwrap();
    updated.save_atomic(&path).unwrap();
    assert_ne!(fs::read(&path).unwrap(),first);
    assert_eq!(G5ProjectArchive::load_file(&path).unwrap(),updated);
    fs::remove_dir_all(root).unwrap();
}
#[test]
fn invalid_file_destination_cannot_overwrite_existing_project(){
    let mut doc=make();
    commit(&mut doc);
    let archive=G5ProjectArchive::from_document(&doc).unwrap();
    let root=std::env::temp_dir().join(format!(
        "rustsolid-g5b-rollback-{}-{:?}",std::process::id(),
        std::thread::current().id()));
    fs::create_dir_all(&root).unwrap();
    let path=root.join("project.json");
    archive.save_atomic(&path).unwrap();
    let stable=fs::read(&path).unwrap();
    assert!(archive.save_atomic(&root).is_err(),
        "renaming file over existing non-empty folder must fail");
    assert_eq!(fs::read(&path).unwrap(),stable);
    assert!(archive.save_atomic(root.join("no-parent").join("project.json")).is_err());
    assert_eq!(fs::read(&path).unwrap(),stable);
    assert_eq!(G5ProjectArchive::load_file(&path).unwrap(),archive);
    fs::remove_dir_all(root).unwrap();
}
#[test]
fn unsupported_source_and_mass_tamper_cannot_be_persisted_as_verified_result(){
    let unsupported=G5BooleanDocument::new("invalid_source",
        BrepBody::from(Solid::block(p(10.,10.,10.),2.,2.,2.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap()),
        tol()).unwrap();
    let archive=G5ProjectArchive::from_document(&unsupported).unwrap();
    assert!(archive.graph.output_bindings.is_empty());
    let restored=archive.reopen().unwrap();
    assert!(restored.model().is_none());
    assert_eq!(restored.undo_depth(),0);
    let mut doc=make();
    commit(&mut doc);
    let mut damaged=G5ProjectArchive::from_document(&doc).unwrap();
    damaged.current.result.as_mut().unwrap().source_upper*=1.5;
    assert!(damaged.reopen().is_err());
}
