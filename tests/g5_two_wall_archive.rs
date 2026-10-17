use rustsolid::{
    BrepBody,CylindricalBrep,G5TwoWallDocument,G5TwoWallOperand,
    G5TwoWallPrepareResult,G5TwoWallProjectArchive,
    G5TwoWallHandleError,GeometryTolerance,Point3,Solid,
};
use std::collections::BTreeSet;
use std::fs;
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn doc()->G5TwoWallDocument {
    G5TwoWallDocument::new("g6_parametric",
        BrepBody::from(Solid::block(q(-1.2,1.,-3.),2.4,3.,6.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(
            q(0.,0.,0.),1.5,5.,tol()).unwrap()),
        tol(),
    ).unwrap()
}
fn commit(doc:&mut G5TwoWallDocument){
    let proposal=doc.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(proposal)=proposal else{
        panic!("G6 certified source not ready");
    };
    doc.commit(proposal).unwrap();
}
fn names(doc:&G5TwoWallDocument)->BTreeSet<String>{
    doc.bindings().into_iter().map(|b|b.name.0).collect()
}
#[test]
fn schema2_archive_restores_two_shells_and_115_logical_producer_links(){
    let mut d=doc();
    commit(&mut d);
    let old=d.bindings()[10].clone();
    let volume=d.model().unwrap().mass.volume;
    let original=names(&d);
    let archive=G5TwoWallProjectArchive::from_document(&d).unwrap();
    assert_eq!(archive.schema_version,2);
    assert_eq!(archive.graph.nodes.len(),3);
    assert_eq!(archive.graph.output_bindings.len(),115);
    assert_eq!(archive.undo.len(),1);
    let text=archive.to_json().unwrap();
    assert!(!text.contains("document_token"));
    assert!(!text.contains("generation"));
    assert!(!text.contains("coedges"));
    let parsed=G5TwoWallProjectArchive::from_json(&text).unwrap();
    assert_eq!(parsed,archive);
    let reopened=parsed.reopen().unwrap();
    assert_eq!(reopened.revision(),1);
    assert_eq!(reopened.model().unwrap().shells.len(),2);
    assert_eq!(reopened.model().unwrap().mass.volume,volume);
    assert_eq!(reopened.undo_depth(),1);
    assert_eq!(names(&reopened),original);
    assert!(matches!(reopened.resolve_handle(old.handle),
        Err(G5TwoWallHandleError::ForeignDocument)));
    reopened.validate().unwrap();
}
#[test]
fn undo_redo_survive_file_reopen_and_rebinding_does_not_reactivate_handles(){
    let mut d=doc();
    commit(&mut d);
    let old=d.bindings()[0].clone();
    d.undo().unwrap();
    assert_eq!(d.redo_depth(),1);
    let archive=G5TwoWallProjectArchive::from_document(&d).unwrap();
    assert_eq!(archive.graph.output_bindings.len(),115);
    let mut open=archive.reopen().unwrap();
    assert!(open.model().is_none());
    assert_eq!(open.redo_depth(),1);
    open.redo().unwrap();
    assert_eq!(open.model().unwrap().shells.len(),2);
    let active=open.resolve_name(&old.name).unwrap();
    assert_ne!(active.handle,old.handle);
    assert_eq!(active.name,old.name);
    assert!(matches!(open.resolve_handle(old.handle),
        Err(G5TwoWallHandleError::ForeignDocument)));
    open.validate().unwrap();
}
#[test]
fn source_changes_across_history_reopen_preserve_both_certified_material_states(){
    let mut d=doc();
    commit(&mut d);
    let old=d.model().unwrap().mass.volume;
    d.replace_source(G5TwoWallOperand::Second,
        BrepBody::from(CylindricalBrep::upright(
            q(0.1,0.,0.),1.5,5.,tol()).unwrap())
    ).unwrap();
    assert!(d.model().is_none());
    let archive=G5TwoWallProjectArchive::from_document(&d).unwrap();
    assert_eq!(archive.undo.len(),2);
    let mut opened=archive.reopen().unwrap();
    assert!(opened.model().is_none());
    opened.undo().unwrap();
    assert_eq!(opened.model().unwrap().mass.volume,old);
    opened.redo().unwrap();
    assert!(opened.model().is_none());
    commit(&mut opened);
    let updated=G5TwoWallProjectArchive::from_document(&opened).unwrap();
    assert_eq!(updated.graph.output_bindings.len(),115);
    let next=updated.reopen().unwrap();
    assert_eq!(next.model().unwrap().shells.len(),2);
    next.validate().unwrap();
}
#[test]
fn archive_fails_closed_when_schema_graph_or_independent_mass_tampered(){
    let mut d=doc();
    commit(&mut d);
    let snapshot=G5TwoWallProjectArchive::from_document(&d).unwrap();
    let mut modified=snapshot.clone();
    modified.schema_version=1; // The old G5b four-shell format is not reused.
    assert!(modified.reopen().is_err());
    let mut modified=snapshot.clone();
    modified.graph.output_bindings[0].producing_node="source_a".into();
    assert!(modified.reopen().is_err());
    let mut modified=snapshot.clone();
    modified.current.result.as_mut().unwrap().flux_volume+=1.;
    assert!(modified.reopen().is_err());
    let mut modified=snapshot.clone();
    modified.current.result.as_mut().unwrap().component_volumes[0]+=0.1;
    assert!(modified.reopen().is_err());
    let mut modified=snapshot.clone();
    modified.document_revision=0;
    assert!(modified.reopen().is_err());
    let mut value:serde_json::Value=
        serde_json::from_str(&snapshot.to_json().unwrap()).unwrap();
    value["unrecognized"] = serde_json::json!(true);
    assert!(G5TwoWallProjectArchive::from_json(&value.to_string()).is_err());
    assert!(G5TwoWallProjectArchive::from_json("{}").is_err());
    d.validate().unwrap();
}
#[test]
fn atomic_file_save_reopens_the_actual_certified_document_without_orphaned_handles(){
    let mut d=doc();
    commit(&mut d);
    let first=G5TwoWallProjectArchive::from_document(&d).unwrap();
    let root=std::env::temp_dir().join(format!(
        "rustsolid-g6d-{}-{:?}",std::process::id(),
        std::thread::current().id()
    ));
    fs::create_dir_all(&root).unwrap();
    let path=root.join("two-wall.gefest-g6.json");
    first.save_atomic(&path).unwrap();
    let original=fs::read(&path).unwrap();
    let actual=G5TwoWallProjectArchive::load_file(&path).unwrap();
    assert_eq!(actual,first);
    let mut reopened=actual.reopen().unwrap();
    reopened.undo().unwrap();
    let modified=G5TwoWallProjectArchive::from_document(&reopened).unwrap();
    modified.save_atomic(&path).unwrap();
    assert_ne!(fs::read(&path).unwrap(),original);
    assert_eq!(G5TwoWallProjectArchive::load_file(&path).unwrap(),modified);
    assert!(modified.save_atomic(root.join("no-parent").join("temp.json")).is_err());
    assert_eq!(G5TwoWallProjectArchive::load_file(&path).unwrap(),modified);
    fs::remove_dir_all(root).unwrap();
}
