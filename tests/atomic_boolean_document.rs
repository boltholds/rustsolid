//! G5: bounded, revision-atomic analytic Boolean document and Undo/Redo.
use rustsolid::{
    BrepBody,CylindricalBrep,GeometryTolerance,G5BooleanDocument,
    G5DocumentOperand,G5Hold,G5PrepareResult,G5TopologyHandleError,
    G5TopologyHandle,G5PreparedBoolean,HistoryDirection,
    Point3,Solid,TopologyEntity,TopologyName,
};
use std::collections::{BTreeMap,BTreeSet};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn source_a()->BrepBody{
    BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
}
fn source_b()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}
fn new_doc(key:&str)->G5BooleanDocument{
    G5BooleanDocument::new(key,source_a(),source_b(),tol()).unwrap()
}
fn prepare(doc:&G5BooleanDocument)->G5PreparedBoolean{
    let result=doc.prepare().unwrap();
    let G5PrepareResult::Ready(v)=result else {
        panic!("expected G5 analytic source-certified proposal: {result:?}");
    };
    v
}
fn invariant(doc:&G5BooleanDocument){
    doc.validate().unwrap();
    if let Some(model)=doc.model(){
        model.validate().unwrap();
        assert_eq!(doc.identity_count(),1+24+36+72+20+20+4);
    }else{
        assert_eq!(doc.identity_count(),0);
    }
}
#[test]
fn first_committed_boolean_is_atomic_real_four_shell_model_with_full_identity_map(){
    let mut doc=new_doc("frame_01");
    invariant(&doc);
    let s=doc.snapshot();
    assert_eq!(s.revision,0);
    assert!(!s.has_result);
    let prepared=prepare(&doc);
    assert!(prepared.volume.volume_certified);
    assert!(!prepared.volume.can_apply);
    let commit=doc.commit(prepared).unwrap();
    assert_eq!(commit.revision_before,0);
    assert_eq!(commit.revision_after,1);
    assert_eq!(commit.shell_components,4);
    assert_eq!(commit.bound_entities,177);
    assert_eq!(commit.changed_names.len(),177);
    assert!(commit.committed_volume>0.);
    assert_eq!(doc.revision(),1);
    assert_eq!(doc.undo_depth(),1);
    assert_eq!(doc.redo_depth(),0);
    assert_eq!(doc.model().unwrap().shells.len(),4);
    assert!(doc.check_snapshot(&s).is_err());
    invariant(&doc);
    let names=doc.bindings();
    assert_eq!(names.len(),177);
    let mut unique=BTreeSet::new();
    for binding in &names{
        assert!(unique.insert(binding.name.clone()));
        assert_eq!(doc.resolve_name(&binding.name),Some(binding.clone()));
        assert_eq!(doc.resolve_handle(binding.handle).unwrap(),*binding);
    }
    assert_eq!(unique.len(),177);
    assert!(names.iter().any(|x|x.name.as_str().contains("/face/cylinder_wall")));
    assert!(names.iter().any(|x|x.name.as_str().contains("/edge/bottom_circle/use/")));
    assert!(names.iter().any(|x|x.entity==TopologyEntity::Body));
    let next=doc.prepare().unwrap();
    assert!(matches!(next,G5PrepareResult::Held(G5Hold::ExistingCommittedResult)));
}
#[test]
fn undo_redo_restore_exact_geometry_but_never_reissue_stale_entity_handles(){
    let mut doc=new_doc("four_corners");
    let commit=doc.commit(prepare(&doc)).unwrap();
    assert_eq!(commit.revision_after,1);
    let original=format!("{:?}",doc.model().unwrap());
    let original_volume=doc.model().unwrap().mass.volume;
    let stable_name=TopologyName(
        "four_corners/difference/corner/x_neg/z_neg/edge/bottom_circle".into()
    );
    let first=doc.resolve_name(&stable_name).unwrap();
    let old=doc.snapshot();
    let undone=doc.undo().unwrap().unwrap();
    assert_eq!(undone.direction,HistoryDirection::Undo);
    assert_eq!(undone.revision_before,1);
    assert_eq!(undone.revision_after,2);
    assert!(!undone.has_result);
    assert!(doc.model().is_none());
    assert_eq!(doc.undo_depth(),0);
    assert_eq!(doc.redo_depth(),1);
    assert_eq!(doc.resolve_name(&stable_name),None);
    assert!(matches!(doc.resolve_handle(first.handle),
        Err(G5TopologyHandleError::StaleHandle)));
    assert!(doc.check_snapshot(&old).is_err());
    invariant(&doc);

    let restored=doc.redo().unwrap().unwrap();
    assert_eq!(restored.direction,HistoryDirection::Redo);
    assert_eq!(restored.revision_before,2);
    assert_eq!(restored.revision_after,3);
    assert!(restored.has_result);
    assert_eq!(doc.undo_depth(),1);
    assert_eq!(doc.redo_depth(),0);
    let active=doc.resolve_name(&stable_name).unwrap();
    assert_eq!(first.name,active.name);
    assert_eq!(first.entity,active.entity);
    assert_ne!(first.handle,active.handle);
    assert!(matches!(doc.resolve_handle(first.handle),
        Err(G5TopologyHandleError::StaleHandle)));
    assert_eq!(doc.resolve_handle(active.handle).unwrap(),active);
    assert_eq!(doc.model().unwrap().mass.volume,original_volume);
    assert_eq!(format!("{:?}",doc.model().unwrap()),original);
    invariant(&doc);
}
#[test]
fn source_edit_after_prepare_forces_revision_conflict_even_if_geometry_is_same(){
    let mut doc=new_doc("source_changes");
    let proposal=prepare(&doc);
    let old=doc.snapshot();
    // Replacing a source with identical geometry is still a logical edit
    // and must not permit reuse of a stale G4c proof (ABA revision check).
    let same=doc.source(G5DocumentOperand::Second).clone();
    assert_eq!(doc.replace_source(G5DocumentOperand::Second,same).unwrap(),1);
    assert_eq!(doc.revision(),1);
    assert_eq!(doc.undo_depth(),1);
    assert_eq!(doc.redo_depth(),0);
    let before=format!("{:?}",doc.snapshot());
    assert!(doc.commit(proposal).is_err());
    assert_eq!(format!("{:?}",doc.snapshot()),before);
    assert!(doc.check_snapshot(&old).is_err());
    invariant(&doc);

    let newer=prepare(&doc);
    assert_eq!(newer.expected_revision(),1);
    let report=doc.commit(newer).unwrap();
    assert_eq!(report.revision_after,2);
    assert_eq!(doc.undo_depth(),2);
    invariant(&doc);
}
#[test]
fn every_rejected_or_tampered_proposal_leaves_sources_result_history_and_names_unchanged(){
    let mut doc=new_doc("rollback");
    let mut broken=prepare(&doc);
    broken.volume.components[0].result_flux.result_volume+=0.2;
    let before=format!("{:?}",doc.snapshot());
    let a=format!("{:?}",doc.source(G5DocumentOperand::First));
    let b=format!("{:?}",doc.source(G5DocumentOperand::Second));
    assert!(doc.commit(broken).is_err());
    assert_eq!(format!("{:?}",doc.snapshot()),before);
    assert_eq!(format!("{:?}",doc.source(G5DocumentOperand::First)),a);
    assert_eq!(format!("{:?}",doc.source(G5DocumentOperand::Second)),b);
    assert_eq!(doc.undo_depth(),0);
    invariant(&doc);

    let mut broken=prepare(&doc);
    broken.candidate.corners[0].cut_arcs[0].pair=rustsolid::FacePairId(u32::MAX);
    assert!(doc.commit(broken).is_err());
    assert_eq!(format!("{:?}",doc.snapshot()),before);

    let good=prepare(&doc);
    let report=doc.commit(good).unwrap();
    assert_eq!(report.bound_entities,177);
    let active=format!("{:?}",doc.snapshot());
    let volume=doc.model().unwrap().mass.volume;
    let old_bindings=doc.bindings();
    // This rejected source is invalid even as a standalone BrepBody.
    let mut invalid=source_a();
    if let BrepBody::Polyhedral(ref mut body)=invalid{
        body.faces.clear();
    }
    assert!(doc.replace_source(G5DocumentOperand::First,invalid).is_err());
    assert_eq!(format!("{:?}",doc.snapshot()),active);
    assert_eq!(doc.model().unwrap().mass.volume,volume);
    assert_eq!(doc.bindings(),old_bindings);
    invariant(&doc);
}
#[test]
fn source_replacement_invalidates_committed_topology_and_undo_restores_new_generation(){
    let mut doc=new_doc("source_replace");
    doc.commit(prepare(&doc)).unwrap();
    let old=doc.bindings();
    let (old_name,old_handle)=(old[5].name.clone(),old[5].handle);
    // Geometry is still in the certified four-corner family but differs
    // from the committed sources.
    let changed=BrepBody::from(CylindricalBrep::upright(
        p(0.03,0.,0.),2.,5.,tol()).unwrap());
    let before=doc.revision();
    let after=doc.replace_source(G5DocumentOperand::Second,changed).unwrap();
    assert_eq!(after,before+1);
    assert_eq!(doc.identity_count(),0);
    assert!(doc.model().is_none());
    assert!(matches!(doc.resolve_handle(old_handle),
        Err(G5TopologyHandleError::StaleHandle)));
    invariant(&doc);
    let transition=doc.undo().unwrap().unwrap();
    assert!(transition.has_result);
    let current=doc.resolve_name(&old_name).unwrap();
    assert_ne!(current.handle,old_handle);
    assert_eq!(current.entity,old[5].entity);
    assert!(matches!(doc.resolve_handle(old_handle),
        Err(G5TopologyHandleError::StaleHandle)));
    invariant(&doc);
    // Redo restores modified source (and retires the resurrected result).
    let redo=doc.redo().unwrap().unwrap();
    assert!(!redo.has_result);
    assert!(doc.model().is_none());
    let next=prepare(&doc);
    let r=doc.commit(next).unwrap();
    assert_eq!(r.shell_components,4);
    assert_eq!(doc.redo_depth(),0);
    invariant(&doc);
}
#[test]
fn handle_from_another_document_is_foreign_even_if_stable_names_match(){
    let mut a=new_doc("matching_key");
    let mut b=new_doc("matching_key");
    let left=prepare(&a);
    let right=prepare(&b);
    let sa=format!("{:?}",a.snapshot());
    let sb=format!("{:?}",b.snapshot());
    assert!(a.commit(right).is_err());
    assert!(b.commit(left.clone()).is_err());
    assert_eq!(format!("{:?}",a.snapshot()),sa);
    assert_eq!(format!("{:?}",b.snapshot()),sb);
    a.commit(left).unwrap();
    b.commit(prepare(&b)).unwrap();
    let first=a.bindings();
    let second=b.bindings();
    assert_eq!(first.len(),second.len());
    let unique_a:BTreeMap<_,_>=first.iter()
        .map(|x|(x.name.clone(),x.handle)).collect();
    let unique_b:BTreeMap<_,_>=second.iter()
        .map(|x|(x.name.clone(),x.handle)).collect();
    assert_eq!(unique_a.keys().collect::<Vec<_>>(),unique_b.keys().collect::<Vec<_>>());
    assert_ne!(unique_a.values().collect::<Vec<_>>(),unique_b.values().collect::<Vec<_>>());
    assert!(matches!(b.resolve_handle(first[0].handle),
        Err(G5TopologyHandleError::ForeignDocument)));
    invariant(&a);invariant(&b);
}
#[test]
fn undo_then_new_commit_clears_redo_without_reusing_old_result_handles(){
    let mut doc=new_doc("redo_branch");
    doc.commit(prepare(&doc)).unwrap();
    let old_handle=doc.bindings()[0].handle;
    doc.undo().unwrap();
    assert_eq!(doc.redo_depth(),1);
    assert_eq!(doc.undo_depth(),0);
    let fresh=prepare(&doc);
    let report=doc.commit(fresh).unwrap();
    assert_eq!(report.revision_before,2);
    assert_eq!(report.revision_after,3);
    assert_eq!(doc.redo_depth(),0);
    assert!(doc.redo().unwrap().is_none());
    assert!(matches!(doc.resolve_handle(old_handle),
        Err(G5TopologyHandleError::StaleHandle)));
    invariant(&doc);
}
#[test]
fn unsupported_sources_fail_closed_during_prepare_without_document_changes(){
    let far=BrepBody::from(Solid::block(p(20.,20.,20.),2.,2.,2.).unwrap());
    let doc=G5BooleanDocument::new("not_intersecting",far,source_b(),tol()).unwrap();
    let before=format!("{:?}",doc.snapshot());
    assert!(matches!(doc.prepare().unwrap(),
        G5PrepareResult::Held(G5Hold::UnsupportedCandidate(_))));
    assert_eq!(format!("{:?}",doc.snapshot()),before);
    invariant(&doc);
}


#[test]
fn read_only_validation_and_snapshot_do_not_change_active_handle_bindings(){
    let mut doc=new_doc("pure_validation");
    doc.commit(prepare(&doc)).unwrap();
    let original_snapshot=format!("{:?}",doc.snapshot());
    let original_bindings=doc.bindings();
    let original_revision=doc.revision();
    for _ in 0..8 {
        doc.validate().unwrap();
        let _=doc.snapshot();
        assert_eq!(doc.bindings(),original_bindings);
        assert_eq!(doc.revision(),original_revision);
    }
    assert_eq!(format!("{:?}",doc.snapshot()),original_snapshot);
}
#[test]
fn stable_semantic_names_are_independent_of_specific_corner_dimensions(){
    let mut symmetric=new_doc("part_family");
    symmetric.commit(prepare(&symmetric)).unwrap();
    let a=BrepBody::from(Solid::block(p(2.,6.,-8.),4.5,2.5,3.5).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        p(4.,1.,-6.2),2.55,10.,tol()).unwrap());
    let mut asymmetric=G5BooleanDocument::new("part_family",a,b,tol()).unwrap();
    asymmetric.commit(prepare(&asymmetric)).unwrap();
    let semantic_symmetric=symmetric.bindings().iter()
        .map(|x|x.name.clone()).collect::<BTreeSet<_>>();
    let semantic_asymmetric=asymmetric.bindings().iter()
        .map(|x|x.name.clone()).collect::<BTreeSet<_>>();
    assert_eq!(semantic_symmetric,semantic_asymmetric);
    assert_ne!(symmetric.model().unwrap().mass.volume,
        asymmetric.model().unwrap().mass.volume);
    // Names are stable logical roles, but process-local handles are not
    // portable across two distinct documents.
    assert!(matches!(symmetric.resolve_handle(asymmetric.bindings()[0].handle),
        Err(G5TopologyHandleError::ForeignDocument)));
    invariant(&symmetric);invariant(&asymmetric);
}
