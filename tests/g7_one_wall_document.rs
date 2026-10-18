//! G7c source-certified G5-style Boolean document, atomicity and Undo/Redo.
use std::collections::BTreeSet;
use rustsolid::{
    BrepBody,CylindricalBrep,G7BooleanDocument,G7Operand,
    G7PrepareResult,G7Hold,G7HandleError,GeometryTolerance,
    HistoryDirection,Point3,Solid,TopologyName,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn sources()->[BrepBody;2]{
    [
        BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),1.5,5.,tol()).unwrap()),
    ]
}
fn doc(key:&str)->G7BooleanDocument{
    let [a,b]=sources();
    G7BooleanDocument::new(key,a,b,tol()).unwrap()
}
fn commit(doc:&mut G7BooleanDocument){
    let result=doc.prepare().unwrap();
    let G7PrepareResult::Ready(proposal)=result else{
        panic!("G7 source-certified difference Held: {result:?}");
    };
    doc.commit(proposal).unwrap();
}
fn names(doc:&G7BooleanDocument)->BTreeSet<TopologyName>{
    doc.bindings().into_iter().map(|x|x.name).collect()
}
#[test]
fn independently_certified_one_wall_material_commits_atomically_with_84_named_entities(){
    let mut d=doc("one_wall_part");
    assert!(d.model().is_none());
    assert_eq!(d.revision(),0);
    let result=d.prepare().unwrap();
    let G7PrepareResult::Ready(proposal)=result else{panic!("G7 Held")};
    assert!(proposal.candidate.closed_shell_proven);
    assert!(proposal.volume.volume_certified);
    assert!(!proposal.candidate.can_apply);
    assert!(!proposal.volume.can_apply);
    let info=d.commit(proposal).unwrap();
    assert_eq!(info.revision_before,0);
    assert_eq!(info.revision_after,1);
    assert_eq!(info.shell_components,1);
    assert!(info.result_volume>0.);
    assert_eq!(info.topology_names.len(),84);
    assert_eq!(d.model().unwrap().vertices.len(),12);
    assert_eq!(d.model().unwrap().edges.len(),18);
    assert_eq!(d.model().unwrap().coedges.len(),36);
    assert_eq!(d.model().unwrap().faces.len(),8);
    assert_eq!(d.model().unwrap().euler_characteristic(),2);
    assert_eq!(d.undo_depth(),1);
    assert_eq!(d.redo_depth(),0);
    assert_eq!(names(&d).len(),84);
    for binding in d.bindings(){
        assert_eq!(d.resolve_name(&binding.name),Some(binding.clone()));
        assert_eq!(d.resolve_handle(binding.handle).unwrap(),binding);
    }
    assert!(names(&d).iter().any(|x|x.as_str().ends_with("original_xmin_upper")));
    assert!(names(&d).iter().any(|x|x.as_str().ends_with("original_xmin_lower")));
    assert!(matches!(d.prepare().unwrap(),G7PrepareResult::Held(G7Hold::AlreadyCommitted)));
    d.validate().unwrap();
}
#[test]
fn persistent_semantic_face_role_rebinds_without_reviving_a_deleted_handle(){
    let mut d=doc("notch");
    commit(&mut d);
    let old_name=TopologyName("notch/difference/notch/face/inner_cylinder".into());
    let old=d.resolve_name(&old_name).unwrap();
    let original=format!("{:?}",d.model().unwrap());
    let undo=d.undo().unwrap().unwrap();
    assert_eq!(undo.direction,HistoryDirection::Undo);
    assert_eq!(undo.revision_before,1);
    assert_eq!(undo.revision_after,2);
    assert!(!undo.has_result);
    assert!(d.model().is_none());
    assert_eq!(d.undo_depth(),0);
    assert_eq!(d.redo_depth(),1);
    assert!(matches!(d.resolve_handle(old.handle),Err(G7HandleError::StaleHandle)));
    let redo=d.redo().unwrap().unwrap();
    assert_eq!(redo.direction,HistoryDirection::Redo);
    assert!(redo.has_result);
    assert_eq!(redo.revision_after,3);
    let active=d.resolve_name(&old_name).unwrap();
    assert_eq!(old.name,active.name);
    assert_eq!(old.entity,active.entity);
    assert_ne!(old.handle,active.handle);
    assert!(matches!(d.resolve_handle(old.handle),Err(G7HandleError::StaleHandle)));
    assert_eq!(format!("{:?}",d.model().unwrap()),original);
    d.validate().unwrap();
}
#[test]
fn stale_preview_rejected_even_for_a_logical_edit_with_identical_geometry(){
    let mut d=doc("conflict");
    let G7PrepareResult::Ready(old)=d.prepare().unwrap()else{panic!("not ready")};
    let same=d.source(G7Operand::Second).clone();
    assert_eq!(d.replace_source(G7Operand::Second,same).unwrap(),1);
    let state=(d.revision(),d.undo_depth(),d.redo_depth(),d.bindings());
    assert!(d.commit(old).is_err());
    assert_eq!((d.revision(),d.undo_depth(),d.redo_depth(),d.bindings()),state);
    commit(&mut d);
    assert_eq!(d.revision(),2);
    d.validate().unwrap();
}
#[test]
fn corrupted_g7a_edge_or_g7b_source_flux_evidence_rejects_the_whole_commit(){
    let mut d=doc("no_partial");
    let G7PrepareResult::Ready(mut altered)=d.prepare().unwrap()else{panic!("not ready")};
    altered.candidate.component.new_cut_arcs[0].pair=rustsolid::FacePairId(u32::MAX);
    assert!(d.commit(altered).is_err());
    assert!(d.model().is_none());
    assert_eq!(d.revision(),0);
    assert_eq!(d.undo_depth(),0);
    let G7PrepareResult::Ready(mut altered)=d.prepare().unwrap()else{panic!("not ready")};
    altered.volume.source.lower*=0.5;
    assert!(d.commit(altered).is_err());
    assert!(d.model().is_none());
    assert_eq!(d.undo_depth(),0);
    assert_eq!(d.bindings().len(),0);
    let G7PrepareResult::Ready(mut altered)=d.prepare().unwrap()else{panic!("not ready")};
    altered.volume.result_flux.volume+=0.2;
    assert!(d.commit(altered).is_err());
    d.validate().unwrap();
}
#[test]
fn foreign_document_handles_and_translated_geometry_do_not_confuse_names(){
    let mut first=doc("stable_roles");
    let mut other=doc("stable_roles");
    commit(&mut first);commit(&mut other);
    let a=first.bindings()[0].clone();
    assert!(matches!(other.resolve_handle(a.handle),
        Err(G7HandleError::ForeignDocument)));
    assert_eq!(names(&first),names(&other));
    let left=BrepBody::from(Solid::block(p(4.,2.,-6.),5.,4.,8.).unwrap());
    let right=BrepBody::from(CylindricalBrep::upright(
        p(5.2,0.,-2.5),1.8,8.,tol()).unwrap());
    let mut shifted=G7BooleanDocument::new(
        "stable_roles",left,right,tol()).unwrap();
    commit(&mut shifted);
    assert_eq!(names(&first),names(&shifted));
    assert_ne!(first.model().unwrap().mass.volume,shifted.model().unwrap().mass.volume);
    shifted.validate().unwrap();
}
#[test]
fn source_change_invalidates_result_handles_and_branching_clears_redo(){
    let mut d=doc("redo_branch");
    commit(&mut d);
    let old=d.bindings()[0].handle;
    let changed=BrepBody::from(CylindricalBrep::upright(
        p(0.1,0.,0.),1.5,5.,tol()).unwrap());
    d.replace_source(G7Operand::Second,changed).unwrap();
    assert_eq!(d.revision(),2);
    assert!(d.model().is_none());
    assert!(matches!(d.resolve_handle(old),Err(G7HandleError::StaleHandle)));
    d.undo().unwrap();
    assert!(d.model().is_some());
    d.redo().unwrap();
    assert!(d.model().is_none());
    commit(&mut d);
    assert_eq!(d.redo_depth(),0);
    assert!(d.redo().unwrap().is_none());
    d.validate().unwrap();
}
#[test]
fn unsupported_two_wall_and_four_corner_sources_fail_closed(){
    for (width,depth,radius) in [(2.4,6.,1.5),(3.,3.,2.)] {
        let x0=-width/2.;
        let a=BrepBody::from(Solid::block(
            p(x0,1.,-depth/2.),width,3.,depth).unwrap());
        let b=BrepBody::from(CylindricalBrep::upright(
            p(0.,0.,0.),radius,5.,tol()).unwrap());
        let d=G7BooleanDocument::new("reject_other_family",a,b,tol()).unwrap();
        assert!(matches!(d.prepare().unwrap(),G7PrepareResult::Held(G7Hold::UnsupportedGeometry(_))));
        assert_eq!(d.revision(),0);
        assert_eq!(d.bindings().len(),0);
        d.validate().unwrap();
    }
}
#[test]
fn read_only_validation_does_not_allocate_new_live_handles(){
    let mut d=doc("pure");
    commit(&mut d);
    let original=d.bindings();
    for _ in 0..6{
        d.validate().unwrap();
        assert_eq!(original,d.bindings());
        assert_eq!(d.revision(),1);
    }
}
