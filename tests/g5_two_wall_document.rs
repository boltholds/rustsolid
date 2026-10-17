use rustsolid::{
    BrepBody,CylindricalBrep,G5TwoWallDocument,
    G5TwoWallPrepareResult,G5TwoWallHold,G5TwoWallOperand,
    G5TwoWallHandleError,GeometryTolerance,
    HistoryDirection,Point3,Solid,TopologyName,
};
use std::collections::BTreeSet;
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn source()->(BrepBody,BrepBody){
    (
        BrepBody::from(Solid::block(q(-1.2,1.,-3.),2.4,3.,6.).unwrap()),
        BrepBody::from(CylindricalBrep::upright(q(0.,0.,0.),1.5,5.,tol()).unwrap()),
    )
}
fn doc(key:&str)->G5TwoWallDocument{
    let (a,b)=source();
    G5TwoWallDocument::new(key,a,b,tol()).unwrap()
}
fn commit(d:&mut G5TwoWallDocument){
    let r=d.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(p)=r else{panic!("G5 two-wall Held: {r:?}")};
    d.commit(p).unwrap();
}
#[test]
fn certified_two_wall_geometry_commits_atomically_with_115_stable_names(){
    let mut d=doc("prismatic_1");
    let before=d.revision();
    assert!(d.model().is_none());
    let p=d.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(p)=p else{panic!("not ready")};
    assert!(p.volume.volume_certified&&!p.candidate.can_apply);
    let result=d.commit(p).unwrap();
    assert_eq!(result.revision_before,before);
    assert_eq!(result.revision_after,1);
    assert_eq!(result.shells,2);
    assert_eq!(result.identities.len(),115);
    assert_eq!(d.bindings().len(),115);
    assert_eq!(d.model().unwrap().shells.len(),2);
    assert_eq!(d.model().unwrap().euler_characteristic(),4);
    let mut unique=BTreeSet::new();
    for b in d.bindings(){
        assert!(unique.insert(b.name.clone()));
        assert_eq!(d.resolve_name(&b.name),Some(b.clone()));
        assert_eq!(d.resolve_handle(b.handle).unwrap(),b);
    }
    assert!(unique.iter().any(|s|s.as_str().ends_with("inner_cylinder_wall")));
    assert_eq!(unique.len(),115);
    assert_eq!(d.undo_depth(),1);
    d.validate().unwrap();
    assert!(matches!(d.prepare().unwrap(),
        G5TwoWallPrepareResult::Held(G5TwoWallHold::AlreadyCommitted)));
}
#[test]
fn undo_redo_restore_analytic_result_with_fresh_handles_and_same_names(){
    let mut d=doc("part_A");
    commit(&mut d);
    let shape=format!("{:?}",d.model().unwrap());
    let reference=TopologyName(
        "part_A/difference/strip/z_pos/edge/bottom_circle".into());
    let old=d.resolve_name(&reference).unwrap();
    let undone=d.undo().unwrap().unwrap();
    assert_eq!(undone.direction,HistoryDirection::Undo);
    assert_eq!(undone.revision_after,2);
    assert!(!undone.has_result);
    assert!(d.model().is_none());
    assert!(matches!(d.resolve_handle(old.handle),
        Err(G5TwoWallHandleError::StaleHandle)));
    let redone=d.redo().unwrap().unwrap();
    assert_eq!(redone.revision_after,3);
    assert!(redone.has_result);
    assert_eq!(format!("{:?}",d.model().unwrap()),shape);
    let active=d.resolve_name(&reference).unwrap();
    assert_eq!(old.name,active.name);
    assert_eq!(old.entity,active.entity);
    assert_ne!(old.handle,active.handle);
    assert!(matches!(d.resolve_handle(old.handle),
        Err(G5TwoWallHandleError::StaleHandle)));
    d.validate().unwrap();
}
#[test]
fn stale_preview_and_invalid_revisions_never_mutate_document(){
    let mut d=doc("source_edit");
    let r=d.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(proposal)=r else{panic!("not ready")};
    let same=d.source(G5TwoWallOperand::Second).clone();
    d.replace_source(G5TwoWallOperand::Second,same).unwrap();
    let before=d.revision();
    assert!(d.commit(proposal).is_err());
    assert_eq!(d.revision(),before);
    assert!(d.model().is_none());
    assert_eq!(d.undo_depth(),1);
    d.validate().unwrap();
    commit(&mut d);
    assert_eq!(d.revision(),2);
}
#[test]
fn malformed_candidate_and_volume_proof_fail_without_partial_changes(){
    let mut d=doc("reject");
    let first=d.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(mut altered)=first else{panic!("not ready")};
    altered.candidate.components[0].analytic_volume+=1.;
    assert!(d.commit(altered).is_err());
    assert!(d.model().is_none());
    assert_eq!(d.undo_depth(),0);
    let next=d.prepare().unwrap();
    let G5TwoWallPrepareResult::Ready(mut altered)=next else{panic!("not ready")};
    altered.volume.components[0].boundary.result_volume+=0.4;
    assert!(d.commit(altered).is_err());
    assert_eq!(d.revision(),0);
    assert_eq!(d.bindings().len(),0);
    d.validate().unwrap();
}
#[test]
fn foreign_document_handles_and_semantic_names_rebind_safely(){
    let mut a=doc("same_part");
    let mut b=doc("same_part");
    commit(&mut a);commit(&mut b);
    let first=a.bindings();
    let second=b.bindings();
    assert_eq!(first.iter().map(|x|&x.name).collect::<Vec<_>>(),
        second.iter().map(|x|&x.name).collect::<Vec<_>>());
    assert_ne!(first[0].handle,second[0].handle);
    assert!(matches!(b.resolve_handle(first[0].handle),
        Err(G5TwoWallHandleError::ForeignDocument)));
    let third=BrepBody::from(Solid::block(q(1.,2.,-5.),3.,4.,9.).unwrap());
    let fourth=BrepBody::from(CylindricalBrep::upright(
        q(2.4,0.,-1.),1.8,8.,tol()).unwrap());
    let mut asym=G5TwoWallDocument::new("same_part",third,fourth,tol()).unwrap();
    commit(&mut asym);
    assert_eq!(first.iter().map(|x|&x.name).collect::<Vec<_>>(),
        asym.bindings().iter().map(|x|&x.name).collect::<Vec<_>>());
    a.validate().unwrap();b.validate().unwrap();asym.validate().unwrap();
}
#[test]
fn source_replace_and_redo_branch_are_atomic(){
    let mut d=doc("edit");
    commit(&mut d);
    let old=d.bindings()[0].handle;
    let replacement=BrepBody::from(CylindricalBrep::upright(
        q(0.1,0.,0.),1.5,5.,tol()).unwrap());
    d.replace_source(G5TwoWallOperand::Second,replacement).unwrap();
    assert!(d.model().is_none());
    assert_eq!(d.revision(),2);
    d.undo().unwrap().unwrap();
    assert!(d.model().is_some());
    assert!(matches!(d.resolve_handle(old),Err(G5TwoWallHandleError::StaleHandle)));
    d.redo().unwrap().unwrap();
    assert!(d.model().is_none());
    commit(&mut d);
    assert_eq!(d.redo_depth(),0);
    assert!(d.redo().unwrap().is_none());
    d.validate().unwrap();
}
#[test]
fn unsupported_originals_remain_held_and_do_not_change_current_state(){
    let a=BrepBody::from(Solid::block(q(-1.5,1.,-1.5),3.,3.,3.).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        q(0.,0.,0.),2.,5.,tol()).unwrap());
    let d=G5TwoWallDocument::new("four_corner_not_two",a,b,tol()).unwrap();
    assert!(matches!(d.prepare().unwrap(),
        G5TwoWallPrepareResult::Held(G5TwoWallHold::UnsupportedGeometry(_))));
    assert!(d.model().is_none());
    assert_eq!(d.revision(),0);
    d.validate().unwrap();
}
