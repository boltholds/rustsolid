//! Independent first-party genus-N structural B-rep / analytic mass witnesses.
//! All geometry is constructed from explicit circular curves and UV loops,
//! not triangulated Boolean approximations.
use rustsolid::{
    BrepStructuralEditor, BrepBody, BrepArena, HoleParameters, HoleEdgeRole,
    HoleTrimRole, StructuralHoleHandle, StructuralHandleError,
    StructuralEdgeRef, StructuralTrimRef, StructuralFaceRef,
    FaceId, EdgeId, CoedgeId, Point3, Solid, LoopRole, HistoryDirection,
    Curve3, Surface3, ThroughHoleBrep,
};
use std::f64::consts::{PI,TAU};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn h(cx:f64,cz:f64,r:f64)->HoleParameters {
    HoleParameters{center_x:cx,center_z:cz,radius:r}
}
fn editor()->BrepStructuralEditor{BrepStructuralEditor::new(block()).unwrap()}
fn eq(a:f64,b:f64){assert!((a-b).abs()<1e-9,"{a} vs {b}");}
fn prove(e:&BrepStructuralEditor,holes:&[HoleParameters]){
    e.validate().unwrap();
    let model=e.model();
    let summary=model.summary();
    let count=holes.len();
    assert_eq!(summary.vertices,8+2*count);
    assert_eq!(summary.edges,12+3*count);
    assert_eq!(summary.coedges,24+6*count);
    assert_eq!(summary.faces,6+count);
    assert_eq!(summary.loops,6+3*count);
    assert_eq!(summary.shells,1);
    assert_eq!(summary.genus,count);
    assert_eq!(summary.euler_characteristic,2-2*count as isize);
    assert_eq!(summary.circular_edges,2*count);
    assert_eq!(summary.cylindrical_faces,count);
    assert_eq!(summary.seam_edges,count);
    assert_eq!(model.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),2*count);
    assert_eq!(model.faces[0].loops.len(),count+1);
    assert_eq!(model.faces[1].loops.len(),count+1);
    let total=holes.iter().map(|h|h.radius*h.radius).sum::<f64>();
    let linear=holes.iter().map(|h|h.radius).sum::<f64>();
    let removed=4.*PI*total;
    let mass=144.-removed;
    eq(model.mass.volume,mass);
    eq(model.mass.surface_area,168.-2.*PI*total+8.*PI*linear);
    let cx=holes.iter().map(|h|4.*PI*h.radius*h.radius*h.center_x).sum::<f64>();
    let cz=holes.iter().map(|h|4.*PI*h.radius*h.radius*h.center_z).sum::<f64>();
    eq(model.mass.centroid.x,-cx/mass);
    eq(model.mass.centroid.z,-cz/mass);
    eq(model.mass.centroid.y,1.);
}

#[test]
fn one_transaction_creates_two_real_cylindrical_faces_and_genus_two(){
    let mut e=editor();
    let before=e.snapshot();
    let mut tx=e.begin_transaction();
    let first=tx.add_through_hole(h(-1.5,-0.5,0.7)).unwrap();
    let second=tx.add_through_hole(h(1.4,0.7,0.9)).unwrap();
    assert_eq!(tx.staged_len(),2);
    let report=tx.commit().unwrap();
    assert_eq!((report.revision_before,report.revision_after),(0,1));
    assert_eq!(report.created.len(),2);
    assert_eq!(report.command_count,2);
    assert!(report.retired.is_empty());
    assert_eq!(e.undo_depth(),1);
    assert_eq!(e.redo_depth(),0);
    prove(&e,&[h(-1.5,-0.5,0.7),h(1.4,0.7,0.9)]);
    assert!(e.check_snapshot(&before).is_err());
    assert_eq!(e.resolve(first).unwrap().topology.side_face,FaceId(6));
    assert_eq!(e.resolve(second).unwrap().topology.side_face,FaceId(7));
    for f in [FaceId(6),FaceId(7)] {
        assert!(matches!(e.model().face_surface(f),Some(Surface3::Cylinder(_))));
    }
}

#[test]
fn one_transaction_unifies_structure_and_carrier_reparameterization(){
    let mut e=editor();
    let initial=e.snapshot();
    let mut tx=e.begin_transaction();
    let first=tx.add_through_hole(h(-1.4,0.0,0.75)).unwrap();
    let second=tx.add_through_hole(h(1.4,0.0,0.75)).unwrap();
    tx.reparameterize_edge(StructuralEdgeRef::Hole{
        hole:first,role:HoleEdgeRole::BottomCircle},1.,0.5).unwrap();
    tx.reparameterize_trim(StructuralTrimRef::Hole{
        hole:second,role:HoleTrimRole::TopCap},1.,TAU).unwrap();
    tx.rebase_face_uv(StructuralFaceRef::Base(FaceId(0)),0.1,-0.2).unwrap();
    tx.rebase_face_uv(StructuralFaceRef::HoleSide(second),0.,0.25).unwrap();
    assert_eq!(tx.staged_len(),6);
    let committed=tx.commit().unwrap();
    assert_eq!(committed.command_count,6);
    assert_eq!(e.revision(),1);
    assert_eq!(e.undo_depth(),1);
    prove(&e,&[h(-1.4,0.,0.75),h(1.4,0.,0.75)]);
    assert_eq!(e.model().summary().genus,2);
    let top=e.resolve(second).unwrap().topology.top_edge;
    assert!(matches!(e.model().edge_curve(top),Some(Curve3::Circle(_))));
    assert!(e.check_snapshot(&initial).is_err());
}

#[test]
fn removing_first_hole_compacts_surviving_topology_without_aba(){
    let mut e=editor();
    let (first,second)={
        let mut tx=e.begin_transaction();
        let a=tx.add_through_hole(h(-1.5,0.,0.6)).unwrap();
        let b=tx.add_through_hole(h(1.5,0.,0.7)).unwrap();
        tx.commit().unwrap();
        (a,b)
    };
    let second_name_before=e.resolve(second).unwrap().topology;
    let mut tx=e.begin_transaction();
    tx.remove_through_hole(first).unwrap();
    let report=tx.commit().unwrap();
    assert_eq!(report.retired.len(),1);
    assert_eq!(report.retired[0].handle,first);
    assert_eq!(report.remapped.len(),1);
    assert_eq!(report.remapped[0].handle,second);
    assert_eq!(report.remapped[0].before,second_name_before);
    assert_eq!(e.resolve(first),Err(StructuralHandleError::StaleHole));
    assert_eq!(e.resolve(second).unwrap().topology.side_face,FaceId(6));
    assert_eq!(e.resolve(second).unwrap().handle,second);
    prove(&e,&[h(1.5,0.,0.7)]);
    let old_first=first;
    let undo=e.undo().unwrap().unwrap();
    assert_eq!(undo.direction,HistoryDirection::Undo);
    prove(&e,&[h(-1.5,0.,0.6),h(1.5,0.,0.7)]);
    let reintroduced=e.handle(first.slot).unwrap();
    assert_ne!(reintroduced,old_first);
    assert_eq!(e.resolve(old_first),Err(StructuralHandleError::StaleHole));
    assert_eq!(e.resolve(second).unwrap().handle,second);
    let redo=e.redo().unwrap().unwrap();
    assert_eq!(redo.direction,HistoryDirection::Redo);
    prove(&e,&[h(1.5,0.,0.7)]);
    assert_eq!(e.resolve(old_first),Err(StructuralHandleError::StaleHole));
}

#[test]
fn removed_hole_slot_is_reused_with_a_new_incarnation(){
    let mut e=editor();
    let first={let mut t=e.begin_transaction();let a=t.add_through_hole(h(-1.5,0.,0.5)).unwrap();t.commit().unwrap();a};
    {let mut t=e.begin_transaction();t.remove_through_hole(first).unwrap();t.commit().unwrap();}
    let next={let mut t=e.begin_transaction();let a=t.add_through_hole(h(1.5,0.,0.5)).unwrap();t.commit().unwrap();a};
    assert_eq!(next.slot,first.slot);
    assert_ne!(next.generation,first.generation);
    assert_eq!(e.resolve(first),Err(StructuralHandleError::StaleHole));
    prove(&e,&[h(1.5,0.,0.5)]);
}
#[test]
fn rejected_stage_poisons_entire_batch_and_preserves_identity(){
    let mut e=editor();
    let original=e.snapshot();
    let mut t=e.begin_transaction();
    t.add_through_hole(h(-1.5,0.,0.5)).unwrap();
    assert!(t.add_through_hole(h(-1.4,0.,0.6)).is_err());
    assert!(t.commit().is_err());
    assert_eq!(e.revision(),0);
    assert_eq!(e.undo_depth(),0);
    assert_eq!(e.bindings().len(),0);
    e.check_snapshot(&original).unwrap();
    prove(&e,&[]);
}
#[test]
fn rollback_on_invalid_late_carrier_edit_and_invalid_geometry(){
    let mut e=editor();
    let before=e.snapshot();
    let mut t=e.begin_transaction();
    let first=t.add_through_hole(h(-1.5,0.,0.5)).unwrap();
    t.reparameterize_edge(StructuralEdgeRef::Hole{
        hole:first,role:HoleEdgeRole::BottomCircle},2.0,0.25).unwrap();
    assert!(t.commit().is_err(),"circle phase scale two is forbidden");
    eq(e.model().mass.volume,144.0);
    e.check_snapshot(&before).unwrap();
    assert_eq!(e.undo_depth(),0);
    let mut t=e.begin_transaction();
    let first=t.add_through_hole(h(-1.5,0.,0.5)).unwrap();
    t.rebase_face_uv(StructuralFaceRef::HoleSide(first),0.125,0.0).unwrap();
    assert!(t.commit().is_err(),"periodic U chart rebase is forbidden");
    assert_eq!(e.revision(),0);
    e.check_snapshot(&before).unwrap();
    prove(&e,&[]);
}
#[test]
fn resize_hole_recomputes_analytic_mass_and_preserves_handle(){
    let mut e=editor();
    let first={let mut t=e.begin_transaction();let a=t.add_through_hole(h(0.,0.,1.)).unwrap();t.commit().unwrap();a};
    let before=e.resolve(first).unwrap();
    let mut t=e.begin_transaction();
    t.resize_through_hole(first,h(0.2,-0.4,1.5)).unwrap();
    let report=t.commit().unwrap();
    assert!(report.remapped.is_empty());
    let after=e.resolve(first).unwrap();
    assert_eq!(before.handle,after.handle);
    assert_eq!(before.topology,after.topology);
    prove(&e,&[h(0.2,-0.4,1.5)]);
    e.undo().unwrap().unwrap();
    prove(&e,&[h(0.,0.,1.)]);
}
#[test]
fn retained_carrier_intent_follows_survivor_after_compaction(){
    let mut e=editor();
    let (first,second)={
        let mut t=e.begin_transaction();
        let a=t.add_through_hole(h(-1.6,0.,0.6)).unwrap();
        let b=t.add_through_hole(h(1.6,0.,0.6)).unwrap();
        t.reparameterize_edge(StructuralEdgeRef::Hole{hole:b,
            role:HoleEdgeRole::BottomCircle},1.,0.4).unwrap();
        t.commit().unwrap();(a,b)
    };
    let before=e.resolve(second).unwrap().topology.bottom_edge;
    assert_eq!(before,EdgeId(15));
    let mut t=e.begin_transaction();
    t.remove_through_hole(first).unwrap();
    t.commit().unwrap();
    assert_eq!(e.resolve(second).unwrap().topology.bottom_edge,EdgeId(12));
    let Curve3::Circle(circle)=e.model().edge_curve(EdgeId(12)).unwrap()
        else {panic!("surviving hole must retain its curved carrier")};
    assert!((circle.frame.x.x-0.4f64.cos()).abs()<1e-10);
    prove(&e,&[h(1.6,0.,0.6)]);
}
#[test]
fn invalid_foreign_and_dropped_staged_hole_refs_remain_stale(){
    let mut a=editor();
    let mut b=editor();
    let abandoned={
        let mut t=a.begin_transaction();
        t.add_through_hole(h(0.,0.,1.)).unwrap()
    };
    let actual={
        let mut t=a.begin_transaction();
        let h=t.add_through_hole(h(0.,0.,1.)).unwrap();
        t.commit().unwrap();h
    };
    assert_eq!(abandoned.slot,actual.slot);
    assert_ne!(abandoned.generation,actual.generation);
    assert_eq!(a.resolve(abandoned),Err(StructuralHandleError::StaleHole));
    let first={
        let mut t=b.begin_transaction();
        let h=t.add_through_hole(h(0.,0.,1.)).unwrap();
        t.commit().unwrap();h
    };
    assert_eq!(a.resolve(first),Err(StructuralHandleError::ForeignEditor));
    let snap=a.snapshot();
    let mut t=a.begin_transaction();
    assert!(t.remove_through_hole(first).is_err());
    assert!(t.commit().is_err());
    a.check_snapshot(&snap).unwrap();
    prove(&a,&[h(0.,0.,1.)]);
}
#[test]
fn genus_growth_to_three_and_back_is_mathematically_certified(){
    let mut e=editor();
    let mut t=e.begin_transaction();
    let a=t.add_through_hole(h(-1.75,-0.5,0.5)).unwrap();
    let b=t.add_through_hole(h(0.,1.25,0.6)).unwrap();
    let c=t.add_through_hole(h(1.75,-0.5,0.4)).unwrap();
    t.commit().unwrap();
    prove(&e,&[h(-1.75,-0.5,0.5),h(0.,1.25,0.6),h(1.75,-0.5,0.4)]);
    let mut t=e.begin_transaction();
    t.remove_through_hole(b).unwrap();
    t.commit().unwrap();
    prove(&e,&[h(-1.75,-0.5,0.5),h(1.75,-0.5,0.4)]);
    assert_eq!(e.resolve(a).unwrap().handle,a);
    assert_eq!(e.resolve(c).unwrap().handle,c);
    assert_eq!(e.resolve(c).unwrap().topology.side_face,FaceId(7));
}
#[test]
fn original_through_hole_can_be_promoted_into_structural_editor(){
    let block=block();
    let donor=ThroughHoleBrep::from_box(&block,h(-1.2,0.,0.75)).unwrap();
    let source_mass=donor.model().mass;
    let mut e=BrepStructuralEditor::from_body(&BrepBody::ThroughHole(donor)).unwrap();
    prove(&e,&[h(-1.2,0.,0.75)]);
    eq(source_mass.volume,e.model().mass.volume);
    let before=e.handle(0).unwrap();
    let mut t=e.begin_transaction();
    t.add_through_hole(h(1.2,0.,0.75)).unwrap();
    t.commit().unwrap();
    assert_eq!(e.resolve(before).unwrap().handle,before);
    prove(&e,&[h(-1.2,0.,0.75),h(1.2,0.,0.75)]);
}
#[test]
fn stage_order_between_structure_and_uv_trim_carriers_remains_atomic(){
    let mut e=editor();
    let mut t=e.begin_transaction();
    t.rebase_face_uv(StructuralFaceRef::Base(FaceId(0)),0.2,-0.4).unwrap();
    let first=t.add_through_hole(h(0.,0.,1.0)).unwrap();
    t.reparameterize_trim(StructuralTrimRef::Hole{hole:first,
        role:HoleTrimRole::BottomCap},1.,TAU).unwrap();
    t.commit().unwrap();
    prove(&e,&[h(0.,0.,1.)]);
    assert_eq!(e.model().faces[0].loops.len(),2);
    assert_eq!(e.undo_depth(),1);
}
#[test]
fn subsequent_commit_after_undo_invalidates_redo_branch(){
    let mut e=editor();
    let first={let mut t=e.begin_transaction();let a=t.add_through_hole(h(0.,0.,0.5)).unwrap();t.commit().unwrap();a};
    let mut t=e.begin_transaction();
    t.resize_through_hole(first,h(0.,0.,0.7)).unwrap();
    t.commit().unwrap();
    e.undo().unwrap().unwrap();
    assert_eq!(e.redo_depth(),1);
    let mut t=e.begin_transaction();
    t.resize_through_hole(first,h(0.,0.,0.9)).unwrap();
    t.commit().unwrap();
    assert_eq!(e.redo_depth(),0);
    assert!(e.redo().unwrap().is_none());
    prove(&e,&[h(0.,0.,0.9)]);
}

#[test]
fn detach_from_source_arena_preserves_original_handles_history_and_mass(){
    let source=BrepArena::from_solid(block()).unwrap();
    let original=source.snapshot().unwrap();
    let original_face=source.handle(rustsolid::TopologyEntity::Face(FaceId(0))).unwrap();
    let mut detached=source.detached_structural_editor().unwrap();
    let mut tx=detached.begin_transaction();
    tx.add_through_hole(h(0.,0.,1.)).unwrap();
    tx.commit().unwrap();
    detached.validate().unwrap();
    prove(&detached,&[h(0.,0.,1.)]);
    source.check_snapshot(&original).unwrap();
    assert_eq!(source.handle(rustsolid::TopologyEntity::Face(FaceId(0))),Some(original_face));
    assert_eq!(source.undo_depth(),0);
    eq(source.model().unwrap().mass.volume,144.);
}
#[test]
fn overlap_in_resize_rejects_without_losing_previous_history(){
    let mut e=editor();
    let (first,second)={
        let mut tx=e.begin_transaction();
        let a=tx.add_through_hole(h(-1.5,0.,0.5)).unwrap();
        let b=tx.add_through_hole(h(1.5,0.,0.5)).unwrap();
        tx.commit().unwrap();(a,b)
    };
    let snapshot=e.snapshot();
    let original_undo=e.undo_depth();
    let mut tx=e.begin_transaction();
    assert!(tx.resize_through_hole(first,h(1.0,0.,0.6)).is_err());
    assert!(tx.resize_through_hole(second,h(1.2,0.,0.4)).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(e.undo_depth(),original_undo);
    assert_eq!(e.revision(),1);
    e.check_snapshot(&snapshot).unwrap();
    prove(&e,&[h(-1.5,0.,0.5),h(1.5,0.,0.5)]);
}
