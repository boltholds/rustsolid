//! First-party mathematical, atomicity, generation and edit-history contracts.
//! A second CAD engine is never the correctness oracle.
use rustsolid::{
    ArenaHandle, ArenaHandleError, ArenaOutput, BrepArena, BrepBody, CylindricalBrep,
    EdgeId, FaceId, GeometryTolerance, IntersectionGraphBuilder, Point3,
    Solid, TopologyDeltaPlanner, TopologyEntity, TopologyName, VertexId,
};
use serde::Deserialize;
use std::collections::BTreeSet;

const CASES:&str=include_str!("fixtures/arena_transactions.v1.json");
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {schema:String,title:String,cases:Vec<Case>}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {id:String,kind:String}

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn cube()->BrepArena {BrepArena::from_solid(Solid::block(p(0.,0.,0.),2.,3.,4.).unwrap()).unwrap()}
fn h(arena:&BrepArena,entity:TopologyEntity)->ArenaHandle{arena.handle(entity).expect("missing original entity")}
fn mass(arena:&BrepArena){
    let model=arena.model().unwrap();
    assert!((model.mass.volume-24.0).abs()<1e-10);
    assert!((model.mass.surface_area-52.0).abs()<1e-10);
    assert_eq!(model.euler_characteristic(),2);
    for (actual,expected) in model.mass.centroid.array().into_iter().zip([1.0,1.5,2.0]) {
        assert!((actual-expected).abs()<1e-10);
    }
    arena.validate().unwrap();
}

fn two_make_commands_one_atomic_commit(){
    let mut arena=cube();
    let before=arena.snapshot().unwrap();
    let old=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let old_face=h(&arena,TopologyEntity::Face(FaceId(0)));
    let face=arena.model().unwrap().faces[0].id;
    let start=h(&arena,TopologyEntity::Vertex(VertexId(0)));
    let end=h(&arena,TopologyEntity::Vertex(VertexId(2)));
    assert_eq!(face,FaceId(0));
    let mut tx=arena.begin_transaction().unwrap();
    tx.split_edge(old,0.25).unwrap();
    tx.split_face(old_face,start,end).unwrap();
    assert_eq!(tx.staged_len(),2);
    let applied=tx.commit().unwrap();
    assert_eq!((applied.report.revision_before,applied.report.revision_after),(0,1));
    assert_eq!(applied.report.changes.len(),2);
    assert_eq!(arena.undo_depth(),1);
    assert_eq!(arena.redo_depth(),0);
    assert_eq!((arena.model().unwrap().vertices.len(),arena.model().unwrap().faces.len()),(9,7));
    assert_eq!(arena.resolve(old).unwrap().entity,TopologyEntity::Edge(EdgeId(0)));
    assert_eq!(arena.resolve(old_face).unwrap().entity,TopologyEntity::Face(FaceId(0)));
    assert!(!applied.created.is_empty());
    let fresh=applied.created.iter().find(|b|matches!(b.entity,TopologyEntity::Vertex(VertexId(8))))
        .expect("MEV must add a named vertex");
    let saved=fresh.handle;
    mass(&arena);
    assert!(arena.check_snapshot(&before).is_err());
    // Undo is ONE modeling step, not two independent commands.
    let undone=arena.undo().unwrap().unwrap();
    assert_eq!((undone.revision_before,undone.revision_after),(1,2));
    assert_eq!(arena.undo_depth(),0);
    assert_eq!((arena.model().unwrap().vertices.len(),arena.model().unwrap().faces.len()),(8,6));
    assert_eq!(arena.resolve(saved),Err(ArenaHandleError::StaleHandle));
    mass(&arena);
    arena.redo().unwrap().unwrap();
    let reincarnated=h(&arena,TopologyEntity::Vertex(VertexId(8)));
    assert_ne!(saved,reincarnated);
    assert_eq!(arena.resolve(saved),Err(ArenaHandleError::StaleHandle));
    assert_eq!((arena.model().unwrap().vertices.len(),arena.model().unwrap().faces.len()),(9,7));
    mass(&arena);
}
fn late_failure_preserves_everything(){
    let mut arena=cube();
    let before=arena.snapshot().unwrap();
    let handles=arena.bindings();
    let original_slots=arena.allocated_slots();
    let mut tx=arena.begin_transaction().unwrap();
    tx.split_edge(before.handle(TopologyEntity::Edge(EdgeId(0))).unwrap(),0.5).unwrap();
    tx.split_face(before.handle(TopologyEntity::Face(FaceId(0))).unwrap(),
        before.handle(TopologyEntity::Vertex(VertexId(0))).unwrap(),
        before.handle(TopologyEntity::Vertex(VertexId(1))).unwrap()).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(arena.revision(),0);
    assert_eq!(arena.undo_depth(),0);
    assert_eq!(arena.allocated_slots(),original_slots);
    assert_eq!(arena.bindings(),handles);
    arena.check_snapshot(&before).unwrap();
    assert_eq!((arena.model().unwrap().vertices.len(),arena.model().unwrap().faces.len()),(8,6));
    mass(&arena);
}
fn dropped_or_mistyped_transaction_cannot_modify_topology(){
    let mut arena=cube();
    let edge=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let face=h(&arena,TopologyEntity::Face(FaceId(0)));
    let vertex=h(&arena,TopologyEntity::Vertex(VertexId(0)));
    let prior=arena.bindings();
    {
        let mut tx=arena.begin_transaction().unwrap();
        tx.split_edge(edge,0.5).unwrap();
        assert!(tx.split_edge(face,0.5).is_err());
        assert!(tx.split_face(edge,vertex,vertex).is_err());
    }
    assert_eq!(arena.revision(),0);
    assert_eq!(prior,arena.bindings());
    assert_eq!(arena.undo_depth(),0);
    mass(&arena);
}
fn kill_compaction_retains_semantic_survivor_handle(){
    let mut arena=cube();
    let edge=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let mut first=arena.begin_transaction().unwrap();
    first.split_edge(edge,0.5).unwrap();first.commit().unwrap();
    let removed=h(&arena,TopologyEntity::Vertex(VertexId(8)));
    let other_edge=h(&arena,TopologyEntity::Edge(EdgeId(5)));
    let mut next=arena.begin_transaction().unwrap();
    next.split_edge(other_edge,0.5).unwrap();next.commit().unwrap();
    let survivor=h(&arena,TopologyEntity::Vertex(VertexId(9)));
    let original_name=arena.resolve(survivor).unwrap().name;
    let mut tx=arena.begin_transaction().unwrap();
    tx.kill_edge_vertex(removed).unwrap();
    let commit=tx.commit().unwrap();
    assert_eq!(arena.resolve(removed),Err(ArenaHandleError::StaleHandle));
    assert_eq!(arena.resolve(survivor).unwrap().entity,TopologyEntity::Vertex(VertexId(8)));
    assert_eq!(arena.handle_name(&original_name),Some(survivor));
    assert!(commit.remapped.iter().any(|item|item.name==original_name
        && item.before==TopologyEntity::Vertex(VertexId(9))
        && item.after==TopologyEntity::Vertex(VertexId(8))));
    assert!(commit.retired.iter().any(|item|item.handle==removed));
    mass(&arena);
    arena.undo().unwrap().unwrap();
    // The same semantic survivor stays stable when the compacted slot changes.
    assert_eq!(arena.resolve(survivor).unwrap().entity,TopologyEntity::Vertex(VertexId(9)));
    arena.redo().unwrap().unwrap();
    assert_eq!(arena.resolve(survivor).unwrap().entity,TopologyEntity::Vertex(VertexId(8)));
    mass(&arena);
}
fn same_slot_reused_after_undo_is_never_aba(){
    let mut arena=cube();
    let e=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let slots=arena.allocated_slots();
    let mut previous:Vec<ArenaHandle>=Vec::new();
    let mut observed_reuse=false;
    for _ in 0..7 {
        let mut tx=arena.begin_transaction().unwrap();
        tx.split_edge(e,0.5).unwrap();
        let result=tx.commit().unwrap();
        let created=result.created.iter().find(|b|b.entity==TopologyEntity::Vertex(VertexId(8))).unwrap().handle;
        for old in &previous {
            assert_eq!(arena.resolve(*old),Err(ArenaHandleError::StaleHandle));
        }
        for newly_allocated in &result.created {
            for old in &previous {
                if newly_allocated.handle.slot==old.slot {
                    observed_reuse=true;
                    assert_ne!(newly_allocated.handle.generation,old.generation);
                }
            }
        }
        previous=result.created.iter().map(|b|b.handle).collect();
        arena.undo().unwrap().unwrap();
        assert_eq!(arena.resolve(created),Err(ArenaHandleError::StaleHandle));
        for old in &previous {assert_eq!(arena.resolve(*old),Err(ArenaHandleError::StaleHandle));}
    }
    assert!(observed_reuse,"allocator should recycle free slots over repeated edit/undo");
    assert_eq!(arena.allocated_slots(),slots+4);
    assert_eq!(arena.live_slots(),slots);
    mass(&arena);
}
fn other_body_handles_never_resolve(){
    let mut a=cube();let b=cube();
    let foreign=h(&b,TopologyEntity::Vertex(VertexId(0)));
    assert_eq!(a.resolve(foreign),Err(ArenaHandleError::ForeignArena));
    let mut tx=a.begin_transaction().unwrap();
    assert!(tx.kill_edge_vertex(foreign).is_err());
    assert_eq!(a.revision(),0);
    mass(&a);
}
fn curved_body_handles_are_queryable_but_not_editable(){
    let cylinder=CylindricalBrep::upright(p(0.,0.,0.),2.,5.,GeometryTolerance::default()).unwrap();
    let mut arena=BrepArena::new(BrepBody::Cylinder(cylinder)).unwrap();
    assert!(!arena.is_editable());
    let snap=arena.snapshot().unwrap();
    assert_eq!(snap.model.summary().circular_edges,2);
    assert_eq!(snap.model.summary().seam_edges,1);
    assert_eq!(snap.model.summary().cylindrical_faces,1);
    let handle=arena.handle(TopologyEntity::Edge(EdgeId(0))).unwrap();
    assert_eq!(arena.resolve(handle).unwrap().entity,TopologyEntity::Edge(EdgeId(0)));
    assert!(arena.begin_transaction().is_err());
    arena.check_snapshot(&snap).unwrap();
    arena.validate().unwrap();
}
fn unclassified_intersection_plan_cannot_be_committed(){
    let arena=BrepArena::from_solid(Solid::block(p(-3.,2.,-3.),6.,1.,6.).unwrap()).unwrap();
    let target=CylindricalBrep::upright(p(0.,0.,0.),2.,5.,GeometryTolerance::default()).unwrap();
    let other=BrepArena::from_cylinder(target).unwrap();
    let snapshot=arena.snapshot().unwrap();
    let second=other.model().unwrap();
    let graph=IntersectionGraphBuilder::new(GeometryTolerance::default()).unwrap()
        .intersect_faces(&snapshot.model,FaceId(0),&second,FaceId(2)).unwrap();
    let plan=TopologyDeltaPlanner::new().plan(&graph,&snapshot.model,FaceId(0),&second,FaceId(2)).unwrap();
    assert!(!plan.is_mutation_ready());
    assert!(arena.require_executable_plan(&plan).is_err());
    arena.check_snapshot(&snapshot).unwrap();
    assert_eq!(arena.model().unwrap().mass.volume,36.0);
    arena.validate().unwrap();
}
fn kill_face_restores_original_analytic_mass_and_face_handles(){
    let mut arena=cube();
    let face=h(&arena,TopologyEntity::Face(FaceId(0)));
    let model=arena.model().unwrap();
    let vertices=model.coedges.iter().filter(|c|c.face==FaceId(0))
        .map(|c|c.start_vertex(&model.edges[c.edge.0 as usize])).collect::<Vec<_>>();
    assert_eq!(vertices.len(),4);
    let from=h(&arena,TopologyEntity::Vertex(vertices[0]));
    let to=h(&arena,TopologyEntity::Vertex(vertices[2]));
    let mut tx=arena.begin_transaction().unwrap();
    tx.split_face(face,from,to).unwrap();
    tx.commit().unwrap();
    let diagonal=h(&arena,TopologyEntity::Edge(EdgeId(12)));
    let removed=h(&arena,TopologyEntity::Face(FaceId(6)));
    let mut tx=arena.begin_transaction().unwrap();
    tx.kill_edge_face(diagonal,removed).unwrap();
    let report=tx.commit().unwrap();
    assert!(report.retired.iter().any(|item|item.handle==removed));
    assert_eq!(arena.resolve(removed),Err(ArenaHandleError::StaleHandle));
    assert_eq!(arena.model().unwrap().faces.len(),6);
    mass(&arena);
    arena.undo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().faces.len(),7);
    arena.redo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().faces.len(),6);
    mass(&arena);
}
fn existing_ids_reject_forged_generation(){
    let arena=cube();
    let mut h=h(&arena,TopologyEntity::Body);
    h.generation+=1;
    assert_eq!(arena.resolve(h),Err(ArenaHandleError::StaleHandle));
}

fn dependent_remainder_edge_chain_is_one_undoable_revision() {
    let mut arena = cube();
    let before = arena.snapshot().unwrap();
    let original = h(&arena, TopologyEntity::Edge(EdgeId(0)));
    let mut tx = arena.begin_transaction().unwrap();
    let first = tx.split_edge_ref(original, 0.25).unwrap();
    assert_eq!(first.index(), 0);
    let second = tx.split_edge_ref(first.output(ArenaOutput::RemainderEdge), 0.5).unwrap();
    assert_eq!(second.index(), 1);
    assert_eq!(tx.staged_len(), 2);
    let commit = tx.commit().unwrap();
    assert_eq!(commit.report.changes.len(), 2);
    assert_eq!((commit.report.revision_before, commit.report.revision_after), (0, 1));
    assert_eq!((arena.undo_depth(), arena.redo_depth()), (1, 0));
    assert_eq!(arena.model().unwrap().vertices.len(), 10);
    let fresh = commit.created.iter().find(|item|item.entity==TopologyEntity::Vertex(VertexId(9)))
        .expect("second MEV vertex").handle;
    arena.check_snapshot(&before).unwrap_err();
    mass(&arena);
    arena.undo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().vertices.len(), 8);
    assert_eq!(arena.resolve(fresh), Err(ArenaHandleError::StaleHandle));
    arena.redo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().vertices.len(), 10);
    assert_eq!(arena.resolve(fresh), Err(ArenaHandleError::StaleHandle));
    mass(&arena);
}

fn dependent_mev_mef_and_diagonal_edge_can_chain() {
    let mut arena = cube();
    let model = arena.model().unwrap();
    let edge = &model.edges[0];
    let incident_face = model.coedges.iter().find(|coedge|coedge.edge==edge.id)
        .expect("edge must be on a face").face;
    let remote = model.coedges.iter().filter(|coedge|coedge.face==incident_face)
        .map(|coedge|coedge.start_vertex(&model.edges[coedge.edge.0 as usize]))
        .find(|vertex|*vertex!=edge.start && *vertex!=edge.end)
        .expect("quadrilateral face has another vertex");
    let edge_handle = h(&arena, TopologyEntity::Edge(edge.id));
    let face_handle = h(&arena, TopologyEntity::Face(incident_face));
    let remote_handle = h(&arena, TopologyEntity::Vertex(remote));
    let mut tx = arena.begin_transaction().unwrap();
    let mev = tx.split_edge_ref(edge_handle, 0.4).unwrap();
    let mef = tx.split_face_ref(
        face_handle, mev.output(ArenaOutput::InsertedVertex), remote_handle
    ).unwrap();
    tx.split_edge_ref(mef.output(ArenaOutput::DiagonalEdge), 0.5).unwrap();
    let committed = tx.commit().unwrap();
    assert_eq!(committed.report.changes.len(), 3);
    assert_eq!(arena.revision(), 1);
    assert_eq!(arena.undo_depth(), 1);
    assert_eq!(arena.model().unwrap().vertices.len(), 10);
    assert_eq!(arena.model().unwrap().faces.len(), 7);
    mass(&arena);
    arena.undo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().faces.len(), 6);
    arena.redo().unwrap().unwrap();
    assert_eq!(arena.model().unwrap().faces.len(), 7);
    mass(&arena);
}

fn kill_created_vertex_in_same_batch_keeps_single_memento() {
    let mut arena = cube();
    let original = arena.snapshot().unwrap();
    let edge = h(&arena, TopologyEntity::Edge(EdgeId(0)));
    let mut tx = arena.begin_transaction().unwrap();
    let mev = tx.split_edge_ref(edge, 0.5).unwrap();
    tx.kill_edge_vertex_ref(mev.output(ArenaOutput::InsertedVertex)).unwrap();
    let committed = tx.commit().unwrap();
    assert_eq!(committed.report.changes.len(), 2);
    assert_eq!(arena.revision(), 1);
    assert_eq!(arena.undo_depth(), 1);
    assert_eq!(arena.model().unwrap().vertices.len(), 8);
    mass(&arena);
    arena.undo().unwrap().unwrap();
    arena.redo().unwrap().unwrap();
    arena.check_snapshot(&original).unwrap_err();
    mass(&arena);
}

fn foreign_transaction_output_rejects_and_poisons_staging() {
    let mut arena = cube();
    let initial = arena.snapshot().unwrap();
    let original = arena.bindings();
    let edge = h(&arena, TopologyEntity::Edge(EdgeId(0)));
    let previous = {
        let mut abandoned = arena.begin_transaction().unwrap();
        abandoned.split_edge_ref(edge, 0.5).unwrap()
    };
    let mut tx = arena.begin_transaction().unwrap();
    tx.split_edge_ref(edge, 0.25).unwrap();
    assert!(tx.kill_edge_vertex_ref(previous.output(ArenaOutput::InsertedVertex)).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(arena.revision(), 0);
    assert_eq!(arena.bindings(), original);
    arena.check_snapshot(&initial).unwrap();
    mass(&arena);
}

fn wrong_output_role_rejects_and_poisons_staging() {
    let mut arena = cube();
    let original = arena.snapshot().unwrap();
    let edge = h(&arena, TopologyEntity::Edge(EdgeId(0)));
    let vertex = h(&arena, TopologyEntity::Vertex(VertexId(0)));
    let mut tx = arena.begin_transaction().unwrap();
    let mev = tx.split_edge_ref(edge, 0.5).unwrap();
    assert!(tx.split_face_ref(
        mev.output(ArenaOutput::CreatedFace), vertex, vertex
    ).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(arena.revision(), 0);
    arena.check_snapshot(&original).unwrap();
    mass(&arena);
}

fn dependent_late_geometry_failure_rolls_back_preflight() {
    let mut arena = cube();
    let initial = arena.snapshot().unwrap();
    let original = arena.bindings();
    let original_slots = arena.allocated_slots();
    let edge = h(&arena, TopologyEntity::Edge(EdgeId(0)));
    let model = arena.model().unwrap();
    let face = model.coedges.iter().find(|coedge|coedge.edge==EdgeId(0))
        .expect("incident face").face;
    let face_handle = h(&arena, TopologyEntity::Face(face));
    let mut tx = arena.begin_transaction().unwrap();
    let mev = tx.split_edge_ref(edge, 0.5).unwrap();
    let vertex = mev.output(ArenaOutput::InsertedVertex);
    tx.split_face_ref(face_handle, vertex, vertex).unwrap();
    assert!(tx.commit().is_err(),"duplicate vertices must reject MEF");
    assert_eq!(arena.revision(), 0);
    assert_eq!(arena.undo_depth(), 0);
    assert_eq!(arena.redo_depth(), 0);
    assert_eq!(arena.bindings(), original);
    assert_eq!(arena.allocated_slots(), original_slots);
    arena.check_snapshot(&initial).unwrap();
    mass(&arena);
}
#[test]
fn pinned_arena_atomicity_identity_and_geometry_witnesses(){
    let m:Manifest=serde_json::from_str(CASES).unwrap();
    assert_eq!(m.schema,"rustsolid.brep_arena.v1");
    assert!(!m.title.is_empty());
    assert_eq!(m.cases.len(),16);
    let mut ids=BTreeSet::new();
    for c in m.cases {
        assert!(ids.insert(c.id.clone()),"duplicate scenario id");
        match c.kind.as_str(){
            "atomic_split_batch"=>two_make_commands_one_atomic_commit(),
            "late_failure_rollback"=>late_failure_preserves_everything(),
            "dropped_and_mistyped"=>dropped_or_mistyped_transaction_cannot_modify_topology(),
            "kill_compaction"=>kill_compaction_retains_semantic_survivor_handle(),
            "generation_aba"=>same_slot_reused_after_undo_is_never_aba(),
            "foreign_body"=>other_body_handles_never_resolve(),
            "curved_read_only"=>curved_body_handles_are_queryable_but_not_editable(),
            "plan_rejected"=>unclassified_intersection_plan_cannot_be_committed(),
            "forged_generation"=>existing_ids_reject_forged_generation(),
            "kill_face_compaction"=>kill_face_restores_original_analytic_mass_and_face_handles(),
            "dependent_edge_chain"=>dependent_remainder_edge_chain_is_one_undoable_revision(),
            "dependent_mev_mef"=>dependent_mev_mef_and_diagonal_edge_can_chain(),
            "dependent_same_batch_inverse"=>kill_created_vertex_in_same_batch_keeps_single_memento(),
            "foreign_transaction_output"=>foreign_transaction_output_rejects_and_poisons_staging(),
            "invalid_output_role"=>wrong_output_role_rejects_and_poisons_staging(),
            "dependent_late_rollback"=>dependent_late_geometry_failure_rolls_back_preflight(),
            other=>panic!("unimplemented executable arena contract {other}"),
        }
    }
}
#[test]
fn permanent_manifest_cannot_accidentally_drop_quality_families(){
    let m:Manifest=serde_json::from_str(CASES).unwrap();
    for required in ["atomic_split_batch","late_failure_rollback","kill_compaction",
        "generation_aba","foreign_body","curved_read_only","plan_rejected","kill_face_compaction",
        "dependent_edge_chain","dependent_mev_mef","dependent_same_batch_inverse",
        "foreign_transaction_output","invalid_output_role","dependent_late_rollback"]{
        assert!(m.cases.iter().any(|c|c.kind==required));
    }
}
#[test]
fn failed_undo_and_stale_snapshot_does_not_mutate_registry(){
    let mut arena=cube();
    let original=arena.snapshot().unwrap();
    assert!(arena.undo().unwrap().is_none());
    assert!(arena.redo().unwrap().is_none());
    arena.check_snapshot(&original).unwrap();
    let edge=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let mut tx=arena.begin_transaction().unwrap();
    tx.split_edge(edge,0.5).unwrap();tx.commit().unwrap();
    assert!(arena.check_snapshot(&original).is_err());
    assert_eq!(arena.revision(),1);
    arena.undo().unwrap().unwrap();
    assert!(arena.check_snapshot(&original).is_err()); // monotonic runtime revision
    mass(&arena);
}

fn failed_staging_poison_prevents_partial_commit(){
    let mut arena=cube();
    let before=arena.snapshot().unwrap();
    let original_bindings=arena.bindings();
    let original_slots=arena.allocated_slots();
    let edge=h(&arena,TopologyEntity::Edge(EdgeId(0)));
    let face=h(&arena,TopologyEntity::Face(FaceId(0)));
    let vertex=h(&arena,TopologyEntity::Vertex(VertexId(0)));
    {
        let mut tx=arena.begin_transaction().unwrap();
        tx.split_edge(edge,0.5).unwrap();
        assert_eq!(tx.staged_len(),1);
        assert!(tx.split_edge(face,0.5).is_err(),"wrong-entity handle must poison transaction");
        assert!(tx.split_face(face,vertex,vertex).is_err(),"subsequent stages are blocked");
        assert!(tx.commit().is_err(),"ignore-Err must not commit previous staged changes");
    }
    assert_eq!(arena.revision(),0);
    assert_eq!(arena.bindings(),original_bindings);
    assert_eq!(arena.allocated_slots(),original_slots);
    assert_eq!((arena.undo_depth(),arena.redo_depth()),(0,0));
    arena.check_snapshot(&before).unwrap();
    mass(&arena);
    {
        let mut tx=arena.begin_transaction().unwrap();
        tx.split_edge(edge,0.5).unwrap();
        assert!(tx.split_edge(edge,f64::NAN).is_err(),"invalid fraction poisons transaction too");
        assert!(tx.commit().is_err());
    }
    assert_eq!(arena.bindings(),original_bindings);
    arena.check_snapshot(&before).unwrap();
    {
        let mut tx=arena.begin_transaction().unwrap();
        tx.split_edge(edge,0.5).unwrap();
        assert!(tx.split_edge(ArenaHandle{arena_token:u64::MAX,..edge},0.5).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(arena.revision(),0);
    arena.check_snapshot(&before).unwrap();
    mass(&arena);
}
#[test]
fn staging_failure_poison_is_atomic_even_when_caller_ignores_error(){
    failed_staging_poison_prevents_partial_commit();
}
