use rustsolid::{
    CommandBatch, CommandHistory, EdgeId, EulerCommand, FaceId, GeometryError,
    HandleError, HistoryDirection, Point3, Solid, TopologyEntity, TopologyName, VertexId,
};

fn cube() -> Solid {
    Solid::block(Point3 { x: 1.0, y: -2.0, z: 7.0 }, 2.0, 3.0, 4.0)
        .unwrap().with_feature_key("assembly/part/feature-1").unwrap()
}

fn name(solid: &Solid, entity: TopologyEntity) -> TopologyName {
    solid.topology_name(entity).unwrap().clone()
}

fn edge_command(solid: &Solid) -> EulerCommand {
    EulerCommand::SplitEdge {
        edge: name(solid, TopologyEntity::Edge(EdgeId(0))), fraction: 0.5,
    }
}

fn face_command(solid: &Solid) -> EulerCommand {
    let face = FaceId(0);
    let corners = &solid.faces[face.0 as usize].boundary;
    EulerCommand::SplitFace {
        face: name(solid, TopologyEntity::Face(face)),
        start: name(solid, TopologyEntity::Vertex(corners[0])),
        end: name(solid, TopologyEntity::Vertex(corners[2])),
    }
}

#[test]
fn command_plus_memento_undo_and_redo_are_reversible() {
    let mut history = CommandHistory::new(cube()).unwrap();
    let retained = history.solid().topology_handle(TopologyEntity::Edge(EdgeId(0))).unwrap();
    let before_mesh = history.solid().mesh.triangles.clone();
    let before_mass = history.solid().mass;
    let report = history.execute(CommandBatch { commands: vec![
        edge_command(history.solid()), face_command(history.solid()),
    ] }).unwrap();
    assert_eq!((report.revision_before, report.revision_after), (0, 1));
    assert_eq!(history.solid().vertices.len(), 9);
    assert_eq!(history.solid().faces.len(), 7);
    assert_eq!(history.solid().euler_characteristic(), 2);
    assert_eq!(history.solid().mass, before_mass);
    assert!(report.journal.frames == 2);
    assert!(report.journal.topology_snapshots < history.solid().coedges.len());
    let new_handle = history.solid().topology_handle(TopologyEntity::Vertex(VertexId(8))).unwrap();
    let new_name = name(history.solid(), TopologyEntity::Vertex(VertexId(8)));
    assert_eq!(history.undo_depth(), 1);

    let undone = history.undo().unwrap().unwrap();
    assert_eq!(undone.direction, HistoryDirection::Undo);
    assert_eq!((undone.revision_before, undone.revision_after), (1, 2));
    assert_eq!((history.solid().vertices.len(),history.solid().faces.len()), (8, 6));
    assert_eq!(history.solid().mesh.triangles, before_mesh);
    assert_eq!(history.solid().mass, before_mass);
    assert!(history.solid().resolve_topology_name(&new_name).is_none());
    assert_eq!(history.solid().resolve_topology_handle(new_handle), Err(HandleError::StaleHandle));
    assert_eq!(history.solid().resolve_topology_handle(retained), Ok(TopologyEntity::Edge(EdgeId(0))));
    assert_eq!(history.command_log().len(), 0);
    history.solid().validate().unwrap();

    let redone = history.redo().unwrap().unwrap();
    assert_eq!(redone.direction, HistoryDirection::Redo);
    assert_eq!((redone.revision_before, redone.revision_after), (2, 3));
    assert_eq!((history.solid().vertices.len(),history.solid().faces.len()), (9, 7));
    assert_eq!(history.solid().mass, before_mass);
    assert_eq!(history.solid().topology_name(TopologyEntity::Vertex(VertexId(8))), Some(&new_name));
    assert_eq!(history.solid().resolve_topology_handle(new_handle), Err(HandleError::StaleHandle));
    let reincarnated = history.solid().topology_handle(TopologyEntity::Vertex(VertexId(8))).unwrap();
    assert_ne!(new_handle, reincarnated);
    history.solid().validate().unwrap();
}

#[test]
fn separate_mementos_support_multiple_undo_redo_cycles() {
    let mut h = CommandHistory::new(cube()).unwrap();
    let edge = edge_command(h.solid());
    let face = face_command(h.solid());
    h.execute(CommandBatch::single(edge)).unwrap();
    h.execute(CommandBatch::single(face)).unwrap();
    assert_eq!((h.solid().revision, h.undo_depth()), (2, 2));
    h.undo().unwrap();
    assert_eq!((h.solid().vertices.len(), h.solid().faces.len()), (9, 6));
    h.undo().unwrap();
    assert_eq!((h.solid().vertices.len(), h.solid().faces.len()), (8, 6));
    assert!(h.undo().unwrap().is_none());
    h.redo().unwrap();
    h.redo().unwrap();
    assert!(h.redo().unwrap().is_none());
    assert_eq!((h.solid().vertices.len(), h.solid().faces.len()), (9, 7));
    assert_eq!(h.solid().revision, 6);
    h.solid().validate().unwrap();
}

#[test]
fn edited_branch_disables_redo_and_failed_command_keeps_history() {
    let mut h=CommandHistory::new(cube()).unwrap();
    let edge=edge_command(h.solid());
    h.execute(CommandBatch::single(edge.clone())).unwrap();
    h.undo().unwrap();
    assert_eq!((h.redo_depth(), h.undo_depth()), (1, 0));
    let failed = EulerCommand::SplitEdge { edge: TopologyName("missing".into()), fraction: 0.5 };
    assert!(matches!(h.execute(CommandBatch::single(failed)), Err(GeometryError::InvalidEdit(_))));
    assert_eq!((h.redo_depth(), h.undo_depth()), (1, 0));
    let revision=h.solid().revision;
    h.execute(CommandBatch::single(edge)).unwrap();
    assert_eq!(h.redo_depth(), 0);
    assert!(h.redo().unwrap().is_none());
    assert_eq!(h.solid().revision, revision+1);
}

#[test]
fn failed_second_command_rolls_back_entire_memento() {
    let mut h=CommandHistory::new(cube()).unwrap();
    let before=format!("{:?}", h.solid());
    let correct=edge_command(h.solid());
    let face_name=name(h.solid(),TopologyEntity::Face(FaceId(0)));
    let a=name(h.solid(),TopologyEntity::Vertex(VertexId(0)));
    let b=name(h.solid(),TopologyEntity::Vertex(VertexId(1)));
    let invalid=EulerCommand::SplitFace { face:face_name, start:a, end:b };
    let err=h.execute(CommandBatch { commands:vec![correct,invalid] }).unwrap_err();
    assert!(matches!(err,GeometryError::InvalidEdit(_)));
    assert_eq!(format!("{:?}", h.solid()), before);
    assert_eq!(h.undo_depth(), 0);
}

#[test]
fn recorded_commands_serialize_and_replay_on_regenerated_primitive() {
    let original=cube();
    let mut h=CommandHistory::new(original.clone()).unwrap();
    h.execute(CommandBatch::single(edge_command(&original))).unwrap();
    h.execute(CommandBatch::single(face_command(&original))).unwrap();
    let log=h.command_log();
    let json=serde_json::to_string(&log).unwrap();
    assert!(json.contains("split_edge"));
    let restored:Vec<CommandBatch>=serde_json::from_str(&json).unwrap();
    let replayed=CommandHistory::replay(cube(), &restored).unwrap();
    assert_eq!(h.solid().mesh.triangles,replayed.solid().mesh.triangles);
    assert_eq!(h.solid().mesh.triangle_faces,replayed.solid().mesh.triangle_faces);
    assert_eq!(h.solid().mass,replayed.solid().mass);
    let entity=TopologyEntity::Vertex(VertexId(8));
    assert_eq!(h.solid().topology_name(entity),replayed.solid().topology_name(entity));
    assert_ne!(h.solid().topology_handle(entity), replayed.solid().topology_handle(entity));
}

#[test]
fn memento_bound_preserves_full_command_log_and_prunes_old_undo() {
    let base=cube();
    let mut h=CommandHistory::with_undo_limit(base.clone(),1).unwrap();
    h.execute(CommandBatch::single(edge_command(&base))).unwrap();
    h.execute(CommandBatch::single(face_command(&base))).unwrap();
    assert_eq!(h.command_log().len(),2);
    assert_eq!((h.archived_depth(),h.undo_depth()),(1,1));
    h.undo().unwrap();
    assert_eq!((h.archived_depth(),h.undo_depth(),h.redo_depth()),(1,0,1));
    assert_eq!(h.command_log().len(),1);
    assert_eq!(h.solid().vertices.len(),9);
    assert!(h.undo().unwrap().is_none());
    let reconstructed=CommandHistory::replay(cube(),&h.command_log()).unwrap();
    assert_eq!(reconstructed.solid().mesh.triangles, h.solid().mesh.triangles);
}

#[test]
fn typed_names_prevent_mixing_edges_and_vertices() {
    let mut h=CommandHistory::new(cube()).unwrap();
    let wrong=name(h.solid(),TopologyEntity::Vertex(VertexId(0)));
    assert!(h.execute(CommandBatch::single(EulerCommand::SplitEdge { edge:wrong, fraction:0.5 })).is_err());
    assert_eq!(h.solid().revision,0);
    assert!(CommandHistory::with_undo_limit(cube(),0).is_err());
    assert!(h.execute(CommandBatch{commands:vec![]}).is_err());
}
