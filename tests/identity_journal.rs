use rustsolid::{
    EdgeId, FaceId, GeometryError, HandleError, Point3, Solid, TopologyEntity,
    TopologyName, VertexId,
};

fn cube() -> Solid {
    Solid::block(Point3 { x: 0.0, y: 0.0, z: 0.0 }, 2.0, 3.0, 4.0).unwrap()
}

#[test]
fn primitive_names_are_deterministic_per_feature_key_but_handles_are_scoped() {
    let a = cube().with_feature_key("part-42/extrude-7").unwrap();
    let b = Solid::block(Point3 { x: 8.0, y: 3.0, z: -5.0 }, 2.0, 3.0, 4.0)
        .unwrap().with_feature_key("part-42/extrude-7").unwrap();
    let face = TopologyEntity::Face(FaceId(0));
    let name = a.topology_name(face).unwrap().clone();
    assert_eq!(name.as_str(), "part-42/extrude-7/face/bottom");
    assert_eq!(b.topology_name(face), Some(&name));
    let a_handle = a.resolve_topology_name(&name).unwrap();
    assert_eq!(a.resolve_topology_handle(a_handle), Ok(face));
    assert_eq!(b.resolve_topology_handle(a_handle), Err(HandleError::ForeignBody));
    assert_ne!(a_handle, b.resolve_topology_name(&name).unwrap());
    assert_eq!(a.topology_parent_name(face), None);
}

#[test]
fn feature_key_rejects_invalid_namespaces_and_changes_after_edit() {
    assert!(cube().with_feature_key("").is_err());
    assert!(cube().with_feature_key("p?x").is_err());
    let mut body = cube().with_feature_key("component/feature").unwrap();
    body.edit_atomic(|tx| { tx.split_edge(EdgeId(0), 0.5)?; Ok(()) }).unwrap();
    let err = body.with_feature_key("changed").unwrap_err();
    assert!(matches!(err, GeometryError::InvalidEdit(_)));
}

#[test]
fn existing_handles_and_names_survive_local_splits() {
    let mut body = cube().with_feature_key("component/extrude").unwrap();
    let edge = TopologyEntity::Edge(EdgeId(0));
    let face = TopologyEntity::Face(FaceId(0));
    let old_edge_handle = body.topology_handle(edge).unwrap();
    let old_face_handle = body.topology_handle(face).unwrap();
    let edge_name = body.topology_name(edge).unwrap().clone();
    let face_name = body.topology_name(face).unwrap().clone();
    let a = body.faces[0].boundary[0];
    let c = body.faces[0].boundary[2];
    let report = body.edit_atomic(|tx| {
        tx.split_edge(EdgeId(0), 0.5)?;
        tx.split_face(FaceId(0), a, c)?;
        Ok(())
    }).unwrap();
    assert_eq!(report.revision_before, 0);
    assert_eq!(report.revision_after, 1);
    assert_eq!(report.named_changes.len(), 9); // 4 from split edge, 5 from split face
    assert_eq!(body.resolve_topology_handle(old_edge_handle), Ok(edge));
    assert_eq!(body.resolve_topology_handle(old_face_handle), Ok(face));
    assert_eq!(body.topology_name(edge), Some(&edge_name));
    assert_eq!(body.topology_name(face), Some(&face_name));
    for change in &report.named_changes {
        assert_eq!(body.resolve_topology_handle(change.handle), Ok(change.entity));
        assert_eq!(body.topology_name(change.entity), Some(&change.name));
        assert_eq!(body.topology_parent_name(change.entity), Some(&change.parent));
        assert_eq!(body.resolve_topology_name(&change.name), Some(change.handle));
    }
    body.validate().unwrap();
}

#[test]
fn newly_generated_names_replay_under_same_feature_history() {
    fn named_split() -> (TopologyName, TopologyName) {
        let mut body = cube().with_feature_key("extrude/1").unwrap();
        let report = body.edit_atomic(|tx| {
            tx.split_edge(EdgeId(0), 0.5)?;
            Ok(())
        }).unwrap();
        (report.named_changes[0].name.clone(), report.named_changes[1].name.clone())
    }
    assert_eq!(named_split(), named_split());
}

#[test]
fn handles_from_aborted_previews_cannot_alias_new_allocations() {
    let mut body = cube();
    let stale;
    {
        let mut tx = body.begin_edit().unwrap();
        tx.split_edge(EdgeId(0), 0.5).unwrap();
        stale = tx.preview().topology_handle(TopologyEntity::Vertex(VertexId(8))).unwrap();
    } // implicit inverse-journal rollback
    assert_eq!(body.vertices.len(), 8);
    assert_eq!(body.resolve_topology_handle(stale), Err(HandleError::StaleHandle));
    body.edit_atomic(|tx| { tx.split_edge(EdgeId(0), 0.5)?; Ok(()) }).unwrap();
    let fresh = body.topology_handle(TopologyEntity::Vertex(VertexId(8))).unwrap();
    assert_ne!(stale.generation, fresh.generation);
    assert_eq!(body.resolve_topology_handle(stale), Err(HandleError::StaleHandle));
    assert_eq!(body.resolve_topology_handle(fresh), Ok(TopologyEntity::Vertex(VertexId(8))));
    body.validate().unwrap();
}

#[test]
fn journals_snapshot_local_changes_not_the_whole_solid() {
    let mut body = cube();
    let mut tx = body.begin_edit().unwrap();
    tx.split_edge(EdgeId(0), 0.5).unwrap();
    let stats = tx.journal_stats();
    assert_eq!(stats.frames, 1);
    assert_eq!(stats.triangle_snapshots, 2);
    assert!(stats.topology_snapshots < tx.preview().coedges.len());
    let report = tx.commit().unwrap();
    assert_eq!(report.journal, stats);
    assert_eq!(report.journal.triangle_snapshots, 2);
    body.validate().unwrap();
}

#[test]
fn face_split_journals_only_owning_triangles() {
    let mut body = cube();
    let a = body.faces[0].boundary[0];
    let c = body.faces[0].boundary[2];
    let mut tx = body.begin_edit().unwrap();
    tx.split_face(FaceId(0), a, c).unwrap();
    assert_eq!(tx.journal_stats().triangle_snapshots, 2);
    assert_eq!(tx.preview().mesh.triangles.len(), 12);
    tx.commit().unwrap();
    body.validate().unwrap();
}

#[test]
fn atomic_failed_batch_rewinds_data_and_naming_sequence() {
    let mut body = cube().with_feature_key("part/extrude").unwrap();
    let original = format!("{body:?}");
    let (a,b) = (body.faces[0].boundary[0], body.faces[0].boundary[1]);
    let stale;
    {
        let mut tx = body.begin_edit().unwrap();
        tx.split_edge(EdgeId(0), 0.5).unwrap();
        stale = tx.preview().topology_handle(TopologyEntity::Vertex(VertexId(8))).unwrap();
        assert!(tx.split_face(FaceId(0), a, b).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(format!("{body:?}"), original);
    assert_eq!(body.resolve_topology_handle(stale), Err(HandleError::StaleHandle));
    let result = body.edit_atomic(|tx| { tx.split_edge(EdgeId(0), 0.5)?; Ok(()) }).unwrap();
    assert!(result.named_changes[0].name.as_str().contains("/edit-1/"));
    body.validate().unwrap();
}

#[test]
fn no_op_transaction_keeps_handles_revision_and_zero_journal() {
    let mut body = cube();
    let handle = body.topology_handle(TopologyEntity::Vertex(VertexId(0))).unwrap();
    let report = body.edit_atomic(|_| Ok(())).unwrap();
    assert_eq!(report.revision_after, 0);
    assert!(report.named_changes.is_empty());
    assert_eq!(report.journal.frames, 0);
    assert_eq!(body.resolve_topology_handle(handle), Ok(TopologyEntity::Vertex(VertexId(0))));
}
