use rustsolid::{CommandBatch,CommandHistory,EdgeId,EulerCommand,FaceId,GeometryError,Point3,Solid,TopologyEntity,TopologyName,VertexId};
fn cube()->Solid{Solid::block(Point3{x:0.0,y:0.0,z:0.0},2.0,3.0,4.0).unwrap().with_feature_key("part-1/feat-0").unwrap()}
fn name(h:&CommandHistory,entity:TopologyEntity)->TopologyName{h.solid().topology_name(entity).unwrap().clone()}
fn split_edge(h:&CommandHistory,edge:EdgeId)->EulerCommand{EulerCommand::MakeEdgeVertex {edge:name(h,TopologyEntity::Edge(edge)),fraction:0.5}}
fn split_face(h:&CommandHistory,face:FaceId)->EulerCommand{
    let ring=&h.solid().faces[face.0 as usize].boundary;
    EulerCommand::MakeEdgeFace{face:name(h,TopologyEntity::Face(face)),start:name(h,TopologyEntity::Vertex(ring[0])),end:name(h,TopologyEntity::Vertex(ring[2]))}
}

#[test]
fn earlier_vertex_can_be_killed_and_undone_independently(){
    let mut h=CommandHistory::new(cube()).unwrap();
    h.execute(CommandBatch::single(split_edge(&h,EdgeId(0)))).unwrap();
    let first_name=name(&h,TopologyEntity::Vertex(VertexId(8)));
    h.execute(CommandBatch::single(split_edge(&h,EdgeId(5)))).unwrap();
    assert_eq!(h.solid().vertices.len(),10);
    h.execute(CommandBatch::single(EulerCommand::KillEdgeVertex {vertex:first_name})).unwrap();
    assert_eq!(h.solid().vertices.len(),9);
    assert_eq!((h.solid().revision,h.undo_depth()),(3,3));
    assert_eq!(h.solid().mass.volume,24.0);
    h.solid().validate().unwrap();
    h.undo().unwrap();
    assert_eq!((h.solid().vertices.len(),h.solid().revision),(10,4));
    h.redo().unwrap();
    assert_eq!((h.solid().vertices.len(),h.solid().revision),(9,5));
    let log=h.command_log();
    assert_eq!(log.len(),3);
    let replay=CommandHistory::replay(cube(),&log).unwrap();
    assert_eq!(replay.solid().mesh.triangles,h.solid().mesh.triangles);
    assert_eq!(replay.solid().mass,h.solid().mass);
}

#[test]
fn non_latest_face_split_is_removable_and_replayable(){
    let mut h=CommandHistory::new(cube()).unwrap();
    h.execute(CommandBatch::single(split_face(&h,FaceId(0)))).unwrap();
    let created_face=name(&h,TopologyEntity::Face(FaceId(6)));
    let diagonal=name(&h,TopologyEntity::Edge(EdgeId(12)));
    h.execute(CommandBatch::single(split_face(&h,FaceId(1)))).unwrap();
    h.execute(CommandBatch::single(EulerCommand::KillEdgeFace {edge:diagonal, removed_face:created_face})).unwrap();
    assert_eq!(h.solid().faces.len(),7);
    h.solid().validate().unwrap();
    let history=serde_json::to_string(&h.command_log()).unwrap();
    let commands:Vec<CommandBatch>=serde_json::from_str(&history).unwrap();
    let replay=CommandHistory::replay(cube(),&commands).unwrap();
    assert_eq!(replay.solid().faces.len(),7);
    h.undo().unwrap();
    assert_eq!(h.solid().faces.len(),8);
    h.redo().unwrap();
    assert_eq!(h.solid().faces.len(),7);
}

#[test]
fn mixed_batch_with_failed_kill_is_atomic(){
    let mut h=CommandHistory::new(cube()).unwrap();
    let original=format!("{:?}",h.solid());
    let command=split_edge(&h,EdgeId(0));
    let bad=EulerCommand::KillEdgeVertex {vertex:name(&h,TopologyEntity::Vertex(VertexId(0)))};
    let err=h.execute(CommandBatch{commands:vec![command,bad]}).unwrap_err();
    assert!(matches!(err,GeometryError::InvalidEdit(_)));
    assert_eq!(format!("{:?}",h.solid()),original);
    assert!(h.command_log().is_empty());
}

#[test]
fn invalid_kill_does_not_clear_redo_branch(){
    let mut h=CommandHistory::new(cube()).unwrap();
    h.execute(CommandBatch::single(split_edge(&h,EdgeId(0)))).unwrap();
    h.undo().unwrap();
    let err=h.execute(CommandBatch::single(EulerCommand::KillEdgeVertex{
        vertex:TopologyName("unknown".into())
    })).unwrap_err();
    assert!(matches!(err,GeometryError::InvalidEdit(_)));
    assert_eq!(h.redo_depth(),1);
}
