use rustsolid::{EdgeId,FaceId,GeometryError,HandleError,Point3,Solid,TopologyEntity,VertexId};
fn cube()->Solid{Solid::block(Point3{x:0.0,y:0.0,z:0.0},2.0,3.0,4.0).unwrap()
    .with_feature_key("part/extrude").unwrap()}

#[test]
fn kev_can_remove_an_earlier_split_after_later_unrelated_split(){
    let mut body=cube();
    let first=body.edit_atomic(|tx|{tx.split_edge(EdgeId(0),0.5)?;Ok(())}).unwrap();
    let first_vertex=match first.changes[0]{rustsolid::EditDelta::SplitEdge(x)=>x.inserted_vertex,_=>panic!()};
    let second=body.edit_atomic(|tx|{tx.split_edge(EdgeId(5),0.5)?;Ok(())}).unwrap();
    let second_vertex=match second.changes[0]{rustsolid::EditDelta::SplitEdge(x)=>x.inserted_vertex,_=>panic!()};
    let second_name=body.topology_name(TopologyEntity::Vertex(second_vertex)).unwrap().clone();
    let second_handle=body.topology_handle(TopologyEntity::Vertex(second_vertex)).unwrap();
    let initial_volume=body.mass.volume;
    let before_revision=body.revision;
    let result=body.kill_edge_vertex(first_vertex).unwrap();
    assert_eq!(result.removed_vertex,first_vertex);
    assert_eq!(body.vertices.len(),9);
    assert_eq!(body.revision,before_revision+1);
    assert_eq!(body.mass.volume,initial_volume);
    body.validate().unwrap();
    let new_handle=body.resolve_topology_name(&second_name).unwrap();
    assert_ne!(new_handle,second_handle);
    assert_eq!(body.resolve_topology_handle(second_handle),Err(HandleError::StaleHandle));
    assert_eq!(body.topology_name(TopologyEntity::Vertex(VertexId(8))),Some(&second_name));
}

#[test]
fn kef_can_merge_an_earlier_face_after_an_unrelated_later_split(){
    let mut body=cube();
    let bottom=body.faces[0].boundary.clone();
    let top=body.faces[1].boundary.clone();
    let first=body.edit_atomic(|tx|{tx.split_face(FaceId(0),bottom[0],bottom[2])?;Ok(())}).unwrap();
    let first=match first.changes[0]{rustsolid::EditDelta::SplitFace(x)=>x,_=>panic!()};
    let second=body.edit_atomic(|tx|{tx.split_face(FaceId(1),top[0],top[2])?;Ok(())}).unwrap();
    let second=match second.changes[0]{rustsolid::EditDelta::SplitFace(x)=>x,_=>panic!()};
    let survivor_name=body.topology_name(TopologyEntity::Face(second.created_face)).unwrap().clone();
    let survivor_handle=body.topology_handle(TopologyEntity::Face(second.created_face)).unwrap();
    let result=body.kill_edge_face(first.diagonal_edge,first.created_face).unwrap();
    assert_eq!(result.removed_face,first.created_face);
    assert_eq!(body.faces.len(),7);
    assert_eq!(body.mass.volume,24.0);
    body.validate().unwrap();
    assert!(body.resolve_topology_name(&survivor_name).is_some());
    assert_eq!(body.resolve_topology_handle(survivor_handle),Err(HandleError::StaleHandle));
}

#[test]
fn kev_rejects_non_valence_two_and_never_mutates_on_error(){
    let mut body=cube();
    let before=format!("{body:?}");
    assert!(matches!(body.kill_edge_vertex(VertexId(0)),Err(GeometryError::InvalidEdit(_))));
    assert_eq!(format!("{body:?}"),before);
}
#[test]
fn kef_rejects_edges_between_non_coplanar_faces(){
    let mut body=cube();
    let before=format!("{body:?}");
    let edge=body.edges[0].clone();
    let removed=body.coedges[edge.coedges[0].0 as usize].face;
    assert!(body.kill_edge_face(edge.id,removed).is_err());
    assert_eq!(format!("{body:?}"),before);
}
#[test]
fn make_then_kill_restores_closed_wound_brep(){
    let mut body=cube();
    let old=body.mass;
    let split=body.edit_atomic(|tx|{tx.split_edge(EdgeId(0),0.25)?;Ok(())}).unwrap();
    let id=match split.changes[0]{rustsolid::EditDelta::SplitEdge(x)=>x.inserted_vertex,_=>panic!()};
    body.kill_edge_vertex(id).unwrap();
    assert_eq!((body.vertices.len(),body.edges.len(),body.coedges.len()),(8,12,24));
    assert_eq!(body.mass,old);
    let original_face=body.faces[0].boundary.clone();
    let split=body.edit_atomic(|tx|{tx.split_face(FaceId(0),original_face[0],original_face[2])?;Ok(())}).unwrap();
    let s=match split.changes[0]{rustsolid::EditDelta::SplitFace(x)=>x,_=>panic!()};
    body.kill_edge_face(s.diagonal_edge,s.created_face).unwrap();
    assert_eq!((body.faces.len(),body.edges.len(),body.loops.len()),(6,12,6));
    body.validate().unwrap();
}
