use rustsolid::{
    CoedgeId, EdgeId, EditDelta, FaceId, FaceRole, GeometryError, Point2, Point3,
    Solid, TopologyEntity, VertexId,
};

fn cube() -> Solid {
    Solid::block(Point3 { x: 3.0, y: -2.0, z: 7.0 }, 2.0, 3.0, 4.0).unwrap()
}
fn approx(a:f64,b:f64) {assert!((a-b).abs()<1e-8,"{a} vs {b}");}

#[test]
fn split_edge_retains_existing_handles_and_produces_valid_mesh() {
    let mut body=cube();
    let before_faces=body.faces.iter().map(|f|f.id).collect::<Vec<_>>();
    let old_edge=body.edges[0].clone();
    let old_id=old_edge.id;
    let old_coedges=old_edge.coedges;
    let before_mass=body.mass;
    let report=body.edit_atomic(|tx| {tx.split_edge(old_id,0.25)?;Ok(())}).unwrap();
    let change=match report.changes[0] { EditDelta::SplitEdge(v)=>v,_=>panic!("wrong event") };
    assert_eq!((report.revision_before,report.revision_after,body.revision),(0,1,1));
    assert_eq!(change.original_edge,old_id);
    assert_eq!(change.inserted_vertex,VertexId(8));
    assert_eq!(change.created_edge,EdgeId(12));
    assert_eq!(body.vertices.len(),9);
    assert_eq!(body.edges.len(),13);
    assert_eq!(body.coedges.len(),26);
    assert_eq!(body.mesh.triangles.len(),14);
    assert_eq!(body.faces.iter().map(|f|f.id).collect::<Vec<_>>(),before_faces);
    assert_eq!(body.edges[old_id.0 as usize].id,old_id);
    assert_eq!(body.coedges[old_coedges[0].0 as usize].id,old_coedges[0]);
    assert_eq!(body.coedges[old_coedges[1].0 as usize].id,old_coedges[1]);
    approx(body.mass.volume,before_mass.volume);
    approx(body.mass.surface_area,before_mass.surface_area);
    assert_eq!(body.mass.centroid,before_mass.centroid);
    assert_eq!(body.euler_characteristic(),2);
    body.validate().unwrap();
}

#[test]
fn split_face_retains_original_face_and_creates_diagonal() {
    let mut body=cube();
    let top=FaceId(1);
    let a=body.faces[1].boundary[0];
    let b=body.faces[1].boundary[2];
    let report=body.edit_atomic(|tx| {tx.split_face(top,a,b)?;Ok(())}).unwrap();
    let change=match report.changes[0] {EditDelta::SplitFace(v)=>v,_=>panic!("wrong event")};
    assert_eq!(change.original_face,top);
    assert_eq!(change.created_face,FaceId(6));
    assert_eq!(change.created_loop.0,6);
    assert_eq!(change.diagonal_edge,EdgeId(12));
    assert_eq!(body.faces[top.0 as usize].role,FaceRole::TopCap);
    assert_eq!(body.faces[change.created_face.0 as usize].role,FaceRole::TopCap);
    assert_eq!((body.vertices.len(),body.edges.len(),body.coedges.len(),body.faces.len(),body.loops.len()),(8,13,26,7,7));
    assert_eq!(body.mesh.triangles.len(),12);
    assert_eq!(body.euler_characteristic(),2);
    approx(body.mass.volume,24.0);
    approx(body.mass.surface_area,52.0);
    let new_faces=body.mesh.triangle_faces.iter().filter(|&&f|f==change.created_face).count();
    assert_eq!(new_faces,1);
    body.validate().unwrap();
}

#[test]
fn transaction_can_compose_edge_and_face_edits() {
    let mut body=cube();
    let a=body.faces[0].boundary[0];
    let c=body.faces[0].boundary[2];
    let result=body.edit_atomic(|tx| {
        let edge=tx.split_edge(EdgeId(0),0.5)?;
        assert_eq!(edge.inserted_vertex,VertexId(8));
        assert_eq!(tx.preview().revision,0);
        tx.split_face(FaceId(0),a,c)?;
        Ok(())
    }).unwrap();
    assert_eq!(result.changes.len(),2);
    assert_eq!(body.revision,1); // one atomic commit, not one per operation
    assert_eq!((body.vertices.len(),body.edges.len(),body.faces.len()),(9,14,7));
    approx(body.mass.volume,24.0);
    body.validate().unwrap();
}

#[test]
fn a_failed_second_edit_rolls_back_the_entire_batch() {
    let mut body=cube();
    let original=format!("{body:?}");
    let (a,b)=(body.faces[0].boundary[0],body.faces[0].boundary[1]);
    let error=body.edit_atomic(|tx| {
        tx.split_edge(EdgeId(0),0.5)?;
        tx.split_face(FaceId(0),a,b)?;
        Ok(())
    }).unwrap_err();
    assert!(matches!(error,GeometryError::InvalidEdit(_)));
    assert_eq!(format!("{body:?}"),original);
    assert_eq!(body.revision,0);
    body.validate().unwrap();
}

#[test]
fn dropping_a_transaction_aborts_and_failed_transactions_cannot_commit() {
    let mut body=cube();
    let original=format!("{body:?}");
    {
        let mut tx=body.begin_edit().unwrap();
        tx.split_edge(EdgeId(0),0.5).unwrap();
        assert_eq!(tx.preview().vertices.len(),9);
    }
    assert_eq!(format!("{body:?}"),original);
    let mut tx=body.begin_edit().unwrap();
    let err=tx.split_edge(EdgeId(0),1e-14).unwrap_err();
    assert!(matches!(err,GeometryError::InvalidEdit(_)));
    assert!(tx.commit().is_err());
    assert_eq!(format!("{body:?}"),original);
}

#[test]
fn face_split_rejects_outside_diagonal_for_concave_face() {
    let profile=[Point2{x:0.0,z:0.0},Point2{x:3.0,z:0.0},
        Point2{x:3.0,z:1.0},Point2{x:1.0,z:1.0},
        Point2{x:1.0,z:3.0},Point2{x:0.0,z:3.0}];
    let mut body=Solid::extrude_xz(&profile,2.0).unwrap();
    let original=format!("{body:?}");
    let (u,v)=(body.faces[0].boundary[2],body.faces[0].boundary[4]);
    let err=body.edit_atomic(|tx|{tx.split_face(FaceId(0),u,v)?;Ok(())}).unwrap_err();
    assert!(format!("{err}").contains("diagonal"));
    assert_eq!(format!("{body:?}"),original);
    body.validate().unwrap();
}

#[test]
fn invalid_or_unrepresentable_edge_split_rejects_without_mutation() {
    let mut body=cube();
    let original=format!("{body:?}");
    for t in [f64::NAN,f64::INFINITY,0.0,1.0,1e-14,1.0-1e-14] {
        assert!(body.edit_atomic(|tx|{tx.split_edge(EdgeId(0),t)?;Ok(())}).is_err());
    }
    assert_eq!(format!("{body:?}"),original);
}

#[test]
fn split_face_handles_orientation_reversal_and_multiple_sessions() {
    let mut body=cube();
    for id in [FaceId(0),FaceId(1)] {
        let ring=body.faces[id.0 as usize].boundary.clone();
        let (u,v)=(ring[2],ring[0]); // input order reversed intentionally
        let report=body.edit_atomic(|tx|{tx.split_face(id,u,v)?;Ok(())}).unwrap();
        assert_eq!(report.revision_after,report.revision_before+1);
        body.validate().unwrap();
    }
    assert_eq!((body.revision,body.faces.len(),body.edges.len()),(2,8,14));
}

#[test]
fn composite_ownership_does_not_follow_twin_or_edge_references() {
    let body=cube();
    let shells=body.topology_children(TopologyEntity::Body);
    assert_eq!(shells.len(),1);
    let faces=body.topology_children(shells[0]);
    assert_eq!(faces.len(),6);
    let loops=body.topology_children(faces[0]);
    assert_eq!(loops.len(),1);
    let coedges=body.topology_children(loops[0]);
    assert_eq!(coedges.len(),4);
    assert!(matches!(coedges[0],TopologyEntity::Coedge(CoedgeId(_))));
    assert!(body.topology_children(coedges[0]).is_empty());
    assert!(body.topology_children(TopologyEntity::Edge(EdgeId(0))).is_empty());
    assert!(body.topology_children(TopologyEntity::Vertex(VertexId(0))).is_empty());
}

#[test]
fn noop_transaction_does_not_advance_revision() {
    let mut body=cube();
    let report=body.edit_atomic(|_|Ok(())).unwrap();
    assert_eq!(report.revision_before,report.revision_after);
    assert_eq!(body.revision,0);
}

#[test]
fn concave_profile_accepts_interior_face_division() {
    let polygon=[Point2{x:0.0,z:0.0},Point2{x:3.0,z:0.0},
        Point2{x:3.0,z:1.0},Point2{x:1.0,z:1.0},
        Point2{x:1.0,z:3.0},Point2{x:0.0,z:3.0}];
    let mut body=Solid::extrude_xz(&polygon,2.0).unwrap();
    let (u,v)=(body.faces[0].boundary[0],body.faces[0].boundary[3]);
    let before=body.mass;
    let change=body.edit_atomic(|tx|{tx.split_face(FaceId(0),u,v)?;Ok(())}).unwrap();
    assert_eq!(change.changes.len(),1);
    approx(body.mass.volume,before.volume);
    approx(body.mass.surface_area,before.surface_area);
    body.validate().unwrap();
}

#[test]
fn far_from_origin_split_keeps_metrics_and_identity() {
    let mut body=Solid::block(Point3{x:1e12,y:1e12,z:-1e12},2.0,3.0,4.0).unwrap();
    let old=body.mass;
    let first=body.edges[0].id;
    let _=body.edit_atomic(|tx|{tx.split_edge(first,0.5)?;Ok(())}).unwrap();
    assert_eq!(body.mass.centroid,old.centroid);
    approx(body.mass.volume,old.volume);
    body.validate().unwrap();
}
