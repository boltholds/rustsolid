use rustsolid::{BrepModel, CylindricalBrep, GeometryTolerance, Point3, QueryLimits,
    Solid, geometry_query};
use serde_json::json;

fn cylinder() -> BrepModel {
    let source=CylindricalBrep::upright(Point3{x:0.0,y:0.0,z:0.0},2.0,5.0,
        GeometryTolerance::default()).unwrap();
    BrepModel::from_cylinder(&source).unwrap()
}
fn box_brep() -> BrepModel {
    let block=Solid::block(Point3{x:1.0,y:2.0,z:3.0},2.0,3.0,4.0).unwrap();
    BrepModel::from_polyhedral(&block,rustsolid::ModelUnits::default()).unwrap()
}
fn run(model:&BrepModel,query:&str)->serde_json::Value {
    geometry_query(model,query,QueryLimits::default()).unwrap()
}

#[test]
fn simple_projection_returns_only_requested_fields() {
    let data=run(&cylinder(),"query { body { kind topology { faceCount seamEdges circularEdges genus } mass { volume centroid { y } } } }");
    assert_eq!(data["body"]["kind"],"analytic_cylinder");
    assert_eq!(data["body"]["topology"],json!({
        "faceCount":3,"seamEdges":1,"circularEdges":2,"genus":0}));
    assert_eq!(data["body"]["mass"]["centroid"],json!({"y":2.5}));
    assert_eq!(data["body"].as_object().unwrap().len(),3);
    assert!(data["body"]["mass"]["volume"].as_f64().unwrap()>60.0);
}

#[test]
fn face_and_edge_filters_and_nested_traversal() {
    let data=run(&cylinder(),r#"{
      body {
        faces(kind:"cylinder",limit:4) {id surface {kind radius} loops {id coedges {id edgeId pcurve {kind}}}}
        edges(seam:true) {id seam closed coedges {id faceId reversed}}
      }
    }"#);
    let faces=data["body"]["faces"].as_array().unwrap();
    assert_eq!(faces.len(),1);
    assert_eq!(faces[0]["id"],2);
    assert_eq!(faces[0]["surface"],json!({"kind":"cylinder","radius":2.0}));
    assert_eq!(faces[0]["loops"][0]["coedges"].as_array().unwrap().len(),4);
    let edges=data["body"]["edges"].as_array().unwrap();
    assert_eq!(edges.len(),1);
    assert_eq!(edges[0]["id"],2);
    assert_eq!(edges[0]["seam"],true);
    assert_eq!(edges[0]["coedges"].as_array().unwrap().len(),2);
}

#[test]
fn polygonal_block_query_has_no_curves_or_seams_and_supports_id_filter() {
    let data=run(&box_brep(),"{body { topology {vertexCount edgeCount faceCount} faces(id:1) {id surface {kind}} edges(limit:2) {id startId endId curve {kind}} }}");
    assert_eq!(data["body"]["topology"],json!({"vertexCount":8,"edgeCount":12,"faceCount":6}));
    assert_eq!(data["body"]["faces"],json!([{"id":1,"surface":{"kind":"plane"}}]));
    let edges=data["body"]["edges"].as_array().unwrap();
    assert_eq!(edges.len(),2);
    assert_eq!(edges[0]["curve"]["kind"],"line");
}

#[test]
fn unsupported_query_features_have_typed_errors() {
    let model=cylinder();
    for (source,expected) in [
        ("mutation { body { kind } }","query.syntax"),
        ("{body { invalid }}","query.unknown_field"),
        ("{body { faces }}","query.missing_selection"),
        ("{body { kind {x} }}","query.scalar_selection"),
        ("{body { faces(limit:500) {id} }}","query.limit"),
        ("{body { faces(kind:banana) {id} }}","query.invalid_argument"),
        ("{body { faces {id} faces {id} }}","query.field"),
        ("{body { edges(seam:123) {id} }}","query.invalid_argument"),
        ("{body { faces {id} }","query.syntax"),
    ] {
        let err=geometry_query(&model,source,QueryLimits::default()).unwrap_err();
        assert_eq!(err.code,expected,"source={source}: {err}");
    }
}

#[test]
fn bounded_query_rejects_excess_bytes_fields_nesting_and_visits() {
    let model=cylinder();
    let default=QueryLimits::default();
    let err=geometry_query(&model,"{body {kind}}",QueryLimits{max_bytes:6,..default}).unwrap_err();
    assert_eq!(err.code,"query.too_large");
    let err=geometry_query(&model,"{body { kind mass {volume} }}",QueryLimits{max_fields:2,..default}).unwrap_err();
    assert_eq!(err.code,"query.fields");
    let err=geometry_query(&model,"{body { edges {coedges {face {loops {id}}}}}}",
        QueryLimits{max_depth:3,..default}).unwrap_err();
    assert_eq!(err.code,"query.depth");
    let err=geometry_query(&model,"{body { edges(limit:10) { coedges { id } } }}",
        QueryLimits{max_visits:2,..default}).unwrap_err();
    assert_eq!(err.code,"query.cost");
}

#[test]
fn geometry_query_does_not_modify_any_geometric_carrier() {
    let model=cylinder();
    let origin_surface=model.geometry.surfaces.clone();
    let origin_curves=model.geometry.curves3.clone();
    for _ in 0..4 {
        run(&model,"{body { shells {id closed faces(limit:3) {id surface {kind}}} mass {volume}}}");
    }
    assert_eq!(model.geometry.surfaces,origin_surface);
    assert_eq!(model.geometry.curves3,origin_curves);
}
