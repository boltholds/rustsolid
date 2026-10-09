//! Stateless one-request JSON stdio adapter. The "edit_solid" request runs an
//! entire batch atomically; source solids cannot be externally mutated halfway.
use rustsolid::{EditDelta, EditReport, EdgeId, FaceId, GeometryError, GeometryTolerance, Point2, Point3, Solid, VertexId};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::io::{self, Read};

#[derive(Debug, Deserialize)]
#[serde(tag = "operation", rename_all = "snake_case")]
enum GeometryRequest {
    ExtrudeProfile { profile: Vec<[f64; 2]>, height: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    Block { origin: [f64; 3], width: f64, height: f64, depth: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    EditSolid { source: SolidSource, edits: Vec<EditRequest> },
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum SolidSource {
    ExtrudeProfile { profile: Vec<[f64; 2]>, height: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    Block { origin: [f64; 3], width: f64, height: f64, depth: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
}
impl SolidSource {
    fn build(self) -> Result<Solid, GeometryError> {
        match self {
            Self::ExtrudeProfile { profile, height, tolerance } => {
                let points: Vec<Point2> = profile.iter().map(|p| Point2 { x: p[0], z: p[1] }).collect();
                Solid::extrude_xz_with_tolerance(&points, height, tolerance.unwrap_or_default())
            }
            Self::Block { origin, width, height, depth, tolerance } => Solid::block_with_tolerance(
                Point3 { x: origin[0], y: origin[1], z: origin[2] }, width, height, depth,
                tolerance.unwrap_or_default(),
            ),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum EditRequest {
    SplitEdge { edge_id: u32, fraction: f64 },
    SplitFace { face_id: u32, start_vertex_id: u32, end_vertex_id: u32 },
}

#[derive(Serialize)]
struct Topology {
    vertices: usize,
    edges: usize,
    coedges: usize,
    loops: usize,
    shells: usize,
    faces: usize,
    euler_characteristic: isize,
}

#[derive(Serialize)]
struct Mass {
    volume: f64,
    surface_area: f64,
    centroid: [f64; 3],
}

#[derive(Serialize)]
struct SolidResponse {
    vertices: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    triangle_face_ids: Vec<u32>,
    topology: Topology,
    mass: Mass,
    bbox: [[f64; 3]; 2],
    /// Included only for the new batch editing operation; geometry.v1 fields
    /// for legacy block and extrusion remain unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    edit_report: Option<Value>,
}

#[derive(Serialize)]
struct Response {
    schema_version: &'static str,
    ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    solid: Option<SolidResponse>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

fn report_json(report: EditReport) -> Value {
    let changes: Vec<Value> = report.changes.into_iter().map(|event| match event {
        EditDelta::SplitEdge(s) => json!({"kind":"split_edge", "original_edge_id":s.original_edge.0,
            "inserted_vertex_id":s.inserted_vertex.0,"created_edge_id":s.created_edge.0,
            "created_coedge_ids":[s.created_coedges[0].0,s.created_coedges[1].0]}),
        EditDelta::SplitFace(s) => json!({"kind":"split_face", "original_face_id":s.original_face.0,
            "created_face_id":s.created_face.0,"created_loop_id":s.created_loop.0,
            "diagonal_edge_id":s.diagonal_edge.0,
            "created_coedge_ids":[s.created_coedges[0].0,s.created_coedges[1].0]}),
    }).collect();
    json!({"revision_before":report.revision_before,"revision_after":report.revision_after,
        "changes":changes})
}

fn process(input: &str) -> Response {
    let parsed: Result<GeometryRequest, _> = serde_json::from_str(input);
    let result: Result<(Solid, Option<EditReport>), GeometryError> = match parsed {
        Ok(GeometryRequest::ExtrudeProfile { profile, height, tolerance }) => {
            SolidSource::ExtrudeProfile { profile, height, tolerance }.build().map(|body|(body,None))
        }
        Ok(GeometryRequest::Block { origin, width, height, depth, tolerance }) => {
            SolidSource::Block { origin, width, height, depth, tolerance }.build().map(|body|(body,None))
        }
        Ok(GeometryRequest::EditSolid { source, edits }) => {
            if edits.len() > 256 {
                Err(GeometryError::InvalidEdit("maximum 256 edits per request".into()))
            } else {
                source.build().and_then(|mut body| {
                    let report=body.edit_atomic(|tx| {
                        for edit in edits {
                            match edit {
                                EditRequest::SplitEdge { edge_id, fraction } => {
                                    tx.split_edge(EdgeId(edge_id),fraction)?;
                                }
                                EditRequest::SplitFace { face_id, start_vertex_id, end_vertex_id } => {
                                    tx.split_face(FaceId(face_id),VertexId(start_vertex_id),VertexId(end_vertex_id))?;
                                }
                            }
                        }
                        Ok(())
                    })?;
                    Ok((body,Some(report)))
                })
            }
        }
        Err(error) => return Response {
            schema_version: "geometry.v1", ok: false, solid: None,
            error: Some(format!("invalid JSON request: {error}")),
        },
    };
    match result {
        Ok((body,report)) => Response {
            schema_version: "geometry.v1", ok: true,
            solid: Some(SolidResponse {
                vertices: body.mesh.vertices.iter().map(|v| v.array()).collect(),
                faces: body.mesh.triangles,
                triangle_face_ids: body.mesh.triangle_faces.iter().map(|id| id.0).collect(),
                topology: Topology {
                    vertices: body.vertices.len(), edges: body.edges.len(),
                    coedges: body.coedges.len(), loops: body.loops.len(), shells: body.shells.len(),
                    faces: body.faces.len(), euler_characteristic: body.euler_characteristic(),
                },
                mass: Mass { volume: body.mass.volume, surface_area: body.mass.surface_area,
                    centroid: body.mass.centroid.array() },
                bbox: [body.bbox.min.array(), body.bbox.max.array()],
                edit_report: report.map(report_json),
            }),
            error: None,
        },
        Err(error) => Response { schema_version: "geometry.v1", ok: false,
            solid: None, error: Some(error.to_string()) },
    }
}

fn main() {
    let mut input = String::new();
    if let Err(error) = io::stdin().read_to_string(&mut input) {
        eprintln!("read failed: {error}");
        std::process::exit(1);
    }
    let response = process(&input);
    match serde_json::to_string(&response) {
        Ok(output) => println!("{output}"),
        Err(error) => { eprintln!("serialization failed: {error}"); std::process::exit(1); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_block_roundtrip() {
        let response = process(r#"{"operation":"block","origin":[1,2,3],"width":2,"height":3,"depth":4}"#);
        assert!(response.ok);
        let solid = response.solid.unwrap();
        assert_eq!(solid.topology.euler_characteristic, 2);
        assert_eq!(solid.topology.coedges, 24);
        assert_eq!(solid.topology.loops, 6);
        assert_eq!(solid.topology.shells, 1);
        assert_eq!(solid.mass.volume, 24.0);
        assert_eq!(solid.triangle_face_ids.len(), solid.faces.len());
        assert!(solid.edit_report.is_none());
    }

    #[test]
    fn malformed_request_reports_error() {
        assert!(!process(r#"{"operation":"boolean"}"#).ok);
        assert!(!process("garbage").ok);
    }

    #[test]
    fn per_request_tolerance_is_accepted_and_validated() {
        let request = r#"{"operation":"block","origin":[0,0,0],"width":0.00000001,"height":1,"depth":1,
            "tolerance":{"absolute_length":1e-12,"relative_length":1e-13}}"#;
        assert!(process(request).ok);
        let invalid = r#"{"operation":"block","origin":[0,0,0],"width":1,"height":1,"depth":1,
            "tolerance":{"absolute_length":-2}}"#;
        let response = process(invalid);
        assert!(!response.ok);
        assert!(response.error.unwrap().contains("invalid tolerance"));
    }

    #[test]
    fn json_edits_apply_atomically_and_return_provenance() {
        let response=process(r#"{"operation":"edit_solid",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "edits":[{"kind":"split_edge","edge_id":0,"fraction":0.5},
                     {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":2}]}"#);
        assert!(response.ok,"{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!(solid.topology.vertices,9);
        assert_eq!(solid.topology.faces,7);
        assert_eq!(solid.mass.volume,24.0);
        let report=solid.edit_report.unwrap();
        assert_eq!(report["revision_before"],0);
        assert_eq!(report["revision_after"],1);
        assert_eq!(report["changes"].as_array().unwrap().len(),2);
        assert_eq!(report["changes"][0]["kind"],"split_edge");
        assert_eq!(report["changes"][1]["kind"],"split_face");
    }

    #[test]
    fn json_failed_batch_does_not_return_partial_solid() {
        let response=process(r#"{"operation":"edit_solid",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "edits":[{"kind":"split_edge","edge_id":0,"fraction":0.5},
                     {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":1}]}"#);
        assert!(!response.ok);
        assert!(response.solid.is_none());
    }
}
