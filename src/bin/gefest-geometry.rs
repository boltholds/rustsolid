//! Stateless one-request JSON stdio adapter. The "edit_solid" request runs an
//! entire batch atomically; source solids cannot be externally mutated halfway.
use rustsolid::{CommandBatch, CommandHistory, CylindricalBrep, EditDelta, EditReport, EdgeId, FaceId, GeometryError, GeometryTolerance, Point2, Point3, Solid, VertexId};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::io::{self, Read};

#[derive(Debug, Deserialize)]
#[serde(tag = "operation", rename_all = "snake_case")]
enum GeometryRequest {
    ExtrudeProfile { profile: Vec<[f64; 2]>, height: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    Block { origin: [f64; 3], width: f64, height: f64, depth: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    AnalyticCylinder {
        origin: [f64;3], radius: f64, height: f64,
        #[serde(default)] segments: Option<usize>,
        #[serde(default)] tolerance: Option<GeometryTolerance>,
    },
    EditSolid { source: SolidSource, edits: Vec<EditRequest> },
    CommandHistory { source: SolidSource, #[serde(default)] feature_key: Option<String>,
        batches: Vec<CommandBatch>, #[serde(default)] undo: usize, #[serde(default)] redo: usize },
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
    #[serde(skip_serializing_if = "Option::is_none")]
    history: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    analytic_brep: Option<Value>,
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
        EditDelta::KillEdgeVertex(k)=>json!({"kind":"kill_edge_vertex",
            "removed_vertex_id":k.removed_vertex.0,"removed_edge_id":k.removed_edge.0,
            "retained_edge_id":k.retained_edge.0}),
        EditDelta::KillEdgeFace(k)=>json!({"kind":"kill_edge_face",
            "removed_edge_id":k.removed_edge.0,"removed_face_id":k.removed_face.0,
            "retained_face_id":k.retained_face.0}),
    }).collect();
    let names = report.named_changes.into_iter().map(|item| json!({
        "entity": format!("{:?}",item.entity),"name":item.name.0,
        "parent_name":item.parent.0
    })).collect::<Vec<_>>();
    json!({"revision_before":report.revision_before,"revision_after":report.revision_after,
        "changes":changes,"name_changes":names,
        "journal":{"frames":report.journal.frames,
            "topology_snapshots":report.journal.topology_snapshots,
            "triangle_snapshots":report.journal.triangle_snapshots}})
}

/// Native periodic-cylinder B-rep returned alongside a display-only mesh.
/// The caller can distinguish true curved geometry through the optional
/// solid.analytic_brep metadata.
fn process_analytic_cylinder(
    origin: [f64; 3], radius: f64, height: f64,
    segments: usize, tolerance: GeometryTolerance,
) -> Response {
    let center = Point3 {x:origin[0],y:origin[1],z:origin[2]};
    let computed = CylindricalBrep::upright(center,radius,height,tolerance)
        .and_then(|body| body.tessellate(segments).map(|mesh|(body,mesh)));
    match computed {
        Err(error) => Response {schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string())},
        Ok((body,mesh)) => {
            let topo = Topology {
                vertices:body.vertices.len(),edges:body.edges.len(),
                coedges:body.coedges.len(),loops:body.loops.len(),
                shells:body.shells.len(),faces:body.faces.len(),
                euler_characteristic:body.euler_characteristic(),
            };
            let seam = body.seam_edge().0;
            let coedges=body.coedges.iter().map(|c| json!({
                "id":c.id.0,"edge_id":c.edge.0,"face_id":c.face.0,
                "loop_id":c.loop_id.0,"twin_id":c.twin.0,
                "next_id":c.next.0,"prev_id":c.prev.0,
                "reversed":c.reversed,"pcurve_id":body.coedge_geometry[c.id.0 as usize].pcurve.0,
            })).collect::<Vec<_>>();
            let analytic=json!({
                "kind":"right_circular_cylinder",
                "surface_carriers":["plane","plane","cylinder"],
                "curve_carriers":["circle","circle","line"],
                "face_surface_ids":[0,1,2],
                "edge_curve_ids":[0,1,2],
                "side_face_id":2,
                "cap_face_ids":[0,1],
                "circular_edge_ids":[0,1],
                "seam_edge_id":seam,
                "seam_coedge_ids":[4,5],
                "u_periodic":true,
                "u_range":[0.0,std::f64::consts::TAU],
                "v_range":[0.0,height],
                "coedges":coedges,
                "facet_segments":segments,
                "mass_is_analytic":true,
            });
            Response {
                schema_version:"geometry.v1",ok:true,
                solid:Some(SolidResponse {
                    vertices:mesh.vertices.iter().map(|v|v.array()).collect(),
                    faces:mesh.triangles,
                    triangle_face_ids:mesh.triangle_faces.iter().map(|f|f.0).collect(),
                    topology:topo,
                    mass:Mass {volume:body.mass.volume,
                        surface_area:body.mass.surface_area,
                        centroid:body.mass.centroid.array()},
                    bbox:[body.bbox.min.array(),body.bbox.max.array()],
                    edit_report:None,history:None,analytic_brep:Some(analytic),
                }),
                error:None,
            }
        }
    }
}

fn process(input: &str) -> Response {
    let parsed: Result<GeometryRequest, _> = serde_json::from_str(input);
    let result: Result<(Solid, Option<EditReport>, Option<Value>), GeometryError> = match parsed {
        Ok(GeometryRequest::ExtrudeProfile { profile, height, tolerance }) => {
            SolidSource::ExtrudeProfile { profile, height, tolerance }.build().map(|body|(body,None,None))
        }
        Ok(GeometryRequest::Block { origin, width, height, depth, tolerance }) => {
            SolidSource::Block { origin, width, height, depth, tolerance }.build().map(|body|(body,None,None))
        }
        Ok(GeometryRequest::AnalyticCylinder {origin,radius,height,segments,tolerance}) => {
            return process_analytic_cylinder(origin,radius,height,
                segments.unwrap_or(64),tolerance.unwrap_or_default());
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
                    Ok((body,Some(report),None))
                })
            }
        }
        Ok(GeometryRequest::CommandHistory { source, feature_key, batches, undo, redo }) => {
            if batches.len() > 256 || undo > 256 || redo > 256 {
                Err(GeometryError::InvalidEdit("history request exceeds 256-item limit".into()))
            } else {
                source.build().and_then(|solid| {
                    let solid = match feature_key {
                        Some(key) => solid.with_feature_key(&key)?,
                        None => solid,
                    };
                    let mut session = CommandHistory::new(solid)?;
                    for batch in batches { session.execute(batch)?; }
                    for _ in 0..undo {
                        if session.undo()?.is_none() {
                            return Err(GeometryError::InvalidEdit("undo depth exceeds history".into()));
                        }
                    }
                    for _ in 0..redo {
                        if session.redo()?.is_none() {
                            return Err(GeometryError::InvalidEdit("redo depth exceeds history".into()));
                        }
                    }
                    let state = json!({
                        "revision": session.solid().revision,
                        "undo_depth": session.undo_depth(), "redo_depth": session.redo_depth(),
                        "archived_batches": session.archived_depth(),
                        "applied_batches": session.command_log(),
                    });
                    Ok((session.into_solid(), None, Some(state)))
                })
            }
        }
        Err(error) => return Response {
            schema_version: "geometry.v1", ok: false, solid: None,
            error: Some(format!("invalid JSON request: {error}")),
        },
    };
    match result {
        Ok((body,report,history_state)) => Response {
            schema_version: "geometry.v1", ok: true,
            solid: Some(SolidResponse {
                vertices: body.mesh.vertices.iter().map(|v| v.array()).collect(),
                faces: body.mesh.triangles,
                triangle_face_ids: body.mesh.triangle_faces.iter().map(|id| id.0).collect(),
                topology: Topology {
                    vertices: body.vertices.len(), edges: body.edges.len(),
                    coedges: body.coedges.len(), loops: body.loops.len(), shells: body.shells.len(),
                    faces: body.faces.len(), euler_characteristic: body.vertices.len() as isize - body.edges.len() as isize + body.faces.len() as isize,
                },
                mass: Mass { volume: body.mass.volume, surface_area: body.mass.surface_area,
                    centroid: body.mass.centroid.array() },
                bbox: [body.bbox.min.array(), body.bbox.max.array()],
                edit_report: report.map(report_json),
                history: history_state,
                analytic_brep: None,
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
    #[test]
    fn json_command_history_undo_redo() {
        let response = process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "feature_key":"part/sample",
            "batches":[
              {"commands":[{"kind":"split_edge","edge":"part/sample/edge/0","fraction":0.5}]},
              {"commands":[{"kind":"split_face","face":"part/sample/face/bottom",
                "start":"part/sample/vertex/0","end":"part/sample/vertex/2"}]}
            ],
            "undo":2,"redo":1
        }"#);
        assert!(response.ok, "{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!((solid.topology.vertices,solid.topology.faces),(9,6));
        let history=solid.history.unwrap();
        assert_eq!(history["revision"],5);
        assert_eq!(history["undo_depth"],1);
        assert_eq!(history["redo_depth"],1);
        assert_eq!(history["applied_batches"].as_array().unwrap().len(),1);
    }

    #[test]
    fn json_command_history_rejects_excess_undo() {
        let response = process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[{"commands":[{"kind":"split_edge","edge":"primitive/edge/0","fraction":0.5}]}],
            "undo":2
        }"#);
        assert!(!response.ok);
        assert!(response.solid.is_none());
    }

    #[test]
    fn json_command_history_independent_kill_and_redo() {
        let result=process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "feature_key":"part/feature",
            "batches":[
                {"commands":[{"kind":"make_edge_vertex","edge":"part/feature/edge/0","fraction":0.5}]},
                {"commands":[{"kind":"kill_edge_vertex","vertex":"part/feature/edge/0/edit-1/vertex"}]}
            ],
            "undo":1, "redo":1
        }"#);
        assert!(result.ok,"{:?}",result.error);
        let solid=result.solid.unwrap();
        assert_eq!(solid.topology.vertices,8);
        assert_eq!(solid.topology.edges,12);
        assert_eq!(solid.mass.volume,24.0);
        let history=solid.history.unwrap();
        assert_eq!(history["revision"],4);
        assert_eq!(history["applied_batches"].as_array().unwrap().len(),2);
    }

    #[test]
    fn json_analytic_cylinder_exposes_exact_brep_and_display_mesh() {
        let result=process(r#"{"operation":"analytic_cylinder",
            "origin":[0,0,0],"radius":2,"height":5,"segments":32}"#);
        assert!(result.ok,"{:?}",result.error);
        let solid=result.solid.unwrap();
        assert_eq!((solid.topology.vertices,solid.topology.edges,solid.topology.coedges,
            solid.topology.loops,solid.topology.faces,solid.topology.shells),(2,3,6,3,3,1));
        assert_eq!(solid.topology.euler_characteristic,2);
        assert_eq!(solid.faces.len(),128);
        assert_eq!(solid.triangle_face_ids.len(),solid.faces.len());
        assert!((solid.mass.volume-std::f64::consts::PI*20.0).abs()<1e-10);
        let analytic=solid.analytic_brep.unwrap();
        assert_eq!(analytic["seam_edge_id"],2);
        assert_eq!(analytic["seam_coedge_ids"],json!([4,5]));
        assert_eq!(analytic["u_periodic"],true);
        assert_eq!(analytic["coedges"].as_array().unwrap().len(),6);
    }

    #[test]
    fn json_analytic_cylinder_rejects_invalid_parameters() {
        for request in [
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":0,"height":1}"#,
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":1,"height":2,"segments":3}"#,
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":1,"height":2,"tolerance":{"absolute_length":-1}}"#,
        ] {
            let response=process(request);
            assert!(!response.ok,"{request}");
            assert!(response.solid.is_none());
        }
    }

}
