//! Strong-atomic source-specific executable through/blind crossing.
//! It does NOT unlock general TopologyDeltaPlan or unrestricted Boolean.
use super::*;
use rustsolid::{ThroughBlindCrossingBrep,MixedFootprintContact};

#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum CrossingSourceEditInput {
    ReplaceThrough {center:[f64;2],radius:f64},
    ReplaceBlind {center:[f64;2],radius:f64,depth:f64},
}

pub(super) fn evaluate(
    source:SolidSource,
    through:mixed_regions_cli::ThroughHoleInput,
    blind:BlindHoleInput,
    batches:Vec<Vec<CrossingSourceEditInput>>,
    undo:usize,
    redo:usize,
    display_segments:Option<usize>,
)->Response {
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if batches.len()>64||undo>64||redo>64{
            return Err(GeometryError::InvalidEdit(
                "crossing source request exceeds 64 batch/history steps".into()));
        }
        let block=source.build()?;
        let t=HoleParameters{
            center_x:through.center[0],center_z:through.center[1],
            radius:through.radius,
        };
        let b=BlindHoleParameters{
            center_x:blind.center[0],center_z:blind.center[1],
            radius:blind.radius,depth:blind.depth,
        };
        let mut body=ThroughBlindCrossingBrep::from_box(&block,t,b)?;
        let mut commits=Vec::new();
        for batch in batches {
            if batch.is_empty()||batch.len()>64{
                return Err(GeometryError::InvalidEdit(
                    "crossing transaction must contain 1..=64 source edits".into()));
            }
            let count=batch.len();
            let mut tx=body.begin_edit();
            for command in batch{
                match command{
                    CrossingSourceEditInput::ReplaceThrough{center,radius}=>{
                        tx.replace_through(HoleParameters{
                            center_x:center[0],center_z:center[1],radius})?;
                    }
                    CrossingSourceEditInput::ReplaceBlind{center,radius,depth}=>{
                        tx.replace_blind(BlindHoleParameters{
                            center_x:center[0],center_z:center[1],radius,depth})?;
                    }
                }
            }
            let revision=tx.commit()?;
            commits.push(json!({"revision":revision,"command_count":count}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo{
            if !body.undo()? {
                return Err(GeometryError::InvalidEdit(
                    "crossing Undo exceeds committed history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":body.revision()}));
        }
        for _ in 0..redo{
            if !body.redo()? {
                return Err(GeometryError::InvalidEdit(
                    "crossing Redo exceeds available history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":body.revision()}));
        }
        body.validate()?;
        let m=body.model();
        let stats=m.summary();
        let plan=body.overlap();
        if plan.contact!=MixedFootprintContact::TwoTransverseCrossings {
            return Err(GeometryError::InvalidEdit(
                "crossing executor lost strict source crossing".into()));
        }
        let display=display_segments.map(|segments|body.tessellate(segments)).transpose()?;
        let quality=display.as_ref().map(|d|json!({
            "segments_per_arc":d.segments_per_arc,
            "top_triangles":d.top_triangles,
            "bottom_triangles":d.bottom_triangles,
            "crescent_triangles":d.crescent_triangles,
            "faceted_volume":d.faceted_volume,
            "faceted_surface_area":d.faceted_surface_area,
            "analytic_volume":d.analytic_volume,
            "faceting_volume_error":d.faceting_volume_error,
            "closed_oriented_manifold_certified":true,
        }));
        let solid=display.map(|d|SolidResponse{
            vertices:d.mesh.vertices.iter().map(|v|v.array()).collect(),
            faces:d.mesh.triangles,
            triangle_face_ids:d.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{
                vertices:stats.vertices,edges:stats.edges,
                coedges:stats.coedges,loops:stats.loops,
                shells:stats.shells,faces:stats.faces,
                euler_characteristic:stats.euler_characteristic,
            },
            mass:Mass{volume:m.mass.volume,surface_area:m.mass.surface_area,
                centroid:m.mass.centroid.array()},
            bbox:[m.bbox.min.array(),m.bbox.max.array()],
            edit_report:None,history:None,brep:None,
            analytic_brep:Some(json!({
                "representation":"analytic_through_blind_crossing",
                "genus":1,"cylindrical_faces":stats.cylindrical_faces,
                "circular_edges":stats.circular_edges,
                "mesh_is_display_only":true,
                "crescent_floor_face_id":8,
                "lower_through_wall_face_id":9,
            })),
        });
        let data=json!({
            "kind":"through_blind_crossing",
            "state":"committed",
            "representation":"analytic_through_blind_crossing",
            "topology_committed":true,
            "generalized_topology_delta_can_apply":false,
            "overlap_area":plan.overlap_area,
            "crossing_points":plan.crossing_points.map(|points|
                points.map(Point3::array)),
            "topology":{
                "vertices":stats.vertices,"edges":stats.edges,
                "coedges":stats.coedges,"loops":stats.loops,
                "faces":stats.faces,"shells":stats.shells,
                "genus":stats.genus,"euler_characteristic":stats.euler_characteristic,
                "inner_loops":m.loops.iter().filter(|l|
                    l.role==rustsolid::LoopRole::Inner).count(),
                "cylindrical_faces":stats.cylindrical_faces,
                "circular_edges":stats.circular_edges,
            },
            "mass":{"volume":m.mass.volume,"surface_area":m.mass.surface_area,
                "centroid":m.mass.centroid.array()},
            "revision":body.revision(),
            "undo_depth":body.undo_depth(),
            "redo_depth":body.redo_depth(),
            "commits":commits,"transitions":transitions,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality,
        });
        Ok((data,solid))
    })();
    match result {
        Ok((data,solid))=>Response{
            schema_version:"geometry.v1",ok:true,solid,error:None,data:Some(data)
        },
        Err(error)=>Response{
            schema_version:"geometry.v1",ok:false,
            solid:None,data:None,error:Some(error.to_string())
        },
    }
}
