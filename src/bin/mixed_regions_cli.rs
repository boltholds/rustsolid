//! Request-scoped strongly atomic mixed through/blind source composition.
//! Never expose a partial B-rep or a facet-only substitute on failure.
use super::*;
use rustsolid::{MixedThroughBlindBrep,MixedThroughBlindParameters};

#[derive(Debug,serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ThroughHoleInput {center:[f64;2],radius:f64}

#[derive(Debug,serde::Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum MixedRegionEditInput {
    ReplaceThrough {index:usize,center:[f64;2],radius:f64},
    ReplaceBlind {index:usize,center:[f64;2],radius:f64,depth:f64},
}
pub(super) fn evaluate(
    source:SolidSource,through:Vec<ThroughHoleInput>,
    blind:Vec<BlindHoleInput>,batches:Vec<Vec<MixedRegionEditInput>>,
    undo:usize,redo:usize,display_segments:Option<usize>,
)->Response{
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if through.len()+blind.len()>16 || batches.len()>64 ||
            undo>64 || redo>64 {
            return Err(GeometryError::InvalidEdit(
                "mixed region request exceeds 16 source cylinders or 64 history/batch steps".into()));
        }
        let parameters=MixedThroughBlindParameters{
            through:through.into_iter().map(|p|HoleParameters{
                center_x:p.center[0],center_z:p.center[1],radius:p.radius
            }).collect(),
            blind:blind.into_iter().map(|p|BlindHoleParameters{
                center_x:p.center[0],center_z:p.center[1],
                radius:p.radius,depth:p.depth
            }).collect(),
        };
        let mut shape=MixedThroughBlindBrep::from_box(&source.build()?,parameters)?;
        let mut commits=Vec::new();
        for batch in batches {
            if batch.is_empty()||batch.len()>256 {
                return Err(GeometryError::InvalidEdit(
                    "mixed region edit batch requires 1..=256 replacement commands".into()));
            }
            let command_count=batch.len();
            let mut tx=shape.begin_edit();
            for cmd in batch{
                match cmd {
                    MixedRegionEditInput::ReplaceThrough{index,center,radius}=>{
                        tx.replace_through(index,HoleParameters{
                            center_x:center[0],center_z:center[1],radius})?;
                    }
                    MixedRegionEditInput::ReplaceBlind{index,center,radius,depth}=>{
                        tx.replace_blind(index,BlindHoleParameters{
                            center_x:center[0],center_z:center[1],radius,depth})?;
                    }
                }
            }
            let revision=tx.commit()?;
            commits.push(json!({"revision":revision,"command_count":command_count}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo{
            if !shape.undo()?{
                return Err(GeometryError::InvalidEdit(
                    "mixed-region Undo exceeds committed history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":shape.revision()}));
        }
        for _ in 0..redo{
            if !shape.redo()?{
                return Err(GeometryError::InvalidEdit(
                    "mixed-region Redo exceeds available history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":shape.revision()}));
        }
        shape.validate()?;
        let display=display_segments.map(|s|shape.tessellate(s)).transpose()?;
        let model=shape.model();let stats=model.summary();
        let quality=display.as_ref().map(|mesh|json!({
            "segments_per_source_circle":mesh.segments_per_circle,
            "top_face_triangles":mesh.top_face_triangles,
            "bottom_face_triangles":mesh.bottom_face_triangles,
            "top_face_loop_count":mesh.top_face_loop_count,
            "bottom_face_loop_count":mesh.bottom_face_loop_count,
            "floor_regions":mesh.floor_regions.iter().map(|f|json!({
                "component_index":f.component_index,
                "face_id":f.floor_face.0,
                "loop_vertices":f.loop_vertices,
                "floor_triangles":f.floor_triangles,
                "polygonal_area":f.polygonal_area,
                "polygonal_perimeter":f.polygonal_perimeter
            })).collect::<Vec<_>>(),
            "faceted_volume":mesh.faceted_volume,
            "faceted_surface_area":mesh.faceted_surface_area,
            "analytic_volume":mesh.analytic_volume,
            "faceting_volume_error":mesh.faceting_volume_error,
            "closed_oriented_manifold_certified":true
        }));
        let solid=display.map(|mesh|SolidResponse{
            vertices:mesh.mesh.vertices.iter().map(|v|v.array()).collect(),
            faces:mesh.mesh.triangles,
            triangle_face_ids:mesh.mesh.triangle_faces.iter()
                .map(|f|f.0).collect(),
            topology:Topology{
                vertices:stats.vertices,edges:stats.edges,
                coedges:stats.coedges,loops:stats.loops,
                shells:stats.shells,faces:stats.faces,
                euler_characteristic:stats.euler_characteristic
            },
            mass:Mass{volume:model.mass.volume,
                surface_area:model.mass.surface_area,
                centroid:model.mass.centroid.array()},
            bbox:[model.bbox.min.array(),model.bbox.max.array()],
            edit_report:None,history:None,
            analytic_brep:Some(json!({
                "representation":"analytic_mixed_through_blind_regions",
                "through_count":shape.through_holes().len(),
                "blind_regions":shape.components().len(),
                "cylindrical_faces":stats.cylindrical_faces,
                "circular_edges":stats.circular_edges,
                "genus":stats.genus,"mesh_is_display_only":true
            })),
            brep:None,
        });
        let data=json!({
            "kind":"mixed_through_blind_regions","state":"committed",
            "representation":"analytic_mixed_through_blind_regions",
            "topology_committed":true,
            "through_count":shape.through_holes().len(),
            "blind_count":shape.blind_holes().len(),
            "blind_regions":shape.components().iter().map(|c|json!({
                "kind":match c.kind{
                    rustsolid::BlindRegionKind::SingleCircle=>"single_circle",
                    rustsolid::BlindRegionKind::TwoCrossingCircles=>"two_crossing_circles"
                },
                "source_indices":c.source_indices,
                "opening_loop_id":c.opening_loop.0,
                "floor_face_id":c.floor_face.0,
                "cylindrical_face_ids":c.cylinder_faces.iter().map(|f|f.0)
                    .collect::<Vec<_>>(),
                "area":c.region_area,"perimeter":c.exposed_perimeter,
                "depth":c.depth
            })).collect::<Vec<_>>(),
            "topology":{
                "vertices":stats.vertices,"edges":stats.edges,
                "coedges":stats.coedges,"loops":stats.loops,
                "faces":stats.faces,"shells":stats.shells,
                "genus":stats.genus,
                "euler_characteristic":stats.euler_characteristic,
                "inner_loops":model.loops.iter().filter(|l|
                    l.role==rustsolid::LoopRole::Inner).count(),
                "cylindrical_faces":stats.cylindrical_faces,
                "circular_edges":stats.circular_edges
            },
            "mass":{
                "volume":model.mass.volume,
                "surface_area":model.mass.surface_area,
                "centroid":model.mass.centroid.array()
            },
            "revision":shape.revision(),
            "undo_depth":shape.undo_depth(),"redo_depth":shape.redo_depth(),
            "commits":commits,"transitions":transitions,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality,
            "generalized_topology_delta_can_apply":false
        });
        Ok((data,solid))
    })();
    match result{
        Ok((data,solid))=>Response{schema_version:"geometry.v1",ok:true,
            solid,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,data:None,error:Some(error.to_string())}
    }
}
