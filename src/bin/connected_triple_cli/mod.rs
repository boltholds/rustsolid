//! Dedicated geometry.v1 constructor for one connected through/blind/blind
//! cylinder chain, with source-atomic edits, two real analytic blind floors,
//! optional closed watertight faceting and exact-in-form source mass.
//! This does NOT enable the generic arbitrary curved Boolean API.
use super::*;
use rustsolid::{ConnectedTripleBrep,ConnectedTripleDisplay,TripleCutParameters};

#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum TripleEditInput {
    ReplaceThrough {center:[f64;2],radius:f64},
    ReplaceBlind {index:usize,center:[f64;2],radius:f64,depth:f64},
}
fn to_through(value:mixed_regions_cli::ThroughHoleInput)->HoleParameters{
    HoleParameters{center_x:value.center[0],center_z:value.center[1],
        radius:value.radius}
}
fn to_blind(value:BlindHoleInput)->BlindHoleParameters{
    BlindHoleParameters{center_x:value.center[0],center_z:value.center[1],
        radius:value.radius,depth:value.depth}
}
pub(super) fn evaluate(
    source:SolidSource,through:mixed_regions_cli::ThroughHoleInput,
    blind:Vec<BlindHoleInput>,batches:Vec<Vec<TripleEditInput>>,
    undo:usize,redo:usize,display_segments:Option<usize>,
)->Response{
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if blind.len()!=2||batches.len()>64||undo>64||redo>64||
            batches.iter().any(|b|b.is_empty()||b.len()>64){
            return Err(GeometryError::InvalidEdit(
                "connected triple accepts exactly two ordered blind bores and bounded source/history batches".into()));
        }
        let box_=source.build()?;
        let mut blinds=blind.into_iter().map(to_blind);
        let params=TripleCutParameters {
            through:to_through(through),
            deep_blind:blinds.next().unwrap(),
            shallow_blind:blinds.next().unwrap(),
        };
        let mut body=ConnectedTripleBrep::from_box(&box_,params)?;
        let mut commits=Vec::new();
        for items in batches {
            let count=items.len();
            let mut tx=body.begin_edit();
            for item in items {
                match item {
                    TripleEditInput::ReplaceThrough{center,radius}=>{
                        tx.replace_through(HoleParameters{
                            center_x:center[0],center_z:center[1],radius})?;
                    }
                    TripleEditInput::ReplaceBlind{index,center,radius,depth}=>{
                        tx.replace_blind(index,BlindHoleParameters{
                            center_x:center[0],center_z:center[1],radius,depth})?;
                    }
                }
            }
            let rev=tx.commit()?;
            commits.push(json!({"revision":rev,"changes":count}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            if !body.undo()?{
                return Err(GeometryError::InvalidEdit(
                    "connected triple undo exceeds available source history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":body.revision()}));
        }
        for _ in 0..redo {
            if !body.redo()?{
                return Err(GeometryError::InvalidEdit(
                    "connected triple redo exceeds available source history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":body.revision()}));
        }
        body.validate()?;
        let model=body.model();let summary=model.summary();
        let expected=body.plan().predicted_remaining_volume
            .ok_or_else(||GeometryError::InvalidEdit("missing source volume witness".into()))?;
        if (model.mass.volume-expected).abs()>model.tolerance.length_at(1.0)*128.0 ||
            summary.genus!=1 || summary.euler_characteristic!=0 ||
            model.faces[0].loops.len()!=2 ||
            model.faces[1].loops.len()!=2 {
            return Err(GeometryError::InvalidTopology(
                "connected triple source mass, Euler or cap ownership mismatch".into()));
        }
        let preview=display_segments.map(|n|body.tessellate(n)).transpose()?;
        let quality=preview.as_ref().map(|d:&ConnectedTripleDisplay|json!({
            "closed_oriented_manifold_certified":true,
            "genus":d.genus,"floor_faces":d.floor_faces,
            "segments_per_analytic_arc":d.segments_per_arc,
            "triangle_count":d.mesh.triangles.len(),
            "triangle_edges":d.triangle_edges,
            "analytic_volume":d.analytic_volume,
            "facet_volume":d.facet_volume,
            "facet_surface_area":d.facet_area,
            "faceting_volume_error":d.faceting_volume_error,
            "independent_polygon_volume":d.independently_polygonized_volume
        }));
        let solid=preview.map(|d|SolidResponse{
            vertices:d.mesh.vertices.iter().map(|p|p.array()).collect(),
            faces:d.mesh.triangles,
            triangle_face_ids:d.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{
                vertices:summary.vertices,edges:summary.edges,
                coedges:summary.coedges,loops:summary.loops,
                shells:summary.shells,faces:summary.faces,
                euler_characteristic:summary.euler_characteristic
            },
            mass:Mass{
                volume:model.mass.volume,
                surface_area:model.mass.surface_area,
                centroid:model.mass.centroid.array()
            },
            bbox:[model.bbox.min.array(),model.bbox.max.array()],
            edit_report:None,history:None,
            analytic_brep:Some(json!({
                "representation":"connected_triple_analytic_brep",
                "through_count":1,"blind_count":2,
                "inner_loops":2,"material_genus":1,
                "shared_coedges_certified":true,
                "analytic_mass_authoritative":true,
                "mesh_is_display_only":true
            })),
            brep:None,
        });
        let params=body.parameters();
        let data=json!({
            "kind":"connected_triple_brep",
            "state":"committed","representation":"analytic_connected_triple",
            "source_constraints":"one through, deep blind crossing through, shallow blind crossing deep blind; through/shallow nonintersecting; distinct floor levels",
            "topology_committed":true,
            "generalized_topology_delta_can_apply":false,
            "source_through":{"center":[params.through.center_x,
                params.through.center_z],"radius":params.through.radius},
            "source_blind":[
                {"center":[params.deep_blind.center_x,params.deep_blind.center_z],
                 "radius":params.deep_blind.radius,"depth":params.deep_blind.depth},
                {"center":[params.shallow_blind.center_x,params.shallow_blind.center_z],
                 "radius":params.shallow_blind.radius,"depth":params.shallow_blind.depth}
            ],
            "topology":{
                "vertices":summary.vertices,"edges":summary.edges,
                "coedges":summary.coedges,"loops":summary.loops,
                "faces":summary.faces,"shells":summary.shells,
                "genus":summary.genus,
                "euler_characteristic":summary.euler_characteristic,
                "inner_cap_loops":2,
                "blind_floor_faces":2,
                "cylindrical_faces":summary.cylindrical_faces,
                "circular_edges":summary.circular_edges
            },
            "mass":{
                "volume":model.mass.volume,
                "surface_area":model.mass.surface_area,
                "centroid":model.mass.centroid.array()
            },
            "slabs":body.plan().slabs.iter().map(|s|json!({
                "y_lower":s.y_lower,"y_upper":s.y_upper,
                "active_disks":s.active_disks,
                "union_area":s.union_area,
                "components":s.connected_union_components,
                "planar_holes":s.planar_holes
            })).collect::<Vec<_>>(),
            "revision":body.revision(),
            "undo_depth":body.undo_depth(),
            "redo_depth":body.redo_depth(),
            "commits":commits,"transitions":transitions,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality
        });
        Ok((data,solid))
    })();
    match result{
        Ok((data,solid))=>Response{
            schema_version:"geometry.v1",ok:true,error:None,solid,data:Some(data)
        },
        Err(err)=>Response{
            schema_version:"geometry.v1",ok:false,
            error:Some(err.to_string()),solid:None,data:None
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const SRC:&str=r#""source":{"kind":"block","origin":[-3,-1,-3],
        "width":6,"height":4,"depth":6}"#;
    const CUTS:&str=r#""through":{"center":[0,0],"radius":0.75},
        "blind":[{"center":[0.7,0],"radius":0.75,"depth":2.0},
                 {"center":[1.4,0.25],"radius":0.6,"depth":1.0}]"#;

    #[test]
    fn cli_exposes_true_connected_analytic_shell_with_optional_watertight_display(){
        let request=format!(r#"{{"operation":"compose_connected_triple",{SRC},{CUTS},
            "display_segments":24}}"#);
        let result=process(&request);
        assert!(result.ok,"{:?}",result.error);
        let data=result.data.unwrap();
        assert_eq!(data["state"],"committed");
        assert_eq!(data["topology_committed"],true);
        assert_eq!(data["generalized_topology_delta_can_apply"],false);
        assert_eq!(data["topology"]["genus"],1);
        assert_eq!(data["topology"]["blind_floor_faces"],2);
        assert_eq!(data["slabs"].as_array().unwrap().len(),3);
        assert!(data["display_quality"]["closed_oriented_manifold_certified"].as_bool().unwrap());
        let mesh=result.solid.unwrap();
        assert_eq!(mesh.faces.len(),mesh.triangle_face_ids.len());
        assert!((data["mass"]["volume"].as_f64().unwrap()-
            data["display_quality"]["analytic_volume"].as_f64().unwrap()).abs()<1e-10);
    }

    #[test]
    fn atomic_parameter_changes_and_history_replay_are_source_certified(){
        let request=format!(r#"{{"operation":"compose_connected_triple",{SRC},{CUTS},
            "batches":[[
                {{"kind":"replace_blind","index":1,"center":[1.43,0.28],
                  "radius":0.58,"depth":1.1}}
            ]],"undo":1,"redo":1}}"#);
        let result=process(&request);
        assert!(result.ok,"{:?}",result.error);
        let data=result.data.unwrap();
        assert_eq!(data["revision"],3);
        assert_eq!(data["commits"][0]["changes"],1);
        assert_eq!(data["transitions"].as_array().unwrap().len(),2);
    }
    #[test]
    fn invalid_connected_triple_never_publishes_partial_geometry(){
        for changed in [
            r#""blind":[{"center":[0.7,0],"radius":0.75,"depth":2},
                        {"center":[0.2,0],"radius":0.6,"depth":1}]"#,
            r#""blind":[{"center":[0.7,0],"radius":0.75,"depth":2},
                        {"center":[1.4,0.25],"radius":0.6,"depth":3}]"#,
        ] {
            let request=format!(r#"{{"operation":"compose_connected_triple",{SRC},
                "through":{{"center":[0,0],"radius":0.75}},{changed}}}"#);
            let result=process(&request);
            assert!(!result.ok);
            assert!(result.data.is_none() && result.solid.is_none());
        }
    }
}
