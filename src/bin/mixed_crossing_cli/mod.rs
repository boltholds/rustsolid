//! Source-certified executed through/blind crossing B-rep through geometry.v1.
//! Only one strict transverse pair is supported; never run the general Boolean
//! planner or substitute facet meshes for analytic geometry.
use super::*;
use rustsolid::{MixedCrossingBrep, MixedFootprintContact};

#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum CrossingEditInput {
    ReplaceThrough {center:[f64;2],radius:f64},
    ReplaceBlind {center:[f64;2],radius:f64,depth:f64},
}
pub(super) fn evaluate(
    source:SolidSource,
    through:mixed_regions_cli::ThroughHoleInput,
    blind:BlindHoleInput,
    batches:Vec<Vec<CrossingEditInput>>,
    undo:usize,redo:usize,display_segments:Option<usize>,
)->Response {
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if batches.len()>64||undo>64||redo>64{
            return Err(GeometryError::InvalidEdit(
                "crossing request exceeds 64 batches or history operations".into()));
        }
        let source=source.build()?;
        let tp=HoleParameters{center_x:through.center[0],center_z:through.center[1],
            radius:through.radius};
        let bp=BlindHoleParameters{center_x:blind.center[0],center_z:blind.center[1],
            radius:blind.radius,depth:blind.depth};
        let mut body=MixedCrossingBrep::from_box(&source,tp,bp)?;
        let mut commits=Vec::new();
        for edits in batches {
            if edits.is_empty()||edits.len()>32 {
                return Err(GeometryError::InvalidEdit(
                    "crossing edit batch requires 1..=32 source changes".into()));
            }
            let len=edits.len();
            let mut tx=body.begin_edit();
            for edit in edits {
                match edit {
                    CrossingEditInput::ReplaceThrough{center,radius} =>
                        tx.replace_through(HoleParameters{
                            center_x:center[0],center_z:center[1],radius
                        })?,
                    CrossingEditInput::ReplaceBlind{center,radius,depth} =>
                        tx.replace_blind(BlindHoleParameters{
                            center_x:center[0],center_z:center[1],radius,depth
                        })?,
                }
            }
            let revision=tx.commit()?;
            commits.push(json!({"revision":revision,"changes":len}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            if !body.undo()? {
                return Err(GeometryError::InvalidEdit(
                    "crossing Undo exceeds committed source history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":body.revision()}));
        }
        for _ in 0..redo {
            if !body.redo()? {
                return Err(GeometryError::InvalidEdit(
                    "crossing Redo exceeds committed source history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":body.revision()}));
        }
        body.validate()?;
        let m=body.model();let s=m.summary();
        let display=display_segments.map(|n|body.tessellate(n)).transpose()?;
        let quality=display.as_ref().map(|d|json!({
            "segments_per_analytic_arc":d.circle_segments_per_arc,
            "closed_oriented_manifold_certified":true,
            "genus":d.genus,
            "facet_triangles":d.mesh.triangles.len(),
            "facet_edges":d.triangle_edges,
            "faceted_volume":d.facet_volume,
            "faceted_surface_area":d.facet_area,
            "analytic_volume":d.analytic_volume,
            "volume_faceting_error":d.facet_volume_error,
            "expected_polygonal_volume":d.expected_polygonal_volume,
            "top_cap_polygon_area":d.top_cap_polygon_area,
            "bottom_cap_polygon_area":d.bottom_cap_polygon_area,
            "floor_polygon_area":d.floor_polygon_area,
        }));
        let solid=display.map(|d|SolidResponse{
            vertices:d.mesh.vertices.iter().map(|p|p.array()).collect(),
            faces:d.mesh.triangles,
            triangle_face_ids:d.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{
                vertices:s.vertices,edges:s.edges,coedges:s.coedges,loops:s.loops,
                shells:s.shells,faces:s.faces,euler_characteristic:s.euler_characteristic,
            },
            mass:Mass{volume:m.mass.volume,surface_area:m.mass.surface_area,
                centroid:m.mass.centroid.array()},
            bbox:[m.bbox.min.array(),m.bbox.max.array()],
            edit_report:None,history:None,analytic_brep:Some(json!({
                "representation":"analytic_crossing_through_blind",
                "inner_loops":2,
                "genus":1,"cylindrical_faces":4,"floor_face_id":10,
                "exact_mass_not_taken_from_mesh":true,
            })),brep:None,
        });
        let p=body.plan();
        let contact=match p.footprint.contact{
            MixedFootprintContact::TwoTransverseCrossings=>"two_transverse_crossings",
            _=>return Err(GeometryError::InvalidEdit(
                "executed crossing must have certified transverse source".into())),
        };
        let data=json!({
            "kind":"mixed_crossing_brep",
            "state":"committed",
            "representation":"analytic_crossing_through_blind",
            "source_contact":contact,
            "revision":body.revision(),
            "undo_depth":body.undo_depth(),"redo_depth":body.redo_depth(),
            "through":{"center":[body.through().center_x,body.through().center_z],
                "radius":body.through().radius},
            "blind":{"center":[body.blind().center_x,body.blind().center_z],
                "radius":body.blind().radius,"depth":body.blind().depth},
            "topology_committed":true,
            "generalized_topology_delta_can_apply":false,
            "source_identity_stable_across_repartition":false,
            "topology":{"vertices":s.vertices,"edges":s.edges,"coedges":s.coedges,
                "loops":s.loops,"faces":s.faces,"shells":s.shells,
                "genus":s.genus,"euler_characteristic":s.euler_characteristic,
                "inner_loops":2,"circular_edges":s.circular_edges,
                "cylindrical_faces":s.cylindrical_faces},
            "topology_roles":{"through_lower_faces":[6,7],"through_upper_face":8,
                "blind_upper_face":9,"blind_floor_face":10,
                "top_aperture_arcs":p.top_boundary_arc_count,
                "bottom_aperture_arcs":p.bottom_boundary_arc_count},
            "mass":{"volume":m.mass.volume,"surface_area":m.mass.surface_area,
                "centroid":m.mass.centroid.array()},
            "source_crossings":p.footprint.crossing_points.map(|items|
                items.into_iter().map(Point3::array).collect::<Vec<_>>()),
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality,
            "commits":commits,"transitions":transitions,
        });
        Ok((data,solid))
    })();
    match result {
        Ok((data,solid))=>Response {
            schema_version:"geometry.v1",ok:true,solid,error:None,data:Some(data),
        },
        Err(err)=>Response{
            schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(err.to_string()),data:None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cli_returns_real_genus_one_brep_and_analytic_mass(){
        let request=r#"{"operation":"compose_mixed_crossing",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5}}"#;
        let result=process(request);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let d=result.data.unwrap();
        assert_eq!(d["state"],"committed");
        assert_eq!(d["representation"],"analytic_crossing_through_blind");
        assert_eq!(d["topology"]["genus"],1);
        assert_eq!(d["topology"]["euler_characteristic"],0);
        assert_eq!(d["topology"]["coedges"],46);
        assert_eq!(d["topology"]["inner_loops"],2);
        assert!(d["topology_committed"].as_bool().unwrap());
        assert_eq!(d["generalized_topology_delta_can_apply"],false);
        assert!((d["mass"]["volume"].as_f64().unwrap()-
            140.14094873829234).abs()<1e-10);
    }
    #[test]
    fn cli_includes_only_certified_mesh_if_requested(){
        let result=process(r#"{"operation":"compose_mixed_crossing",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5},
            "display_segments":32}"#);
        assert!(result.ok,"{:?}",result.error);
        let d=result.data.unwrap();
        let solid=result.solid.unwrap();
        assert_eq!(d["display_quality"]["genus"],1);
        assert_eq!(d["display_quality"]["closed_oriented_manifold_certified"],true);
        assert_eq!(solid.faces.len(),solid.triangle_face_ids.len());
        assert_eq!(solid.analytic_brep.unwrap()["genus"],1);
    }
    #[test]
    fn cli_parameter_changes_are_atomic_and_history_is_replayable(){
        let result=process(r#"{"operation":"compose_mixed_crossing",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5},
            "batches":[[
                {"kind":"replace_through","center":[0,0],"radius":0.55},
                {"kind":"replace_blind","center":[0.55,0],"radius":0.5,"depth":1.4}
            ]],"undo":1,"redo":1}"#);
        assert!(result.ok,"{:?}",result.error);
        let d=result.data.unwrap();
        assert_eq!(d["revision"],3);
        assert_eq!(d["commits"][0]["changes"],2);
        assert_eq!(d["transitions"].as_array().unwrap().len(),2);
        assert_eq!(d["topology"]["genus"],1);
    }
    #[test]
    fn cli_rejects_tangency_late_commit_and_invalid_mesh_without_partial_output(){
        for request in [
            r#"{"operation":"compose_mixed_crossing",
                "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
                "through":{"center":[0,0],"radius":0.5},
                "blind":{"center":[1,0],"radius":0.5,"depth":1.5}}"#,
            r#"{"operation":"compose_mixed_crossing",
                "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
                "through":{"center":[0,0],"radius":0.5},
                "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5},
                "display_segments":7}"#,
            r#"{"operation":"compose_mixed_crossing",
                "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
                "through":{"center":[0,0],"radius":0.5},
                "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5},
                "batches":[[
                    {"kind":"replace_through","center":[0,0],"radius":0.55},
                    {"kind":"replace_blind","center":[1.05,0],"radius":0.5,"depth":1.5}
                ]]}"#,
        ]{
            let r=process(request);
            assert!(!r.ok);
            assert!(r.solid.is_none()&&r.data.is_none());
        }
    }
}
