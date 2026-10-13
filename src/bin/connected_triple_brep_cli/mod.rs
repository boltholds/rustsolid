//! Dedicated and strictly bounded executed T-A-B analytic triple B-rep.
//! Differs from the read-only connected_triple_stitch_plan: this action
//! returns one validated canonical shell and, optionally, certified facets.
//! No general curved Boolean, triple junction or B-rep topology delta unlock.
use super::*;
use rustsolid::{
    ConnectedTripleBrep,ConnectedTripleDisplay,TripleChainParameters,
};
#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct TripleInput {
    through:mixed_regions_cli::ThroughHoleInput,
    deeper_blind:BlindHoleInput,
    shallower_blind:BlindHoleInput,
}
impl TripleInput{
    fn geometry(self)->TripleChainParameters{
        TripleChainParameters{
            through:HoleParameters{
                center_x:self.through.center[0],
                center_z:self.through.center[1],radius:self.through.radius
            },
            deeper_blind:BlindHoleParameters{
                center_x:self.deeper_blind.center[0],
                center_z:self.deeper_blind.center[1],
                radius:self.deeper_blind.radius,depth:self.deeper_blind.depth
            },
            shallower_blind:BlindHoleParameters{
                center_x:self.shallower_blind.center[0],
                center_z:self.shallower_blind.center[1],
                radius:self.shallower_blind.radius,
                depth:self.shallower_blind.depth
            }
        }
    }
}
pub(super) fn evaluate(
    source:SolidSource,through:mixed_regions_cli::ThroughHoleInput,
    deeper_blind:BlindHoleInput,shallower_blind:BlindHoleInput,
    batches:Vec<Vec<TripleInput>>,undo:usize,redo:usize,
    display_segments:Option<usize>,
)->Response{
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if batches.len()>64||undo>64||redo>64||
            batches.iter().any(|b|b.is_empty()||b.len()>64){
            return Err(GeometryError::InvalidEdit(
                "triple source batches and undo/redo exceed bounded request limits".into()));
        }
        let box_=source.build()?;
        let first=TripleInput{through,deeper_blind,shallower_blind}.geometry();
        let mut body=ConnectedTripleBrep::from_box(&box_,first)?;
        let mut revisions=Vec::new();
        for edits in batches{
            let len=edits.len();
            let mut tx=body.begin_edit();
            for edit in edits{
                tx.replace_sources(edit.geometry())?;
            }
            let revision=tx.commit()?;
            revisions.push(json!({"revision":revision,"commands":len}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo{
            if !body.undo()?{
                return Err(GeometryError::InvalidEdit(
                    "triple undo exceeds committed history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":body.revision()}));
        }
        for _ in 0..redo{
            if !body.redo()?{
                return Err(GeometryError::InvalidEdit(
                    "triple redo exceeds committed history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":body.revision()}));
        }
        body.validate()?;
        let m=body.model();
        let stats=m.summary();
        if stats.vertices!=22||stats.edges!=38||stats.coedges!=76||
            stats.faces!=18||stats.euler_characteristic!=0||stats.genus!=1{
            return Err(GeometryError::InvalidTopology(
                "connected triple B-rep does not satisfy canonical source topology".into()));
        }
        let quality=display_segments.map(|n|body.tessellate(n)).transpose()?;
        let quality_meta=quality.as_ref().map(|d:&ConnectedTripleDisplay|json!({
            "closed_oriented_manifold_certified":true,
            "segments_per_analytic_arc":d.segments_per_source_arc,
            "genus":d.genus,
            "triangle_edges":d.triangle_edges,
            "triangle_count":d.mesh.triangles.len(),
            "faceted_volume":d.faceted_volume,
            "faceted_surface_area":d.faceted_surface_area,
            "analytical_volume":d.analytical_volume,
            "faceting_volume_error":d.faceting_volume_error,
            "arc_chord_volume_error_bound":d.arc_chord_volume_error_bound,
            "independent_polygonized_volume":d.predicted_polygonized_volume,
            "independent_polygonized_area":d.predicted_polygonized_area,
            "floor_polygon_areas":d.floor_polygon_areas
        }));
        let solid=quality.map(|d|SolidResponse{
            vertices:d.mesh.vertices.iter().map(|p|p.array()).collect(),
            faces:d.mesh.triangles,
            triangle_face_ids:d.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{
                vertices:stats.vertices,edges:stats.edges,coedges:stats.coedges,
                loops:stats.loops,shells:stats.shells,faces:stats.faces,
                euler_characteristic:stats.euler_characteristic
            },
            mass:Mass{
                volume:m.mass.volume,surface_area:m.mass.surface_area,
                centroid:m.mass.centroid.array()
            },
            bbox:[m.bbox.min.array(),m.bbox.max.array()],
            edit_report:None,history:None,
            analytic_brep:Some(json!({
                "representation":"analytic_connected_triple",
                "closed":true,"genus":1,
                "floor_face_ids":[16,17],
                "inner_loops":2,
                "source_analytic_mass_authoritative":true,
                "display_mesh_is_derived":true
            })),
            brep:None,
        });
        let plan=body.stitch();
        let inputs=body.inputs();
        let data=json!({
            "kind":"connected_triple_brep",
            "state":"committed",
            "representation":"analytic_connected_triple",
            "topology_committed":true,
            "generalized_topology_delta_can_apply":false,
            "source_stitch_plan_can_apply":false,
            "source":{
                "through":{"center":[inputs.through.center_x,inputs.through.center_z],
                    "radius":inputs.through.radius},
                "deeper_blind":{
                    "center":[inputs.deeper_blind.center_x,inputs.deeper_blind.center_z],
                    "radius":inputs.deeper_blind.radius,
                    "depth":inputs.deeper_blind.depth},
                "shallower_blind":{
                    "center":[inputs.shallower_blind.center_x,inputs.shallower_blind.center_z],
                    "radius":inputs.shallower_blind.radius,
                    "depth":inputs.shallower_blind.depth}
            },
            "topology":{"vertices":stats.vertices,"edges":stats.edges,
                "coedges":stats.coedges,"faces":stats.faces,
                "loops":stats.loops,"shells":stats.shells,
                "inner_loops":2,"genus":stats.genus,
                "euler_characteristic":stats.euler_characteristic,
                "circular_edges":stats.circular_edges,
                "cylindrical_faces":stats.cylindrical_faces},
            "mass":{"volume":m.mass.volume,"surface_area":m.mass.surface_area,
                "centroid":m.mass.centroid.array()},
            "height_levels":plan.levels,
            "global_source_junctions":plan.junctions.len(),
            "global_atomic_arcs":plan.atomic_arcs.len(),
            "depth_floor_count":plan.floors.len(),
            "depth_floor_face_ids":[16,17],
            "slab_areas":plan.slabs.iter().map(|s|s.union_area).collect::<Vec<_>>(),
            "revision":body.revision(),
            "undo_depth":body.undo_depth(),"redo_depth":body.redo_depth(),
            "commits":revisions,"transitions":transitions,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality_meta,
            "scope":"One connected T-A-B chain; T/B disjoint, T/A and A/B strictly transverse, distinct blind depth levels. No arbitrary triple junctions."
        });
        Ok((data,solid))
    })();
    match result{
        Ok((data,solid))=>Response{schema_version:"geometry.v1",ok:true,
            solid,error:None,data:Some(data)},
        Err(e)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(e.to_string()),data:None}
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    const ROOT:&str=r#""source":{"kind":"block",
        "origin":[-4,-2,-3],"width":8,"height":5,"depth":6}"#;
    const PARAMS:&str=r#""through":{"center":[-1.15,0],"radius":0.8},
        "deeper_blind":{"center":[0,0],"radius":0.8,"depth":2.8},
        "shallower_blind":{"center":[1.15,0],"radius":0.8,"depth":1.2}"#;
    #[test]
    fn cli_executes_connected_three_source_canonical_brep_with_two_real_floors(){
        let req=format!(r#"{{"operation":"compose_connected_triple",{ROOT},{PARAMS},
            "display_segments":32}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        let d=r.data.unwrap();
        assert_eq!(d["state"],"committed");
        assert_eq!(d["topology_committed"],true);
        assert_eq!(d["topology"]["vertices"],22);
        assert_eq!(d["topology"]["edges"],38);
        assert_eq!(d["topology"]["coedges"],76);
        assert_eq!(d["topology"]["genus"],1);
        assert_eq!(d["depth_floor_face_ids"],json!([16,17]));
        assert_eq!(d["display_quality"]["closed_oriented_manifold_certified"],true);
        assert!(d["display_quality"]["faceting_volume_error"].as_f64().unwrap()>=-2e-6);
        let mesh=r.solid.unwrap();
        assert_eq!(mesh.faces.len(),mesh.triangle_face_ids.len());
        assert_eq!(mesh.analytic_brep.unwrap()["genus"],1);
    }
    #[test]
    fn triple_source_transaction_can_commit_undo_and_redo(){
        let req=format!(r#"{{"operation":"compose_connected_triple",{ROOT},{PARAMS},
            "batches":[[{{"through":{{"center":[-1.2,0],"radius":0.8}},
                "deeper_blind":{{"center":[0,0],"radius":0.8,"depth":2.8}},
                "shallower_blind":{{"center":[1.15,0],"radius":0.8,"depth":1.35}}}}]],
            "undo":1,"redo":1}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        let d=r.data.unwrap();
        assert_eq!(d["revision"],3);
        assert_eq!(d["commits"][0]["commands"],1);
        assert_eq!(d["transitions"].as_array().unwrap().len(),2);
        assert_eq!(d["topology"]["genus"],1);
    }
    #[test]
    fn connected_triple_bad_depth_or_through_shallow_overlap_fails_closed(){
        for bad in [
            r#""through":{"center":[-1.15,0],"radius":0.8},
               "deeper_blind":{"center":[0,0],"radius":0.8,"depth":2.8},
               "shallower_blind":{"center":[-0.4,0],"radius":0.8,"depth":1.2}"#,
            r#""through":{"center":[-1.15,0],"radius":0.8},
               "deeper_blind":{"center":[0,0],"radius":0.8,"depth":2.8},
               "shallower_blind":{"center":[1.15,0],"radius":0.8,"depth":5.0}"#,
        ]{
            let req=format!(r#"{{"operation":"compose_connected_triple",{ROOT},{bad},
                "display_segments":32}}"#);
            let r=process(&req);
            assert!(!r.ok);
            assert!(r.solid.is_none()&&r.data.is_none());
        }
    }
}
