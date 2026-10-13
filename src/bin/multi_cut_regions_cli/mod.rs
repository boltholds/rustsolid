//! Read-only geometry.v1 API for multi-cut disk arrangements and depth slabs.
//! Every exposed arc belongs to an actual cutter circle, not a polygon facet.
//! Does not return B-rep faces, vertices or executable Boolean transactions.
use super::*;
use rustsolid::{MultiCutRegionPlanner,MultiCutPlanState,MultiCutHold};

pub(super) fn evaluate(
    source:SolidSource,
    through:Vec<mixed_regions_cli::ThroughHoleInput>,
    blind:Vec<BlindHoleInput>,
)->Response{
    let result=(|| -> Result<Value,GeometryError> {
        let solid=source.build()?;
        let through=through.into_iter().map(|h|HoleParameters{
            center_x:h.center[0],center_z:h.center[1],radius:h.radius
        }).collect::<Vec<_>>();
        let blind=blind.into_iter().map(|h|BlindHoleParameters{
            center_x:h.center[0],center_z:h.center[1],
            radius:h.radius,depth:h.depth
        }).collect::<Vec<_>>();
        let plan=MultiCutRegionPlanner::plan(&solid,&through,&blind)?;
        plan.validate_against(&solid,&through,&blind)?;
        let (state,hold)=match plan.state{
            MultiCutPlanState::Ready=>("ready",None),
            MultiCutPlanState::Held(reason)=>("held",Some(match reason{
                MultiCutHold::NearTangencyOrCoincidence=>"near_tangency_or_coincidence",
                MultiCutHold::NearTripleJunction=>"near_triple_junction",
                MultiCutHold::UnresolvedArc=>"unresolved_arc",
                MultiCutHold::IncompleteBoundaryCycles=>"incomplete_boundary_cycles",
                MultiCutHold::UnresolvedDepthEvent=>"unresolved_depth_event",
                MultiCutHold::NumericalAreaAmbiguity=>"numerical_area_ambiguity",
            }))
        };
        let slabs=plan.slabs.iter().map(|s|json!({
            "y_lower":s.y_lower,"y_upper":s.y_upper,
            "active_disks":s.active_disks,
            "union_area":s.union_area,
            "removed_volume":s.removed_volume,
            "area_resolution_budget":s.area_resolution_budget,
            "boundary_vertex_count":s.vertices.len(),
            "connected_union_components":s.connected_union_components,
            "planar_holes":s.planar_holes,
            "boundary_loops":s.boundary_loops.iter().map(|lp|json!({
                "is_hole":lp.is_hole,"signed_area":lp.signed_area,
                "arcs":lp.arcs.iter().map(|arc|json!({
                    "disk_index":arc.disk_index,
                    "start_vertex":arc.start_vertex,
                    "end_vertex":arc.end_vertex,
                    "start_angle":arc.start_angle,
                    "end_angle":arc.end_angle,
                    "midpoint":arc.midpoint.array()
                })).collect::<Vec<_>>()
            })).collect::<Vec<_>>()
        })).collect::<Vec<_>>();
        Ok(json!({
            "kind":"multi_cut_region_plan",
            "state":state,"hold_reason":hold,
            "through_count":through.len(),"blind_count":blind.len(),
            "source_cutters":through.iter().enumerate().map(|(i,t)|json!({
                "index":i,"kind":"through","center":[t.center_x,t.center_z],
                "radius":t.radius
            })).chain(blind.iter().enumerate().map(|(i,b)|json!({
                "index":i+through.len(),"kind":"blind",
                "center":[b.center_x,b.center_z],
                "radius":b.radius,"depth":b.depth
            }))).collect::<Vec<_>>(),
            "slabs":slabs,
            "removed_volume":plan.removed_volume,
            "predicted_remaining_volume":plan.predicted_remaining_volume,
            "length_budget":plan.length_budget,
            "source_revision":plan.source_revision,
            "can_apply":false,
            "topology_committed":false,
            "note":"Source-bound analytical multi-disk union, not a committed curved Boolean B-rep. Numerical slack is not a formal interval guarantee."
        }))
    })();
    match result {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(err)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(err.to_string()),data:None},
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const ROOT:&str=r#""source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6}"#;

    #[test]
    fn cli_three_cutters_with_two_blind_depths_return_three_full_boundary_graphs(){
        let req=format!(r#"{{"operation":"multi_cut_region_plan",{ROOT},
            "through":[{{"center":[0,0],"radius":0.75}}],
            "blind":[{{"center":[0.7,0],"radius":0.75,"depth":2}},
                     {{"center":[1.4,0.25],"radius":0.6,"depth":1}}]}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"ready");
        assert_eq!(d["through_count"],1);
        assert_eq!(d["blind_count"],2);
        assert_eq!(d["slabs"].as_array().unwrap().len(),3);
        assert_eq!(d["slabs"][0]["active_disks"],json!([0]));
        assert_eq!(d["slabs"][1]["active_disks"],json!([0,1]));
        assert_eq!(d["slabs"][2]["active_disks"],json!([0,1,2]));
        assert!(d["predicted_remaining_volume"].as_f64().unwrap()>0.0);
        assert_eq!(d["can_apply"],false);
        assert_eq!(d["topology_committed"],false);
    }

    #[test]
    fn cli_detects_exact_tangency_and_returns_no_mass_and_no_solid(){
        let req=format!(r#"{{"operation":"multi_cut_region_plan",{ROOT},
            "through":[{{"center":[0,0],"radius":0.7}}],
            "blind":[{{"center":[1.4,0],"radius":0.7,"depth":1.5}}]}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"held");
        assert_eq!(d["hold_reason"],"near_tangency_or_coincidence");
        assert!(d["slabs"].as_array().unwrap().is_empty());
        assert!(d["predicted_remaining_volume"].is_null());
        assert_eq!(d["can_apply"],false);
    }

    #[test]
    fn cli_invalid_source_rejects_without_partial_data(){
        let req=format!(r#"{{"operation":"multi_cut_region_plan",{ROOT},
            "through":[{{"center":[0,0],"radius":0.7}}],
            "blind":[{{"center":[0.7,0],"radius":0.7,"depth":4.0}}]}}"#);
        let r=process(&req);
        assert!(!r.ok);
        assert!(r.solid.is_none()&&r.data.is_none());
    }
}
