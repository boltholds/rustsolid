//! Read-only three-cylinder connected source chain: global 3D halfedge
//! incidence, seam atomization and two source-bound finite floor interfaces.
use super::*;
use rustsolid::{
    ConnectedTripleStitchPlanner,TripleChainParameters,
    TripleStitchState,TripleStitchHold,TriplePhysicalEdgeKey
};

pub(super) fn evaluate(source:SolidSource,
    through:mixed_regions_cli::ThroughHoleInput,
    deeper_blind:BlindHoleInput,shallower_blind:BlindHoleInput)->Response {
    let result=(|| -> Result<Value,GeometryError> {
        let box_=source.build()?;
        let p=TripleChainParameters {
            through:HoleParameters {
                center_x:through.center[0],center_z:through.center[1],
                radius:through.radius
            },
            deeper_blind:BlindHoleParameters {
                center_x:deeper_blind.center[0],
                center_z:deeper_blind.center[1],
                radius:deeper_blind.radius,depth:deeper_blind.depth
            },
            shallower_blind:BlindHoleParameters {
                center_x:shallower_blind.center[0],
                center_z:shallower_blind.center[1],
                radius:shallower_blind.radius,depth:shallower_blind.depth
            }
        };
        let plan=ConnectedTripleStitchPlanner::plan(&box_,p)?;
        plan.validate_against(&box_,p)?;
        let (state,hold)=match plan.state{
            TripleStitchState::Ready=>("ready",None),
            TripleStitchState::Held(reason)=>("held",Some(match reason{
                TripleStitchHold::UnresolvedSourceArrangement=>"unresolved_source_arrangement",
                TripleStitchHold::UnsupportedCircleContact=>"unsupported_circle_contact",
                TripleStitchHold::NearTangency=>"near_tangency",
                TripleStitchHold::AmbiguousArcExposure=>"ambiguous_arc_exposure",
                TripleStitchHold::OpenPlanarContour=>"open_planar_contour",
                TripleStitchHold::NonmanifoldPhysicalEdge=>"nonmanifold_physical_edge",
                TripleStitchHold::InconsistentAnalyticArea=>"inconsistent_analytic_area",
                TripleStitchHold::InconsistentEulerCharacteristic=>"inconsistent_euler"
            }))
        };
        let slabs=plan.slabs.iter().map(|s|json!({
            "y_lower":s.y_lower,"y_upper":s.y_upper,
            "active_source_indices":s.active_sources,
            "union_area":s.union_area,
            "removed_volume":s.removed_volume,
            "contour_arc_ids":s.exposed_contour.parts.iter()
                .map(|a|a.arc).collect::<Vec<_>>(),
            "closed":true
        })).collect::<Vec<_>>();
        let floors=plan.floors.iter().map(|s|json!({
            "y":s.y,
            "area":s.floor_area,
            "closed":true,
            "old_source_boundary_arcs":s.incoming_lower_arcs,
            "new_source_boundary_arcs":s.born_upper_arcs,
            "oriented_uses":s.contour.parts.iter().map(|x|json!({
                "atom":x.arc,"reversed_for_ccw_area":x.reversed
            })).collect::<Vec<_>>()
        })).collect::<Vec<_>>();
        Ok(json!({
            "kind":"connected_triple_stitch_plan",
            "state":state,"hold_reason":hold,
            "can_apply":false,"topology_committed":false,
            "source_revision":plan.source_revision,
            "source":{"through":{"center":[p.through.center_x,p.through.center_z],
                "radius":p.through.radius},
                "deeper_blind":{"center":[p.deeper_blind.center_x,p.deeper_blind.center_z],
                    "radius":p.deeper_blind.radius,"depth":p.deeper_blind.depth},
                "shallower_blind":{"center":[p.shallower_blind.center_x,p.shallower_blind.center_z],
                    "radius":p.shallower_blind.radius,"depth":p.shallower_blind.depth}},
            "levels":plan.levels,
            "junctions":plan.junctions.iter().map(|j|json!({
                "id":j.id,"source_pair":j.source_pair,"point_xz":j.position_xz
            })).collect::<Vec<_>>(),
            "atoms":plan.atomic_arcs.iter().map(|a|json!({
                "id":a.id,"source":a.source,
                "start_junction":a.start_junction,
                "end_junction":a.end_junction,
                "start_angle":a.start_angle,"end_angle":a.end_angle,
                "length":a.arc_length,
                "exposed_in_slabs":a.visible_in_slab
            })).collect::<Vec<_>>(),
            "slabs":slabs,"floor_interfaces":floors,
            "source_wall_patches":plan.walls.iter().map(|w|json!({
                "source":w.source,"atom":w.arc,"slab":w.slab,
                "y_lower":w.y_lower,"y_upper":w.y_upper,
                "uv_bounds":w.uv_bounds,
                "corners":w.corners.map(Point3::array),
                "analytic_area":w.analytic_area,
                "material_orientation":"inward_radial"
            })).collect::<Vec<_>>(),
            "edge_incidence":plan.shared_physical_edges.iter().map(|e|json!({
                "key":match e.key{
                    TriplePhysicalEdgeKey::CircleArc{arc,level}=>
                        json!({"kind":"circular_arc","atom":arc,"level":level}),
                    TriplePhysicalEdgeKey::VerticalTrace{junction,slab}=>
                        json!({"kind":"vertical_trace","junction":junction,"slab":slab})
                },
                "use_orientation_signs":e.opposite_use_signs
            })).collect::<Vec<_>>(),
            "predicted_shell":{
                "vertices":plan.predicted_shell_vertices,
                "edges":plan.predicted_shell_edges,
                "faces":plan.predicted_shell_faces,
                "euler":plan.predicted_euler,
                "genus":plan.predicted_genus
            },
            "analytic_surface_area":plan.analytic_surface_area,
            "predicted_remaining_volume":plan.predicted_remaining_volume,
            "length_budget":plan.length_budget,
            "note":"Halfedge incidence and analytical UV patch geometry are a PREVIEW, not a constructed canonical B-rep or certified watertight display."
        }))
    })();
    match result{
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(e)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(e.to_string()),data:None}
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    const SOURCE:&str=r#""source":{"kind":"block",
        "origin":[-4,-2,-3],"width":8,"height":5,"depth":6}"#;
    const FIXTURE:&str=r#""through":{"center":[-1.15,0],"radius":0.8},
        "deeper_blind":{"center":[0,0],"radius":0.8,"depth":2.8},
        "shallower_blind":{"center":[1.15,0],"radius":0.8,"depth":1.2}"#;
    #[test]
    fn three_connected_sources_return_global_physical_incidence_without_solid(){
        let req=format!(r#"{{"operation":"connected_triple_stitch_plan",{SOURCE},{FIXTURE}}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"ready");
        assert_eq!(d["junctions"].as_array().unwrap().len(),4);
        assert_eq!(d["atoms"].as_array().unwrap().len(),8);
        assert_eq!(d["slabs"].as_array().unwrap().len(),3);
        assert_eq!(d["floor_interfaces"].as_array().unwrap().len(),2);
        assert_eq!(d["predicted_shell"]["genus"],1);
        assert_eq!(d["predicted_shell"]["euler"],0);
        assert_eq!(d["predicted_shell"]["edges"],38);
        assert_eq!(d["edge_incidence"].as_array().unwrap().len(),26);
        assert_eq!(d["can_apply"],false);
        assert_eq!(d["topology_committed"],false);
    }
    #[test]
    fn unsupported_third_circle_contact_returns_held_without_false_topology(){
        let req=format!(r#"{{"operation":"connected_triple_stitch_plan",{SOURCE},
            "through":{{"center":[-1.15,0],"radius":0.8}},
            "deeper_blind":{{"center":[0,0],"radius":0.8,"depth":2.8}},
            "shallower_blind":{{"center":[-0.4,0],"radius":0.8,"depth":1.2}}
        }}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"held");
        assert!(d["predicted_remaining_volume"].is_null());
        assert!(d["floor_interfaces"].as_array().unwrap().is_empty());
        assert_eq!(d["can_apply"],false);
    }
    #[test]
    fn invalid_source_rejects_without_partial_result(){
        let req=format!(r#"{{"operation":"connected_triple_stitch_plan",{SOURCE},
            "through":{{"center":[-1.15,0],"radius":0.8}},
            "deeper_blind":{{"center":[0,0],"radius":0.8,"depth":2.8}},
            "shallower_blind":{{"center":[1.15,0],"radius":0.8,"depth":5}}
        }}"#);
        let r=process(&req);
        assert!(!r.ok);
        assert!(r.data.is_none()&&r.solid.is_none());
    }
}
