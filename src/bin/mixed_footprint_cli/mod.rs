//! Stateless, read-only JSON bridge for the independent mixed footprint plan.
use super::*;
use rustsolid::{MixedFootprintOverlapPlanner, MixedFootprintContact, FootprintSource};

pub(super) fn evaluate(
    source: SolidSource,
    through: mixed_regions_cli::ThroughHoleInput,
    blind: BlindHoleInput,
) -> Response {
    let result = (|| -> Result<Value, GeometryError> {
        let box_ = source.build()?;
        let t = HoleParameters {
            center_x: through.center[0], center_z: through.center[1],
            radius: through.radius,
        };
        let b = BlindHoleParameters {
            center_x: blind.center[0], center_z: blind.center[1],
            radius: blind.radius, depth: blind.depth,
        };
        let plan = MixedFootprintOverlapPlanner::plan(&box_, t, b)?;
        plan.validate_against(&box_,t,b)?;
        let contact = match plan.contact {
            MixedFootprintContact::Disjoint => "disjoint",
            MixedFootprintContact::TwoTransverseCrossings => "two_transverse_crossings",
            MixedFootprintContact::BlindInsideThrough => "blind_inside_through",
            MixedFootprintContact::ThroughInsideBlind => "through_inside_blind",
            MixedFootprintContact::Coincident => "coincident",
            MixedFootprintContact::HeldNearTangency => "held_near_tangency",
        };
        let points = plan.crossing_points.map(|array| {
            array.into_iter().map(Point3::array).collect::<Vec<_>>()
        });
        let arcs = plan.exposed_arcs.map(|array| {
            array.into_iter().map(|a| json!({
                "source": match a.source {
                    FootprintSource::Through => "through",
                    FootprintSource::Blind => "blind",
                },
                "start_crossing": a.start_crossing,
                "end_crossing": a.end_crossing,
                "start_angle": a.start_angle,
                "sweep_radians": a.sweep_radians,
                "midpoint": a.midpoint.array(),
            })).collect::<Vec<_>>()
        });
        Ok(json!({
            "kind":"mixed_footprint_overlap_plan",
            "state":contact,
            "through":{"center":[t.center_x,t.center_z],"radius":t.radius},
            "blind":{"center":[b.center_x,b.center_z],
                     "radius":b.radius,"depth":b.depth},
            "overlap_area":plan.overlap_area,
            "predicted_remaining_volume":plan.predicted_remaining_volume,
            "length_budget":plan.length_budget,
            "crossings":points,
            "exposed_arcs":arcs,
            "can_apply":false,
            "topology_committed":false,
            "source_revision":plan.source_revision,
            "note":"Read-only analytic prediction. No B-rep, Boolean commit or vendor-code dependency."
        }))
    })();
    match result {
        Ok(data) => Response {
            schema_version:"geometry.v1",ok:true,solid:None,error:None,data:Some(data),
        },
        Err(err) => Response {
            schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(err.to_string()),data:None,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn cli_predicts_crossing_volume_without_publishing_solid() {
        let r=process(r#"{"operation":"mixed_footprint_overlap_plan",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5}}"#);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"two_transverse_crossings");
        assert_eq!(d["topology_committed"],false);
        assert_eq!(d["can_apply"],false);
        assert_eq!(d["crossings"].as_array().unwrap().len(),2);
        assert_eq!(d["exposed_arcs"].as_array().unwrap().len(),2);
        assert!((d["predicted_remaining_volume"].as_f64().unwrap()-
            140.14094873829234).abs()<1e-9);
        assert!((d["overlap_area"].as_f64().unwrap()-
            0.3070924246521893).abs()<1e-10);
    }

    #[test]
    fn cli_near_tangent_returns_held_with_no_volume_or_geometry() {
        let r=process(r#"{"operation":"mixed_footprint_overlap_plan",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[1,0],"radius":0.5,"depth":1.5}}"#);
        assert!(r.ok,"{:?}",r.error);
        assert!(r.solid.is_none());
        let d=r.data.unwrap();
        assert_eq!(d["state"],"held_near_tangency");
        assert!(d["predicted_remaining_volume"].is_null());
        assert!(d["crossings"].is_null());
        assert_eq!(d["can_apply"],false);
    }

    #[test]
    fn cli_invalid_source_has_no_partial_result() {
        let r=process(r#"{"operation":"mixed_footprint_overlap_plan",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.5},
            "blind":{"center":[0.5,0],"radius":0.5,"depth":4.0}}"#);
        assert!(!r.ok);
        assert!(r.solid.is_none() && r.data.is_none());
    }

    #[test]
    fn cli_nested_containment_uses_distinct_depths() {
        let r=process(r#"{"operation":"mixed_footprint_overlap_plan",
            "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
            "through":{"center":[0,0],"radius":0.9},
            "blind":{"center":[0.1,0],"radius":0.3,"depth":2}}"#);
        assert!(r.ok,"{:?}",r.error);
        let d=r.data.unwrap();
        assert_eq!(d["state"],"blind_inside_through");
        assert!((d["predicted_remaining_volume"].as_f64().unwrap()-
            (144.0-PI*0.9*0.9*4.0)).abs()<1e-9);
        assert!(d["exposed_arcs"].is_null());
        assert_eq!(d["topology_committed"],false);
    }
}
