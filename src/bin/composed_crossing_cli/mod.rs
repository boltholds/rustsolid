//! Statless geometry.v1 adapter for source-certified composed crossing clusters.
//! Multiple separated source-bound analytical components share ONE closed shell;
//! connected triples or inter-cluster touching remain unsupported.
use super::*;
use rustsolid::{
    CertifiedCrossingPair,ComposedCrossingBrep,
    ComposedCrossingParameters,ComposedCrossingDisplay,
};

#[derive(Debug,Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct PairInput{
    through:mixed_regions_cli::ThroughHoleInput,
    blind:BlindHoleInput,
}
fn pair(p:PairInput)->CertifiedCrossingPair{
    CertifiedCrossingPair{
        through:HoleParameters{center_x:p.through.center[0],
            center_z:p.through.center[1],radius:p.through.radius},
        blind:BlindHoleParameters{center_x:p.blind.center[0],
            center_z:p.blind.center[1],radius:p.blind.radius,depth:p.blind.depth}
    }
}
fn blind(h:BlindHoleInput)->BlindHoleParameters{
    BlindHoleParameters{center_x:h.center[0],center_z:h.center[1],
        radius:h.radius,depth:h.depth}
}
#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum ComposedEditInput {
    ReplacePair {index:usize,through:mixed_regions_cli::ThroughHoleInput,blind:BlindHoleInput},
    ReplaceIsolatedBlind {index:usize,center:[f64;2],radius:f64,depth:f64},
}
pub(super) fn evaluate(
    source:SolidSource,pairs:Vec<PairInput>,isolated_blind:Vec<BlindHoleInput>,
    batches:Vec<Vec<ComposedEditInput>>,undo:usize,redo:usize,
    display_segments:Option<usize>,
)->Response{
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if batches.len()>64||undo>64||redo>64||
            batches.iter().any(|b|b.is_empty()||b.len()>64){
            return Err(GeometryError::InvalidEdit(
                "composed source edits exceed the bounded batch/history budget".into()));
        }
        let box_=source.build()?;
        let parameters=ComposedCrossingParameters{
            pairs:pairs.into_iter().map(pair).collect(),
            isolated_blind:isolated_blind.into_iter().map(blind).collect(),
        };
        let mut body=ComposedCrossingBrep::from_box(&box_,parameters)?;
        let mut commits=Vec::new();
        for steps in batches{
            let len=steps.len();
            let mut tx=body.begin_edit();
            for item in steps {
                match item {
                    ComposedEditInput::ReplacePair{index,through,blind}=>{
                        tx.replace_pair(index,pair(PairInput{through,blind}))?;
                    }
                    ComposedEditInput::ReplaceIsolatedBlind{index,center,radius,depth}=>{
                        tx.replace_isolated_blind(index,BlindHoleParameters{
                            center_x:center[0],center_z:center[1],radius,depth
                        })?;
                    }
                }
            }
            let rev=tx.commit()?;
            commits.push(json!({"revision":rev,"command_count":len}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            if !body.undo()? {
                return Err(GeometryError::InvalidEdit(
                    "composed source undo exceeds stored history".into()));
            }
            transitions.push(json!({"kind":"undo","revision":body.revision()}));
        }
        for _ in 0..redo {
            if !body.redo()? {
                return Err(GeometryError::InvalidEdit(
                    "composed source redo exceeds stored history".into()));
            }
            transitions.push(json!({"kind":"redo","revision":body.revision()}));
        }
        body.validate()?;
        let m=body.model();
        let st=m.summary();
        let preview=display_segments.map(|n|body.tessellate(n)).transpose()?;
        let quality=preview.as_ref().map(|d:&ComposedCrossingDisplay|json!({
            "closed_oriented_manifold_certified":true,
            "genus":d.genus,
            "segments_per_analytic_arc":d.segments_per_source_arc,
            "triangle_edges":d.triangle_edges,
            "triangle_count":d.mesh.triangles.len(),
            "facet_volume":d.facet_volume,
            "facet_area":d.facet_surface_area,
            "analytic_volume":d.analytic_volume,
            "volume_faceting_error":d.volume_faceting_error,
            "independent_polygon_volume":d.expected_polygonal_volume
        }));
        let solid=preview.map(|display|SolidResponse{
            vertices:display.mesh.vertices.iter().map(|p|p.array()).collect(),
            faces:display.mesh.triangles,
            triangle_face_ids:display.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{
                vertices:st.vertices,edges:st.edges,coedges:st.coedges,
                loops:st.loops,shells:st.shells,faces:st.faces,
                euler_characteristic:st.euler_characteristic
            },
            mass:Mass{volume:m.mass.volume,
                surface_area:m.mass.surface_area,
                centroid:m.mass.centroid.array()},
            bbox:[m.bbox.min.array(),m.bbox.max.array()],
            edit_report:None,history:None,
            analytic_brep:Some(json!({
                "representation":"composed_crossing_shell",
                "source_certified":true,
                "genus":st.genus,
                "inner_loops":2*body.parameters().pairs.len()+
                    body.parameters().isolated_blind.len(),
                "mesh_is_display_only":true
            })),
            brep:None,
        });
        let sl=body.slab_plan();
        let data=json!({
            "kind":"composed_crossing_shell",
            "representation":"analytic_composed_crossing_shell",
            "state":"committed",
            "topology_committed":true,
            "generalized_topology_delta_can_apply":false,
            "source_pair_count":body.parameters().pairs.len(),
            "isolated_blind_count":body.parameters().isolated_blind.len(),
            "source_pairs":body.parameters().pairs.iter().map(|p|json!({
                "through":{"center":[p.through.center_x,p.through.center_z],
                            "radius":p.through.radius},
                "blind":{"center":[p.blind.center_x,p.blind.center_z],
                          "radius":p.blind.radius,"depth":p.blind.depth}
            })).collect::<Vec<_>>(),
            "topology":{"vertices":st.vertices,"edges":st.edges,
                "coedges":st.coedges,"loops":st.loops,"faces":st.faces,
                "shells":st.shells,"genus":st.genus,
                "euler_characteristic":st.euler_characteristic,
                "inner_loops":2*body.parameters().pairs.len()+
                    body.parameters().isolated_blind.len()},
            "mass":{"volume":m.mass.volume,
                "surface_area":m.mass.surface_area,
                "centroid":m.mass.centroid.array()},
            "source_slab_count":sl.slabs.len(),
            "slab_analytic_volume":sl.predicted_remaining_volume,
            "slabs":sl.slabs.iter().map(|s|json!({
                "y_lower":s.y_lower,"y_upper":s.y_upper,
                "active_disk_indices":s.active_disks,
                "boundary_components":s.connected_union_components,
                "planar_holes":s.planar_holes,
                "union_area":s.union_area
            })).collect::<Vec<_>>(),
            "revision":body.revision(),
            "undo_depth":body.undo_depth(),"redo_depth":body.redo_depth(),
            "commits":commits,"transitions":transitions,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":quality,
            "limitation":"Crossing pairs are physically separated from every other group. Connected triple intersections and arbitrary curved Booleans are not supported."
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
mod tests {
    use super::*;
    const SOURCE:&str=r#""source":{"kind":"block","origin":[-5,-1,-4],"width":10,"height":4,"depth":8}"#;
    const PAIRS:&str=r#""pairs":[
        {"through":{"center":[-2,0],"radius":0.45},
         "blind":{"center":[-1.55,0],"radius":0.45,"depth":1.5}},
        {"through":{"center":[1.5,0],"radius":0.5},
         "blind":{"center":[2,0],"radius":0.5,"depth":2.5}}
    ]"#;
    #[test]
    fn composed_crossings_produce_one_genus_two_source_shell(){
        let req=format!(r#"{{"operation":"compose_crossing_clusters",{SOURCE},{PAIRS},
            "display_segments":24}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        let data=r.data.unwrap();
        assert_eq!(data["source_pair_count"],2);
        assert_eq!(data["topology"]["genus"],2);
        assert_eq!(data["topology"]["coedges"],68);
        assert_eq!(data["source_slab_count"],3);
        assert_eq!(data["topology_committed"],true);
        assert_eq!(data["generalized_topology_delta_can_apply"],false);
        assert_eq!(data["display_quality"]["closed_oriented_manifold_certified"],true);
        assert!((data["slab_analytic_volume"].as_f64().unwrap()-
            data["mass"]["volume"].as_f64().unwrap()).abs()<1e-9);
        let mesh=r.solid.unwrap();
        assert_eq!(mesh.faces.len(),mesh.triangle_face_ids.len());
    }
    #[test]
    fn composed_source_undo_redo_is_atomic(){
        let req=format!(r#"{{"operation":"compose_crossing_clusters",{SOURCE},{PAIRS},
            "batches":[[
                {{"kind":"replace_pair","index":0,
                  "through":{{"center":[-2.1,0],"radius":0.45}},
                  "blind":{{"center":[-1.65,0],"radius":0.45,"depth":1.7}}}}
            ]],"undo":1,"redo":1}}"#);
        let r=process(&req);
        assert!(r.ok,"{:?}",r.error);
        let d=r.data.unwrap();
        assert_eq!(d["revision"],3);
        assert_eq!(d["transitions"].as_array().unwrap().len(),2);
        assert_eq!(d["topology"]["genus"],2);
    }
    #[test]
    fn connected_triplet_and_invalid_mesh_do_not_publish_partial_shape(){
        for extra in [
            r#""isolated_blind":[{"center":[-2,0.85],"radius":0.4,"depth":1}]"#,
            r#""display_segments":5"#,
        ]{
            let req=format!(r#"{{"operation":"compose_crossing_clusters",{SOURCE},{PAIRS},{extra}}}"#);
            let r=process(&req);
            assert!(!r.ok);
            assert!(r.solid.is_none()&&r.data.is_none());
        }
    }
}
