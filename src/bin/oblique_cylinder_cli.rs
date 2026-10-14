//! Opt-in geometry.v1 bridge for ONE source-certified cylinder vs global-Y
//! plane halfspace. The plane is horizontal in world coordinates but can be
//! oblique to a freely oriented cylinder axis. No arbitrary body-body Boolean.
use super::{Mass, Response, SolidResponse, Topology};
use rustsolid::{
    BrepModel, ClosedCylinderImprint, CylindricalBrep, FaceId, Frame3,
    GeometryError, GeometryTolerance, GraphState, ImprintResult, Mesh,
    ObliqueClipOutcome, ObliqueHalfspaceEditor, ObliqueHalfspaceSide,
    Point3, Solid, Surface3, TopologyDeltaPlanner, IntersectionGraphBuilder,
};
use serde::Deserialize;
use serde_json::{json, Value};

fn invalid(s: impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn p(v:[f64;3])->Point3{Point3{x:v[0],y:v[1],z:v[2]}}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn cross(a:Point3,b:Point3)->Point3{
    Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}
}
fn norm(a:Point3)->f64{a.x.hypot(a.y).hypot(a.z)}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct CylinderInput {
    origin:[f64;3],
    /// Cylinder's physical +V axis, not a world-space Euler angle.
    axis:[f64;3],
    /// Initial angular U=0 direction, perpendicular to 'axis'.
    radial_x:[f64;3],
    radius:f64,
    height:f64,
}
impl CylinderInput {
    fn build(self,tol:GeometryTolerance)->Result<CylindricalBrep,GeometryError>{
        let axis=p(self.axis);let radial=p(self.radial_x);
        let a=norm(axis);let r=norm(radial);
        if !a.is_finite() || !r.is_finite() ||
            a<=tol.absolute_length || r<=tol.absolute_length {
            return Err(invalid("cylinder axis and radial_x must be finite nonzero vectors"));
        }
        let n=Point3{x:axis.x/a,y:axis.y/a,z:axis.z/a};
        let x=Point3{x:radial.x/r,y:radial.y/r,z:radial.z/r};
        if dot(n,x).abs()>tol.angular {
            return Err(invalid("cylinder radial_x must be perpendicular to its axis"));
        }
        let y=cross(n,x);
        let frame=Frame3::from_axes(p(self.origin),x,y,tol)?;
        CylindricalBrep::new(frame,self.radius,self.height,tol)
    }
}

#[derive(Debug,Clone,Copy,Deserialize)]
#[serde(rename_all="snake_case")]
pub(super) enum KeepSide{Negative,Positive}
impl KeepSide{
    fn convert(self)->ObliqueHalfspaceSide{
        match self{Self::Negative=>ObliqueHalfspaceSide::Negative,
                   Self::Positive=>ObliqueHalfspaceSide::Positive}
    }
    fn as_str(self)->&'static str{
        match self{Self::Negative=>"negative",Self::Positive=>"positive"}
    }
}

fn failure(err:GeometryError)->Response{
    Response{schema_version:"geometry.v1",ok:false,solid:None,
        error:Some(err.to_string()),data:None}
}
fn held(stage:&'static str,reason:String)->Response{
    Response{schema_version:"geometry.v1",ok:false,solid:None,
        error:Some(format!("{stage} held: {reason}")),data:Some(json!({
            "kind":"oblique_cylinder_clip","state":"held","stage":stage,
            "reason":reason,"topology_committed":false,
            "general_body_boolean":false,
        }))}
}

fn actual_plane(source:&CylindricalBrep,height:f64)
    ->Result<(BrepModel,FaceId),GeometryError>
{
    if !height.is_finite(){return Err(invalid("plane_height must be finite"));}
    let bounds=source.bbox;
    let span=(bounds.max.x-bounds.min.x)
        .hypot(bounds.max.y-bounds.min.y)
        .hypot(bounds.max.z-bounds.min.z);
    if !span.is_finite(){return Err(invalid("source cylinder has invalid bounds"));}
    let pad=2.0*span.max(1.0);
    let min_x=bounds.min.x-pad;
    let min_z=bounds.min.z-pad;
    let width=bounds.max.x-bounds.min.x+2.0*pad;
    let depth=bounds.max.z-bounds.min.z+2.0*pad;
    if ![min_x,min_z,width,depth].iter().all(|v|v.is_finite()){
        return Err(invalid("finite plane face bounds overflow"));
    }
    let block=Solid::block_with_tolerance(
        Point3{x:min_x,y:height,z:min_z},
        width,pad,depth,source.tolerance,
    )?;
    let model=BrepModel::from_polyhedral(&block,source.units)?;
    let face=model.faces.iter().find_map(|face|{
        let Surface3::Plane(s)=model.face_surface(face.id)? else{return None};
        let eps=source.tolerance.length_at(span.max(1.0));
        ((s.frame.origin.y-height).abs()<=eps &&
         s.frame.normal.y.abs()>1.0-source.tolerance.angular).then_some(face.id)
    }).ok_or_else(||invalid("unable to select physical bounded plane face"))?;
    Ok((model,face))
}

fn face_kind(surface:Option<Surface3>)->&'static str{
    match surface {
        Some(Surface3::Plane(_))=>"plane",Some(Surface3::Cylinder(_))=>"cylinder",
        Some(Surface3::Sphere(_))=>"sphere",None=>"missing",
    }
}
fn edge_kind(curve:Option<rustsolid::Curve3>)->&'static str{
    match curve {
        Some(rustsolid::Curve3::Line(_))=>"line",
        Some(rustsolid::Curve3::Circle(_))=>"circle",
        Some(rustsolid::Curve3::Harmonic(_))=>"ellipse",
        None=>"missing",
    }
}
fn pcurve_kind(curve:Option<rustsolid::Curve2>)->&'static str{
    match curve {
        Some(rustsolid::Curve2::Line(_))=>"line",
        Some(rustsolid::Curve2::Circle(_)) |
        Some(rustsolid::Curve2::OrientedCircle(_))=>"circle",
        Some(rustsolid::Curve2::Harmonic(_))=>"harmonic",
        Some(rustsolid::Curve2::CylinderSection(_))=>"cylinder_section",
        None=>"missing",
    }
}
fn map_solid(model:&BrepModel,mesh:&Mesh,kind:&str,
    segments:usize,history:Value,quality:Value,
)->SolidResponse{
    let t=model.summary();
    let coedges=model.coedges.iter().map(|c|{
        let curve=model.geometry.coedges[c.id.0 as usize];
        json!({"id":c.id.0,"edge_id":c.edge.0,"face_id":c.face.0,
            "loop_id":c.loop_id.0,"twin_id":c.twin.0,
            "next_id":c.next.0,"prev_id":c.prev.0,"reversed":c.reversed,
            "pcurve_id":curve.pcurve.0,"pcurve_kind":pcurve_kind(model.coedge_pcurve(c.id))})
    }).collect::<Vec<_>>();
    SolidResponse{
        vertices:mesh.vertices.iter().map(|p|p.array()).collect(),
        faces:mesh.triangles.clone(),
        triangle_face_ids:mesh.triangle_faces.iter().map(|face|face.0).collect(),
        topology:Topology{
            vertices:t.vertices,edges:t.edges,coedges:t.coedges,
            loops:t.loops,shells:t.shells,faces:t.faces,
            euler_characteristic:t.euler_characteristic,
        },
        mass:Mass{volume:model.mass.volume,surface_area:model.mass.surface_area,
            centroid:model.mass.centroid.array()},
        bbox:[model.bbox.min.array(),model.bbox.max.array()],
        edit_report:None,
        history:Some(history),
        analytic_brep:Some(json!({
            "kind":kind,
            "analytic_mass":true,
            "display_mesh_only":true,
            "face_surface_kinds":model.faces.iter()
                .map(|f|face_kind(model.face_surface(f.id))).collect::<Vec<_>>(),
            "edge_carrier_kinds":model.edges.iter()
                .map(|e|edge_kind(model.edge_curve(e.id))).collect::<Vec<_>>(),
            "coedges":coedges,
            "new_elliptic_edge":if kind=="oblique_halfspace_clipped_cylinder" {Some(1)}else{Some(4)},
            "facet_segments":segments,
            "mesh_quality":quality,
        })),
        brep:None,
    }
}

pub(super) fn evaluate(
    cylinder:CylinderInput,plane_height:f64,retain:KeepSide,
    display_segments:usize,undo:usize,redo:usize,
    tolerance:GeometryTolerance,
)->Response{
    let result=(||->Result<Response,GeometryError>{
        tolerance.validate()?;
        if !(12..=4096).contains(&display_segments){
            return Err(invalid("display_segments must be in 12..=4096"));
        }
        // Exactly one bounded edit is performed; nonsensical or unbounded
        // histories are an error, never a quietly shortened replay.
        if undo>1 || redo>undo {
            return Err(invalid("one cut supports undo=0/1 and redo<=undo"));
        }
        let source=cylinder.build(tolerance)?;
        // Negative/positive sides refer to global Y for this endpoint;
        // a downward or almost-horizontal cylinder reverses / destabilizes
        // that convention and is unsupported here.
        if source.frame.normal.y<=tolerance.angular.sin() {
            return Err(invalid("cylinder axis must have a positive world-Y component"));
        }
        let original=BrepModel::from_cylinder(&source)?;
        let (plane,face)=actual_plane(&source,plane_height)?;
        let graph=IntersectionGraphBuilder::new(tolerance)?
            .intersect_faces(&plane,face,&original,FaceId(2))?;
        if graph.state!=GraphState::Ready{
            return Ok(held("intersection",format!("{:?}",graph.state)));
        }
        let plan=TopologyDeltaPlanner::new().plan(
            &graph,&plane,face,&original,FaceId(2))?;
        let mut imprint=ClosedCylinderImprint::new(source)?;
        match imprint.imprint(&graph,&plan,&plane,&original)?{
            ImprintResult::Committed(_)=>{},
            ImprintResult::Held(reason)=>return Ok(held("imprint",format!("{reason:?}"))),
        }
        let mut clip=ObliqueHalfspaceEditor::new(&imprint)?;
        let region=clip.classify()?;
        match clip.clip(retain.convert())? {
            ObliqueClipOutcome::Committed(_)=>{},
            ObliqueClipOutcome::Held(reason)=>return Ok(held("clip",format!("{reason:?}"))),
        }
        let mut transitions=vec![json!({"kind":"clip","revision":clip.revision()})];
        for _ in 0..undo{
            if !clip.undo()?{
                return Err(invalid("cut Undo exceeds stored history"));
            }
            transitions.push(json!({"kind":"undo","revision":clip.revision()}));
        }
        for _ in 0..redo{
            if !clip.redo()?{
                return Err(invalid("cut Redo exceeds stored history"));
            }
            transitions.push(json!({"kind":"redo","revision":clip.revision()}));
        }
        clip.model().validate()?;
        let applied=clip.applied().is_some();
        let (mesh,quality)=if applied {
            let display=clip.tessellate(display_segments)?;
            clip.validate_display(&display)?;
            (display.mesh,json!({
                "closed_oriented_manifold_certified":true,
                "mesh_euler":display.certificate.mesh_euler,
                "facet_volume":display.certificate.facet_volume,
                "source_polygon_volume":display.certificate.source_polygon_volume,
                "face_triangle_counts":display.certificate.triangle_faces,
                "source_revision":display.certificate.revision,
            }))
        } else {
            // Undo returns the imprinted full cylinder, not the original
            // three-face source and certainly not a falsely clipped body.
            let display=imprint.tessellate(display_segments)?;
            imprint.validate_display(&display)?;
            (display.mesh,json!({
                "closed_oriented_manifold_certified":true,
                "mesh_euler":display.certificate.mesh_euler,
                "facet_volume":display.certificate.facet_volume,
                "source_polygon_volume":display.certificate.inscribed_prism_volume,
                "face_triangle_counts":[display_segments,display_segments,
                    2*display_segments,2*display_segments],
                "source_revision":imprint.revision(),
            }))
        };
        let history=json!({
            "revision":clip.revision(),"undo_depth":clip.undo_depth(),
            "redo_depth":clip.redo_depth(),"transitions":transitions,
        });
        let state=if applied {"committed"}else{"imprinted_after_undo"};
        let kind=if applied {"oblique_halfspace_clipped_cylinder"}
            else {"imprinted_cylinder"};
        let solid=map_solid(clip.model(),&mesh,kind,display_segments,history.clone(),quality.clone());
        Ok(Response{schema_version:"geometry.v1",ok:true,solid:Some(solid),
            error:None,data:Some(json!({
                "kind":"oblique_cylinder_clip","state":state,
                "retained_halfspace":if applied {Some(retain.as_str())}else{None},
                "plane":{"kind":"world_y","height":plane_height},
                "source_kind":"analytic_cylinder","general_body_boolean":false,
                "topology_committed":applied,
                "analytic_cell_classification":{
                    "source_revision":region.source_revision,
                    "axial_mean":region.axial_mean,
                    "axial_amplitude":region.axial_amplitude,
                    "face_sides":["negative","positive","negative","positive"],
                    "elliptic_edge_id":region.cut_edge.0,
                },
                "history":history,"mesh_certificate":quality,
            }))})
    })();
    match result {Ok(value)=>value,Err(e)=>failure(e)}
}
