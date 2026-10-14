//! Source-certified bounded Boolean: retain one side of the oblique
//! transverse section through a capped analytical cylinder. The halfspace
//! itself is mathematically unbounded. This is NOT arbitrary body/body CSG.
//! The cut is constructed from first-party harmonic surface geometry, not
//! a sampled Boolean or proprietary geometric-kernel implementation.
use std::f64::consts::{PI,TAU};
use crate::{
    BrepFace,BrepModel,BrepOrigin,ClosedCylinderImprint,Coedge,CoedgeGeometry,
    CoedgeId,Curve2,Curve2Id,Curve3,Curve3Id,CylindricalBrep,
    Edge,EdgeGeometry,EdgeId,FaceGeometry,FaceId,FaceSense,Frame3,
    GeometryError,HarmonicCurve2,HarmonicCurve3,Line2,Loop,LoopId,
    LoopRole,MassProperties,ParameterDomain2,ParameterRange,PlaneSurface,
    Point2Param,Point3,Shell,ShellId,Surface3,SurfaceId,Vertex,VertexId,
    BoundingBox,
};

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
fn n(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn dist(a:Point3,b:Point3)->f64{n(a.sub(b))}
fn bounded(a:f64,b:f64)->ParameterRange{
    ParameterRange::Bounded{start:a,end:b}
}
fn bounded_uv(umax:f64,vmax:f64)->ParameterDomain2{
    ParameterDomain2{
        u:bounded(-umax,umax),
        v:bounded(-vmax,vmax),
        u_periodic:false,v_periodic:false,
    }
}
fn range_extrema(a:Point3,b:Point3,c:Point3)->(Point3,Point3){
    let amplitude=|x:f64,y:f64|x.hypot(y);
    let radius=Point3{x:amplitude(b.x,c.x),
        y:amplitude(b.y,c.y),z:amplitude(b.z,c.z)};
    (a.sub(radius),a.add(radius))
}
fn envelope(a:(Point3,Point3),b:(Point3,Point3))->BoundingBox{
    BoundingBox{
        min:Point3{x:a.0.x.min(b.0.x),y:a.0.y.min(b.0.y),
            z:a.0.z.min(b.0.z)},
        max:Point3{x:a.1.x.max(b.1.x),y:a.1.y.max(b.1.y),
            z:a.1.z.max(b.1.z)},
    }
}

/// A strict planar halfspace relative to the source ellipse orientation.
/// Negative retains material below the section along the source axial frame;
/// Positive retains material above it.
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum ObliqueHalfspaceSide {Negative,Positive}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum ObliqueCutHold {AlreadyClipped,IncompleteImprint,NearDegenerateSection}

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ObliqueRegionClassification {
    pub source_revision:u64,
    pub cut_edge:EdgeId,
    pub elliptic_center:Point3,
    pub upward_normal:Point3,
    pub axial_mean:f64,
    pub axial_amplitude:f64,
    /// [original lower circle, upper circle, lower side, upper side]
    pub face_sides:[ObliqueHalfspaceSide;4],
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ObliqueClipCommit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub retained_side:ObliqueHalfspaceSide,
    pub analytic_mass:MassProperties,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum ObliqueClipOutcome {Committed(ObliqueClipCommit), Held(ObliqueCutHold)}

#[derive(Debug,Clone,Copy)]
struct Section {
    ellipse:HarmonicCurve3,
    mean:f64,
    cos:f64,
    sin:f64,
    amp:f64,
    oriented_normal:Point3,
}

fn section_info(source:&CylindricalBrep,imprint:&BrepModel)->Result<Section,GeometryError>{
    imprint.validate()?;
    let Curve3::Harmonic(ellipse)=imprint.edge_curve(EdgeId(4))
        .ok_or_else(||invalid("imprint lacks the physical elliptical cut edge"))?
        else{return Err(invalid("imprinted edge has no harmonic carrier"));};
    let a=source.frame.normal;
    let mean=ellipse.center.sub(source.frame.origin).dot(a);
    let cos=ellipse.cosine.dot(a);
    let sin=ellipse.sine.dot(a);
    let span=cos.hypot(sin);
    let source_span=source.height.max(source.radius);
    let eps=source.tolerance.length_at(source_span)*64.0;
    let expected_center=source.frame.origin.add(a.scale(mean));
    let expected_cos=source.frame.x.scale(source.radius).add(a.scale(cos));
    let expected_sin=source.frame.y.scale(source.radius).add(a.scale(sin));
    if dist(ellipse.center,expected_center)>eps ||
        n(ellipse.cosine.sub(expected_cos))>eps ||
        n(ellipse.sine.sub(expected_sin))>eps ||
        ![mean,cos,sin,span].iter().all(|x|x.is_finite()) ||
        mean-span<=eps || mean+span>=source.height-eps {
        return Err(invalid("source harmonic cut does not separate two strict regions"));
    }
    let raw=ellipse.cosine.cross(ellipse.sine);
    let len=n(raw);
    if !len.is_finite() || len<=source.radius*eps {
        return Err(invalid("oblique cut plane is singular"));
    }
    let oriented_normal=raw.scale(1.0/len);
    if oriented_normal.dot(a)<=source.tolerance.angular {
        return Err(invalid("oblique cut has inconsistent axial orientation"));
    }
    Ok(Section{ellipse,mean,cos,sin,amp:span,oriented_normal})
}

/// Precisely bounded analytic face-cell classification, not point-sampling
/// heuristics: every lower cell lies on the negative side, every upper
/// cell on the positive side because 0<v_section(theta)<height holds
/// for every theta, and the plane equation is affine along each ruling.
pub fn classify_imprinted_cylinder(imprint:&ClosedCylinderImprint)
    ->Result<ObliqueRegionClassification,GeometryError>
{
    let model=imprint.model();
    if model.vertices.len()!=3 || model.edges.len()!=5 ||
        model.faces.len()!=4 || model.loops.len()!=4 ||
        model.coedges.len()!=10 {
        return Err(invalid("oblique halfspace classification requires one closed imprint"));
    }
    let section=section_info(imprint.source_cylinder(),model)?;
    Ok(ObliqueRegionClassification{
        source_revision:imprint.revision(),
        cut_edge:EdgeId(4),
        elliptic_center:section.ellipse.center,
        upward_normal:section.oriented_normal,
        axial_mean:section.mean,
        axial_amplitude:section.amp,
        face_sides:[
            ObliqueHalfspaceSide::Negative,
            ObliqueHalfspaceSide::Positive,
            ObliqueHalfspaceSide::Negative,
            ObliqueHalfspaceSide::Positive,
        ],
    })
}

fn model_mass(source:&CylindricalBrep,section:Section,side:ObliqueHalfspaceSide)
    ->Result<MassProperties,GeometryError>
{
    let h=source.height;
    let radius=source.radius;
    let circle_area=PI*radius*radius;
    let slope_elliptical_area=n(section.ellipse.cosine.cross(section.ellipse.sine))*PI;
    let (axial_height,radial_direction,axial_centroid)=match side {
        ObliqueHalfspaceSide::Negative=>(
            section.mean,1.0,
            section.mean*0.5+section.amp*section.amp/(8.0*section.mean)),
        ObliqueHalfspaceSide::Positive=>(
            h-section.mean,-1.0,
            (h+section.mean)*0.5-
                section.amp*section.amp/(8.0*(h-section.mean))),
    };
    let radial=source.frame.x.scale(radial_direction*section.cos*radius/(4.0*axial_height))
        .add(source.frame.y.scale(radial_direction*section.sin*radius/(4.0*axial_height)));
    let centroid=source.frame.origin.add(source.frame.normal.scale(axial_centroid)).add(radial);
    let mass=MassProperties{
        volume:circle_area*axial_height,
        surface_area:circle_area+slope_elliptical_area+
            TAU*radius*axial_height,
        centroid,
    };
    if mass.volume<=0.0 || !mass.volume.is_finite() ||
       mass.surface_area<=0.0 || !mass.surface_area.is_finite() ||
       !mass.centroid.is_finite() {
        return Err(invalid("oblique cut analytic mass became invalid"));
    }
    Ok(mass)
}

fn split_bbox(source:&CylindricalBrep,section:Section,side:ObliqueHalfspaceSide)->BoundingBox{
    let axis=source.frame.normal;
    let r=source.radius;
    let ellipse=range_extrema(section.ellipse.center,
        section.ellipse.cosine,section.ellipse.sine);
    let original_base=match side{
        ObliqueHalfspaceSide::Negative=>source.frame.origin,
        ObliqueHalfspaceSide::Positive=>
            source.frame.origin.add(axis.scale(source.height)),
    };
    let horizontal=range_extrema(original_base,
        source.frame.x.scale(r),source.frame.y.scale(r));
    envelope(horizontal,ellipse)
}

fn circle_uv(model:&BrepModel,id:CoedgeId,face:u32,loop_index:u32,
    new_curves:&mut Vec<Curve2>,
)->Result<CoedgeGeometry,GeometryError>{
    let curve=model.coedge_pcurve(id).ok_or_else(||invalid("source coedge has no trim"))?;
    let id_new=Curve2Id(u32::try_from(new_curves.len())
        .map_err(|_|invalid("UV catalog too large"))?);
    new_curves.push(curve);
    let domain=model.geometry.coedges[id.0 as usize].domain;
    Ok(CoedgeGeometry{pcurve:id_new,domain,
        face:FaceId(face),loop_id:LoopId(loop_index)})
}
fn new_curve(trim:Curve2,range:ParameterRange,face:u32,loop_index:u32,
    curves:&mut Vec<Curve2>)->Result<CoedgeGeometry,GeometryError>{
    let id=Curve2Id(u32::try_from(curves.len())
        .map_err(|_|invalid("UV catalog too large"))?);
    curves.push(trim);
    Ok(CoedgeGeometry{pcurve:id,domain:range,
        face:FaceId(face),loop_id:LoopId(loop_index)})
}

fn build(source:&CylindricalBrep,imprint:&BrepModel,section:Section,
    retained:ObliqueHalfspaceSide,
)->Result<BrepModel,GeometryError>{
    let below=retained==ObliqueHalfspaceSide::Negative;
    let r=source.radius;
    let h=source.height;
    let cut_norm=if below {section.oriented_normal}
        else{section.oriented_normal.scale(-1.0)};
    let plane_x=section.ellipse.cosine.scale(1.0/n(section.ellipse.cosine));
    let plane_y=cut_norm.cross(plane_x);
    let frame=Frame3::from_axes(section.ellipse.center,plane_x,plane_y,source.tolerance)?;
    let face_cut=Surface3::Plane(PlaneSurface{frame});
    let uv_center=frame.project(section.ellipse.center)?;
    let uv_cos=frame.project(section.ellipse.center.add(section.ellipse.cosine))?;
    let uv_sin=frame.project(section.ellipse.center.add(section.ellipse.sine))?;
    let reverse=if below {1.0}else{-1.0};
    let ellipse_trim=Curve2::Harmonic(HarmonicCurve2{
        center:uv_center,
        cosine:Point2Param{u:uv_cos.u-uv_center.u,v:uv_cos.v-uv_center.v},
        sine:Point2Param{
            u:reverse*(uv_sin.u-uv_center.u),
            v:reverse*(uv_sin.v-uv_center.v)},
    });
    // Use exact harmonic coefficient amplitudes, independent of the
    // shifted-phase chosen samples.
    let uv_extent_u=(uv_cos.u-uv_center.u).hypot(uv_sin.u-uv_center.u);
    let uv_extent_v=(uv_cos.v-uv_center.v).hypot(uv_sin.v-uv_center.v);
    if uv_extent_u<=source.tolerance.absolute_length ||
        uv_extent_v<=source.tolerance.absolute_length{
        return Err(invalid("elliptic cap has degenerate planar UV span"));
    }
    let (endpoint_a,endpoint_b,original_circle_coedges,original_seam,
        original_ellipse,original_side,original_circle_edge,original_seam_edge)=
    if below {
        (VertexId(0),VertexId(2),[CoedgeId(0),CoedgeId(1)],
            [CoedgeId(4),CoedgeId(5)],CoedgeId(9),FaceId(2),
            EdgeId(0),EdgeId(2))
    }else{
        (VertexId(2),VertexId(1),[CoedgeId(2),CoedgeId(3)],
            [CoedgeId(6),CoedgeId(7)],CoedgeId(8),FaceId(3),
            EdgeId(1),EdgeId(3))
    };
    let vertex_a=imprint.vertices[endpoint_a.0 as usize].position;
    let vertex_b=imprint.vertices[endpoint_b.0 as usize].position;
    let circle=imprint.edge_curve(original_circle_edge)
        .ok_or_else(||invalid("missing retained rim curve"))?;
    let ellipse=imprint.edge_curve(EdgeId(4))
        .ok_or_else(||invalid("missing elliptic rim curve"))?;
    let seam=imprint.edge_curve(original_seam_edge)
        .ok_or_else(||invalid("missing retained seam curve"))?;
    let circle_geometry=imprint.geometry.edges[original_circle_edge.0 as usize];
    let ellipse_geometry=imprint.geometry.edges[4];
    let seam_geometry=imprint.geometry.edges[original_seam_edge.0 as usize];
    let edges=vec![
        Edge{id:EdgeId(0),start:VertexId(0),end:VertexId(0),
            coedges:[CoedgeId(0),CoedgeId(1)]},
        Edge{id:EdgeId(1),start:VertexId(1),end:VertexId(1),
            coedges:[CoedgeId(2),CoedgeId(3)]},
        Edge{id:EdgeId(2),start:VertexId(0),end:VertexId(1),
            coedges:[CoedgeId(4),CoedgeId(5)]},
    ];
    let vertices=vec![
        Vertex{id:VertexId(0),position:vertex_a},
        Vertex{id:VertexId(1),position:vertex_b},
    ];
    let coedges=vec![
        Coedge{id:CoedgeId(0),edge:EdgeId(0),face:FaceId(2),loop_id:LoopId(2),
            reversed:false,next:CoedgeId(4),prev:CoedgeId(5),twin:CoedgeId(1)},
        Coedge{id:CoedgeId(1),edge:EdgeId(0),face:FaceId(0),
            loop_id:LoopId(0),reversed:true,
            next:CoedgeId(1),prev:CoedgeId(1),twin:CoedgeId(0)},
        Coedge{id:CoedgeId(2),edge:EdgeId(1),face:FaceId(2),loop_id:LoopId(2),
            reversed:true,next:CoedgeId(5),prev:CoedgeId(4),twin:CoedgeId(3)},
        Coedge{id:CoedgeId(3),edge:EdgeId(1),face:FaceId(1),
            loop_id:LoopId(1),reversed:false,
            next:CoedgeId(3),prev:CoedgeId(3),twin:CoedgeId(2)},
        Coedge{id:CoedgeId(4),edge:EdgeId(2),face:FaceId(2),loop_id:LoopId(2),
            reversed:false,next:CoedgeId(2),prev:CoedgeId(0),twin:CoedgeId(5)},
        Coedge{id:CoedgeId(5),edge:EdgeId(2),face:FaceId(2),loop_id:LoopId(2),
            reversed:true,next:CoedgeId(0),prev:CoedgeId(2),twin:CoedgeId(4)},
    ];
    // Below: lower circle->seam->ellipse(rev)->seam.
    // Above: ellipse->seam->upper circle(rev)->seam.
    let loops=vec![
        Loop{id:LoopId(0),face:FaceId(0),role:LoopRole::Outer,
            first_coedge:CoedgeId(1)},
        Loop{id:LoopId(1),face:FaceId(1),role:LoopRole::Outer,
            first_coedge:CoedgeId(3)},
        Loop{id:LoopId(2),face:FaceId(2),role:LoopRole::Outer,
            first_coedge:CoedgeId(0)},
    ];
    let surfaces=if below {
        vec![imprint.face_surface(FaceId(0))
                .ok_or_else(||invalid("missing lower planar cap support"))?,
            face_cut,
            imprint.face_surface(FaceId(2))
                .ok_or_else(||invalid("missing cylindrical side support"))?]
    }else{
        vec![face_cut,
            imprint.face_surface(FaceId(1))
                .ok_or_else(||invalid("missing upper planar cap support"))?,
            imprint.face_surface(FaceId(3))
                .ok_or_else(||invalid("missing upper cylindrical support"))?]
    };
    let circ_bound=bounded(-r,r);
    let side_bound=imprint.faces[original_side.0 as usize].geometry.bounds;
    let face_geometry=|s:u32,lp:u32,bounds:ParameterDomain2|FaceGeometry{
        surface:SurfaceId(s),bounds,outer_loop:LoopId(lp),tolerance:source.tolerance,
    };
    let faces=vec![
        BrepFace{id:FaceId(0),loops:vec![LoopId(0)],
            geometry:face_geometry(0,0,if below {
                ParameterDomain2{u:circ_bound,v:circ_bound,
                    u_periodic:false,v_periodic:false}
            }else{bounded_uv(uv_extent_u,uv_extent_v)}),
            sense:FaceSense::Forward},
        BrepFace{id:FaceId(1),loops:vec![LoopId(1)],
            geometry:face_geometry(1,1,if below {
                bounded_uv(uv_extent_u,uv_extent_v)
            }else{
                ParameterDomain2{u:circ_bound,v:circ_bound,
                    u_periodic:false,v_periodic:false}}),
            sense:FaceSense::Forward},
        BrepFace{id:FaceId(2),loops:vec![LoopId(2)],
            geometry:face_geometry(2,2,side_bound),sense:FaceSense::Forward},
    ];
    let mut curves2=Vec::with_capacity(6);
    let mut coedge_geom=Vec::with_capacity(6);
    // The owned pcurves must follow the *new* six coedge IDs in order:
    // lower = circle side, circle cap, ellipse side, ellipse cap, seams.
    // upper = ellipse side, ellipse cap, circle side, circle cap, seams.
    if below {
        coedge_geom.push(circle_uv(imprint,original_circle_coedges[0],2,2,&mut curves2)?);
        coedge_geom.push(circle_uv(imprint,original_circle_coedges[1],0,0,&mut curves2)?);
        coedge_geom.push(circle_uv(imprint,original_ellipse,2,2,&mut curves2)?);
        coedge_geom.push(new_curve(ellipse_trim,bounded(0.0,TAU),1,1,&mut curves2)?);
    }else{
        coedge_geom.push(circle_uv(imprint,original_ellipse,2,2,&mut curves2)?);
        coedge_geom.push(new_curve(ellipse_trim,bounded(0.0,TAU),0,0,&mut curves2)?);
        coedge_geom.push(circle_uv(imprint,original_circle_coedges[0],2,2,&mut curves2)?);
        coedge_geom.push(circle_uv(imprint,original_circle_coedges[1],1,1,&mut curves2)?);
    }
    coedge_geom.push(circle_uv(imprint,original_seam[0],2,2,&mut curves2)?);
    coedge_geom.push(circle_uv(imprint,original_seam[1],2,2,&mut curves2)?);
    if coedge_geom.len()!=6 {return Err(invalid("cut failed to create 6 coedge trims"));}
    let curves3=if below {vec![circle,ellipse,seam]}else{vec![ellipse,circle,seam]};
    let edge_geom=if below {
        vec![circle_geometry,ellipse_geometry,seam_geometry]
    }else{vec![ellipse_geometry,circle_geometry,seam_geometry]};
    // Deduplicated 3D curve catalog: each physical edge has a distinct
    // source curve, while both coedges reference that *same* physical edge.
    let edge_geom=edge_geom.into_iter().enumerate().map(|(i,g)|
        EdgeGeometry{curve:Curve3Id(i as u32),domain:g.domain,tolerance:g.tolerance}
    ).collect();
    let body=BrepModel{
        origin:BrepOrigin::AnalyticObliqueClip,
        vertices,edges,coedges,loops,faces,
        shells:vec![Shell{id:ShellId(0),
            faces:vec![FaceId(0),FaceId(1),FaceId(2)],closed:true}],
        geometry:crate::BrepGeometry{
            curves3,curves2,surfaces,edges:edge_geom,coedges:coedge_geom},
        mass:model_mass(source,section,retained)?,
        bbox:split_bbox(source,section,retained),
        tolerance:source.tolerance,
        units:source.units,
    };
    body.validate()?;
    let expected=PI*source.radius.powi(2)*if below{
        section.mean
    }else{source.height-section.mean};
    if (body.mass.volume-expected).abs()>source.tolerance.length_at(source.height)
        *source.radius.powi(2)*32.0 {
        return Err(invalid("oblique halfspace volume failed analytic disk integral"));
    }
    Ok(body)
}

/// Atomic owner of a *single restricted volume-changing clip*, initialized
/// from a previously source-certified and closed cylinder imprint.
pub struct ObliqueHalfspaceEditor {
    source:CylindricalBrep,
    original:BrepModel,
    current:BrepModel,
    section:Section,
    imprint_revision:u64,
    revision:u64,
    applied:Option<ObliqueHalfspaceSide>,
    undo:Vec<(BrepModel,Option<ObliqueHalfspaceSide>)>,
    redo:Vec<(BrepModel,Option<ObliqueHalfspaceSide>)>,
}
impl ObliqueHalfspaceEditor {
    pub fn new(imprint:&ClosedCylinderImprint)->Result<Self,GeometryError>{
        let classification=classify_imprinted_cylinder(imprint)?;
        let source=imprint.source_cylinder().clone();
        let original=imprint.model().clone();
        let section=section_info(&source,&original)?;
        if classification.cut_edge!=EdgeId(4) {
            return Err(invalid("imprint region classifier source changed"));
        }
        Ok(Self{
            source,original:original.clone(),current:original,
            section,imprint_revision:imprint.revision(),
            revision:0,applied:None,undo:vec![],redo:vec![],
        })
    }
    pub fn model(&self)->&BrepModel{&self.current}
    pub fn source(&self)->&CylindricalBrep{&self.source}
    pub fn revision(&self)->u64{self.revision}
    pub fn applied(&self)->Option<ObliqueHalfspaceSide>{self.applied}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}

    pub fn classify(&self)->Result<ObliqueRegionClassification,GeometryError>{
        let oriented_normal=self.section.oriented_normal;
        self.original.validate()?;
        Ok(ObliqueRegionClassification{
            source_revision:self.imprint_revision,cut_edge:EdgeId(4),
            elliptic_center:self.section.ellipse.center,
            upward_normal:oriented_normal,axial_mean:self.section.mean,
            axial_amplitude:self.section.amp,
            face_sides:[
                ObliqueHalfspaceSide::Negative,ObliqueHalfspaceSide::Positive,
                ObliqueHalfspaceSide::Negative,ObliqueHalfspaceSide::Positive,
            ],
        })
    }

    pub fn clip(&mut self,retain:ObliqueHalfspaceSide)
        ->Result<ObliqueClipOutcome,GeometryError>
    {
        if self.applied.is_some(){
            return Ok(ObliqueClipOutcome::Held(ObliqueCutHold::AlreadyClipped));
        }
        self.original.validate()?;
        self.current.validate()?;
        // All writes are scoped to a private candidate. An error from
        // B-rep topology or mass/orientation validation discards it.
        let candidate=build(&self.source,&self.original,self.section,retain)?;
        // Independent polygonal watertightness, ownership and finite-volume
        // tests are mandatory *before* the candidate becomes observable.
        crate::oblique_clip_mesh::preflight_clip_display(
            &self.source,&candidate,retain)?;
        let before=self.revision;
        let after=before.checked_add(1).ok_or_else(||invalid("oblique clip revision exhausted"))?;
        let old=std::mem::replace(&mut self.current,candidate);
        if self.undo.len()>=128{self.undo.remove(0);}
        self.undo.push((old,self.applied));
        self.redo.clear();
        self.applied=Some(retain);
        self.revision=after;
        Ok(ObliqueClipOutcome::Committed(ObliqueClipCommit{
            revision_before:before,revision_after:after,
            retained_side:retain,analytic_mass:self.current.mass,
        }))
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{
        let Some((candidate,selection))=self.undo.last().cloned() else{return Ok(false)};
        candidate.validate()?;
        let rev=self.revision.checked_add(1).ok_or_else(||invalid("oblique cut revision overflow"))?;
        let old=std::mem::replace(&mut self.current,candidate);
        let old_selection=self.applied;
        self.undo.pop();self.redo.push((old,old_selection));
        self.applied=selection;self.revision=rev;
        Ok(true)
    }
    pub fn redo(&mut self)->Result<bool,GeometryError>{
        let Some((candidate,selection))=self.redo.last().cloned() else{return Ok(false)};
        candidate.validate()?;
        let rev=self.revision.checked_add(1).ok_or_else(||invalid("oblique cut revision overflow"))?;
        let old=std::mem::replace(&mut self.current,candidate);
        let old_selection=self.applied;
        self.redo.pop();if self.undo.len()>=128{self.undo.remove(0);}
        self.undo.push((old,old_selection));
        self.applied=selection;self.revision=rev;
        Ok(true)
    }
}
