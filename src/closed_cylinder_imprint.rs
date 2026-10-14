//! Narrow, source-preserving execution of one closed oblique intersection
//! imprint on a capped analytic cylinder. This changes the *partitioning* of
//! the lateral surface, not the material or mass of the original cylinder.
//! A generic Boolean TopologyDeltaPlan remains preview-only.
use std::f64::consts::TAU;

use crate::{
    BrepModel, BrepOrigin, Coedge, CoedgeGeometry,
    CoedgeId, Curve2, Curve2Id, Curve3, Curve3Id, CylinderSection2,
    CylindricalBrep, Edge, EdgeGeometry, EdgeId, FaceGeometry, FaceId,
    FaceSense, GeometryError, GraphNodeKind, GraphState,
    IntersectionCurve, IntersectionGraph, Line2, Loop, LoopId, LoopRole,
    ParameterDomain2, ParameterRange, PlanArcRole, PlanState, Point2Param,
    Surface3, TopologyDeltaPlan, Vertex, VertexId,
};

const MAX_HISTORY:usize=128;

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn uv(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn line(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:uv(b.u-a.u,b.v-a.v)})
}
fn domain(lo:f64,hi:f64)->ParameterRange{
    ParameterRange::Bounded{start:lo,end:hi}
}
fn length(p:crate::Point3)->f64 {p.x.hypot(p.y).hypot(p.z)}
fn distance(a:crate::Point3,b:crate::Point3)->f64 {length(a.sub(b))}

/// Exact structural identity of the bounded canonical cylinder operand:
/// compare actual topology and geometry fields instead of serializing
/// Debug output or treating matching dimensions as a trusted source.
fn same_cylinder_source(a:&BrepModel,b:&BrepModel)->bool{
    a.origin==b.origin && a.units==b.units &&
    a.tolerance==b.tolerance && a.mass==b.mass && a.bbox==b.bbox &&
    a.geometry.curves3==b.geometry.curves3 &&
    a.geometry.curves2==b.geometry.curves2 &&
    a.geometry.surfaces==b.geometry.surfaces &&
    a.geometry.edges==b.geometry.edges &&
    a.geometry.coedges==b.geometry.coedges &&
    a.vertices.len()==b.vertices.len() &&
    a.vertices.iter().zip(&b.vertices).all(|(x,y)|
        x.id==y.id && x.position==y.position) &&
    a.edges.len()==b.edges.len() &&
    a.edges.iter().zip(&b.edges).all(|(x,y)|
        x.id==y.id && x.start==y.start && x.end==y.end && x.coedges==y.coedges) &&
    a.coedges.len()==b.coedges.len() &&
    a.coedges.iter().zip(&b.coedges).all(|(x,y)|
        x.id==y.id && x.edge==y.edge && x.face==y.face &&
        x.loop_id==y.loop_id && x.reversed==y.reversed &&
        x.next==y.next && x.prev==y.prev && x.twin==y.twin) &&
    a.loops.len()==b.loops.len() &&
    a.loops.iter().zip(&b.loops).all(|(x,y)|
        x.id==y.id && x.face==y.face && x.role==y.role &&
        x.first_coedge==y.first_coedge) &&
    a.faces.len()==b.faces.len() &&
    a.faces.iter().zip(&b.faces).all(|(x,y)|
        x.id==y.id && x.loops==y.loops &&
        x.geometry==y.geometry && x.sense==y.sense) &&
    a.shells.len()==b.shells.len() &&
    a.shells.iter().zip(&b.shells).all(|(x,y)|
        x.id==y.id && x.closed==y.closed && x.faces==y.faces)
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum ImprintHold {
    AlreadyImprinted,
    NotClosedTransverseSection,
    NearCapOrDegenerateSection,
    UnsupportedSource,
    AmbiguousInput,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct ImprintCommit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub new_vertex:VertexId,
    pub new_ellipse_edge:EdgeId,
    pub new_upper_face:FaceId,
    pub ellipse_coedges:[CoedgeId;2],
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum ImprintResult { Committed(ImprintCommit), Held(ImprintHold) }

/// Owns a *validated* canonical B-rep model, with full-snapshot rollback.
/// Only the source cylinder's own side face can be imprinted. No arbitrary
/// face splitting, body Boolean, multi-shell edit or topological naming is
/// authorized by this operator.
pub struct ClosedCylinderImprint {
    source:CylindricalBrep,
    current:BrepModel,
    revision:u64,
    undo:Vec<BrepModel>,
    redo:Vec<BrepModel>,
}
impl ClosedCylinderImprint {
    pub fn new(source:CylindricalBrep)->Result<Self,GeometryError>{
        source.validate()?;
        let current=BrepModel::from_cylinder(&source)?;
        Ok(Self{source,current,revision:0,undo:vec![],redo:vec![]})
    }
    pub fn model(&self)->&BrepModel {&self.current}
    pub fn source_cylinder(&self)->&CylindricalBrep {&self.source}
    pub fn revision(&self)->u64 {self.revision}
    pub fn undo_depth(&self)->usize {self.undo.len()}
    pub fn redo_depth(&self)->usize {self.redo.len()}

    /// A narrowly scoped executor for the closed, non-tangent,
    /// plane/cylinder topology-intent branch. The general delta remains
    /// can_apply=false; this operator *only* imprints the original
    /// cylinder's lateral surface, leaving both input volumes untouched.
    pub fn imprint(&mut self,graph:&IntersectionGraph,
        plan:&TopologyDeltaPlan,first:&BrepModel,second:&BrepModel,
    )->Result<ImprintResult,GeometryError>{
        graph.validate_against(first,graph.source_faces[0].face,
            second,graph.source_faces[1].face)?;
        plan.validate_against(graph,first,graph.source_faces[0].face,
            second,graph.source_faces[1].face)?;
        if self.current.faces.len()!=3 {
            return Ok(ImprintResult::Held(ImprintHold::AlreadyImprinted));
        }
        if graph.state!=GraphState::Ready {
            return Ok(ImprintResult::Held(ImprintHold::AmbiguousInput));
        }
        if plan.state!=PlanState::PreviewOnly || plan.can_apply ||
            plan.arcs.len()!=1 || graph.arcs.len()!=1 ||
            graph.nodes.len()!=1 || graph.nodes[0].kind!=GraphNodeKind::PeriodicAnchor ||
            plan.arcs[0].role!=PlanArcRole::ProposeNewCarrier ||
            !graph.arcs[0].closed || graph.arcs[0].contact!=crate::TrimContact::Crossing ||
            graph.arcs[0].from!=graph.arcs[0].to ||
            (graph.arcs[0].parameter_start).abs()>1e-9 ||
            (graph.arcs[0].parameter_end-TAU).abs()>1e-9 {
            return Ok(ImprintResult::Held(ImprintHold::NotClosedTransverseSection));
        }
        let (plane,cylinder)=match (
            first.face_surface(graph.source_faces[0].face),
            second.face_surface(graph.source_faces[1].face),
        ){
            (Some(Surface3::Plane(_)),Some(Surface3::Cylinder(c)))=>
                (first,c),
            _=>return Ok(ImprintResult::Held(ImprintHold::UnsupportedSource)),
        };
        if plane.faces.is_empty() ||
            graph.source_faces[1].face!=FaceId(2) ||
            cylinder.frame!=self.source.frame ||
            cylinder.radius!=self.source.radius {
            return Ok(ImprintResult::Held(ImprintHold::UnsupportedSource));
        }

        // Model identity is checked against the unmodified source, not a
        // plausibly similar cylinder or a forged allocation of physical IDs.
        let expected=BrepModel::from_cylinder(&self.source)?;
        if !same_cylinder_source(second,&expected) {
            return Ok(ImprintResult::Held(ImprintHold::UnsupportedSource));
        }
        let IntersectionCurve::PlaneCylinderSection(section)=graph.arcs[0].curve else {
            return Ok(ImprintResult::Held(ImprintHold::NotClosedTransverseSection));
        };
        if section.cylinder()!=cylinder {
            return Ok(ImprintResult::Held(ImprintHold::UnsupportedSource));
        }
        let [mean,cosine,sine]=section.axial_coefficients()?;
        let span=cosine.hypot(sine);
        let margin=self.source.tolerance.length_at(
            self.source.height.max(self.source.radius))*64.0;
        if !(mean-span>margin && mean+span<self.source.height-margin) {
            return Ok(ImprintResult::Held(ImprintHold::NearCapOrDegenerateSection));
        }
        let ellipse=section.harmonic_curve()?;
        let position=ellipse.evaluate(0.0)?;
        let seam_point=self.source.frame.origin
            .add(self.source.frame.x.scale(self.source.radius))
            .add(self.source.frame.normal.scale(mean+cosine));
        if distance(position,seam_point)>margin {
            return Err(invalid("analytical section does not meet the original seam"));
        }
        let checkpoint=self.revision.checked_add(1)
            .ok_or_else(||invalid("imprint revision overflow"))?;
        let candidate=build_closed_imprint(&self.source,&self.current,ellipse,
            CylinderSection2{angle_offset:0.0,direction:1.0,
                axial_center:mean,axial_cosine:cosine,axial_sine:sine},
            position,mean+cosine,mean-span,mean+span)?;
        candidate.validate()?;
        if candidate.mass!=expected.mass || candidate.bbox!=expected.bbox ||
           candidate.euler_characteristic()!=2 {
            return Err(invalid("imprint modified source mass, bounds or Euler invariant"));
        }
        let old=std::mem::replace(&mut self.current,candidate);
        if self.undo.len()==MAX_HISTORY {self.undo.remove(0);}
        self.undo.push(old);
        self.redo.clear();
        let before=self.revision;self.revision=checkpoint;
        Ok(ImprintResult::Committed(ImprintCommit{
            revision_before:before,revision_after:checkpoint,
            new_vertex:VertexId(2),new_ellipse_edge:EdgeId(4),
            new_upper_face:FaceId(3),ellipse_coedges:[CoedgeId(8),CoedgeId(9)],
        }))
    }

    pub fn undo(&mut self)->Result<bool,GeometryError>{
        if self.undo.is_empty(){return Ok(false);}
        let next=self.revision.checked_add(1).ok_or_else(||invalid("imprint revision overflow"))?;
        let restored=self.undo.last().expect("nonempty history").clone();
        restored.validate()?;
        let old=std::mem::replace(&mut self.current,restored);
        self.undo.pop();self.redo.push(old);self.revision=next;
        Ok(true)
    }
    pub fn redo(&mut self)->Result<bool,GeometryError>{
        if self.redo.is_empty(){return Ok(false);}
        let next=self.revision.checked_add(1).ok_or_else(||invalid("imprint revision overflow"))?;
        let restored=self.redo.last().expect("nonempty redo").clone();
        restored.validate()?;
        let old=std::mem::replace(&mut self.current,restored);
        self.redo.pop();if self.undo.len()==MAX_HISTORY{self.undo.remove(0);}
        self.undo.push(old);self.revision=next;
        Ok(true)
    }
}

fn build_closed_imprint(source:&CylindricalBrep,model:&BrepModel,
    ellipse:crate::HarmonicCurve3,section_uv:CylinderSection2,
    junction:crate::Point3,seam_v:f64,min_v:f64,max_v:f64,
)->Result<BrepModel,GeometryError>{
    if model.origin!=BrepOrigin::AnalyticCylinder ||
        model.vertices.len()!=2 || model.edges.len()!=3 ||
        model.coedges.len()!=6 || model.faces.len()!=3 ||
        model.loops.len()!=3 || model.shells.len()!=1 {
        return Err(invalid("closed imprint requires a pristine canonical cylinder"));
    }
    let t=source.tolerance;
    let mut out=model.clone();
    out.vertices.push(Vertex{id:VertexId(2),position:junction});
    out.edges[2].end=VertexId(2);
    out.geometry.edges[2].domain=domain(0.0,seam_v);
    out.edges.push(Edge{id:EdgeId(3),start:VertexId(2),end:VertexId(1),
        coedges:[CoedgeId(6),CoedgeId(7)]});
    out.geometry.edges.push(EdgeGeometry{
        curve:out.geometry.edges[2].curve,domain:domain(seam_v,source.height),tolerance:t,
    });
    let ellipse_id=Curve3Id(u32::try_from(out.geometry.curves3.len())
        .map_err(|_|invalid("spatial curve catalog exhausted"))?);
    out.geometry.curves3.push(Curve3::Harmonic(ellipse));
    out.edges.push(Edge{id:EdgeId(4),start:VertexId(2),end:VertexId(2),
        coedges:[CoedgeId(8),CoedgeId(9)]});
    out.geometry.edges.push(EdgeGeometry{
        curve:ellipse_id,domain:domain(0.0,TAU),tolerance:t});

    // Lower patch cycle: bottom circle, U=2PI seam, ellipse backwards,
    // U=0 seam. Upper patch: ellipse forward, U=2PI seam, top circle
    // backwards, U=0 seam. Both same-face seam uses have opposite direction.
    out.coedges[4].next=CoedgeId(9);
    out.coedges[5].prev=CoedgeId(9);
    out.coedges[2].face=FaceId(3);
    out.coedges[2].loop_id=LoopId(3);
    out.coedges[2].next=CoedgeId(7);
    out.coedges[2].prev=CoedgeId(6);
    out.coedges.extend([
        Coedge{id:CoedgeId(6),edge:EdgeId(3),face:FaceId(3),loop_id:LoopId(3),
            reversed:false,next:CoedgeId(2),prev:CoedgeId(8),twin:CoedgeId(7)},
        Coedge{id:CoedgeId(7),edge:EdgeId(3),face:FaceId(3),loop_id:LoopId(3),
            reversed:true,next:CoedgeId(8),prev:CoedgeId(2),twin:CoedgeId(6)},
        Coedge{id:CoedgeId(8),edge:EdgeId(4),face:FaceId(3),loop_id:LoopId(3),
            reversed:false,next:CoedgeId(6),prev:CoedgeId(7),twin:CoedgeId(9)},
        Coedge{id:CoedgeId(9),edge:EdgeId(4),face:FaceId(2),loop_id:LoopId(2),
            reversed:true,next:CoedgeId(5),prev:CoedgeId(4),twin:CoedgeId(8)},
    ]);
    out.loops.push(Loop{id:LoopId(3),face:FaceId(3),role:LoopRole::Outer,
        first_coedge:CoedgeId(8)});
    out.faces[2].geometry.bounds.v=domain(0.0,max_v);
    let mut upper_face=out.faces[2].clone();
    upper_face.id=FaceId(3);
    upper_face.loops=vec![LoopId(3)];
    upper_face.geometry=FaceGeometry{
        surface:upper_face.geometry.surface,outer_loop:LoopId(3),
        bounds:ParameterDomain2{v:domain(min_v,source.height),
            ..upper_face.geometry.bounds},tolerance:t,
    };
    upper_face.sense=FaceSense::Forward;
    out.faces.push(upper_face);
    out.shells[0].faces.push(FaceId(3));

    // Keep the original face's circle trim carriers unaltered and append
    // COW trims for each physically split seam and the new ellipse.
    out.geometry.coedges[2].face=FaceId(3);
    out.geometry.coedges[2].loop_id=LoopId(3);
    let mut bind=|coedge:u32,curve:Curve2,interval:ParameterRange,face:u32,lp:u32|->Result<(),GeometryError>{
        let id=Curve2Id(u32::try_from(out.geometry.curves2.len())
            .map_err(|_|invalid("UV catalog exhausted"))?);
        out.geometry.curves2.push(curve);
        let entry=CoedgeGeometry{pcurve:id,domain:interval,
            face:FaceId(face),loop_id:LoopId(lp)};
        if coedge<(out.geometry.coedges.len() as u32){
            out.geometry.coedges[coedge as usize]=entry;
        }else if coedge==out.geometry.coedges.len() as u32{
            out.geometry.coedges.push(entry);
        }else{return Err(invalid("coedge geometry order diverged"));}
        Ok(())
    };
    let unit=domain(0.0,1.0);
    bind(4,line(uv(TAU,0.0),uv(TAU,seam_v)),unit,2,2)?;
    bind(5,line(uv(0.0,seam_v),uv(0.0,0.0)),unit,2,2)?;
    bind(6,line(uv(TAU,seam_v),uv(TAU,source.height)),unit,3,3)?;
    bind(7,line(uv(0.0,source.height),uv(0.0,seam_v)),unit,3,3)?;
    bind(8,Curve2::CylinderSection(section_uv),domain(0.0,TAU),3,3)?;
    bind(9,Curve2::CylinderSection(CylinderSection2{
        angle_offset:TAU,direction:-1.0,..section_uv}),domain(0.0,TAU),2,2)?;
    Ok(out)
}
