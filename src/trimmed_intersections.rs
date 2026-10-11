//! First-party clipping of analytic supporting-surface intersections against
//! *real* bounded B-rep faces. This module never infers exactness from facets.
//!
//! Supported face trims: a simple single-loop planar polygon, a circular disk
//! on a plane, or the full periodic lateral rectangle of a right cylinder.
//! It deliberately fails closed on coincidence, holes, unrelated periodic
//! trims, unsupported surfaces, and numerically unresolved events.
use crate::{
    BrepModel, Curve2, FaceId, GeometryError, GeometryTolerance, IntersectionAmbiguity,
    IntersectionCertificate, IntersectionContact, IntersectionCurve, ParameterRange,
    Point2Param, Point3, Surface3, SurfaceIntersection, SurfaceIntersectionEngine,
};
use std::f64::consts::{PI, TAU};

fn invalid(msg: impl Into<String>) -> GeometryError {
    GeometryError::InvalidEdit(msg.into())
}
fn uncertain(msg: &str) -> GeometryError {
    invalid(format!("uncertain trim event: {msg}"))
}
fn event_result(events:Result<Vec<f64>,GeometryError>)
    -> Result<Option<Vec<f64>>,GeometryError>
{
    match events {
        Err(GeometryError::InvalidEdit(message)) if message.starts_with("uncertain trim event:")=>Ok(None),
        Err(error)=>Err(error),
        Ok(v)=>Ok(Some(v)),
    }
}
fn mag(v:Point3)->f64 {v.x.hypot(v.y).hypot(v.z)}
fn dist(a:Point3,b:Point3)->f64 {mag(a.sub(b))}
fn cross(a:[f64;2],b:[f64;2])->f64 {a[0]*b[1]-a[1]*b[0]}
fn wrap(u:f64)->f64 {u.rem_euclid(TAU)}
fn param(range:ParameterRange,f:f64)->Result<f64,GeometryError>{
    match range {
        ParameterRange::Bounded{start,end} if start.is_finite() && end.is_finite() && end>start =>
            Ok(start+(end-start)*f),
        _=>Err(invalid("trim requires a finite parameter domain")),
    }
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TrimPosition {Outside, Interior, Boundary}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TrimContact {Crossing, Tangent, BoundaryOverlap}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TrimAmbiguity {NearParallel, NearTangency, CoincidentSupports, UnresolvedBoundaryEvent}

/// An analytic support curve restricted to a finite interval; its parameter
/// belongs to the untrimmed IntersectionCurve, not to any tessellation.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TrimmedBranch {
    pub carrier:IntersectionCurve,
    pub parameter_start:f64,
    pub parameter_end:f64,
    pub start:Point3,
    pub end:Point3,
    pub contact:TrimContact,
    pub closed:bool,
}
impl TrimmedBranch {
    pub fn evaluate(&self,fraction:f64)->Result<Point3,GeometryError>{
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction){
            return Err(invalid("trimmed parameter fraction must lie in [0,1]"));
        }
        self.carrier.evaluate(self.parameter_start+(self.parameter_end-self.parameter_start)*fraction)
    }
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TrimmedPoint {pub point:Point3,pub contact:TrimContact}
#[derive(Debug,Clone,PartialEq)]
pub enum TrimmedIntersection {
    Empty,
    Ambiguous {reason:TrimAmbiguity},
    Geometry {
        branches:Vec<TrimmedBranch>,
        contacts:Vec<TrimmedPoint>,
        /// Sample-based residual on both infinite supporting surfaces, NOT
        /// a proof of interval completeness or Boolean classification.
        certificate:IntersectionCertificate,
    },
}
impl TrimmedIntersection {
    pub fn kind(&self)->&'static str{
        match self {
            Self::Empty=>"empty",Self::Ambiguous{..}=>"ambiguous",
            Self::Geometry{branches,contacts,..} if !branches.is_empty() && !contacts.is_empty()=>"mixed",
            Self::Geometry{branches,..} if !branches.is_empty()=>"segments",
            Self::Geometry{..}=>"points",
        }
    }
    pub fn branches(&self)->&[TrimmedBranch]{
        if let Self::Geometry{branches,..}=self {branches} else {&[]}
    }
    pub fn contacts(&self)->&[TrimmedPoint]{
        if let Self::Geometry{contacts,..}=self {contacts} else {&[]}
    }
}

#[derive(Debug,Clone)]
enum TrimDomain {
    Polygon {frame:crate::Frame3, vertices:Vec<Point2Param>},
    Disk {frame:crate::Frame3, center:Point2Param, radius:f64},
    PeriodicSide {cylinder:crate::CylinderSurface, lower:f64, upper:f64},
}
impl TrimDomain {
    fn from_face(body:&BrepModel,id:FaceId,tol:GeometryTolerance)->Result<Self,GeometryError>{
        let face=body.faces.get(id.0 as usize).filter(|f|f.id==id)
            .ok_or_else(||invalid(format!("unknown trim face {:?}",id)))?;
        if face.loops.len()!=1 {return Err(invalid("trimmed intersection currently requires exactly one outer loop"));}
        let surface=body.face_surface(id).ok_or_else(||invalid("face has no support surface"))?;
        match surface {
            Surface3::Plane(plane)=>{
                let ids=loop_coedges(body,face.loops[0])?;
                if ids.len()==1 {
                    let co=body.coedge_pcurve(ids[0]).ok_or_else(||invalid("missing circular disk pcurve"))?;
                    if let Curve2::Circle(circle)=co {
                        if circle.radius<=tol.absolute_length {return Err(invalid("disk radius below tolerance"));}
                        return Ok(Self::Disk{frame:plane.frame,center:circle.center,radius:circle.radius});
                    }
                }
                if ids.len()<3 {return Err(invalid("planar polygon needs at least 3 boundary coedges"));}
                let mut vertices=Vec::with_capacity(ids.len());
                for id in ids {
                    let co=body.coedge_pcurve(id).ok_or_else(||invalid("missing polygon pcurve"))?;
                    let Curve2::Line(line)=co else {
                        return Err(invalid("mixed planar circle/line loops are not supported by the first trim slice"));
                    };
                    let g=&body.geometry.coedges[id.0 as usize];
                    vertices.push(line.evaluate(param(g.domain,0.0)?)?);
                }
                Ok(Self::Polygon{frame:plane.frame,vertices})
            }
            Surface3::Cylinder(cylinder)=>{
                // Only the full native periodic side rectangle is currently
                // provable. The two u boundaries are a seam, not a physical
                // edge of the face. Validated generic B-rep UV trims must not
                // be replaced by a bounding-box guess for arbitrary geometry.
                let ids=loop_coedges(body,face.loops[0])?;
                if ids.len()!=4 || !face.geometry.bounds.u_periodic ||
                   face.geometry.bounds.v_periodic ||
                   !matches!(face.geometry.bounds.u,ParameterRange::Bounded{start,end}
                       if start.abs()<=tol.angular && (end-TAU).abs()<=tol.angular) {
                    return Err(invalid("only an uncut periodic cylindrical side is supported"));
                }
                let ParameterRange::Bounded{start, end}=face.geometry.bounds.v else {
                    return Err(invalid("cylindrical side needs finite axial bounds"));
                };
                if start>=end {return Err(invalid("invalid cylindrical side axial bounds"));}
                Ok(Self::PeriodicSide{cylinder,lower:start,upper:end})
            }
            Surface3::Sphere(_)=>Err(invalid("sphere-face clipping is not implemented")),
        }
    }
    fn classify(&self,p:Point3,tol:GeometryTolerance)->Result<TrimPosition,GeometryError>{
        match self {
            Self::Polygon{frame,vertices}=>{
                let uv=frame.project(p)?;
                let eps=tol.absolute_length;
                let mut inside=false;
                for k in 0..vertices.len() {
                    let a=vertices[k];let b=vertices[(k+1)%vertices.len()];
                    let dx=b.u-a.u;let dy=b.v-a.v;
                    let vx=uv.u-a.u;let vy=uv.v-a.v;
                    let length=dx.hypot(dy);
                    if length<=eps{return Err(invalid("degenerate planar trim edge"));}
                    let along=(vx*dx+vy*dy)/(length*length);
                    if along>=-eps/length && along<=1.0+eps/length &&
                       cross([dx,dy],[vx,vy]).abs()/length<=eps {
                        return Ok(TrimPosition::Boundary);
                    }
                    if (a.v>uv.v)!=(b.v>uv.v) {
                        let orient=crate::orient2d([a.u,a.v],[b.u,b.v],[uv.u,uv.v])?;
                        let to_right=if b.v>a.v {orient==crate::Orientation::Positive}
                            else {orient==crate::Orientation::Negative};
                        if to_right {inside=!inside;}
                    }
                }
                Ok(if inside {TrimPosition::Interior}else{TrimPosition::Outside})
            }
            Self::Disk{frame,center,radius}=>{
                let uv=frame.project(p)?;
                let d=(uv.u-center.u).hypot(uv.v-center.v);
                Ok(if (d-radius).abs()<=tol.absolute_length {TrimPosition::Boundary}
                   else if d<*radius {TrimPosition::Interior}
                   else {TrimPosition::Outside})
            }
            Self::PeriodicSide{cylinder,lower,upper}=>{
                let local=p.sub(cylinder.frame.origin);
                let axial=local.dot(cylinder.frame.normal);
                let eps=tol.length_at(cylinder.radius.max(upper-lower));
                Ok(if axial<lower-eps || axial>upper+eps {TrimPosition::Outside}
                   else if (axial-lower).abs()<=eps || (axial-upper).abs()<=eps {
                       TrimPosition::Boundary
                   } else {TrimPosition::Interior})
            }
        }
    }

    fn events(&self,curve:IntersectionCurve,cylinder:Option<crate::CylinderSurface>,tol:GeometryTolerance)->Result<Vec<f64>,GeometryError>{
        match curve {
            IntersectionCurve::Line(line)=>self.line_events(line,tol),
            IntersectionCurve::PlaneCylinderSection(section)=>self.section_events(section,
                cylinder.ok_or_else(||invalid("missing generating cylinder"))?,tol),
        }
    }
    fn line_events(&self,line:crate::Line3,tol:GeometryTolerance)->Result<Vec<f64>,GeometryError>{
        match self {
            Self::Polygon{frame,vertices}=>{
                let p=frame.project(line.origin)?;
                let direction=[line.direction.dot(frame.x),line.direction.dot(frame.y)];
                let len=direction[0].hypot(direction[1]);
                if len<0.5{return Err(invalid("line direction lies outside face support plane"));}
                let mut out=Vec::new();
                for k in 0..vertices.len() {
                    let a=vertices[k];let b=vertices[(k+1)%vertices.len()];
                    let segment=[b.u-a.u,b.v-a.v];
                    let q=[a.u-p.u,a.v-p.v];
                    let den=cross(direction,segment);
                    let size=segment[0].hypot(segment[1]);
                    if den.abs()<=tol.angular*size {
                        // Collinear boundary: clip its two endpoints normally,
                        // but don't invent a transverse crossing.
                        let side=cross(direction,q).abs()/len;
                        if side<=tol.absolute_length {
                            out.push((q[0]*direction[0]+q[1]*direction[1])/(len*len));
                            let to_b=[b.u-p.u,b.v-p.v];
                            out.push((to_b[0]*direction[0]+to_b[1]*direction[1])/(len*len));
                        }
                        continue;
                    }
                    let f=cross(q,direction)/den;
                    if f>=-tol.absolute_length/size && f<=1.0+tol.absolute_length/size {
                        out.push(cross(q,segment)/den);
                    }
                }
                Ok(out)
            }
            Self::Disk{frame,center,radius}=>{
                let p=frame.project(line.origin)?;
                let x=[p.u-center.u,p.v-center.v];
                let d=[line.direction.dot(frame.x),line.direction.dot(frame.y)];
                let a=d[0]*d[0]+d[1]*d[1];
                if a<0.25 {return Err(invalid("disk support line cannot be projected robustly"));}
                let b=x[0]*d[0]+x[1]*d[1];
                let disc=b*b-a*(x[0]*x[0]+x[1]*x[1]-radius*radius);
                let scale=a*radius*radius;
                if disc < -tol.length_at(*radius)*radius {return Ok(Vec::new());}
                if disc.abs()<=tol.length_at(*radius)*radius && disc.abs()>64.0*f64::EPSILON*scale {
                    return Err(uncertain("near-tangent disk clipping"));
                }
                let root=disc.max(0.0).sqrt();
                if root==0.0 {Ok(vec![-b/a])} else {Ok(vec![(-b-root)/a,(-b+root)/a])}
            }
            Self::PeriodicSide{cylinder,lower,upper}=>{
                let t0=line.origin.sub(cylinder.frame.origin).dot(cylinder.frame.normal);
                let dir=line.direction.dot(cylinder.frame.normal);
                if dir.abs()<=tol.angular {
                    if t0>=*lower && t0<=*upper {return Err(uncertain("generator parallel to cylindrical trim bounds"));}
                    return Ok(Vec::new());
                }
                Ok(vec![(lower-t0)/dir,(upper-t0)/dir])
            }
        }
    }
    fn section_events(&self,section:crate::PlaneCylinderSection,section_cylinder:crate::CylinderSurface,tol:GeometryTolerance)
        ->Result<Vec<f64>,GeometryError>
    {
        match self {
            Self::PeriodicSide{cylinder,lower,upper}=>{
                // For the periodic section S(u,v(u)), solve its exact axial
                // boundary equations A cos u + B sin u = target.
                let plane_height0=section.axial_parameter(0.0)?;
                let plane_height_pi=section.axial_parameter(PI)?;
                let mean=(plane_height0+plane_height_pi)*0.5;
                let cos_coeff=(plane_height0-plane_height_pi)*0.5;
                let sin_coeff=section.axial_parameter(PI*0.5)?-mean;
                let amp=cos_coeff.hypot(sin_coeff);
                let mut out=Vec::new();
                if amp<=tol.absolute_length {return Ok(out);}
                let phase=sin_coeff.atan2(cos_coeff);
                for target in [*lower,*upper] {
                    let val=(target-mean)/amp;
                    if val.abs()>1.0+tol.absolute_length/amp {continue;}
                    if (val.abs()-1.0).abs()<=tol.absolute_length/amp &&
                        (val.abs()-1.0).abs()>32.0*f64::EPSILON {
                        return Err(uncertain("axial section near cylinder axial trim"));
                    }
                    if val.abs()<=1.0 {
                        let delta=val.acos();
                        out.push(wrap(phase+delta));out.push(wrap(phase-delta));
                    }
                }
                // Only valid if this is the actual section's cylinder.
                // Residual validation on both supports is performed afterward.
                if !cylinder.radius.is_finite(){return Err(invalid("invalid side radius"));}
                Ok(out)
            }
            Self::Polygon{vertices,..}=>{
                let mut out=Vec::new();
                // Planar polygon edges intersect the generating cylinder by a
                // quadratic. These are exact algebraic event equations in
                // real arithmetic, not sampling-based guesses.
                // Edge carrier endpoints are recovered in world coordinates
                // via the surface frame below.
                let Self::Polygon{frame,..}=self else{unreachable!()};
                for k in 0..vertices.len(){
                    let a=vertices[k];let b=vertices[(k+1)%vertices.len()];
                    let p=frame.position(a.u,a.v)?;let q=frame.position(b.u,b.v)?;
                    out.extend(segment_cylinder_angles(p,q,section_cylinder,tol)?);
                }
                Ok(out)
            }
            Self::Disk{frame,center,radius}=>{
                // The full native cap is a circular disk. A horizontal
                // cylinder section is concentric with that disk; it is either
                // entirely inside, outside, or lies on its circular boundary.
                // Non-horizontal elliptic sections need a conic/conic solver.
                let low=section.axial_parameter(0.0)?;
                let high=section.axial_parameter(PI)?;
                let side=section.axial_parameter(PI/2.0)?;
                if (low-high).abs()>tol.absolute_length || (low-side).abs()>tol.absolute_length {
                    return Err(invalid("oblique conic/circular disk boundary not supported"));
                }
                let center_world=section_cylinder.frame.origin
                    .add(section_cylinder.frame.normal.scale(low));
                let coords=frame.project(center_world)?;
                let offset=(coords.u-center.u).hypot(coords.v-center.v);
                if offset>tol.absolute_length ||
                   (frame.normal.dot(section_cylinder.frame.normal).abs()-1.0).abs()>tol.angular ||
                   (radius-section_cylinder.radius).abs()>tol.length_at(*radius) &&
                   *radius<section_cylinder.radius {
                    return Err(invalid("nonconcentric circular disk/section requires circle intersection"));
                }
                Ok(Vec::new())
            }
        }
    }
}


/// Classify a point against a validated face's actual UV trimming domain.
/// Internal first-party graph validation uses this exact classifier instead of
/// inventing face membership from bounding boxes or tessellation.
pub(crate) fn classify_trim_point(
    model:&BrepModel,face:FaceId,point:Point3,tolerance:GeometryTolerance,
)->Result<TrimPosition,GeometryError>{
    let domain=TrimDomain::from_face(model,face,tolerance)?;
    domain.classify(point,tolerance)
}

fn segment_cylinder_angles(start:Point3,end:Point3,cylinder:crate::CylinderSurface,tol:GeometryTolerance)
    ->Result<Vec<f64>,GeometryError>
{
    let p=start.sub(cylinder.frame.origin);
    let d=end.sub(start);
    let radial=|v:Point3|v.sub(cylinder.frame.normal.scale(v.dot(cylinder.frame.normal)));
    let a=radial(d);let b=radial(p);
    let qa=a.dot(a);let qb=2.0*a.dot(b);let qc=b.dot(b)-cylinder.radius*cylinder.radius;
    if qa<=tol.absolute_length*tol.absolute_length {
        if qc.abs()<=tol.length_at(cylinder.radius)*cylinder.radius {
            return Err(uncertain("polygon boundary overlaps cylinder generator"));
        }
        return Ok(Vec::new());
    }
    let discriminant=qb*qb-4.0*qa*qc;
    let budget=tol.length_at(cylinder.radius)*cylinder.radius*qa;
    if discriminant < -budget{return Ok(Vec::new());}
    if discriminant.abs()<=budget && discriminant.abs()>64.0*f64::EPSILON*qa*cylinder.radius*cylinder.radius {
        return Err(uncertain("near-tangent polygon boundary/cylinder intersection"));
    }
    let root=discriminant.max(0.0).sqrt();
    let candidates=[(-qb-root)/(2.0*qa),(-qb+root)/(2.0*qa)];
    let mut result=Vec::new();
    for fraction in candidates {
        if fraction< -tol.absolute_length || fraction>1.0+tol.absolute_length {continue;}
        let q=start.add(d.scale(fraction.clamp(0.0,1.0)));
        let relative=q.sub(cylinder.frame.origin);
        result.push(wrap(relative.dot(cylinder.frame.y).atan2(relative.dot(cylinder.frame.x))));
    }
    Ok(result)
}
fn loop_coedges(body:&BrepModel,loop_id:crate::LoopId)->Result<Vec<crate::CoedgeId>,GeometryError>{
    let start=body.loops.get(loop_id.0 as usize)
        .ok_or_else(||invalid("missing outer trim loop"))?.first_coedge;
    let mut out=Vec::new();let mut current=start;
    for _ in 0..body.coedges.len(){
        out.push(current);
        current=body.next_in_loop(current).ok_or_else(||invalid("broken trim loop next"))?;
        if current==start {return Ok(out);}
    }
    Err(invalid("unbounded trim loop traversal"))
}

#[derive(Debug,Clone,Copy)]
pub struct TrimmedIntersectionEngine {pub tolerance:GeometryTolerance}
impl TrimmedIntersectionEngine {
    pub fn new(tolerance:GeometryTolerance)->Result<Self,GeometryError>{tolerance.validate()?;Ok(Self{tolerance})}
    pub fn intersect_faces(&self,a:&BrepModel,face_a:FaceId,b:&BrepModel,face_b:FaceId)
        ->Result<TrimmedIntersection,GeometryError>
    {
        a.validate()?;b.validate()?;
        let sa=a.face_surface(face_a).ok_or_else(||invalid("missing first face"))?;
        let sb=b.face_surface(face_b).ok_or_else(||invalid("missing second face"))?;
        let ta=TrimDomain::from_face(a,face_a,self.tolerance)?;
        let tb=TrimDomain::from_face(b,face_b,self.tolerance)?;
        let support=SurfaceIntersectionEngine::new(self.tolerance)?.intersect(sa,sb)?;
        let (contact,curves)=match support {
            SurfaceIntersection::Empty=>return Ok(TrimmedIntersection::Empty),
            SurfaceIntersection::Coincident=>return Ok(TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::CoincidentSupports}),
            SurfaceIntersection::Ambiguous{reason,..}=>return Ok(TrimmedIntersection::Ambiguous{reason:match reason{
                IntersectionAmbiguity::NearParallel=>TrimAmbiguity::NearParallel,
                IntersectionAmbiguity::NearTangency=>TrimAmbiguity::NearTangency,
            }}),
            SurfaceIntersection::Curves{contact,curves}=>(contact,curves),
        };
        let certificate=SurfaceIntersection::Curves{contact,curves:curves.clone()}
            .certify(sa,sb,self.tolerance)?;
        let mut branches=Vec::new();let mut contacts=Vec::new();
        for curve in curves {
            let periodic=matches!(curve.parameter_domain(),ParameterRange::Bounded{..});
            let generating_cylinder=match (sa,sb) {
                (Surface3::Cylinder(c),_)|(_,Surface3::Cylinder(c))=>Some(c),
                _=>None,
            };
            let Some(mut marks)=event_result(ta.events(curve,generating_cylinder,self.tolerance))? else {
                return Ok(TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::UnresolvedBoundaryEvent});
            };
            let Some(next)=event_result(tb.events(curve,generating_cylinder,self.tolerance))? else {
                return Ok(TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::UnresolvedBoundaryEvent});
            };
            marks.extend(next);
            if periodic {marks.push(0.0);marks.push(TAU);}
            marks.retain(|x|x.is_finite());
            marks.sort_by(f64::total_cmp);
            let scale=a.bbox.max.sub(a.bbox.min);
            let span=mag(scale).max(mag(b.bbox.max.sub(b.bbox.min)));
            let eps=if periodic {self.tolerance.angular.max(self.tolerance.length_at(span)/span.max(1.0))}
                else {self.tolerance.length_at(span)};
            marks.dedup_by(|x,y|(*x-*y).abs()<=eps);
            if marks.len()>256 {return Err(invalid("too many trim boundary events"));}
            for pair in marks.windows(2) {
                let lo=pair[0];let hi=pair[1];
                if hi-lo<=eps {continue;}
                let midpoint=curve.evaluate(0.5*(lo+hi))?;
                let pa=ta.classify(midpoint,self.tolerance)?;
                let pb=tb.classify(midpoint,self.tolerance)?;
                if pa==TrimPosition::Outside || pb==TrimPosition::Outside {continue;}
                let start=curve.evaluate(lo)?;let end=curve.evaluate(hi)?;
                let ca=ta.classify(start,self.tolerance)?;
                let cb=tb.classify(start,self.tolerance)?;
                let da=ta.classify(end,self.tolerance)?;
                let db=tb.classify(end,self.tolerance)?;
                if [ca,cb,da,db].iter().any(|p|*p==TrimPosition::Outside){
                    return Ok(TrimmedIntersection::Ambiguous{reason:TrimAmbiguity::UnresolvedBoundaryEvent});
                }
                let kind=if contact==IntersectionContact::Tangent {TrimContact::Tangent}
                    else if pa==TrimPosition::Boundary || pb==TrimPosition::Boundary {TrimContact::BoundaryOverlap}
                    else {TrimContact::Crossing};
                branches.push(TrimmedBranch{carrier:curve,parameter_start:lo,parameter_end:hi,
                    start,end,contact:kind,closed:periodic && (hi-lo-TAU).abs()<=eps});
            }
            for (index,&t) in marks.iter().enumerate() {
                if periodic && (index==0 || index+1==marks.len()) {continue;}
                let point=curve.evaluate(t)?;
                let pa=ta.classify(point,self.tolerance)?;let pb=tb.classify(point,self.tolerance)?;
                if pa==TrimPosition::Outside || pb==TrimPosition::Outside ||
                    (pa!=TrimPosition::Boundary && pb!=TrimPosition::Boundary) {continue;}
                if branches.iter().any(|branch|dist(branch.start,point)<=eps || dist(branch.end,point)<=eps){continue;}
                contacts.push(TrimmedPoint{point,contact:if contact==IntersectionContact::Tangent{
                    TrimContact::Tangent}else{TrimContact::BoundaryOverlap}});
            }
        }
        if branches.is_empty() && contacts.is_empty(){return Ok(TrimmedIntersection::Empty);}
        Ok(TrimmedIntersection::Geometry{branches,contacts,certificate})
    }
}
