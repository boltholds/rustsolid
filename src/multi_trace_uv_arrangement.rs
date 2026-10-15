//! G2: source-checked, read-only planar UV arrangement of several transverse
//! analytic LINE traces with real circular source holes.
//!
//! Limitations: strictly convex line-loop outer boundary, disjoint circular
//! single-coedge holes, complete material-spanning line traces, finite limits.
//! Tangencies, missing cut intervals, overlapping traces, and ambiguous chart
//! events fail closed. This DOES NOT commit topology or classify solid volume.
use crate::{
    BrepModel, CoedgeId, Curve2, Curve3, EdgeId, FaceId, FacePairGraph, FacePairId,
    GeometryError, GeometryTolerance, GlobalArcRef, GraphOperand, GraphState,
    IntersectionCurve, LoopId, LoopRole, MultiFaceGraphState,
    MultiFaceIntersectionGraph, ParameterRange, Point2Param, Surface3,
    TrimContact, UvCircularSpan, UvCycleKind,
};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

const MAX_LINES: usize = 8;
const MAX_HOLES: usize = 64;
const MAX_PAIRS: usize = 512;
const MAX_PIECES: usize = 2048;
const AMBIGUITY_FACTOR: f64 = 8.0;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn p(u: f64, v: f64) -> Point2Param { Point2Param { u, v } }
fn add(a: Point2Param, b: [f64; 2], t: f64) -> Point2Param {
    p(a.u + b[0]*t, a.v + b[1]*t)
}
fn sub(a: Point2Param, b: Point2Param) -> [f64; 2] { [a.u-b.u, a.v-b.v] }
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 { a[0]*b[0] + a[1]*b[1] }
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 { a[0]*b[1]-a[1]*b[0] }
fn norm(a: [f64; 2]) -> f64 { a[0].hypot(a[1]) }
fn distance(a: Point2Param, b: Point2Param) -> f64 { norm(sub(a,b)) }
fn interpolate(a: Point2Param, b: Point2Param, f: f64) -> Point2Param {
    p(a.u+(b.u-a.u)*f, a.v+(b.v-a.v)*f)
}
fn area_line(a: Point2Param, b: Point2Param) -> f64 {
    0.5*(a.u*b.v-a.v*b.u)
}
fn length3(model: &BrepModel) -> f64 {
    let v=model.bbox.max.sub(model.bbox.min);
    v.x.hypot(v.y).hypot(v.z).max(1.0)
}
fn bounded(range: ParameterRange) -> Result<(f64,f64), GeometryError> {
    match range {
        ParameterRange::Bounded {start,end} if start.is_finite() &&
            end.is_finite() && start < end => Ok((start,end)),
        _ => Err(invalid("UV arrangement requires a finite source coedge domain")),
    }
}
fn loop_ids(model: &BrepModel, loop_id: LoopId) -> Result<Vec<CoedgeId>,GeometryError> {
    let mut out=Vec::new();
    let first=model.loops.get(loop_id.0 as usize)
        .ok_or_else(||invalid("missing source UV loop"))?.first_coedge;
    let mut id=first;
    for _ in 0..=model.coedges.len() {
        let co=model.coedges.get(id.0 as usize)
            .filter(|co|co.id==id && co.loop_id==loop_id)
            .ok_or_else(||invalid("invalid source loop coedge"))?;
        out.push(id);
        id=co.next;
        if id==first{return Ok(out);}
    }
    Err(invalid("unbounded source trim loop"))
}
fn polygon_area(points: &[Point2Param]) -> f64 {
    let origin=points[0];
    points.iter().enumerate().map(|(i,&x)|
        0.5*cross(sub(x,origin),
            sub(points[(i+1)%points.len()],origin))).sum()
}
fn anchored_span_area(span:MultiUvBoundarySpan,origin:Point2Param)->f64 {
    match span.circular {
        Some(circle)=>{
            let ds=circle.angle_end.sin()-circle.angle_start.sin();
            let dc=circle.angle_end.cos()-circle.angle_start.cos();
            0.5*(circle.radius*
                ((circle.center.u-origin.u)*ds-(circle.center.v-origin.v)*dc)
                +circle.radius*circle.radius*
                (circle.angle_end-circle.angle_start))
        }
        None=>0.5*cross(sub(span.start,origin),sub(span.end,origin)),
    }
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiTraceUvHold {
    UpstreamAmbiguity,
    MissingTransverseCut,
    UnsupportedSurface,
    UnsupportedTrim,
    UnsupportedCurve,
    BudgetExceeded,
    TangencyOrVertexEvent,
    OverlappingTraces,
    IncompletePhysicalTrace,
    AmbiguousEventLocation,
    NonManifoldIncidence,
    DisconnectedCell,
    AmbiguousHoleContainment,
    AreaMismatch,
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub enum MultiUvBoundaryProvenance {
    Source {
        edge: EdgeId, coedge: CoedgeId, loop_id: LoopId,
        parameter_start: f64, parameter_end: f64,
    },
    Cut { segment: u32, source: GlobalArcRef },
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct MultiUvBoundarySpan {
    pub start: Point2Param,
    pub end: Point2Param,
    pub provenance: MultiUvBoundaryProvenance,
    pub circular: Option<UvCircularSpan>,
}
impl MultiUvBoundarySpan {
    pub fn signed_area(self) -> f64 {
        if let Some(circle)=self.circular {
            let ds=circle.angle_end.sin()-circle.angle_start.sin();
            let dc=circle.angle_end.cos()-circle.angle_start.cos();
            0.5*(circle.radius*(circle.center.u*ds-circle.center.v*dc)
                +circle.radius*circle.radius*(circle.angle_end-circle.angle_start))
        } else {
            area_line(self.start,self.end)
        }
    }
    pub fn midpoint(self) -> Point2Param {
        match self.circular {
            Some(c)=>{
                let t=(c.angle_start+c.angle_end)*0.5;
                p(c.center.u+c.radius*t.cos(),c.center.v+c.radius*t.sin())
            },
            None=>interpolate(self.start,self.end,0.5),
        }
    }
    fn tangent_start(self) -> [f64;2] {
        match self.circular {
            Some(c) if c.angle_start>c.angle_end =>
                [c.angle_start.sin(),-c.angle_start.cos()],
            Some(c)=>[-c.angle_start.sin(),c.angle_start.cos()],
            None=>sub(self.end,self.start),
        }
    }
    fn tangent_end(self) -> [f64;2] {
        match self.circular {
            Some(c) if c.angle_start>c.angle_end =>
                [c.angle_end.sin(),-c.angle_end.cos()],
            Some(c)=>[-c.angle_end.sin(),c.angle_end.cos()],
            None=>sub(self.end,self.start),
        }
    }
}

#[derive(Debug,Clone,PartialEq)]
pub struct MultiUvBoundaryCycle {
    pub kind: UvCycleKind,
    pub spans: Vec<MultiUvBoundarySpan>,
    pub signed_area: f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct MultiUvFaceCell {
    pub cycles: Vec<MultiUvBoundaryCycle>,
    pub area: f64,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct MultiUvCut {
    /// Multiple pieces may originate from one physical intersection arc.
    pub segment: u32,
    pub source: GlobalArcRef,
    pub start: Point2Param,
    pub end: Point2Param,
}
#[derive(Debug,Clone,PartialEq)]
pub struct MultiTraceUvPlan {
    pub face: FaceId,
    pub operand: GraphOperand,
    pub source_pairs: Vec<FacePairId>,
    pub cells: Vec<MultiUvFaceCell>,
    pub cuts: Vec<MultiUvCut>,
    pub original_face_area: f64,
    /// Still needs G3 whole-solid classification and G4/G5 B-rep commit.
    pub can_apply: bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum MultiTraceUvResult {
    Ready(MultiTraceUvPlan),
    Held(MultiTraceUvHold),
}

#[derive(Debug,Clone,Copy)]
struct SourceSegment {
    start: Point2Param,
    end: Point2Param,
    coedge: CoedgeId,
    edge: EdgeId,
    loop_id: LoopId,
    parameter_start: f64,
    parameter_end: f64,
}
impl SourceSegment {
    fn span(self, a: f64, b: f64) -> MultiUvBoundarySpan {
        MultiUvBoundarySpan {
            start: interpolate(self.start,self.end,a),
            end: interpolate(self.start,self.end,b),
            provenance: MultiUvBoundaryProvenance::Source {
                edge:self.edge,coedge:self.coedge,loop_id:self.loop_id,
                parameter_start:self.parameter_start+
                    (self.parameter_end-self.parameter_start)*a,
                parameter_end:self.parameter_start+
                    (self.parameter_end-self.parameter_start)*b,
            },
            circular:None,
        }
    }
}
#[derive(Debug,Clone,Copy)]
struct SourceCircle {
    center: Point2Param,
    radius: f64,
    edge: EdgeId,
    coedge: CoedgeId,
    loop_id: LoopId,
}
impl SourceCircle {
    fn point(self, angle: f64) -> Point2Param {
        p(self.center.u+self.radius*angle.cos(),
          self.center.v+self.radius*angle.sin())
    }
    fn span(self, start: f64, end: f64) -> MultiUvBoundarySpan {
        MultiUvBoundarySpan {
            start:self.point(start),end:self.point(end),
            provenance:MultiUvBoundaryProvenance::Source {
                edge:self.edge,coedge:self.coedge,loop_id:self.loop_id,
                parameter_start:start,parameter_end:end,
            },
            circular:Some(UvCircularSpan {
                center:self.center,radius:self.radius,
                angle_start:start,angle_end:end,
            }),
        }
    }
}
struct SourceFace {
    outer: Vec<SourceSegment>,
    holes: Vec<SourceCircle>,
    original_area: f64,
    extent: f64,
}
fn source_face(model: &BrepModel, face: FaceId, eps: f64)
    -> Result<Result<SourceFace,MultiTraceUvHold>,GeometryError>
{
    let selected=model.faces.get(face.0 as usize)
        .filter(|candidate|candidate.id==face)
        .ok_or_else(||invalid("unknown planar UV source face"))?;
    if !matches!(model.face_surface(face),Some(Surface3::Plane(_))) {
        return Ok(Err(MultiTraceUvHold::UnsupportedSurface));
    }
    if selected.loops.is_empty() || selected.loops.len()>MAX_HOLES+1 {
        return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
    }
    let ids=loop_ids(model,selected.loops[0])?;
    if ids.len()<3 || ids.len()>MAX_PIECES/4 {
        return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
    }
    let mut raw=Vec::new();
    for id in ids {
        let co=model.coedges.get(id.0 as usize)
            .ok_or_else(||invalid("source coedge is absent"))?;
        let Some(Curve2::Line(line))=model.coedge_pcurve(id) else {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        };
        let (start,end)=bounded(model.geometry.coedges[id.0 as usize].domain)?;
        let a=line.evaluate(start)?;
        let b=line.evaluate(end)?;
        if distance(a,b)<=eps*AMBIGUITY_FACTOR {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        }
        raw.push(SourceSegment{
            start:a,end:b,coedge:id,edge:co.edge,loop_id:co.loop_id,
            parameter_start:start,parameter_end:end,
        });
    }
    let points:Vec<_>=raw.iter().map(|x|x.start).collect();
    let signed=polygon_area(&points);
    if signed.abs()<=eps*eps {
        return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
    }
    let ccw=signed>0.0;
    let extent=length3(model);
    for i in 0..raw.len() {
        let a=points[i];
        let b=points[(i+1)%points.len()];
        let c=points[(i+2)%points.len()];
        let rotation=cross(sub(b,a),sub(c,b));
        if distance(raw[i].end,b)>eps*AMBIGUITY_FACTOR ||
            rotation.abs()<=eps*extent ||
            (rotation>0.0)!=ccw {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        }
    }
    if !ccw {
        for segment in &mut raw {
            std::mem::swap(&mut segment.start,&mut segment.end);
            std::mem::swap(&mut segment.parameter_start,&mut segment.parameter_end);
        }
    }
    let mut holes=Vec::new();
    let mut holes_area=0.0;
    for &lid in selected.loops.iter().skip(1) {
        let inner=model.loops.get(lid.0 as usize)
            .ok_or_else(||invalid("missing inner UV loop"))?;
        if inner.role!=LoopRole::Inner || inner.face!=face {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        }
        let edge_ids=loop_ids(model,lid)?;
        if edge_ids.len()!=1 {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        }
        let id=edge_ids[0];
        let Some(Curve2::Circle(circle))=model.coedge_pcurve(id) else {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        };
        if !circle.radius.is_finite() || circle.radius<=eps*AMBIGUITY_FACTOR {
            return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
        }
        holes_area+=PI*circle.radius*circle.radius;
        holes.push(SourceCircle{
            center:circle.center,radius:circle.radius,
            edge:model.coedges[id.0 as usize].edge,
            coedge:id,loop_id:lid,
        });
    }
    // Guard against touching, overlapping or unbounded source circles.
    // These are intentionally not repaired by numerical snapping.
    for (i,hole) in holes.iter().enumerate() {
        for segment in &raw {
            let margin=cross(sub(segment.end,segment.start),
                sub(hole.center,segment.start))/
                distance(segment.end,segment.start);
            if margin<=hole.radius+eps*AMBIGUITY_FACTOR {
                return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
            }
        }
        for other in holes.iter().skip(i+1) {
            if distance(hole.center,other.center)<=
                hole.radius+other.radius+eps*AMBIGUITY_FACTOR {
                return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
            }
        }
    }
    let area=signed.abs()-holes_area;
    if !area.is_finite() || area<=eps*eps {
        return Ok(Err(MultiTraceUvHold::UnsupportedTrim));
    }
    Ok(Ok(SourceFace{outer:raw,holes,original_area:area,extent}))
}

#[derive(Debug,Clone,Copy)]
struct SourceInterval {
    start:f64,end:f64,
    source:GlobalArcRef,
}
#[derive(Debug,Clone)]
struct TraceLine {
    origin:Point2Param,
    direction:[f64;2],
    intervals:Vec<SourceInterval>,
}
impl TraceLine {
    fn along(&self,q:Point2Param)->f64 {
        dot(sub(q,self.origin),self.direction)
    }
    fn signed_distance(&self,q:Point2Param)->f64 {
        cross(self.direction,sub(q,self.origin))
    }
    fn point(&self,t:f64)->Point2Param {
        add(self.origin,self.direction,t)
    }
}

fn line_segment_event(line:&TraceLine,segment:SourceSegment,eps:f64)
    ->Result<Option<(f64,f64)>,MultiTraceUvHold>
{
    let delta=sub(segment.end,segment.start);
    let denom=cross(line.direction,delta);
    let signed_a=line.signed_distance(segment.start);
    let signed_b=line.signed_distance(segment.end);
    if denom.abs()<=1e-12*norm(delta) {
        if signed_a.abs()<=eps*AMBIGUITY_FACTOR ||
            signed_b.abs()<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::OverlappingTraces);
        }
        return Ok(None);
    }
    let s= -signed_a/denom;
    if s< -eps/norm(delta) || s>1.0+eps/norm(delta) {
        return Ok(None);
    }
    if s<=eps*AMBIGUITY_FACTOR/norm(delta) ||
        s>=1.0-eps*AMBIGUITY_FACTOR/norm(delta) {
        return Err(MultiTraceUvHold::TangencyOrVertexEvent);
    }
    let t=line.along(interpolate(segment.start,segment.end,s));
    Ok(Some((t,s)))
}

fn circle_line_events(line:&TraceLine,hole:SourceCircle,eps:f64)
    ->Result<Vec<(f64,f64)>,MultiTraceUvHold>
{
    let center_vector=sub(hole.center,line.origin);
    let axial=dot(center_vector,line.direction);
    let perpendicular=cross(line.direction,center_vector);
    let residual=hole.radius-perpendicular.abs();
    if residual.abs()<=eps*AMBIGUITY_FACTOR {
        return Err(MultiTraceUvHold::TangencyOrVertexEvent);
    }
    if residual<0.0 {return Ok(Vec::new());}
    let half=(hole.radius*hole.radius-perpendicular*perpendicular).sqrt();
    let mut out=Vec::with_capacity(2);
    for t in [axial-half,axial+half] {
        let q=line.point(t);
        let angle=(q.v-hole.center.v).atan2(q.u-hole.center.u).rem_euclid(TAU);
        out.push((t,angle));
    }
    Ok(out)
}

fn line_crossing(a:&TraceLine,b:&TraceLine,eps:f64)
    ->Result<Option<(f64,f64)>,MultiTraceUvHold>
{
    let determinant=cross(a.direction,b.direction);
    let between=sub(b.origin,a.origin);
    if determinant.abs()<=1e-10 {
        if a.signed_distance(b.origin).abs()<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::OverlappingTraces);
        }
        return Ok(None);
    }
    let ta=cross(between,b.direction)/determinant;
    let tb=cross(between,a.direction)/determinant;
    if !ta.is_finite()||!tb.is_finite() {
        return Err(MultiTraceUvHold::AmbiguousEventLocation);
    }
    Ok(Some((ta,tb)))
}
fn sort_knots(knots:&mut Vec<f64>,eps:f64)
    ->Result<(),MultiTraceUvHold>
{
    if knots.iter().any(|x|!x.is_finite()) || knots.len()>MAX_PIECES {
        return Err(MultiTraceUvHold::BudgetExceeded);
    }
    knots.sort_by(f64::total_cmp);
    let mut unique=Vec::with_capacity(knots.len());
    for &t in knots.iter() {
        if let Some(&last)=unique.last() {
            let gap=t-last;
            // A graph endpoint, a line/circle root and a trim witness may
            // describe the SAME event; keep one parameter anchor.
            if gap<=eps {continue;}
            // Distinct roots closer than the wider ambiguity budget cannot
            // be treated as reliably separate vertices.
            if gap<=eps*AMBIGUITY_FACTOR {
                return Err(MultiTraceUvHold::AmbiguousEventLocation);
            }
        }
        unique.push(t);
    }
    *knots=unique;
    Ok(())
}
fn inside_face(face:&SourceFace,q:Point2Param,eps:f64)
    ->Result<bool,MultiTraceUvHold>
{
    for segment in &face.outer {
        let side=cross(sub(segment.end,segment.start),sub(q,segment.start)) /
            distance(segment.end,segment.start);
        if side.abs()<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::AmbiguousEventLocation);
        }
        if side<0.0 {return Ok(false);}
    }
    for hole in &face.holes {
        let radial=distance(q,hole.center)-hole.radius;
        if radial.abs()<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::AmbiguousEventLocation);
        }
        if radial<0.0 {return Ok(false);}
    }
    Ok(true)
}
fn add_span(spans:&mut Vec<MultiUvBoundarySpan>,span:MultiUvBoundarySpan)
    ->Result<(),MultiTraceUvHold>
{
    if spans.len()>=MAX_PIECES {
        return Err(MultiTraceUvHold::BudgetExceeded);
    }
    spans.push(span);
    Ok(())
}

#[derive(Debug,Clone,Copy)]
struct ArrangedEdge {
    from:usize,
    to:usize,
    span:MultiUvBoundarySpan,
}
fn locate(nodes:&mut Vec<Point2Param>,q:Point2Param,eps:f64)
    ->Result<usize,MultiTraceUvHold>
{
    let mut found=None;
    for (i,&candidate) in nodes.iter().enumerate() {
        let d=distance(candidate,q);
        if d<=eps {
            if found.is_some() {
                return Err(MultiTraceUvHold::AmbiguousEventLocation);
            }
            found=Some(i);
        } else if d<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::AmbiguousEventLocation);
        }
    }
    if let Some(i)=found {return Ok(i);}
    nodes.push(q);
    Ok(nodes.len()-1)
}
fn angle(v:[f64;2])->Result<f64,MultiTraceUvHold>{
    if !norm(v).is_finite()||norm(v)<1e-14 {
        return Err(MultiTraceUvHold::AmbiguousEventLocation);
    }
    Ok(v[1].atan2(v[0]))
}
fn extract_cycles(spans:Vec<MultiUvBoundarySpan>,eps:f64)
    ->Result<Vec<MultiUvBoundaryCycle>,MultiTraceUvHold>
{
    let mut nodes=Vec::new();
    let mut edges=Vec::new();
    for span in spans {
        if span.circular.is_none() && distance(span.start,span.end)<=eps*AMBIGUITY_FACTOR {
            return Err(MultiTraceUvHold::AmbiguousEventLocation);
        }
        let from=locate(&mut nodes,span.start,eps)?;
        let to=locate(&mut nodes,span.end,eps)?;
        if span.circular.is_none() && from==to {
            return Err(MultiTraceUvHold::NonManifoldIncidence);
        }
        edges.push(ArrangedEdge{from,to,span});
    }
    let mut incoming=vec![Vec::new();nodes.len()];
    let mut outgoing=vec![Vec::new();nodes.len()];
    for (i,e) in edges.iter().enumerate() {
        incoming[e.to].push(i);
        outgoing[e.from].push(i);
    }
    if incoming.iter().zip(outgoing.iter())
        .any(|(a,b)|a.is_empty()||a.len()!=b.len()||b.len()>4) {
        // Two transverse lines can have four outgoing halfedges.
        // Higher-order intersections need a separate certified event solver.
        return Err(MultiTraceUvHold::NonManifoldIncidence);
    }
    for uses in &outgoing {
        for (i,&left) in uses.iter().enumerate() {
            let a=angle(edges[left].span.tangent_start())?;
            for &right in uses.iter().skip(i+1) {
                let b=angle(edges[right].span.tangent_start())?;
                let delta=(a-b).rem_euclid(TAU);
                if delta<=1e-10 || (TAU-delta)<=1e-10 {
                    return Err(MultiTraceUvHold::NonManifoldIncidence);
                }
            }
        }
    }
    let mut next=vec![0;edges.len()];
    for (i,edge) in edges.iter().enumerate() {
        let tangent=edge.span.tangent_end();
        let back=angle([-tangent[0],-tangent[1]])?;
        let mut best=None::<(usize,f64)>;
        for &index in &outgoing[edge.to] {
            let forward=angle(edges[index].span.tangent_start())?;
            let clockwise=(back-forward).rem_euclid(TAU);
            // The exact opposite halfedge is the same physical cut seen
            // backwards. Never choose it as the continuation of the face to
            // the left; an endpoint with no other exit is a dangling cut.
            if clockwise<=1e-10 {continue;}
            if let Some((_,score))=best {
                if (score-clockwise).abs()<=1e-10 {
                    return Err(MultiTraceUvHold::NonManifoldIncidence);
                }
                if clockwise<score {best=Some((index,clockwise));}
            } else {
                best=Some((index,clockwise));
            }
        }
        next[i]=best.ok_or(MultiTraceUvHold::NonManifoldIncidence)?.0;
    }
    let mut visited=vec![false;edges.len()];
    let mut cycles=Vec::new();
    for first in 0..edges.len() {
        if visited[first]{continue;}
        let mut cursor=first;
        let mut path=Vec::new();
        loop {
            if visited[cursor] {
                if cursor==first {break;}
                return Err(MultiTraceUvHold::NonManifoldIncidence);
            }
            visited[cursor]=true;
            path.push(edges[cursor].span);
            cursor=next[cursor];
            if cursor==first{break;}
            if path.len()>edges.len(){
                return Err(MultiTraceUvHold::NonManifoldIncidence);
            }
        }
        // Integrate in a cycle-local coordinate frame: a cap translated
        // millions of units must not lose its small true area through
        // cancellation of world-coordinate shoelace terms.
        let anchor=path[0].start;
        let signed_area=path.iter().map(|span|
            anchored_span_area(*span,anchor)).sum::<f64>();
        if !signed_area.is_finite()||signed_area.abs()<=eps*eps {
            return Err(MultiTraceUvHold::NonManifoldIncidence);
        }
        cycles.push(MultiUvBoundaryCycle {
            kind:if signed_area>0.0 {UvCycleKind::Outer}
                else{UvCycleKind::Inner},
            spans:path,signed_area,
        });
    }
    Ok(cycles)
}
fn angle_on_arc(angle:f64,start:f64,end:f64,slack:f64)->Option<bool>{
    let low=start.min(end);
    let high=start.max(end);
    let mut matched=false;
    for k in -2..=2 {
        let value=angle+(k as f64)*TAU;
        if value>=low-slack && value<=high+slack {
            if value<=low+slack || value>=high-slack {
                return None;
            }
            matched=true;
        }
    }
    Some(matched)
}
fn ray_parity(cycle:&MultiUvBoundaryCycle,q:Point2Param,eps:f64,ray:[f64;2])
    ->Option<bool>
{
    let mut count=0usize;
    for span in &cycle.spans {
        if let Some(circle)=span.circular {
            let center=sub(circle.center,q);
            let projected=dot(center,ray);
            let perpendicular=cross(ray,center);
            let discr=circle.radius*circle.radius-perpendicular*perpendicular;
            if discr < -eps*circle.radius {continue;}
            if discr.abs()<=eps*circle.radius {return None;}
            let root=discr.sqrt();
            for t in [projected-root,projected+root] {
                if t<=eps*AMBIGUITY_FACTOR {if t>=-eps*AMBIGUITY_FACTOR{return None;}continue;}
                let hit=add(q,ray,t);
                let theta=(hit.v-circle.center.v).atan2(hit.u-circle.center.u);
                let inside=angle_on_arc(theta,circle.angle_start,circle.angle_end,
                    eps*AMBIGUITY_FACTOR/circle.radius)?;
                if inside{count+=1;}
            }
        } else {
            let direction=sub(span.end,span.start);
            let denominator=cross(ray,direction);
            let offset=sub(span.start,q);
            if denominator.abs()<=1e-12*norm(direction) {
                if cross(ray,offset).abs()<=eps*AMBIGUITY_FACTOR {return None;}
                continue;
            }
            let t=cross(offset,direction)/denominator;
            let f=cross(offset,ray)/denominator;
            if t<=eps*AMBIGUITY_FACTOR {
                if t>=-eps*AMBIGUITY_FACTOR && f>=0.0 && f<=1.0{return None;}
                continue;
            }
            if f>=-eps/norm(direction) && f<=1.0+eps/norm(direction) {
                if f<=eps*AMBIGUITY_FACTOR/norm(direction) ||
                    f>=1.0-eps*AMBIGUITY_FACTOR/norm(direction) {
                    return None;
                }
                count+=1;
            }
        }
    }
    Some(count%2==1)
}
pub(crate) fn contains(cycle:&MultiUvBoundaryCycle,q:Point2Param,eps:f64)
    ->Result<bool,MultiTraceUvHold>
{
    for direction in [0.173_f64,0.523,1.117,1.937,2.773] {
        let ray=[direction.cos(),direction.sin()];
        if let Some(inside)=ray_parity(cycle,q,eps,ray) {
            return Ok(inside);
        }
    }
    Err(MultiTraceUvHold::AmbiguousHoleContainment)
}

fn arrange_planar(
    model:&BrepModel, face_id:FaceId, operand:GraphOperand,
    pair_ids:Vec<FacePairId>, mut lines:Vec<TraceLine>, eps:f64,
)->Result<MultiTraceUvResult,GeometryError>
{
    if lines.is_empty() {return Ok(MultiTraceUvResult::Held(
        MultiTraceUvHold::MissingTransverseCut));}
    if lines.len()>MAX_LINES {return Ok(MultiTraceUvResult::Held(
        MultiTraceUvHold::BudgetExceeded));}
    let face=match source_face(model,face_id,eps)? {
        Ok(data)=>data,
        Err(reason)=>return Ok(MultiTraceUvResult::Held(reason)),
    };

    // Each line must cross the convex outer trim at two nonvertex points.
    // Line/hole and line/line intersections, as well as arc endpoints,
    // become exact-in-form event knots before any directed edge is created.
    let mut outer_events=vec![Vec::<(usize,f64,f64)>::new();lines.len()];
    let mut hole_events=vec![Vec::<Vec<(f64,f64)>>::new();lines.len()];
    let mut line_events=vec![Vec::<f64>::new();lines.len()];
    for (li,line) in lines.iter().enumerate() {
        for (edge_index,&edge) in face.outer.iter().enumerate() {
            match line_segment_event(line,edge,eps) {
                Ok(Some((t,f)))=>outer_events[li].push((edge_index,t,f)),
                Ok(None)=>{},
                Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
            }
        }
        if outer_events[li].len()!=2 {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::IncompletePhysicalTrace));
        }
        let mut hole_roots=Vec::new();
        for &hole in &face.holes {
            match circle_line_events(line,hole,eps) {
                Ok(v)=>hole_roots.push(v),
                Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
            }
        }
        hole_events[li]=hole_roots;
    }
    for i in 0..lines.len() {
        for j in i+1..lines.len() {
            match line_crossing(&lines[i],&lines[j],eps) {
                Ok(Some((ti,tj)))=>{
                    line_events[i].push(ti);
                    line_events[j].push(tj);
                },
                Ok(None)=>{},
                Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
            }
        }
    }

    let mut spans=Vec::new();
    // A source outer segment is cut at every intersection; source parameter
    // values are preserved without creating physical B-rep vertices.
    for (ei,&edge) in face.outer.iter().enumerate() {
        let mut fractions=vec![0.0,1.0];
        for events in &outer_events {
            for &(index,_,fraction) in events {
                if index==ei {fractions.push(fraction);}
            }
        }
        if let Err(why)=sort_knots(&mut fractions,eps/distance(edge.start,edge.end)) {
            return Ok(MultiTraceUvResult::Held(why));
        }
        for pair in fractions.windows(2) {
            if let Err(why)=add_span(&mut spans,edge.span(pair[0],pair[1])) {
                return Ok(MultiTraceUvResult::Held(why));
            }
        }
    }

    // Preserve real circle arcs: no chordization, even where multiple cuts
    // partition one hole boundary into many spans.
    for (hole_index,&hole) in face.holes.iter().enumerate() {
        let mut angles=Vec::new();
        for per_line in &hole_events {
            for &(_,angle) in &per_line[hole_index] {
                angles.push(angle);
            }
        }
        angles.sort_by(f64::total_cmp);
        for pair in angles.windows(2) {
            if pair[1]-pair[0]<=eps*AMBIGUITY_FACTOR/hole.radius {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::TangencyOrVertexEvent));
            }
        }
        if angles.len()>1 && angles[0]+TAU-angles[angles.len()-1]
            <=eps*AMBIGUITY_FACTOR/hole.radius {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::TangencyOrVertexEvent));
        }
        if angles.is_empty() {
            if let Err(why)=add_span(&mut spans,hole.span(0.0,-TAU)) {
                return Ok(MultiTraceUvResult::Held(why));
            }
        } else {
            for i in 0..angles.len() {
                let start=angles[i];
                let end=if i==0 {angles[angles.len()-1]-TAU} else{angles[i-1]};
                if let Err(why)=add_span(&mut spans,hole.span(start,end)) {
                    return Ok(MultiTraceUvResult::Held(why));
                }
            }
        }
    }

    let mut cuts=Vec::new();
    for (li,line) in lines.iter_mut().enumerate() {
        let mut boundary:Vec<f64>=outer_events[li].iter().map(|(_,t,_)|*t).collect();
        boundary.sort_by(f64::total_cmp);
        let lo=boundary[0];let hi=boundary[1];
        if hi-lo<=eps*AMBIGUITY_FACTOR {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::IncompletePhysicalTrace));
        }
        let mut knots=vec![lo,hi];
        for per_hole in &hole_events[li]{
            for &(t,_) in per_hole {
                if t>lo && t<hi {knots.push(t);}
            }
        }
        for &event in &line_events[li] {
            if event>lo && event<hi {knots.push(event);}
            else if (event-lo).abs()<=eps*AMBIGUITY_FACTOR ||
                    (event-hi).abs()<=eps*AMBIGUITY_FACTOR {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::TangencyOrVertexEvent));
            }
        }
        for interval in &line.intervals {
            for t in [interval.start,interval.end] {
                if t>lo && t<hi {knots.push(t);}
                else if (t-lo).abs()<=eps*AMBIGUITY_FACTOR ||
                        (t-hi).abs()<=eps*AMBIGUITY_FACTOR {
                    // Boundary alignment is normal, do not duplicate the knot.
                } else if t<lo-eps*AMBIGUITY_FACTOR ||
                    t>hi+eps*AMBIGUITY_FACTOR {
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::IncompletePhysicalTrace));
                }
            }
        }
        // Interior line-line roots outside the planar source are not added.
        // A root strictly between outer bounds must lie outside all holes
        // unless the cut really crosses that hole (already event-split).
        if let Err(why)=sort_knots(&mut knots,eps) {
            return Ok(MultiTraceUvResult::Held(why));
        }
        let mut covered=BTreeSet::new();
        for interval in knots.windows(2) {
            let midpoint=line.point((interval[0]+interval[1])*0.5);
            let material=match inside_face(&face,midpoint,eps) {
                Ok(result)=>result,
                Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
            };
            if !material {continue;}
            let candidate=line.intervals.iter()
                .enumerate()
                .filter(|(_,arc)| {
                    arc.start<=interval[0]+eps && arc.end>=interval[1]-eps
                }).collect::<Vec<_>>();
            if candidate.len()!=1 {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::IncompletePhysicalTrace));
            }
            let (source_index,source)=candidate[0];
            covered.insert(source_index);
            let start=line.point(interval[0]);
            let end=line.point(interval[1]);
            if cuts.len()>=MAX_PIECES/2 {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::BudgetExceeded));
            }
            let id=cuts.len() as u32;
            cuts.push(MultiUvCut {
                segment:id,source:source.source,start,end,
            });
            for (a,b) in [(start,end),(end,start)] {
                if let Err(why)=add_span(&mut spans,MultiUvBoundarySpan {
                    start:a,end:b,circular:None,
                    provenance:MultiUvBoundaryProvenance::Cut {
                        segment:id,source:source.source,
                    },
                }) {
                    return Ok(MultiTraceUvResult::Held(why));
                }
            }
        }
        if covered.len()!=line.intervals.len() {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::IncompletePhysicalTrace));
        }
    }

    let cycles=match extract_cycles(spans,eps) {
        Ok(cycles)=>cycles,
        Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
    };
    let mut cells=Vec::<MultiUvFaceCell>::new();
    let mut negative=Vec::new();
    for cycle in cycles {
        match cycle.kind {
            UvCycleKind::Outer=>cells.push(MultiUvFaceCell{
                area:cycle.signed_area,cycles:vec![cycle],
            }),
            UvCycleKind::Inner=>negative.push(cycle),
        }
    }
    if cells.is_empty() {
        return Ok(MultiTraceUvResult::Held(
            MultiTraceUvHold::DisconnectedCell));
    }
    for cycle in negative {
        // Negative-area boundary cycles in this bounded slice must be
        // standalone, uncut original circular holes. A cut-derived negative
        // cycle would require full topological nesting classification.
        let mut associated=None;
        for span in &cycle.spans {
            let (Some(circle),MultiUvBoundaryProvenance::Source{loop_id,..})=
                (span.circular,span.provenance) else {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::AmbiguousHoleContainment));
            };
            if associated.is_some_and(|(id,center):(LoopId,Point2Param)|
                id!=loop_id || distance(center,circle.center)>eps) {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::AmbiguousHoleContainment));
            }
            associated=Some((loop_id,circle.center));
        }
        let Some((_,center))=associated else {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::AmbiguousHoleContainment));
        };
        let mut owner=None;
        for (i,cell) in cells.iter().enumerate() {
            match contains(&cell.cycles[0],center,eps) {
                Ok(true) if owner.is_none()=>owner=Some(i),
                Ok(true)=>return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::AmbiguousHoleContainment)),
                Ok(false)=>{},
                Err(why)=>return Ok(MultiTraceUvResult::Held(why)),
            }
        }
        let Some(i)=owner else {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::AmbiguousHoleContainment));
        };
        cells[i].area+=cycle.signed_area;
        cells[i].cycles.push(cycle);
    }

    if cells.iter().any(|cell|!cell.area.is_finite()||cell.area<=eps*eps) {
        return Ok(MultiTraceUvResult::Held(MultiTraceUvHold::AreaMismatch));
    }
    let total=cells.iter().map(|cell|cell.area).sum::<f64>();
    if !total.is_finite() ||
        (total-face.original_area).abs()>eps*face.extent*512.0 {
        return Ok(MultiTraceUvResult::Held(MultiTraceUvHold::AreaMismatch));
    }
    // Every analytic cut piece is used twice, in opposite directions, by
    // the reconstructed face cycles. This checks incidence but never
    // allocates an EdgeId or commits a B-rep face.
    for cut in &cuts {
        let mut forward=0;
        let mut reverse=0;
        for span in cells.iter().flat_map(|cell|&cell.cycles)
            .flat_map(|cycle|&cycle.spans) {
            if span.provenance != (MultiUvBoundaryProvenance::Cut {
                segment:cut.segment,source:cut.source
            }) {continue;}
            if distance(span.start,cut.start)<=eps &&
                distance(span.end,cut.end)<=eps {forward+=1;}
            else if distance(span.start,cut.end)<=eps &&
                distance(span.end,cut.start)<=eps {reverse+=1;}
            else {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::NonManifoldIncidence));
            }
        }
        if forward!=1||reverse!=1 {
            return Ok(MultiTraceUvResult::Held(
                MultiTraceUvHold::NonManifoldIncidence));
        }
    }
    cells.sort_by(|a,b|{
        let anchor=|cell:&MultiUvFaceCell|cell.cycles[0].spans[0].start;
        let a=anchor(a);let b=anchor(b);
        a.u.total_cmp(&b.u).then(a.v.total_cmp(&b.v))
    });
    Ok(MultiTraceUvResult::Ready(MultiTraceUvPlan {
        face:face_id,operand,source_pairs:pair_ids,cells,cuts,
        original_face_area:face.original_area,can_apply:false,
    }))
}

/// A face-pair intersection which lives wholly on an *existing physical
/// boundary coedge* is not a new interior face split. G2 may skip it only
/// after original node witnesses agree on ONE CoedgeId and the entire
/// 3D segment is inside that coedge's bounded analytic line carrier.
/// A tangent circle, an interior coincident segment, or a synthetic chart
/// seam is never silently dropped by this gate.
fn certified_original_boundary_overlap(
    model:&BrepModel,
    source_face:FaceId,
    source_operand:GraphOperand,
    graph:&crate::IntersectionGraph,
    arc:&crate::GraphArc,
    eps:f64,
)->Result<bool,GeometryError>{
    if arc.closed || arc.contact!=TrimContact::BoundaryOverlap {
        return Ok(false);
    }
    let Some(start)=graph.nodes.get(arc.from.0 as usize)
        .filter(|node|node.id==arc.from) else{return Ok(false)};
    let Some(end)=graph.nodes.get(arc.to.0 as usize)
        .filter(|node|node.id==arc.to) else{return Ok(false)};
    for source in start.boundary_sources.iter().filter(|w|
        w.operand==source_operand && w.face==source_face) {
        if !end.boundary_sources.contains(source) {continue;}
        let Some(coedge)=model.coedges.get(source.coedge.0 as usize)
            .filter(|coedge|coedge.id==source.coedge &&
                coedge.edge==source.edge &&
                coedge.face==source_face) else{continue};
        let _=coedge;
        let Some(Curve3::Line(line))=model.edge_curve(source.edge) else{continue};
        let Some(edge)=model.geometry.edges.get(source.edge.0 as usize) else{continue};
        let ParameterRange::Bounded{start:lo,end:hi}=edge.domain else{continue};
        let speed_squared=line.direction.length_squared();
        if !speed_squared.is_finite()||speed_squared<=0.0 {continue;}
        let slack=eps/speed_squared.sqrt();
        let mut valid=true;
        for fraction in [0.0,0.125,0.25,0.5,0.75,0.875,1.0] {
            let point=arc.evaluate(fraction)?;
            let parameter=point.sub(line.origin).dot(line.direction)/speed_squared;
            if !parameter.is_finite()||parameter<lo-slack||
                parameter>hi+slack ||
                point.sub(line.evaluate(parameter)?).length_squared().sqrt()>eps
            {
                valid=false;
                break;
            }
        }
        if valid{return Ok(true);}
    }
    Ok(false)
}

impl MultiTraceUvPlan {
    pub fn validate_against(
        &self, first:&BrepModel, second:&BrepModel,
        pairs:&[FacePairGraph], tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        let rebuilt=Self::build_selected(
            first,second,pairs,self.operand,self.face,tolerance)?;
        if rebuilt!=MultiTraceUvResult::Ready(self.clone()) {
            return Err(invalid("G2 arrangement differs from source-bound reconstruction"));
        }
        Ok(())
    }

    pub fn build(
        first:&BrepModel,second:&BrepModel,
        graph:&MultiFaceIntersectionGraph,operand:GraphOperand,
        face:FaceId,tolerance:GeometryTolerance,
    )->Result<MultiTraceUvResult,GeometryError>{
        graph.validate_against(first,second)?;
        if graph.state==MultiFaceGraphState::Held {
            return Ok(MultiTraceUvResult::Held(MultiTraceUvHold::UpstreamAmbiguity));
        }
        Self::build_selected(first,second,&graph.pairs,operand,face,tolerance)
    }

    /// Explicitly selected source graph pairs are allowed for a bounded
    /// local face preview. This does not establish global Boolean completeness.
    pub fn build_selected(
        first:&BrepModel,second:&BrepModel,pairs:&[FacePairGraph],
        operand:GraphOperand,face_id:FaceId,tolerance:GeometryTolerance,
    )->Result<MultiTraceUvResult,GeometryError>{
        tolerance.validate()?;
        first.validate()?;second.validate()?;
        if pairs.len()>MAX_PAIRS {
            return Ok(MultiTraceUvResult::Held(MultiTraceUvHold::BudgetExceeded));
        }
        let model=if operand==GraphOperand::First{first}else{second};
        let Some(Surface3::Plane(plane))=model.face_surface(face_id) else{
            return Ok(MultiTraceUvResult::Held(MultiTraceUvHold::UnsupportedSurface));
        };
        let eps=tolerance.length_at(length3(model));
        let mut unique=BTreeSet::new();
        let mut pair_ids=Vec::new();
        let mut lines=Vec::<TraceLine>::new();
        for pair in pairs {
            if !unique.insert(pair.id) {
                return Err(invalid("duplicate pair IDs in selected G2 graph"));
            }
            pair.graph.validate_against(
                first,pair.first_face,second,pair.second_face)?;
            // The local graph validator does not prove that an Empty pair
            // contains every possible trimmed intersection. Recompute
            // each selected pair, including empty ones, with the declared
            // construction policy before accepting cut completeness.
            let recomputed=crate::IntersectionGraphBuilder::new(tolerance)?
                .intersect_faces(first,pair.first_face,second,pair.second_face)?;
            if pair.graph!=recomputed {
                return Err(invalid("G2 selected pair differs from source geometry"));
            }
            let selected=if operand==GraphOperand::First {
                pair.first_face==face_id
            }else{
                pair.second_face==face_id
            };
            if !selected {continue;}
            if matches!(pair.graph.state,GraphState::Ambiguous(_)) {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::UpstreamAmbiguity));
            }
            if pair.graph.state==GraphState::Empty {continue;}
            if pair.graph.arcs.is_empty() {
                return Ok(MultiTraceUvResult::Held(
                    MultiTraceUvHold::MissingTransverseCut));
            }
            let mut has_real_cut=false;
            for arc in &pair.graph.arcs {
                if certified_original_boundary_overlap(
                    model,face_id,operand,&pair.graph,arc,eps)?
                {
                    // The original physical boundary is preserved by G4's
                    // source-edge pairing, never treated as a new UV cut.
                    continue;
                }
                if arc.closed||arc.contact!=TrimContact::Crossing {
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::TangencyOrVertexEvent));
                }
                let IntersectionCurve::Line(line)=arc.curve else {
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::UnsupportedCurve));
                };
                let a=plane.frame.project(arc.evaluate(0.0)?)?;
                let b=plane.frame.project(arc.evaluate(1.0)?)?;
                let mut direction=[
                    line.direction.dot(plane.frame.x),
                    line.direction.dot(plane.frame.y),
                ];
                let magnitude=norm(direction);
                if !magnitude.is_finite()||magnitude<0.5 {
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::UnsupportedCurve));
                }
                direction=[direction[0]/magnitude,direction[1]/magnitude];
                if direction[0]< -1e-12 ||
                    (direction[0].abs()<=1e-12&&direction[1]<0.0) {
                    direction=[-direction[0],-direction[1]];
                }
                let index=lines.iter().position(|existing|{
                    cross(existing.direction,direction).abs()<1e-10 &&
                    existing.signed_distance(a).abs()<=eps
                });
                if index.is_none() && lines.iter().any(|existing|{
                    cross(existing.direction,direction).abs()<1e-10 &&
                    existing.signed_distance(a).abs()<=eps*AMBIGUITY_FACTOR
                }){
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::OverlappingTraces));
                }
                // Adjacent cuts from different opposite faces cannot be
                // stitched by projected 2D closeness alone. They require
                // G1 source-topology reconciliation at their endpoints.
                if let Some(i)=index {
                    if lines[i].intervals.iter().any(|segment|
                        segment.source.pair!=pair.id) {
                        return Ok(MultiTraceUvResult::Held(
                            MultiTraceUvHold::OverlappingTraces));
                    }
                }
                let line_index=match index {
                    Some(i)=>i,
                    None=>{
                        if lines.len()>=MAX_LINES {
                            return Ok(MultiTraceUvResult::Held(
                                MultiTraceUvHold::BudgetExceeded));
                        }
                        lines.push(TraceLine {
                            origin:a,direction,intervals:Vec::new(),
                        });
                        lines.len()-1
                    }
                };
                let trace=&mut lines[line_index];
                let ta=trace.along(a);
                let tb=trace.along(b);
                if !ta.is_finite()||!tb.is_finite()||
                    (ta-tb).abs()<=eps*AMBIGUITY_FACTOR {
                    return Ok(MultiTraceUvResult::Held(
                        MultiTraceUvHold::AmbiguousEventLocation));
                }
                trace.intervals.push(SourceInterval {
                    start:ta.min(tb),end:ta.max(tb),
                    source:GlobalArcRef {pair:pair.id,local:arc.id},
                });
                has_real_cut=true;
            }
            if has_real_cut {pair_ids.push(pair.id);}
        }
        pair_ids.sort();
        arrange_planar(model,face_id,operand,pair_ids,lines,eps)
    }
}

#[cfg(test)]
mod arrangement_tests {
    use super::*;
    use crate::{
        BrepBody, FacePairId, GlobalArcRef, GraphArcId,
        HoleParameters, Point3, Solid, ThroughHoleBrep,
    };

    fn q(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
    fn block()->BrepModel {
        BrepBody::from(Solid::block(q(-3.,0.,-3.),6.,5.,6.)
            .unwrap()).shared().unwrap()
    }
    fn pierced()->BrepModel {
        let body=Solid::block(q(-3.,0.,-3.),6.,5.,6.).unwrap();
        BrepBody::from(ThroughHoleBrep::from_box(&body,HoleParameters{
            center_x:0.,center_z:0.,radius:1.,
        }).unwrap()).shared().unwrap()
    }
    fn trace(
        model:&BrepModel,start:Point3,end:Point3,
        intervals:&[(Point3,Point3)],pair_id:u32,
    )->TraceLine {
        let Surface3::Plane(plane)=model.face_surface(FaceId(1)).unwrap() else {
            panic!("synthetic support must be planar");
        };
        let a=plane.frame.project(start).unwrap();
        let b=plane.frame.project(end).unwrap();
        let dv=sub(b,a);
        let n=norm(dv);
        let direction=[dv[0]/n,dv[1]/n];
        let mut line=TraceLine {
            origin:a,direction,intervals:Vec::new(),
        };
        for (i,&(from,to)) in intervals.iter().enumerate() {
            let x=line.along(plane.frame.project(from).unwrap());
            let y=line.along(plane.frame.project(to).unwrap());
            line.intervals.push(SourceInterval{
                start:x.min(y),end:x.max(y),
                source:GlobalArcRef{
                    pair:FacePairId(pair_id),local:GraphArcId(i as u32),
                },
            });
        }
        line
    }
    fn assert_four(plan:MultiTraceUvPlan,area:f64) {
        assert_eq!(plan.cells.len(),4);
        assert_eq!(plan.cuts.len(),4);
        assert!(!plan.can_apply);
        for cell in &plan.cells {
            assert!((cell.area-area).abs()<1e-7,
                "unexpected quarter area: {} vs {area}",cell.area);
        }
        assert!((plan.cells.iter().map(|c|c.area).sum::<f64>()-4.*area)
            .abs()<1e-7);
    }
    #[test]
    fn orthogonal_crossings_form_four_true_face_cycles() {
        let model=block();
        let x=trace(&model,q(0.,5.,-3.),q(0.,5.,3.),
            &[(q(0.,5.,-3.),q(0.,5.,3.))],0);
        let z=trace(&model,q(-3.,5.,0.),q(3.,5.,0.),
            &[(q(-3.,5.,0.),q(3.,5.,0.))],1);
        let result=arrange_planar(&model,FaceId(1),GraphOperand::First,
            vec![FacePairId(0),FacePairId(1)],vec![x,z],1e-9).unwrap();
        let MultiTraceUvResult::Ready(plan)=result else {
            panic!("orthogonal crossing must form four cells: {result:?}");
        };
        assert_four(plan,9.0);
    }
    #[test]
    fn crossing_inside_a_true_hole_keeps_four_quarter_circle_arcs() {
        let model=pierced();
        let x=trace(&model,q(0.,5.,-3.),q(0.,5.,3.),&[
            (q(0.,5.,-3.),q(0.,5.,-1.)),
            (q(0.,5.,1.),q(0.,5.,3.)),
        ],0);
        let z=trace(&model,q(-3.,5.,0.),q(3.,5.,0.),&[
            (q(-3.,5.,0.),q(-1.,5.,0.)),
            (q(1.,5.,0.),q(3.,5.,0.)),
        ],1);
        let result=arrange_planar(&model,FaceId(1),GraphOperand::First,
            vec![FacePairId(0),FacePairId(1)],vec![x,z],1e-9).unwrap();
        let MultiTraceUvResult::Ready(plan)=result else {
            panic!("quarter-annulus cells must preserve circular arcs: {result:?}");
        };
        assert_eq!(plan.cells.iter().flat_map(|c|&c.cycles)
            .flat_map(|c|&c.spans).filter(|s|s.circular.is_some()).count(),4);
        assert_four(plan,9.0-PI/4.0);
    }
    #[test]
    fn identical_lines_are_not_accepted_as_two_independent_cuts() {
        let model=block();
        let line=trace(&model,q(0.,5.,-3.),q(0.,5.,3.),
            &[(q(0.,5.,-3.),q(0.,5.,3.))],0);
        let mut duplicate=line.clone();
        duplicate.intervals[0].source.pair=FacePairId(1);
        let result=arrange_planar(&model,FaceId(1),GraphOperand::First,
            vec![FacePairId(0),FacePairId(1)],
            vec![line,duplicate],1e-9).unwrap();
        assert_eq!(result,MultiTraceUvResult::Held(
            MultiTraceUvHold::OverlappingTraces));
    }
}
