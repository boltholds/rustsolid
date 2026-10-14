//! Exact-in-form 2D arrangement of ONE transverse straight intersection
//! trace across a convex planar B-rep face with disjoint circular inner loops.
//!
//! This creates complete oriented UV *cell boundaries* (including true
//! circular arcs), groups their holes, and reconciles both directed uses of
//! every new physical cut segment. It is a preview, not a topology commit.
//! No polygonal faceting of the source circles or proprietary algorithms.
use std::f64::consts::{PI, TAU};
use crate::{
    BrepModel, CoedgeId, CoedgeReconciliationPlan, EdgeId, FaceId,
    GeometryError, GeometryTolerance,
    GraphArcId, GraphOperand, GraphState, IntersectionCurve, IntersectionGraph,
    LoopId, LoopRole, ParameterRange, Point2Param, Point3, Surface3, TrimContact,
};

const MAX_HOLES:usize=64;
const MAX_FRAGMENTS:usize=512;
const SAMPLE_EPS_MULTIPLIER:f64=8.0;

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn p2(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn dist3(a:Point3,b:Point3)->f64 {
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

fn delta(a:Point2Param,b:Point2Param)->[f64;2]{[b.u-a.u,b.v-a.v]}
fn d(a:Point2Param,b:Point2Param)->f64{(a.u-b.u).hypot(a.v-b.v)}
fn cross(a:[f64;2],b:[f64;2])->f64{a[0]*b[1]-a[1]*b[0]}
fn cross_points(a:Point2Param,b:Point2Param)->f64{a.u*b.v-a.v*b.u}
fn uv_distance_along(origin:Point2Param,dir:[f64;2],point:Point2Param)->f64{
    cross(dir,[point.u-origin.u,point.v-origin.v])
}
fn lerp(a:Point2Param,b:Point2Param,t:f64)->Point2Param{
    p2(a.u+(b.u-a.u)*t,a.v+(b.v-a.v)*t)
}
fn param(range:ParameterRange,f:f64)->Result<f64,GeometryError>{
    match range {
        ParameterRange::Bounded{start,end} if start<end &&
            start.is_finite() && end.is_finite()=>Ok(start+(end-start)*f),
        _=>Err(invalid("UV arrangement requires bounded source coedge carrier")),
    }
}
fn loop_ids(model:&BrepModel,id:LoopId)->Result<Vec<CoedgeId>,GeometryError>{
    let first=model.loops.get(id.0 as usize)
        .ok_or_else(||invalid("missing trim loop"))?.first_coedge;
    let mut current=first;let mut out=Vec::new();
    for _ in 0..model.coedges.len(){
        out.push(current);
        current=model.next_in_loop(current)
            .ok_or_else(||invalid("broken source loop traversal"))?;
        if current==first{return Ok(out);}
    }
    Err(invalid("unbounded source loop traversal"))
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum UvArrangementHold{
    UnresolvedIntersection,
    MissingTransverseCut,
    UnsupportedSurface,
    UnsupportedSourceTrim,
    NonConvexOuterLoop,
    TangentOrVertexEvent,
    NonManifoldUvIncidence,
    DisconnectedRegion,
    InconsistentSourceArea,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum UvCutSide{Positive,Negative}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum UvCycleKind{Outer,Inner}

/// Source coedge fragments preserve the original UV parameter values.
/// Circular wraps may use an unwrapped negative angle; they refer to
/// the SAME source circular edge, not a new topological seam edge.
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum UvBoundaryProvenance{
    Source{edge:EdgeId,coedge:CoedgeId,loop_id:LoopId,
        parameter_start:f64,parameter_end:f64},
    Cut{arc:GraphArcId},
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct UvCircularSpan{
    pub center:Point2Param,
    pub radius:f64,
    pub angle_start:f64,
    pub angle_end:f64,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct UvBoundarySpan{
    pub start:Point2Param,
    pub end:Point2Param,
    pub provenance:UvBoundaryProvenance,
    pub circular:Option<UvCircularSpan>,
}
impl UvBoundarySpan{
    pub fn signed_area(self)->f64{
        if let Some(c)=self.circular{
            let ds=c.angle_end.sin()-c.angle_start.sin();
            let dc=c.angle_end.cos()-c.angle_start.cos();
            0.5*(c.radius*(c.center.u*ds-c.center.v*dc)
                +c.radius*c.radius*(c.angle_end-c.angle_start))
        }else{0.5*cross_points(self.start,self.end)}
    }
    pub fn midpoint(self)->Point2Param{
        if let Some(c)=self.circular{
            let angle=(c.angle_start+c.angle_end)*0.5;
            p2(c.center.u+c.radius*angle.cos(),c.center.v+c.radius*angle.sin())
        }else{lerp(self.start,self.end,0.5)}
    }
}
#[derive(Debug,Clone,PartialEq)]
pub struct UvBoundaryCycle{
    pub kind:UvCycleKind,
    pub spans:Vec<UvBoundarySpan>,
    pub signed_area:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct UvFaceCell{
    pub side:UvCutSide,
    pub cycles:Vec<UvBoundaryCycle>,
    pub area:f64,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct UvCutTwin{
    pub arc:GraphArcId,
    /// A genuine physical cut segment has one directed use in each UV cell.
    pub positive_start:Point2Param,
    pub positive_end:Point2Param,
    pub negative_start:Point2Param,
    pub negative_end:Point2Param,
}
#[derive(Debug,Clone,PartialEq)]
pub struct UvFaceRegionPlan{
    pub face:FaceId,
    pub operand:GraphOperand,
    pub cells:[UvFaceCell;2],
    pub cut_twins:Vec<UvCutTwin>,
    pub original_face_area:f64,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum UvFaceArrangement{
    Ready(UvFaceRegionPlan),
    Held(UvArrangementHold),
}
#[derive(Debug,Clone,Copy)]
struct Draft{
    side:UvCutSide,
    span:UvBoundarySpan,
}

struct TopologyAssembler{
    items:Vec<Draft>,
    eps:f64,
}
impl TopologyAssembler{
    fn new(eps:f64)->Self{Self{items:Vec::new(),eps}}
    fn add(&mut self,side:UvCutSide,span:UvBoundarySpan)->Result<(),GeometryError>{
        if self.items.len()>=MAX_FRAGMENTS {
            return Err(invalid("too many planar UV arrangement fragments"));
        }
        if !span.start.u.is_finite()||!span.start.v.is_finite()||
            !span.end.u.is_finite()||!span.end.v.is_finite(){
            return Err(invalid("nonfinite planar arrangement fragment"));
        }
        if d(span.start,span.end)<=self.eps && span.circular.is_none(){
            return Err(invalid("noncircular zero-length UV edge"));
        }
        self.items.push(Draft{side,span});
        Ok(())
    }
    fn classify_fragment(&mut self,span:UvBoundarySpan,
        origin:Point2Param,dir:[f64;2])->Result<(),UvArrangementHold>
    {
        let midpoint=span.midpoint();
        let signed=uv_distance_along(origin,dir,midpoint);
        if signed.abs()<=self.eps*SAMPLE_EPS_MULTIPLIER {
            return Err(UvArrangementHold::TangentOrVertexEvent);
        }
        self.add(if signed>0.0{UvCutSide::Positive}else{UvCutSide::Negative},span)
            .map_err(|_|UvArrangementHold::NonManifoldUvIncidence)
    }
    fn cycles(&self,side:UvCutSide)->Result<Vec<UvBoundaryCycle>,UvArrangementHold>{
        let filtered=self.items.iter().filter(|item|item.side==side)
            .map(|item|item.span).collect::<Vec<_>>();
        let mut vertices=Vec::<Point2Param>::new();
        let mut ids=Vec::<(usize,usize)>::new();
        for span in &filtered{
            let mut locate=|point:Point2Param|{
                if let Some((i,_))=vertices.iter().enumerate().find(|(_,x)|
                    d(**x,point)<=self.eps){
                    i
                }else{
                    let i=vertices.len();vertices.push(point);i
                }
            };
            let a=locate(span.start);let b=locate(span.end);
            ids.push((a,b));
        }
        let mut incoming=vec![Vec::<usize>::new();vertices.len()];
        let mut outgoing=vec![Vec::<usize>::new();vertices.len()];
        for (i,&(a,b)) in ids.iter().enumerate(){
            outgoing[a].push(i);incoming[b].push(i);
        }
        if incoming.iter().zip(outgoing.iter()).any(|(a,b)|a.len()!=1||b.len()!=1){
            return Err(UvArrangementHold::NonManifoldUvIncidence);
        }
        let mut visited=vec![false;filtered.len()];
        let mut out=Vec::new();
        for start in 0..filtered.len(){
            if visited[start]{continue;}
            let mut cycle=Vec::new();let mut cursor=start;
            loop {
                if visited[cursor]{
                    if cursor==start{break;}
                    return Err(UvArrangementHold::NonManifoldUvIncidence);
                }
                visited[cursor]=true;
                cycle.push(filtered[cursor]);
                let next=outgoing[ids[cursor].1][0];
                cursor=next;
                if cursor==start{break;}
                if cycle.len()>filtered.len(){
                    return Err(UvArrangementHold::NonManifoldUvIncidence);
                }
            }
            let signed_area=cycle.iter().map(|c|c.signed_area()).sum::<f64>();
            if !signed_area.is_finite()||signed_area.abs()<=self.eps*self.eps{
                return Err(UvArrangementHold::NonManifoldUvIncidence);
            }
            out.push(UvBoundaryCycle{
                kind:if signed_area>0.0{UvCycleKind::Outer}else{UvCycleKind::Inner},
                spans:cycle,signed_area,
            });
        }
        Ok(out)
    }
}

fn push_split_segment(assembler:&mut TopologyAssembler,
    provenance:UvBoundaryProvenance, a:Point2Param,b:Point2Param,
    origin:Point2Param,dir:[f64;2],
)->Result<usize,UvArrangementHold>{
    let sa=uv_distance_along(origin,dir,a);
    let sb=uv_distance_along(origin,dir,b);
    if sa.abs()<=assembler.eps*SAMPLE_EPS_MULTIPLIER ||
       sb.abs()<=assembler.eps*SAMPLE_EPS_MULTIPLIER{
        return Err(UvArrangementHold::TangentOrVertexEvent);
    }
    let mut intervals=vec![(0.0,1.0)];
    let crosses=sa.signum()!=sb.signum();
    if crosses{
        let root=sa/(sa-sb);
        if root<=assembler.eps/d(a,b) || root>=1.0-assembler.eps/d(a,b){
            return Err(UvArrangementHold::TangentOrVertexEvent);
        }
        intervals=vec![(0.0,root),(root,1.0)];
    }
    for (x,y) in intervals{
        let sub=UvBoundarySpan{
            start:lerp(a,b,x),end:lerp(a,b,y),
            provenance:match provenance{
                UvBoundaryProvenance::Source{edge,coedge,loop_id,parameter_start,parameter_end}=>{
                    UvBoundaryProvenance::Source{
                        edge,coedge,loop_id,
                        parameter_start:parameter_start+(parameter_end-parameter_start)*x,
                        parameter_end:parameter_start+(parameter_end-parameter_start)*y,
                    }
                }
                x=>x,
            },
            circular:None,
        };
        assembler.classify_fragment(sub,origin,dir)?;
    }
    Ok(usize::from(crosses))
}

fn polygon_area(vertices:&[Point2Param])->f64{
    vertices.iter().enumerate().map(|(i,p)|
        cross_points(*p,vertices[(i+1)%vertices.len()])).sum::<f64>()*0.5
}

fn arrange(first:&BrepModel,second:&BrepModel,graph:&IntersectionGraph,
    operand:GraphOperand,tol:GeometryTolerance,
)->Result<UvFaceArrangement,GeometryError>{
    tol.validate()?;
    let (a,b)=(graph.source_faces[0].face,graph.source_faces[1].face);
    graph.validate_against(first,a,second,b)?;
    if graph.state!=GraphState::Ready || graph.arcs.is_empty(){
        return Ok(UvFaceArrangement::Held(UvArrangementHold::UnresolvedIntersection));
    }
    let (model,face_id)=match operand{
        GraphOperand::First=>(first,a),
        GraphOperand::Second=>(second,b),
    };
    let source=model.faces.get(face_id.0 as usize)
        .ok_or_else(||invalid("missing source face for planar arrangement"))?;
    let Surface3::Plane(plane)=model.face_surface(face_id)
        .ok_or_else(||invalid("source face has no surface"))? else{
        return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSurface));
    };
    if source.loops.is_empty() || source.loops.len()>MAX_HOLES+1 {
        return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
    }
    let IntersectionCurve::Line(line)=graph.arcs[0].curve else{
        return Ok(UvFaceArrangement::Held(UvArrangementHold::MissingTransverseCut));
    };
    if graph.arcs.iter().any(|arc|arc.curve!=graph.arcs[0].curve||
        arc.closed||arc.contact!=TrimContact::Crossing){
        return Ok(UvFaceArrangement::Held(UvArrangementHold::MissingTransverseCut));
    }
    let origin=plane.frame.project(line.origin)?;
    let direct=[line.direction.dot(plane.frame.x),line.direction.dot(plane.frame.y)];
    let mag=direct[0].hypot(direct[1]);
    if !mag.is_finite() || mag<0.5{
        return Ok(UvFaceArrangement::Held(UvArrangementHold::MissingTransverseCut));
    }
    let dir=[direct[0]/mag,direct[1]/mag];
    let scale=model.bbox.max.sub(model.bbox.min);
    let extent=scale.x.hypot(scale.y).hypot(scale.z).max(1.0);
    let eps=tol.length_at(extent).max(graph.tolerance.length);
    let mut assembler=TopologyAssembler::new(eps);

    // The real outer loop: each original coedge and its exact parameter
    // interval is retained, but the fragments are oriented counterclockwise
    // in this planar UV chart so face material is always on their left.
    let outer=loop_ids(model,source.loops[0])?;
    if outer.len()<3||outer.len()>MAX_FRAGMENTS/2{
        return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
    }
    let mut edges=Vec::new();let mut points=Vec::new();
    for id in outer{
        let Some(crate::Curve2::Line(c))=model.coedge_pcurve(id) else{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        };
        let g=model.geometry.coedges[id.0 as usize];
        let ta=param(g.domain,0.0)?;let tb=param(g.domain,1.0)?;
        let pa=c.evaluate(ta)?;let pb=c.evaluate(tb)?;
        if d(pa,pb)<=eps*SAMPLE_EPS_MULTIPLIER{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        }
        points.push(pa);
        edges.push((id,pa,pb,ta,tb));
    }
    let sign=polygon_area(&points);
    if sign.abs()<=eps*eps{
        return Ok(UvFaceArrangement::Held(UvArrangementHold::NonConvexOuterLoop));
    }
    let ccw=sign>0.0;
    // Strict convexity is necessary for the two-cell certification below.
    let n=points.len();
    for i in 0..n{
        let a=points[i];let b=points[(i+1)%n];let c=points[(i+2)%n];
        let turn=cross(delta(a,b),delta(b,c));
        if !turn.is_finite()||turn.abs()<=eps*extent ||
            (turn>0.0)!=ccw{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::NonConvexOuterLoop));
        }
        if d(edges[i].2,edges[(i+1)%n].1)>eps*SAMPLE_EPS_MULTIPLIER{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        }
    }
    let original_outer=sign.abs();
    let mut outer_crossings=0;
    for (id,pa,pb,ta,tb) in edges{
        let c=&model.coedges[id.0 as usize];
        let span=UvBoundaryProvenance::Source{
            edge:c.edge,coedge:id,loop_id:c.loop_id,
            parameter_start:if ccw{ta}else{tb},
            parameter_end:if ccw{tb}else{ta},
        };
        outer_crossings+=push_split_segment(&mut assembler,span,
            if ccw{pa}else{pb},if ccw{pb}else{pa},origin,dir)
            .map_err(|_|invalid("outer boundary crosses cut at uncertain vertex/tangency"))?;
    }
    if outer_crossings!=2 {
        return Ok(UvFaceArrangement::Held(UvArrangementHold::MissingTransverseCut));
    }

    let mut total_holes_area=0.0;
    for &loop_id in source.loops.iter().skip(1){
        let Some(lp)=model.loops.get(loop_id.0 as usize) else{
            return Err(invalid("missing inner UV loop"));
        };
        if lp.role!=LoopRole::Inner || lp.face!=face_id {
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        }
        let ids=loop_ids(model,loop_id)?;
        if ids.len()!=1 {
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        }
        let id=ids[0];
        let Some(crate::Curve2::Circle(circle))=model.coedge_pcurve(id) else{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        };
        if circle.radius<=eps*SAMPLE_EPS_MULTIPLIER || !circle.radius.is_finite(){
            return Ok(UvFaceArrangement::Held(UvArrangementHold::UnsupportedSourceTrim));
        }
        total_holes_area+=PI*circle.radius*circle.radius;
        let c=&model.coedges[id.0 as usize];
        let provenance=|x:f64,y:f64|UvBoundaryProvenance::Source{
            edge:c.edge,coedge:id,loop_id,parameter_start:x,parameter_end:y};
        let normal=[-dir[1],dir[0]];
        let centre_offset=(circle.center.u-origin.u)*normal[0]
            +(circle.center.v-origin.v)*normal[1];
        let mut pairs=Vec::new();
        let target=-centre_offset/circle.radius;
        if (target.abs()-1.0).abs()<=eps*SAMPLE_EPS_MULTIPLIER/circle.radius {
            return Ok(UvFaceArrangement::Held(UvArrangementHold::TangentOrVertexEvent));
        }
        if target.abs()<1.0 {
            let phase=normal[1].atan2(normal[0]);
            let angle=target.acos();
            let mut roots=[
                (phase+angle).rem_euclid(TAU),
                (phase-angle).rem_euclid(TAU),
            ];
            roots.sort_by(f64::total_cmp);
            if roots.iter().any(|t|t.min(TAU-*t)<=eps/circle.radius){
                return Ok(UvFaceArrangement::Held(UvArrangementHold::TangentOrVertexEvent));
            }
            pairs.push((roots[1],roots[0]));
            pairs.push((roots[0],roots[1]-TAU));
        }else {
            pairs.push((0.0,-TAU));
        }
        for (start,end) in pairs{
            let pos=|angle:f64|p2(
                circle.center.u+circle.radius*angle.cos(),
                circle.center.v+circle.radius*angle.sin());
            let span=UvBoundarySpan{
                start:pos(start),end:pos(end),
                provenance:provenance(start,end),
                circular:Some(UvCircularSpan{
                    center:circle.center,radius:circle.radius,
                    angle_start:start,angle_end:end}),
            };
            assembler.classify_fragment(span,origin,dir)
                .map_err(|_|invalid("inner circular trim is tangent/uncertain"))?;
        }
    }

    let mut cut_twins=Vec::new();
    for arc in &graph.arcs{
        let p=plane.frame.project(arc.evaluate(0.0)?)?;
        let q=plane.frame.project(arc.evaluate(1.0)?)?;
        let sa=uv_distance_along(origin,dir,p);
        let sb=uv_distance_along(origin,dir,q);
        if sa.abs()>eps*SAMPLE_EPS_MULTIPLIER ||
           sb.abs()>eps*SAMPLE_EPS_MULTIPLIER ||
           d(p,q)<=eps*SAMPLE_EPS_MULTIPLIER{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::MissingTransverseCut));
        }
        let on_axis=|point:Point2Param| (point.u-origin.u)*dir[0]+(point.v-origin.v)*dir[1];
        let (start,end)=if on_axis(p)<on_axis(q){(p,q)}else{(q,p)};
        let pos=UvBoundarySpan{start,end,
            provenance:UvBoundaryProvenance::Cut{arc:arc.id},circular:None};
        let neg=UvBoundarySpan{start:end,end:start,
            provenance:UvBoundaryProvenance::Cut{arc:arc.id},circular:None};
        assembler.add(UvCutSide::Positive,pos)?;
        assembler.add(UvCutSide::Negative,neg)?;
        cut_twins.push(UvCutTwin{
            arc:arc.id,positive_start:start,positive_end:end,
            negative_start:end,negative_end:start,
        });
    }

    let mut cells=Vec::new();
    for side in [UvCutSide::Positive,UvCutSide::Negative]{
        let cycles=match assembler.cycles(side){
            Ok(v)=>v,
            Err(why)=>return Ok(UvFaceArrangement::Held(why)),
        };
        if cycles.iter().filter(|c|c.kind==UvCycleKind::Outer).count()!=1 {
            return Ok(UvFaceArrangement::Held(UvArrangementHold::DisconnectedRegion));
        }
        let area=cycles.iter().map(|c|c.signed_area).sum::<f64>();
        if !area.is_finite() || area<=eps*eps{
            return Ok(UvFaceArrangement::Held(UvArrangementHold::InconsistentSourceArea));
        }
        cells.push(UvFaceCell{side,cycles,area});
    }
    let area=original_outer-total_holes_area;
    let sum=cells[0].area+cells[1].area;
    if !area.is_finite()||area<=0.0 ||
        (area-sum).abs()>eps*extent*256.0 {
        return Ok(UvFaceArrangement::Held(UvArrangementHold::InconsistentSourceArea));
    }
    Ok(UvFaceArrangement::Ready(UvFaceRegionPlan{
        face:face_id,operand,
        cells:[cells.remove(0),cells.remove(0)],
        cut_twins,original_face_area:area,can_apply:false,
    }))
}

impl UvFaceRegionPlan {
    pub fn build(first:&BrepModel,second:&BrepModel,
        graph:&IntersectionGraph,operand:GraphOperand,
        tolerance:GeometryTolerance,
    )->Result<UvFaceArrangement,GeometryError>{
        arrange(first,second,graph,operand,tolerance)
    }
    pub fn certify_physical_twins(&self,
        first:&BrepModel,second:&BrepModel,
        graph:&IntersectionGraph,peers:&CoedgeReconciliationPlan,
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>
    {
        self.validate_against(first,second,graph,tolerance)?;
        peers.validate_against(first,second,graph,tolerance)?;
        if self.cut_twins.len()!=peers.edges.len() {
            return Err(invalid("UV cells and physical edge reconciliation disagree on arc count"));
        }
        let model=if self.operand==GraphOperand::First{first}else{second};
        let face=model.face_surface(self.face)
            .ok_or_else(||invalid("UV cell face lost its source support"))?;
        let bbox=model.bbox.max.sub(model.bbox.min);
        let extent=bbox.x.hypot(bbox.y).hypot(bbox.z).max(1.0);
        let eps=tolerance.length_at(extent).max(graph.tolerance.evaluation_budget)*8.0;
        for twin in &self.cut_twins {
            let edge=peers.edges.iter().find(|x|x.source_arc==twin.arc)
                .ok_or_else(||invalid("UV cut arc lacks a physical 3D edge candidate"))?;
            let endpoints=[face.evaluate(twin.positive_start.u,
                twin.positive_start.v)?.point,
                face.evaluate(twin.positive_end.u,twin.positive_end.v)?.point];
            let natural=dist3(endpoints[0],edge.start)<=eps &&
                dist3(endpoints[1],edge.end)<=eps;
            let reversed=dist3(endpoints[0],edge.end)<=eps &&
                dist3(endpoints[1],edge.start)<=eps;
            if !natural && !reversed {
                return Err(invalid("planar cut cell segment not co-located with physical 3D edge"));
            }
            if d(twin.positive_start,twin.negative_end)>eps ||
               d(twin.positive_end,twin.negative_start)>eps {
                return Err(invalid("opposite intra-face cut coedges do not share endpoints"));
            }
            let own=if self.operand==GraphOperand::First{&edge.uses[0]}
                else{&edge.uses[1]};
            if own.source.face!=self.face{
                return Err(invalid("reconciled physical edge refers to a different planar face"));
            }
            let uv_start=own.evaluate(if natural==own.reversed{1.0}else{0.0})?;
            let uv_end=own.evaluate(if natural==own.reversed{0.0}else{1.0})?;
            if d(uv_start,twin.positive_start)>eps ||
               d(uv_end,twin.positive_end)>eps {
                return Err(invalid("planar UV cut differs from the corresponding owned pcurve"));
            }
            for (side,required_start,required_end) in [
                (UvCutSide::Positive,twin.positive_start,twin.positive_end),
                (UvCutSide::Negative,twin.negative_start,twin.negative_end),
            ]{
                let cell=self.cells.iter().find(|x|x.side==side)
                    .ok_or_else(||invalid("UV cut has no material-side cell"))?;
                let matching=cell.cycles.iter().flat_map(|c|&c.spans)
                    .filter(|x|x.provenance==UvBoundaryProvenance::Cut{arc:twin.arc})
                    .collect::<Vec<_>>();
                if matching.len()!=1||
                   d(matching[0].start,required_start)>eps||
                   d(matching[0].end,required_end)>eps {
                    return Err(invalid("physical cut edge must have exactly one opposite use in each UV cell"));
                }
            }
        }
        Ok(())
    }
    pub fn validate_against(&self,first:&BrepModel,second:&BrepModel,
        graph:&IntersectionGraph,tolerance:GeometryTolerance)->Result<(),GeometryError>
    {
        let regenerated=arrange(first,second,graph,self.operand,tolerance)?;
        if regenerated!=UvFaceArrangement::Ready(self.clone()){
            return Err(invalid("UV face-cell arrangement changed or sources were tampered"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod distance_regression {
    use super::*;

    #[test]
    fn physical_endpoints_use_euclidean_distance_in_all_three_axes() {
        let a=Point3{x:10.0,y:-7.0,z:3.0};
        let b=Point3{x:8.0,y:-7.0,z:3.0};
        assert!((dist3(a,b)-2.0).abs()<1e-13);
        assert!((dist3(b,a)-2.0).abs()<1e-13);
        let c=Point3{x:8.0,y:-4.0,z:9.0};
        assert!((dist3(a,c)-7.0).abs()<1e-13);
        assert!((dist3(c,a)-7.0).abs()<1e-13);
        assert_eq!(dist3(a,a),0.0);
    }
}
