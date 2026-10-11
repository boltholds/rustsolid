//! Immutable, first-party intersection event graph over validated B-rep faces.
//!
//! The graph records analytic arcs and contacts with typed provenance for both
//! operands, boundary-edge witnesses, and local tolerance evidence. It never
//! changes either input B-rep and is not a Boolean/TopologyDelta executor.
//! Closed periodic intersections retain one SYNTHETIC parameter anchor and a
//! self-loop; the periodic seam is not promoted to a physical cut.
use crate::{
    BrepModel, CoedgeId, Curve3, EdgeId, FaceId, GeometryError,
    GeometryTolerance, IntersectionCertificate, IntersectionCurve, ParameterRange,
    Point3, Surface3, SurfaceId, TrimAmbiguity, TrimContact, TrimPosition,
    TrimmedBranch, TrimmedIntersection, TrimmedIntersectionEngine,
};
use std::collections::BTreeSet;

fn invalid(msg: impl Into<String>) -> GeometryError { GeometryError::InvalidTopology(msg.into()) }
fn length(v: Point3) -> f64 { v.x.hypot(v.y).hypot(v.z) }
fn distance(a: Point3, b: Point3) -> f64 { length(a.sub(b)) }
fn extent(model: &BrepModel) -> f64 { length(model.bbox.max.sub(model.bbox.min)) }

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GraphOperand { First, Second }
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GraphNodeId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GraphArcId(pub u32);

/// Model-local face/surface identity. `GraphOperand` is essential: FaceId(0)
/// in one body MUST NOT be confused with FaceId(0) in another body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphFaceSource {
    pub operand: GraphOperand,
    pub face: FaceId,
    pub surface: SurfaceId,
    pub tolerance: GeometryTolerance,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GraphBoundarySource {
    pub operand: GraphOperand,
    pub face: FaceId,
    pub edge: EdgeId,
    pub coedge: CoedgeId,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphTolerance {
    /// Engineering/geometric merging tolerance determined from local feature
    /// extents and the engine + both face policies (never global origin size).
    pub length: f64,
    pub angular: f64,
    /// Floating point representability and sampling allowance; NEVER used
    /// to silently fuse distinct geometric nodes.
    pub evaluation_budget: f64,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphNodeKind {
    Endpoint,
    IsolatedContact,
    /// A closed curve needs a graph anchor, but its periodic seam is not a cut.
    PeriodicAnchor,
}
#[derive(Debug, Clone, PartialEq)]
pub struct GraphNode {
    pub id: GraphNodeId,
    pub point: Point3,
    pub kind: GraphNodeKind,
    /// Contact kind is retained for isolated touch points (and possible
    /// overlapping events merged with another branch endpoint).
    pub contact: Option<TrimContact>,
    /// [first face, second face] at this point (both must not be Outside).
    pub membership: [TrimPosition; 2],
    pub boundary_sources: Vec<GraphBoundarySource>,
    pub tolerance: GraphTolerance,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphArc {
    pub id: GraphArcId,
    pub from: GraphNodeId,
    pub to: GraphNodeId,
    pub curve: IntersectionCurve,
    pub parameter_start: f64,
    pub parameter_end: f64,
    pub closed: bool,
    pub contact: TrimContact,
    pub source_faces: [GraphFaceSource; 2],
    pub tolerance: GraphTolerance,
}
impl GraphArc {
    pub fn evaluate(self, fraction: f64) -> Result<Point3, GeometryError> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("graph arc evaluation fraction outside [0,1]"));
        }
        self.curve.evaluate(self.parameter_start+(self.parameter_end-self.parameter_start)*fraction)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphState { Empty, Ready, Ambiguous(TrimAmbiguity) }
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionGraph {
    pub source_faces: [GraphFaceSource; 2],
    pub state: GraphState,
    pub nodes: Vec<GraphNode>,
    pub arcs: Vec<GraphArc>,
    /// Sample certificate for the underlying unbounded support intersection;
    /// not a proof of trimmed-interval completeness.
    pub support_certificate: Option<IntersectionCertificate>,
    pub tolerance: GraphTolerance,
}
impl IntersectionGraph {
    pub fn kind(&self) -> &'static str {
        match self.state {GraphState::Empty=>"empty", GraphState::Ready=>"ready",
            GraphState::Ambiguous(_)=>"ambiguous"}
    }
    pub fn incident_arcs(&self, node: GraphNodeId) -> Vec<GraphArcId> {
        self.arcs.iter().filter(|a|a.from==node || a.to==node).map(|a|a.id).collect()
    }
    pub fn degree(&self, node: GraphNodeId) -> usize {
        self.arcs.iter().map(|a|usize::from(a.from==node)+usize::from(a.to==node)).sum()
    }
    /// Verify the graph against the actual source geometry, not merely its own
    /// stored coordinates. Invalid/tampered arc, provenance or trim fails closed.
    pub fn validate_against(&self, a:&BrepModel, fa:FaceId, b:&BrepModel, fb:FaceId)
        -> Result<(),GeometryError>
    {
        a.validate()?;b.validate()?;
        let sources=[face_source(a,fa,GraphOperand::First)?,
            face_source(b,fb,GraphOperand::Second)?];
        if self.source_faces!=sources {return Err(invalid("graph source face identity or tolerance changed"));}
        if !(self.tolerance.length.is_finite() && self.tolerance.length>0.0 &&
             self.tolerance.evaluation_budget.is_finite() &&
             self.tolerance.evaluation_budget>=self.tolerance.length) {
            return Err(invalid("invalid graph tolerance policy"));
        }
        match self.state {
            GraphState::Empty|GraphState::Ambiguous(_) => {
                if !self.nodes.is_empty() || !self.arcs.is_empty() || self.support_certificate.is_some(){
                    return Err(invalid("unresolved/empty intersections must not expose cut topology"));
                }
                return Ok(());
            }
            GraphState::Ready=>{
                if self.nodes.is_empty() || self.support_certificate.is_none(){
                    return Err(invalid("ready graph lacks contact topology or support evidence"));
                }
            }
        }
        let mut seen=BTreeSet::new();
        for (i,node) in self.nodes.iter().enumerate() {
            if node.id.0 as usize!=i || !node.point.is_finite() ||
                node.tolerance!=self.tolerance ||
                (node.kind==GraphNodeKind::IsolatedContact && node.contact.is_none()) {
                return Err(invalid("invalid graph vertex identity or position"));
            }
            for existing in &self.nodes[..i] {
                if distance(node.point,existing.point)<=self.tolerance.length {
                    return Err(invalid("duplicate intersection vertex within geometric tolerance"));
                }
            }
            let membership=classify_both(a,fa,b,fb,node.point,self.tolerance)?;
            if membership!=node.membership || membership.contains(&TrimPosition::Outside) {
                return Err(invalid("graph node no longer lies within both trimmed faces"));
            }
            let expected=boundary_witnesses(a,fa,b,fb,node.point,membership,self.tolerance)?;
            if node.boundary_sources!=expected || !seen.insert(node.id){
                return Err(invalid("graph boundary-edge provenance is inconsistent"));
            }
            if node.kind==GraphNodeKind::PeriodicAnchor && self.degree(node.id)!=2 {
                return Err(invalid("periodic anchor is not a closed analytic loop"));
            }
            certify_point(a,fa,b,fb,node.point,self.tolerance.evaluation_budget)?;
        }
        for (i,arc) in self.arcs.iter().enumerate() {
            if arc.id.0 as usize!=i || arc.source_faces!=sources ||
                arc.tolerance!=self.tolerance ||
                !arc.parameter_start.is_finite() || !arc.parameter_end.is_finite() ||
                arc.parameter_start>=arc.parameter_end {
                return Err(invalid("invalid graph arc identity, provenance or interval"));
            }
            let start=self.nodes.get(arc.from.0 as usize).filter(|n|n.id==arc.from)
                .ok_or_else(||invalid("graph arc start missing"))?;
            let end=self.nodes.get(arc.to.0 as usize).filter(|n|n.id==arc.to)
                .ok_or_else(||invalid("graph arc end missing"))?;
            if arc.closed != (arc.from==arc.to) ||
                (arc.closed && start.kind!=GraphNodeKind::PeriodicAnchor) {
                return Err(invalid("invalid periodic self-loop topology"));
            }
            if !arc.closed && distance(start.point,end.point)<=self.tolerance.length {
                return Err(invalid("open graph arc collapses within tolerance"));
            }
            for f in [0.0,0.125,0.25,0.5,0.75,0.875,1.0] {
                let p=arc.evaluate(f)?;
                if (f==0.0 && distance(p,start.point)>self.tolerance.evaluation_budget) ||
                   (f==1.0 && distance(p,end.point)>self.tolerance.evaluation_budget) {
                    return Err(invalid("graph curve does not meet its stored endpoints"));
                }
                let positions=classify_both(a,fa,b,fb,p,self.tolerance)?;
                if positions.contains(&TrimPosition::Outside) {
                    return Err(invalid("graph arc leaves selected B-rep trim"));
                }
                certify_point(a,fa,b,fb,p,self.tolerance.evaluation_budget)?;
                if f>0.0 && f<1.0 && arc.contact==TrimContact::BoundaryOverlap &&
                    positions.iter().all(|x|*x!=TrimPosition::Boundary) {
                    return Err(invalid("boundary-overlap arc is not on a face boundary"));
                }
            }
        }
        for node in &self.nodes {
            if node.kind==GraphNodeKind::IsolatedContact && self.degree(node.id)!=0 {
                return Err(invalid("isolated point contact was wired to a curve"));
            }
            if node.kind==GraphNodeKind::Endpoint && self.degree(node.id)==0 {
                return Err(invalid("unreferenced branch endpoint"));
            }
        }
        Ok(())
    }
}

fn face_source(model:&BrepModel,face:FaceId,operand:GraphOperand)
    ->Result<GraphFaceSource,GeometryError>
{
    let src=model.faces.get(face.0 as usize).filter(|f|f.id==face)
        .ok_or_else(||invalid("graph references unknown B-rep face"))?;
    src.geometry.tolerance.validate()?;
    Ok(GraphFaceSource{operand,face,surface:src.geometry.surface,tolerance:src.geometry.tolerance})
}
fn tolerance_policy(a:&BrepModel, fa:FaceId,b:&BrepModel,fb:FaceId,
    engine:GeometryTolerance)->Result<GraphTolerance,GeometryError>
{
    let sa=face_source(a,fa,GraphOperand::First)?;
    let sb=face_source(b,fb,GraphOperand::Second)?;
    let feature_extent=extent(a).max(extent(b));
    if !feature_extent.is_finite(){return Err(invalid("B-rep extent overflow"));}
    let linear=engine.length_at(feature_extent)
        .max(sa.tolerance.length_at(extent(a))).max(sb.tolerance.length_at(extent(b)));
    let max_coord=[a.bbox.min,a.bbox.max,b.bbox.min,b.bbox.max].into_iter()
        .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()]).fold(1.0f64,f64::max);
    let roundoff=64.0*f64::EPSILON*max_coord;
    if !linear.is_finite() || !roundoff.is_finite(){return Err(invalid("invalid graph precision"));}
    Ok(GraphTolerance{length:linear,angular:engine.angular.max(sa.tolerance.angular)
        .max(sb.tolerance.angular),evaluation_budget:linear.max(roundoff)})
}
fn classify_both(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId,
    point:Point3,tol:GraphTolerance)->Result<[TrimPosition;2],GeometryError>
{
    let policy=GeometryTolerance{absolute_length:tol.length,relative_length:0.0,angular:tol.angular};
    Ok([crate::trimmed_intersections::classify_trim_point(a,fa,point,policy)?,
        crate::trimmed_intersections::classify_trim_point(b,fb,point,policy)?])
}
fn on_edge(model:&BrepModel, edge:EdgeId,p:Point3,epsilon:f64)->Result<bool,GeometryError>{
    let curve=model.edge_curve(edge).ok_or_else(||invalid("boundary edge lacks 3D curve"))?;
    let ge=&model.geometry.edges[edge.0 as usize];
    let ParameterRange::Bounded{start,end}=ge.domain else {return Err(invalid("unbounded B-rep edge carrier"))};
    match curve {
        Curve3::Line(line)=>{
            let d=line.direction;let dd=d.dot(d);
            if dd<=0.0 || !dd.is_finite(){return Err(invalid("invalid edge line direction"));}
            let t=p.sub(line.origin).dot(d)/dd;
            let delta=epsilon/dd.sqrt();
            if t<start-delta || t>end+delta {return Ok(false);}
            Ok(distance(line.evaluate(t.clamp(start,end))?,p)<=epsilon)
        }
        Curve3::Circle(circle)=>{
            let v=p.sub(circle.frame.origin);
            let normal=v.dot(circle.frame.normal);
            if normal.abs()>epsilon{return Ok(false);}
            let x=v.dot(circle.frame.x);let y=v.dot(circle.frame.y);
            if (x.hypot(y)-circle.radius).abs()>epsilon {return Ok(false);}
            let angle=y.atan2(x).rem_euclid(std::f64::consts::TAU);
            let angular_slack=epsilon/circle.radius;
            Ok(angle>=start-angular_slack && angle<=end+angular_slack)
        }
    }
}
fn boundary_witnesses(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId,
    point:Point3,membership:[TrimPosition;2],tol:GraphTolerance)
    ->Result<Vec<GraphBoundarySource>,GeometryError>
{
    let mut out=BTreeSet::new();
    for (model,face,operand,position) in [(a,fa,GraphOperand::First,membership[0]),
                                          (b,fb,GraphOperand::Second,membership[1])] {
        if position!=TrimPosition::Boundary {continue;}
        let face_object=&model.faces[face.0 as usize];
        for loop_id in &face_object.loops {
            let start=model.loops[loop_id.0 as usize].first_coedge;
            let mut co=start;
            for _ in 0..=model.coedges.len() {
                let use_=&model.coedges[co.0 as usize];
                let e=&model.edges[use_.edge.0 as usize];
                let twins=e.coedges;
                // A periodic seam is a UV chart artifact, not a physical
                // boundary of the selected cylindrical side face.
                let seam=model.coedges[twins[0].0 as usize].face==
                         model.coedges[twins[1].0 as usize].face;
                if !seam && on_edge(model,e.id,point,tol.evaluation_budget)? {
                    out.insert(GraphBoundarySource{operand,face,edge:e.id,coedge:co});
                }
                co=use_.next;
                if co==start{break;}
            }
        }
    }
    Ok(out.into_iter().collect())
}
fn implicit_distance(surface:Surface3,p:Point3)->Result<f64,GeometryError>{
    let d=match surface {
        Surface3::Plane(plane)=>plane.frame.normal.dot(p.sub(plane.frame.origin)).abs(),
        Surface3::Cylinder(c)=>{
            let v=p.sub(c.frame.origin);let z=v.dot(c.frame.normal);
            (length(v.sub(c.frame.normal.scale(z)))-c.radius).abs()
        }
        Surface3::Sphere(_)=>return Err(invalid("sphere intersection graph is not implemented")),
    };
    if d.is_finite(){Ok(d)}else{Err(invalid("non-finite implicit surface residual"))}
}
fn certify_point(a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId,
    point:Point3,budget:f64)->Result<(),GeometryError>
{
    let first=a.face_surface(fa).ok_or_else(||invalid("missing first support"))?;
    let second=b.face_surface(fb).ok_or_else(||invalid("missing second support"))?;
    if implicit_distance(first,point)? > budget || implicit_distance(second,point)? > budget {
        return Err(invalid("intersection graph point does not lie on both surfaces"));
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub struct IntersectionGraphBuilder {pub tolerance:GeometryTolerance}
impl IntersectionGraphBuilder {
    pub fn new(tolerance:GeometryTolerance)->Result<Self,GeometryError>{
        tolerance.validate()?;Ok(Self{tolerance})
    }
    pub fn intersect_faces(&self,a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId)
        ->Result<IntersectionGraph,GeometryError>
    {
        let engine=TrimmedIntersectionEngine::new(self.tolerance)?;
        let intersection=engine.intersect_faces(a,fa,b,fb)?;
        let sources=[face_source(a,fa,GraphOperand::First)?,face_source(b,fb,GraphOperand::Second)?];
        let precision=tolerance_policy(a,fa,b,fb,self.tolerance)?;
        let mut graph=IntersectionGraph{source_faces:sources,state:GraphState::Empty,
            nodes:vec![],arcs:vec![],support_certificate:None,tolerance:precision};
        match intersection {
            TrimmedIntersection::Empty=>return Ok(graph),
            TrimmedIntersection::Ambiguous{reason}=>{
                graph.state=GraphState::Ambiguous(reason);return Ok(graph);
            }
            TrimmedIntersection::Geometry{branches,contacts,certificate}=>{
                graph.state=GraphState::Ready;graph.support_certificate=Some(certificate);
                let branches=coalesce_periodic_seam(branches,precision.length);
                for br in branches {
                    if graph.arcs.len()>=256 {return Err(invalid("intersection graph exceeds arc budget"));}
                    if !br.closed && distance(br.start,br.end)<=precision.length {
                        // Never promote a numerically collapsed segment to
                        // a topological cut edge without further refinement.
                        graph.state=GraphState::Ambiguous(TrimAmbiguity::UnresolvedBoundaryEvent);
                        graph.nodes.clear();graph.arcs.clear();graph.support_certificate=None;
                        return Ok(graph);
                    }
                    let kind=if br.closed {GraphNodeKind::PeriodicAnchor}else{GraphNodeKind::Endpoint};
                    let start=insert_node(&mut graph,a,fa,b,fb,br.start,kind,None)?;
                    let end=if br.closed {start}else{insert_node(&mut graph,a,fa,b,fb,br.end,kind,None)?};
                    let id=GraphArcId(graph.arcs.len() as u32);
                    graph.arcs.push(GraphArc{id,from:start,to:end,curve:br.carrier,
                        parameter_start:br.parameter_start,parameter_end:br.parameter_end,
                        closed:br.closed,contact:br.contact,source_faces:sources,tolerance:precision});
                }
                for c in contacts {
                    insert_node(&mut graph,a,fa,b,fb,c.point,GraphNodeKind::IsolatedContact,Some(c.contact))?;
                }
            }
        }
        graph.validate_against(a,fa,b,fb)?;
        Ok(graph)
    }
}
fn insert_node(graph:&mut IntersectionGraph,a:&BrepModel,fa:FaceId,b:&BrepModel,fb:FaceId,
    point:Point3,kind:GraphNodeKind,contact:Option<TrimContact>)->Result<GraphNodeId,GeometryError>
{
    if !point.is_finite(){return Err(invalid("non-finite graph endpoint"));}
    for node in &mut graph.nodes {
        if distance(point,node.point)<=graph.tolerance.length {
            if kind==GraphNodeKind::PeriodicAnchor || node.kind==GraphNodeKind::PeriodicAnchor {
                return Err(invalid("periodic anchor coincides with another unrelated event"));
            }
            if kind==GraphNodeKind::Endpoint && node.kind==GraphNodeKind::IsolatedContact {
                node.kind=GraphNodeKind::Endpoint;
            }
            if node.contact.is_none() {node.contact=contact;}
            return Ok(node.id);
        }
    }
    if graph.nodes.len()>=512 {return Err(invalid("intersection graph exceeds vertex budget"));}
    let membership=classify_both(a,fa,b,fb,point,graph.tolerance)?;
    if membership.contains(&TrimPosition::Outside){
        return Err(invalid("graph vertex is outside one source trim"));
    }
    certify_point(a,fa,b,fb,point,graph.tolerance.evaluation_budget)?;
    let boundaries=boundary_witnesses(a,fa,b,fb,point,membership,graph.tolerance)?;
    let id=GraphNodeId(graph.nodes.len() as u32);
    graph.nodes.push(GraphNode{id,point,kind,contact,membership,
        boundary_sources:boundaries,tolerance:graph.tolerance});
    Ok(id)
}

/// Periodic UV charts cut one physical arc at u=0/TAU. Merge first and last
/// fragments into a single unwrapped analytic parameter interval, rather than
/// introducing an artificial topological split at the chart seam.
fn coalesce_periodic_seam(mut branches:Vec<TrimmedBranch>,epsilon:f64)->Vec<TrimmedBranch>{
    if branches.len()<2 {return branches;}
    let period=std::f64::consts::TAU;
    let first=branches[0];
    let last=*branches.last().expect("at least two branches");
    if !matches!(first.carrier,IntersectionCurve::PlaneCylinderSection(_)) ||
        first.carrier!=last.carrier || first.closed || last.closed ||
        first.contact!=last.contact ||
        first.parameter_start.abs()>1.0e-9 ||
        (last.parameter_end-period).abs()>1.0e-9 ||
        distance(first.start,last.end)>epsilon {return branches;}
    let replacement=TrimmedBranch{
        carrier:first.carrier, parameter_start:last.parameter_start,
        parameter_end:first.parameter_end+period, start:last.start,
        end:first.end, contact:first.contact, closed:false,
    };
    branches.pop();
    branches.remove(0);
    branches.push(replacement);
    branches
}
