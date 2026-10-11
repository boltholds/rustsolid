//! Deterministic, immutable, source-bound topology change *intent*.
//!
//! Converts a validated analytic IntersectionGraph into symbolic vertex, edge
//! and face work orders. Nothing in this module changes a B-rep, allocates a
//! persistent entity ID, asserts region-in/out ownership or commits a Boolean.
//! Every retained reference is operand-namespaced. Preview plans are checked
//! by regenerating them from the currently validated geometric graph.
use crate::{
    BrepModel, CoedgeId, Curve3, EdgeId, FaceId, GeometryError, GraphArcId,
    GraphBoundarySource, GraphFaceSource, GraphNode, GraphNodeId, GraphNodeKind,
    GraphOperand, GraphState, IntersectionCurve, IntersectionGraph, ParameterRange,
    Point3, TrimAmbiguity, TrimContact, TrimPosition, VertexId,
};
use std::collections::BTreeSet;
use std::f64::consts::TAU;

fn invalid(s:impl Into<String>)->GeometryError {GeometryError::InvalidTopology(s.into())}
fn length(p:Point3)->f64 {p.x.hypot(p.y).hypot(p.z)}
fn distance(a:Point3,b:Point3)->f64 {length(a.sub(b))}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum PlanState {
    /// Empty intersection or isolated contact: no topological cut is proposed.
    NoChange,
    /// Geometrically checked candidates; never permission to mutate geometry.
    PreviewOnly,
    /// Geometric or topological degeneracy requires a specialized resolver.
    Held,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord)]
pub enum PlanRequirement {
    NumericalRefinement,
    GeneralizedBrepEditor,
    RegionClassification,
    InnerLoopConstruction,
    OpenTraceCompletion,
    CurvedEdgeSubdivision,
    MultipleSplitsSameEdge,
    TangencyClassification,
    BoundaryReconciliation,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum PlanArcRole {ProposeNewCarrier, ReconcileBoundary, ContactOnly}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum PlanFaceEffect {
    /// Endpoints belong to the actual physical boundary of this face.
    CandidateBoundarySplit,
    /// An endpoint is interior, so an isolated open arc cannot split the face.
    PendingOpenTrace,
    /// Closed curve can introduce an inner trim loop or a new region.
    CandidateInnerLoop,
    /// Intersection runs along an original face boundary: keep topology until reconciliation.
    ReconcileExistingBoundary,
    /// Tangent curves are not automatically classified as cutting edges.
    RecordTangency,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum PlanVertexPlacement {
    /// No physical vertex is introduced at the parameter origin of a closed loop.
    VirtualPeriodicAnchor,
    /// An isolated contact is evidence, not a request to split geometry.
    ContactOnly,
    /// Reuse an original vertex on this operand. No new vertex should be minted.
    ReuseVertex {vertex:VertexId},
    /// Proposed insertion on an original 3D edge carrier, before local rewiring.
    SplitEdge {edge:EdgeId,coedge:CoedgeId,parameter:f64,fraction:f64},
    /// A cut endpoint not yet connected to this face's physical boundary.
    CreateFaceVertex {face:FaceId},
}
#[derive(Debug,Clone,PartialEq)]
pub struct PlannedVertex {
    pub node:GraphNodeId,
    pub operand:GraphOperand,
    pub face:FaceId,
    pub point:Point3,
    pub placement:PlanVertexPlacement,
    pub physical_boundary_sources:Vec<GraphBoundarySource>,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct PlannedArc {
    pub arc:GraphArcId,
    pub from:GraphNodeId,
    pub to:GraphNodeId,
    pub contact:TrimContact,
    pub closed:bool,
    pub analytic_curve:IntersectionCurve,
    pub parameter_start:f64,
    pub parameter_end:f64,
    pub source_faces:[GraphFaceSource;2],
    pub role:PlanArcRole,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct PlannedFaceChange {
    pub arc:GraphArcId,
    pub operand:GraphOperand,
    pub face:FaceId,
    pub effect:PlanFaceEffect,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TopologyDeltaPlan {
    pub state:PlanState,
    pub source_faces:[GraphFaceSource;2],
    pub original_graph_state:GraphState,
    pub vertices:Vec<PlannedVertex>,
    pub arcs:Vec<PlannedArc>,
    pub face_changes:Vec<PlannedFaceChange>,
    pub requirements:Vec<PlanRequirement>,
    /// No persistent entity ids are generated, no elements are retired and
    /// nothing can be applied while the common B-rep is still read-only.
    pub can_apply:bool,
}
impl TopologyDeltaPlan {
    /// Validation rederives all placements, sources and requirements from the
    /// original mathematically checked graph and *current* source B-reps.
    /// It also catches accidental or malicious mutation of the preview plan.
    pub fn validate_against(&self, graph:&IntersectionGraph,
        first:&BrepModel,first_face:FaceId,second:&BrepModel,second_face:FaceId)
        -> Result<(),GeometryError>
    {
        let expected=TopologyDeltaPlanner::new().plan(graph,first,first_face,second,second_face)?;
        if self!=&expected {return Err(invalid("topology delta plan differs from regenerated source-bound intent"));}
        Ok(())
    }
    pub fn split_edge_requests(&self)->usize {
        self.vertices.iter().filter(|v|matches!(v.placement,PlanVertexPlacement::SplitEdge{..})).count()
    }
    pub fn is_mutation_ready(&self)->bool {self.can_apply}
}

#[derive(Debug,Clone,Copy,Default)]
pub struct TopologyDeltaPlanner;
impl TopologyDeltaPlanner {
    pub fn new()->Self {Self}

    pub fn plan(&self,graph:&IntersectionGraph,
        first:&BrepModel,first_face:FaceId,second:&BrepModel,second_face:FaceId)
        -> Result<TopologyDeltaPlan,GeometryError>
    {
        graph.validate_against(first,first_face,second,second_face)?;
        let mut plan=TopologyDeltaPlan{
            state:PlanState::NoChange,
            source_faces:graph.source_faces,
            original_graph_state:graph.state,
            vertices:vec![],arcs:vec![],face_changes:vec![],
            requirements:vec![],can_apply:false,
        };
        match graph.state {
            GraphState::Empty=>return Ok(plan),
            GraphState::Ambiguous(_)=>{
                plan.state=PlanState::Held;
                plan.requirements.push(PlanRequirement::NumericalRefinement);
                return Ok(plan);
            }
            GraphState::Ready=>{},
        }
        // Do not generate independent source vertices for mere touch points.
        // A full curve tangent to a surface is also not assumed to cut it.
        for node in &graph.nodes {
            let incident=graph.incident_arcs(node.id);
            let has_cut=incident.iter().any(|id|graph.arcs[id.0 as usize].contact==TrimContact::Crossing);
            for (source_index,(body,face,operand)) in [
                (first,first_face,GraphOperand::First),
                (second,second_face,GraphOperand::Second),
            ].into_iter().enumerate() {
                let witnesses=node.boundary_sources.iter().copied()
                    .filter(|w|w.operand==operand && w.face==face).collect::<Vec<_>>();
                let placement=if node.kind==GraphNodeKind::PeriodicAnchor {
                    PlanVertexPlacement::VirtualPeriodicAnchor
                } else if node.kind==GraphNodeKind::IsolatedContact || !has_cut {
                    PlanVertexPlacement::ContactOnly
                } else {
                    determine_placement(body,face,node,operand,
                        graph.tolerance.length,&witnesses)?
                };
                if node.membership[source_index]==TrimPosition::Boundary && witnesses.is_empty(){
                    return Err(invalid("graph boundary node lacks physical source witnesses"));
                }
                plan.vertices.push(PlannedVertex{node:node.id,operand,face,
                    point:node.point,placement,physical_boundary_sources:witnesses});
            }
        }
        let mut needs=BTreeSet::<PlanRequirement>::new();
        let mut serious=false;
        for arc in &graph.arcs {
            let role=match arc.contact {
                TrimContact::Crossing=>PlanArcRole::ProposeNewCarrier,
                TrimContact::BoundaryOverlap=>{
                    needs.insert(PlanRequirement::BoundaryReconciliation);
                    serious=true;PlanArcRole::ReconcileBoundary
                }
                TrimContact::Tangent=>{
                    needs.insert(PlanRequirement::TangencyClassification);
                    serious=true;PlanArcRole::ContactOnly
                }
            };
            plan.arcs.push(PlannedArc{arc:arc.id,from:arc.from,to:arc.to,
                contact:arc.contact,closed:arc.closed,analytic_curve:arc.curve,
                parameter_start:arc.parameter_start,parameter_end:arc.parameter_end,
                source_faces:arc.source_faces,role});
            for (idx,source) in arc.source_faces.iter().enumerate() {
                let effect=match role {
                    PlanArcRole::ReconcileBoundary=>PlanFaceEffect::ReconcileExistingBoundary,
                    PlanArcRole::ContactOnly=>PlanFaceEffect::RecordTangency,
                    PlanArcRole::ProposeNewCarrier=>{
                        needs.insert(PlanRequirement::GeneralizedBrepEditor);
                        needs.insert(PlanRequirement::RegionClassification);
                        if arc.closed {
                            needs.insert(PlanRequirement::InnerLoopConstruction);
                            PlanFaceEffect::CandidateInnerLoop
                        } else {
                            let from=&graph.nodes[arc.from.0 as usize];
                            let to=&graph.nodes[arc.to.0 as usize];
                            if from.membership[idx]==TrimPosition::Boundary &&
                               to.membership[idx]==TrimPosition::Boundary {
                                PlanFaceEffect::CandidateBoundarySplit
                            }else{
                                needs.insert(PlanRequirement::OpenTraceCompletion);
                                PlanFaceEffect::PendingOpenTrace
                            }
                        }
                    }
                };
                plan.face_changes.push(PlannedFaceChange{arc:arc.id,
                    operand:source.operand,face:source.face,effect});
            }
        }
        let mut same_edge=BTreeSet::new();
        for v in &plan.vertices {
            if let PlanVertexPlacement::SplitEdge{edge,..}=v.placement {
                if !same_edge.insert((v.operand,edge)){
                    needs.insert(PlanRequirement::MultipleSplitsSameEdge);
                }
                let model=if v.operand==GraphOperand::First {first}else{second};
                if matches!(model.edge_curve(edge),Some(Curve3::Circle(_))) {
                    needs.insert(PlanRequirement::CurvedEdgeSubdivision);
                }
            }
        }
        if serious {plan.state=PlanState::Held;}
        else if !plan.arcs.is_empty() {plan.state=PlanState::PreviewOnly;}
        plan.requirements=needs.into_iter().collect();
        Ok(plan)
    }
}

/// A relative parameter along the actual 3D carrier of a physical edge.
/// Constructed values are f64 proposals and cannot be committed without a
/// dedicated local tolerance/canonical-reparameterization check.
fn edge_parameter(model:&BrepModel,edge:EdgeId,p:Point3,epsilon:f64)
    ->Result<(f64,f64),GeometryError>
{
    let e=model.edges.get(edge.0 as usize).filter(|e|e.id==edge)
        .ok_or_else(||invalid("plan refers to missing source edge"))?;
    let binding=model.geometry.edges.get(edge.0 as usize)
        .ok_or_else(||invalid("source edge lacks analytic binding"))?;
    let ParameterRange::Bounded{start,end}=binding.domain else {
        return Err(invalid("cannot split an unbounded physical edge"));
    };
    let curve=model.edge_curve(edge).ok_or_else(||invalid("missing source edge carrier"))?;
    let t=match curve {
        Curve3::Line(line)=>{
            let denom=line.direction.dot(line.direction);
            if denom<=0.0 || !denom.is_finite(){return Err(invalid("non-finite line carrier"));}
            line.direction.dot(p.sub(line.origin))/denom
        }
        Curve3::Circle(circle)=>{
            let v=p.sub(circle.frame.origin);
            let x=v.dot(circle.frame.x);let y=v.dot(circle.frame.y);
            let angle=y.atan2(x).rem_euclid(TAU);
            // Return the equivalent parameter inside the original bounded
            // domain, not an arbitrary angle measured from global origin.
            let candidate=[angle,angle+TAU,angle-TAU].into_iter()
                .find(|v|*v>=start-epsilon/circle.radius &&
                         *v<=end+epsilon/circle.radius)
                .ok_or_else(||invalid("circle parameter falls outside carrier range"))?;
            candidate
        }
    };
    let fraction=(t-start)/(end-start);
    if !t.is_finite() || !fraction.is_finite() || fraction<0.0 || fraction>1.0 {
        return Err(invalid("edge split parameter is outside its carrier"));
    }
    let q=curve.evaluate(t)?;
    if distance(q,p)>epsilon {return Err(invalid("edge split does not lie on original analytic carrier"));}
    if e.start!=e.end && (fraction==0.0 || fraction==1.0) {
        return Err(invalid("new split lies on an existing endpoint"));
    }
    Ok((t,fraction))
}

fn determine_placement(model:&BrepModel,face:FaceId,node:&GraphNode,
    operand:GraphOperand,epsilon:f64,boundaries:&[GraphBoundarySource])
    ->Result<PlanVertexPlacement,GeometryError>
{
    if boundaries.is_empty(){return Ok(PlanVertexPlacement::CreateFaceVertex{face});}
    let mut unique_edges=BTreeSet::new();
    let mut original_vertices=BTreeSet::new();
    for w in boundaries {
        if w.operand!=operand || w.face!=face {return Err(invalid("cross-operand witness in node placement"));}
        unique_edges.insert(w.edge);
        let edge=model.edges.get(w.edge.0 as usize).filter(|e|e.id==w.edge)
            .ok_or_else(||invalid("missing witnessed edge"))?;
        for vertex in [edge.start,edge.end] {
            let position=model.vertices.get(vertex.0 as usize)
                .ok_or_else(||invalid("witness references unknown source vertex"))?.position;
            if distance(position,node.point)<=epsilon {original_vertices.insert(vertex);}
        }
    }
    if original_vertices.len()>1 {
        return Err(invalid("multiple source vertices lie within proposed merge tolerance"));
    }
    if let Some(vertex)=original_vertices.into_iter().next(){
        return Ok(PlanVertexPlacement::ReuseVertex{vertex});
    }
    if unique_edges.len()!=1 {
        return Err(invalid("source boundary junction without a unique edge or existing vertex"));
    }
    let edge=*unique_edges.iter().next().expect("single source edge");
    let witness=boundaries.iter().find(|w|w.edge==edge).expect("edge from witness");
    let (parameter,fraction)=edge_parameter(model,edge,node.point,epsilon)?;
    Ok(PlanVertexPlacement::SplitEdge{edge,coedge:witness.coedge,parameter,fraction})
}
