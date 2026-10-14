//! Reconcile a REAL spatial intersection edge with its TWO independently
//! evaluable source-face UV pcurves, without publishing any topological IDs.
//! The halfedge directions are opposite in 3D. Their eventual face-loop
//! orientation/material use remains an explicit future Boolean decision.
use crate::{
    BrepModel, Curve2, Curve3, GeometryError,
    GeometryTolerance, GraphArcId, GraphFaceSource, GraphNodeId,
    GraphNodeKind, GraphState, IntersectionCurve,
    IntersectionGraph, IntersectionPcurve, ParameterRange,
    Point2Param, Point3, TrimContact,
};
fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn dist(a:Point3,b:Point3)->f64 {
    let d=Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z};
    d.x.hypot(d.y).hypot(d.z)
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum EdgeReconciliationHold {
    UnresolvedGraph,
    ContactNotTransverse,
    DuplicateOrOverlappingEdge,
    MissingPhysicalBoundaryWitness,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ReconciledCoedgeUse{
    pub source:GraphFaceSource,
    /// The common physical 3D edge runs in the natural parameter direction
    /// for the first coedge, and in the reverse direction for the second.
    pub reversed:bool,
    pub uv:IntersectionPcurve,
    pub curve2:Curve2,
}
impl ReconciledCoedgeUse{
    pub fn evaluate(&self,along_use:f64)->Result<Point2Param,GeometryError>{
        if !along_use.is_finite() || !(0.0..=1.0).contains(&along_use){
            return Err(invalid("reconciled coedge fraction must be inside [0,1]"));
        }
        self.uv.evaluate(if self.reversed{1.0-along_use}else{along_use})
    }
    pub fn point(&self,along_use:f64)->Result<Point3,GeometryError>{
        let uv=self.evaluate(along_use)?;
        Ok(self.uv.surface.evaluate(uv.u,uv.v)?.point)
    }
}
#[derive(Debug,Clone,PartialEq)]
pub struct ReconciledPhysicalEdge{
    pub source_arc:GraphArcId,
    pub origin_node:GraphNodeId,
    pub destination_node:GraphNodeId,
    pub curve:Curve3,
    pub parameter_domain:ParameterRange,
    pub start:Point3,
    pub end:Point3,
    pub uses:[ReconciledCoedgeUse;2],
    /// If true, the common endpoint is only a parameter anchor and MUST NOT
    /// become a separate material vertex merely due to a periodic chart.
    pub virtual_periodic_anchor:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub struct CoedgeReconciliationPlan{
    pub source_faces:[GraphFaceSource;2],
    pub edges:Vec<ReconciledPhysicalEdge>,
    /// An independent arrangement and topological commit are still required.
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum EdgeReconciliationResult{
    Ready(CoedgeReconciliationPlan),
    Held(EdgeReconciliationHold),
}
impl CoedgeReconciliationPlan{
    pub fn build(first:&BrepModel,second:&BrepModel,graph:&IntersectionGraph,
        tolerance:GeometryTolerance)->Result<EdgeReconciliationResult,GeometryError>
    {
        tolerance.validate()?;
        graph.validate_against(
            first,graph.source_faces[0].face,
            second,graph.source_faces[1].face)?;
        if graph.state!=GraphState::Ready || graph.arcs.is_empty() {
            return Ok(EdgeReconciliationResult::Held(
                EdgeReconciliationHold::UnresolvedGraph));
        }
        let mut edges:Vec<ReconciledPhysicalEdge>=Vec::new();
        for arc in &graph.arcs {
            if arc.contact!=TrimContact::Crossing {
                return Ok(EdgeReconciliationResult::Held(
                    EdgeReconciliationHold::ContactNotTransverse));
            }
            let from=graph.nodes.get(arc.from.0 as usize)
                .ok_or_else(||invalid("reconciled source vertex absent"))?;
            let to=graph.nodes.get(arc.to.0 as usize)
                .ok_or_else(||invalid("reconciled destination vertex absent"))?;
            let virtual_periodic_anchor=arc.closed &&
                arc.from==arc.to && from.kind==GraphNodeKind::PeriodicAnchor;
            if !arc.closed && from.kind==GraphNodeKind::PeriodicAnchor {
                return Ok(EdgeReconciliationResult::Held(
                    EdgeReconciliationHold::MissingPhysicalBoundaryWitness));
            }
            if !arc.closed && to.kind==GraphNodeKind::PeriodicAnchor {
                return Ok(EdgeReconciliationResult::Held(
                    EdgeReconciliationHold::MissingPhysicalBoundaryWitness));
            }
            if arc.closed && !virtual_periodic_anchor {
                return Ok(EdgeReconciliationResult::Held(
                    EdgeReconciliationHold::UnresolvedGraph));
            }
            if !virtual_periodic_anchor &&
                (from.boundary_sources.is_empty()||to.boundary_sources.is_empty()){
                return Ok(EdgeReconciliationResult::Held(
                    EdgeReconciliationHold::MissingPhysicalBoundaryWitness));
            }
            let start=arc.evaluate(0.0)?;let end=arc.evaluate(1.0)?;
            let eps=arc.tolerance.evaluation_budget.max(
                tolerance.length_at(dist(start,end).max(1.0)));
            if dist(start,from.point)>eps || dist(end,to.point)>eps{
                return Err(invalid("graph edge endpoints are not tied to source nodes"));
            }
            let [p0,p1]=arc.pcurves(first,second)?;
            let curve=match arc.curve{
                IntersectionCurve::Line(l)=>Curve3::Line(l),
                IntersectionCurve::PlaneCylinderSection(s)=>Curve3::Harmonic(
                    s.harmonic_curve()?),
            };
            let edge=ReconciledPhysicalEdge{
                source_arc:arc.id,
                origin_node:arc.from,destination_node:arc.to,
                curve,parameter_domain:ParameterRange::bounded(
                    arc.parameter_start,arc.parameter_end)?,
                start,end,virtual_periodic_anchor,
                uses:[
                    ReconciledCoedgeUse{
                        source:arc.source_faces[0],reversed:false,
                        curve2:p0.to_curve2()?,uv:p0,
                    },
                    ReconciledCoedgeUse{
                        source:arc.source_faces[1],reversed:true,
                        curve2:p1.to_curve2()?,uv:p1,
                    },
                ],
            };
            for f in [0.0,0.08,0.25,0.5,0.75,0.92,1.0]{
                let expected=arc.evaluate(f)?;
                let edge_t=arc.parameter_start+(arc.parameter_end-arc.parameter_start)*f;
                let from_curve=curve.evaluate(edge_t)?;
                if dist(from_curve,expected)>eps {
                    return Err(invalid("3D physical edge carrier differs from graph arc"));
                }
                if dist(edge.uses[0].point(f)?,expected)>eps ||
                   dist(edge.uses[1].point(1.0-f)?,expected)>eps {
                    return Err(invalid("two source pcurves disagree on the shared 3D edge"));
                }
            }
            let midpoint=curve.evaluate(0.5*(arc.parameter_start+arc.parameter_end))?;
            for old in &edges{
                if dist(old.start,start)<=eps && dist(old.end,end)<=eps {
                    let old_t=match old.parameter_domain{
                        ParameterRange::Bounded{start,end}=>(start+end)*0.5,
                        ParameterRange::Unbounded=>return Err(invalid(
                            "reconciled physical edge has unbounded source domain")),
                    };
                    if dist(old.curve.evaluate(old_t)?,midpoint)<=eps {
                        return Ok(EdgeReconciliationResult::Held(
                            EdgeReconciliationHold::DuplicateOrOverlappingEdge));
                    }
                }
            }
            edges.push(edge);
        }
        Ok(EdgeReconciliationResult::Ready(Self{
            source_faces:graph.source_faces,
            edges,can_apply:false,
        }))
    }
    pub fn validate_against(&self,first:&BrepModel,second:&BrepModel,
        graph:&IntersectionGraph,tolerance:GeometryTolerance)
        ->Result<(),GeometryError>
    {
        let computed=Self::build(first,second,graph,tolerance)?;
        if computed!=EdgeReconciliationResult::Ready(self.clone()){
            return Err(invalid("reconciled coedges no longer match source graph"));
        }
        Ok(())
    }
}
