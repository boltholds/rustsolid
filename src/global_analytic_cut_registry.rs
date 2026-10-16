//! G4b2b2: source-bound, full face-pair analytic cut carrier registry.
//!
//! This is the bridge from G1 globally reconciled physical endpoints to two
//! independently evaluable G2 source pcurves for every real crossing arc.
//! No proximity-based vertex merging, persistent result EdgeId or closed shell
//! is inferred. A transverse partial cylinder intersection is a valid
//! **preview**, never permission to commit a partial Boolean B-rep.

use crate::{
    BrepBody, Curve2, Curve3, FaceId, GeometryError, GeometryTolerance,
    GlobalArcRef, GlobalEventId, GlobalEventKind, GlobalNodeReconciliation,
    GlobalNodeReconciliationPlan, GlobalNodeRef, GraphNodeKind, GraphOperand,
    IntersectionCurve, MultiFaceGraphState, MultiFaceIntersectionGraph,
    MultiFaceIntersectionGraphBuilder, ParameterRange, Point3,
    TrimContact,
};

fn invalid(reason:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(reason.into())
}
fn dist(a:Point3,b:Point3)->f64{
    let d=a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}
fn revision(body:&BrepBody)->Option<u64>{
    match body{
        BrepBody::Polyhedral(s)=>Some(s.revision),
        BrepBody::ThroughHole(s)=>Some(s.revision()),
        BrepBody::BlindHole(s)=>Some(s.revision()),
        _=>None,
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticCutRegistryHold{
    UpstreamGraphAmbiguity,
    UnreconciledPhysicalVertex,
    ContactNotTransverse,
    NoCertifiedCrossing,
    UnsupportedAnalyticPcurve,
    TooManyArcs,
    InvalidSourceEvent,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct SourceCertifiedCutArc{
    /// Pair+local arc identity is never replaced by a flat ambiguous ID.
    pub source:GlobalArcRef,
    pub first_face:FaceId,
    pub second_face:FaceId,
    pub carrier:Curve3,
    pub parameter_domain:ParameterRange,
    /// These are two independent, source-face **Curve2** carriers,
    /// not generated coedges or synthetic mesh intersections.
    pub face_pcurves:[Curve2;2],
    /// Closed periodic traces retain no actual physical corner events:
    /// the graph's chart origin is virtual and must not be welded.
    pub physical_endpoints:Option<[GlobalEventId;2]>,
    pub virtual_periodic_anchor:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub struct GlobalAnalyticCutPlan {
    pub arcs:Vec<SourceCertifiedCutArc>,
    pub global_nodes:GlobalNodeReconciliationPlan,
    /// Source graph includes ALL A/B face pairs, including Empty pairs.
    pub audited_face_pairs:usize,
    /// Identifies only the source event space. NOT resulting B-rep vertices.
    pub physical_event_count:usize,
    /// The source graph may also contain non-cut contacts, never promoted.
    pub discarded_original_boundary_contacts:usize,
    pub source_revisions:[Option<u64>;2],
    pub construction_tolerance:GeometryTolerance,
    pub can_apply:bool,
    pub closed_shell_proven:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum GlobalAnalyticCutResult{
    Ready(GlobalAnalyticCutPlan),
    Held(AnalyticCutRegistryHold),
}

impl GlobalAnalyticCutPlan {
    pub fn validate_against(&self,first:&BrepBody,second:&BrepBody,
        tolerance:GeometryTolerance)->Result<(),GeometryError>
    {
        if self.source_revisions!=[revision(first),revision(second)]||
            self.construction_tolerance!=tolerance {
            return Err(invalid("G4b2b2 source revision or tolerance policy changed"));
        }
        let rebuilt=Self::build(first,second,tolerance)?;
        if rebuilt!=GlobalAnalyticCutResult::Ready(self.clone()){
            return Err(invalid("G4b2b2 global analytic cut registry changed"));
        }
        Ok(())
    }

    pub fn build(first:&BrepBody,second:&BrepBody,
        tolerance:GeometryTolerance)->Result<GlobalAnalyticCutResult,GeometryError>
    {
        tolerance.validate()?;first.validate()?;second.validate()?;
        let a=first.shared()?;let b=second.shared()?;
        let graph=MultiFaceIntersectionGraphBuilder::new(tolerance)?
            .build(&a,&b)?;
        // Includes full face coverage, checks every Empty state and
        // recomputes every certified source pair.
        graph.validate_against(&a,&b)?;
        if graph.state==MultiFaceGraphState::Held {
            return Ok(GlobalAnalyticCutResult::Held(
                AnalyticCutRegistryHold::UpstreamGraphAmbiguity));
        }
        let global=GlobalNodeReconciliation::build(&graph,&a,&b)?;
        global.validate_against(&graph,&a,&b)?;
        let GlobalNodeReconciliation::Ready(global)=global else{
            return Ok(GlobalAnalyticCutResult::Held(
                AnalyticCutRegistryHold::UnreconciledPhysicalVertex));
        };
        if graph.pairs.len()>512 || graph.pairs.iter()
            .map(|pair|pair.graph.arcs.len()).sum::<usize>()>2048 {
            return Ok(GlobalAnalyticCutResult::Held(
                AnalyticCutRegistryHold::TooManyArcs));
        }
        let mut arcs=Vec::new();
        let mut contacts=0usize;
        for pair in &graph.pairs {
            for arc in &pair.graph.arcs {
                match arc.contact {
                    TrimContact::BoundaryOverlap=>{
                        // Existing source boundary, not a transverse newly
                        // material-separating Boolean edge.
                        contacts+=1;
                        continue;
                    },
                    TrimContact::Tangent=>return Ok(GlobalAnalyticCutResult::Held(
                        AnalyticCutRegistryHold::ContactNotTransverse)),
                    TrimContact::Crossing=>{},
                }
                let endpoints=if arc.closed {
                    let Some(node)=pair.graph.nodes.get(arc.from.0 as usize)else{
                        return Ok(GlobalAnalyticCutResult::Held(
                            AnalyticCutRegistryHold::InvalidSourceEvent));
                    };
                    if arc.from!=arc.to||node.kind!=GraphNodeKind::PeriodicAnchor{
                        return Ok(GlobalAnalyticCutResult::Held(
                            AnalyticCutRegistryHold::InvalidSourceEvent));
                    }
                    None
                }else{
                    let resolve=|id:crate::GraphNodeId|->Option<GlobalEventId>{
                        let local=GlobalNodeRef{pair:pair.id,local:id};
                        let event=global.event_for(local)?;
                        if event.kind==GlobalEventKind::BoundaryEndpoint {
                            Some(event.id)
                        }else{None}
                    };
                    let (Some(start),Some(end))=(resolve(arc.from),resolve(arc.to))else{
                        return Ok(GlobalAnalyticCutResult::Held(
                            AnalyticCutRegistryHold::InvalidSourceEvent));
                    };
                    if start==end{return Ok(GlobalAnalyticCutResult::Held(
                        AnalyticCutRegistryHold::InvalidSourceEvent));}
                    Some([start,end])
                };
                let [first_uv,second_uv]=arc.pcurves(&a,&b)?;
                let c2=[first_uv.to_curve2()?,second_uv.to_curve2()?];
                let c3=match arc.curve {
                    IntersectionCurve::Line(line)=>Curve3::Line(line),
                    IntersectionCurve::PlaneCylinderSection(section)=>
                        Curve3::Harmonic(section.harmonic_curve()?),
                };
                let domain=ParameterRange::bounded(
                    arc.parameter_start,arc.parameter_end)?;
                let extent=a.bbox.max.sub(a.bbox.min).length_squared().sqrt()
                    .max(b.bbox.max.sub(b.bbox.min).length_squared().sqrt());
                let eps=tolerance.length_at(extent.max(1.0))*32.0+
                    arc.tolerance.evaluation_budget;
                let Some(surface_a)=a.face_surface(pair.first_face) else{
                    return Err(invalid("G4b2b2 missing source A surface"));
                };
                let Some(surface_b)=b.face_surface(pair.second_face) else{
                    return Err(invalid("G4b2b2 missing source B surface"));
                };
                for f in [0.0_f64,0.125,0.25,0.5,0.75,0.875,1.0] {
                    let t=arc.parameter_start+
                        (arc.parameter_end-arc.parameter_start)*f;
                    let xyz=arc.evaluate(f)?;
                    if dist(c3.evaluate(t)?,xyz)>eps ||
                        dist(surface_a.evaluate(
                            c2[0].evaluate(t)?.u,c2[0].evaluate(t)?.v
                        )?.point,xyz)>eps ||
                        dist(surface_b.evaluate(
                            c2[1].evaluate(t)?.u,c2[1].evaluate(t)?.v
                        )?.point,xyz)>eps
                    {
                        return Ok(GlobalAnalyticCutResult::Held(
                            AnalyticCutRegistryHold::UnsupportedAnalyticPcurve));
                    }
                }
                arcs.push(SourceCertifiedCutArc{
                    source:GlobalArcRef{pair:pair.id,local:arc.id},
                    first_face:pair.first_face,second_face:pair.second_face,
                    carrier:c3,parameter_domain:domain,
                    face_pcurves:c2,physical_endpoints:endpoints,
                    virtual_periodic_anchor:arc.closed,
                });
            }
        }
        if arcs.is_empty(){
            return Ok(GlobalAnalyticCutResult::Held(
                AnalyticCutRegistryHold::NoCertifiedCrossing));
        }
        Ok(GlobalAnalyticCutResult::Ready(Self{
            arcs,physical_event_count:global.events.iter()
                .filter(|event|event.kind==GlobalEventKind::BoundaryEndpoint)
                .count(),
            global_nodes:global,audited_face_pairs:graph.pairs.len(),
            discarded_original_boundary_contacts:contacts,
            source_revisions:[revision(first),revision(second)],
            construction_tolerance:tolerance,
            can_apply:false,closed_shell_proven:false,
        }))
    }
}
