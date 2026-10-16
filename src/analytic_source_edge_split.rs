//! G4b2b3a: genuine, source-certified physical splits of original B-rep edges.
//!
//! A new G1 boundary event becomes ONE vertex per source operand, using its
//! original source-edge witness, never geometric-nearness vertex welding.
//! Edges and BOTH original coedges are partitioned with their actual analytic
//! Curve3/Curve2 parameter subdomains, then the source B-rep is revalidated.
//! The two resulting models are subdivided SOURCE bodies, NOT a Boolean shell.
//! Curved original edge splitting, material-region selection, curved patch
//! assembly, independent volume certification and G5 remain Held.

use std::collections::{BTreeMap, BTreeSet};
use crate::{
    AnalyticCutRegistryHold, BrepBody, BrepModel, Coedge, CoedgeGeometry,
    CoedgeId, Edge, EdgeGeometry, EdgeId, GeometryError, GeometryTolerance,
    GlobalAnalyticCutPlan, GlobalAnalyticCutResult, GlobalEventId,
    GlobalEventKind, GraphOperand, MultiFaceIntersectionGraphBuilder,
    ParameterRange, Point3, Vertex, VertexId,
};

const MAX_NEW_VERTICES: usize = 2048;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn distance(a: Point3, b: Point3) -> f64 {
    let d = a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}
fn bounded(domain: ParameterRange) -> Result<[f64; 2], GeometryError> {
    match domain {
        ParameterRange::Bounded { start, end } => Ok([start, end]),
        _ => Err(invalid("G4b2b3a requires a bounded source edge/coedge")),
    }
}
fn interval(parent: [f64; 2], lo: f64, hi: f64)
    -> Result<ParameterRange, GeometryError>
{
    ParameterRange::bounded(
        parent[0] + (parent[1] - parent[0])*lo,
        parent[0] + (parent[1] - parent[0])*hi,
    )
}
fn oriented_interval(parent: [f64; 2], lo: f64, hi: f64, reversed: bool)
    -> Result<ParameterRange, GeometryError>
{
    if reversed { interval(parent, 1.0-hi, 1.0-lo) }
    else { interval(parent, lo, hi) }
}
fn original_vertex(model: &BrepModel, id: VertexId) -> Result<Point3, GeometryError> {
    model.vertices.get(id.0 as usize)
        .filter(|v| v.id == id).map(|v| v.position)
        .ok_or_else(|| invalid("G4b2b3a missing original source vertex"))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalyticSourceEdgeSplitHold {
    UpstreamRegistry(AnalyticCutRegistryHold),
    MissingSourceEdgeWitness,
    UnsupportedOriginalEdgeCurve,
    UncertainSourceEvent,
    DuplicatePhysicalEvent,
    EventIdentityConflict,
    MissingInteriorSplit,
    SplitBudgetExceeded,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SourceEdgePiece {
    /// Original EdgeId is kept as lineage, never confused across operands.
    pub original_edge: EdgeId,
    pub result_edge: EdgeId,
    pub result_coedges: [CoedgeId; 2],
    pub vertices: [VertexId; 2],
    /// Increasing fractions along the ORIGINAL physical source edge.
    pub source_fractions: [f64; 2],
    pub curve_domain: ParameterRange,
    /// Each original coedge's own UV parameter domain, in directed-use order.
    pub coedge_domains: [ParameterRange; 2],
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnalyticSourceEdgeSplit {
    pub operand: GraphOperand,
    pub original_edge: EdgeId,
    pub original_coedges: [CoedgeId; 2],
    /// Only G1 events STRICTLY interior to this physical edge create vertices.
    pub interior_events: Vec<GlobalEventId>,
    pub pieces: Vec<SourceEdgePiece>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnalyticSplitEventBinding {
    pub operand: GraphOperand,
    pub event: GlobalEventId,
    pub source_edges: Vec<EdgeId>,
    /// One result vertex even when many face pairs witnessed this G1 event.
    pub result_vertex: VertexId,
    pub created_vertex: bool,
}

#[derive(Debug, Clone)]
pub struct GlobalAnalyticEdgeSplitPlan {
    pub registry: GlobalAnalyticCutPlan,
    /// Fully validated source B-reps with real vertex, edge and coedge splits.
    /// These are NOT oriented retained patches of the Boolean operation.
    pub split_sources: [BrepModel; 2],
    pub edges: Vec<AnalyticSourceEdgeSplit>,
    pub event_bindings: Vec<AnalyticSplitEventBinding>,
    pub new_source_vertices: [usize; 2],
    pub can_apply: bool,
    pub boolean_shell_proven: bool,
}
#[derive(Debug, Clone)]
pub enum GlobalAnalyticEdgeSplitResult {
    Ready(GlobalAnalyticEdgeSplitPlan),
    Held(AnalyticSourceEdgeSplitHold),
}

#[derive(Debug, Clone, Copy)]
struct EdgeKnot {
    fraction: f64,
    event: GlobalEventId,
    snapped_to_edge: Point3,
}
#[derive(Debug, Clone)]
struct EventAssignment {
    vertex: VertexId,
    created: bool,
    source_edges: BTreeSet<EdgeId>,
}

impl GlobalAnalyticEdgeSplitPlan {
    pub fn validate_against(
        &self, first: &BrepBody, second: &BrepBody,
        tolerance: GeometryTolerance,
    ) -> Result<(), GeometryError> {
        self.registry.validate_against(first, second, tolerance)?;
        let GlobalAnalyticEdgeSplitResult::Ready(rebuilt) =
            Self::build(first, second, tolerance)? else {
            return Err(invalid("G4b2b3a source edge split is no longer supported"));
        };
        if self.edges != rebuilt.edges
            || self.event_bindings != rebuilt.event_bindings
            || self.new_source_vertices != rebuilt.new_source_vertices
            || self.registry != rebuilt.registry
            || self.can_apply != rebuilt.can_apply
            || self.boolean_shell_proven != rebuilt.boolean_shell_proven
            || format!("{:?}", self.split_sources) != format!("{:?}", rebuilt.split_sources)
        {
            return Err(invalid("G4b2b3a physical source-edge split differs from sources"));
        }
        self.split_sources[0].validate()?;
        self.split_sources[1].validate()?;
        Ok(())
    }

    pub fn build(
        first: &BrepBody, second: &BrepBody,
        tolerance: GeometryTolerance,
    ) -> Result<GlobalAnalyticEdgeSplitResult, GeometryError> {
        tolerance.validate()?;
        let registry = match GlobalAnalyticCutPlan::build(first, second, tolerance)? {
            GlobalAnalyticCutResult::Ready(plan) => plan,
            GlobalAnalyticCutResult::Held(reason) =>
                return Ok(GlobalAnalyticEdgeSplitResult::Held(
                    AnalyticSourceEdgeSplitHold::UpstreamRegistry(reason))),
        };
        let originals = [first.shared()?, second.shared()?];
        let graph = MultiFaceIntersectionGraphBuilder::new(tolerance)?
            .build(&originals[0], &originals[1])?;
        graph.validate_against(&originals[0], &originals[1])?;
        if registry.audited_face_pairs != graph.pairs.len() {
            return Err(invalid("G4b2b3a source graph and certified registry disagree"));
        }
        let used_events: BTreeSet<GlobalEventId> = registry.arcs.iter()
            .filter_map(|arc| arc.physical_endpoints)
            .flat_map(|endpoints| endpoints.into_iter()).collect();
        let mut witnesses: BTreeMap<(GraphOperand, EdgeId), BTreeMap<GlobalEventId,Point3>> =
            BTreeMap::new();
        for event in &registry.global_nodes.events {
            if !used_events.contains(&event.id) { continue; }
            if event.kind != GlobalEventKind::BoundaryEndpoint {
                return Ok(GlobalAnalyticEdgeSplitResult::Held(
                    AnalyticSourceEdgeSplitHold::MissingSourceEdgeWitness));
            }
            let mut has_source = false;
            for member in &event.members {
                let pair = graph.pairs.iter().find(|pair| pair.id == member.pair)
                    .ok_or_else(||invalid("G4b2b3a missing source face-pair graph"))?;
                let node = pair.graph.nodes.get(member.local.0 as usize)
                    .filter(|node|node.id == member.local)
                    .ok_or_else(||invalid("G4b2b3a missing source graph node"))?;
                for witness in &node.boundary_sources {
                    let source = match witness.operand {
                        GraphOperand::First => &originals[0],
                        GraphOperand::Second => &originals[1],
                    };
                    let edge = source.edges.get(witness.edge.0 as usize)
                        .filter(|edge| edge.id == witness.edge)
                        .ok_or_else(||invalid("G4b2b3a unknown witness EdgeId"))?;
                    let coedge = source.coedges.get(witness.coedge.0 as usize)
                        .filter(|co|co.id == witness.coedge && co.face == witness.face
                            && co.edge == edge.id)
                        .ok_or_else(||invalid("G4b2b3a invalid original coedge witness"))?;
                    if !edge.coedges.contains(&coedge.id) {
                        return Err(invalid("G4b2b3a witness not incident to physical source edge"));
                    }
                    witnesses.entry((witness.operand, edge.id))
                        .or_default().insert(event.id, event.point);
                    has_source = true;
                }
            }
            if !has_source {
                return Ok(GlobalAnalyticEdgeSplitResult::Held(
                    AnalyticSourceEdgeSplitHold::MissingSourceEdgeWitness));
            }
        }

        let mut split_sources = originals.clone();
        let mut edges = Vec::new();
        let mut event_assignments: BTreeMap<(GraphOperand,GlobalEventId),EventAssignment> =
            BTreeMap::new();
        let mut new_source_vertices = [0usize;2];
        for (operand_index, operand) in
            [(0usize,GraphOperand::First),(1usize,GraphOperand::Second)]
        {
            let source = &originals[operand_index];
            let model = &mut split_sources[operand_index];
            let extent = source.bbox.max.sub(source.bbox.min)
                .length_squared().sqrt().max(1.0);
            let coordinate_scale = source.vertices.iter()
                .flat_map(|v| [v.position.x.abs(),v.position.y.abs(),v.position.z.abs()])
                .fold(1.0_f64,f64::max);
            let eps = tolerance.length_at(extent)*16.0 +
                coordinate_scale*32.0*f64::EPSILON;
            for ((owner, edge_id), edge_events) in witnesses.iter() {
                if *owner != operand { continue; }
                let original_edge = source.edges.get(edge_id.0 as usize)
                    .filter(|e|e.id == *edge_id)
                    .ok_or_else(||invalid("G4b2b3a source EdgeId vanished"))?;
                // A period-closed source circle needs angular split policy;
                // never infer an arbitrary periodic vertex from XYZ proximity.
                if original_edge.start == original_edge.end {
                    return Ok(GlobalAnalyticEdgeSplitResult::Held(
                        AnalyticSourceEdgeSplitHold::UnsupportedOriginalEdgeCurve));
                }
                let Some(crate::Curve3::Line(line)) = source.edge_curve(*edge_id)
                else {
                    return Ok(GlobalAnalyticEdgeSplitResult::Held(
                        AnalyticSourceEdgeSplitHold::UnsupportedOriginalEdgeCurve));
                };
                let a = original_vertex(source, original_edge.start)?;
                let b = original_vertex(source, original_edge.end)?;
                let delta = b.sub(a);
                let length_squared = delta.length_squared();
                if length_squared <= eps*eps {
                    return Ok(GlobalAnalyticEdgeSplitResult::Held(
                        AnalyticSourceEdgeSplitHold::UncertainSourceEvent));
                }
                let length = length_squared.sqrt();
                let edge_geom = *source.geometry.edges.get(edge_id.0 as usize)
                    .ok_or_else(||invalid("G4b2b3a missing parent edge domain"))?;
                let edge_domain = bounded(edge_geom.domain)?;
                if distance(line.evaluate(edge_domain[0])?,a)>eps
                    || distance(line.evaluate(edge_domain[1])?,b)>eps
                {
                    return Err(invalid("G4b2b3a source edge parameterization not aligned"));
                }
                let mut knots = Vec::<EdgeKnot>::new();
                for (&event,&point) in edge_events {
                    let projected_fraction = point.sub(a).dot(delta)/length_squared;
                    let projected = a.add(delta.scale(projected_fraction));
                    if !projected_fraction.is_finite()
                        || projected_fraction < -eps/length
                        || projected_fraction > 1.0+eps/length
                        || distance(projected,point)>eps
                    {
                        return Ok(GlobalAnalyticEdgeSplitResult::Held(
                            AnalyticSourceEdgeSplitHold::UncertainSourceEvent));
                    }
                    let fraction = projected_fraction.clamp(0.0,1.0);
                    let on_curve = line.evaluate(
                        edge_domain[0]+(edge_domain[1]-edge_domain[0])*fraction)?;
                    if distance(on_curve,point)>eps {
                        return Ok(GlobalAnalyticEdgeSplitResult::Held(
                            AnalyticSourceEdgeSplitHold::UncertainSourceEvent));
                    }
                    knots.push(EdgeKnot {fraction,event,snapped_to_edge:on_curve});
                }
                knots.sort_by(|a,b|a.fraction.total_cmp(&b.fraction)
                    .then(a.event.cmp(&b.event)));
                for window in knots.windows(2) {
                    if (window[1].fraction-window[0].fraction)*length<=eps*2.0 {
                        return Ok(GlobalAnalyticEdgeSplitResult::Held(
                            AnalyticSourceEdgeSplitHold::DuplicatePhysicalEvent));
                    }
                }
                let mut interior = Vec::<EdgeKnot>::new();
                for knot in knots {
                    let at_start = knot.fraction*length<=eps;
                    let at_end = (1.0-knot.fraction)*length<=eps;
                    let existing = if at_start { Some(original_edge.start) }
                        else if at_end { Some(original_edge.end) } else { None };
                    let assignment_key = (operand,knot.event);
                    let vertex = if let Some(id) = existing {
                        id
                    } else {
                        if new_source_vertices[operand_index]>=MAX_NEW_VERTICES {
                            return Ok(GlobalAnalyticEdgeSplitResult::Held(
                                AnalyticSourceEdgeSplitHold::SplitBudgetExceeded));
                        }
                        let id=VertexId(u32::try_from(model.vertices.len())
                            .map_err(|_|invalid("G4b2b3a vertex ID overflow"))?);
                        model.vertices.push(Vertex {
                            id,position:knot.snapped_to_edge,
                        });
                        new_source_vertices[operand_index]+=1;
                        interior.push(knot);
                        id
                    };
                    use std::collections::btree_map::Entry;
                    match event_assignments.entry(assignment_key) {
                        Entry::Occupied(mut entry) => {
                            if entry.get().vertex!=vertex {
                                return Ok(GlobalAnalyticEdgeSplitResult::Held(
                                    AnalyticSourceEdgeSplitHold::EventIdentityConflict));
                            }
                            entry.get_mut().source_edges.insert(*edge_id);
                        },
                        Entry::Vacant(entry) => {
                            let mut sources=BTreeSet::new(); sources.insert(*edge_id);
                            entry.insert(EventAssignment {
                                vertex,created:existing.is_none(),source_edges:sources,
                            });
                        },
                    }
                }
                if interior.is_empty() { continue; }
                let mut fractions = vec![0.0_f64];
                let mut vertices = vec![original_edge.start];
                let mut interior_events = Vec::new();
                for knot in &interior {
                    fractions.push(knot.fraction);
                    let vertex = event_assignments[&(operand,knot.event)].vertex;
                    vertices.push(vertex);
                    interior_events.push(knot.event);
                }
                fractions.push(1.0);
                vertices.push(original_edge.end);
                let n = fractions.len()-1;
                let old_coedges = original_edge.coedges;
                let source_coedges: [Coedge; 2] = [
                    source.coedges[old_coedges[0].0 as usize].clone(),
                    source.coedges[old_coedges[1].0 as usize].clone(),
                ];
                let source_coedge_geometry: [CoedgeGeometry;2] = [
                    source.geometry.coedges[old_coedges[0].0 as usize],
                    source.geometry.coedges[old_coedges[1].0 as usize],
                ];
                // Keep the original edge ID for its first physical piece.
                // All new source edges are appended; original names survive.
                let mut result_edge_ids = vec![original_edge.id];
                for _ in 1..n {
                    let id=EdgeId(u32::try_from(model.edges.len())
                        .map_err(|_|invalid("G4b2b3a edge ID overflow"))?);
                    model.edges.push(original_edge.clone());
                    model.geometry.edges.push(edge_geom);
                    result_edge_ids.push(id);
                }
                let mut uses = [Vec::<CoedgeId>::new(),Vec::<CoedgeId>::new()];
                for side in 0..2 {
                    let directed_first = if source_coedges[side].reversed {n-1}else{0};
                    for segment in 0..n {
                        let id = if segment==directed_first {
                            source_coedges[side].id
                        } else {
                            let id=CoedgeId(u32::try_from(model.coedges.len())
                                .map_err(|_|invalid("G4b2b3a coedge ID overflow"))?);
                            model.coedges.push(source_coedges[side].clone());
                            model.geometry.coedges.push(source_coedge_geometry[side]);
                            id
                        };
                        uses[side].push(id);
                    }
                }
                let mut pieces = Vec::new();
                for segment in 0..n {
                    let result_edge = result_edge_ids[segment];
                    let coedge_ids = [uses[0][segment],uses[1][segment]];
                    let domain = interval(edge_domain,fractions[segment],fractions[segment+1])?;
                    model.edges[result_edge.0 as usize] = Edge {
                        id:result_edge,start:vertices[segment],
                        end:vertices[segment+1],coedges:coedge_ids,
                    };
                    model.geometry.edges[result_edge.0 as usize] = EdgeGeometry {
                        domain,..edge_geom
                    };
                    let mut uv_domains = [ParameterRange::Unbounded;2];
                    for side in 0..2 {
                        let coedge = &source_coedges[side];
                        let id=coedge_ids[side];
                        let uv_domain=oriented_interval(
                            bounded(source_coedge_geometry[side].domain)?,
                            fractions[segment],fractions[segment+1],
                            coedge.reversed)?;
                        uv_domains[side]=uv_domain;
                        // A previous source-edge split may already have
                        // spliced a neighbouring loop use. Preserve its
                        // CURRENT links when reusing an original coedge ID.
                        let current=model.coedges[id.0 as usize].clone();
                        model.coedges[id.0 as usize] = Coedge {
                            id,edge:result_edge,face:coedge.face,
                            loop_id:coedge.loop_id,reversed:coedge.reversed,
                            next:current.next,prev:current.prev,
                            twin:coedge_ids[1-side],
                        };
                        model.geometry.coedges[id.0 as usize] = CoedgeGeometry {
                            domain:uv_domain,..source_coedge_geometry[side]
                        };
                    }
                    pieces.push(SourceEdgePiece {
                        original_edge:*edge_id,result_edge,result_coedges:coedge_ids,
                        vertices:[vertices[segment],vertices[segment+1]],
                        source_fractions:[fractions[segment],fractions[segment+1]],
                        curve_domain:domain,coedge_domains:uv_domains,
                    });
                }
                // Replace each original coedge in its ORIGINAL directed loop
                // by a chain, respecting its sense and periodic chart use.
                for side in 0..2 {
                    let old = &source_coedges[side];
                    let chain:Vec<CoedgeId> = if old.reversed {
                        uses[side].iter().copied().rev().collect()
                    } else {
                        uses[side].clone()
                    };
                    let prev=model.coedges[old.id.0 as usize].prev;
                    let next=model.coedges[old.id.0 as usize].next;
                    model.coedges[prev.0 as usize].next=chain[0];
                    model.coedges[next.0 as usize].prev=*chain.last().unwrap();
                    for (index,&id) in chain.iter().enumerate() {
                        model.coedges[id.0 as usize].prev =
                            if index==0 {prev}else{chain[index-1]};
                        model.coedges[id.0 as usize].next =
                            if index+1==chain.len() {next}else{chain[index+1]};
                    }
                }
                edges.push(AnalyticSourceEdgeSplit {
                    operand,original_edge:*edge_id,original_coedges:old_coedges,
                    interior_events,pieces,
                });
            }
            model.validate()?;
            if model.euler_characteristic()!=source.euler_characteristic()
                || model.faces.len()!=source.faces.len()
                || model.loops.len()!=source.loops.len()
                || model.mass!=source.mass
            {
                return Err(invalid("G4b2b3a source subdivision changed shell invariants"));
            }
        }
        if new_source_vertices==[0,0] {
            return Ok(GlobalAnalyticEdgeSplitResult::Held(
                AnalyticSourceEdgeSplitHold::MissingInteriorSplit));
        }
        let event_bindings = event_assignments.into_iter().map(
            |((operand,event),assignment)| AnalyticSplitEventBinding {
                operand,event,source_edges:assignment.source_edges.into_iter().collect(),
                result_vertex:assignment.vertex,created_vertex:assignment.created,
            }).collect();
        Ok(GlobalAnalyticEdgeSplitResult::Ready(Self {
            registry,split_sources,edges,event_bindings,new_source_vertices,
            can_apply:false,boolean_shell_proven:false,
        }))
    }
}
