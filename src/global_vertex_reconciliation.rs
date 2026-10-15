//! Read-only G1 event reconciliation across independently certified face-pair graphs.
//!
//! Source graphs keep their original operand/face/edge/coedge identities.
//! Endpoint events are combined only when adjacent source-face topology supplies
//! a common physical edge or vertex on EACH changed operand. Pure proximity is
//! never proof of identity. Open curve crossings at intermediate parameters are
//! not split here; that belongs to global UV arrangement (G2).
//!
//! No persistent B-rep vertices, face splits, or Boolean edits are produced.

use crate::{
    BrepModel, EdgeId, FaceId, GeometryError, GlobalNodeRef, GraphNodeKind,
    GraphOperand, MultiFaceGraphState, MultiFaceIntersectionGraph, Point3, VertexId,
};
use std::collections::{BTreeMap, BTreeSet};

const MAX_PAIRS: usize = 512;
const MAX_EVENTS: usize = 2048;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn distance(a: Point3, b: Point3) -> f64 {
    let delta = a.sub(b);
    delta.x.hypot(delta.y).hypot(delta.z)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GlobalEventId(pub u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlobalEventKind {
    /// An endpoint of one or more analytically intersecting, open branches.
    BoundaryEndpoint,
    /// A point-only contact. It must not silently become a cut endpoint.
    IsolatedContact,
    /// Merely the parameter origin of a closed branch, not a physical vertex.
    VirtualPeriodicAnchor,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReconciledGlobalEvent {
    pub id: GlobalEventId,
    /// Representative from the first local event, not a snapped B-rep vertex.
    pub point: Point3,
    pub kind: GlobalEventKind,
    /// Ordered by face pair and original local node identity.
    pub members: Vec<GlobalNodeRef>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GlobalNodeReconciliationPlan {
    pub events: Vec<ReconciledGlobalEvent>,
    pub source_local_nodes: usize,
    /// G1 is read-only; a topology mutation still requires G2--G5.
    pub can_apply: bool,
}

impl GlobalNodeReconciliationPlan {
    /// Stable lookup from a certified local graph event to its G1 cluster.
    /// The returned ID is a preview identity, never a persistent B-rep ID.
    pub fn event_for(&self, local: GlobalNodeRef) -> Option<&ReconciledGlobalEvent> {
        self.events.iter().find(|event| event.members.contains(&local))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlobalNodeReconciliationHold {
    UpstreamAmbiguity,
    BudgetExceeded,
    /// Two events are closer than the resolution policy, but do not have
    /// source-topology evidence that they are the same physical point.
    MissingSharedTopologyWitness,
    /// Different face-pair tolerances leave identity numerically undecidable.
    UnresolvedNearCoincidence,
    /// A transitive group is not a mutually consistent geometric clique.
    NonCliqueOrDuplicatePair,
}

#[derive(Debug, Clone, PartialEq)]
pub enum GlobalNodeReconciliation {
    Ready(GlobalNodeReconciliationPlan),
    Held(GlobalNodeReconciliationHold),
}

#[derive(Debug, Clone)]
struct PreparedEvent {
    reference: GlobalNodeRef,
    point: Point3,
    kind: GraphNodeKind,
    faces: [FaceId; 2],
    resolution: f64,
    edges: BTreeSet<(GraphOperand, EdgeId)>,
    vertices: BTreeSet<(GraphOperand, VertexId)>,
}

impl PreparedEvent {
    fn common_feature(&self, other: &Self, operand: GraphOperand) -> bool {
        self.edges.iter().any(|feature| {
            feature.0 == operand && other.edges.contains(feature)
        }) || self.vertices.iter().any(|feature| {
            feature.0 == operand && other.vertices.contains(feature)
        })
    }

    /// Equal source faces are automatically compatible; different source
    /// faces require proof of adjacency through an identical real B-rep
    /// edge/vertex on that operand. A periodic chart origin cannot be welded.
    fn weldable(&self, other: &Self) -> bool {
        if self.reference.pair == other.reference.pair ||
            self.kind != GraphNodeKind::Endpoint ||
            other.kind != GraphNodeKind::Endpoint {
            return false;
        }
        [(GraphOperand::First, 0), (GraphOperand::Second, 1)]
            .iter()
            .all(|(operand, index)| {
                self.faces[*index] == other.faces[*index] ||
                    self.common_feature(other, *operand)
            })
    }
}

fn root(parents: &[usize], mut index: usize) -> usize {
    while parents[index] != index {
        index = parents[index];
    }
    index
}

fn reconcile_prepared(events: &[PreparedEvent]) -> GlobalNodeReconciliation {
    if events.len() > MAX_EVENTS {
        return GlobalNodeReconciliation::Held(
            GlobalNodeReconciliationHold::BudgetExceeded,
        );
    }
    let mut parents: Vec<usize> = (0..events.len()).collect();
    for i in 0..events.len() {
        for j in i + 1..events.len() {
            if events[i].reference.pair == events[j].reference.pair {
                continue;
            }
            let d = distance(events[i].point, events[j].point);
            // Both policies must resolve the point as one event. The wider
            // envelope is treated as uncertain rather than silently merged.
            let close = events[i].resolution.min(events[j].resolution);
            let uncertain = events[i].resolution.max(events[j].resolution);
            if d > uncertain {
                continue;
            }
            if d > close {
                return GlobalNodeReconciliation::Held(
                    GlobalNodeReconciliationHold::UnresolvedNearCoincidence,
                );
            }
            if !events[i].weldable(&events[j]) {
                return GlobalNodeReconciliation::Held(
                    GlobalNodeReconciliationHold::MissingSharedTopologyWitness,
                );
            }
            let a = root(&parents, i);
            let b = root(&parents, j);
            if a != b {
                let smaller = a.min(b);
                parents[a] = smaller;
                parents[b] = smaller;
            }
        }
    }

    let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..events.len() {
        clusters.entry(root(&parents, i)).or_default().push(i);
    }
    let mut output = Vec::with_capacity(clusters.len());
    for members in clusters.into_values() {
        let mut seen_pairs = BTreeSet::new();
        // Pairwise checks disallow tolerance-chain welding (A~B and B~C
        // but A not~C) and two distinct events from one source pair.
        for (i, &left) in members.iter().enumerate() {
            if !seen_pairs.insert(events[left].reference.pair) {
                return GlobalNodeReconciliation::Held(
                    GlobalNodeReconciliationHold::NonCliqueOrDuplicatePair,
                );
            }
            for &right in members.iter().skip(i + 1) {
                if distance(events[left].point, events[right].point)
                    > events[left].resolution.min(events[right].resolution)
                    || !events[left].weldable(&events[right])
                {
                    return GlobalNodeReconciliation::Held(
                        GlobalNodeReconciliationHold::NonCliqueOrDuplicatePair,
                    );
                }
            }
        }
        let representative = &events[members[0]];
        let kind = match representative.kind {
            GraphNodeKind::Endpoint => GlobalEventKind::BoundaryEndpoint,
            GraphNodeKind::IsolatedContact => GlobalEventKind::IsolatedContact,
            GraphNodeKind::PeriodicAnchor => GlobalEventKind::VirtualPeriodicAnchor,
        };
        output.push(ReconciledGlobalEvent {
            id: GlobalEventId(output.len() as u32),
            point: representative.point,
            kind,
            members: members.iter().map(|&i| events[i].reference).collect(),
        });
    }
    GlobalNodeReconciliation::Ready(GlobalNodeReconciliationPlan {
        events: output,
        source_local_nodes: events.len(),
        can_apply: false,
    })
}

impl GlobalNodeReconciliation {
    /// Certifies the exact stored plan against freshly inspected source graphs.
    /// A modified plan cannot be applied or mistaken for a valid preview.
    pub fn validate_against(
        &self,
        graph: &MultiFaceIntersectionGraph,
        first: &BrepModel,
        second: &BrepModel,
    ) -> Result<(), GeometryError> {
        let expected = Self::build(graph, first, second)?;
        if self != &expected {
            return Err(invalid("global-node plan differs from source-bound reconstruction"));
        }
        Ok(())
    }

    pub fn build(
        graph: &MultiFaceIntersectionGraph,
        first: &BrepModel,
        second: &BrepModel,
    ) -> Result<Self, GeometryError> {
        graph.validate_against(first, second)?;
        if graph.state == MultiFaceGraphState::Held {
            return Ok(Self::Held(GlobalNodeReconciliationHold::UpstreamAmbiguity));
        }
        if graph.pairs.len() > MAX_PAIRS {
            return Ok(Self::Held(GlobalNodeReconciliationHold::BudgetExceeded));
        }

        let mut prepared = Vec::new();
        for pair in &graph.pairs {
            for node in &pair.graph.nodes {
                if prepared.len() >= MAX_EVENTS {
                    return Ok(Self::Held(GlobalNodeReconciliationHold::BudgetExceeded));
                }
                let mut edges = BTreeSet::new();
                let mut vertices = BTreeSet::new();
                for source in &node.boundary_sources {
                    let model = match source.operand {
                        GraphOperand::First => first,
                        GraphOperand::Second => second,
                    };
                    let physical_edge = model.edges.get(source.edge.0 as usize)
                        .filter(|edge| edge.id == source.edge)
                        .ok_or_else(|| invalid("missing source physical edge"))?;
                    edges.insert((source.operand, source.edge));
                    // Different edge IDs meeting at an original material
                    // vertex are a valid shared-topology witness. Nearness
                    // to an unrelated vertex is never considered here.
                    for id in [physical_edge.start, physical_edge.end] {
                        let vertex = model.vertices.get(id.0 as usize)
                            .filter(|vertex| vertex.id == id)
                            .ok_or_else(|| invalid("missing source material vertex"))?;
                        if distance(node.point, vertex.position) <= node.tolerance.length {
                            vertices.insert((source.operand, id));
                        }
                    }
                }
                prepared.push(PreparedEvent {
                    reference: GlobalNodeRef { pair: pair.id, local: node.id },
                    point: node.point,
                    kind: node.kind,
                    faces: [pair.first_face, pair.second_face],
                    resolution: node.tolerance.length,
                    edges,
                    vertices,
                });
            }
        }
        Ok(reconcile_prepared(&prepared))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FacePairId, GraphNodeId};

    fn point(x: f64) -> Point3 {
        Point3 { x, y: 0.0, z: 0.0 }
    }
    fn endpoint(pair: u32, x: f64, first: u32, second: u32, witness: EdgeId) -> PreparedEvent {
        PreparedEvent {
            reference: GlobalNodeRef {
                pair: FacePairId(pair),
                local: GraphNodeId(0),
            },
            point: point(x),
            kind: GraphNodeKind::Endpoint,
            faces: [FaceId(first), FaceId(second)],
            resolution: 1.0e-8,
            edges: [(GraphOperand::Second, witness)].into_iter().collect(),
            vertices: BTreeSet::new(),
        }
    }

    #[test]
    fn adjacent_faces_with_one_shared_physical_edge_weld() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(1, 2.0e-10, 0, 2, EdgeId(9));
        let GlobalNodeReconciliation::Ready(plan) = reconcile_prepared(&[a,b]) else {
            panic!("shared physical edge should reconcile");
        };
        assert_eq!(plan.events.len(), 1);
        assert_eq!(plan.events[0].members.len(), 2);
        assert_eq!(plan.events[0].kind, GlobalEventKind::BoundaryEndpoint);
        assert!(!plan.can_apply);
    }

    #[test]
    fn identical_coordinates_without_shared_source_topology_are_held() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(1, 0.0, 0, 2, EdgeId(10));
        assert!(matches!(
            reconcile_prepared(&[a, b]),
            GlobalNodeReconciliation::Held(
                GlobalNodeReconciliationHold::MissingSharedTopologyWitness)
        ));
    }

    #[test]
    fn tolerance_chain_fails_closed_even_with_edge_provenance() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(1, 0.6e-8, 0, 2, EdgeId(9));
        let c = endpoint(2, 1.2e-8, 0, 3, EdgeId(9));
        assert!(matches!(
            reconcile_prepared(&[a, b, c]),
            GlobalNodeReconciliation::Held(
                GlobalNodeReconciliationHold::NonCliqueOrDuplicatePair)
        ));
    }

    #[test]
    fn periodic_chart_anchors_never_become_physical_vertices() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let mut b = endpoint(1, 0.0, 0, 2, EdgeId(9));
        b.kind = GraphNodeKind::PeriodicAnchor;
        assert!(matches!(
            reconcile_prepared(&[a,b]),
            GlobalNodeReconciliation::Held(
                GlobalNodeReconciliationHold::MissingSharedTopologyWitness)
        ));
    }

    #[test]
    fn disparate_precision_budgets_are_held_in_uncertain_band() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let mut b = endpoint(1, 1.5e-8, 0, 2, EdgeId(9));
        b.resolution = 2.0e-8;
        assert!(matches!(
            reconcile_prepared(&[a,b]),
            GlobalNodeReconciliation::Held(
                GlobalNodeReconciliationHold::UnresolvedNearCoincidence)
        ));
    }

    #[test]
    fn material_vertex_from_different_incident_edges_is_shared() {
        let mut a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let mut b = endpoint(1, 1.0e-10, 0, 2, EdgeId(10));
        a.vertices.insert((GraphOperand::Second, VertexId(4)));
        b.vertices.insert((GraphOperand::Second, VertexId(4)));
        let GlobalNodeReconciliation::Ready(plan) = reconcile_prepared(&[a,b]) else {
            panic!("one original material vertex must reconcile incident edges");
        };
        assert_eq!(plan.events.len(), 1);
        assert_eq!(plan.events[0].members.len(), 2);
    }

    #[test]
    fn changing_both_faces_requires_witnesses_on_both_operands() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(1, 0.0, 3, 2, EdgeId(9));
        assert!(matches!(
            reconcile_prepared(&[a,b]),
            GlobalNodeReconciliation::Held(
                GlobalNodeReconciliationHold::MissingSharedTopologyWitness)
        ));
    }

    #[test]
    fn sufficiently_separated_events_stay_independent() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(1, 2.0e-8, 0, 2, EdgeId(9));
        let GlobalNodeReconciliation::Ready(plan) = reconcile_prepared(&[a,b]) else {
            panic!("separate events outside both tolerance envelopes should remain separate");
        };
        assert_eq!(plan.events.len(),2);
    }

    #[test]
    fn solitary_periodic_anchor_is_virtual_and_stays_unwelded() {
        let mut a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        a.kind = GraphNodeKind::PeriodicAnchor;
        let GlobalNodeReconciliation::Ready(plan) = reconcile_prepared(&[a]) else {
            panic!("a periodic chart anchor is a permitted nonphysical preview event");
        };
        assert_eq!(plan.events[0].kind,GlobalEventKind::VirtualPeriodicAnchor);
        assert!(!plan.can_apply);
    }

    #[test]
    fn same_pair_events_must_not_weld() {
        let a = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let b = endpoint(0, 0.0, 0, 1, EdgeId(9));
        let GlobalNodeReconciliation::Ready(plan) = reconcile_prepared(&[a,b]) else {
            panic!("per-pair node identity is already certified upstream");
        };
        assert_eq!(plan.events.len(), 2);
    }
}
