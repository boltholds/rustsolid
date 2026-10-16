# RustSolid v0.41 — G4b2b3a: source-bound physical edge subdivision

**Status: real topology-preserving subdivision of both original B-reps.**
This is an intermediate carrier for G4b2b3, not the Difference result,
not a selected curved face-patch shell, and not G4c/G5.

The preceding G4b2b2 GlobalAnalyticCutPlan has exact source-certified
intersection curves and G1 physical boundary endpoint identities, but
has no material edge, vertex or coedge objects for partially intersected
solids. This slice turns the certified boundary events into a new
source-copy topological subdivision, without changing the originals.

GlobalAnalyticEdgeSplitPlan::build(left,right,tolerance):

1. Reconstructs and independently validates the complete G4b2b2
   source cut registry and all face-pair graphs.
2. Takes only endpoints of registered transverse open arcs, expands
   their real G1 member-node source coedge witnesses, and groups by
   (GraphOperand, original EdgeId, GlobalEventId). No XYZ-only weld is
   permitted. Periodic synthetic anchors do not create vertices.
3. For each originally straight, nonclosed physical source edge,
   analytically projects G1 events onto its original bounded Curve3
   and rejects off-edge, unresolved near-end or duplicate events.
   An event at an existing source vertex retains that VertexId; a
   strictly interior event creates one new source-local VertexId.
4. Splits each physical edge into ordered Edges and splits BOTH
   directed coedges with independent original Curve2 domains,
   respecting reversed coedge orientation.
5. Splices new directed uses into affected Loops with reciprocal
   next/prev links and two opposite twin coedges per physical edge.
6. Calls BrepModel::validate() on BOTH split source models, proving
   closed source incidence, UV/3D carrier agreement, unchanged Euler
   characteristic and unchanged source mass metadata.
7. validate_against reconstructs the split from current sources and
   compares the exact IDs, segments, domains and resulting B-rep models.

The source operands remain unchanged. All new IDs belong to a
throwaway preview, not persistent document topology.

## Fail-closed limits

- A periodic original circular edge with an interior event is Held:
  it needs an angular split policy and independent recharting.
- Near-duplicate G1 events and conflicting original-vertex IDs are Held.
- An open arc endpoint without original B-rep edge witness is Held.
- Intersection curves still need retained material UV regions and
  true cut-edge coedge twins on both different source faces.
- These are complete, subdivided ORIGINAL source shells. They are not
  a Difference result, so boolean_shell_proven=false and can_apply=false.
- G4b2b3b must finish retained curved face patches, physical cut-edge
  twins, radial vertex fans and multiple connected components. G4c/G5
  remain downstream.

## Regression cases

tests/analytic_source_edge_split.rs: partial box-cylinder source
subdivisions, true coedge domains, event identity reuse, Euler/mass
invariants, tamper and source revision invalidation, disjoint cases.

Run with cargo test --test analytic_source_edge_split, then
cargo test --all-targets.
