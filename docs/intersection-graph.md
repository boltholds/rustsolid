# RustSolid v0.15 — analytic IntersectionGraph and provenance

The read-only `IntersectionGraphBuilder` converts the **bounded** v0.14
`TrimmedIntersectionEngine` result into a typed event graph. It does not edit
either source body, create Boolean results, commit topology deltas or claim
certified completeness for NURBS intersections.

## Contract

- `IntersectionGraph::state` is `Empty`, `Ready` or `Ambiguous(reason)`.
  Empty and unresolved results expose **no candidate cut topology**. Unsupported
  trimmed face configurations return an explicit error.
- `GraphFaceSource` namespaces local `FaceId` and `SurfaceId` by
  `GraphOperand::{First,Second}` and retains their individual tolerance policies.
  A face with the same numeric ID in another body remains distinct.
- Every `GraphArc` stores its analytic `IntersectionCurve`, bounded parameter
  interval, endpoints by typed `GraphNodeId`, `TrimContact` class and **both**
  parent face sources. `GraphBoundarySource` records actual coedge/edge
  incidences where each endpoint lies on a physical face boundary.
- `GraphNode::membership` classifies each vertex against **both actual UV
  trim domains** as `Interior` or `Boundary`. `Outside` is never accepted.
  Point-only contacts retain their `TrimContact` in the node.
- `GraphNodeKind::PeriodicAnchor` marks an artificial parameter origin on a
  **closed self-loop** (degree two); it is not an instruction to split the
  physical circular edge at the cylinder's UV seam.
- When a partial periodic arc crosses the `0 ↔ 2π` chart boundary, its first
  and last fragments are combined into **one unwrapped analytic interval**.
  No chart seam becomes a physical graph junction.
- `GraphTolerance::length` combines the engine tolerance and each *face's*
  tolerance relative to local feature extent. It governs potential vertex
  merging and detecting collapsed spans. `evaluation_budget` additionally
  allows representability at far-away global coordinates but **never governs
  topology merging**. Numerical ambiguity near topology-changing contacts is
  preserved for downstream refinement.
- `IntersectionGraph::validate_against` checks source face identities, the
  provenance policy, graph references, edge/node adjacency, analytic endpoint
  and interior samples on both supporting surfaces, **membership in both
  face trims** and the actual source edges near boundary points. Corrupting
  any of these fails the validation test.

## Independent correctness and regression policy

`tests/fixtures/intersection_graph.v1.json` pins 12 original B-rep scenarios:
line segments, isolated corner contact, boundary overlap, disjoint face trims,
coincident supports, complete and half-circle, two/tangent cylinder generators,
near-tangent uncertainty, periodic cap/side contact, and heights outside the
cylinder. `tests/intersection_graph.rs` has independent analytic expectations
for coordinates and arc lengths, operand swap, translation,
UV reparameterization, provenance tampering and uncertainty. The suite uses
no Parasolid or OpenCascade outputs as truth.

A finite set of samples only establishes measured support/trim consistency,
not a general mathematical proof of completeness of the analytic intersection.
For a future `TopologyDelta`, require both a geometry-based classification
and a separate topological validity proof before mutating a B-rep.

## Example

```rust
use rustsolid::{IntersectionGraphBuilder, GeometryTolerance};

let graph = IntersectionGraphBuilder::new(GeometryTolerance::default())?
    .intersect_faces(&first_brep, first_face, &second_brep, second_face)?;
if graph.kind() == "ready" {
    graph.validate_against(&first_brep, first_face, &second_brep, second_face)?;
    for arc in &graph.arcs {
        let geometric_midpoint = arc.evaluate(0.5)?;
        // Use provenance and classification; do not directly mutate topology.
    }
}
```

The `geometry.v1` CLI `intersection_graph` operation exposes read-only
metadata, nodes, provenance, analytic arcs and a limited number of sample
coordinates. It has no edit side effects.

## Next independently gated stage

Convert a validated `IntersectionGraph` into a **preview-only, typed
`TopologyDeltaPlan`**: vertices to insert, edges/faces to split, face-region
classification, retained/created/retired carrier mapping and explicit uncertainty.
The plan must not be committed until sparse generational arenas, local
tolerance propagation, nontrivial shell topology and transaction journaling
are fully supported.
