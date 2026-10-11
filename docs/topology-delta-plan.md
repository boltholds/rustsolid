# TopologyDeltaPlan v0.16 — immutable, source-bound geometric intent

This is **not** a Boolean implementation or a mutation API. A verified
`IntersectionGraph` proposes where a later generalized B-rep editor *might*
need vertices, edge subdivision, analytic carriers, and face-region boundaries.
A single pair of intersecting faces does not determine which side of a closed
solid should be kept. No topology is modified and `can_apply` is **always false**.

## Three layers

1. **SurfaceIntersectionEngine** computes exact-in-form analytic supporting
   curves with approximate numerical evaluations.
2. **TrimmedIntersectionEngine / IntersectionGraph** classifies those curves
   against both actual source face domains and records physical boundary
   witnesses, contact types, source face identity, and local tolerances.
3. **TopologyDeltaPlanner** maps verified graph nodes to *operand-namespaced*
   symbolic `PlanVertexPlacement`s, assigns a `PlanArcRole` to every graph arc,
   and proposes per-face `PlanFaceEffect` values and typed prerequisites.
   `TopologyDeltaPlan::validate_against` regenerates a plan from the still
   validated sources and rejects missing, tampered or stale intents.

### Per-node decisions

- `ReuseVertex` identifies an already existing source vertex (never creates
  two vertices for the same physical endpoint within one operand).
- `SplitEdge` records a source edge, its directed coedge, and a bounded
  floating-point carrier parameter/fraction, **without executing a split**.
- `CreateFaceVertex` signals a cut endpoint in the interior of one source
  face. It is *not* a valid isolated Euler split.
- `VirtualPeriodicAnchor` represents a closed loop's arbitrary angular
  parameter origin. It is not a physical cut, even on a cylinder seam.
- `ContactOnly` signals isolated point contacts or unclassified tangencies.

Both operands are kept separate using `GraphOperand`; `FaceId(0)` or
`EdgeId(0)` from unrelated bodies never silently refer to one entity.

### Per-face decisions

| Intent | Precondition |
|---|---|
| `CandidateBoundarySplit` | An open *crossing* arc terminates on the physical boundary at both endpoints of this face. |
| `PendingOpenTrace` | At least one endpoint is in the face interior; further curves / region classification are necessary. |
| `CandidateInnerLoop` | A closed crossing curve requires a new trim/region; a periodic anchor must not become an actual split edge. |
| `ReconcileExistingBoundary` | A curve overlaps a physical face boundary. Treat as a boundary matching problem, not an unconditional Boolean cut. |
| `RecordTangency` | Tangency requires dedicated topology classification and may need no cut at all. |

The plan also records why it **cannot be committed**: generalized writable
B-rep arenas, region classification, inner trim-loop creation, open trace
completion, multiple splits on the same edge, curved-edge subdivision,
boundary reconciliation, tangency or numerical refinement.

An `Empty` graph is `NoChange`; an ambiguous graph is `Held`. Isolated point
contacts are `NoChange` with preserved contact evidence. Crossing arcs give
`PreviewOnly`, while tangency and overlapping boundary arcs give `Held`.
None of these state labels authorizes model mutation.

## Independent regression cases

`tests/fixtures/topology_delta_plan.v1.json` pins 12 analytic input scenarios;
`tests/topology_delta_plan.rs` constructs the actual B-reps and checks
classifications, topology actions, source identity, tolerance requirements,
closed-loop representation, mathematical coordinates, operand symmetry,
reproducibility, and tampering detection. No third-party CAD engine is used
as a correctness oracle.

Run:

```bash
cargo test --test topology_delta_plan
cargo test --all-targets
```

## JSON API

`geometry.v1` adds `topology_delta_plan`, accepting the same `first` and
`second` face operands as `intersection_graph`. It returns an immutable plan:
node bindings (including carrier parameters), analytic arc intentions,
per-face changes, outstanding requirements and `can_apply: false`.
`intersection_graph`, `intersect_faces`, `command_history` and legacy commands
retain their request/response shapes.

## Next gate, NOT covered by v0.16

Build a generalized mutable B-rep arena with stable generational IDs, locally
scoped tolerance propagation, robust face-region classification, provenance
for retained/new/retired edges and surfaces, and atomic journaled transactions.
Every real commit must validate a candidate body, and a failed operation must
leave the old model byte-for-byte equivalent. This preview plan is deliberately
not used to perform unsafe partial Boolean edits.
