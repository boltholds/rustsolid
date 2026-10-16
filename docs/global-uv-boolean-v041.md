# RustSolid v0.41 — Global arrangement / volume-changing Boolean gate

Status: implementation specification; **not a completed Boolean**.

## Baseline inspected
- Cargo.toml currently declares v0.40.0.
- src/planar_uv_arrangement.rs supports one straight cut of one convex planar face with circular holes; its plan is a preview.
- src/coedge_reconciliation.rs binds a spatial intersection edge to two face-specific pcurves; its plan is not a commit.
- src/topology_delta_plan.rs explicitly forbids applying its preview to a B-rep.
- Existing worktree contains staged/uncommitted v0.40 changes; do not reset, commit or overwrite these changes as part of v0.41.

## Scope of first executable vertical slice
Operands: watertight oriented solids A and B, with planar faces and cylindrical cutters in a bounded supported configuration. A must contain a pre-existing circular through-hole; B is an overlapping cylinder. Difference A\\B must remove strictly positive material and intersect the original hole boundary transversely. Any unresolved coincidence, tangency, nonmanifold vertex, self-intersection, missing loop or invalid shell yields Held and preserves A and B byte-for-byte.

## Required pipeline
1. Collect all face-pair intersections, retaining operand/face namespaces and analytic 3D carriers.
2. Cluster intersection events with explicit tolerance bounds; never merge distinct intersections merely because their UV projections coincide on different charts.
3. Build an arrangement **per affected face** from all its intersecting curves, source trims and periodic seams. Split source trims and curves at every event; preserve analytic circles/arcs.
4. Trace closed directed cycles, determine nesting/holes, and construct material regions. Require signed-area conservation within a documented tolerance.
5. Classify representative points of each open region against the **whole opposite solid**, with boundary ambiguity reported as Held. Confirm orientation using local normals; do not classify only against one opposing face.
6. Select regions for Difference (A outside B; B inside A with reversed orientation), stitch shared physical edges and face-specific pcurves. Reject nonmanifold incidence.
7. Build candidate B-rep in temporary storage, assign persistent IDs only at commit, validate and atomically swap into the document.

## Invariants and independent checks
- Each ordinary manifold edge has exactly two oppositely oriented coedge uses; seam uses must be explicitly paired.
- Face loops close in UV and evaluate onto the same 3D carrier at multiple samples and endpoints.
- Closed shell is orientable and has no unpaired boundary; compute Euler characteristic per shell/component and compare to its genus, not to a hard-coded 2.
- Analytic volume for the chosen fixture must be independently derived (e.g. inclusion-exclusion with controlled quadrature), not reused from the same Boolean implementation.
- Result must satisfy 0 <= volume(A\\B) <= volume(A), and strictly less than volume(A) for the transverse fixture.
- On Held or validation error, original operands, revision and undo stack are unchanged.
- Undo and redo must reproduce identical geometry and topology, within stable ID policy.

## Proposed contracts
```rust
pub enum GlobalBooleanResult {
    Ready(ValidatedCandidateBrep),
    Held(GlobalBooleanHold),
}
pub enum GlobalBooleanHold {
    UnsupportedSurface,
    UnresolvedIntersection,
    AmbiguousRegionClassification,
    InconsistentUvArrangement,
    NonManifoldIncidence,
    InvalidCandidateShell,
    VolumeMismatch,
}
```
A Ready candidate is **not** permission to mutate the document until an atomic transaction validates the source revision and swaps the candidate.

## Acceptance tests
1. Existing hole + overlapping cylindrical subtractor: positive removal, correct number of components, watertight.
2. Disjoint subtractor: unchanged volume and topology.
3. Fully contained subtractor: additional cavity/shell as appropriate.
4. Tangency: Held (until exact contact policy exists).
5. Coincident cylinder boundaries: Held.
6. Two cuts on same physical edge: correctly ordered splits or Held.
7. Periodic seam crossing: one physical edge, chart-specific pcurves.
8. Reverse input orientation: reject invalid operands or normalize before operation.
9. Corrupted source coedge: Held/error, no mutation.
10. Failed candidate validation: no revision or undo history change.
11. Undo/redo: topology and volume restored.

## Delivery gates
G1: read-only multi-face intersection graph with tests.
G2: multi-trace UV arrangement, area/loop invariants.
G3: whole-solid region classifier, ambiguity tests.
G4: candidate B-rep builder and independent volume oracle.
G5: atomic transaction and Undo/Redo.

Do not label v0.41 complete until G1–G5 execute and pass on the repository test runner.

## G1 implementation checkpoint (2026-10-14)
- `src/multi_face_intersection_graph.rs` collects pairwise immutable graphs with operand-specific face/arc/node references. Its validator now regenerates every face-pair intersection using the stored construction tolerance; fabricated Empty pairs therefore cannot conceal actual crossings.
- `src/global_vertex_reconciliation.rs` combines only **open-arc endpoint** events from distinct pairs. On every changed operand, adjacency must be backed by the same physical source edge or the same original material vertex (derived from certified boundary witnesses). Coincident XYZ coordinates alone cannot establish equivalence.
- Mixed precision envelopes, incompatible topology witnesses, non-clique proximity chains, and synthetic periodic anchors near physical endpoints fail closed with `Held`.
- Point-only contacts and periodic chart anchors remain non-material events. Interior arc crossings are not yet split; they require a G2 multi-trace arrangement.
- Rebuilt/source-checked preview plans keep `can_apply = false`. No persistent B-rep IDs are allocated.
- Coverage: `tests/multi_face_intersection_graph.rs`, `tests/global_vertex_reconciliation.rs`, plus internal reconciliation unit tests.
- **Validation blocker:** `cargo_test` returns `Permission denied (os error 13)` while opening `/workspace/projects/rustsolid/Cargo.toml`. Code has been written, but compilation and runtime pass status are **not verified**. Do not mark G1 complete until the runner's file permissions are repaired and tests succeed.

## G2 implementation checkpoint (2026-10-14)
- Added `src/multi_trace_uv_arrangement.rs` with bounded line/line and line/circle UV event construction, original circular arc preservation, cycle stitching, nesting of uncut circular holes, area and directed incidence checks.
- Added 8 integration and 3 internal test scenarios. The two-line/three-hole geometric orientation and area were independently smoke-checked in a mathematical prototype, not in the Rust test runner.
- The plan is a read-only source-reconstructed preview; cut pieces retain pair-scoped `GlobalArcRef`; `can_apply=false`.
- Coincident lines, unresolved tangencies, truncated cutter arcs, unverified multiway intersections and topology outside the bounded slice return Held.
- Neither default nor nested `cargo_test` could read its manifest due to workspace permissions; G2 is **not yet compiled/tested**.

## G3 implementation checkpoint (2026-10-15)

- Added `src/whole_solid_region_classification.rs`: source-validated whole-body classification for each complete planar G2 UV cell, using analytic boundary-extrema witnesses rather than a single representative point.
- Bounded opposite solids: validated AABB polyhedron and finite analytic cylinder. Unsupported B-reps or unresolved bounds fail closed.
- The selected graph pair list must include every face of the opposite whole body, including Empty pairs.
- For axis-aligned box boundaries, one-f64-ULP UV→3D overshoot is accepted **only when the G2 cut is physically witnessed on that exact opposing face**; unmatched near-face slivers remain Held. `WholeSolidCellCertificate` retains the witness bitmap.
- Difference/Union/Intersection side selection is explicit; difference retains opposite-body Inside patches with reversed normal. `can_apply=false`.
- Ran `cargo test --test whole_solid_region_classification` on the source-matched temporary Linux Docker copy: **5 passed, 0 failed**. G2 tests: **8 passed, 0 failed** after correcting an exact-tangency expectation; geometric construction still Held on tangency.
- Full suite revealed a preexisting G1 weakness: disjoint finite bodies could produce Held because untrimmed supporting planes coincide. Added a strict, validated body-bbox separation early-out in `IntersectionGraphBuilder` (uncertain gaps still follow analytic intersection). Targeted `multi_face_intersection_graph` tests now pass **3/3**.
- Full suite must be rerun with the G1 correction before G3 can be marked regression-complete. AI Workbench's own cargo runner continues to lack access to Cargo.toml.

The RustSolid package version remains v0.40 until G4/G5 are shipped and a release is intentionally cut.

## Full regression verification (2026-10-15)

The temporary source-matched Linux Docker checkout completed `cargo test --quiet` with **exit code 0**: **530 passed, 0 failed, 62 test groups** (including G1/G2/G3 and preexisting RustSolid suites).

The verification checkout used the same `06ea004cb8e63bb6f53e6aeaeec83b5806258cfe` base revision plus the AI Workbench staged/unstaged Rust source and integration test files. It was an isolated copy; no commit/push/reset was made in the authoritative project.

AI Workbench's own `cargo_test` tool still reports a permissions error opening its mounted Cargo manifest, so this pass is specifically attributed to the copied checkout's real Rust 1.88 Linux test execution.

G4 topology reconstruction and G5 atomic Boolean commit remain future work. `can_apply` is still false on the G3 result.


## G4a implementation checkpoint (2026-10-15)

- Added a **read-only, face-local** `FacePatchAssemblyPlan` with source
  coedge provenance, analytic 3D graph line witnesses, directed retained
  loops, and source-/cross-operand cut references.
- Selected G3 material cells are required to retain precisely one side
  of each local cut; the two G2 uses are **not** assigned as physical
  B-rep radial twins. Source trims are independently reconstructed from
  their original 2D pcurves.
- The effective outward shell orientation combines original face sense
  with the reversed subtractor patches used for Difference.
- Regression suite: `tests/face_patch_assembly.rs`; design and
  limitations: `docs/face-patch-assembly-g4a-v041.md`.
- **G4 is not complete.** External face-to-face coedge pairing, complete
  solid-face scheduling, temporary manifold shell construction,
  independent volume oracle and G5 atomic commit/Undo/Redo remain.
  The G4a `can_apply` flag is always false.

## G4b1 implementation checkpoint (2026-10-15)

- `src/global_edge_stitch.rs` pairs selected source-backed G4a coedges
  from multiple faces by actual original coedge twins or by the same
  certified A/B face-pair analytic graph arc. Its temporary edge IDs
  are never assigned as persistent `EdgeId`.
- Physical intersection twins must be on different source operands,
  have opposing 3D directions, identical original bounded graph arc
  provenance and physical endpoint witnesses. Missing or incompatible
  uses are recorded or Held; no synthetic shell faces are invented.
- G2 now ignores a `BoundaryOverlap` as an interior cut **only** when
  both endpoint graph witnesses refer to the same original source
  coedge AND every checked point is on the bounded 3D line carrier.
  Actual tangencies and fabricated source evidence still fail closed.
- `src/radial_shell_audit.rs` proves closed **combinatorial** radial
  incidence, oriented vertex fans and integer Euler/genus for source-
  complete candidate patch sets. Synthetic tetrahedra and disjoint
  shells exercise the positive audit; open loops, self-paired edges
  and coincident unlinked vertices are rejected.
- Current real G4a fixture builds and source-validates a physical
  cross-operand edge using distinct planar source faces. Their many
  other unchanged/partially trimmed faces are not yet G4a-classifiable,
  so `missing_source_faces` remains nonempty, the radial shell audit
  does **not** claim closure, and `can_apply=false`.
- Next G4b2: retained original/unchanged face scheduling, global
  cut-vertex reconciliation and canonical closed `BrepModel` candidate.
  G4c/G5 volume oracle and atomic transaction remain outstanding.

Tests: `tests/global_edge_stitch.rs` (6 source-backed integration
cases), plus five `radial_shell_audit` unit checks. See
`docs/global-edge-stitch-g4b1-v041.md`. The CI suite succeeds on
`feat/g4b-cross-face-edge-stitch-20261015` without modifying `main`.

## G4b2a implementation checkpoint (2026-10-16)

- Added `src/canonical_box_candidate.rs`: an actual immutable
  `BrepModel` candidate for the deliberately limited **two exact
  axis-aligned box** Boolean class. The existing source-certified
  `BooleanKernel::apply_solids` performs the planar arrangement and
  actual coedge construction; G4b2a independently audits its result.
- Every one of the twelve original source faces is partitioned on
  both source boxes' EXACT supporting-plane grid, classified by
  Boolean occupancy on both sides of every tile, and matched to
  exactly one result `FaceId` when exposed. The face coverage ledger
  distinguishes untouched, discarded and partitioned faces with
  conservative surface area accounting and right-Difference senses.
- Every resulting `VertexId` is tied to a unique triple of exact
  original XYZ plane indices, without proximity welding. The actual
  B-rep's coedges are independently traversed by the G4b radial
  vertex-fan, Euler and genus auditor. Source revision and tolerance
  are bound to preview validation.
- A tiny box-specific inclusion/exclusion cross-check verifies
  volume against the original bounds; it is not a general curved
  G4c mass oracle.
- This is the first **real validated closed canonical candidate**
  within this strict class: `closed_shell_proven=true` only for
  a single genus-zero AABB-CSG shell, but `can_apply=false`.
- Unlike the new source-independent G4b1 edge-matching interface,
  this specialized path does not yet schedule/cut arbitrary NURBS
  or cylindrical surface patches. It **reuses** the existing planar
  Boolean engine rather than claiming G1/G2/G3 generality.
- Remaining: source-bound global carrier assembly for arbitrary
  planes/cylinders, periodic edge reconciliation and robust 3D
  cut-vertex identity, independent curved mass oracle, then G5
  revision-bound atomic commit with Undo/Redo.

Tests: `tests/canonical_box_candidate.rs` and
`docs/canonical-box-candidate-g4b2a-v041.md`. This is an experimental
feature branch; do not label v0.41 complete or merge to stable CAD
until general G4/G5 are implemented and reviewed.

## G4b2b source-certified analytic carrier checkpoint (2026-10-16)

The bounded curved source-preview and full graph carrier layers are now
composable, but this is **NOT yet a general solid/solid Boolean**.

- `AnalyticBoreCandidate` independently certifies the existing
  **true closed** through-hole (genus 1) and blind-bore (genus 0)
  source constructors for a strict box minus an upright cylindrical
  cutter. It verifies both physical circular rim coedges, the
  period-`2π` same-face seam, all nine original source-face
  dispositions, original box carriers, cut-floor ownership,
  and correct material normals; every cap-to-cylinder closed graph
  retains a virtual, nonmaterial periodic anchor.
- `AnalyticSourcePcurvePair` preserves the independent original
  **plane Harmonic Curve2** and **cylinder CylinderSection Curve2**
  for each closed source graph and validates that they map to the
  same 3D analytic curve, even under cutter-frame rotation.
- `GlobalAnalyticCutPlan` consumes ALL source face-pair graphs,
  validates the G1 global event clustering against physical source
  edges/vertices, and publishes every transverse source-owned
  analytic arc with two UV carriers and source-checked global endpoint
  event identities. It works for the actual **partially intersecting
  box/cylinder** fixture spanning multiple source faces. Local
  `GlobalArcRef` remains pair-scoped; periodic anchors never
  become physical vertices.
- The through/blind results have real canonical closed B-rep
  topology and strict source-specific analytic mass checks. The
  partial multi-face registry has only verified **cut curves and
  endpoint identities**, **not** material-region coverage or a
  closed B-rep. Unsupported tangencies and insufficient sources
  return Held; `can_apply=false` across the G4b2b previews.
- Next G4b2b3: use these source-derived analytic curves and
  G1 events to split the *real* affected source trims globally,
  assemble oriented multi-face B-rep regions, reconcile periodic
  coedges, and certify every radial vertex fan with the curved
  shell/genus checks. Only then can general G4c volume and G5
  source-revision atomic commit/Undo/Redo be claimed.

Tests: `tests/analytic_bore_candidate.rs` (9 integration),
`tests/global_analytic_cut_registry.rs` (5 integration),
`docs/analytic-bore-candidate-g4b2b-v041.md`, and
`docs/global-analytic-cut-registry-g4b2b-v041.md`.
