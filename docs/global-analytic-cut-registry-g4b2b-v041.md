# RustSolid v0.41 — G4b2b2: global analytic cut carrier registry

**Status: globally source-certified cutting-curve preview. Not a completed
curved Boolean and NOT a closed candidate B-rep.**

This stage complements the *real closed* through-hole and blind-hole
`AnalyticBoreCandidate` in G4b2b1. It covers the harder case where one
cylinder intersects several faces of a different source body and cutting
curves must share genuine material vertices across different face-pair
intersection graphs.

`GlobalAnalyticCutPlan::build(left,right,tolerance)` operates on **all
faces** of both source `BrepBody` operands, not an arbitrary selected
subset.

It composes existing first-party kernels:

1. `MultiFaceIntersectionGraphBuilder` evaluates every source
   `FaceId(A) × FaceId(B)` pair, including ones reporting Empty. Its
   source-bound validator independently recomputes every pair.
2. `GlobalNodeReconciliation` builds a G1 identity for each true
   intersection endpoint, merging events across face-pair graphs only
   when original material source edges/vertices prove their shared
   identity. Distinct nearby XYZ coordinates are **never** enough.
   Non-clique tolerance chains, unsupported adjacency and ambiguous
   contacts return Held.
3. Each `TrimContact::Crossing` arc retains a pair-scoped
   `GlobalArcRef` and both original source face IDs. Its independent
   `IntersectionPcurve::to_curve2` maps the exact analytic 3D carrier
   to **two different source UV charts** (including harmonic plane
   sections, cylindrical generators and unwrapped cylindrical U).
4. At seven locations, `Curve3` and both new `Curve2` carriers are
   independently evaluated on their real source `Surface3` to ensure
   the same 3D point. No triangulation, extrapolated line or stitched
   render mesh supplies the curve geometry.
5. Open-arc endpoints refer to certified `GlobalEventId`s from G1.
   A periodic closed curve has **no physical endpoint ID**: its chart
   anchor stays virtual. BoundaryOverlap is recorded separately and
   never promoted to a transverse new Boolean edge.
6. `validate_against` rederives all source face graphs, G1 endpoint
   identities and 3D/UV curve carriers, and checks source revision
   stamps and geometric tolerance.

The result `SourceCertifiedCutArc` is a source-owned cut *intent*,
not a persistent `EdgeId` or `CoedgeId`. Real topology construction
must still choose retained material sides and link these verified
intersections with G2/G3 source trim regions.

## Geometry and failure policy

The acceptance fixture is a box with faces `x,z=±1.5` crossing a
finite radius-2 upright cylinder over `y=1..4`. Multiple different
face-pair graphs share genuine physical endpoints where source box
edges intersect cylinder sections. Tests require that G1 recognizes
those shared events *using source topology*, and all source plane/
cylinder pcurve witnesses reconstruct their analytic 3D carriers.

Disjoint solids return `Held(NoCertifiedCrossing)`. An uncertain
near-tangent cylinder returns `Held(UpstreamGraphAmbiguity)`.
Changed body revisions, moved cutter geometry, edited graph arc IDs
and modified source event coordinates invalidate a saved plan.

`closed_shell_proven=false` and `can_apply=false` for all registry
results. G4b2b3 must still split partially intersected source faces
at all G1 events, weld only evidence-backed vertices, resolve period
seams, and build the complete oriented manifold `BrepModel`.
G4c volume/general integration and G5 atomic Undo/Redo are downstream.

## Files and tests

- `src/global_analytic_cut_registry.rs`
- `tests/global_analytic_cut_registry.rs`
- `src/analytic_bore_candidate.rs`
- `tests/analytic_bore_candidate.rs`

```bash
cargo test --test global_analytic_cut_registry
cargo test --test analytic_bore_candidate
cargo test --all-targets
```
