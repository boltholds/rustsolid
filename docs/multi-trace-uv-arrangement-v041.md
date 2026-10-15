# RustSolid v0.41 / G2 — Multi-Trace UV Arrangement

Status: implementation written in the AI Workbench `rustsolid` worktree, **compiled and regression-tested on a source-matched isolated Linux/Docker checkout**. G2 remains preview-only, not a Boolean commit.

## Public source

- `src/multi_trace_uv_arrangement.rs`
- `MultiTraceUvPlan::build(first, second, full_graph, operand, face, tolerance)`
- `MultiTraceUvPlan::build_selected(first, second, source_pairs, operand, face, tolerance)`
- `MultiTraceUvPlan::validate_against(...)`
- `tests/multi_trace_uv_arrangement.rs`

All outputs are immutable previews: `can_apply = false`. They do not allocate persistent B-rep IDs, classify volumetric material, or perform Boolean commits.

## Supported bounded geometry

- Strictly convex outer planar face bounded by original `Curve2::Line` coedges.
- Up to eight distinct transverse supporting lines.
- Up to 64 disjoint circular inner loops represented as original `Curve2::Circle` coedges.
- Each source line trace must cover all **material** intervals where its supporting line crosses the entire planar face. The true graph arcs may be interrupted by circular holes.
- Exact-in-form circle and line boundary spans; no circle faceting.
- Cut segments carry a `GlobalArcRef` (face-pair and local arc ID) and a G2-local cut-piece index. One analytic 3D graph arc may produce multiple G2 pieces due to crossings.
- Selected face-pair graphs are rederived against the source solids using the supplied tolerance, including pairs claiming to be Empty.
- Adjacent collinear trace intervals from different face pairs are intentionally **not** welded without G1 source-topology confirmation.

## Algorithm

1. Validate source B-reps and exact selected `IntersectionGraph` witnesses.
2. Read the original planar UV trimming loops; enforce convexity, circle separation and circular containment.
3. Gather and group line arcs, retaining pair-specific provenance.
4. Find analytic line/outer-segment, line/circle and line/line events.
5. Sort/snap only repeated references to one event; reject numerically unresolved near pairs.
6. Split original outer lines, real inner circular loops and supported cut intervals at all events.
7. Add oppositely directed halfedges for each cut interval; source loops are oriented with material to the left.
8. Trace face cycles using the immediately clockwise outgoing halfedge from the *backwards* incoming tangent; do not backtrack along a cut twin.
9. Attach uncut circular inner cycles to the unique containing positive cycle using analytic ray intersections.
10. Check signed area conservation in an anchor-relative UV frame (for translation stability) and verify exactly one forward and one reverse use per G2 cut piece.
11. Return a `Ready` preview or typed `Held`, never a partial topology commit.

## Deliberate Held cases

Coincident or near-parallel overlapping lines, missing physical trace coverage, source vertex contact, circle tangent contact, close unresolved events, nonmanifold radial incidence, overlap between source holes, unverifiable nesting, unsupported nonconvex/curved support, multiway line intersections, unclosed cycles and area mismatch.

The G2 overlay of two crossing lines can reconstruct planar cells even where the 3D triple-intersection event is not yet independently certified. **This is not permission to create a physical B-rep vertex.** G3/G4 must establish material classification and full 3D edge/vertex incidence.

## Test cases authored

- Two parallel cutter faces: three regions, one hole divided by the first line, exact expected areas.
- Circular hole entirely inside the middle region: preserved negative-area `UvCycleKind::Inner`.
- Finite source cutter face that ends inside the target face: `Held(IncompletePhysicalTrace)`.
- Two selected graphs with the same supporting line: `Held` rather than duplicated cuts.
- Every cut piece has its original `GlobalArcRef` and exactly two opposite halfedge uses.
- Tangency on hole boundary remains upstream `Held`.
- Three real holes and two transverse cutter traces preserve all circular fragments and signed area.
- A forged Empty face-pair graph is rejected by source rederivation.
- Module tests exercise orthogonal line intersections with four square cells; an intersection inside a circular hole preserves four quarter-circle arcs; coincident synthetic lines are rejected.

## Build validation history

Both ordinary `cargo_test` and a separately written narrow validation manifest at `validation-g2/Cargo.toml` fail before compilation with:

```
error: failed to read `/workspace/projects/rustsolid/Cargo.toml`
Caused by: Permission denied (os error 13)
```

The alternate manifest produces the identical permissions error for `/workspace/projects/rustsolid/validation-g2/Cargo.toml`. This remains an AI Workbench runner/mount issue. On an isolated checkout matched to the current Rust sources, `cargo test --test multi_trace_uv_arrangement` passed **8/8**, and `cargo test --quiet` passed **530/530 across 62 test groups** with exit code 0 (Rust 1.88, Linux Docker). The exact tangent fixture's upstream graph can be Ready while G2 correctly holds its non-transverse vertex event; the expectation was repaired in the test.

## G2 completion criteria

- Fix AI Workbench's native Cargo runner permissions (the independent checkout already provides executable validation).
- Verified: `cargo test --test multi_trace_uv_arrangement` and full 530-test suite pass on the source-matched isolated runner.
- Cross-check multiple holes and 3D pairing with G1 reconciliation.
- Manually verify B-rep UV orientation and chart seam behavior under translated, rotated and scaled fixtures.
- No topological mutation until G3-G5.
