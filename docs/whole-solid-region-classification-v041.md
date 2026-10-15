# RustSolid v0.41 — G3: Whole-Solid Region Classification

Status: implementation written in AI Workbench. **Not an executable Boolean**; G4/G5 remain outstanding.

## Source layout

- `src/whole_solid_region_classification.rs`
- `tests/whole_solid_region_classification.rs`
- Independent analytic UV line/circle boundary representation: G2 `MultiTraceUvPlan`
- Original whole-body point membership: `BooleanRegionClassifier::classify_point`

The G3 entrypoint `WholeSolidRegionPlan::build(left, right, pairs, uv, operation, tolerance)` classifies each **complete** G2 face cell against the **whole opposing solid**, rather than against a single intersecting face. `validate_against` reconstructs the G2 plan and repeats classification to reject tampered or stale previews.

## Strict prerequisites

1. Both source B-rep bodies validate.
2. A G2 `Ready` analytic planar UV arrangement remains identical to source reconstruction.
3. The supplied source-face pairs include **every opposing face**, including those with Empty intersection graphs. An arrangement constructed from a selected subset of crossing pairs may be geometrically useful, but is insufficient for G3 and returns `Held(IncompleteOppositeFaceCoverage)`.
4. No unresolved G2 topology or ambiguous contact.
5. Opposing whole solid is a certified axis-aligned polyhedral box, or a finite upright/rotated analytic cylinder. Other B-rep body variants return `Held(UnsupportedOppositeSolid)` until their bounding predicates are independently implemented.
6. No persistent B-rep topology or transaction state is changed.

## Whole-cell evidence (stronger than point sampling)

G3 evaluates the extrema of each relevant **linear functional** over every directed boundary span of each UV region, including exact stationary angles of circular arcs; it does not facettize circles into chords.

For a certified AABB, the whole cell is **Inside** only when every world-coordinate interval lies within the corresponding box slab. It is **Outside** only when a single axis is entirely on the exterior side of a slab. Otherwise it remains `Held`.

For an analytic cylinder, the cell's axial and cylinder-frame transverse coordinate intervals produce safe bounds on radial distance. The entire cell is **Inside** when its axial bounds lie between the two caps and its conservative radial upper bound is within the radius. It is **Outside** when an axial slab or radial lower bound separates the entire region. Otherwise `Held`.

Floating-point tolerance cannot silently enlarge Inside/Outside regions. Comparisons preserve strict bounds, except when the G2 cell has a **physically validated cut edge** on the specific opposing AABB face and the discrepancy is within a representability-sized (`64 * f64::EPSILON * coordinate_scale`) rounding budget. Each such face-supported allowance is retained in `WholeSolidCellCertificate::AxisAlignedBox.cut_plane_support`. Unwitnessed near-boundary margins return `Held(UnresolvedNumericBoundary)` or `Held(NonUniformCellBounds)`.

After the whole-cell bound proof, a separate strictly interior point is chosen geometrically from the analytic directed boundary. The preexisting `BooleanRegionClassifier` independently checks that point against the **entire opposing BrepBody**. `Boundary`, missing interior witnesses, and contradictions never become Ready.

## Boolean selection policy (preview only)

| Operation | Source A, other volume | Source B, other volume |
| --- | --- | --- |
| Union | Keep Outside | Keep Outside |
| Intersection | Keep Inside | Keep Inside |
| Difference A \\ B | Keep Outside | Keep Inside, reverse orientation |

The `can_apply` flag stays **false** for all Ready classifications. Real G4 reconstruction still needs to pair analytic cut edges and solve material orientation consistently across all faces.

## Authored regression coverage

- True annulus and a box cutter, three regions labeled Outside / Inside / Outside.
- Union, Intersection, and Difference selection consistency.
- Missing otherwise-empty opposing source-face pairs: Held.
- Analytic cylinder crossed by a planar vertical face at two generators, with radial proof.
- Source bodies unchanged after classification/validation, tampered evidence rejected.
- Module-level mathematical checks for arc extrema, negatively unwrapped whole circles, second-operand Difference selection and translated coordinates.

## Validation notes

- Initial AI Workbench `cargo_test` failed before reading the manifest: `Permission denied (os error 13)`.
- A temporary validation copy of the project was created on the connected development device at the exact source commit and overlaid with the worktree source files. It is **not** the authoritative RustSolid repository.
- Standalone `rustfmt` syntax checking passed for G2/G3 modules and their integration tests after correcting a G2 negative-number comparison spacing bug.
- On the isolated, source-matched temporary clone (Rust 1.88, Linux Docker): `cargo test --test whole_solid_region_classification` passed **5/5** and `cargo test --test multi_trace_uv_arrangement` passed **8/8** after correcting its overly strict exact-tangency expectation. These are real compiled tests of the copied worktree, not tests on the AI Workbench runner.
- A precision probe found `1.7000000000000002` reconstructed for an AABB face at `1.7`. G3 now only permits the representability-sized adjustment if that box face is confirmed by the original analytic G2 physical cut.
- The complete regression suite is tracked separately; this document does not claim a full pass until its exit status is confirmed.

## G3 acceptance gates

1. Full source build with all v0.40 + v0.41 modules.
2. All G2/G3 unit and integration tests pass, and no preexisting regressions.
3. Nonconvex, near-tangent, ambiguous radial boundaries and incomplete source-face coverage fail closed.
4. Every G3 Ready witness is source-revalidated and `can_apply=false`.
5. Full material-surface classification across all source faces belongs to G4 scheduling; the current entrypoint certifies one bounded planar source-face arrangement at a time.

## Full regression verification (2026-10-15)

The temporary source-matched Linux Docker checkout completed `cargo test --quiet` with **exit code 0**: **530 passed, 0 failed, 62 test groups** (including G1/G2/G3 and preexisting RustSolid suites).

The verification checkout used the same `06ea004cb8e63bb6f53e6aeaeec83b5806258cfe` base revision plus the AI Workbench staged/unstaged Rust source and integration test files. It was an isolated copy; no commit/push/reset was made in the authoritative project.

AI Workbench's own `cargo_test` tool still reports a permissions error opening its mounted Cargo manifest, so this pass is specifically attributed to the copied checkout's real Rust 1.88 Linux test execution.

G4 topology reconstruction and G5 atomic Boolean commit remain future work. `can_apply` is still false on the G3 result.
