# RustSolid v0.36 — dual analytic UV witnesses for oblique plane/cylinder trims

RustSolid v0.36 establishes a **read-only, mathematically defined
intersection-preview milestone** on top of the existing v0.35 source-certified
curved B-rep. It does not create a generalized curved Boolean executor.

## What is new

`IntersectionPcurve` represents one *analytic* surface-UV mapping of an
existing 3D `IntersectionCurve`, restricted to a finite parameter interval.
For each crossing `IntersectionGraph::arc`:

```rust
let [uv_first, uv_second] = arc.pcurves(&first_brep, &second_brep)?;
let position_3d = arc.evaluate(0.5)?;
let position_from_first = uv_first.reconstruct(0.5)?;
let position_from_second = uv_second.reconstruct(0.5)?;
```

For a plane/cylinder section, the supporting curve is

```text
C(u) = S_cyl(u, v(u))
v(u) = (n_plane · (O_plane - O_cyl)
        - r[n_plane·X_cyl cos(u) + n_plane·Y_cyl sin(u)])
        / (n_plane·A_cyl)
```

The cylinder pcurve is `(u, v(u))`, with **unwrapped** angular `u`; the
plane pcurve is the analytic inverse projection of `C(u)` into its
right-handed surface frame. Both use the *same* global intersection parameter
even when a bounded physical arc traverses `u=0/2π`. A face-side seam is
never treated as a new material boundary.

For plane/plane support intersections, the pcurve is the linear projection
into each support frame. For generator-line/cylinder intersections, the
cylindrical pcurve has fixed `u` and the line's varying axial `v`.

The graph validator samples endpoints and interior parameters of **both**
surface pcurves and refuses a graph when either reconstructed point disagrees
with the corresponding 3D curve beyond the source-model error budget. These
checks are numerical witnesses; they are not continuous interval proofs.

## Independent executable fixtures

`tests/oblique_plane_cylinder_preview.rs` has five test cases:

1. A horizontal finite box face intersecting a tilted finite cylinder:
   one endpoint on the physical face boundary `z=1`, the other on the
   cylinder's upper circular edge `v=2.4`; the arc crosses the parameter
   seam between them and the graph creates exactly one continuous physical
   arc with two endpoint nodes, not an extra seam node.
2. Both analytic pcurves reconstruct the same 3D crossing, with strictly
   increasing unwrapped cylinder U and checked operand swap.
3. Three independent local side probes per operand classify the curved
   crossing for Boolean Difference; the graph is a `ClassifiedPreview`,
   with meaningful orientation labels and no topology commit.
4. A physically *inclined* finite planar side face from an extruded,
   independently defined sloping XZ profile. Its bounded intersection with
   the tilted cylinder again crosses the seam, ends on two distinct real
   boundaries, and passes surface-UV/3D reconstruction and region preview.
5. A near-grazing axial boundary event within the numerical tolerance
   returns `TrimAmbiguity::UnresolvedBoundaryEvent`. The graph is
   `Ambiguous`, Boolean region is `Held`, and the symbolic topology delta
   is also `Held`, with no false physical edges.

Test command:

```sh
cargo test --test oblique_plane_cylinder_preview
cargo test --all-targets
cargo check --release --all-targets
```

All new mathematical fixtures use first-party equations and source geometry,
not proprietary/decompiled CAD algorithms or binary output.

## Precise implementation boundary

- The new `IntersectionPcurve` is a **preview evaluator**, not yet a
  persisted `Curve2` in a canonical coedge geometry catalog.
- Existing clipping supports a single outer polygon, a planar disk, or
  the original full periodic cylindrical side; multiple inner loops,
  arbitrary trimmed NURBS and conic/conic boundary intersection remain
  unsupported and are *not* silently approximated by mesh polygons.
- `TopologyDeltaPlan` records the pending physical edge splits and face
  work orders, but `can_apply` remains false.
- `BrepModel` is still a read-only shared interchange representation,
  and `BrepArena` cannot commit general curved Boolean transactions.
  The dedicated v0.35 connected triple-cylinder B-rep remains separately
  executable and fully source-bound.
- Tangencies, ambiguous boundary events and missing numerical witnesses
  remain `Held` rather than being treated as valid material cuts.

Next independent milestone: turn preview pcurve witnesses into
coedge-owned `Curve2` support for a deliberately bounded executor, add
full face-cell classification and then commit atomic topology only after
manifold, orientation, source-mass, watertight display and Undo/Redo
acceptance gates pass.

## Verification (2026-10-14)

The isolated WSL test checkout was cloned at the same parent revision
`06ea004cb8e63bb6f53e6aeaeec83b5806258cfe`, and the new Rust source,
new integration test and updated Cargo manifest were copied from this
AI Workbench worktree for validation.

- `cargo test --test oblique_plane_cylinder_preview -- --nocapture`:
  **5 passed, 0 failed**.
- `cargo test --all-targets`: **passed**.
- `cargo check --release --all-targets`: **passed**.
- Existing preexisting unused-import/dead-code compiler warnings remain;
  no runtime execution of any proprietary DLL took place.

Evidence from fixtures establishes behavior only for their tested analytic
families and numerical scales. Production CAD parity and complete
continuous-surface bounds are not claimed.
