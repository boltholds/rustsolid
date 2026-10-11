# RustSolid v0.19 — first executable Boolean B-rep slice

This is a **strictly scoped planar Boolean executor**, not a general curved,
NURBS or Parasolid-compatible Boolean engine. It accepts **two independently
validated, exactly axis-aligned six-face, eight-vertex cuboids**. An arbitrary
concave extrusion, displaced prism, faceted cylinder, analytical cylinder, or
other B-rep is never silently substituted by its bounding box.

## Supported operations and result representation

`BooleanKernel::apply(&BrepBody, &BrepBody, BooleanOperation)`
and `BooleanKernel::apply_solids(&Solid, &Solid, BooleanOperation)` expose
`Union`, `Intersection`, and ordered `Difference` (A − B). A read-only
`BrepArena::boolean(&other, op)` delegates to the same implementation without
mutating either operand or either undo/redo stack.

`BooleanResult.outcome` is either:

- `BooleanOutcome::Solid(Solid)`: a **new independent canonical polyhedral
  body** with a fresh identity incarnation, manifold B-rep, face-to-mesh
  ownership, complete analytic mass, validated signed volume and model-space
  tolerance. All faces are exact planar quads (not an approximation).
- `BooleanOutcome::Empty`: the resulting 3D set has zero volume. There is
  no invented degenerate `Solid` object.

`BooleanResult.face_sources` contains one `BooleanFaceSource` per output
B-rep face, indicating `result_face`, which input `operand` supplied its
plane, the `source_face` and whether the source orientation is `reversed`.
For A − B, internal cut faces sourced from B have reversed orientation.
This is per-face lineage, **not** persistent topological naming through
arbitrary parametric regeneration. Newly minted `TopologyName` values
belong to the fresh `boolean/` feature namespace.

### Exact arrangement algorithm

The six bounding planes of two *certified* AABBs divide each Cartesian
coordinate into at most three nonzero intervals. Their product yields at most
27 open 3D cells, with exact source-coordinate boundaries. Each cell is
classified by membership of its entire interval in A and B. The Boolean
truth table selects which cells belong to the result. Only quadrilateral
patches separating a selected cell from a nonselected neighbor become output
faces. Adjacent patches share the same grid vertices, preventing mismatched
T-junctions. All retained patches have outward winding.

Volume and centroid use closed-form integrals over occupied cells, and area
sums only exposed patches; these quantities are **not** inferred from a
display mesh. The newly built B-rep is then independently checked against
the strict `Solid::validate` contract, including canonical coedge incidence,
fan connectivity, Euler characteristic, triangle face ownership, tessellation
closure, signed volume and area consistency.

### Unsupported topologies are errors

Current `Solid` validates **exactly one connected closed shell with genus 0
and one outer loop per face**. Therefore:

- A union yielding two separated solids cannot be represented yet.
- Removing an enclosed inner box creates a cavity and a second closed shell.
- Removing a through-slab may create disconnected components.
- Edge-only and point-only contact can create non-manifold unions.
- Through-holes produce genus 1, which current `Solid` rejects.

Such requests return `BooleanError::UnsupportedResultTopology`, never a
partial or superficially plausible body. The `Empty` outcome covers a
zero-volume 3D result; tangential contact does not create a solid.
Nearly coincident but unequal source plane coordinates closer than the
combined positional resolution return `BooleanError::AmbiguousBoundary`.
No implicit snapping or epsilon-driven reinterpretation is performed.

A curved operand or non-box polyhedral operand yields
`BooleanError::UnsupportedOperand`.

## Rust API example

```rust
use rustsolid::{BooleanKernel, BooleanOperation, BooleanOutcome, Point3, Solid};
let a = Solid::block(Point3{x:0.0,y:0.0,z:0.0},2.0,3.0,4.0)?;
let b = Solid::block(Point3{x:1.0,y:1.0,z:1.0},2.0,3.0,4.0)?;
let result = BooleanKernel::apply_solids(&a, &b, BooleanOperation::Intersection)
    .map_err(|e| rustsolid::GeometryError::InvalidEdit(e.to_string()))?;
if let BooleanOutcome::Solid(solid) = result.outcome {
    assert_eq!(solid.mass.volume,6.0);
    solid.validate()?;
}
```

## CLI geometry.v1

```json
{
  "operation": "boolean",
  "mode": "difference",
  "left": {"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
  "right":{"kind":"block","origin":[1,1,1],"width":2,"height":3,"depth":4}
}
```

A solid Boolean returns `ok:true`, a regular `solid` payload (vertices,
triangles and ownership, topology, mass, bbox), and `data` with
`kind:"boolean"`, `mode`, `state:"solid"`, a representation label, and
`face_sources`. Empty output returns `ok:true`, `solid:null`/omitted,
and `data.state:"empty"`. Unsupported cases return `ok:false`, a detailed
`error` string, and **no partial solid**.

## Deliberate architecture boundary

`TopologyDeltaPlan` from v0.16 stays preview-only, and
`BrepArena::require_executable_plan` still refuses all generalized plans.
This box-specific executor does not mutate analytic curves, parameter-space
trims, inner face loops or multiple shell topology. It is an acceptance
foundation for the next stages:

1. Multiple shells and inner-loop representation, with a validator capable
   of distinguishing cavities, components, orientation, and genus.
2. A common mutable geometry/trim registry for line, circle, cylinder and
   eventually NURBS carriers, reusing generational identities.
3. Face-region classification from the intersection graph, with verified
   inside/on/outside predicates and tangent/coincident classifications.
4. Atomic generalized surface splitting, trim loop reconstruction and
   provenance-aware Boolean commits into a stable modeling feature tree.
5. Failure injection, differential independent geometry checks, and
   tolerance/metamorphic regressions before general-purpose CAD use.

The first-party witnesses live in
[`tests/boolean_box.rs`](../tests/boolean_box.rs). GitHub Actions runs
this integration test explicitly as a required quality gate.
