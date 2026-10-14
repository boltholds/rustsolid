# RustSolid multi-face region preview — v0.39 basis, v0.40 extension

This component gathers **source-owned local intersection and Boolean
region evidence**. It does not perform general B-rep mutations.

## API

```rust
let full = MultiFaceRegionPlan::plan(
    &left, &right, BooleanOperation::Difference,
    GeometryTolerance::default(),
)?;

// Restrict analysis to one genuine face pair without inventing
// global completeness or ignoring unrelated source geometry.
let selected = MultiFaceRegionPlan::plan_selected(
    &left, &right, BooleanOperation::Difference,
    GeometryTolerance::default(),
    &[(FaceId(1), FaceId(2))],
)?;
selected.validate_against(&left, &right, GeometryTolerance::default())?;
assert!(!selected.can_apply);
```

The full `plan` evaluates all pairs of faces (up to 64 source faces per
operand and 512 pairs). The selected variant evaluates only an
explicitly named, validated subset. The selected variant **does not**
claim knowledge of the other source faces.

For every classified pair, v0.40 retains:

- `IntersectionGraph`: first-party bounded analytic 3D curves and
  their endpoint/face provenance;
- `GraphRegionPlan`: three local side samples per arc per operand,
  not necessarily complete material cells;
- `CoedgeReconciliationPlan`: one unique physical 3D carrier per
  source graph arc with two opposite independently evaluable
  owner-face UV pcurves;
- `UvFaceRegionPlan` when one participating *planar* source face has
  a strictly convex straight outer contour and supported native
  circular inner loops, crossed by one transverse straight line.

In particular, **native circular holes are no longer rejected
globally**. Their actual line/circle intersections and source loops
participate in clipping and face-cell construction. Unsupported curved
inner contours, nonconvex multi-cell arrangements, and other
uncertain trimming configurations remain held.

## Strict all-or-nothing Held

Any unsupported face pair, numerical ambiguity, tangent or coincident
boundary, nonseparating source region, unresolved pcurve ownership, or
failed required hole-cell construction places the entire selected plan
into a typed `Held` with zero retained partial geometric witnesses.
No partial mesh or body is ever committed.

The same rules apply to the full Cartesian face-pair plan. An
`Empty` result has no proposed cut. A `ClassifiedPreview` is still
**not** permission to execute a generalized Boolean.

The source B-rep stays immutable. All stored values can be checked by
`validate_against`, which reconstructs the selected witnesses from
the original operands and requested tolerance; forged counts or
pcurve/edge witnesses are rejected.

## Tests / limitations

`tests/multi_face_region_preview.rs` verifies a disjoint pair, at
least four real face-pair intersections, a real native through-hole
source, near-tangency, and a selected annular cap cut producing
exactly two analytic arcs, two physical edge candidates and two
complete planar UV cells.

For the full mathematical contract, area integrals, source coedge
fragment spans, seam/periodic treatment and the boundary between
face-cell readiness and general Boolean proof, see
[planar UV cells and coedge reconciliation](planar-uv-cells-and-coedge-reconciliation.md).

**Not yet supported**: all global face-cell arrangements on arbitrary
NURBS, curved inner loops, multiple mutually intersecting cuts per
face, boundary coincidence, disconnected nonconvex source regions,
generalized source-independent B-rep topology execution and global
material-region volume integration. `can_apply` is always false.
