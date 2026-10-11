# RustSolid v0.14 — bounded B-rep face intersections

This is the first **bounded geometric intersection** slice. Unlike `SurfaceIntersectionEngine`, which intersects infinite supporting surfaces, `TrimmedIntersectionEngine::intersect_faces` accepts a `BrepModel` and a `FaceId` for **each** operand. It clips analytic support branches against the actual UV trim data, without substituting render facets for geometry.

## Supported first slice

- One simple, closed **planar polygon** with straight `Curve2::Line` coedges; supports both convex and concave outlines.
- One **circular disk** bounded by a single closed `Curve2::Circle` coedge.
- One **untrimmed periodic cylindrical side**, bounded axially and represented with two distinct seam p-curves. The seam is *not* a physical face boundary.
- Plane/plane support intersections yield finite line segments or isolated vertex contacts. A line running exactly along a planar face edge is classified as `BoundaryOverlap`.
- Plane/cylinder support intersections yield finite angular arcs or generator segments. Closed circular intersections, partial arcs, parallel-axis secants and tangent generators are supported.
- Boundary events are determined analytically: 2D line/segment intersections and line/circle roots on planar trims; 3D segment/cylinder quadratic roots; roots of the sinusoidal axial boundary equation on periodic cylindrical trims. This is **event-based clipping**, not resolution-dependent sampling.
- A clipped result contains the original analytic `IntersectionCurve`, a bounded parameter interval, geometric endpoints, contact type and a sample-based two-surface residual certificate.

### Deliberate safety restrictions

- The existing `BrepModel` currently supports one closed two-use manifold shell. The clipping layer handles **one outer loop** per selected face; holes, multiple loops and arbitrary UV domains are rejected rather than silently treated as rectangular bounds.
- Coincident supporting surfaces return `Ambiguous(CoincidentSupports)`, because a complete face/face overlap operation is separate from line clipping.
- Near-tangent or near-parallel supporting surfaces retain the existing explicit ambiguity classification. Numerically unresolved trim-boundary events produce `Ambiguous(UnresolvedBoundaryEvent)` instead of speculative cut edges.
- The only supported conic/circular-disk pairing is the provably concentric horizontal cylinder section. Oblique conic/disk boundaries need an independent conic/conic clipping solver and are currently rejected.
- `IntersectionCertificate` reports sampled residuals on the infinite supporting surfaces; it does **not** prove the absence of an unobserved boundary crossing or certify a Boolean solid.
- Numerical positions and event roots are `f64` approximations. Robust orientation signs are used for polygon containment, but interval-certified intersection construction is future work.

## API

```rust
let a = BrepModel::from_polyhedral(&block, ModelUnits::default())?;
let b = BrepModel::from_cylinder(&cylinder)?;
let trimmed = TrimmedIntersectionEngine::new(GeometryTolerance::default())?
    .intersect_faces(&a, FaceId(0), &b, FaceId(2))?;
for branch in trimmed.branches() {
    let point = branch.evaluate(0.5)?;
    // This remains a mathematical intersection carrier, not a Boolean edge.
}
```

The stateless `geometry.v1` CLI supports `intersect_faces`, taking `{ "source": ..., "face_id": N }` operands and optional `samples` and `tolerance` fields. The response is a read-only `data` object with `kind`, branch parameters, geometry samples, boundary contacts and a residual certificate. Existing requests and responses are unchanged.

## Independent regression contracts

`tests/fixtures/trimmed_intersections.v1.json` permanently pins 12 face/face scenarios. `tests/trimmed_intersections.rs` checks analytically known endpoints, circle and generator equations, angular arc measure, tangent/boundary classification, symmetry under operand swap, integrity rejection, and separation of infinite supports from bounded trims. `cargo test --test trimmed_intersections` is a required CI step. No Parasolid/OCCT result is treated as a mathematical oracle.

## Next milestones

1. Face trimming on oblique conic sections by circular/elliptic UV edges, along with multi-loop and hole classification.
2. Stable interval-certified event isolation and explicit per-vertex/per-edge error bounds.
3. Mutable shared B-rep topology with generational slots, transaction-safe subdivision and provenance.
4. Assemble Boolean cut loops, classify inside/outside volumes, and support sewing/validation — **not included** in this slice.
