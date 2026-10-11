# RustSolid v0.13 — analytical surface intersections and independent contracts

This phase adds `SurfaceIntersectionEngine` for **untrimmed** supporting surfaces:
plane/plane and plane/right circular cylinder. An intersection curve is
mathematical geometry; it is not yet a trimmed B-rep edge, a Boolean result or
a proof that two bounded faces intersect.

## Classification

- Plane/plane: empty, coincident, one crossing line, or unresolved near-parallel.
- Plane/cylinder with nonzero `n · axis`: a closed analytic conic section,
  parameterized by cylinder angle `u` with axial coordinate derived exactly
  from the plane equation in the ideal real-arithmetic model.
- Plane/cylinder parallel to the cylinder axis: zero, one tangent or two
  straight generators, determined by the intersection of a line and circle
  in the cylinder transverse plane.
- A separation comparable with engineering tolerance around tangent status,
  or an angle comparable with the angular tolerance around axis parallelism,
  is classified `Ambiguous`; no artificial Boolean topology is constructed.

For a cylinder frame `(O,X,Y,A)`, radius `r`, plane `(P,n)`:

```
S(u,v) = O + r(X cos u + Y sin u) + A v
n · (S(u,v) - P) = 0
v(u) = [n · (P - O) - r(n · X cos u + n · Y sin u)] / (n · A)
```

The last equation is only used when the denominator is comfortably above the
angular conditioning boundary. If `n · A = 0`, the circle equation
`r[(n · X) cos u + (n · Y) sin u] = n · (P - O)` determines the generators.

## Independent acceptance contracts

`tests/fixtures/intersections.v1.json` pins 14 analytic, metamorphic and
boundary witnesses. `tests/surface_intersection_contracts.rs` enforces
coverage for BOTH supported pairs; it rejects missing categories, duplicates
and forged classifications. Additional focused tests are in
`tests/intersections.rs`.

Each constructed branch is evaluated at multiple parameters. The
`IntersectionCertificate` stores the measured maximum residual against each
implicit supporting surface and the applied numerical tolerance. It **does
not** prove completeness, nonintersection, or exactness of all points.

Tests use known equations (line of two coordinate planes, circular and
oblique cylinder sections), operand swap, translation, cylinder reparameterization
and deliberate modifications to the expected or computed geometry. Neither
Parasolid nor OpenCascade is used as an oracle.

## Read-only CLI

```json
{
  "operation": "intersect_surfaces",
  "first": {
    "kind": "plane", "origin": [0,3,0],
    "x_axis": [1,0,0], "y_axis": [0,0,-1]
  },
  "second": {
    "kind": "cylinder", "origin": [0,0,0], "radius": 2,
    "x_axis": [1,0,0], "y_axis": [0,0,-1]
  },
  "samples": 9
}
```

The response uses the existing `geometry.v1` envelope and returns a `data`
object with result kind, branch sample coordinates and a residual certificate.
Requests with unsupported surface types, invalid frames or excessive samples
fail without mutating a stored model. Existing CLI operations remain intact.

## Follow-up gates

1. Develop bounded face/face clipping by UV trim loops and distinguish open
   crossing segments, boundary contacts, overlap, and true disjointness.
2. Extend exact predicates to intersection construction and separately
   track *certified* intervals versus f64 sampling.
3. Integrate intersection curves into a writable generational B-rep arena
   with Command/Memento and face provenance.
4. Add more surface pairs only with independent mathematical oracles,
   metamorphic properties and near-degeneracy classifications.
