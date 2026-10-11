# RustSolid v0.25 — analytic top-entry blind cylindrical Boolean

RustSolid now executes a second independently validated curved subtraction:
**a vertical, finite-depth cylindrical bore opening on the upper plane of
a certified axis-aligned box**. Unlike the v0.21 full-height through-hole,
this result contains an actual *new planar circular floor*, plus a true
cylindrical side wall and a circular inner trim loop on the existing
top face. No triangles are used to define the Boolean boundary.

## Exact source admissibility

The executable Boolean is `BooleanKernel::apply(&box,&cylinder,
BooleanOperation::Difference)`, with:

- Left: a fully validated **six-face, eight-vertex axis-aligned box**. A
  generic convex or concave shape is not silently replaced by its AABB.
- Right: a validated `CylindricalBrep` with a strictly **+Y** axis
  (machine-representability, not angular-model-tolerance snapping), its
  mathematical upper cap on the box's top (+Y) face, and lower cap
  **strictly inside** the material above the box bottom.
- The entire circular cross-section is strictly inside all four box
  vertical walls, with local engineering-tolerance clearance.
- The cylinder depth is finite and strictly between zero and the box
  vertical height. Cylinder side surfaces intersect no other box sides.
- The real top box face and cutter cylindrical side must produce **one
  complete, source-validated, closed, transverse analytic
  `IntersectionGraph` arc** at the opening boundary.

The last condition is verified through `IntersectionGraphBuilder`,
which rechecks its mathematical curve against both actual input B-reps.
The closed curve is a physical top-face trim and source carrier, not
a circle inferred from a triangulated cylinder. Ambiguous, coincident,
partial, fragmented and tangent contact graphs fail closed.

**Plane coincidence has a stronger rule than physical tolerance snapping:**
the cutter upper cap must meet the top source plane to the floating-point
representability budget; an unequal near-coincident plane is not silently
moved. Likewise the existing v0.21 through-bore no longer tolerates a
physically distinct near-bottom cutter being reinterpreted as a
full-height through cut.

General curved bodies, rotated/oblique bore axes, cylinders contacting
the four lateral box walls, multiple overlapping bores, disjoint shells,
and arbitrary intersections remain unsupported. Unsupported operands
return errors, not faceted Boolean stand-ins.

## Analytic region and shell topology

The result is `BooleanOutcome::BlindHole(BlindHoleBrep)`. A boxed
top-entry blind cut has:

| Canonical entity | Count |
| --- | ---: |
| Vertices | 10 |
| Edges | 15 |
| Coedges | 30 |
| Faces | 8 |
| Loops | 9 |
| Inner face loops | 1 |
| Closed shells | 1 |
| Circular edges | 2 |
| Cylindrical faces | 1 |
| Periodic seam edges | 1 |
| Euler characteristic | 2 |
| Genus | 0 |

The original bottom box face is **not** punctured. The original top
face now has a circular inner loop, shared with the top circle edge of
the newly introduced cylindrical side. Its lower circle is shared with
a distinct planar disk face at the bore floor. Its physical periodic
seam retains separate UV trims for the two directed coedge uses.
All faces, loops, coedges and 3D/2D geometric carrier references pass
the shared manifold and UV/3D validation.

The model additionally records an explicit `FaceSense` on each
`BrepFace`, independent of its `Surface3` UV frame. A standalone
cylinder's radial support has normal pointing **outward**, whereas
a subtractive cylindrical hole has a face normal pointing **into the
empty bore**: `FaceSense::Reversed`. The top-entry blind bore's
new circular floor is physically **+Y**, although its canonical UV
circle shares a **−Y plane support parameterization**; its face also
has `FaceSense::Reversed`. The derived
`BrepModel::oriented_face_normal(face,u,v)` evaluates the actual
material-oriented normal. The common B-rep validator also checks that
faces produced by each existing canonical source have the expected
material sense; changes that invert only the face orientation fail closed. This is represented in the normal B-rep query
schema as `faces { sense }` and in `inspect_brep` metadata.

The finite-depth shape is **genus zero**, since only one cap has an
inner loop. The Euler count follows the cell-complex-with-holes
invariant `10-15+8-1=2`.

## Exact mass vs. display faceting

Let `V0`, `S0` and `C0` be the original box's analytic volume, area
and centroid, `r` the bore radius, `d` its depth and `topY` the box's
upper Y coordinate. Let the removed cylindrical volume be
`Vr=πr²d` and its centroid `Cr=(cx, topY-d/2, cz)`.

```
V = V0 - Vr
S = S0 + 2πrd
C = C0 + (C0 - Cr) * Vr / V
```

The bore's removed circular top-cap area is exactly replaced by its
circular planar floor. Consequently their areas cancel and only the
new cylindrical lateral area contributes to `S`. Mass calculations
use box-relative subtraction moments, reducing cancellation for large
absolute world positions.

`BlindHoleBrep::tessellate(12..=512)` produces a **separate viewport
representation**. It triangulates the upper rectangular cap with one
circular inner loop through the ISC-licensed earcutr display triangulator,
retains all original bottom/side box triangles, creates a topologically
shared floor disk fan and a cylinder wall from the same circular vertex
indices, and independently certifies:

- One two-use, oppositely directed incidence per mesh edge.
- Valid face lineage and no degenerate facets.
- Plane and cylinder support distance and **material** normal signs.
- The planar top annulus area against a separate polygon area oracle.
- The polyhedral signed tetrahedral volume against
  `Vfacet=V0-Apolygon*d`.
- The polyhedral surface area against
  `Sfacet=S0+Ppolygon*d`, and the expected convergence to analytic mass
  as the number of circle segments increases.

Facet volume and surface area remain separate from analytic mass. A
near/exact circle radius, topological genus and 3D carrier parameters
never depend on triangle count.

## Parametric edits and arena lifetime

`BlindHoleBrep` owns a certified original box, parameters
`BlindHoleParameters {center_x,center_z,radius,depth}`,
canonical analytic B-rep, monotonic revision, and bounded Undo/Redo
parameter history.

`BlindHoleTransaction` stages changes to radius, center and depth in
one isolated candidate. The entire prospective B-rep and analytic mass
are rebuilt and validated before publication. A poisoned stage,
unsupported depth or lateral contact leaves body, revision and history
untouched. Undo and Redo similarly reconstruct validated candidates
before exposing them.

`BrepArena::from_blind_hole()` and
`BrepArena::update_blind_hole()` integrate the feature into the
existing sparse generational identity registry. A completely valid
candidate and reconciled registry are published together: surviving
named faces and circle edges keep handles across radius/depth edits.
A failed update never invalidates previously valid topology handles
or source snapshots.

`BrepCarrierEditor::from_body(&BrepBody::BlindHole(...))` also supports
source-preserving analytic curve/UV chart reparameterization for this
shape without changing its exact physical mass. Generic structural
operations such as adding a second blind bore to the same body are
**not** yet executable through `BrepStructuralEditor`.

## Rust and JSON APIs

Rust:

```rust
use rustsolid::{
    BlindHoleBrep, BlindHoleParameters, BrepBody, BooleanKernel,
    BooleanOperation, BooleanOutcome, Point3, Solid, CylindricalBrep,
    GeometryTolerance
};
let box_=Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let cutter=CylindricalBrep::upright(
    Point3{x:0.,y:1.,z:0.},1.,2.,GeometryTolerance::default())?;
let result=BooleanKernel::apply(&BrepBody::from(box_),
    &BrepBody::from(cutter),BooleanOperation::Difference)?;
if let BooleanOutcome::BlindHole(blind)=result.outcome {
    blind.validate()?;
    assert_eq!(blind.model().summary().genus,0);
    let display=blind.tessellate(64)?;
    assert!(display.faceting_volume_error>0.);
}
```

Existing `geometry.v1` `boolean` command:

```json
{
  "operation":"boolean","mode":"difference",
  "left":{"kind":"block","origin":[-3,-1,-3],
          "width":6,"height":4,"depth":6},
  "right":{"kind":"analytic_cylinder","origin":[0,1,0],
           "radius":1,"height":2}
}
```

Successful results return `state:"solid"`,
`representation:"analytic_cylindrical_blind_hole"`,
`source_opening_trace:{state:"ready",arcs:1,closed:true,...}`,
`opening_face_id:1`, `floor_face_id:7`, `side_face_id:6`,
`genus:0`, `topology_committed:true`, verified display faceting
metrics, and per-face Boolean source lineage. The normal `solid`
payload contains derived triangles and exact mass/topological counts.

The independent `inspect_brep`, `geometry_query` and
`carrier_transaction` source `kind:"analytic_blind_hole"`
takes a box `origin,width,height,depth`, a circle
`center:[cx,cz],radius` and separate `bore_depth`.
This avoids confusing the block Z-depth with the bore Y-depth.

## Independent acceptance and future scope

The **22 independent** `tests/blind_hole.rs` regressions cover exact Euler/genus/loop ownership,
analytical mass and centroid, source-bound opening traces and tampering
protection, true cavity material classification, closed display mesh
orientation, convergence, translated-world robustness, atomic radius/
center/depth editing, arena Undo/Redo and topology handle preservation,
and strict refusal of ambiguous or unsupported Boolean cuts. The JSON
CLI and mandatory GitHub Actions independently repeat critical
through/blind/no-result gates.

The original analytic through-hole, exact axis-aligned box CSG and
genus-N multi-through-hole structural editor are preserved.

**Still not implemented:** A generalized mutable face-region arrangement
for arbitrary plane/cylinder crossings, oblique/blind intersections,
multiple cavities, arbitrarily oriented ruled surfaces, multi-shell
topology, analytic NURBS, fillets, and a general
`TopologyDeltaPlan::can_apply` executor. The new verified operation
is a strict, useful first instance of **finite-depth face-region
reconstruction**, not a universal Boolean engine.
