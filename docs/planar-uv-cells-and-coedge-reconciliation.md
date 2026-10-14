# RustSolid v0.40 — circular-hole UV cells and cross-face edge reconciliation

RustSolid v0.40 adds a deliberately bounded **geometric/topological
preview** for one transverse *straight* cut across a convex planar
B-rep face containing disjoint native circular inner loops. This is a
reusable mathematical step toward generalized curved B-rep Boolean:
it is **not** a Boolean executor and never mutates an input body.

## 1. Real planar inner-loop clipping

`TrimmedIntersectionEngine` now recognizes a planar cap whose B-rep
`face.loops` contains one outer polygon or disk and up to 64 inner
`Curve2::Circle` loops. Every hole is checked for strict containment
in its outer trim boundary and strict separation from every other
hole using *exact circle geometry* (not faceted display).

Membership is derived from the actual oriented UV domain:

```text
Inside(face) = Inside(outer) ∧ not Inside(any inner hole)
Boundary(face) = Boundary(outer) or Boundary(any inner hole)
```

A line-carrier's original boundary events are augmented with the exact
line/circle quadratic roots for *each* hole. This produces genuine
truncated 3D intersection arcs on either side of holes without
inventing geometry inside material-free openings. Near tangency
continues to produce `TrimmedIntersection::Ambiguous` and cannot
be mistaken for a transverse Boolean edge. Unsupported curved
inner-loop carriers/conic-conic hole events remain explicitly
unsupported.

## 2. Whole planar face-cell arrangement (bounded analytic family)

`UvFaceRegionPlan::build` uses the already source-validated
`IntersectionGraph` and a chosen source operand to construct
**complete closed oriented UV loops** for both regions on the
two sides of a single transverse supporting *line*.

Within the supported family:
- outer loop: strictly convex finite polygon, with native linear
  coedges and two transverse crossings of the cut line;
- inner loops: up to 64 strictly separated native circles; some may be
  entirely on one side, while others are cut into true circular arcs;
- each cell has one positively oriented outer loop and zero or more
  negatively oriented *inner* loops;
- source fragments retain `EdgeId`, `CoedgeId`, `LoopId` and original
  carrier parameter endpoints; circle parameters may unwrap across
  `2π` without creating a physical seam vertex;
- every retained intersection arc has one new cut fragment in each
  cell with opposite orientation and the **same** source `GraphArcId`.

Exact-in-form Green integrals are used for each straight or circular
fragment, not polygonal approximations of holes:

```text
A_line(a,b) = (a.u*b.v - a.v*b.u)/2

2*A_circle_arc = r*Cx*(sin(end)-sin(start))
               - r*Cy*(cos(end)-cos(start))
               + r*r*(end-start)
```

The signed area of every resulting cycle is computed directly from
its actual physical carriers. All cell areas are positive, and their
sum must reproduce the independently computed source-face area:

```text
A_positive + A_negative
    = |A_outer_polygon| - π * Σ(radius_inner²)
```

Unresolved cut-through vertices, tangent circles, coincident edges,
nonconvex/disconnected domains, nonmanifold UV incidences and
mismatched source area are *Held* (or explicitly rejected), with
no published partial cells.

This is a complete 2D planar arrangement for **the declared family and
one cutting line**, not full general UV arrangements for several
intersecting curves, NURBS pcurves or multi-face cell complexes.

## 3. Two truly source-owned coedges for a physical 3D intersection

`CoedgeReconciliationPlan` assigns **no permanent topology IDs**.
It reconstructs each validated `GraphArc` as one original analytic
3D curve plus two independently evaluable `Curve2` pcurve carriers:

```text
IntersectionGraph.arc
  ↳ 3D Curve3::Line / Curve3::Harmonic  (one physical carrier)
  ↳ first face:  Curve2::Line/Harmonic/CylinderSection (forward)
  ↳ second face: Curve2::Line/Harmonic/CylinderSection (reversed)
```

For actual physical arcs, the two uses have opposing 3D orientation,
the same geometric endpoints, and pointwise-equal 3D reconstructions.
The natural unwrapped cylinder U is preserved across `0/2π`. A
periodic graph anchor on a **closed** ellipse is explicitly
`virtual_periodic_anchor=true`: it must not become a new material
vertex solely due to chart periodicity.

`UvFaceRegionPlan::certify_physical_twins` now cross-checks the
**intra-face pair** of opposite UV-cell cut fragments against these
**cross-face** pcurve owners. For each cut arc:
- the physical 3D endpoints must coincide with the source-face UV
  reconstruction (either allowable parameter direction);
- each cell contains **exactly one** use of that graph arc;
- the two uses have opposite endpoints/orientation;
- each UV span is aligned with the reconciled source-face pcurve.

Any mismatch prevents readiness of the multi-face preview.

## 4. Multi-face preview update

`MultiFaceRegionPlan::plan_selected` can now isolate one bounded set
of actual source-face pairs. For each supported pair it retains:

```text
    source face pair
      -> 3D+UV IntersectionGraph
      -> pair-local BooleanRegionClassifier samples
      -> CoedgeReconciliationPlan
      -> UvFaceRegionPlan (when a planar face admits certified cells)
```

The original `plan()` still checks the full Cartesian face-pair set.
Either entry point is **all or nothing**. Any unsupported geometry,
uncertain contact, hole domain beyond the supported circular form,
failed UV-cell construction where holes require it, or mismatch in
the physical edge reconciliation yields `Held` with zero retained
partial pairs/arcs/cells.

`MultiFaceRegionPlan::can_apply` remains **false** even when every
local face pair is classified and every UV cell passes area and twin
validation. Local side samples are not a proof that **all** 3D
material regions across all source faces have been assigned
consistently.

## 5. Executable first-party numerical fixtures

The v0.40 integration tests use original independently constructed
RustSolid sources, not a proprietary kernel.

- One through-bored box cap vs a real perpendicular face: the
  intersection line is clipped into **two** physical arcs by the
  circular hole. Both new endpoints are owned by the original
  circular edge and neither point lies inside the empty opening.
- A rectangular cap without holes: two closed cells of equal area.
- A square cap with one centred circle: two closed cells of area
  `18 - π/2` each.
- A circle entirely to one side: the containing cell retains a
  single negative-area native circular loop; no fake seam split.
- **Three disjoint through-holes:** one intersected circle is
  partitioned into exact arcs and two untouched holes retain their
  complete circular inner loops. The two cell areas match independent
  `πr²` arithmetic.
- A crossing pair of planar source faces: one true 3D edge per arc,
  with two opposite source-owned UV carriers and a verified
  per-cell cut-twin ledger.
- A transverse closed elliptic plane/cylinder curve: a single
  physical elliptic edge, two owned analytic pcurves, one *virtual*
  periodic anchor.
- An oblique section crossing the cylinder's periodic seam: one
  physical curve, no false seam vertex.
- Near-tangent hole events, unsupported nonconvex faces, disjoint
  supports and tampered source-bound plans are rejected/held.

The tests confirm geometry and topology of their finite analytic
fixtures. They do not imply formal interval certification or safety
for arbitrary floating-point scales.

## 6. Why generalized Boolean commit stays disabled

No general mutable `BrepModel` operator is activated by this stage.
A complete 3D Boolean still requires, at minimum:

1. global arrangement of **all** intersections of multiple curves
   on every affected face (including NURBS, coincident regions,
   holes, and open/closed shells);
2. globally consistent ownership/inside-outside classification of
   every bounded 3D material cell, not just arc-side probes;
3. one global physical `EdgeId` per unique intersection branch,
   paired coedges, correct face loop rewiring and source provenance;
4. atomic candidate construction, genus/Euler and manifold validation,
   a watertight face-owned mesh and independently integrated mass
   and centroid, followed by Undo/Redo with no partial publish.

Existing v0.35–v0.39 source-certified Boolean/imprint transactions
remain executable for **their** restricted source families. This
milestone is a generalizable geometric contract, not another
hidden special-case arbitrary Boolean.

## Commands

```sh
cargo test --test planar_uv_hole_clipping
cargo test --test planar_uv_arrangement
cargo test --test coedge_reconciliation
cargo test --test multi_face_region_preview
cargo test --all-targets
cargo check --release --all-targets
```

No decompiled Parasolid/Siemens algorithms, machine code, opaque
private data structures or proprietary SDKs are included.
