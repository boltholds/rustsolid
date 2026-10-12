# RustSolid v0.27 — executable crossing analytic circle-arc arrangement

The first *intersecting* curved Boolean region implemented by RustSolid is
the subtraction of the **union of two partially overlapping, top-entry,
equal-depth upright cylinders** from an exact certified axis-aligned box.

This is a genuine analytical B-rep with crossing-generated vertices, open
circular arcs, shared vertical edges, UV trim loops, a common floor face,
per-surface material orientation and independently integrated mass. It is
not based on triangle clipping and is not the general curved
`TopologyDeltaPlan` executor.

## Source proof and rejection policy

Both supplied `BlindHoleParameters` undergo the v0.25 independent
top-entry cylinder and source `IntersectionGraph` validation. The input
block must have actual axis-aligned planar faces; treating an arbitrary
polyhedron's AABB as occupied material is prohibited.

The two planar circular footprints must have exactly **two transverse
crossings**:

`|r1-r2| + tolerance < center_distance < r1+r2 - tolerance`.

Disjoint circles are still supported by the v0.26
[`CompositeBlindBrep`](composite-blind.md) with separate floor faces.
Contained, coincident, tangent, nearly tangent, improperly positioned,
insufficiently separated from box walls and out-of-bounds circles
**fail closed** in the crossed arrangement builder.

Both cylinders must have the **same** finite depth, with equality
checked to a floating-point representability budget. Differing depths
would create a stepped floor and horizontal split curves; this more
general arrangement remains future work.

Only +Y cylinders entering the **top** of the box are admitted. Oblique
cuts, three-way crossing intersections and independent components are
not silently approximated.

## Building a real analytic arrangement

`CrossedDiskRegion::from_box` solves the two circle equations in the
box's XZ chart. Let `D` be the distance between centers, `e` their
unit separation vector, and `x1=(D²+r1²-r2²)/(2D)`. The two
intersections are `C1+x1*e ± h*perp(e)`, with
`h²=r1²-x1²`. Both points are checked against the two mathematical
source circles. The near-tangent case `h≈0` is rejected.

For each circle the planner considers the two angular paths between
crossings and retains **only the oriented arc on the exterior of the
other disk**. Midpoints and additional fractions certify its position
outside the other circle. The two surviving arcs must form one closed
oriented loop. Their native `Curve3::Circle` parameters can cross
`2π` without splitting a mathematical arc at an artificial seam.

The resulting B-rep uses exactly:

| Property | Count |
| --- | ---: |
| Vertices | 12 |
| Edges | 18 |
| Coedges | 36 |
| Faces | 9 |
| Loops | 10 |
| Closed shells | 1 |
| Circular arc edges | 4 |
| Shared vertical cylinder crossing edges | 2 |
| Cylinder faces | 2 |
| Planar union-floor faces | 1 |
| Top inner-loop contours | 1 |
| Euler characteristic | 2 |
| Genus | 0 |

The top planar box face obtains a *single inner loop consisting of
two actual circular arc coedges*. The independent upward-facing floor
has a single outer loop consisting of the reverse uses of the two lower
circle arc edges. Each cylindrical side is a truly trimmed rectangular
analytic parameter patch with a bottom arc, vertical crossing edge,
reversed top arc, and opposite vertical crossing edge. Both shared
vertical edges have **two opposing coedges on different cylindrical
surfaces**, not fictitious periodic same-face seams. No cylindrical
patch uses a closed full circle outside its true visible domain.

The two material-facing cylindrical faces and the floor use the
required reversed `FaceSense`; source box faces stay forward.
`BrepModel::validate` independently checks coedge twins, next/prev
cycles, vertex coincidence, UV-to-3D correspondence, one connected
closed shell, orientation and Euler characteristic.

`CrossedDiskRegion` retains the two crossing witnesses, exposed arc
domains and source indices, exact disk-union area, exposed perimeter,
removed-area moments and cut depth. It cannot be confused with
the independent preview-only `TopologyDeltaPlan`.

## Mathematical mass

The circle intersection lens is split at its common intersection
chord into two exact circular segments. For radius `r` and signed
chord coordinate `x`, each segment has area:

`Aseg = r² acos(x/r) - x sqrt(r²-x²)`

The union area is
`Aunion=πr1²+πr2²-Aseg1-Aseg2`, and its perimeter is the sum
of the two surviving circular arc lengths.

For the common depth `d`, original box volume `V0` and surface
area `S0`:

`V = V0 - d*Aunion`

`S = S0 + d*(arc_length1+arc_length2)`

The removed cylinder union has constant XZ cross-section and its axial
centroid is `y = box_top_y-d/2`. XZ moments are derived analytically
from the circle disks minus both chord-bounded lens segments, using
box-relative coordinates to reduce cancellation. No tessellation
resolution influences these mass values.

## Atomic committed source edits

`IntersectingBlindBrep` owns the original certified box, its two
`BlindHoleParameters`, the separately reconstructed
`CrossedDiskRegion`, a complete canonical `BrepModel`, its revision
and a bounded Undo/Redo history.

`IntersectingBlindTransaction::replace(index, parameters)` stages new
center/radius/depth values on an isolated copy. A commit rebuilds and
validates the *entire* crossed-circle arrangement and exact analytic
mass before replacing one published state. A malformed stage poisons
the transaction even if its error was ignored. A late unsupported
topological transition (e.g., tangency or different depths) rolls back
all staged changes, revision and history. Undo/Redo restore validated
candidates and increase revision monotonically.

No source `Solid` or previously constructed B-rep is mutated by the
builder. The common `geometry_query` API can inspect committed
`BrepModel` instances, including exact curved topology and face sense.

## Common kernel ownership and material classification

`BrepBody::IntersectingBlind(IntersectingBlindBrep)` implements the normal
`BrepSource` contract for canonical analytical queries and derived
viewport meshing. The three-way `BooleanRegionClassifier` distinguishes
material from the union of the two missing disk columns, the common
floor plane, both visible cylindrical walls, the top opening, and the
unchanged box envelope. Intersections hidden within the *other* removed
disk are **not** misidentified as exterior cylinder faces.

`BrepArena::from_intersecting_blind` creates a sparse
owner-generational topology registry for these nine faces and the
split crossing edges. `update_intersecting_blind([params_a,params_b])`
isolates the entire reassembly, commits valid candidates together with
identity registry reconciliation and preserves semantic handles for
surviving face/edge owners. Failed topology transitions leave the
source arena, revisions and handles unchanged. Arena Undo/Redo and
`BrepCarrierEditor::from_body` support this shape as well.

## Display representation

`IntersectingBlindBrep::tessellate(segments_per_arc)` constructs a
derived **watertight facet mesh** with analytical arc samples from the
same `Curve3` carriers, one top annulus and a compound nonconvex floor.
The two cylindrical walls reuse **exactly the same mesh ring indices**
at both real cylinder crossing edges.

A complete certificate checks planar support and correct top/floor
normals, inward-facing analytical cylinder facets, two opposing triangle
uses of every undirected mesh edge, projected material area, positive
signed tetrahedral volume and an independently reconstructed polyhedral
area/mass oracle:

`Vfacet = Vbox - d*Apolygon_union`

`Sfacet = Sbox + d*Ppolygon_union`

Display volume is an approximation and must converge toward the
independent analytical volume. It is never written into `BrepModel.mass`.

## CLI

The specialized `geometry.v1` operation:

```json
{
  "operation":"compose_intersecting_blind",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "holes":[
    {"center":[-0.7,0],"radius":1,"depth":2},
    {"center":[0.7,0],"radius":1,"depth":2}
  ],
  "batches":[[
    {"kind":"replace","index":0,"center":[-0.8,0.1],
     "radius":1.1,"depth":2.25},
    {"kind":"replace","index":1,"center":[0.7,-0.1],
     "radius":1,"depth":2.25}
  ]],
  "undo":1,
  "redo":1
}
```

An optional `display_segments` integer (12–512 **per exposed arc**)
requests the certified closed display mesh through the normal `solid`
vertices/triangles/`triangle_face_ids` contract. The
`data.display_quality` certificate distinguishes faceted and exact
analytical mass. When `display_segments` is omitted, no mesh is
computed or returned, and the old stateless response remains a canonical
analytic topology result. Invalid resolutions and failed display
certificates return no partial `data` or `solid` payload.

The result exposes source-validated crossing coordinates, the two arc
domains and source cutter indices, shared vertical crossing edge IDs,
single top and floor face IDs, exact mass/centroid, genus and topology,
one-step edit reports, and Undo/Redo transitions. The old v0.26
`compose_blind` remains the executor for **disjoint** bore regions.

Only this specialized crossing operation authorizes its own canonical
reassembly; `TopologyDeltaPlan::can_apply` remains false globally.
There is no promise of generalized NURBS/plane-cylinder arrangements or
multi-shell Boolean operators yet.

## Tests

First-party `tests/intersecting_blind.rs` pins exact Euler/incidence,
true split `Curve3` arcs, shared cylinder-crossing edge positions,
equal and unequal circle radii, rotated center lines, source permutation,
explicit tangency/containment/differing-depth refusal, lens-area oracle, independent Simpson-integration control for unequal
rotated circle unions, box-relative centroid, three-way void membership,
Arena generational naming/Undo, canonical querying and strong-atomic
edit/Undo.
`tests/intersecting_blind_mesh.rs` separately tests manifold display
closure, nonconvex floor geometry, 12/48/192 sampling convergence,
material triangle lineage, meshing after committed edits and moderate
world translation. Explicit GitHub Actions gates cover these suites
and independent geometry.v1 positive/negative smoke tests.
