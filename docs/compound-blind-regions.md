# RustSolid v0.28 — composed analytical blind-region components

RustSolid can now reconstruct and edit a **single connected canonical B-rep**
containing several independent top-entry blind-cut regions. Each independent
region is either (a) one circular bore, or (b) the union of **exactly two
partially overlapping, nonnested, equal-depth circular bores**. Both kinds
may coexist on the same box's top face, provided different regions are
strictly separated.

This is a useful executable **component arrangement** built from the
original v0.25/v0.27 analytical cut witnesses, not an implementation of
universal arbitrary-face CSG or of `TopologyDeltaPlan::can_apply`.

## Bounded source and arrangement contract

`BlindRegionComposer::build(&Solid, &[BlindHoleParameters])` requires a
certified six-face axis-aligned box and **1–16 top-entry cylindrical
holes**, each independently passing v0.25 source-bound closed opening
intersection certification. The input parameters give center X/Z,
radius, and finite positive depth from the box upper (+Y) plane.

The planner examines every pair of circle footprints in model-local XZ:

- Distinct, separated circles form different region components.
- Two partially overlapping, transversely crossing circles form **one**
  composite component only when their lower floor planes coincide
  within floating-point representability.
- Exact/near tangency, containment, coincidence and unresolved near
  contacts are rejected.
- Any circle intersecting **two other circles** is rejected, including
  chained or triple crossings. More complex planar arrangements require
  a general arc-intersection scheduler not claimed by this release.

An independent single circle may have a different depth from a crossing
pair or another single circle: their footprints do not interact.
No oblique, side-entry, intersecting-through, multi-shell or cross-depth
overlapping cuts are silently approximated.

The connected components are ordered deterministically by their
smallest source-hole input index. `BlindRegionComponent` retains
`source_indices`, `kind`, true circular/arc edge ownership, top loop
ID, cylindrical side face IDs, floor face ID, analytical region area,
exposed perimeter and depth.

## True topology construction

Each component is independently constructed by the corresponding
validated v0.25 `BlindHoleBrep` or v0.27 `IntersectingBlindBrep`.
The compositor extracts **only** that donor's non-box vertices, edges,
coedges, loops, supporting surfaces and 2D/3D analytic carriers. Every
cross-reference is remapped to the final dense destination catalog,
while the six original box faces and their carrier IDs stay untouched.

The source box's **single top planar face** gains exactly one distinct
inner UV loop for each independent aperture region:

- For an isolated bore, that loop is a real analytic closed circle.
- For a two-circle crossing region, that loop is a composed two-arc
  analytical boundary. Both cylinder patches share the true two
  vertical crossing edges and their opposing coedges.

Each independent region also gets its own **upward-facing analytical
planar floor**, sharing lower analytic circle or crossing-arc boundaries
with its cylinder walls. Material face normals use the typed
`FaceSense::Reversed` on every interior cylinder and floor face,
distinct from the UV support's parameterized orientation.

Let `S` be the number of singleton regions and `P` the number
of two-circle crossing regions. The canonical B-rep has:

| Property | Count |
| --- | ---: |
| Vertices | 8 + 2S + 4P |
| Edges | 12 + 3S + 6P |
| Coedges | 24 + 6S + 12P |
| Faces | 6 + 2S + 3P |
| Loops | 6 + 3S + 4P |
| Top inner loops | S + P |
| Independent planar floors | S + P |
| Analytic cylindrical faces | S + 2P |
| Shells | 1 |
| Euler characteristic | 2 |
| Genus | 0 |

Every completed body is validated against the common canonical
`BrepModel` manifold, twin/next/prev, UV-to-3D, material sense,
closed-shell and Euler/genus checks. A failed attachment or regional
carrier mismatch aborts the **entire candidate**.

The operation is allowed to renumber component-specific dense IDs when
sources are reclassified or reordered; **those numbers are not durable
CAD feature IDs**.

## Exact mass and centroid

The original box has volume `V0`, area `A0`, centroid `C0`.
For each aperture component `i`, take its **analytic** XZ union area
`Ai`, exposed analytical circle-arc perimeter `Pi` and depth `di`.

```
V = V0 − Σ(Ai*di)
A = A0 + Σ(Pi*di)
```

For singleton regions `Ai=πri²` and `Pi=2πri`. For crossed pairs,
the v0.27 source-verified `CrossedDiskRegion` computes the exact
two-disk union area, exposed perimeter and area first moments using
circular segment integrals.

The actual removed volume centroid Y-coordinate is
`topY−di/2`, while X/Z moments are evaluated relative to the original
box centroid. These independent moments combine into the final exact
analytic `MassProperties`. Faceting resolution cannot affect the
committed volume, surface area or centroid.

## Strong-atomic source parameter updates

`CompoundBlindBrep` owns the original certified box, an ordered
vector of source `BlindHoleParameters`, the currently connected
analytical components, the final canonical `BrepModel`, and bounded
Undo/Redo parameter snapshots.

`CompoundBlindTransaction::replace(source_index, new_parameters)`
stages one or more source modifications. On commit the full planner,
component donors, analytical topology graft and mass calculator run
on an **isolated candidate**. Only a completely validated result can
replace the current model in one revision/Undo step.

A group may **split** into two independent holes, or two separate holes
may **join** as a crossing pair, provided every new configuration meets
the same strict shape and collision contract. Source indices remain
stable within the fixed-length parameter vector, but dense topological
names are regenerated. The API does not pretend to preserve arbitrary
cross-regeneration face handles.

An invalid staging operation poisons the transaction, and a valid
sequence followed by an invalid late geometric result leaves the
published geometry, revision, mass, Undo/Redo branches and original
source unchanged. Undo/Redo revalidate entire candidates and increase
revision monotonically.

`BrepBody::CompoundBlind` supports the common analytical B-rep view,
read-only `BrepArena::from_compound_blind`, detached shape-preserving
carrier editing, `geometry_query`, and correct material/boundary
classification. Mutation of a compound body **through the arena's
dense-ID semantic handle registry is not offered** until an explicit
stable component-to-entity lineage journal exists.

## Derived mesh and independent certificate

`CompoundBlindBrep::tessellate(segments_per_circle)` samples the
current true 3D analytic circular carriers, projects their complete
boundary rings into a stable box-local XZ chart and reconstructs:

- One *shared top planar face* triangulated with as many distinct inner
  aperture loops as there are region components.
- One analytically owned floor per component. Nonconvex two-arc
  crossing floors are sampled separately, not approximated as two
  complete overlapping disks.
- Inward-facing cylinder side facets. Both crossing cylinders share
  exactly the same ring indices at each common vertical edge.
- The original lower cap and four lateral box face triangles, with
  their original B-rep source face IDs.

Earcut's polygon triangulation is a **proposal** for viewport display
only, independently rejected if face area, material containment, per
triangle source/normal/carrier distance, edge-pair closure, positive
signed volume or faceted area disagree with analytical source geometry.

The polygonized component aperture area `Apoly_i` and perimeter
`Ppoly_i` supply independent display oracles:

```
Vfacet = Vbox − Σ(Apoly_i * depth_i)
Sfacet = Sbox + Σ(Ppoly_i * depth_i)
```

The display result `CompoundBlindDisplay` reports the mesh, separate
polygonized floor region certificates, top face loop count, derived
volume/area, true analytic volume, and faceting error. It never modifies
canonical `BrepModel.mass`, shape topology or revision.

Resolution is constrained to 12–512 segments per source circle, and
8,192 sampled source-circle points in one display request.

## Rust

```rust
use rustsolid::{
    BlindHoleParameters,CompoundBlindBrep,Point3,Solid
};
let block=Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let holes=vec![
    BlindHoleParameters{center_x:-1.5,center_z:-0.5,radius:0.75,depth:1.8},
    BlindHoleParameters{center_x:-0.4,center_z:-0.5,radius:0.70,depth:1.8},
    BlindHoleParameters{center_x: 1.7,center_z: 1.2,radius:0.55,depth:2.6},
];
let mut body=CompoundBlindBrep::from_box(&block,holes)?;
assert_eq!(body.components().len(),2);
assert_eq!(body.model().summary().genus,0);
let display=body.tessellate(48)?;
assert_eq!(display.floor_regions.len(),2);
```

## geometry.v1

The specialized executable operation is `compose_blind_regions`:

```json
{
  "operation":"compose_blind_regions",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "holes":[
    {"center":[-1.5,-0.5],"radius":0.75,"depth":1.8},
    {"center":[-0.4,-0.5],"radius":0.7,"depth":1.8},
    {"center":[1.7,1.2],"radius":0.55,"depth":2.6}
  ],
  "display_segments":48
}
```

The `data` result includes topology, exact mass, each component's
source indices, area, perimeter, loop and floor/side face IDs. If
`display_segments` is specified, a separate `solid` viewport mesh
includes real triangle-to-B-rep-face ownership and a
`data.display_quality` certificate.

Optional `batches:[[{kind:"replace",index,center,radius,depth},...]]`
and `undo` / `redo` form one source-verified structural
repartitioning transaction per batch. Missing display parameters
continue to produce analytical results without unnecessary faceting.

A normal `inspect_brep` or `geometry_query` source may use
`kind:"analytic_compound_blind"` with the same block, `holes`, and
an optional `segments`.

## Deliberate limits and next gates

This release **does not** create arbitrary planar arrangements when
three or more source circles overlap, handle depth-changing stepped
floors, build multiple material shells, implement general NURBS,
or turn global `TopologyDeltaPlan` into an executable Boolean journal.
General topological name persistence across changing region connectivity
is likewise not claimed.

The next rigorous extension is a proper planar **circle-arc arrangement
graph** that supports 3+ mutually intersecting contours (including
junction valence and robust same-face trim splits), then a heterogeneous
step-floor/oblique-surface region scheduler with independently proven
mass and manifold continuity.
