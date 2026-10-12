# RustSolid v0.29 — one exact analytical through/blind B-rep

RustSolid v0.29 composes **full-height through-cylinder cuts and
independent finite-depth blind regions in the same connected analytical
B-rep**. The canonical solid preserves real 3D circular edges, cylindrical
surfaces, 2D UV trims, coedge incidence and material-facing orientation.

This milestone closes a practical feature gap in the v0.23 through-hole
and v0.28 compound-blind tools. It does **not** execute arbitrary curved
Booleans, intersect a through-bore with another bore, or apply the
preview-only general `TopologyDeltaPlan`.

## Geometry admitted by the executor

`MixedThroughBlindBrep::from_box(&Solid, MixedThroughBlindParameters)`
starts from one source-certified six-face axis-aligned block. The
`through` list has at least one upright, full-height cylindrical
cut. The `blind` list also has at least one top-entry, finite-depth
circular cutter. Together they contain at most **16 source cylinders**.

The source-proven through-hole builder checks every cutter's radius,
clearance from four side walls and pairwise through-hole separation.
The v0.28 blind-region planner verifies every blind cylinder using its
source-bound top-opening `IntersectionGraph`; it groups blind holes
into connected regions of either one circular cutter or precisely two
transversely crossing equal-depth cutters. Distinct blind regions must
be strictly separated. The existing rejection rules for chained/
three-way crossings, containment, near tangency and different depths
inside one crossing component remain in force.

**Every through footprint must be strictly separated from every blind
footprint.** Touching, overlapping or nearly contacting
through/blind boundaries are rejected regardless of whether the blind
bore ends above the lower box plane. Those heterogeneous intersections
would require new side-wall trim and shared-edge construction.

There is no silent projection, faceted substitution or approximation of
unsupported shape combinations.

## Topology and the correct genus

The B-rep is assembled from complete independent analytical donors:

1. `multi_through_hole::build_multiple` constructs a genuinely genus-T
   shell with T full-height cylindrical walls. Both top and bottom box
   faces have one inner loop per through-bore.
2. The v0.28 `plan_groups` and `graft_region` insert complete blind
   source regions into that shell, remapping vertices, 3D curves,
   edges, coedges, UV pcurves, loops, surfaces and faces.
3. Every blind component contributes an **additional** inner loop
   on the shared top cap, one true upward-facing planar floor, and
   one or two inward-facing analytical cylinder patches. Blind floors
   do **not** create inner loops on the box bottom.
4. All added cylinder and floor faces use `FaceSense::Reversed` relative
   to their mathematical parameterization. Original box faces remain
   forward. The entire candidate passes `BrepModel::validate()`.

For T through holes, S independent singleton blind regions and P
two-circle crossing blind regions:

| Quantity | Verified formula |
| --- | ---: |
| Vertices | 8 + 2T + 2S + 4P |
| Edges | 12 + 3T + 3S + 6P |
| Coedges | twice the edge count |
| Faces | 6 + T + 2S + 3P |
| Loops | 6 + 3T + 3S + 4P |
| Top cap loops | 1 + T + S + P |
| Bottom cap loops | 1 + T |
| Inner UV loops across all faces | 2T + S + P |
| Shells | 1 |
| Euler characteristic | 2 − 2T |
| **Genus** | **T** |

The crucial test is that adding a blind floor changes neither the
number of through tunnels nor the genus of the closed shell. A
crossing blind pair has one compound aperture and a single united
floor, not two phantom overlapping planar floors.

Dense edge/face/coedge IDs of freshly assembled components are
**not durable feature IDs**. The mixed body can be exported as a
read-only `BrepArena` snapshot; source-specific parametric edits
remain inside the `MixedThroughBlindBrep` owner rather than
falsely preserving names across topology repartition.

## Exact mass, independent from tessellation

Let the uncut box have volume `V0`, area `A0`, centroid `C0`
and height `H`.

For each through hole `t`, its removed area is `At = πrt²`,
removed volume `Vt = At H`, side area `2πrt H`. Because it
passes through the box the upper and lower circular cap areas are
**both removed**; its net area change is
`2πrt H − 2 At`.

For each independent blind region `b` (singleton or crossing pair),
take its exact planar area `Ab`, exposed perimeter `Pb` and
depth `db`. Its net area change is `Pb db`; the top circular
opening is replaced by its independent floor.

Consequently:

```
V = V0 − Σ(At * H) − Σ(Ab * db)
A = A0 + Σ(2πrt*H − 2At) + Σ(Pb * db)
```

All XZ mass moments use the original box centroid as origin. Through
cylinder removal has unchanged axial Y centroid, whereas a blind
region's centroid lies at `topY − db/2`. For crossing blind pairs,
the source-proven exact union-of-two-circles moments from v0.27
are reused. The final centroid is the original centroid minus all
relative removed-volume moments divided by the remaining volume.

No triangle count enters any of these formulas.

## Atomic parameter edits

`MixedThroughBlindBrep` stores the verified original box, source
parameter vectors, actual connected blind-component descriptors,
a complete analytical model, revision and bounded Undo/Redo history.

One `MixedThroughBlindTransaction` may stage
`replace_through(index, HoleParameters)` and
`replace_blind(index, BlindHoleParameters)` together. A single
commit must revalidate *every* source hole, the cross-family spacing,
the blind crossing grouping, all analytical graph links and exact mass
before any published state changes. A malformed staging command
poisons the transaction. Late failure due to invalid overlap,
tangency or unsatisfied depth constraints leaves body, history and
revision unchanged. Undo/Redo use independent validated rebuilds
with monotonically increasing runtime revision.

The lists and group topology may change when a blind pair joins or
splits, but this milestone changes **parameters of existing source
cylinders**, not their source-list length.

## Derived closed display mesh

`MixedThroughBlindBrep::tessellate(segments)` uses the analytical
circles/arc segments already bound to the B-rep and derives a
closed, oriented triangle representation.

The top cap is one plane with **all** through and blind apertures.
The bottom cap is a different plane with **only through apertures**.
Both triangulations come from their own polygon-with-holes charts
through `earcutr`, with explicit area and exclusion checks.

Every blind region separately triangulates its genuine floor,
including the nonconvex two-arc crossing floor. Every through or
blind cylinder side shares the *same* 3D ring vertex IDs with its
corresponding top/bottom cap or blind floor, including the two
crossing-edge endpoints of paired blind cylinders.

The mesher checks:
- all mesh edges have exactly two oppositely directed triangle uses;
- all planar facets lie on their actual analytical planes;
- cylindrical facets lie on their mathematical radius and face the
  empty bore according to `FaceSense::Reversed`;
- upper cap and lower cap projected material areas match their
  **different** aperture sets;
- signed tetrahedral volume agrees with an independent polygonal
  swept-region volume oracle;
- facet surface area agrees with a polygonal perimeter/area oracle;
- facet volume does not materially lie below the exact
  analytical volume, and converges toward it as segmentation rises.

`MixedHoleDisplay` reports both cap triangle/loop counts, per-blind-floor
certificates, faceted volume, surface area, exact analytical volume and
faceting volume error. It never overwrites `BrepModel.mass`.

## Rust API

```rust
use rustsolid::{
    MixedThroughBlindBrep,MixedThroughBlindParameters,
    HoleParameters,BlindHoleParameters,Point3,Solid,
};
let source=Solid::block(
    Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let mut body=MixedThroughBlindBrep::from_box(
    &source,
    MixedThroughBlindParameters {
        through:vec![HoleParameters {
            center_x:-1.7,center_z:0.0,radius:0.5
        }],
        blind:vec![BlindHoleParameters {
            center_x:1.5,center_z:0.0,radius:0.75,depth:2.0
        }],
    },
)?;
assert_eq!(body.model().summary().genus,1);
body.validate()?;
let display=body.tessellate(48)?;
assert_eq!(display.bottom_face_loop_count,2);
assert_eq!(display.top_face_loop_count,3);
```

## geometry.v1

The `compose_mixed_regions` request is independent of the old
`compose_blind_regions` and `structural_transaction` APIs:

```json
{
  "operation":"compose_mixed_regions",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "through":[{"center":[-1.7,0],"radius":0.5}],
  "blind":[{"center":[1.5,0],"radius":0.75,"depth":2}],
  "display_segments":48
}
```

The result includes `data.kind="mixed_through_blind_regions"`,
actual topology/genus, blind-region source lineage, independent floor
face IDs, exact mathematical mass, optional `solid` for display and
`data.display_quality` with independent closure and volume witnesses.

If `display_segments` is omitted the response contains analytical
topology and mass **without** a rendered `solid` payload.
Unsupported cutters or invalid display requests return
`ok:false` with neither data nor partial mesh.

The same analytical body is inspectable with
`inspect_brep`, `geometry_query` and the shape-preserving
`carrier_transaction` using B-rep source kind
`analytic_mixed_regions` and fields
`origin,width,height,depth,through,blind`.

## Known limits and future gate

This is a **bounded heterogenous feature family**, not a universal
analytical Boolean kernel. The following remain unsupported:

- a through-cylinder crossing or touching any blind-cylinder or its
  compound trim (even if their Y-depth ranges could avoid each other);
- depth changes across two **overlapping** blind circles, three-way or
  chained blind crossings and tangential contacts;
- oblique/side-entry bores, curved support surfaces beyond supported
  cylinders, NURBS, non-box source bodies and multi-shell cavities;
- generic intersection-driven splitting and rebuilding of arbitrary
  analytical faces through `TopologyDeltaPlan`;
- permanent CAD feature naming across unrelated regeneration.

The kernel continues to refuse those cases rather than fabricating
faceted Boolean results. The next research/implementation gate should
be analytical **stepped blind-floor intersections at differing
depths**, including vertical patch split curves and planar shelf
regions, with genuinely independent mass verification.
