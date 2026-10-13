# RustSolid v0.31 — certified one-through/one-blind crossing B-rep

This stage constructs an **actual analytical B-rep**, not merely the read-only
v0.30 footprint plan. It supports exactly one full-height upright cylindrical
through hole and one top-entry finite-depth cylindrical blind hole inside an
independently certified axis-aligned block. Both footprints must intersect
at exactly two well-separated transverse crossing points.

The implementation is first-party mathematical geometry, not a port of
private Siemens algorithms, source code, headers or proprietary data layouts.

## Actual topology

One connected closed orientable shell contains:

| Physical structure | Result |
|---|---:|
| Vertices | 14 |
| Edges | 23 |
| Opposite directed coedges | 46 |
| Faces | 11 |
| Loops | 13 |
| Inner face trim loops | 2 |
| Shells | 1 |
| Euler characteristic (counting face holes) | 0 |
| Genus | 1 |

Face IDs 0 and 1 are the source box bottom and top caps. Both contain one
*actual inner UV loop*. Faces 2–5 are the original box exterior sides.
The five new analytical faces are: two split lower through-cylinder walls
(6–7), one upper through-cylinder wall (8), one exposed upper blind-cylinder
wall (9), and the source-certified blind floor crescent (10).

All new physical circular and vertical edges have exactly two oppositely
oriented coedges, reciprocal next/previous loop links, and analytic 3D
carriers. Each coedge has its own validated UV p-curve on the owning
plane/cylinder support. The floor has reversed support orientation, giving
the correct +Y material-facing normal. The upper aperture is a two-arc
union contour and the lower aperture is the through circle. The blind floor
uses the exterior blind arc plus the complementary through arc.

The original v0.30 area-decomposition planner has seven mathematical
surface-area regions; the canonical B-rep has eleven topological faces
because there are four box exterior faces and two lower through-wall faces.

The connected-shell Euler witness is:

    chi = V - E + F - inner_loops = 14 - 23 + 11 - 2 = 0
    genus = (2 - chi)/2 = 1

The shared canonical BrepModel validator additionally checks
coedge ownership, opposite twins, face-loop cycles, face adjacency and
sampled UV-to-3D carrier agreement. This is not an interval proof for
all floating-point parameters or arbitrary non-manifold B-reps.

## Exact-in-form analytic mass

Let V0 be the block volume, H its height, d the blind depth, At and Ab
the through/blind disk areas, and Ai their source-certified overlap:

    Vremaining = V0 - At*H - (Ab - Ai)*d

The same source geometry supplies first-party analytical surface area and
relative centroid moments. No faceted triangle integration is substituted
for analytic mass.

For the regression fixture (box 6×4×6, radii 0.5 and 0.5, center spacing
0.5 and depth 1.5), Vremaining = 140.14094873829234.

## Strong atomic source edits

MixedCrossingBrep is a dedicated source owner, independent of the
existing MixedThroughBlindBrep multi-cut composer. Its typed transaction
replaces source through/blind parameters, builds a complete isolated
candidate, validates it and publishes it atomically in one revision.
Undo and Redo reassemble independently validated source states. Failed
or poisoned transactions do not change the current mass or topology.

IDs are local to the current B-rep reconstruction; durable cross-regeneration
topological naming is NOT claimed.

## Derived watertight display

MixedCrossingMesher triangulates the exact source-bound circle arcs,
the original box sides, both cap openings, and the crescent floor,
sharing the physical boundary vertex IDs across neighboring faces.
The display is rejected unless:

- Every triangle edge has two oppositely directed usages.
- All triangles have owning analytic face IDs.
- Signed mesh volume is positive and agrees with independently computed
  polygonized cap-section volume.
- Faceted area agrees with polygonized cap/side/cylinder area.
- Faceted volume cannot undershoot the analytic volume beyond a bounded
  numeric slack.
- Display Euler characteristic is zero (one through tunnel).

The segment budget is 12..=512 per arc. Increasing display resolution
reduces facet error without modifying the analytic B-rep or its mass.

## Public APIs

Rust library:

~~~rust
use rustsolid::{Solid,Point3,HoleParameters,BlindHoleParameters,MixedCrossingBrep};
let block=Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let through=HoleParameters{center_x:0.,center_z:0.,radius:0.5};
let blind=BlindHoleParameters{center_x:0.5,center_z:0.,radius:0.5,depth:1.5};
let body=MixedCrossingBrep::from_box(&block,through,blind)?;
body.validate()?;
assert_eq!(body.model().summary().genus,1);
let display=body.tessellate(32)?;
assert_eq!(display.genus,1);
# Ok::<(),rustsolid::GeometryError>(())
~~~

geometry.v1 JSON:

~~~json
{
  "operation":"compose_mixed_crossing",
  "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
  "through":{"center":[0,0],"radius":0.5},
  "blind":{"center":[0.5,0],"radius":0.5,"depth":1.5},
  "display_segments":32
}
~~~

The endpoint returns the committed *dedicated* analytic B-rep status,
topology/genus, analytic volume and material orientation roles. A derived
triangle mesh and faceting-quality metadata appear only when requested.
It supports source-parameter batches and Undo/Redo, with invalid inputs
returning neither partial solid nor partial data. The generic Boolean and
TopologyDeltaPlan paths stay unchanged.

## Verify

~~~sh
cargo test --test mixed_footprint
cargo test --test mixed_crossing_surfaces
cargo test --test mixed_crossing_brep
cargo test --test mixed_crossing_mesh
cargo test --bin gefest-geometry mixed_crossing_cli
cargo test --all-targets
cargo check --release --all-targets
~~~

## Explicit unsupported operations

No intersecting third cylinder, multiple mutually interacting holes,
oblique cutter, open sheet, arbitrary curved Boolean, full tolerant CAD
healing, NURBS, X_T/STEP exchange, persistent names across connectivity
changes, or production-certified extreme-precision operation is implied.

The old multi-hole MixedThroughBlindBrep builder still rejects
through/blind footprint contacts. General TopologyDeltaPlan remains
preview-only and cannot be committed.
