# RustSolid v0.35 — a canonical analytical B-rep for a connected three-cylinder chain

This is a **real executable B-rep**, not a read-only region plan.
The initial scope is intentionally narrow: one full-height vertical through
cylinder intersects a deeper blind cylinder; that blind cylinder intersects
a shallower blind cylinder. The through cutter and the shallow cutter
must stay strictly separated in XZ. Both finite blind depths must be
strictly distinct and the disks entirely inside a certified axis-aligned box.

The connected three-cylinder chain is NOT a general simultaneous triple
junction. Near tangency, internal/external contact, union holes, invalid
depth order and unrelated source solids are rejected before publication.
The generic curved Boolean planner remains preview-only.

## What is physically constructed?

Each source cylinder's angular chart is divided at the union of **all**
intersection event angles from each horizontal depth slab. This means
the same physical junction has one canonical vertex ID on each plane,
even when a cylinder becomes exposed or hidden at the next blind floor.

Each exposed circular sector, over its physical height interval, has:
- its own analytical cylindrical face with a true CylinderSurface carrier;
- two shared circular edges (one on each depth level);
- two shared vertical line edges;
- one closed four-coedge loop with source-matched 2D UV trims.

Both original source box caps own a true inner circular/compound trim
loop. The two different blind-floor heights own two separate planar floor
faces. Their physical loop boundaries use the exact horizontal edges
that terminate one slab or emerge in the next. Floor faces have upward
material-facing normals; subtractive cylinder surfaces have inward
material-facing normals.

The resulting B-rep contains one orientable connected closed shell.
Its genus is one and its Euler characteristic is zero. These values,
and the matching coedge incidence, loop cycles and UV/3D trim agreement,
are checked by the canonical BrepModel validator and source replay.
The surface patches may contain genuine source-aligned subdivision seams;
there are no faceted CAD faces or invented intersection edges.

## Exact-in-form source volume, area, centroid

The independent v0.32/v0.34 source planners report three constant-active
height slabs. Source material volume equals the original block volume minus
the sum of union-disk area multiplied by each slab's height.

Analytical centroid moments are computed from first-party Green integrals
of the actual exposed circles, not derived from mesh triangles.
For each arc, the planar first moments follow:

    Mx = (1/2) integral(x^2 dz)
    Mz = -(1/2) integral(z^2 dx)

These are integrated over the physical height intervals. An independent
test oracle uses cosine-smoothed Simpson integration over the union of
2D vertical disk intervals, including the first moments of that union.
Tests check source-analytic volume and all centroid coordinates against
that different algorithm.

Surface area sums the actual exposed cylinder sectors, original side
walls, both pierced caps and two material-facing floors.

## Derived watertight mesh

ConnectedTripleMesher samples **each canonical B-rep edge only once**.
All adjacent cylinder, cap and floor faces reuse the same physical
vertex indices. Each planar face is triangulated with its own real
outer and inner loops; cylinder patches are sampled between their paired
upper/lower circular edges.

A candidate is rejected unless every undirected triangle edge has
exactly two opposite directed uses, all face owners exist, signed volume
is positive, mesh Euler characteristic is zero and facet volume agrees
with an independent polygonized horizontal-section volume calculation.

The analytic B-rep and mass never change when display resolution changes.
Supported sampling is 12 to 128 segments per source arc. Floating-point
budgets are explicit but are not a mathematically formal interval proof
for all real coordinates.

## API and example

Rust source: ConnectedTripleBrep::from_box(source, TripleCutParameters).
The builder supports source-bound replace_through/replace_blind
transactions, atomic commit, Undo and Redo. Invalid candidates publish
no partial topology.

JSON geometry.v1 operation: compose_connected_triple

~~~json
{
  "operation":"compose_connected_triple",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "through":{"center":[0,0],"radius":0.75},
  "blind":[
    {"center":[0.7,0],"radius":0.75,"depth":2.0},
    {"center":[1.4,0.25],"radius":0.6,"depth":1.0}
  ],
  "display_segments":24
}
~~~

The blind cutters must be ordered deep first, shallow second.
A successful result returns committed analytic topology, two floor
faces, exact-in-form mass, three verified depth slabs and optional
watertight faceted display. Failed operations return no partial
solid/data; general TopologyDeltaPlan remains non-executable.

## Acceptance commands

~~~sh
cargo test --test connected_triple_brep
cargo test --bin gefest-geometry connected_triple_cli
cargo test --all-targets
cargo check --release --all-targets
~~~

The implementation uses independently derived computational geometry;
it does not embed Siemens/Parasolid proprietary code, decompilation,
private headers or undisclosed SDK layouts.
