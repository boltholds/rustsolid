# RustSolid v0.32 — multi-cylinder planar region graph and depth-slab volume planner

RustSolid v0.32 adds a **source-bound, read-only analytical planner** for
several interacting upright cylindrical cuts in one certified axis-aligned
block. It handles up to 12 source cutters, with at least one full-height
through cylinder and one top-entry blind cylinder. Blind cylinders can have
different finite depths and their circular footprints may overlap with one
another and with one or more through cylinders.

This is NOT yet a new executable multi-cylinder Boolean or an editable
B-rep. The already certified v0.31 one-through/one-blind producer stays
available separately. General curved B-rep topology commits remain refused.

## Mathematical construction

Let y0 and y1 be consecutive levels of a sorted event sequence consisting
of the block bottom, each independent blind floor, and the block top.

Within each open slab (y0, y1), each through-cylinder footprint is always
active. A blind cylinder of depth d is active only above (topY-d).
Therefore the active cutter set is constant in that slab.

For each active 2D circle, all strict pair intersections are calculated in
source coordinates. Their angular parameters split the circular boundary
into candidate intervals. Intervals fully covered by another source disk
are discarded; each remaining arc has one exact-in-form circular carrier
and a midpoint that lies outside every other active disk. Arcs reference
the same source-bound pair-intersection vertex IDs on both incident circles.

The directed arc graph must satisfy precisely one outgoing and one incoming
arc at each exposed boundary vertex. Its arcs are assembled into complete
oriented boundary cycles. Positive signed cycles are exterior boundaries of
connected disk-union components; negative cycles delimit genuine holes in
the 2D disk union. These planar hole counts are NOT interpreted as a 3D
solid genus before constructing and validating actual 3D B-rep incidence.

Each directed arc contributes its exact-in-form Green line integral:

    area_arc = 1/2 * (r² * (theta_end-theta_start)
                   + r*(cx-xref)*(sin(theta_end)-sin(theta_start))
                   + r*(cz-zref)*(cos(theta_start)-cos(theta_end)))

The reference origin is taken from a local cutter, avoiding unnecessary
world-coordinate cancellation. Outer and hole cycle integrals are summed
with a compensated accumulator. No triangulation or rasterization is used
to construct the analytical source area.

The removed material volume is obtained by summing the active planar union
area in each constant-active interval:

    V_removed  = Σ [A_union(active cutters at slab k) * (y_k+1 - y_k)]
    V_remaining = V_source - V_removed

This also handles connected three-circle chains, three pairwise crossings,
isolated components, nested hidden disks, different blind depths and
multiple closed planar boundary loops, without assuming pairwise disk
overlap areas can simply be added.

## Fail-closed geometric events

The planner produces MultiCutPlanState::Held(reason) with NO volume and
NO partial slab collection for unresolved events such as:

- source circles touching externally/internally, near coincidence;
- circle pair crossings lying on or too close to a third circle support;
- sub-resolution angular intervals or incomplete directed boundary cycles;
- distinct blind floor levels too close to separate reliably;
- uncertain signed area or physically impossible area bounds.

Invalid inputs (non-box source, nonfinite cutters, missing source family,
out-of-range blind depth/radius, lateral wall collisions or more than
12 cutters) return an explicit error. Numerical resolution slack is
reported but is NOT a rigorous IEEE-754 interval error certificate.

When Ready, the plan exposes source-indexed arcs, physical intersection
points, oriented planar loops, connected disk-union component counts, hole
counts, each slab's active cutter list and predicted volume. A stored
plan is validated by independent source replay of all witnesses.
The original source B-rep, its revisions, editor handles and history
are never changed. can_apply() is ALWAYS false.

## Independent numerical cross-check

Tests do NOT use the Green contour formula as their truth oracle. A
different implementation integrates the union of 2D vertical disk
intervals along x with composite Simpson quadrature and cosine endpoint
smoothing. The tests compare each slab union area and the final swept
volume for three interacting cylinders at two different blind depths,
an all-pairwise-crossing triple and a three-disk annular arrangement with
one planar hole. They also test area monotonicity, translation,
permutation, nested/disjoint source hiding, near triple junction holds,
malformed input refusal, and tampered-plan rejection.

Because independent numerical quadrature has discretization error,
these tests use explicit finite tolerances. CI success does not imply
formal verification or production-level robustness for all CAD scales.

## Rust interface

~~~rust
use rustsolid::{
    Solid,Point3,HoleParameters,BlindHoleParameters,
    MultiCutRegionPlanner,MultiCutPlanState
};
let block=Solid::block(
    Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.
)?;
let through=[HoleParameters{center_x:0.,center_z:0.,radius:0.75}];
let blind=[
    BlindHoleParameters{center_x:0.7,center_z:0.,radius:0.75,depth:2.0},
    BlindHoleParameters{center_x:1.4,center_z:0.25,radius:0.6,depth:1.0},
];
let plan=MultiCutRegionPlanner::plan(&block,&through,&blind)?;
assert_eq!(plan.state,MultiCutPlanState::Ready);
assert_eq!(plan.slabs.len(),3);
assert!(!plan.can_apply());
plan.validate_against(&block,&through,&blind)?;
# Ok::<(),rustsolid::GeometryError>(())
~~~

## JSON geometry.v1

~~~json
{
  "operation": "multi_cut_region_plan",
  "source": {"kind":"block","origin":[-3,-1,-3],
             "width":6,"height":4,"depth":6},
  "through": [{"center":[0,0],"radius":0.75}],
  "blind": [
    {"center":[0.7,0],"radius":0.75,"depth":2},
    {"center":[1.4,0.25],"radius":0.6,"depth":1}
  ]
}
~~~

Successful responses have data.kind="multi_cut_region_plan",
data.state="ready", source indices, per-slab boundary loop/arc witnesses,
union areas and predicted_remaining_volume. When held,
data.state="held" and hold_reason explain the refusal, with no area/
volume claim. All responses have can_apply=false, topology_committed=false
and omit the solid result.

The existing "compose_mixed_crossing" operation remains the only
dedicated through/blind overlap B-rep executor in this release.

## What would unlock an executable generalization?

Before any multi-cylinder topology commit, construct globally consistent
3D edge/face carrier ownership across every blind floor event, splice
layer-changing topological contours, establish mutually opposite coedges,
physical shell adjacency, UV/3D trim equality, material orientation,
multi-shell/void component semantics, Euler/genus proof, stable transaction
lineage and an independently certified watertight display. Exact-in-form
volume and closed planar arcs ALONE are insufficient for that claim.

RustSolid remains a clean-room first-party kernel: proprietary vendor
binary implementation details, disassembly, SDK internals and source
fragments are not imported into this project.
