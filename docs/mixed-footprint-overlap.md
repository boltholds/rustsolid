# RustSolid v0.30 — source-bound through/blind footprint overlap plan

The first v0.30 capability is deliberately **read-only mathematical planning**:
one full-height upright cylindrical through cut and one top-entry finite-depth
cylindrical blind cut, both strictly contained by a certified axis-aligned box.
It does not modify, approximate, or construct a curved B-rep. The existing
MixedThroughBlindBrep builder still refuses **all** overlapping through/blind
footprints. General TopologyDeltaPlan continues to report can_apply=false.

This is independently developed geometry; no proprietary kernel code,
decompiler output, private ABI layouts or third-party CAD oracle are used.

## Exact source and classifications

MixedFootprintOverlapPlanner::plan accepts a real validated Solid (certified
as a six-face axis-aligned block), a HoleParameters through cutter and a
BlindHoleParameters finite-depth cutter. Both horizontal circular footprints
must have positive radius and strict side-wall clearance. The blind depth
must lie strictly between zero and full box height. The circle classification
is one of:

- Disjoint
- TwoTransverseCrossings
- BlindInsideThrough
- ThroughInsideBlind
- Coincident (exactly identical circle center and radius)
- HeldNearTangency (near external/internal tangency, near coincidence or
  insufficient geometric representability)

Held has no overlap-area or volume prediction. It never invents cut edges.

For two transverse circle crossings the plan contains two points on both
actual source circles at the box top plane, and one *counterclockwise exposed
outer arc* per source circle. Both arcs connect at the crossing witnesses and
their midpoints lie outside the other circle. These are mathematical witness
arcs only, not B-rep edges or topology deltas.

## Independent mass law

Let V0 be box volume, H the box height, d the blind depth, At and Ab the
through/blind footprint areas and Ai their intersection area. Then

    V_remaining = V0 - At * H - (Ab - Ai) * d

The through cutter has already removed all height H of its footprint.
The blind cutter additionally removes material outside the through footprint,
but only over depth d. This gives the disjoint, strictly crossing, contained
and coincident cases without treating generated triangles as geometry.

For a pair of strictly crossing disks with center distance D and radii rt/rb:

    x1 = (D² + rt² - rb²) / (2D)
    x2 = D - x1
    h  = sqrt(rt² - x1²)
    Ai = rt² * atan2(h, x1) + rb² * atan2(h, x2) - D*h

The implementation requires positive finite dimensions and numerical
certificates on both crossing points; near tangent contacts are held before
computing a suspect lens or topology. A model-space length budget combines
source tolerance and a floating-point representability allowance for
large translated world coordinates.

With a box of dimensions 6×4×6, rt=rb=0.5, D=0.5, d=1.5:

    Ai          = 0.3070924246521893
    V_remaining = 140.14094873829234

The tests independently integrate the **minimum disk strip width**, using a
cosine change of variable and composite Simpson integration. This verification
does not reuse the closed-form lens formula.

## Rust library

Example (read-only, no topology edits):

~~~rust
use rustsolid::{
    Solid, Point3, HoleParameters, BlindHoleParameters,
    MixedFootprintOverlapPlanner, MixedFootprintContact,
};

let block = Solid::block(Point3{x:-3.,y:-1.,z:-3.}, 6., 4., 6.)?;
let through = HoleParameters{center_x:0.,center_z:0.,radius:0.5};
let blind = BlindHoleParameters{
    center_x:0.5,center_z:0.,radius:0.5,depth:1.5,
};
let plan = MixedFootprintOverlapPlanner::plan(&block,through,blind)?;
assert_eq!(plan.contact,MixedFootprintContact::TwoTransverseCrossings);
assert!(!plan.can_apply());
plan.validate_against(&block,through,blind)?;
assert!((plan.predicted_remaining_volume.unwrap()-140.14094873829234).abs()<1e-9);
# Ok::<(),rustsolid::GeometryError>(())
~~~

The plan carries the source bounds, analytic source volume, source revision,
tolerance and hole parameters, and validation recomputes every witness. It
does not create any new face, loop, coedge, handle or Undo/Redo transaction.

## JSON geometry.v1

~~~json
{
  "operation": "mixed_footprint_overlap_plan",
  "source": {"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
  "through": {"center":[0,0],"radius":0.5},
  "blind": {"center":[0.5,0],"radius":0.5,"depth":1.5}
}
~~~

For strict crossing, the response has ok=true, solid omitted, data.kind set
to mixed_footprint_overlap_plan, state=two_transverse_crossings, exact-in-form
overlap_area and predicted_remaining_volume, two crossing points and two
exposed arcs. For a held near-tangent, both areas and witnesses are null,
can_apply=false, topology_committed=false. Invalid box or dimensions fail
with ok=false and no data/solid.

## Acceptance tests and boundaries

Run:

~~~sh
cargo test --test mixed_footprint
cargo test --all-targets
~~~

Tests check the independent strip integral, unequal radii and permutations,
all contact classes, source-circle residuals, exposed-arc closure, mass
inclusion/exclusion, depth monotonicity, coordinate translation, refusal at
tangency, source-plan tampering, and legacy Boolean refusal.

**Not implemented:** joined top opening, split-height through walls, blind
floor trimmed by an existing through void, cylinder face intersection edges,
global curved Boolean, manifold/genus proof or triangulation of the new
mixed intersection. These must appear in a separate geometry/topology PR with
independent mass, Euler, UV/3D trim, material-sense and watertight tests.

The preexisting docs/parasolid-api-crosswalk.md reports historical v0.14
scope and is not a current v0.30 feature matrix.
