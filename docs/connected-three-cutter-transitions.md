# RustSolid v0.34 — connected three-cylinder cross-layer transition witnesses

This first step beyond the v0.33 composition of physically **disjoint**
crossing pairs studies **one connected three-source arrangement**: a
full-height through cylinder and two top-entry blind cylinders with
possibly different depths.

The resulting implementation is a first-party **read-only, source-replay
validated face-splitting plan**. It does *not* yet emit a canonical B-rep
with globally joined UV trim loops. It does not authorize a Boolean commit.

The v0.32 MultiCutRegionPlanner already returns complete exterior planar
circle-union arcs at each active-depth slab and computes an independently
verified material volume. The new ConnectedCutTransitionPlanner compares
the actual circle-owned angular intervals between neighboring slabs.

At every blind-floor transition, it partitions the circle arcs at every
existing lower/upper endpoint and labels each subinterval as:

- Continuous: physical cylinder side is exposed immediately above and
  below the floor.
- Terminating: the exposed lower cylinder boundary is covered by an
  newly active blind cutter immediately above the floor.
- Emerging: an upper exposed cylinder boundary first appears above the
  floor (including the new cutter).

This is a geometric transition witness, not a coedge allocated on the
B-rep. Every witness endpoint is a point on a source circle at the
physical floor height. It remains source-indexed, with a scalar radius,
world-space endpoints and angular interval.

The new void area at a depth transition is the difference between exact
Green-integral planar union areas in the upper and lower slabs. That area
is also the material-facing horizontal blind-floor area at the event.
All events require positive finite area and resolution above a disclosed
numerical budget.

A connected three-source boundary is required by the acceptance gate:
three different disk sources must appear on the boundary of a single
connected upper void. Separated two-circle groups plus an unrelated
circle do not pass as evidence for the new scope.

The plan always has can_apply=false. This is crucial: 3D surface splitting
will require true shared physical edge IDs, opposite coedges and
source-specific UV parameter curves on both incident faces. It must also
handle contour changes at the floor without crossing improperly between
circular patches. Until a canonical candidate passes manifold connectivity,
face sense, Euler/genus, mass and watertight triangulation checks, the
new plan remains non-executable.

## First-party regression fixture

For a 6×4×6 source block: one through disk centered (0,0), radius 0.75;
blind cutters centered (0.7,0), radius 0.75, depth 2.0 and (1.4,0.25),
radius 0.6, depth 1.0. All sources are upright and strictly interior.

The v0.32 planner returns 3 slabs and independent analytic volume;
v0.34 supplies two floor events with source-indexed emerging,
terminating or continuing circular wall sectors. Regression tests
validate each interval's endpoint against its circle, source-height
agreement, source replay integrity and refusal of a disconnected
three-source group or a tangency.

## Explicit next gate

The next separate PR should construct the entire 3D connected B-rep for
this one-through/two-blind bounded scenario. It must share vertices at
real floor intersection junctions, split the through and blind cylinders
along the global arrangement, allocate each edge exactly twice with
opposite coedges and closed loop cycles, and prove material orientation
at every face. Only after independent watertight display and source mass
checks should it allow source-atomic Undo/Redo.

This independent implementation uses no proprietary Siemens/Parasolid
code, hidden SDK structures or decompiled algorithms.
