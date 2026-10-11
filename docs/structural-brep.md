# RustSolid v0.23 — exact genus-N structural B-rep transactions

This milestone executes a bounded family of **true topology-changing**
analytic operations within the same transaction engine used to edit curves
and UV trims. A certified six-face axis-aligned box may receive, remove and
resize up to 32 independent, fully penetrating circular holes with +Y axes.
The result is a single connected, orientable closed B-rep of genus N, where
N is the number of remaining holes.

It does **not** authorize arbitrary curved `TopologyDeltaPlan` application.
It is the first independently testable mechanism for creating and retiring
analytic topological entities and recomputing mass while retaining
source-validated, atomic carrier reparameterization.

## Canonical multi-hole B-rep construction

Each real through-hole is an analytic subtractive feature, carrying:

- Two real circle edges, each with opposing coedges shared by one cap's
  inner loop and the hole's cylindrical side.
- A third seam edge with two oppositely directed coedges on the *same*
  periodic cylindrical face, distinguished by separate UV pcurves at
  `u=0` and `u=2π`.
- Two new inner circular UV trim loops on the top and bottom planar caps.
- One new periodic cylindrical supporting surface/face with its own
  `v=[0,height]` interval and physical 3D circle and seam carriers.
- Two new topological vertices at the circle's chart seam endpoints.

The original box retains six faces and its own 3D/UV carriers. A complete
candidate for N holes has:

| Canonical topological quantity | Formula |
| --- | --- |
| Vertices | 8 + 2N |
| Edges | 12 + 3N |
| Coedges | 24 + 6N |
| Faces | 6 + N |
| Loops | 6 + 3N |
| Inner trim loops | 2N |
| Shells | 1 |
| Euler characteristic | 2 − 2N |
| Genus | N |

All holes must be strictly inside the four vertical box sides, remain
separated by more than the combined radii plus model-space tolerance,
and span exactly the full box height. Intersecting, tangent, touching,
oblique, partial-depth, blind and sub-tolerance holes are refused. The
32-hole limit bounds reconstruction work and input complexity.

### Independent analytic mass contract

Let the original box have volume `V0`, area `A0`, centroid `C0`,
height `h`, and each hole have radius `ri` and cross-sectional center
`(xi,zi)`. Then:

`V = V0 − πh Σ(ri²)`

`A = A0 − 2π Σ(ri²) + 2πh Σ(ri)`

The centroid is calculated from the original center and *relative removed
moments*: `C = C0 + (Σ Vi (C0 − Ci))/V`, with
`Vi=π ri² h` and `Ci=(xi, C0.y, zi)`. This avoids subtracting two
huge absolute moments for parts far from the world origin.

These properties are recomputed on each complete structural rebuild.
The result's geometric topology, genus and UV/3D carrier correspondence
are validated by the shared `BrepModel::validate` contract. Triangle
meshing is not part of the mathematical Boolean boundary.

## Unified atomic editor

`BrepStructuralEditor::new(Solid)` or
`BrepStructuralEditor::from_through_hole(&ThroughHoleBrep)` starts a
private candidate model with process-local generational hole slots.
`BrepArena::detached_structural_editor()` is an adapter for a verified
box or single-hole arena; it deliberately leaves the arena's original
Euler and feature history untouched.

`BrepStructuralTransaction` supports:

- `add_through_hole(HoleParameters) -> StructuralHoleHandle`
- `remove_through_hole(StructuralHoleHandle)`
- `resize_through_hole(StructuralHoleHandle, HoleParameters)`
- `reparameterize_edge(StructuralEdgeRef, scale, offset)`
- `reparameterize_trim(StructuralTrimRef, scale, offset)`
- `rebase_face_uv(StructuralFaceRef, du, dv)`

Edge/trim/face targets may address either a stable base-box topology ID
or an existing/newly staged hole handle and an explicit carrier role.
The hole reference identifies its feature across dense topological
renumbering; it is not an unsafe numeric approximation to a 3D curve ID.

At commit the editor:

1. Checks all current holes and strict pairwise separation.
2. Reconstructs the complete canonical B-rep, **creating or removing**
   real topological and analytic carrier entities.
3. Replays surviving shape-preserving carrier/UV intents through the v0.22
   `BrepCarrierEditor`. Hole carrier references are rebound using their
   stable feature slot; a removed hole's associated intents are retired
   with it. Replay is bounded to 2048 intents, applied in private,
   256-command carrier batches.
4. Verifies all manifold, UV, seam, face-role, Euler/genus and carrier
   invariants; calculates the new analytic mass and centroid.
5. Publishes the complete `State` in **one revision**, producing one
   Undo/Redo Memento, or rejects it with **no externally visible changes**.

Ignoring a failed stage poisons its entire transaction. A late invalid
curve parameterization after valid hole insertions similarly cancels the
whole batch. The old candidate, revision, handles, Undo/Redo and cached
geometry remain unchanged on failure.

## Topological identities and undo safety

A `StructuralHoleHandle` has `editor_token`, sparse `slot`, and
monotonically minted `generation`. Each binding reports the current
`StructuralHoleTopology`: cylindrical face, two circular edges, seam edge,
and three loops. A surviving handle retains the same incarnation even if
another feature is removed and its dense indices compact.

Abandoned staged handles and removed handles are stale. Undo/Redo never
resurrect previously retired generations: a restored feature gets a fresh
generation, while features surviving the transition retain theirs.
The commit reports created, retired and remapped feature bindings.
A revision-bound `StructuralSnapshot` is invalid against later revisions.
Undo history is bounded to 128 committed steps.

Source `BrepArena` handles are not reused in the detached structural
editor, and neither the public process-local feature slots nor carrier
generations are durable serialized CAD feature IDs.

## Rust usage

```rust
use rustsolid::{
    BrepStructuralEditor, HoleParameters, HoleEdgeRole, StructuralEdgeRef,
    Point3, Solid
};
let block=Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let mut editor=BrepStructuralEditor::new(block)?;
let mut tx=editor.begin_transaction();
let hole=tx.add_through_hole(HoleParameters {
    center_x:-1.5, center_z:0.0, radius:0.5
})?;
tx.add_through_hole(HoleParameters {
    center_x:1.5, center_z:0.0, radius:0.7
})?;
tx.reparameterize_edge(
    StructuralEdgeRef::Hole{hole,role:HoleEdgeRole::BottomCircle},
    1.0,0.4
)?;
tx.commit()?;
assert_eq!(editor.model().summary().genus,2);
editor.validate()?;
```

## geometry.v1: structural_transaction

```json
{
  "operation": "structural_transaction",
  "source": {"kind":"block","origin":[-3,-1,-3],
             "width":6,"height":4,"depth":6},
  "batches": [[
    {"kind":"add_through_hole","center":[-1.5,0],"radius":0.5},
    {"kind":"add_through_hole","center":[1.5,0],"radius":0.7},
    {"kind":"reparameterize_hole_edge","slot":1,
     "role":"bottom_circle","scale":1,"offset":0.4}
  ]],
  "undo": 1,
  "redo": 1
}
```

The response returns `data.kind="structural_brep"`,
`topology_committed:true`, current genus, exact analytic volume/area/
centroid, curve and trim carrier counts, created/retired/remapped hole
bindings, revision and history transitions. All operations in each inner
`batches` array commit atomically; slots in JSON identify features
within the current request/session.

The multi-hole result is a fully validated, queryable canonical analytic
`BrepModel` through the Rust API. **General multi-hole display
tessellation is not available in this milestone**, and the JSON operation
sets `display_mesh_available:false` without inventing a partial mesh or
discarding true circular boundaries. Direct `inspect_brep` of a multi-hole
body is not yet a supported standalone `BrepInput` representation.

The existing `boolean` endpoint continues to execute v0.19 planar box
CSG and the v0.21 *single* analytic through-hole. The structural editor
is the explicit path to more than one hole, not a general Boolean kernel.

## Remaining engineering gates

1. A standalone heterogeneous mutable owner store with generic insert,
   split, merge and retire operators for arbitrary nonconcentric analytic
   carriers; no box-only source assumption.
2. General face-region construction from **all** pairwise analytic
   intersection graphs, including intersections of multiple curved
   surfaces and complete determination of inside/outside/touching regions.
3. Multi-shell and multi-component B-rep topology, cavity ownership,
   arbitrary genus and rigorous oriented-face/inner-loop validation.
4. Closed display tessellation of arbitrary multi-hole planar caps and
   curved side walls as a derived cache, followed by STEP-level exchange.
5. Independent mass integrals for arbitrary trimmed analytic surfaces,
   failure-injected generalized Boolean commit and genuine persistent
   feature lineage across regeneration.

The global `TopologyDeltaPlan::can_apply == false` guard remains.
For engineering production, the existing Gefest OpenCascade path remains
the recommended default until equivalent broad coverage is demonstrated.

Eighteen independent first-party regressions include insertion-order
metamorphism, genus and mass invariants, different source/hole counts,
failed partial commits, ABA on Undo/Redo, repeat-stage budget exhaustion,
and references surviving dense topology compaction.

Tests: `tests/structural_editor.rs` and the dedicated mandatory
`geometry.v1` CI smoke/atomic-refusal gates.
