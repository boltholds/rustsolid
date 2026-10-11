# RustSolid v0.26 — composition of analytical through and blind cuts

RustSolid now builds a **single true analytic B-rep** containing a mixed
sequence of strictly disjoint vertical cylindrical subtractive features.
Each cut may penetrate the full box or stop on an independent internal
circular planar floor. Composition uses a common topology/trim carrier
assembly layer, exact mass moments and one closed shared shell.

The model **does not** subtract polygonized triangle meshes, nor quietly
allow all arbitrary analytical Booleans. This is a bounded, executable
step toward general face-region reconstruction, with explicit refusal of
intersecting/touching cuts or unsupported source geometries.

## Geometric input contract

Supported left operand: a fully validated, exact six-face axis-aligned
box. Supported right-operand family: up to **32** upright circular
cylinders, each of which independently certifies either:

- **Through:** exactly spans the box between its lower and upper
  horizontal end faces and lies strictly within the four side planes.
- **Blind:** opens precisely on the top (+Y) plane and terminates on a
  horizontal planar floor strictly inside the material above the
  original box bottom.

All circular openings must be pairwise separated in the top XZ sketch
by more than their combined radii plus local tolerance. Thus two
cutters cannot overlap or touch each other, even if their depths differ.
This avoids an unsupported cylindrical-cylinder intersection, merged
or broken contour graph and material-boundary ambiguity.

For **each blind cut**, the source block top plane and the original
cylindrical cutter side must generate one independently revalidated,
closed transverse `IntersectionGraph` opening. The same proof is
performed when cuts are supplied as typed `AnalyticCut` parameters
rather than pre-built `CylindricalBrep` source operands.
Through-cut plane alignment follows the existing strict v0.25
float-representability rules; no arbitrary rotation is snapped into a
vertical shape.

## Canonical common-topology assembly

`AnalyticCut` distinguishes the existing exact, validated
`HoleParameters` (through) and `BlindHoleParameters` (blind)
without flattening their semantic difference to a mesh.

The v0.26 B-rep constructor:

1. Certifies the exact source box, validates all input features and
   checks their pairwise spatial separation.
2. Builds an independent canonical **single-cut** analytical B-rep
   through the already validated v0.21/v0.25 fragment producers.
3. Rebinds each fragment's true vertices, circular and seam edges,
   two directed uses per edge, coedge UV curves, inner/outer loops,
   cylinder support and optional floor support into one shared
   `BrepModel` geometry catalog.
4. Keeps the **same original six box faces** and their original
   edge/vertex owners. Each through-cut attaches an inner loop to
   *both* caps. Each blind-cut attaches an inner loop only to the
   upper cap, introduces a new separate circular floor face and
   leaves the original bottom cap untouched.
5. Remaps all fragment-owned coedge `next/prev/twin` pointers and
   loops to new dense IDs; preserves the original support-surface
   parameterization and explicit `FaceSense::Reversed` on every
   subtractive cylindrical wall/floor.
6. Recomputes the **entire resulting body's** exact volume, surface
   area and centroid from independently derived removed-cylinder
   moments. Rejects any failed manifold, UV/3D carrier, Euler,
   face-owner, face orientation or analytical mass check.

Every cut adds **two** actual seam vertices, **three** analytic edges,
**six** directed coedges and **three** topological loops. A blind cut
has two additional faces (cylindrical side and floor), while a
through-cut has one (cylindrical side). For `t` through and `b`
blind cuts:

| Quantity | Exactly |
| --- | --- |
| Vertices | `8 + 2(t+b)` |
| Edges | `12 + 3(t+b)` |
| Coedges | `24 + 6(t+b)` |
| Faces | `6 + t + 2b` |
| Loops | `6 + 3(t+b)` |
| Inner circular face loops | `2t + b` |
| Shells | `1` |
| Euler characteristic | `2 - 2t` |
| Genus | `t` |

This makes a critical distinction: **blind holes do not increase
genus**. They introduce floor disks; through-holes create tunnels.

### Exact analytical mass

For a source box of height `H`, with volume `V0`, area `S0`,
centroid `C0`, and a cut `i` of radius `ri` and removed depth
`di`, the removed volume is `Vi = πri²di`.

```
V = V0 - Σ Vi
S = S0 + Σ through (2π ri H - 2π ri²)
       + Σ blind   (2π ri di)
C = C0 + (Σ Vi (C0 - Ci)) / V
```

Here the removed volume centroid `Ci` lies at the cutter's
`(cx, cz)` and either the box height-midpoint (through) or
`topY - di/2` (blind). Calculating relative moments improves
large-world-coordinate stability.

This exact mass is **never** derived from the display triangles.

## Mixed cap regions and certified viewport tessellation

`CompositeRegionMesher` reconstructs two *different* planar regions
from the exact canonical 3D edge carriers:

- **Upper cap:** one exterior rectangular loop minus **all**
  circular openings.
- **Lower cap:** the exterior rectangular loop minus **only
  through-cut** circular openings.
- **For each blind cut:** an additional circular floor disk at the
  true cutter termination plane.
- **For each cut:** the analytically supported inward-facing cylinder
  wall, sampled as facets sharing both ring vertex arrays with the
  cap and/or floor.

Only the display face regions are triangulated with ISC-licensed
`earcutr`. Each mesh is checked against the current B-rep for
face ownership, correct analytical support planes/radii, normal
directions and material `FaceSense`; each undirected facet edge
must appear exactly twice with opposite directions. The mesher checks
cap area, tetrahedral signed volume and surface area against a
**separate polygonal source oracle**:

```
Vfacet = V0 - Σ Apolygon_i di
Sfacet = S0 + Σ through (Ppolygon_i H - 2 Apolygon_i)
            + Σ blind   (Ppolygon_i di)
```

The inscribed polygon removes slightly less than its true analytical
circle, so the faceted material volume should approach the exact
analytic volume from above as segment resolution increases.

The display supports 12–512 segments per circle and at most 16,384
sampled circumference vertices. Segmentation changes neither the
source B-rep, its genus, real `Curve3::Circle` edges nor mass.

## Revisioned strongly atomic editing and identity

`BrepCompositeEditor` owns the exact source box, typed cuts,
source-bound canonical B-rep, monotonic revision and undo/redo history.
Its transaction can `add`, `remove` or `replace` a cut's type,
radius, center and/or blind depth, all in one validated **atomic
commit**. Invalid staging poisons the entire batch even if its error
was ignored. All geometry construction occurs on a private candidate;
failures cannot publish partial faces, carriers or mass.

`CompositeCutHandle` has an editor token, sparse feature slot and
generational incarnation. A surviving cut keeps the same handle if
an earlier cut disappears and its dense topology indices compact.
Removed or abandoned incarnations cannot be used again; Undo/Redo
cannot accidentally revive them, even when the same slot is replaced.
Cloning an editor creates an independent namespace: a handle belonging
to its original must be rejected as foreign in the clone.

The commit reports created, retired and remapped feature bindings
containing the currently owned cylinder face, optional floor face,
circle/seam edges and loops.

A `BrepCompositeEditor` is also a `BrepSource`, accepted through
`BrepBody::Composite` for common B-rep queries, full-region
membership and meshing. `BrepArena::from_composite` permits an
independent read-only Arena view, and
`BrepArena::detached_composite_editor` creates an editable candidate
from supported original box/through/blind/composite sources without
modifying the originating Arena.

The public serializable stable cross-session CAD feature identity
problem **remains open**. Process-local handles are not durable source
identifiers.

## Rust Boolean composition

```rust
use rustsolid::{
    BrepBody, BooleanKernel, BooleanOperation, BooleanOutcome,
    CylindricalBrep, GeometryTolerance, Point3, Solid,
};
let box_=BrepBody::from(Solid::block(
    Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?);
let cutters=vec![
    BrepBody::from(CylindricalBrep::upright(
        Point3{x:-1.5,y:-1.,z:0.},0.6,4.,
        GeometryTolerance::default())?),
    BrepBody::from(CylindricalBrep::upright(
        Point3{x:1.4,y:1.5,z:0.3},0.75,1.5,
        GeometryTolerance::default())?),
];
let result=BooleanKernel::apply_sequence(
    &box_,&cutters,BooleanOperation::Difference)?;
if let BooleanOutcome::Composite(editor)=result.outcome {
    assert_eq!(editor.model().summary().genus,1);
    let viewport=editor.tessellate(64)?;
    assert!(viewport.faceting_volume_error>=0.0);
}
```

`apply_sequence` only admits `BooleanOperation::Difference`, a
certified exact box on the left, and one to 32 fully validated
cylindrical B-rep cutters. It does **not** turn arbitrary solids into
cylinders by inspecting their bounding boxes. Each cutter is first
certified to be either a through-hole or a top-entry blind hole, then
inserted into the shared canonical geometry. The old one-cutter
`BooleanKernel::apply` behavior remains unchanged.

`BooleanResult.face_sources` records source operand left/right,
source face (cylinder side=2, floor cap=0), result face and reversed
orientation. The composite editor's individual bindings give the
complete cutter-to-feature correspondence.

## geometry.v1

Execute a **one-shot** analytical Boolean from real cylinder inputs:

```json
{
  "operation":"boolean_sequence",
  "mode":"difference",
  "left":{"kind":"block","origin":[-3,-1,-3],
          "width":6,"height":4,"depth":6},
  "right":[
    {"kind":"analytic_cylinder","origin":[-1.5,-1,0],
     "radius":0.6,"height":4},
    {"kind":"analytic_cylinder","origin":[1.4,1.5,0.3],
     "radius":0.75,"height":1.5}
  ],
  "segments":48
}
```

Or use an **atomic feature history** session with typed cuts:

```json
{
  "operation":"composite_transaction",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "batches":[[
    {"kind":"add","cut":{"kind":"through",
      "center":[-1.5,0],"radius":0.6}},
    {"kind":"add","cut":{"kind":"blind",
      "center":[1.4,0.3],"radius":0.75,"depth":1.5}}
  ]],
  "display_segments":48
}
```

Additional commands are `{"kind":"remove","slot":0}` and
`{"kind":"replace","slot":1,"cut":{...}}`. The session also permits
`undo`/`redo` counts. `data.cut_count`, `through_cuts`,
`blind_cuts`, genus, exact mass, topology, stable feature generation,
per-cut side/floor lineage, created/retired/remapped histories and
`display_quality` are returned. When the optional
`display_segments` is omitted, the transaction returns *no*
`solid` payload, preserving the explicit distinction between
analytic topology and display tessellation.

Both `inspect_brep` and `geometry_query` now accept the input
`kind:"analytic_composite_cuts"` with box dimensions,
`cuts:[{kind:"through",...},{kind:"blind",...}]` and optional
tessellation segments.

An unsupported or invalid cut sequence, failing UV/3D topology,
inconsistent analytic mass, intersecting cut or failed mesh
certificate produces `ok:false` and **no partial solid or result
data**.

## Acceptance and important limitations

`tests/composite_cut.rs` pins 22 independent regressions: mixed
genus/mass/centroid, input-order metamorphism, true face/UV ownership,
planar floor loops, signed/closed display manifold and its volume
convergence, obsolete/foreign/abandoned handles, clone isolation,
atomic batch rollback, query/Arena compatibility, 32 mixed-cut stress,
source-bound Boolean operand lineage, 32-feature stress, deterministic\nmesh replay, and 33rd-feature refusal. Mandatory GitHub Actions
gates run both this module and independent CLI positive/failure cases,
alongside all earlier v0.25 tests.

**Remaining limits**: This version does **not** execute a general
`TopologyDeltaPlan`, splice crossing cylinders, merge intersecting
openings, create lateral or tilted blind holes, build non-AABB
multi-shell bodies or calculate NURBS intersections. The composition
mechanism is reusable across these two supported cut types, but the
fragment construction still relies on certified **canonical box +
circular vertical cylinder** source families.

The next architectural gate is a general oriented intersection-region
arrangement and mutable face/edge/UV topology on arbitrary supported
analytical operands, with new independent mass integration and
source-persistent feature naming.
