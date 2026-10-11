# RustSolid v0.26 — exact composition of separate blind-cut face regions

This milestone provides the first bounded multi-cut assembly of analytic
plane/cylinder B-reps: a certified axis-aligned block minus **1–16 pairwise
separated, strictly vertical top-entry blind cylindrical bores**, with
different radii, depths and cross-sectional centers.

It is a true analytic B-rep **builder** rather than a general-purpose
`TopologyDeltaPlan` executor. Each source blind bore independently
passes the v0.25 closed-opening `IntersectionGraph` witness and
comprehensive canonical B-rep validation before any topology is grafted.
No display mesh is used to define the cut boundaries.

A candidate is built privately from the first certified bore. Each later
bore contributes its two circular edges, seam edge, six coedges,
three loops, two seam vertices, analytic cylinder side face, and planar
floor face. Its top circular cap loop is linked into the **same existing
top block face** as a separate inner contour. The graph retains distinct
`Curve3`, `Curve2` and `Surface3` carriers, including each periodic
cylinder's twin coedges and material-facing reversed `FaceSense`.

The source boxes must be exact six-face blocks, the bore openings must
remain strictly inside all four box walls, the individual cutting cylinders
must end inside the material, and all hole footprints must be separated by
more than the two radii plus local tolerance. Intersecting, touching and
near-tangent circles are rejected even if their bore depths differ.

For N holes the canonical graph has:

| Topological property | Count |
| --- | ---: |
| Vertices | 8 + 2N |
| Edges | 12 + 3N |
| Coedges | 24 + 6N |
| Faces | 6 + 2N |
| Loops | 6 + 3N |
| Top cap inner loops | N |
| Independent planar floors | N |
| Connected closed shells | 1 |
| Euler characteristic | 2 |
| Genus | 0 |

The exact analytical mass is independently recomputed using removed
cylinder volumes `Vi=π ri² di`, side-wall areas `2π ri di`,
and centroids `Ci=(xi, topY-di/2, zi)`:

`V = Vbox - Σ Vi`

`S = Sbox + Σ (2π ri di)`

`C = Cbox + Σ [Vi (Cbox - Ci)] / V`

Every resulting `BrepModel` is fully validated, including topology,
genus, per-coedge UV/3D correspondence, analytic carrier ownership,
and the canonical physical face orientation contract.

## API

Rust:

```rust
let model = CompositeBlindBrep::build(&block, &[
    BlindHoleParameters {center_x:-1.5,center_z:0.0,radius:0.6,depth:1.2},
    BlindHoleParameters {center_x: 1.5,center_z:0.0,radius:0.7,depth:2.4},
])?;
model.validate()?;
assert_eq!(model.summary().faces,10);
assert_eq!(model.summary().genus,0);
```

JSON `geometry.v1`:

```json
{
  "operation":"compose_blind",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "holes":[
    {"center":[-1.5,0],"radius":0.6,"depth":1.2},
    {"center":[1.5,0],"radius":0.7,"depth":2.4}
  ]
}
```

The API returns `data.kind="composite_blind"`, a verified analytic
topology summary, exact mass and a list of floor face IDs (for this
example, `[7,9]`). Errors are fail-closed and never return partial
B-rep or mass results.

## Honest scope

- No generic intersecting or oblique cylinders and no multi-shell cavities.
- No persistent parametric feature history/Undo/Redo for the entire
  composite object yet: the result is a newly owned **immutable analytic
  B-rep**. The v0.25 *single*-blind bore retains its separate parameter
  transaction engine.
- No composite display tessellation yet; the JSON operation returns
  topology/mass data, not a pretend mesh. The source-verified single-blind
  and independent multi-through display meshes remain available.
- Huge absolute coordinate systems may be rejected by the existing
  v0.25 source intersection tolerance/graph validator. The composition
  does not bypass that certificate or advertise unrestricted
  large-coordinate support.
- General `TopologyDeltaPlan::can_apply` remains false, and the
  `BooleanKernel::apply` entry point is unchanged for generic multiple
  cuts. This new operation is explicitly `CompositeBlindBrep::build`.

The next correct step is a general heterogenous region arrangement that
can insert multiple blind and through cuts in a shared mutable owner graph,
with persistent face lineage, meshing of several circular floors, and
fully source-certified edge splitting at crossing curves.
