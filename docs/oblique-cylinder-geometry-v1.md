# RustSolid v0.39 — opt-in `geometry.v1` oblique-cylinder clip

This milestone exposes the already source-certified v0.38 bounded
cylinder/plane-halfspace Boolean through the existing stateless JSON
stdio bridge. It does **not** modify Gefest's production OCC path or
enable arbitrary B-rep/body Boolean operations.

## Request

```json
{
  "operation": "oblique_cylinder_clip",
  "cylinder": {
    "origin": [0.0, 0.0, 0.0],
    "axis": [-0.29552020666133955, 0.955336489125606, 0.0],
    "radial_x": [0.955336489125606, 0.29552020666133955, 0.0],
    "radius": 2.0,
    "height": 5.0
  },
  "plane_height": 2.5,
  "retain": "negative",
  "display_segments": 64,
  "undo": 0,
  "redo": 0
}
```

- `cylinder.axis` and `radial_x` are nonzero vectors that must be
  perpendicular to the model's angular tolerance. They determine a
  right-handed cylindrical chart; neither a tessellation nor Euler
  angles are used for the analytic source.
- `axis` must project **positively** onto world +Y. The finite
  plane is **horizontal in world space** at `y=plane_height`. It may
  therefore be oblique to the cylinder axis. This endpoint does not
  accept arbitrary plane normals or a second solid as its cutter.
- `retain="negative"` keeps the material below the plane;
  `retain="positive"` keeps the material above it.
- `display_segments` must be 12..=4096, default 64.
- `undo` may be 0 or 1, with `redo<=undo`. Both replay on a
  **per-request in-memory** editor and never persist between commands.
- The optional `tolerance` uses the existing `GeometryTolerance`
  fields (absolute_length, relative_length, angular) in model units.

A wide validated rectangular face on the cutting plane is generated
internally as a **temporary analytic trim witness**. The certified
pipeline is unchanged:

```text
CylindricalBrep source
    + finite plane witness
    -> IntersectionGraphBuilder (analytic 3D + dual UV)
    -> TopologyDeltaPlanner (source-bound preview)
    -> ClosedCylinderImprint (physical ellipse/coedge pcurves)
    -> ObliqueHalfspaceEditor (whole-face classification + selected side)
    -> CertifiedObliqueClipDisplay (manifold mesh and analytic mass)
    -> geometry.v1 JSON
```

The global `TopologyDeltaPlan::can_apply` remains false. Only the
specialized source-verified executor is allowed to commit its one
certified topology.

## Successful response

`ok=true` contains `solid` with `vertices`, `faces`,
`triangle_face_ids`, analytical `mass`, `bbox`, `topology`
and `analytic_brep`. The `analytic_brep.kind` is
`oblique_halfspace_clipped_cylinder`, with the ellipse carrier,
coedge owners, UV pcurve kinds, and a display certificate.

`data.kind="oblique_cylinder_clip"`; `data.state="committed"`;
`data.retained_halfspace` is the selected side. The 3-face
topology is a real, validated analytical B-rep. Display triangles
carry physical FaceId ownership; at `n` segments the mesh has
`4n` triangles and `2n+2` vertices. The facet volume differs
from the exact analytical volume and must not replace it.

A request with `undo=1,redo=0` returns
`state="imprinted_after_undo"`, four faces, `6n` display triangles
and the original source volume, with `topology_committed=false`.
Undo followed by Redo returns the clipped solid with revision 3.

## Fail-closed responses

Unresolvable intersections, tangency, near-cap events, coincident
supports and unsupported source geometry never publish a partial
`solid`. Such a response has `ok=false` and, for an expected
numerical ambiguity, `data.state="held"`,
`data.topology_committed=false`, and a typed stage/reason. Invalid
request geometry or budgets produce an ordinary `ok=false` error.

Every request is stateless from Gefest's perspective: no pre-existing
CAD scene, B-rep or persistent Undo/Redo history can be partially
modified by a failed native execution.

## Verification

```sh
cargo test --test oblique_cylinder_cli
cargo test --all-targets
cargo check --release --all-targets
cargo build --release --bin gefest-geometry
```

Four executable CLI tests check:
- positive/negative clips and independent source-volume/centroid sums;
- B-rep shells/vertices/faces and triangle owner multiplicities;
- undo/redo final state versus clipped/imprinted topology;
- invalid parameters, unsupported orientation, near-cap and
  out-of-bounds resolution returning no partial mesh.

The third-party Gefest repo exposes this operation **only** through a
separate opt-in backend route; the existing OpenCascade sketch/extrude
workflow remains the default. No Parasolid binaries or decompiled code
are used.

## Next gate

Cross-face region arrangements are still only locally sampled previews.
For any future multi-face or inner-loop Boolean, the engine must
construct complete face cells, classify each against the other source,
reconcile shared physical edges/coedges, integrate volume
independently, and pass atomic topology/mesh/Undo-Redo certification.
A multi-face preview cannot be promoted to `can_apply=true` on the
strength of a few point samples.
