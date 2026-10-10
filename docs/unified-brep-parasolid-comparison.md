# RustSolid v0.9 — unified topology/geometry graph

This first **read-only shared B-rep** is an explicit transitional layer over the
existing editable polyhedral `Solid` and the exact `CylindricalBrep`. A `BrepSource`
trait (with a convenience `BrepBody` enum) maps both into `BrepModel`, preserving
native 3D curves, surfaces, and their separate 2D face trims. No faceted
approximation is substituted for an analytic cylinder in the geometry graph.

## Shared model

- `Vertex`, `Edge`, `Coedge`, `Loop` and `Shell` keep their existing typed IDs.
- `BrepFace` owns loop references and a `FaceGeometry { surface, bounds, ... }`.
- `BrepGeometry` is the shared catalog of `Curve3`, `Surface3`, `Curve2`, with
  edge/face/coedge carrier bindings and explicit parameter ranges.
- Read-only queries cover `face_surface`, `edge_curve`, `coedge_pcurve`,
  `next_in_loop`, `next_of_edge`, `coedges_of_edge`, `loops_of_face`.
- One manifold edge has two directed uses. A circular edge may have coincident
  start/end vertices. Two twin coedges may belong to **the same periodic face**,
  in which case distinct UV representations are required.
- The shared checker verifies indices/ownership, closed cycles, reciprocal
  coedge links, radial twins, connected shell, Euler parity, and sample-based
  agreement of a 3D edge curve with each `Surface3(Curve2(t))` trim.
- The number of visualization triangles is independent of the analytic graph.
  Geometry/evaluation and display meshing have separate interfaces.

## Scope and remaining gaps

The common model **does not** yet replace source editing or own an independent
transactional arena. `Solid::edit_atomic()` and the Euler/Command/Memento path
still operate on polyhedra. The curved cylinder has its independent producer and
validator. `BrepModel` is created afresh from its source; it does not serve as
a durable cache or as a new persistent-topology-naming implementation.

Only one closed, connected orientable manifold shell with two coedges per edge
is supported. The Euler genus check is a *necessary* global invariant, not a
complete proof against pinched vertices, self-intersection or ambiguous trim
containment. General radial cycles, open-sheet bodies, inner trim loops, holes,
exact NURBS, curved Euler editing, face intersections and Boolean are future
work. UV/3D consistency is checked at sampled parameters and is not a formal
certificate for arbitrary parametric curves.

## Functional architecture comparison with Parasolid 38.01.207

The user-supplied decompilation is a proprietary binary-derived **analysis
reference**, not public source. These are observed *API capabilities* and
subsystem boundaries, not copied implementation algorithms:

| Parasolid API signal | RustSolid current equivalent | Gap |
|---|---|---|
| `PK_FIN_ask_next_in_loop`, `PK_FIN_ask_next_of_edge` | common directed loop traversal and two-use edge traversal | arbitrary radial fin cycles |
| `PK_FACE_ask_oriented_surf`, `PK_EDGE_ask_curve` | face/edge carrier IDs in shared graph | general surface/edge representation and orientation |
| `PK_SURF_eval`, `PK_CURVE_eval` | analytic `Surface3`/`Curve3::evaluate` | derivatives beyond first order, conics, NURBS |
| `PK_FIN_ask_curve` / geometry | per-coedge `Curve2` trim and parameter range | robust complete trim topology and reparameterization |
| `PK_EDGE_ask_precision` | `GeometryTolerance` per body, edge bindings | tolerance propagation and certification |
| `PK_SESSION_ask_general_topology` | closed two-use manifold shell | sheets, generalized non-manifold edges |
| `PK_SURF_intersect_surf` / `PK_BODY_boolean` | only plane/plane support and isolated Euler operations | clipped surface intersections, classification, Booleans |
| `PK_MARK_*`, `PK_PARTITION_*` | Command/Memento + per-body journal | partition-level transactions and incremental persistent I/O |

Inspection of representative wrapper pseudocode shows the API layer performs
substantial validation, error/status handling and delegation into internal
helpers. Function lengths, or the presence of a symbol, cannot prove the
internal intersection or Boolean method is fully recoverable.

The underlying standard mathematics of curves/surfaces, parameterization and
computational topology can be independently implemented from open references.
Do not include Siemens source fragments, type layouts, internal helper code,
or large decompiler-derived translations in this repository. New algorithms
should have reproducible correctness tests and separately documented provenance.
