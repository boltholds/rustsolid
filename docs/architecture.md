# RustSolid architecture

RustSolid is an independent geometric modeling kernel for Gefest CAD, distinct
from the Rust constraint-solving workspace. The production CAD path still uses
OpenCascade; the RustSolid JSON stdio interface is opt-in and backward
compatible with `geometry.v1`.

```text
Gefest sketches -- constrained by solver-core
         |
   geometry.v1 stdio / Rust library
         |
   RustSolid (model-space tolerance policy)
      |                  |
 Oriented B-rep      Polygon geometry
 Shell                Extrusion / Block
  Face                  |
   Loop             Triangulation + analytic mass
    Coedge             |
     Edge--Vertex    Selection mesh -> FaceId
```

## v0.2: oriented B-rep topology

Canonical incidence is expressed by typed handles: `ShellId`, `FaceId`,
`LoopId`, `CoedgeId`, `EdgeId`, `VertexId`. Collections use deterministic,
zero-based, model-local IDs. These are **not persistent topological naming IDs**
for re-generated solids. A `Face` stores canonical `loops` and a checked legacy
`boundary` vertex view. `Coedge` represents one directed edge usage by a face:
`edge`, `face`, `loop_id`, `reversed`, `next`, `prev`, and `twin`. An `Edge`
references exactly two oppositely oriented coedges. A `Loop` identifies its
first coedge and whether it is an outer or inner loop; the **current** builder
and validator support only one outer loop per face. A `Shell` owns its faces.

`Solid::check_topology()` returns the first `TopologyIssue`, with stable `code`
and typed `TopologyEntity` witness. It checks:

- ID/index agreement and in-range references;
- next/prev reciprocity, closed directed cycles and face boundary consistency;
- exactly two coedges per edge, reciprocal twins, opposite directions, distinct faces;
- one shell owning each face exactly once; a connected face-adjacency graph;
- a single incident coedge fan at each vertex (no pinched vertices);
- face planarity, nondegenerate area, and non-self-intersecting boundaries;
- genus-zero Euler characteristic `V-E+F=2`.

`Solid::validate()` additionally checks triangle/face ownership, mesh/B-rep
vertex coordinate agreement, triangle normals, manifold tessellation, signed
volume against analytic mass, and finite bounds. It uses the stored tolerance.
Mesh and B-rep share the same vertices in this slice; more general per-surface
tessellation will require its own traceability contract.

## Tolerance semantics

`GeometryTolerance` has three fields:

| Field | Default | Meaning |
| --- | ---: | --- |
| `absolute_length` | `1e-9` | Absolute model-space distance (in caller's length units) |
| `relative_length` | `1e-12` | Relative distance wrt local feature extent |
| `angular` | `1e-8` | Maximum angular deviation, radians |

The effective linear tolerance is `max(absolute_length, relative_length *
feature_extent)`; area tests use effective linear tolerance times feature
extent. Both tolerances and selected input dimensions are validated; a feature
smaller than its effective tolerance is rejected, not silently merged. A part
moved to large absolute coordinates still uses its *local* feature extent,
though sub-ULP detail cannot be recovered from floating-point inputs.

Call `Solid::extrude_xz_with_tolerance` or `Solid::block_with_tolerance` for
custom tolerances, or include an optional `tolerance` object with
`geometry.v1` input. Omitted fields receive the defaults. This is a **numerical
validation policy**, not an exact-predicate framework, tolerance propagation
system, or production-grade geometric robustness proof. Angular tests currently
apply to triangle/face normal agreement; more operations will extend its use.

## v0.3: topology edit transactions

- `Solid::begin_edit()` opens a single-writer transaction over a private clone.
  `Solid::edit_atomic(|tx| ...)` runs multiple changes atomically, revalidating
  the working B-rep and its tessellation after each successful operation.
- `EditTransaction::split_edge(edge_id, fraction)` inserts one vertex along
  `Edge.start → Edge.end`, subdivides the two opposite coedge usages and their
  triangle boundaries, retaining the original edge ID for the `Edge.start`
  segment and both original coedge IDs.
- `EditTransaction::split_face(face_id, a, b)` requires a simple planar face
  and two distinct, nonadjacent boundary vertices joined by an interior
  diagonal. It adds an edge, two coedges, one loop and one face. Existing
  boundary coedges and the original face/loop IDs remain in place; the two
  planar subfaces are retessellated. Exterior, touching or degenerate
  diagonals are rejected.
- `commit()` verifies all invariants before replacing the body and increments
  `revision` **once** for a nonempty batch. Failed edits poison their
  transaction; dropping without commit is an implicit rollback. The original
  body remains byte-for-byte equivalent at the Rust data-model level.
- `EditReport` carries typed `SplitEdge` and `SplitFace` provenance with source
  and created IDs; the stateless `geometry.v1` CLI accepts `edit_solid` and
  emits `edit_report` without changing legacy responses.
- Mesh validation additionally requires exact topological edge-to-triangle
  ownership and surface-area conservation. Mass/bounds remain unchanged by
  these geometry-preserving subdivisions.

### Composite versus graph

The *ownership* hierarchy `Solid → Shell → Face → Loop → Coedge` is
Composite-like. `Solid::topology_children()` returns only ownership children.
Edges and vertices can be referenced by many coedges, and `twin/next/prev`
form cycles: these are non-owning graph links addressed by typed handles.
Treating the full B-rep as a recursive object-tree would duplicate shared
geometry or cause cyclic ownership. `Solid` acts as the aggregate root for
transactions, not a universal polymorphic `Component` class.

ID preservation is local to these append-only edit operations. Rebuilding
solids from feature history has no persistent naming algorithm yet.
Transactions presently clone the full solid and use O(n) validation; an
incremental journal/arena and revisioned handles are planned after correctness
and provenance fixtures are stable.

## Unsupported and safety limits

Currently only one closed, connected genus-zero polyhedral shell with one
outer wire per face is supported: no internal wires/holes, multiple shells,
shared non-manifold edges, analytic/NURBS surfaces, Boolean, fillet/chamfer,
STEP/X_T or persistent naming after edits. `FaceRole` carries extrusion
provenance only. The geometric validity checker does not yet prove arbitrary
3D face-face nonintersection. Never feed unvalidated solids into manufacturing.

## Planned gates

1. Extend the topology transaction editor with deletion, local face rewiring,
   robust invalid-handle diagnostics and partial journal-based undo.
2. Planes, cylinders, conics, parametric surface domains and trim curves.
3. Robust adaptive/exact orientation and intersection predicates with
   differential tests against license-compatible geometry libraries.
4. Half-edge subdivision operations, intersection traces, face splitting,
   Boolean operations, fillets and topological provenance.
5. Tessellation with per-face parametric traceability, STEP import/export,
   benchmarks and integration with Gefest's feature history.

## IP boundary

No proprietary Siemens/Parasolid code, binaries, headers or decompiler outputs
may be committed to RustSolid. Implementation is independent and based on
standard computational geometry concepts; any future third-party code must
have a compatible license and documented provenance.
