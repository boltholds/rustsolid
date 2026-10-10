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

## v0.4: feature naming and incremental edit rollback

- `TopologyName`: deterministic feature-key + construction/operation path for reproducible, identically indexed feature histories. `NameChange` records parent-to-child provenance and the original name stays with the retained portion.
- `TopologyHandle`: entity ID, body token and unique process-local generation; prevents a handle from an aborted preview matching a new entity allocated at the same numeric slot. These handles cannot be persisted across sessions.
- `EditTransaction`: works on an exclusive borrow of Solid. Each operation records the pre-change local edge/coedge/face/loop/shell state and affected triangle slots; appended vectors are reverted by truncation. Drop, closure error or failed validation restores frames in reverse order without cloning the entire body. Structural validation still traverses the complete model per edit.
- Naming across arbitrary parametric regeneration, face merging/splitting matching, deletion of committed IDs, distributed edit conflicts and robust Boolean operators remain future capabilities. `TopologyName` alone does not prove geometric identity after topology changes.


## v0.5: Command + Memento geometry history

`CommandHistory` holds a validated `Solid`, a full applied `CommandBatch` log,
a bounded stack of inverse `EditMemento` values, and a redo stack of commands.
Commands are immutable and ordered; their selections reference `TopologyName`
rather than volatile numeric slots or process-local generational handles.

- `EulerCommand::SplitEdge` adds one vertex and one edge (`ΔV = ΔE = +1`).
- `EulerCommand::SplitFace` adds one edge and one face (`ΔE = ΔF = +1`).
- Undo performs the strictly scoped inverse deletion by replaying the **local
  inverse journal frames** in reverse. Each batch is an atomic checkpoint.
- Redo re-executes the saved commands and gives recreated entities *fresh
  handle generations*. Stale handles must not silently resurrect.
- Undo and Redo both increase the runtime revision monotonically. Semantic
  topological names remain deterministic if feature key, construction and
  command sequence remain identical.
- A new command after Undo invalidates the redo branch **only if successful**.
- `undo_limit` bounds in-memory Mementos; pruned entries remain in the command
  log as non-undoable history. `command_log()` excludes undone commands and
  can be serialized through serde and replayed against an equivalent primitive.

The `geometry.v1` stdio API includes an optional `command_history` request with
`source`, optional `feature_key`, `batches`, and optional `undo`/`redo` counts.
Its `solid.history` response includes `revision`, `undo_depth`, `redo_depth`,
and `applied_batches`. Older `block`, `extrude_profile` and `edit_solid` request
and response shapes remain backward compatible.

**Limits:** The two operations are Euler-characteristic-preserving splits, not
an arbitrary solid-building or independent deletion algebra. Their inverse
removals are currently valid only for the generated elements in an undoable
LIFO command history. Multi-shell/holes, robust general KEV/KEF, Boolean and
surface intersection remain unimplemented. Mementos and runtime handles are
process-local; durable state requires a regenerated source + serialized command
log (with future schema migration support). Full-model validation still runs
per operation; inverse journals are incremental in stored *data*, not yet in
validation complexity.

## v0.6: independent inverse Euler edits and analytic plane predicates

The first independent inverse operations are `Solid::kill_edge_vertex(vertex_id)`
(valence-two collinear vertex collapse) and `Solid::kill_edge_face(edge_id,
removed_face)` (merge two coplanar adjacent faces). They can remove a supported
split without reversing subsequent **unrelated** history commands. Live Euler
operators still work only for one closed genus-zero polyhedral shell; this is
not yet a general MEV/KEV/MEF/KEF calculus for holes and arbitrary shells.

Unlike the append-only `split_*` path, a kill operation currently rebuilds the
contiguous vertex/face topology and triangulation on a candidate body, validates
it, rebinds names with original provenance, and swaps it in on success. This is
an intentionally slower, whole-body transactional fallback until reusable
sparse slot arenas and local Euler edit proofs are implemented. Moving an ID
invalidates its old `TopologyHandle`; resolve its `TopologyName` to get a new
runtime handle. Untouched elements at unchanged IDs retain their handle. Killed
IDs are no longer addressable by name. No mesh or B-rep will be partially
mutated on a failed kill.

`EulerCommand` accepts make (split) and kill variants. If a batch contains a
kill operator, `CommandHistory` uses a guarded full-Solid Memento and executes
all commands on a working copy. Incremental inverse frames remain the default
for make-only batches. Undo of a kill restores the former geometry; redo
replays the named command. The command log remains the source of truth for
rebuilding after serialization. This path is intentionally not low-memory.

`predicates.rs` uses `robust` 1.2 (MIT OR Apache-2.0) for adaptive exact-sign
orientation tests with finite IEEE-754 coordinates, and provides independently
unit-tested plane/plane and segment/plane classifications. Constructed
intersection points are approximate f64; a line between infinite supporting
planes is **not** a trimmed face intersection. Surface-surface intersection
curves, NURBS, exact intersection constructions and arbitrary Boolean solid
splitting remain future work.

## v0.7: Geometry Foundation (first vertical slice)

The independent `geometry` module introduces explicit `Curve2`, `Curve3`,
`Surface3`, `Frame3`, `ParameterRange`, `ParameterDomain2`, typed carrier IDs,
UV face domains, and 2D trimming curves. Evaluated 3D curves include lines and
circles; evaluated analytic surfaces include planes, cylinders, and spheres.
Surface evaluation returns position and partial derivatives. Domains declare
periodicity (cylinder/sphere `u` is periodic), while 2D coedge trims have their
own bounded parameter interval.

**B-rep integration:** `GeometryStore::from_solid(solid, units)` (or
`solid.geometry_store(units)`) maps every *existing polygonal* B-rep edge to an
analytic line, every planar face to a planar supporting surface with projected
UV bounds, and each directed coedge to an oriented 2D trimming line on its
owning face's surface. It verifies edge endpoints and coedge UV-to-3D
reconstruction against scoped tolerances. Carrier IDs and bindings belong to a
revision- and body-token-bound **derived view**, not yet an owned, persistent
`GeometryStore` in `Solid`. After editing, rebuild the view; `check_current`
rejects stale revisions or bodies. Curved Surface3 evaluation is available,
but cylinder/sphere surfaces are **not yet bindable as curved B-rep faces** and
there is no true curve/surface intersection or Boolean.

`GeometryStore` supports per-edge and per-vertex tolerance overrides, validated
against each bound edge length. Underlying B-rep edit validation still uses
its model-wide `GeometryTolerance` until the full tolerance propagation layer
is implemented. Length units are explicit metadata (`Millimeter`, `Meter`,
`Inch`); analytic curve/surface parameters use radians and kernel model units.
No implicit rescaling is performed when selecting a unit; convert values
explicitly before constructing a body.

The `PredicateKernel` centralizes `robust` adaptive exact-sign 2D orientation
predicates and a separate model-tolerance classification. The legacy profile
builder, B-rep planar loop validator, and face editing now use the same robust
2D determinant function instead of separate naive cross products. Area
thresholds are still tolerance-based and there are remaining non-predicate
geometric operations to harden. Exact orientation signs do not imply exact
intersection coordinates, stable floating-point construction, or full 3D
intersection certificates.

Acceptance tests cover planar edge/face/trim associations, UV reconstruction,
analytic surface derivatives, scoped tolerances, stale revision detection,
large-coordinate translations and numerically challenging orientations.
The pre-existing `geometry.v1` JSON API remains unchanged.

Remaining work: permanent geometry ownership in sparse arenas, native curved
B-rep trims, arc/conic and rational B-spline support, bounded-surface
intersections, tolerance propagation through edits, stable geometry naming,
per-surface meshing, and face classification for Boolean modeling.

## v0.8 preview: first native curved B-rep (right circular cylinder)

The **separate** `CylindricalBrep` geometry type models an exact analytic right
circular cylinder, independently of tessellation. This first implementation
intentionally does **not** pretend to be a subtype of the older polygonal
`Solid`, whose face model requires >=3 distinct boundary vertices and exactly
two incident coedges on **different** faces. Those invariants would incorrectly
reject circles and periodic seam edges. Both types reuse the same typed vertex,
edge, coedge, loop and shell data structures, as well as `Curve3`, `Curve2`,
`Surface3` and their geometry bindings.

A complete capped cylinder has 2 vertices, 3 edges, 6 coedges, 3 loops,
3 faces and 1 closed shell (`V−E+F=2`). The two circular edges are closed and
have a single vertex each. The seam edge joins the two vertices and is used
twice by the same lateral face. The lateral UV trim is the rectangle
`[0, 2π] × [0, height]`: its two vertical sides map to one 3D seam edge via
**different** coedge pcurves at `u=0` and `u=2π`. Each cap has a single
closed circle coedge. The bottom cap frame is reversed, keeping all three
outward surface orientations consistent.

`CylindricalBrep::validate()` traverses directed coedge loops, checks twins,
incidence, UV/3D mapping at endpoints and interior samples, analytic normals,
periodicity, mass and axis-aligned bounds. The scoped validator is not a
mathematical proof against every possible geometric defect.

`CylindricalBrep::tessellate(segments)` produces a separate watertight indexed
mesh with `triangle_faces` preserving exact face ownership. It includes no
physical seam duplication, so vertices are welded at the UV seam in 3D. The
analytic cylinder radius, curvature and volume remain unchanged by facet
resolution. Facet volume is smaller than analytic volume; never use mesh
volume as the solid's authoritative volume.

The opt-in `geometry.v1` CLI operation `analytic_cylinder` returns the same
existing `solid` mesh fields, plus `solid.analytic_brep` describing the analytic
carrier relationships and periodic seam. All legacy requests remain unchanged.
This new curved model is not yet the canonical `GeometryStore` for other
solids and has no native Boolean, fillet, hole, serialization, history or
Euler edit operations.

### Mathematical definition

With a right-handed frame `(O, X, Y, N=X×Y)` and `r>0`, `h>0`:

- `S(u,v) = O + r (X cos u + Y sin u) + v N`
- `∂S/∂u = r (-X sin u + Y cos u)`; `∂S/∂v = N`
- Normal `normalize(∂S/∂u × ∂S/∂v) = X cos u + Y sin u`
- `u∈[0,2π]` periodic; `v∈[0,h]`
- `V=πr²h`; `A=2πr²+2πrh`

This is an independent implementation of public differential geometry, not
code or proprietary numerical internals copied from Parasolid.
