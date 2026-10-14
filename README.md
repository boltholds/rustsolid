# RustSolid

RustSolid is an **independent Rust solid geometry kernel** being developed for [Gefest CAD](https://github.com/boltholds/Gefest-CAD).

> **Status: v0.39 experimental kernel:** v0.38's certified oblique cylinder/halfspace Boolean is now exposed through a narrowly typed `geometry.v1` JSON operation with analytic face ownership, history replay and fail-closed errors. The v0.35 connected triple-cylinder Boolean remains available separately. **General Boolean between arbitrary B-rep bodies remains unsupported.**

## Current capabilities

- **v0.39: all-or-nothing multi-face region preview.** A source-namespaced
  `MultiFaceRegionPlan` aggregates bounded geometric intersection graphs and
  local arc-side Boolean evidence across multiple finite face pairs. Any
  ambiguous pair, near tangency, or unsupported inner UV loop holds the
  *entire* plan with no partial region report; general Boolean mutation
  remains disabled. [Multi-face preview contract](docs/multi-face-region-preview.md).

- **v0.39: opt-in `oblique_cylinder_clip` JSON endpoint.** A tilted analytic cylinder is cut by a bounded world-horizontal plane witness. Two retained halfspace results produce genuine elliptic-cap B-reps and watertight, per-face-owned display meshes with analytic mass/centroid. One-shot Undo/Redo, typed `Held` and tolerance budgets are exposed without changing the default Gefest/OpenCascade backend. [v0.39 integration contract](docs/oblique-cylinder-geometry-v1.md).

- **v0.38: certified oblique imprint display and bounded halfspace Boolean.**
  A single canonical physical elliptical edge is sampled once, shared by
  upper/lower cylinder-side faces, and validated as a four-face watertight
  display with exact physical face ownership. The closed intersection divides
  the source cylinder into analytically classified upper/lower regions. A
  narrow first-party volume-changing operator retains either side and builds
  a closed three-face body with its own elliptical planar cap, exact-in-form
  volume, first moments, bounds and certified face-owned mesh. Transactions
  fail closed with atomic Undo/Redo and independent disk-quadrature tests.
  This is a cylinder/halfspace contract, **not generalized body/body CSG**.
  See [v0.38 analytic clipped-cylinder milestone](docs/oblique-imprint-mesh-halfspace-clip.md).

- **v0.37: executable closed oblique cylinder-side imprint.** Converts a verified
  closed plane/cylinder trace into a true physical spatial elliptical edge,
  two opposite coedges with independently stored analytical 2D cylinder
  trims and a two-face partition of the cylindrical wall. One source-bound
  topology-intent plan is checked before committing; invalid/near-cap cases
  are Held. The model is atomically validated with monotonic Undo/Redo and
  unchanged source mass. This is a **shape-preserving imprint**, not a
  volume-changing Boolean, and the original display mesh does not yet carry
  the two new face owner IDs. See [v0.37 closed oblique cylinder imprint](docs/closed-oblique-cylinder-imprint.md).

- **v0.36: paired analytic pcurve witnesses for oblique plane/cylinder
  face intersections.** A physically inclined finite planar face and a
  finite tilted cylinder yield shared 3D intersection arcs and two
  reconstructible UV parameter curves. Validated seam crossing does not
  introduce a physical B-rep vertex; cut endpoints belong to actual face
  boundaries. Preview-only Boolean region probes are classified; near-grazing
  trim events remain Held and no general topology delta can commit.
  See [v0.36 oblique UV intersection contract](docs/oblique-plane-cylinder-pcurves.md).


- Build a **connected three-cylinder crossing chain** in one analytic B-rep: the through and deeper blind cutters intersect, and the deeper and shallower blind cutters intersect, while the outside pair stays separated. Source-backed angular events define shared physical line/circle edges, cylinder UV trims, two upward blind floors and one closed genus-one shell. Analytic source mass/centroid are independently checked against union-of-slices quadrature; display triangles reuse shared B-rep edges. Exposed in the `compose_connected_triple` JSON API. See [v0.35 connected triple B-rep](docs/connected-triple-brep.md).

- Plan **connected through/two-blind circular boundary transitions across different floor depths**, partitioning cylinder sectors into continuing, terminating and emerging source-owned arcs. Validate floor-area growth and source replay, without premature multi-cut topology commits. See [v0.34 transition contracts](docs/connected-three-cutter-transitions.md).

- Build **one analytical B-rep from multiple strictly separated, individually interacting through/blind crossing clusters** with different blind depths, shared upper/lower cap loops, full coedge/UV remapping, source-atomic Undo/Redo and a certified global watertight mesh. The new `compose_crossing_clusters` API also permits strictly disjoint stand-alone blind cuts. See [v0.33 composed crossing shells](docs/composed-crossing-shells.md).

- Plan multiple **interacting through and blind cylindrical cutters** (up to 12) at different finite blind depths using exact-in-form exposed circular arcs, closed source-indexed contour loops, connected-component and planar-hole classification, and per-depth-slab volume integration. Tangencies, unresolved triple junctions and numerical ambiguity produce a typed Held state, with no topology commit. See [v0.32 multi-cut layered regions](docs/multi-cut-layered-regions.md).

- Construct **one true crossing through/blind analytical solid** in a certified axis-aligned block: shared 3D edges, two-use coedges, UV trim loops, one closed genus-one shell, independent analytic mass, source-atomic Undo/Redo and certified watertight display triangulation. Available through the Rust API and `compose_mixed_crossing` geometry.v1 command. See [v0.31 mixed crossing B-rep](docs/mixed-crossing-brep.md).

- Read-only geometric planning for one **overlapping** upright through-hole and a finite-depth top-entry blind-hole: exact-in-form circular lens area, independent volume inclusion/exclusion, source-validated crossing points and exposed arc witnesses, typed near-tangency holds, and no topology commit. See [v0.30 overlap contract](docs/mixed-footprint-overlap.md).

- Extrude a simple, possibly concave, planar XZ profile along +Y, or create a block.
- Maintain explicit oriented B-rep: shells, material-facing face sense independent of surface UV frames, faces, loops, coedges (twin/next/prev), edges and vertices.
- Edit the polyhedral B-rep atomically: split edges and faces, maintain twin/loop links, preserve unaffected local handles, validate and rollback on failure.
- Chain Euler operations using typed outputs created earlier in the same transaction, with one Undo/Redo step for the batch.
- Compute exact planar B-rep Boolean union/intersection/difference for validated axis-aligned boxes and bounded analytical cylinder cuts, including one closed B-rep with both disjoint full-height through holes and independent single/crossed blind regions.
- Classify source volumes for validated polyhedral solids, cylinders, through-holes and finite-depth blind bores; construct source-validated plane/cylinder opening traces and a two-circle transverse crossing arrangement. General curved Boolean plans beyond these certified cases remain preview-only.
- Atomically edit radius and center of an analytic through-hole, and radius/center/depth of a blind bore, with Undo/Redo, independent B-rep reconstruction and stable semantic handles through BrepArena.
- Edit the canonical shared B-rep's 3D curve carriers, per-coedge 2D UV trims, and support-surface charts through source-preserving transactions, copy-on-write catalogs, owner-bound generational handles and Undo/Redo.
- Create, remove and resize up to 32 independent full-height analytic cylindrical holes in a box, including true inner UV loops/coedges and genus/analytic mass changes; combine structural work and carrier reparameterization in one atomic transaction.
- Derive watertight, oriented display meshes for analytical multi-hole bodies, including top and bottom caps with **different** UV-loop sets for mixed through/blind cuts, nonconvex blind floors, signed-volume certificates and source face ownership.
- Validate closed manifold incidence, vertex fans and planar non-self-intersecting loops; the restricted polyhedral editor has Euler characteristic 2, while the multi-hole analytic B-rep has Euler characteristic 2−2N.
- Configure model-space absolute/relative/angular tolerances for each solid or JSON request.
- Triangulate polygon caps and sides, retaining per-triangle source face IDs for CAD selection.
- Compute prism volume, surface area, centroid and axis-aligned bounds; translate bodies.
- Ray-pick a triangle and map the hit back to the owning topological face.
- Accept `geometry.v1` JSON commands through the `gefest-geometry` stdio executable.

This implementation is **original code** based on common computational geometry methods. It contains no Parasolid binaries, Siemens headers or decompiled implementation. The historical `PK_*` API inventory, if used for planning, is not source code or an implementation specification.

## Build and test

For v0.35 connected triple analytical source and watertight display, run `cargo test --test connected_triple_brep && cargo test --bin gefest-geometry connected_triple_cli`.

For v0.33 composed crossing-shell topology and global mesh, run `cargo test --test composed_crossing_brep`.
For the independent v0.32 multi-cut region graph planner, run `cargo test --test multi_cut_regions`.
For the executable v0.31 one-pair crossing topology and display, run `cargo test --test mixed_crossing_brep && cargo test --test mixed_crossing_mesh`.
For the earlier v0.30 read-only overlap planner, run `cargo test --test mixed_footprint`.

Requires a stable Rust toolchain; the minimum supported Rust version has not been revalidated for the v0.24 triangulation dependency:

```sh
cargo test --all-targets
cargo build --release
```

Example JSON request:

```sh
printf '%s' '{"operation":"extrude_profile","profile":[[0,0],[2,0],[2,3],[0,3]],"height":4}' | target/release/gefest-geometry
```

A batch of edits can also be sent as one atomic request:

```sh
printf '%s' '{"operation":"edit_solid","source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},"edits":[{"kind":"split_edge","edge_id":0,"fraction":0.5},{"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":2}]}' | target/release/gefest-geometry
```

The `edit_report` includes one revision increment and a list of typed provenance events, with newly allocated IDs. If any operation fails, no partial solid is returned. The old construction operations remain compatible.

The response carries `schema_version: "geometry.v1"`, `ok`, and either a `solid` or an `error`. `solid` includes `vertices`, indexed triangle `faces`, `triangle_face_ids`, counts of vertices/edges/coedges/loops/shells/faces and Euler characteristic under `topology`, analytic `mass` and `bbox`. The request accepts an optional `tolerance` object; omitted fields use the defaults. The existing triangle field named `faces` is retained for API compatibility. `geometry.v1` is a provisional bridge contract, not a STEP or Parasolid serialization format.

Example library use:

```rust
use rustsolid::{Point2, Solid};

fn main() -> Result<(), rustsolid::GeometryError> {
let profile = [
    Point2 { x: 0.0, z: 0.0 },
    Point2 { x: 10.0, z: 0.0 },
    Point2 { x: 10.0, z: 5.0 },
    Point2 { x: 0.0, z: 5.0 },
];
let solid = Solid::extrude_xz(&profile, 8.0)?;
assert_eq!(solid.euler_characteristic(), 2);
assert_eq!(solid.mass.volume, 400.0);
    Ok(())
}
```

## Library transactions

```rust
use rustsolid::{Solid, Point3, EdgeId, FaceId};

fn main() -> Result<(), rustsolid::GeometryError> {
let mut body = Solid::block(Point3 { x:0.0, y:0.0, z:0.0 }, 2.0, 3.0, 4.0)?;
let corners = body.faces[0].boundary.clone();
let report = body.edit_atomic(|tx| {
    tx.split_edge(EdgeId(0), 0.5)?; // fraction along oriented canonical edge
    tx.split_face(FaceId(0), corners[0], corners[2])?;
    Ok(())
})?;
assert_eq!(report.revision_after, 1);
    Ok(())
}
```

The edit closure runs against a cloned working model, validates after every edit,
and swaps it into the body only on a successful commit. Dropping a transaction
or returning an error leaves the original body untouched. Existing `VertexId`,
`FaceId`, `EdgeId`, `CoedgeId` and `LoopId` survive supported local edits;
new elements get appended IDs. Regenerating a body from an edited sketch **does
not yet preserve IDs**, and these are not persistent topological names.

The [Composite pattern](https://en.wikipedia.org/wiki/Composite_pattern) is
appropriate for *ownership* (`Solid → Shell → Face → Loop → Coedge`).
`Coedge.twin`, shared `Edge` and `Vertex` references instead form a graph;
they must not be recursively owned by multiple parents. The API exposes
`topology_children()` for ownership traversal only.

## Gefest integration

The existing Gefest backend uses CadQuery/OpenCascade, and its geometric constraint solving is handled by a separate Rust workspace (`solver-core`). RustSolid is deliberately independent of both. When using the experimental Gefest Python adapter, set `GEFEST_GEOMETRY_BINARY` to the compiled `target/release/gefest-geometry`; opt in with `GEFEST_GEOMETRY_BACKEND=rust`. The OpenCascade path should remain the production default until parity, reproducibility and regression gates pass.

## Known limits

The polyhedral `Solid` editor still represents **one connected genus-zero shell with one outer loop per face**. The shared analytic B-rep supports multiple disjoint full-height cylindrical holes and one closed genus-N shell in the bounded structural editor. The original Boolean endpoint handles exact axis-aligned box CSG, a single full-height through-cylinder or a certified top-entry blind cylinder; additional through-holes use the explicit `structural_transaction` path. Complete global Boolean region classification, arbitrary cross-depth or oblique intersecting cuts, connected regions of three or more crossing cylinders, generalized through/blind **intersections beyond one certified transverse pair and all tangencies**, arbitrary structural cuts, multi-component bodies, NURBS, general cavities, source-independent meshing and stable topological naming across connectivity changes, STEP/X_T, fillets and production-grade numerical guarantees remain future work. Not suitable for safety-critical or precision manufacturing use.

The validated Euler command history supports atomic MEV/MEF/KEV/KEF batches with Undo/Redo. Box Boolean evaluation creates an independent B-rep rather than mutating one operand's history. General Boolean splitting of curved topology, vertex move, face extrusion and fillets remain unsupported. Arena transaction candidates still clone their command history (at least O(n) memory work).

See [docs/architecture.md](docs/architecture.md) for invariant details, tolerance semantics, limits, and the next gates.

## License

Mozilla Public License 2.0 (`MPL-2.0`). See [LICENSE](LICENSE).

## v0.4: topology identity and inverse journals

For deterministic topology names call `solid.with_feature_key("part/feature")` before editing. `topology_name`, `resolve_topology_name`, `topology_handle` and `resolve_topology_handle` separate history-relative names from process-local generational handles. Edited descendants have parent provenance and are included in `EditReport.named_changes`. Transactions now mutate through an exclusive borrow with operation-local inverse snapshots and rollback on any failure. `EditReport.journal` exposes snapshot counts. General geometric face matching on arbitrary regeneration and persistent serialized handles remain out of scope.


## Command + Memento history (v0.5)

The typed `EulerCommand` records ordered, replayable CAD operations using
feature-relative `TopologyName` selectors. `CommandHistory` executes a batch
atomically and retains an opaque, bounded inverse-journal Memento for Undo.
Redo replays the command and regenerates fresh handles, preserving the ABA
protection established in v0.4. Undo/Redo revisions increase monotonically.

```rust
use rustsolid::{
    CommandBatch, CommandHistory, EdgeId, EulerCommand, Point3, Solid,
    TopologyEntity,
};
let body = Solid::block(Point3 { x: 0., y: 0., z: 0. }, 2., 3., 4.)?
    .with_feature_key("part-1/extrude-1")?;
let mut history = CommandHistory::new(body)?;
let edge = history.solid().topology_name(TopologyEntity::Edge(EdgeId(0)))
    .unwrap().clone();
history.execute(CommandBatch::single(EulerCommand::SplitEdge {
    edge, fraction: 0.5,
}))?;
history.undo()?;
history.redo()?;
# Ok::<(), rustsolid::GeometryError>(())
```

For a one-request replay/Undo/Redo demonstration, send a JSON
`{"operation":"command_history", "source":{"kind":"block", ...},
 "feature_key":"part-1/extrude-1", "batches":[{"commands":[...]}],
 "undo":1, "redo":1}` request to `gefest-geometry`.

The present split-edge / split-face operators preserve Euler characteristic;
undo removes *their own* generated topology using incremental inverse frames.
Independent make/kill Euler operations and arbitrary topology reconciliation
are future work.

## v0.6 scoped independent inverse Euler operators

`Solid::kill_edge_vertex(vertex)` collapses a straight degree-two vertex, and `Solid::kill_edge_face(edge, removed_face)` merges adjacent coplanar faces. These inverse Euler edits can target compatible earlier splits even after unrelated modeling edits. The kill path currently rebuilds and compacts the closed polygonal B-rep; named surviving entities are rebound, and any changed numeric ID invalidates its former runtime generational handle. This path has a full-body copy and is deliberately separate from the incremental make-only transaction path. These are scoped inverse operators, not unrestricted solid topology algebra.

`predicates` provides adaptive exact-sign `orient2d`/`orient3d` (using the MIT/Apache-licensed `robust` crate), and tolerance-aware plane/plane and segment/plane classifications. Geometric intersection coordinates are approximate floating point; arbitrary trimmed-surface Boolean and NURBS editing are not implemented. v0.19 provides a separate exact-planar Boolean path for two certified axis-aligned boxes.

### Geometry Foundation v0.7

The new `geometry` module separates analytic carriers from B-rep topology:
`Curve2` (2D p-curves), `Curve3` (lines, circles), `Surface3` (planes, cylinders,
spheres), explicit parameter domains, UV frames and evaluated derivatives.
`GeometryStore::from_solid` connects the existing polyhedral B-rep to planar
support surfaces, line edge carriers and oriented UV trim lines. Each view is
bound to one body incarnation and revision and must be rebuilt after an edit.

```rust
use rustsolid::{GeometryStore, ModelUnits, Point3, Solid, FaceId, Surface3};

let body = Solid::block(Point3 { x:0.0, y:0.0, z:0.0 },2.0,3.0,4.0)?;
let geometry = GeometryStore::from_solid(&body, ModelUnits::default())?;
geometry.validate_bindings(&body)?;
assert!(matches!(geometry.face_surface(FaceId(0)), Some(Surface3::Plane(_))));
# Ok::<(), rustsolid::GeometryError>(())
```

Per-edge and per-vertex geometry-view tolerance overrides are supported.
`PredicateKernel` distinguishes adaptive exact-sign orientation from
engineering-scale "within tolerance" classification, and shared robust 2D
predicates now back the core polygon operations. This remains an independent
polyhedral modeling prototype: surfaces are evaluable mathematically, but
curved surfaces are not yet bound as native curved topology or exportable STEP.


## Native analytic B-rep cylinder (preview)

```rust
use rustsolid::{CylindricalBrep, GeometryTolerance, Point3, FaceId, Surface3};
let cylinder = CylindricalBrep::upright(
    Point3 { x:0.0, y:0.0, z:0.0 }, 5.0, 20.0, GeometryTolerance::default())?;
assert!(matches!(cylinder.face_surface(FaceId(2)), Some(Surface3::Cylinder(_))));
assert_eq!(cylinder.coedges[4].face, cylinder.coedges[5].face); // UV seam
let render_mesh = cylinder.tessellate(96)?;
assert_eq!(render_mesh.triangle_faces.len(), render_mesh.triangles.len());
# Ok::<(), rustsolid::GeometryError>(())
```

The topology contains exact circular edges and an analytic cylindrical side
with distinct periodic UV seam trims. The displayed triangles are a
resolution-dependent approximation; the underlying curves, surfaces and
analytic mass do not depend on facet count. Supports the new `analytic_cylinder`
JSON CLI request (default 64 facets). Native curved Boolean and Euler operations
are not yet supported.


## v0.9 — common B-rep interface

`BrepModel` is a single, geometry-backed *read-only* topology graph for both
`Solid` and `CylindricalBrep`. Build it using the extensible `BrepSource` trait,
`BrepModel::from_polyhedral`, `BrepModel::from_cylinder`, or the convenience
`BrepBody::shared()` method. The same methods navigate faces, loops, coedges
and edge uses, and query their `Surface3`, `Curve3` and `Curve2` geometry.

Closed circular edges and same-face periodic UV seams remain analytic; they are
not converted into faceted B-rep edges. `BrepBody::tessellate(segments)` produces
a separate display mesh. Call `BrepModel::validate()` before using untrusted
or mutated geometry graph data.

A new `inspect_brep` JSON request can query either a polyhedral block/profile
or an analytic cylinder. It returns the original `geometry.v1` mesh envelope
with an additional `solid.brep` summary. Editing remains source-specific. This
stage does **not** provide Boolean modeling, arbitrary radial edge incidence,
interior holes, or general curved face editing.

See [docs/unified-brep-parasolid-comparison.md](docs/unified-brep-parasolid-comparison.md).

## v0.10: boundary accuracy probes and Geometry Query Language

RustSolid now exposes a small, typed, read-only [GraphQL-inspired B-rep query language](docs/geometry-query-and-probes.md) and a seeded [boundary-focused model generator](docs/geometry-query-and-probes.md) for analytic and metamorphic precision checks. Queries project nested fields and filter geometric carriers with strict depth/field/output budgets. The probe suite reports independent analytic errors separately from expected tessellation losses.

CLI operations: `geometry_query` (requires `source` and a selection `query`) and `probe_kernel` (optional `seed`). Both return a machine-readable `data` object under `geometry.v1`. These are the first testing and inspection tools, not a general GraphQL implementation or a proof of industrial CAD robustness.

## v0.11: first-party regression corpus

No external CAD engine is treated as a correctness oracle. See [self-contained regression fixtures](docs/regression-corpus.md): pinned independent analytical golden data, metamorphic relations, and reproducible counterexample shrinking, plus a bounded seed-based discovery suite. CI checks the permanent 20-fixture corpus. The `regression_corpus` JSON request produces an outcome report without importing OpenCascade or Parasolid.

## v0.12: mandatory operator quality contracts

The [modeling operator contract gate](docs/operator-contracts.md) covers each registered primitive constructor and every `EulerCommand` variant with an independent analytic oracle, a metamorphic property and at least two invalid-boundary scenarios. The executable, pinned inventory lives in [tests/fixtures/operator_contracts.v1.json](tests/fixtures/operator_contracts.v1.json), while [tests/operator_contracts.rs](tests/operator_contracts.rs) exercises real modeling operations and intentionally malformed manifests. Adding an Euler command without updating its exhaustive operation-kind match cannot compile; omitting its quality contract fails CI. No external CAD backend serves as a truth oracle.

## v0.13: analytical support-surface intersections

[The intersection engine](docs/intersections.md) constructs untrimmed analytic curves for plane/plane and plane/cylinder pairs, without deriving geometry from tessellation. Results distinguish transverse crossings, tangency, disjoint/coincident supports, and **numerically ambiguous** near-tangent/near-parallel cases. Sample residual certificates compare every constructed branch with both source surface equations. The pinned [intersection fixtures](tests/fixtures/intersections.v1.json) have their own compulsory analytic/metamorphic/boundary test gate in CI. Neither OpenCascade nor Parasolid serves as a truth oracle. This is not a trimmed face Boolean operation.

## v0.14 — bounded B-rep face intersections

The [trimmed intersection engine](docs/trimmed-intersections.md) clips analytic plane/plane and plane/cylinder support curves to **actual UV contours of two selected B-rep faces**. It supports simple planar polygons, circular disks and full periodic cylindrical sides. Result branches carry analytic curves and finite parameter ranges; isolated point contacts, boundary overlaps and near-degenerate ambiguity remain distinct. It has [12 pinned, independent fixtures](tests/fixtures/trimmed_intersections.v1.json) in CI. This is not yet a Boolean operation or complete arbitrary NURBS/UV-trim support.

## Evidence-based Parasolid capability comparison

The versioned [official Siemens ↔ previously observed PK API ↔ RustSolid capability map](docs/parasolid-api-crosswalk.md) records release-specific, publicly documented behavior (including v38.0 versus v39.0 differences), current kernel gaps and independently falsifiable acceptance tests. Earlier decompiler-derived symbol names are clearly labeled **not reverified in this audit**; private Siemens code and inferred internal algorithms are never incorporated. A first-party Python CI test ensures every ledger item remains sourced, versioned and paired with a non-vendor correctness oracle.

## v0.15: typed IntersectionGraph candidate topology

[IntersectionGraph](docs/intersection-graph.md) consumes the verified v0.14 UV-clipped intersection branches and builds an **immutable** graph with typed operand-namespaced face provenance, local face tolerance records, physical boundary edge/coedge witnesses, analytic intervals, contact classification, and read-only validation against both input B-reps. Full periodic circles are graph self-loops with a synthetic non-cut seam anchor; partial arcs crossing u=0/2π are coalesced into one unwrapped physical branch. The pinned [12-scenario graph regression corpus](tests/fixtures/intersection_graph.v1.json) is checked by CI along with invariance/tampering tests. The `geometry.v1` API gains read-only `intersection_graph` queries; no Boolean cut or B-rep mutation occurs.

## v0.16 — immutable TopologyDeltaPlan

The [IntersectionGraph](docs/intersection-graph.md) can now be converted into a preview-only [TopologyDeltaPlan](docs/topology-delta-plan.md). Each plan records source vertex reuse, proposed edge splitting, interior face vertices, virtual periodic anchors, original source curves, per-face cut intentions, and explicit prerequisites. It does not edit any B-rep or decide which regions a Boolean would keep. Every plan reports can_apply=false. Twelve pinned mathematical scenarios and tampering tests are enforced by CI.


## v0.29: one exact analytical body with through tunnels and blind floors

The new [heterogeneous mixed-region builder](docs/mixed-through-blind-regions.md)
combines one or more actual full-height analytical through-bores with
separately verified top-entry blind components in the **same connected
closed B-rep shell**. A blind component may be a circular disk or a
union of exactly two partially crossing, equal-depth disks.
Cross-family source footprints are strictly disjoint: contact,
tangency or overlap fails closed.

Both source families share the original top box face. The bottom face
contains only the through-hole inner loops; blind components end in
separate upward-facing analytical floor faces. The shell's **genus is
exactly the number of through holes**, independent of the number of
blind floors and crossing arcs. The complete source-generated candidate
rechecks all 3D/UV carrier references, face senses, Euler characteristic,
analytical surface area and exact removed-volume moments.

`MixedThroughBlindBrep` supports isolated, strong-atomic
parameter changes for existing through and blind sources, including
blind-pair connectivity repartition, with monotonically revisioned
Undo/Redo. `BrepBody::MixedRegions`, common volume queries,
`geometry_query`, `inspect_brep` and request-local
`geometry.v1 compose_mixed_regions` expose the validated shape.
Its optional `display_segments` derives a watertight face-owned
mesh with separate cap triangulation and analytic-vs-faceted
mass convergence witnesses. The underlying exact B-rep is not
replaced by triangles. Structural source-list insertion/deletion,
intersecting through/blind side trims and arbitrary curved Boolean
operators remain future work.

## v0.28: mixed connected blind-cut region components

The new [compound analytical blind-region compositor](docs/compound-blind-regions.md)
builds a **single exact genus-zero B-rep** from up to 16 strictly certified
top-entry blind cylindrical cuts, partitioned into independent connected
footprint regions. Each region is either a full isolated circle or a
**two-circle partially overlapping pair sharing one floor depth**.
Independent regions can have different depths.

The compositor reuses the source-certified v0.25 circular bores and v0.27
true split circular-arc crossings, then remaps their actual analytic
vertices, edges, coedges, UV loops, surfaces and independent floor faces
into one manifold shell. All apertures share the box's original top
face, with one inner loop per connected region. Volume, area and
centroid derive from exact source area/perimeter/moments, not facets.

`CompoundBlindBrep` supports **atomic changes to the connected region
classification**: a crossing pair can split into separate single holes
and be restored through Undo/Redo with one revision per entire batch.
The shared `BrepBody` query and read-only arena inspect this topology;
generic arena mutation with persistent feature names after renumbering
is deliberately not claimed.

`geometry.v1 compose_blind_regions` returns component-to-source indexes,
real floor/side topology IDs, analytical mass, revision/history and
optional `display_segments` with a **watertight multi-floor display mesh**
and an independent polygonal mass certificate. The ordinary
`analytic_compound_blind` B-rep source kind also supports inspection
and geometry queries. Tests cover two pairs, mixed pair/single,
atomic repartition, fail-closed three-way crossings and facet
convergence.

Global `TopologyDeltaPlan` is still preview-only, and three-way
circle arrangements or overlapping cuts with stepped floors are not
silently approximated.

## v0.27: executable analytic crossing-circle region reassembly

The new [intersecting-blind analytic Boolean](docs/intersecting-blind.md)
supports subtracting **two partially overlapping, nonnested, equal-depth**
top-entry cylinders from a validated block. Exact circle-circle
intersection witnesses split both cylinder boundaries into visible
`Curve3::Circle` arcs. Two real shared vertical crossing edges join
the cylindrical patches. The existing box top face gets one compound
two-arc inner UV loop, and a new upward-facing floor has a matching
two-arc outer loop.

This is a **genuine new genus-zero analytical B-rep**, with 12
vertices, 18 edges, 36 coedges, 10 loops and 9 faces. The union-of-two-disks
area, exposed perimeter and analytic removed moments determine exact
volume, area and centroid, independently of display segmentation.
`IntersectingBlindBrep` supports fully atomic two-source parameter edits
with one-step Undo/Redo. True tangencies, containment, differing floor
depths, tilted axes and other unsupported cuts fail closed.

A new display mesher independently certifies a watertight top opening,
two inward-facing cylinder strips and one nonconvex united floor. The
specialized `geometry.v1 compose_intersecting_blind` operation provides
source arc/crossing evidence, exact mass, revisions and optional
`display_segments` with face-owned triangles. It does not make the
general `TopologyDeltaPlan` executable.

## v0.26: composed independent analytic blind-cut regions

[CompositeBlindBrep](docs/composite-blind.md) builds one exact, validated
genus-zero analytical B-rep from a certified axis-aligned block and up to
16 pairwise-disjoint top-entry blind cylindrical cuts. Each cut has its
own certified closed IntersectionGraph source opening, planar disk floor,
true circular 3D edges, cylindrical side and UV trims. All openings share
the original top planar face with one inner loop per bore.

The composition rechecks full manifold topology and recomputes exact
volume, area and centroid using disjoint removed cylinder moments. It is
available via Rust `CompositeBlindBrep::build` and the new
`geometry.v1 compose_blind` operation. It intentionally returns analytic
topology and mass **without inventing a display mesh**, and it does not
claim to implement general intersecting curved Boolean arrangements or
mutable composite history. Previously available single-hole editing and
mesh paths remain unaffected.

## v0.25: source-certified top-entry analytic blind cylindrical Boolean

[The new blind-bore builder](docs/analytic-blind-hole.md) subtracts an
upright cylinder whose top cap exactly meets the top face of an axis-aligned
block while its lower cap stays strictly inside the material. The actual
top face gets an inner circular UV loop, a new independent planar **floor**
closes the blind cavity, and the cylindrical wall uses the original
analytic radius, 3D circular edges, and periodic seam carriers.

An independent `IntersectionGraph` check certifies that the two source
surfaces produce exactly one full closed top-opening trace. The resulting
B-rep has 8 faces, 9 loops, one inner loop, and genus 0. Volume, surface
area, and centroid use analytical removed-cylinder moments; the display
mesh is separately certified to be watertight, correctly oriented and
convergent with increasing circumference resolution.

A new typed `FaceSense` distinguishes the physical material-facing normal
from the surface's parameterization normal. This represents the inward
wall and upward-facing circular floor without synthesizing false UV
geometry. The `geometry_query` face `sense` field and `inspect_brep` face
sense metadata expose this orientation to Gefest CAD.

`BooleanOutcome::BlindHole(BlindHoleBrep)`,
`BrepArena::from_blind_hole`/`update_blind_hole`, and source kind
`analytic_blind_hole` expose the result and its strongly atomic
radius/center/depth editing. Existing full through-hole and box CSG
paths retain their contracts. Independent source-controlled regressions
and mandatory CI Boolean/CLI gates cover the opening, floor, mass,
orientation, rejected cutter alignment and Undo/Redo.

**Not a universal analytic Boolean kernel:** oblique cylinders,
intersecting cuts, arbitrary trimmed surfaces, multiple cavities and
general `TopologyDeltaPlan` execution remain unsupported.

## v0.24: certified multi-loop analytic display meshing

The new [FaceRegionMesher](docs/face-region-meshing.md) samples the *current*
analytic circular B-rep edges to triangulate both planar box caps with zero
or more circular inner UV loops. It builds inward-facing cylindrical wall
facets and reuses the original box wall triangles, sharing mesh vertex IDs
at all contours. Every displayed result is checked for correct per-face
ownership, outward winding, opposite directed edge uses, full manifold
closure, face-region area, and independent signed-volume/area formulas.

The resulting `CertifiedDisplayMesh` reports faceted volume and area,
exact analytical mass and volume faceting error separately. Tessellation
resolution never changes the source B-rep, its genus or editing revision.
The existing `structural_transaction` JSON API accepts an **optional**
`display_segments:12..=512` and then returns the display `solid`
payload and detailed quality certificate. Without this option, its previous
no-mesh response remains compatible. First-party
[`multi_loop_mesh` tests](tests/multi_loop_mesh.rs) and mandatory CLI CI
witnesses cover two/three holes, convergence, remote coordinates and
component orientation.

Generalized non-box Boolean operations, editable intersecting curved
trim arrangements, and arbitrary NURBS remain unsupported.

## v0.23: unified genus-changing structural and carrier transactions

[The structural B-rep editor](docs/structural-brep.md) supports a verified
axis-aligned box with up to 32 pairwise-disjoint, fully penetrating +Y
cylindrical holes. Each insertion creates actual analytic circular edges,
paired coedges, inner UV loops and a periodic cylindrical face. Deletion
retires those entities, reindexes surviving dense topology and recomputes
volume, surface area and centroid from exact source geometry.

`BrepStructuralTransaction` combines insertion, removal, radius/center
changes and the v0.22 carrier reparameterizations in **one atomic edit**.
Generational feature handles survive dense compaction, while retired and
abandoned handles never revive. Undo/Redo treats the entire batch as one
step. The editor can be initialized from an existing single-hole feature
or detached from a validated `BrepArena` without modifying that arena.

The `geometry.v1` API adds `structural_transaction`, returning a validated
analytic B-rep topology summary, source feature bindings, exact mass and
commits with created/retired/remapped identities. Multi-hole tessellation
is not implemented yet; the JSON response does not fabricate a display
mesh. [First-party structural regressions](tests/structural_editor.rs) and
mandatory CLI/CI gates verify genus, mathematical mass, atomicity, ABA,
and source history preservation.

Only this bounded analytical feature family is executable; general
`TopologyDeltaPlan` application remains preview-only.

## v0.22: canonical transactional carrier and UV-trim editor

The new [BrepCarrierEditor](docs/carrier-editor.md) owns a mutable validated
`BrepModel`, including 3D curve, 2D coedge trim, and supporting surface
catalogs. Copy-on-write replacements and dense-catalog compaction preserve
typed generational handles bound to their owning edges, coedges and faces,
rather than to the changing carrier-array indices.

A `BrepCarrierTransaction` can replace compatible analytic carriers,
reparameterize line/circle edges and UV trim lines, shift circular trim
parameters by complete turns, and rebase planar or cylindrical face UV
coordinates together with all owned coedge trims. Only shape-preserving
edits pass commit: original physical paths and support surfaces, analytic mass,
topology and model bounds are checked before candidate publication. Ignored
staging errors poison the entire batch; Undo/Redo retain a single modeling step
with monotonic revisions and stale-snapshot detection.

A separate `BrepArena::detached_carrier_editor()` can start from an existing
arena snapshot. It does not silently mutate that arena's Euler or curved
parametric feature history. The `geometry.v1` CLI exposes
`carrier_transaction` for safe request-local reparameterizations; direct
typed replacement is provided by the Rust API. Independent regressions
in [tests/carrier_editor.rs](tests/carrier_editor.rs) and mandatory GitHub
Actions gates cover genus-one trims, periodic seams, copy-on-write aliasing,
transaction failures and history.

This slice **does not** create or delete topological edges/loops or enable
generalized curved Boolean application of `TopologyDeltaPlan`.

## v0.21: executable genus-one analytic through-hole

[The first real curved Boolean](docs/analytic-through-hole.md) supports subtracting a
vertical, fully penetrating analytic cylinder from a certified axis-aligned
block. The resulting `ThroughHoleBrep` contains 7 analytic faces, 9 UV loops
(including 2 inner circular loops), 2 circular edges, a periodic cylindrical
seam and one genus-one shell with Euler characteristic 0. Exact mass and
boundary geometry are independent from the separate display tessellation.

`BooleanKernel::apply` and `BrepArena::boolean` return
`BooleanOutcome::Analytic` for the supported subtraction. A
`ThroughHoleTransaction` stages radius/center modifications on a private
candidate, validates before publication, and supports Undo/Redo.
`BrepArena::update_hole` additionally reconciles the generational topology
registry before publishing the changed geometry.

The `geometry.v1` `boolean` CLI accepts box minus a full analytic cylinder,
returning genus, trims and source-face lineage metadata. A direct
`analytic_through_hole` source can also be inspected and queried.
General curved Boolean application of `TopologyDeltaPlan` is still refused.

## v0.20: analytic cut region classification

[BooleanRegionClassifier](docs/boolean-regions.md) determines whether a
spatial point lies strictly inside, outside or on a validated polyhedral B-rep
or analytic right circular cylinder. Polyhedral membership uses oriented
solid-angle winding of the actual triangle boundary; curved cylinders use
their exact local-frame radial and axial equations.

For a validated two-face `IntersectionGraph`, the classifier samples both
sides of each analytic crossing arc on **each** source face, evaluates
the opposite body's volume membership, and records the retained and discarded
regions for Union, Intersection and ordered Difference. It handles periodic
circles and open plane-cylinder traces, preserves original face provenance,
and refuses unsupported tangencies or ambiguous near-boundary probes.
The witness plan can be revalidated against both current source bodies.
This is *read-only local region evidence*: face-region reconstruction,
inner UV loops, multiple shells, and general curved Boolean commits
remain future gates. `TopologyDeltaPlan` cannot yet be applied.

The `geometry.v1` CLI exposes `boolean_region_plan`, with `first` and
`second` face operands, a `mode`, and explicit `can_apply:false`.
First-party tests and the dedicated CI gate verify these classifications.

## v0.19: exact planar box Boolean B-rep

[The first Boolean executor](docs/boolean-box.md) performs `union`, `intersection`
and ordered `difference` on two independently certified axis-aligned
parallelepipeds. It partitions the exact source face coordinates into Cartesian
3D regions, classifies occupied cells, builds outward planar B-rep patches and
validates every result through `Solid::validate`, including manifold incidence,
Euler characteristic, mesh ownership, volume and area. This uses exact planar
source positions, **not approximate faceting of curved input**.

`BooleanResult` explicitly distinguishes an `Empty` result from a validated
`Solid`, and reports per-result-face lineage back to left/right source
faces (including reversed cut faces for subtraction). Non-manifold results,
multi-body outputs, enclosed cavities, near-coincident distinct planes,
general extrusions and analytic cylinders are rejected with explicit errors.
`BrepArena::boolean` leaves both source arenas untouched and creates a new
independent result, rather than falsely preserving historical handles.

The `geometry.v1` CLI exposes `{"operation":"boolean","mode":"union|intersection|difference",
"left":{...},"right":{...}}` with standard solid output and face lineage
metadata. [Independent Boolean regressions](tests/boolean_box.rs) and CLI smoke
tests are pinned by CI. See the [Boolean scope and next gates](docs/boolean-box.md).

## v0.18: dependency-aware B-rep transactions

The [BrepArena](docs/brep-arena.md) now supports typed references to geometry
created earlier in the *same* batch. Use `split_edge_ref` /
`split_face_ref` to obtain a transaction-scoped `ArenaStep` and pass
`step.output(ArenaOutput::InsertedVertex)`,
`step.output(ArenaOutput::RemainderEdge)`,
`step.output(ArenaOutput::CreatedFace)` or
`step.output(ArenaOutput::DiagonalEdge)` to later `*_ref` methods.
The original handle-only methods remain supported.

A CLI example uses the first split's newly created remainder edge as the
second split's input. The output reference is scoped to the current batch
and has a zero-based `step` index:

```sh
printf '%s' '{"operation":"arena_transaction","source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},"batches":[[{"kind":"split_edge","edge_id":0,"fraction":0.25},{"kind":"split_edge","edge_ref":{"step":0,"output":"remainder_edge"},"fraction":0.5}]],"undo":1,"redo":1}' | target/release/gefest-geometry
```

Forward references, foreign transaction outputs, invalid output roles, retired
elements and failed geometry validations reject the whole batch. The kernel
compiles dependencies from actual provisional Euler outputs into replayable
semantic names; it does not guess future numeric B-rep IDs. Each successful
batch remains a single revision and Memento, including mixed MEV/MEF/KEV/KEF
operations. Nine new pinned first-party regression cases bring the arena suite
to 19 scenarios. Dependent preflight currently clones and validates a private
working body; it is not yet the optimized sparse mutable geometry store.

The analytical cylinder remains read-only. This release does not implement
general Boolean application of `TopologyDeltaPlan` or a Parasolid-compatible
B-rep editor.

## v0.17: sparse generational B-rep arena and atomic edits

The first [BrepArena](docs/brep-arena.md) manages topology through named sparse generational handles, providing read-only geometry snapshots for polyhedral solids and analytical cylinders. It executes type-checked MEV/MEF/KEV/KEF commands as one strongly atomic `CommandHistory` batch with a shared Undo/Redo Memento and named created/retired/remapped entity reports. Candidate edits and the identity registry are reconciled before publishing, so a late failed command changes nothing. Every retired slot increases its generation before reuse. The analytic cylinder is intentionally read-only, and no preview-only `TopologyDeltaPlan` can bypass that restriction. The new `geometry.v1` `arena_transaction` CLI workflow is stateless and returns transient handles valid only within the current process. [Pinned first-party regressions](tests/fixtures/arena_transactions.v1.json) and CI check generations, rollback, reindexing, Undo/Redo, and analytic mass invariance. Candidate cloning is correctness-first; persistent O(changed-elements) editing is future work.
