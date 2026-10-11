# RustSolid

RustSolid is an **independent Rust solid geometry kernel** being developed for [Gefest CAD](https://github.com/boltholds/Gefest-CAD).

> **Status: v0.18 experimental kernel: dependency-aware Euler batches, sparse generational identity arena, atomic polyhedral editing and read-only analytic B-rep snapshots.** It is not currently a replacement for Parasolid or OpenCascade. Its API does not claim compatibility with proprietary kernels.

## Current capabilities

- Extrude a simple, possibly concave, planar XZ profile along +Y, or create a block.
- Maintain explicit oriented B-rep: shells, faces, loops, coedges (twin/next/prev), edges and vertices.
- Edit the polyhedral B-rep atomically: split edges and faces, maintain twin/loop links, preserve unaffected local handles, validate and rollback on failure.
- Chain Euler operations using typed outputs created earlier in the same transaction, with one Undo/Redo step for the batch.
- Validate closed manifold incidence, vertex fans, planar non-self-intersecting loops and Euler characteristic 2.
- Configure model-space absolute/relative/angular tolerances for each solid or JSON request.
- Triangulate polygon caps and sides, retaining per-triangle source face IDs for CAD selection.
- Compute prism volume, surface area, centroid and axis-aligned bounds; translate bodies.
- Ray-pick a triangle and map the hit back to the owning topological face.
- Accept `geometry.v1` JSON commands through the `gefest-geometry` stdio executable.

This implementation is **original code** based on common computational geometry methods. It contains no Parasolid binaries, Siemens headers or decompiled implementation. The historical `PK_*` API inventory, if used for planning, is not source code or an implementation specification.

## Build and test

Requires a stable Rust toolchain (Rust 1.70+):

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

The `v0.4` model supports only a **single planar boundary loop without holes**; its B-rep is polygonal, not an exact analytic solid. There is no Boolean engine, fillets, NURBS, analytic trim curves/surfaces, per-edge tolerances, shells with holes, persistent naming across model regeneration, STEP/X_T export, or production-grade robust predicates. Extrusion uses XZ/+Y coordinates. Not suitable for safety-critical or precision manufacturing use.

Only geometry-preserving edge subdivision and internal planar face division are transactional today; deletion, vertex move, face extrusion, Boolean splitting and fillets are not implemented yet. Transactions currently clone the entire body (O(n) copy cost).

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

`predicates` provides adaptive exact-sign `orient2d`/`orient3d` (using the MIT/Apache-licensed `robust` crate), and tolerance-aware plane/plane and segment/plane classifications. Geometric intersection coordinates are approximate floating point; trimmed surface intersection, NURBS, and solid Boolean are not implemented.

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
