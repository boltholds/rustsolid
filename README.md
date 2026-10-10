# RustSolid

RustSolid is an **independent Rust solid geometry kernel** being developed for [Gefest CAD](https://github.com/boltholds/Gefest-CAD).

> **Status: v0.4 experimental transactional topology slice — planar, genus-zero polyhedral solids only.** It is not currently a replacement for Parasolid or OpenCascade. Its API does not claim compatibility with proprietary kernels.

## Current capabilities

- Extrude a simple, possibly concave, planar XZ profile along +Y, or create a block.
- Maintain explicit oriented B-rep: shells, faces, loops, coedges (twin/next/prev), edges and vertices.
- Edit the polyhedral B-rep atomically: split edges and faces, maintain twin/loop links, preserve unaffected local handles, validate and rollback on failure.
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
