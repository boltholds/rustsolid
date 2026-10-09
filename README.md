# RustSolid

RustSolid is an **independent Rust solid geometry kernel** being developed for [Gefest CAD](https://github.com/boltholds/Gefest-CAD).

> **Status: v0.1 prototype — planar, genus-zero polyhedral solids only.** It is not currently a replacement for Parasolid or OpenCascade. Its API does not claim compatibility with proprietary kernels.

## Current capabilities

- Extrude a simple, possibly concave, planar XZ profile along +Y, or create a block.
- Maintain explicit vertex/edge/face topology with oriented face boundaries.
- Check closed, consistently wound manifold edge incidences and Euler characteristic 2.
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

The response carries `schema_version: "geometry.v1"`, `ok`, and either a `solid` or an `error`. `solid` includes `vertices`, indexed triangle `faces`, `triangle_face_ids`, counts and Euler characteristic under `topology`, analytic `mass` and `bbox`. `geometry.v1` is a provisional bridge contract, not a STEP or Parasolid serialization format.

Example library use:

```rust
use rustsolid::{Point2, Solid};

let profile = [
    Point2 { x: 0.0, z: 0.0 },
    Point2 { x: 10.0, z: 0.0 },
    Point2 { x: 10.0, z: 5.0 },
    Point2 { x: 0.0, z: 5.0 },
];
let solid = Solid::extrude_xz(&profile, 8.0)?;
assert_eq!(solid.euler_characteristic(), 2);
assert_eq!(solid.mass.volume, 400.0);
# Ok::<(), rustsolid::GeometryError>(())
```

## Gefest integration

The existing Gefest backend uses CadQuery/OpenCascade, and its geometric constraint solving is handled by a separate Rust workspace (`solver-core`). RustSolid is deliberately independent of both. When using the experimental Gefest Python adapter, set `GEFEST_GEOMETRY_BINARY` to the compiled `target/release/gefest-geometry`; opt in with `GEFEST_GEOMETRY_BACKEND=rust`. The OpenCascade path should remain the production default until parity, reproducibility and regression gates pass.

## Known limits

The `v0.1` model supports only a **single planar boundary loop without holes**; its B-rep is polygonal, not an exact analytic solid. There is no Boolean engine, fillets, NURBS, analytic trim curves/surfaces, per-edge tolerances, shells with holes, persistent naming across changing topologies, STEP/X_T export, or production-grade robust predicates. Extrusion uses XZ/+Y coordinates. Not suitable for safety-critical or precision manufacturing use.

See [docs/architecture.md](docs/architecture.md) for module boundaries and an incremental roadmap.

## License

Mozilla Public License 2.0 (`MPL-2.0`). See [LICENSE](LICENSE).
