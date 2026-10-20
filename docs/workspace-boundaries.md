# RustSolid workspace boundaries — RFC 0001

Status: **in progress**. This branch starts the restructuring on top of
`feat/g8c-tilted-one-wall-analytic-brep-20261019`.

## Breaking-change policy

We do not preserve the old public `rustsolid` API for Gefest or third-party
consumers. Change callers rather than introduce versioned aliases, facades,
old field names, or a permanent `legacy` crate. Existing G4–G8 experiments
are retained as behavioral and geometric evidence, not as frozen Rust API
contracts. Do not silently delete the *underlying algorithms* when an API is
removed.

## Crate dependency DAG (target, not yet all implemented)

```text
foundation
   ├── predicates                  # filtered/exact-sign numeric decisions
   └── geometry                    # parameterized curves and surfaces
           ├── uv                   # Curve2, polygonal adapter and arrangement
           ├── intersect            # curve/curve, curve/surface, SSI
           └── topology             # B-rep incidence and generation-aware IDs
                ├── validation      # independent geometry/topology oracles
                ├── meshing         # display/CAE facets, face provenance
                └── boolean         # depends also on intersect, uv, validation
                        └── document # transactions, history, persistence
                                      └── gefest-geometry API adapter
```

Dependencies must go downward only: in particular, `geometry` must never
depend on `topology`, `topology` must never depend on `boolean`, and
`validation` must not depend on `boolean`. An independent mass/flux oracle
cannot read the candidate Boolean implementation it is meant to check.
No crate should import or re-export an unrelated high-level subsystem to
work around this constraint.

## Implemented in this branch

- `rustsolid-foundation`: existing point types, error and tolerance policy.
- `rustsolid-predicates`: adaptive orientation, local-scale classification,
  analytic supporting-plane calculations. It has no Solid/B-rep dependency.
- `rustsolid-geometry`: analytic curves and surfaces, parameter domains, frames
  and bounded geometry IDs. The topology-backed `GeometryStore` stays
  in a higher layer as an adapter until its ownership is redesigned.
- `rustsolid-topology`: vertex/edge/face/coedge/loop/shell incidence records,
  typed topology IDs and issue types. Face-to-coedge assembly now runs on a
  scratch copy and commits face loops only when it succeeds. Solid-owned Euler
  operators, arena mutations and validators remain for a subsequent cut.
- `rustsolid-intersect`: untrimmed analytic supporting-surface intersections.
  Does not depend on Solid, coedges or Boolean; trimmed intersections remain
  in the topology-aware layer. The pinned 14-case analytic, metamorphic and
  boundary-intersection test corpus now lives inside this owning crate.
- `rustsolid-uv`: polyline-only polygonal Boolean through `geo`, using
  `rustsolid-geometry::Point2Param` instead of a duplicate UV point type; validates
  input and output polygons. It neither stores analytic curves nor determines
  final 3D topology. Its output MUST NOT be committed as certified B-rep
  without source-provenance and geometric-error reconciliation.
- Root `rustsolid`: temporarily holds G4–G8 behavior and still owns the
  topology-aware `face_support_plane` helper. It is a migration workspace
  member, NOT a commitment to the old public API.

The old root `Point2` uses `(x,z)` because that is the current local
coordinate convention. New UV code uses explicitly named `(u,v)` coordinates.
This avoids a hidden axis swap when integrating `geo`.

## Migration ownership

| Existing module(s) | Target crate |
| --- | --- |
| `geometry.rs`, analytic `Curve2/3`, surfaces | `geometry` — analytic part moved; `GeometryStore` adapter remains |
| `intersections.rs`, `trimmed_intersections.rs` | `intersect` — pure supporting-surface part moved; trim adapter pending |
| `planar_uv_arrangement.rs`, `multi_trace_uv_arrangement.rs` | `uv` |
| `brep.rs`, `unified_brep.rs`, `coedge_reconciliation.rs`, `brep_arena.rs` | `topology` — entity records moved; validation/assembly remain to migrate |
| `global_vertex_reconciliation.rs`, `global_edge_stitch.rs` | `topology` / `boolean` boundary |
| `boolean.rs`, `boolean_region.rs`, cut plans | `boolean` |
| `independent_*volume*.rs`, `radial_shell_audit.rs` | `validation` |
| `face_region_mesher.rs`, `*_mesh.rs` | `meshing` |
| `atomic_boolean_document.rs`, `g5*`, `g7*archive.rs` | `document` |
| `src/bin/gefest-geometry.rs` | `apps/gefest-geometry` |
| Regression fixtures & performance experiments | `rustsolid-bench` |

Some modules need to be split, not moved wholesale. In particular the old
`geometry.rs` reads Solid data; separate the pure analytic types from the
Solid-to-geometry projection. The latter belongs to a higher layer.

## Dependency audit

- Already present: `robust` for orient2d/orient3d, and `earcutr` for display
  triangulation. The newly isolated predicates crate now owns `robust`.
- Newly trialed: `geo`, in a limited polygon-only package. Do not replace
  circular-arc source geometry with sampled polylines.
- Evaluate `rstar`, `truck-geometry`, `truck-topology`, `truck-shapeops`
  against behavior on G7/G8 before introducing them in production code.
  No dependency is accepted merely because it has a similarly named API.

## Verification gates

Run the CI-enforced layer dependency policy before tests. A new crate must be
registered explicitly in `tests/check_workspace_boundaries.py`.

```bash
python3 tests/check_workspace_boundaries.py
cargo test -p rustsolid-foundation --all-targets
cargo test -p rustsolid-predicates --all-targets
cargo test -p rustsolid-geometry --all-targets
cargo test -p rustsolid-intersect --all-targets
cargo test -p rustsolid-topology --all-targets
cargo test -p rustsolid-intersect --test surface_intersection_contracts
cargo test -p rustsolid-uv --all-targets
cargo test --workspace --all-targets
```

For every transfer of G4–G8 functionality, assert:
1. Exact source curve and pcurve coedge ownership are preserved or explicitly
   unavailable; no polygonized canonical 3D geometry.
2. Each closed shell is manifold with correct orientation.
3. Face/edge provenance is unambiguous or operation is rejected.
4. Independent volume/flux oracle encloses the computed mass.
5. A failed operation leaves the current document unchanged.
6. Time, allocations and compile/test time are measured per crate.

CI independently tests the low-level crates. The root regression corpus
remains active until the corresponding tests are rewritten against their
new owning crate. A working G8c result is required before removing the
old implementation.
