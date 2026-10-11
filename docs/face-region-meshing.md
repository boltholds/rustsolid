# RustSolid v0.24 — certified display tessellation of multi-loop analytic faces

This release closes a **display**, not a general Boolean, prerequisite:
the v0.23 canonical box with 0–32 strictly disjoint full-height analytical
cylindrical holes can now produce a watertight, consistently oriented,
face-owned triangle mesh, while preserving its exact analytic B-rep
and mathematical mass.

The result does NOT replace circles with polygons in `BrepModel`, nor
enable `TopologyDeltaPlan` execution. Genus N, the true 3D circle
and seam edges, planar caps with N inner trim loops, and all cylindrical
carrier / UV references remain the authoritative geometry.

## Planar region reconstruction

`FaceRegionMesher::tessellate(&BrepStructuralEditor, segments)` begins
with a full independent `BrepStructuralEditor::validate()`, including
source hole separation, canonical topology and genus, exact analytic mass,
and the v0.22 carrier/UV transformation equivalence gates.

For each cap, the mesher reconstructs one polygon-with-holes:

1. Read the actual four B-rep outer cap corners, preserving shared source
   box vertex indices.
2. Walk the structural editor's live feature bindings in dense B-rep order.
   For each one, sample its **current 3D circular edge carrier and bounded
   analytic parameter interval**, not cached material from an earlier
   parameterization. The bottom and top circles must correspond to one
   physical, vertical edge chain. Surface displacement, wrong radius or
   inconsistent top/bottom angular registration are refused.
3. Project all contour positions to a **box-local XZ plane**, not absolute
   world coordinates. Each circular inner contour is clockwise in this chart;
   the outer box contour is counterclockwise.
4. Triangulate the verified two-dimensional region using
   [`earcutr` 0.5.0](https://docs.rs/earcutr/0.5.0/earcutr/).
   This ISC-licensed library is used **only to propose display triangles**;
   neither its output nor any ear clipping result serves as a CAD Boolean
   truth oracle. The implementation was independently verified against
   the crate's documented hole-index API.
5. Reject triangles with invalid indices, zero/sub-tolerance projected
   area, wrongly oriented winding, centroids outside the outer region or
   inside an excluded polygonal hole, or an accumulated area inconsistent
   with the difference between the outer contour and all sampled holes.
6. Build the opposite cap from the same physical vertices with reversed
   winding. Both caps share exactly the ring vertex IDs that the cylinder
   wall triangulation uses.

The cap region `FaceRegionCertificate` reports its face ID, number of
loops, triangle count and *faceted* material area. It does not claim to
have exactly integrated curved cap trim boundaries.

## Closed boundary and independent mathematical certificates

The outer box walls reuse the original validated triangles and source
face IDs. Each analytic cylindrical wall uses the identical paired
bottom/top ring indices already used by the cap meshes. Cylinder
triangles are wound **inward**, toward the physical hole void.

Before publishing the mesh the kernel checks:

- Every undirected triangle edge occurs **exactly twice**, with opposite
  directed uses; watertightness is structural, not a shaded visual guess.
- Every triangle references a live face and finite distinct vertices.
  Planar patches lie on their actual support planes with outward winding.
  Cylindrical facets have all vertices on their analytic radius at the
  corresponding hole and face inward toward the removed material.
- The two cap triangle areas equal the independently computed polygon
  region area and have positive physical material measure.
- The signed tetrahedral volume of the entire oriented mesh is positive
  and matches a separate **polygonized hole volume oracle**:
  `Vfacet = Vbox − h Σ Apolygon_i`.
- The sum of triangle areas matches the independent display area formula:
  `Sfacet = Sbox − 2 Σ Apolygon_i + h Σ Ppolygon_i`.
- Display volume cannot lie materially below the exact analytic volume:
  the inscribed circle polygons remove less material than their true
  cylindrical voids.

The independently computed **analytic** mass remains
`V = Vbox − πh Σ(r_i²)` and is never overwritten by
`Vfacet`. As segments increase, the difference between polygonal and
analytical volume must converge toward zero; this is pinned by tests.

The certified container `CertifiedDisplayMesh` reports `mesh`,
`bottom_region`, `top_region`, `faceted_volume`,
`faceted_surface_area`, `analytic_volume`,
`faceting_volume_error`, and the number of circle segments.

## Geometry API

Rust:

```rust
use rustsolid::{BrepStructuralEditor,HoleParameters,Point3,Solid};
let source=Solid::block(
    Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.)?;
let mut editor=BrepStructuralEditor::new(source)?;
let mut tx=editor.begin_transaction();
tx.add_through_hole(HoleParameters{
    center_x:-1.5,center_z:0.,radius:0.5
})?;
tx.add_through_hole(HoleParameters{
    center_x:1.5,center_z:0.,radius:0.7
})?;
tx.commit()?;
let certified=editor.tessellate(64)?;
assert_eq!(certified.bottom_region.loop_count,3);
assert_eq!(editor.model().summary().genus,2);
assert!(certified.faceting_volume_error>0.);
```

Existing `geometry.v1` request `structural_transaction` gains an
**optional** `display_segments` integer:

```json
{
  "operation":"structural_transaction",
  "source":{"kind":"block","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6},
  "batches":[[
    {"kind":"add_through_hole","center":[-1.5,0],"radius":0.5},
    {"kind":"add_through_hole","center":[1.5,0],"radius":0.7}
  ]],
  "display_segments":64
}
```

When supplied, the response includes the regular `solid` rendering
payload of vertices, triangles, triangle-to-B-rep face ownership,
**true analytic topology** and exact mass. `data.display_quality`
contains the triangulation certificate, including polygonal area,
facet count, signed volume, faceting error and loop counts.
`data.display_mesh_available:true`, `display_mesh_included:true`,
`topology_committed:true` and `generalized_topology_delta_can_apply:false`
remain explicit.

When `display_segments` is **omitted**, the exact v0.23 response shape
is retained: no `solid` payload, no extra mesh computation;
`display_mesh_included:false` signals no mesh requested.
`display_mesh_available:true` now means the capability exists, not
that a mesh was requested.

Segmentation is bounded to **12–512 segments per circle**, and
`hole_count * segments <= 16,384` caps total circular tessellation
work. Invalid values and unsupported/certification-failing geometries
return `ok:false` with no partial `solid` or result data. Triangle
counts and positions are deterministic for a fixed input and segmentation,
but their dense indices are not durable CAD feature identifiers.

## Independent regression and remaining scope

`tests/multi_loop_mesh.rs` pins **12 independent cases**:

- Zero-hole box, two/three asymmetric multi-hole combinations.
- Exact loop count, owner-face/triangle association and watertight
  opposite-edge usage.
- 12/48/192-segment convergence without changing analytic mass or revision.
- Shared edge sampling after structural compaction, circle phase
  reparameterization and cap/cylinder UV frame rebasing.
- Insertion-order metamorphism and detached Arena invariance.
- Segmentation errors and near/far coordinate translations.
- Several radii, offsets and tessellation resolutions.
- The 32-hole bound with 33 planar cap loops and deterministic repeated output.

GitHub Actions runs this entire integration module and an explicit
geometry.v1 smoke/negative case in addition to all previous suites.

**Still unsupported:** arbitrary planar polygon + arbitrary curved face
trims, arbitrary body Boolean intersection regions, general shape
deformation, multiple shell components/cavities, oblique or partial-depth
holes, and generalized `TopologyDeltaPlan` application. This release
provides a **certified display triangulation for the v0.23 canonical
structural family**, not a production-quality universal mesher.
