# RustSolid v0.33 — one canonical B-rep from separately certified crossing clusters

RustSolid v0.33 implements a *bounded executable composition*, not a universal
multi-cylinder Boolean. It combines between **one and four independent
strictly transverse through/blind crossing pairs** and optionally additional
isolated blind holes within a genuine certified axis-aligned source box.
The total source budget is 12 cylinders. Blind depths may differ between
components. Every disk of one cluster must be strictly separated from every
disk of every other cluster; near tangent, connected triple overlap, or
disconnected topology ambiguities are explicitly refused.

The result is ONE validated analytic B-rep, not a pile of overlapping closed
objects, and not a polygonized replacement for exact face geometry.

## Cross-layer topological contract

Each v0.31 crossing pair is independently built as a source-certified
closed genus-one B-rep. One isolated blind hole is independently built as a
closed genus-zero donor. Donors share exactly the same eight original box
vertices, twelve edges, twenty-four directed uses and six support surfaces.
Only the donors' truly new geometry and incidence are remapped into
a single canonical shell:

- donor 3D line/circle carriers are assigned fresh stable-in-this-snapshot IDs;
- source plane/cylinder surfaces and every UV p-curve are rebased;
- all new vertex/edge/coedge/loop/face IDs are remapped, including twin,
  next/prev, trim face/loop owner and edge curve references;
- original bottom cap face 0 is given ONE genuine inner loop per through pair;
- original top cap face 1 receives ONE inner loop per crossing pair or blind
  component; blind floors attach at the correct independent world Y levels;
- every newly allocated face is material-oriented, opposite every physical
  twin use, and belongs to the *same* outer shell.

The dedicated v0.31 crossing donor already owns the split-height cylindrical
walls where through and blind materials interact. The composer does not
manufacture fake seams across disconnected clusters.

For N crossing pairs and B isolated blind holes, the number of inner face
trim loops is 2N+B, the genus is N, and the Euler characteristic is
2-2N. Example: **two distinct crossing pairs**, with two finite blind depths,
share a shell of:

- 20 vertices, 34 edges, 68 coedges, 20 loops, 16 faces, one shell;
- four true inner cap loops;
- Euler characteristic -2, genus 2.

An additional isolated blind hole adds two source-proven faces (wall/floor),
one top-cap inner loop and no additional tunnel genus.

## Independent depth-slab mass gate

Before allowing topology assembly, the v0.32 MultiCutRegionPlanner
independently reconstructs source-indexed 2D disk-union contours at each
blind floor level and integrates their analytical areas across open height
slabs. The complete resulting *global* volume must equal, within an explicit
floating-point budget, the sum of independently rebuilt donor analytical
removed volumes. Both methods derive from the same source inputs but use
different formulas and combinatorial decompositions.

Mass and centroid remain analytical and source-relative. The result
recalculates and verifies face/edge/coedge carrier bindings and source
geometry on replay. Undo/Redo operates on isolated parameter candidates;
failed edits publish no partially modified shell.

A completed v0.32 region plan with all arcs closed is NOT sufficient to
authorize a connected three-way crossing. Distinct cluster separation is an
additional acceptance precondition.

## One source-global display mesh

The display mesher does NOT concatenate independently closed donor solids.
Instead it consumes only each donor's *actual* new curved wall/floor faces,
extracts the hole loops from their already certified top/bottom facets, and
maps all donor vertices to one global index space while preserving the
source box's original eight vertices. It triangulates **one** top cap
and **one** bottom cap with all their holes together, then appends the source
box's exterior walls just once.

It then requires:

1. Exactly two oppositely oriented uses for every triangle edge across the
   entire result; no duplicate source-box caps/walls.
2. Triangle source face IDs lie inside the one canonical composed model.
3. Mesh Euler characteristic agrees with 2-2N and the expected shell genus.
4. Positive signed facet volume agrees with *independently composed*
   donor facet volumes, and facet surface area agrees with their independently
   composed area budget.
5. The displayed volume does not undershoot analytic mass by more than the
   disclosed numerical slack; angular refinement approaches analytic mass.
6. The underlying B-rep is unmodified by display generation.

The display is a certified **derived facet representation**. Analytical
cylinders, p-curves and mass stay authoritative. This does not constitute a
formal all-real-number numerical proof; the validators check finite samples
and explicit floating-point budgets.

## Rust API

~~~rust
use rustsolid::{
    Solid,Point3,HoleParameters,BlindHoleParameters,CertifiedCrossingPair,
    ComposedCrossingParameters,ComposedCrossingBrep
};

let box_=Solid::block(Point3{x:-5.,y:-1.,z:-4.},10.,4.,8.)?;
let pair=|x:f64,r:f64,depth:f64|CertifiedCrossingPair{
    through:HoleParameters{center_x:x,center_z:0.,radius:r},
    blind:BlindHoleParameters{
        center_x:x+r,center_z:0.,radius:r,depth
    }
};
let body=ComposedCrossingBrep::from_box(&box_,
    ComposedCrossingParameters{
        pairs:vec![pair(-2.,0.45,1.5),pair(1.5,0.5,2.5)],
        isolated_blind:vec![]
    })?;
body.validate()?;
assert_eq!(body.model().summary().genus,2);
let mesh=body.tessellate(32)?;
assert_eq!(mesh.genus,2);
# Ok::<(),rustsolid::GeometryError>(())
~~~

## JSON geometry.v1

~~~json
{
  "operation": "compose_crossing_clusters",
  "source": {"kind":"block","origin":[-5,-1,-4],
             "width":10,"height":4,"depth":8},
  "pairs": [
    {
      "through":{"center":[-2,0],"radius":0.45},
      "blind":{"center":[-1.55,0],"radius":0.45,"depth":1.5}
    },
    {
      "through":{"center":[1.5,0],"radius":0.5},
      "blind":{"center":[2,0],"radius":0.5,"depth":2.5}
    }
  ],
  "isolated_blind": [],
  "display_segments": 32
}
~~~

Output is a committed *dedicated* bounded shell, with analytic mass, source
slab summary, coedge/genus data and optional derived facet mesh and quality.
A batch can replace a certified pair or standalone blind source, followed
by request-local Undo/Redo. Invalid late operations, near inter-cluster
contact, or invalid segment budgets return no partial solid or result data.
This API does not enable the general Boolean entrypoint or the generic
TopologyDeltaPlan executor.

## Run the independent gates

~~~sh
cargo test --test composed_crossing_brep
cargo test --bin gefest-geometry composed_crossing_cli
cargo test --all-targets
cargo check --release --all-targets
~~~

## Not yet implemented

Directly interacting three-cutter sources (one through plus two overlapping
blinds sharing the same connected void), crossing chains across clusters,
triple junctions, nonmanifold joins, arbitrary-depth regional transitions
inside one connected cluster, oblique geometry, general NURBS or source
independent Boolean classification are still **explicitly unsupported**.

Solving them requires a shared global intersection graph with source-level
region splitting and verified 3D face ownership, rather than assembling
separated certified components. The v0.32 read-only global plan remains the
appropriate exploration and classification stage for those future cases.
