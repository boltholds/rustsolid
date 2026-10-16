# RustSolid v0.41 — G4b2b1: source-certified analytic bore B-rep candidate

**Status:** executable, immutable, bounded curved-solid candidate for the
two supported analytical Boolean Difference families:

1. Certified axis-aligned box minus a strictly interior **full-height
   upright cylinder**, producing a real closed genus-one through-hole.
2. Certified axis-aligned box minus a strictly interior **top-entry
   upright cylinder of partial depth**, producing a real closed genus-zero
   blind bore with a separate analytic planar floor.

The actual analytic topology constructor is **reused**, not copied or
approximated. `ThroughHoleBrep` and `BlindHoleBrep` already own exact
original cylinder/plane carriers, original-source physical edges, two
face-specific analytic UV pcurves on each circular rim, the periodic
cylinder seam, and their canonical validators. The new
`AnalyticBoreCandidate::build` is a source-bound *independent verifier*
of those real closed candidates. This is not arbitrary B-rep Boolean.

## Validation pipeline

- Validate both original operands and their geometric tolerance policy.
  Reject non-AABB inputs, non-upright cylinders, tangency, near-coincident
  endpoints, mismatched caps and unsupported partial contact without
  substituting meshes or bounding boxes.
- Invoke the already-certified analytic Boolean kernel, accepting only its
  true `AnalyticThroughHole` or `AnalyticBlindHole` B-rep output.
- Independently build source-bound **face-pair intersection graphs**
  for both original box caps against the original cylinder side in the
  through case, or for the upper opening in the blind case. Each graph
  must contain one true closed analytic arc and one *virtual* periodic
  chart anchor. Such anchors are not welded across different face pairs
  into fictitious physical material vertices.
- Persist the independently derived, source-owned
  `AnalyticSourcePcurvePair` for each closed face-pair graph. These store
  a `Curve2::Harmonic` on the original box cap and a
  `Curve2::CylinderSection` on the original cutter's periodic wall, with
  their source graph parameter interval. Both curves re-evaluate to
  the same 3D arc at multiple certified parameter values, even when
  the cutter's original angular UV frame is rotated. These are
  provenance witnesses, not newly allocated coedges.
- Verify the generated physical lower and upper circular edges
  (`EdgeId(12)`, `EdgeId(13)`) against the source cylinder radius and
  correct source cap or newly constructed blind floor plane. On each
  edge, evaluate BOTH independent pcurves (a `Circle2` in the planar
  face and a `Line2` on the cylindrical wall) and the physical
  `Circle3` carrier at multiple angles; validate source cylinder
  graph arcs in their own original angular charts.
- Validate the same-face periodic seam `EdgeId(14)` as one *physical*
  carrier with two coedges in distinct UV charts (`u=0` and
  `u=2π`), never as two disconnected radial edges. Source graph
  periodic anchor identities do not become stable result IDs.
- Compare all **12 source box edges** and **24 original coedge** carriers
  with the produced model. There are **9 explicit source face ownership
  records**: 6 original box faces, 1 resulting cylindrical wall from
  source right FaceId(2), and the two original cylinder caps (one
  becomes the blind floor when applicable; others are discarded).
  Opened box faces have true circular inner loops; face senses of the
  resulting cylinder wall and blind floor are reversed to face the void.
- Independently recompute constrained analytic **volume, surface area
  and centroid** from original box dimensions and the actual cylinder
  radius, axis and length. This is the mathematical oracle for the
  bounded upright-bore geometry, not general G4c integration.
- Check canonical B-rep `validate()`, shell ownership, closed loop
  incidence and genus: through-hole V=10 E=15 F=7,
  Euler characteristic 0, genus 1, 2 planar inner loops; blind-hole
  V=10 E=15 F=8, Euler characteristic 2, genus 0, one inner loop and
  a planar circular floor.

## Source revalidation

`validate_against` verifies the original left source edit revision,
the complete right cylinder frame (including angular UV chart),
radius/height and the construction tolerance. It then independently
recreates the candidate and checks geometric carriers, source
provenance, each closed face-pair graph, circle/rim witnesses and the
periodic seam.

Both original source bodies, their revisions and their undo/redo
histories remain untouched; `can_apply=false` for every candidate.

## Limitations — next G4b2b2 / G4c / G5

This **is not yet a global analytic B-rep Boolean executor** for
two arbitrary solids or multiple intersecting cylinders.

- The planar box/cylinder-through and box/cylinder-blind families are
  already implemented by existing dedicated RustSolid producers. We
  verify them independently rather than claiming these constructors
  were written afresh.
- A partial cylinder that intersects a box *side* or crosses another
  cutter is held, not clipped by facets.
- Source face-pair graphs of the two rim intersections keep their own
  virtual periodic anchors; true physical vertex reconciliation across
  multiple *open* cut arcs still requires G1 event/corner witnesses
  and globally scheduled G2/G3 patches.
- A general B-rep candidate builder must support arbitrary partially
  trimmed analytic support faces, physical edge splits at globally
  certified events, periodic seam chart migration, nested loops and
  full manifold vertex-fan incidence.
- A general G4c volume oracle and G5 atomic document transaction
  with stable naming and Undo/Redo remain explicitly outstanding.

## Verification

`tests/analytic_bore_candidate.rs` covers full-height and strict
partial-depth cutters, true circle carriers, both result-face and
original-source independent pcurves, source graph evidence, complete original face lineage, a
nonmaterial periodic seam, off-center mass and centroid, UV frame
rotation, failures for tangent/unsupported cuts, and stale source
or tolerance policies.

```bash
cargo test --test analytic_bore_candidate
cargo test --all-targets
```

This is a stacked feature branch and intentionally does not alter
RustSolid `main`.
