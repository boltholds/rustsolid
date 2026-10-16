# RustSolid v0.41 G4b2b3b: bounded analytic curved patch partition

Status: source-certified preview of a real four-shell curved B-rep for
canonical AABB minus a full-height upright cylinder whose disk crosses all
four box walls and leaves four separated corner regions.

Each corner has 6 physical vertices, 9 edges, 18 coedges, 5 analytic faces
(two planar caps, two planar walls, one inward-oriented cylindrical wall),
one connected genus-zero closed shell. Totals: V=24, E=36, F=20,
4 components, Euler characteristic 8, total genus 0.

Evidence chain:
- G4b2b2 recomputes the entire source face-pair graph, real G1 physical
  vertex events and two source-owned analytic UV pcurves per crossing arc.
- G4b2b3a subdivides both original source B-reps at source-certified events.
- All retained original material edges match source-subdivided physical
  edges with original source-face ownership.
- All 16 new cut edges match unique global analytic source arcs, with
  carrier samples and real source event identities.
- Every retained patch has a positive retained-material side witness and
  a negative removed-material side witness, using strict analytic point
  membership. Ambiguous near-boundary conditions are Held.
- Every resulting edge has two opposite physical coedges; each material
  vertex fan is audited independently and all 4 shell components are
  validated without geometric proximity welding.

Only the new AnalyticPartialBoxCylinderDifference origin permits four
shells in BrepModel::validate(). Previous producers retain the one-shell
contract. Multi-shell validation checks per-shell connectedness, face
ownership, shared-vertex/edge exclusions, Euler and genus.

The source-specific analytical mass/area/centroid values are preview
metadata. Their comparison against independent numerical quadrature in
tests does not substitute for the general independent G4c volume oracle.

Unsupported cases include tangencies, tilted or recharted cylinders,
general curved UV arrangements, multiple interacting cutters, and
arbitrary trimmed surfaces. The result is a real bounded-family geometric
candidate, NOT a general curved Boolean. G5 atomic source-revision
transaction and Undo/Redo are not implemented.
closed_shell_proven=true only after all source checks; can_apply=false.

Tests: cargo test --test analytic_face_patch_partition, then cargo test
--all-targets. Keep PR draft until G4c/G5 are implemented and reviewed.
