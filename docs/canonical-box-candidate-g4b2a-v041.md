# RustSolid G4b2a — Closed canonical candidate for exact AABB Boolean

**Scope:** executable, source-audited closed result B-rep only for two
certified axis-aligned rectangular solids with distinct, well-separated
supporting planes. This is a specialized continuation of the G4b1
topological work, **not** a general multi-surface Boolean or a G5 commit.

## Architecture and reuse

The existing `BooleanKernel::apply_solids` already implements an exact
Cartesian halfspace arrangement for two original AABBs and constructs a
separate planar `Solid` with physical EdgeId/CoedgeId twins and
`BooleanFaceSource` lineage. G4b2a deliberately reuses that **existing
source-proven planar builder** as its geometry executor instead of
pretending the current G2 planar cut algorithm can classify every
cylindrical/partial source surface.

The independent `CanonicalBoxCandidate::build` gates this executor:

1. Both operands are validated and proven to be **actual AABBs**, not
   arbitrary meshes whose bounding boxes were substituted as solids.
2. Every source supporting plane coordinate from BOTH bodies is
   separated from every other plane on the same axis by more than the
   declared resolution. No source boundaries are merged by proximity.
3. The independent face-coverage auditor enumerates all twelve original
   planar faces. For each, the Cartesian arrangement partitions its
   UV rectangle using only exact source-plane coordinates. Both sides
   of every open rectangle are classified by the two box halfspaces.
   The operation's truth table determines whether the tile is an
   exposed physical result face and its outward orientation.
4. EVERY produced result quad must map **exactly once** to the expected
   source face/tile/normal. Duplicates, missing pieces, extra surfaces,
   unsupported orientations and incorrect Difference-B reversal fail
   the candidate. No omitted original face is implicitly treated as
   fully discarded.
5. Every new result VertexId is validated against the original exact
   three-axis source-coordinate lattice. A unique tuple of three
   coordinate-plane indices binds every material vertex; equal XYZ
   derived from different topological fans are not welded by epsilon.
6. A **real** canonical `BrepModel::from_polyhedral` is constructed,
   including planar Surface3 and Curve2/Curve3 stores, physical
   shared edges, distinct face-specific coedges and closed Shell.
   `BrepModel::validate` checks UV-to-3D reconstruction and manifold
   incidence.
7. The preexisting `certify_radial_shell` independently traverses the
   ACTUAL candidate's coedges and face loops. It proves exactly two
   radial uses for every real edge, closed oriented vertex fans, no
   coincident disconnected topological vertices, one component,
   Euler characteristic 2 and genus 0.
8. The result volume is cross-checked with the elementary *AABB-only*
   inclusion-exclusion expression derived from original box extents.
   This is **not** a curved-solid or general G4c mass oracle.

The result is a genuine immutable `CanonicalBoxCandidate` owning a
validated `BrepModel`, typed `SourceFaceCoverage` for ALL original
faces, exact `CandidateVertexLattice`, and a
`RadialTopologyCertificate`. The `closed_shell_proven` predicate is
true for this **restricted planar fixture only**.

`can_apply=false`: neither the user's document nor the original
operands, revision/identity incarnation, or Undo/Redo stack is touched.
`validate_against` regenerates the full candidate from the current
sources and rejects changed result geometry, provenance or classification.

## Source face dispositions

Every original face is classified as:

- `RetainedUncut`: the entire original face area remains a result
  boundary, even if the exact Cartesian grid subdivides it into quads;
- `Discarded`: no physical boundary area is retained;
- `Partitioned`: a strict positive subset survives. Every tile is
  independently accounted for and area conservation is checked.

The actual Boolean surface uses must have the source `FaceId` and
`BooleanOperand` of that tile, including the correct reversed
outward sense for subtractor faces in Difference.

## Deliberate Held cases and remaining G4b2 work

- Non-AABB polyhedra, cylindrical sources and generic NURBS carriers
  return `Held(UnsupportedOperand)`, never fake boxes derived from a
  loose bounding box.
- Exact-coincident or resolution-close source supporting planes
  return `Held(NearCoincidentSourcePlanes)` in this strict first slice.
- An empty intersection has no closed nonempty candidate and returns
  `Held(EmptyResult)`. Multiple disjoint or nonmanifold result shells
  remain unsupported until independently validated.
- The existing box engine's unsupported cavities/components and
  unsupported geometrical conditions remain Held, not partial bodies.

This adapter is not yet driven by **all** G1/G2/G3/G4a plans for each
surface. It demonstrates that G4b's topological audit can certify an
actual complete B-rep without weakening any of the existing analytical
provenance gates. Generalizing from complete AABB tiles to B-rep
multi-trace curved patches requires source-face scheduling, consistent
global cut-vertex identities across all G1 graph pairs, cylinder
pcurves, periodic seam ownership, robust candidate mass integration,
then G5 atomic document commit.

## Tests

`tests/canonical_box_candidate.rs` checks three Boolean operations on
strict overlapping boxes, unchanged/partitioned/discarded original
faces, twelve-face coverage, exact coordinate-grid vertex identities,
closed manifold topology, Difference-B surface senses, tamper/rebuild,
non-box refusal, near/coincident plane refusal and empty intersection.

```bash
cargo test --test canonical_box_candidate
cargo test --all-targets
```

All work remains stacked on the experimental G4b1 PR, without
modifying the primary `main` branch.
