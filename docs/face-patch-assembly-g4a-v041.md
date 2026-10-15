# RustSolid v0.41 — G4a: source-owned face patch / coedge draft

**Status: experimental, read-only face-local topology preparation.**
This is **not** G4 complete and is **not** a volume-changing general Boolean.
All outputs preserve `can_apply=false`; no persistent B-rep topology is
created or committed.

## Boundary between G3 and G4

G1 verifies individual face-pair intersection graphs and reconciles material
endpoints only when physical source topology agrees. G2 constructs oriented
analytic UV cycles for one selected source face, including circular holes
without polygonizing them. G3 proves each complete UV cell uniformly inside
or outside the opposite **whole solid**, and chooses the retained cells for
Union/Intersection/Difference.

`FacePatchAssemblyPlan::build(left,right,pairs,uv,g3,tolerance)` then:

1. Source-validates both operands, G2, and G3 before inspecting patch uses.
2. Enumerates **every** directed source/cut span into `FacePatchDraftUse`,
   retaining original operand, FaceId, source EdgeId/CoedgeId/LoopId or
   `GlobalArcRef`, the face-local cell and cycle identities, and a kept
   flag. These indices are **not** persistent topological names.
3. Certifies source trim segments against the original coedge UV carrier at
   both parameter endpoints. It never invents source edge identity by
   proximity.
4. Reconstructs cut spans from the source plane in XYZ, verifies their
   bounded support **on the original analytic graph arc**, and rejects
   unsupported nonlinear intersection carriers. Each G2 cut requires two
   original uses with opposite UV orientation.
5. Requires the G3 result to retain **exactly one** of those two uses.
   If both are kept or both are discarded, it refuses the patch rather
   than converting an ambiguous material boundary into a physical edge.
6. Applies the actual original face sense and G3's Difference-B orientation
   reversal to the eventual outward-facing coedge direction and loop order.
7. Records the counterparty operand and FaceId for every cut. This is the
   **required external twin**, not a made-up local link between two
   coedges from the same source face.

The resulting plan reports `open_cut_uses`, `open_source_uses`,
`retained_loops`, `cut_pairs`, and the analytic (not triangulated)
`kept_area_uv`. `validate_against` rebuilds G2 and G3 from the
current source and rejects changed bodies, classifications, operations or
tampered draft coordinates.

## Deliberate limitations and subsequent work

The two local uses of a G2 cut are *opposite sides of a subdivision in
ONE source face*. They are **not** opposite radial coedges of a physical
edge in the final B-rep. Their provisional identity must not be treated as
a physical EdgeId or topological twin.

**G4b** must schedule and classify **every affected face** of both
operands, split graph arcs consistently across chart boundaries,
reconcile same-XYZ vertices only under G1 source evidence, and pair
matching physical cut halves across the *two different source faces*.
It must also stitch unchanged source coedges, close all inner loops,
validate orientation and manifold radial valence, and produce a
fully checked candidate B-rep.

The current G2 supports bounded full-face planar line traces. It may
return Held for a trimmed intersection that ends halfway along the
opposing face, so the availability of one G4a source patch does **not**
authorize a full global candidate. Curved source surface patches and
tangencies also remain unsupported here.

**G4c** needs an independent analytic volume oracle and shell/genus
checks for the first original hole + crossing cutter fixture. **G5**
introduces atomic revision-bound commit and Undo/Redo only after G4b/c
produce a certified candidate.

No production release may claim G4/G5 or general Boolean complete
based on the existence of this G4a module.

## Source and regression

- `src/face_patch_assembly.rs`
- `tests/face_patch_assembly.rs`

```bash
cargo test --test face_patch_assembly
cargo test --all-targets
```

Tests exercise the existing analytic through-hole + transverse box
cutter for **both Difference operand orientations**, whole-face
retention, directed closed loop reconstruction, source/graph carrier
consistency, modified operations, stale source bodies, cut-ID forgery,
and the non-applicable `can_apply=false` boundary.

The G1/G2/G3 predecessors remain source-checked and immutable.
