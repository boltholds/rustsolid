# RustSolid v0.41 — G4b1: source-certified cross-face edge stitching

**Status: guarded global incidence preview, not a Boolean B-rep commit.**

This step consumes several already certified, read-only G4a
`FacePatchAssemblyPlan` objects and their own source `FacePairGraph`
sets. It never creates persistent `EdgeId`, modifies an original
`BrepBody`, assigns any mass properties, or treats G2's two local cut
sides as a radial pair.

## Physical coedge matching

`GlobalEdgeStitchPlan::build(left,right,inputs,tolerance)` reconstructs
and validates G2/G3/G4a for every source patch before pairing any use.
Duplicate input faces and inconsistent operations are rejected.

Two matching families are supported:

**Original source edge**: both selected uses come from the SAME operand
and SAME source `EdgeId`. Their source coedges must be actual
`coedge.twin` of one another in the validated canonical `BrepModel`,
with opposite directed XYZ endpoints. Distinct source EdgeIds are not
merged because they happen to be geometrically coincident. Matching
subspans may later become separate result edges, but retain their
original edge/coedge provenance.

**New intersection edge**: two selected uses originate from DIFFERENT
operands. They must name the SAME actual pair of source faces, the SAME
analytical graph arc when regenerated against those sources, and
opposite directed XYZ endpoints. Each source UV span remains separate.
Unlike an arbitrary geometric proximity test, an arc's face-pair
identity and bounded original 3D carrier are mandatory. The first
bounded gate supports an entire certified `IntersectionCurve::Line`
graph arc; an unverified mid-arc subdivision cannot silently create a
new material vertex.

The two endpoints of such a new physical edge must also have original
boundary witnesses on BOTH operands in the same graph. Uncertified
endpoint sharing, partial arcs and same-direction uses return
`Held`, not a fabricated edge.

Each successful pair is recorded as `StitchedEdgeDraft` with two
different `StitchFaceOwner` sources, two distinct original UV
parametric boundary spans, a common XYZ carrier and a **provisional
draft ID**. Nothing published here is a persistent kernel `EdgeId`.

## Radial incidence is not shell closure

`GlobalEdgeStitchPlan` audits directed source loops, counts each
retained face coedge exactly once, verifies continuity, and records
all paired physical edge witnesses. It separately reports:

- unmatched `open_uses`, with the missing source coedge or missing
  opposing certified cut patch,
- `missing_source_faces` for every source body face not proven by a
  G4a plan,
- `radially_paired` (all selected uses have exactly one paired use).

The scope of the G2/G3 face-local classifier means that real parts
will often be only partially certified: G2 requires a full transverse
planar cut, and cannot yet publish unaffected or partially trimmed
faces. These are **explicit outstanding faces**, never presumed
discarded. No shell is labeled closed based on radial matching alone.

To claim a valid closed B-rep, G4b2 must also provide certified
unaffected-face selection, split intermediate graph arcs with a globally
reconciled vertex identity, verify single connected manifold vertex
fans, close periodic surface seams and all trim loops, and construct
and validate the actual canonical `BrepModel` with shell/genus checks.
Then G4c must supply an independent analytic volume oracle. Only G5
may atomically commit a revision-bound candidate with Undo/Redo.

Until every gate is implemented and executed,
`closed_shell_proven=false` and `can_apply=false`.

## Tests

`tests/global_edge_stitch.rs` covers:

- one valid G4a face cannot masquerade as a closed global shell,
- duplicate source faces and forged source coordinates fail closed,
- both-operand certified face crossing, including exact opposite
  edge orientation and distinct face-local chart witnesses,
- empty input cannot produce an empty success candidate.

```bash
cargo test --test global_edge_stitch
cargo test --all-targets
```

The downstream RustSolid G4a interface is unchanged, and neither
operation touches the project's primary `main` branch.
