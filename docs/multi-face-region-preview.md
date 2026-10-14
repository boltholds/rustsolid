# RustSolid v0.39 — all-or-nothing multi-face region preview

This is a read-only generalization of the **source-provenance and local
region classification** pipeline. It is not a global face-cell
partitioner, topology executor, or generalized Boolean.

## API

```rust
let preview = MultiFaceRegionPlan::plan(
    &left, &right, BooleanOperation::Difference,
    GeometryTolerance::default(),
)?;
preview.validate_against(&left, &right, GeometryTolerance::default())?;
assert!(!preview.can_apply);
```

The planner validates both closed B-rep operands, evaluates every
supported pair of finite source faces through
`IntersectionGraphBuilder`, and classifies the local arc-side samples
via the existing `BooleanRegionClassifier`. A successful
`ClassifiedPreview` records each source-face pair, its mathematically
grounded `IntersectionGraph`, and its local `GraphRegionPlan`.
Source face IDs remain **operand-namespaced**. Empty graph pairs do not
invent boundary edges.

To remain bounded, one call supports at most 64 faces per source body
and at most 512 face pairs. These limits are on the mathematical
preview, not a guarantee of eventual Boolean complexity.

## All-or-nothing Held policy

The planner returns a typed `Held` with **no partial face pairs, arcs
or classification samples** when any of the following occurs:

- an input face has an inner UV loop or a non-single-loop trim domain;
- a face pair is unsupported by current analytic intersection families;
- a pair has unresolved tangency, overlap or numerical ambiguity;
- an isolated physical contact would otherwise look like a real cut;
- a local side-probe cannot classify the separating region;
- face/pair work exceeds the fixed budget.

Existing source B-reps are immutable throughout. `can_apply` is
**always false**; a separate topology editor cannot interpret this
preview as permission to change material.

The existing `TrimmedIntersectionEngine` supports only selected
planar polygons/disks and full periodic cylinder-side rectangles.
Proper inner-loop clipping, multi-arc arrangements, and global face
cells are **not** implemented by this milestone.

## Verified fixtures

`tests/multi_face_region_preview.rs` provides:

1. truly disjoint box and cylinder => no false cut;
2. a finite box intersecting a cylinder across several real faces =>
   multiple per-face arc/classification witnesses;
3. a through-hole B-rep with genuine inner cap trim loops => typed
   `Held(UnsupportedInnerLoop)` and zero leaked partial reports;
4. near-tangent input => typed `Held` and zero leaked partial reports.

The validated plan is deterministic and recomputable against the same
source bodies. A forged count breaks `validate_against`.

## The remaining hard problem

Pair-local arc-side truth values are insufficient to choose **all**
material regions on a complex trimmed surface. Next steps require
complete 2D face arrangements, embedded intersection arcs, winding
rules for all outer/inner trim loops, explicit region-cell adjacency,
and globally consistent labeling across multiple source faces.
Only then can a general `TopologyDeltaPlan` executor construct and
glue shared edges/coedges and independently certify the resulting
closed B-rep and mass. No proprietary kernel code is used.
