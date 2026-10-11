# RustSolid v0.20 — source-bound Boolean region classification

This milestone moves the generalized surface Boolean pipeline from *where
would the cut be?* to *which side of the cut belongs to each result?* It is
**read-only evidence for future topology edits**. Exact box Boolean execution
from v0.19 continues to work unchanged.

## 1. Three-valued containment of complete solids

`BooleanRegionClassifier::classify_point(&BrepBody, point)` returns one of:

- `Inside` — strictly within the validated body.
- `Outside` — strictly outside the validated body.
- `Boundary` — within the local positional resolution of its boundary.

The result cannot silently reclassify an on-boundary point as inside or
outside. Invalid inputs and numerically undecidable solid-angle sums return an
error.

For a closed planar-polyhedral `Solid`, signed solid angles of **all oriented
triangles** provide a volumetric winding decision after testing actual
triangle faces for boundary proximity. This supports concave closed
polyhedra; the bounding box is used solely as a conservative quick rejection,
never as a stand-in for the shape.

For a `CylindricalBrep`, membership is determined analytically in its
verified local orthonormal frame, using the radial inequality
`hypot(x,y) <= radius` and axial bounds `0 <= v <= height`. The method also
works when the cylinder axis is rotated or the source is translated.

The polyhedral winding is numerical, not adaptive exact volume orientation
for all possible conditions. Unresolved cases are errors rather than invented
answers. Only the two existing `BrepBody` source variants are supported.

## 2. Local source-face decisions along intersection arcs

`BooleanRegionClassifier::classify_graph(graph, left, left_face,
right, right_face, operation)`:

1. Converts both actual bodies to their shared B-rep representation and
   **revalidates the intersection graph** against its original face sources.
2. On each crossing arc, evaluates its exact-in-form analytic curve and
   differential tangent at fractions 0.2, 0.5, and 0.8.
3. On **each** source face, constructs a pair of small probes across the
   intersection direction but tangent to the source surface. Both probes are
   projected onto the actual plane or cylinder carrier and must lie in the
   source face's existing trim domain.
4. Classifies both probes against the **opposite full-volume body**.
   Each genuine transverse cut must produce one inside and one outside
   witness. Touches, ambiguous near-tangents and unavailable trim clearance
   produce a `Held` plan, never guessed cut topology.
5. Applies a Boolean truth table to each probe:

| Operation | Portion of A to retain | Portion of B to retain |
| --- | --- | --- |
| Union | Outside B | Outside A |
| Intersection | Inside B | Inside A |
| Difference (A−B) | Outside B | Inside A, with orientation reversed |

These decisions are emitted as `ArcRegionSample` with operand-local face
identity, curve position, both source-surface probe coordinates, opposite
volume classifications, keep/discard flags and orientation-change intent.

`GraphRegionPlan` has three states:

- `NoCut` for an empty graph or isolated zero-dimensional contact;
- `ClassifiedPreview` for supported transverse curve samples;
- `Held` for unresolved geometric classification, tangency, face-boundary
  overlap, or insufficient separation.

It always reports `can_apply() == false`.

## 3. Semantic trust and limits

`GraphRegionPlan::validate_against` regenerates the entire plan and
rechecks the original `IntersectionGraph` against both current B-reps.
Changing keep flags, source face identifiers or witness coordinates invalidates
the proof. No source body, arena registry, Undo/Redo history, or
`TopologyDeltaPlan` is modified during classification.

**This is LOCAL region evidence for one face pair.** Three samples per
intersection arc cannot certify all patches of a complex multi-face Boolean,
resolve hidden intersections with unrelated faces, complete open traces or
guarantee fully labeled 2D region arrangements. Even a
`ClassifiedPreview` does **not** mean a complete solid Boolean is executable.

This version does not introduce inner UV trim loops, multi-shell body
validation, mutable analytic B-rep carriers, Boolean feature persistence or
the generalized `TopologyDeltaPlan` executor. Those remain required gates.

## 4. geometry.v1 read-only CLI

```json
{
  "operation": "boolean_region_plan",
  "mode": "difference",
  "first": {
    "source": {"kind":"block","origin":[-3,2,-3],"width":6,"height":1,"depth":6},
    "face_id": 0
  },
  "second": {
    "source": {"kind":"analytic_cylinder","origin":[0,0,0],"radius":2,"height":5},
    "face_id": 2
  }
}
```

Returns `ok:true`, `data.kind="boolean_region_plan"`,
`data.state="classified_preview"`, `intersection_arcs:1`, and six
per-face-region witness samples. The response also contains the independent
source-validated `topology_intent`: symbolic vertices/arcs/face-change
counts and **remaining full Boolean requirements**, including
`region_classification` because sampling a single face pair does not certify
global region-arrangement completeness.
`topology_committed:false`, `can_apply:false`, and **no solid result**
are mandatory. `Held` results carry a typed `hold_reason`.

Run:

```sh
cargo test --test boolean_region
cargo test --all-targets
cargo check --release --all-targets
```

Tests pin concave winding against false bbox classifications, boundary
semantics, exact cylindrical membership in rotated frames, two open
plane/cylinder generators, one periodic closed intersection circle,
Union/Intersection/Difference dual-side decisions, tangency and near-tangent
holding, source-evidence tampering, and large-world translation robustness.
GitHub Actions includes a dedicated test gate and a JSON API smoke test.

## 5. Next executable topology gates

1. Formal multi-loop/multi-shell B-rep validation, including hole orientation,
   genus, connected components, containment and manifold fan checks.
2. A mutable shared carrier store for curved 3D edges, per-coedge UV curves
   and face surface orientation, with scoped generation-based handles.
3. Region arrangement completeness across **all** intersecting face pairs,
   not merely local arc samples. This must handle tangencies, coincident
   supports, overlapping boundaries, open traces and periodic seam charts.
4. Strong-atomic construction of faces, loops, edges and shells with source
   lineage, independent mathematical mass/closure witnesses, Undo/Redo,
   and failure-injection regressions.
5. Actual subtraction of a fully contained cylindrical through-hole from a
   block as a first constrained analytic Boolean **only after** the above
   topology requirements have independently passed.

The conservative v0.16 `TopologyDeltaPlan::can_apply == false` guard
remains in force throughout this phase.
