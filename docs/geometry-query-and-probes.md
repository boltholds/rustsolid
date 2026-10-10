# RustSolid: boundary-case accuracy probes and geometry query language

## Why this split exists

The kernel must be judged by **independent oracles**, not simply by its own
`validate()` method returning `Ok`. Exhaustive Cartesian products of model
parameters are prohibitively large and mostly uninformative. The first probe
suite targets failure boundaries: resolution thresholds, high aspect ratios,
self-intersection, concavity, reversed loops, large translations, periodic
surfaces, facet convergence, and geometry-preserving topology edits.

Two tools are deliberately separate:

- `accuracy.rs`: seeded, fixed-budget boundary corpus with exact analytic
  oracles and metamorphic conservation properties, returning measured errors.
- `query.rs`: a bounded, read-only, typed **GraphQL-inspired** query DSL over the
  shared `BrepModel`. It introspects surfaces, curves, loops, coedges, and mass
  independently of render tessellation.

Neither implementation uses proprietary Parasolid code. The inspiration for
query traversal is the public concept of an entity graph; the mathematical
oracles are conventional volume, surface area and inscribed-polygon formulae.

## Geometry Query Language (GQL subset)

Example:

```graphql
query ProbeCylinder {
  body {
    kind
    topology { faceCount edgeCount seamEdges eulerCharacteristic }
    mass { volume centroid { x y z } }
    faces(kind: "cylinder") {
      id
      surface { kind radius }
      loops { id coedges(limit: 8) { id edgeId twinId pcurve { kind } } }
    }
    edges(seam: true) { id closed coedges { id faceId reversed } }
  }
}
```

The only supported top-level field is `body`. Fields are **typed**, selections
are nested, and list fields accept a bounded `limit` and optional `id` filter.
Faces accept `kind: "plane" | "cylinder" | "sphere"`; edges accept
`curve: "line" | "circle"`, `seam: bool` and `closed: bool`.

This is intentionally **not spec-compliant GraphQL**. It does not accept
mutations, variables, aliases, fragments, subscriptions, schema introspection,
custom resolvers or arbitrary execution. Unsupported fields/arguments are
rejected with an offset and error code. Both the parser and resolver limit
byte length, syntax nesting, selection count, list length, visits and total
response bytes. Default list limit is 32, max 128.

Library:

```rust
let model = BrepModel::from_cylinder(&cylinder)?;
let value = geometry_query(&model,
    "{ body { faces(kind: \"cylinder\") { id surface { kind radius } } } }",
    QueryLimits::default())?;
```

CLI (`geometry.v1`):

```json
{
  "operation": "geometry_query",
  "source": {
    "kind": "analytic_cylinder", "origin": [0,0,0],
    "radius": 2, "height": 5
  },
  "query": "{ body { topology { faceCount seamEdges } mass { volume } } }"
}
```

Returns `{"schema_version":"geometry.v1","ok":true,"data":{"body":...}}`.
Existing `block`, `analytic_cylinder`, `edit_solid`, `inspect_brep` and
`command_history` outputs retain their previous fields.

## Deterministic boundary probe suite

CLI:

```json
{"operation":"probe_kernel","seed":13}
```

The `data` response includes counts, every named case and its serializable
reproducer, plus error metrics. The suite currently exercises blocks,
concave polygon prisms, cylinders with true circular carriers, and
geometry-preserving `split_edge`/`split_face` edits.

For valid models, measure:

- analytic volume vs independent volume formula (relative error);
- analytic area vs independent surface-area formula (relative error);
- centroid and bounding box vs independent expected values;
- consistency of `Curve3(t)` with `Surface3(Curve2(t))` across interior samples;
- signed mesh volume vs exact volume (**faceting**, not analytic-kernel error);
- the mesh-volume deficit vs an independently derived formula for an
  inscribed regular polygon: `1 − sin(2π/n)/(2π/n)`.

For expected-invalid models, explicit rejection is a pass. A case that
succeeds unexpectedly, fails validity unexpectedly, or exceeds an oracle
budget is a failed probe with a reproducer and seed. Seed variation affects
only bounded perturbations of selected border cases; the corpus is not an
uncontrolled random mesh soup or an exhaustive combinatorial search.

**Limitations:** this corpus is a starter fitness measure, not proof of
Parasolid-level robustness. Mass formula agreement is not evidence of correct
surface/surface intersections or Booleans. The next step is a separate
comparison runner against license-compatible OCCT outputs, with tolerances,
golden STEP fixtures, operation sequences and shrinking of failing examples.
The query layer remains read-only and operates on a validated snapshot, not a
persistent shared geometry database.
