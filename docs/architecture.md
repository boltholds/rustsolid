# RustSolid architecture and next steps

## Intent

RustSolid is the exact/approximate geometric modeling subsystem for Gefest CAD. Keep the existing Gefest `solver-core` for constraint solving; the constraint solver determines sketch parameters, RustSolid constructs and queries solid geometry. Retain the OpenCascade reference path until independently implemented functionality meets parity benchmarks.

```text
Gefest sketch + constraint solver
           |
     stable geometry API
           |
  RustSolid modeling kernel
   |         |           |
 Topology  Geometry     Operators
 (B-rep)   (analytic)   (extrude, ...)
    \       /              |
   Topological naming + history
            |
   Tessellation / picking / formats
```

## v0.1 (implemented)

- Polygonal, one-shell orientable genus-zero boundary topology.
- Consistently oriented closed extrusion along +Y in XZ convention.
- Ear-clipping triangulation of simple polygonal caps without holes.
- Mesh-to-face ownership, bounds, analytic prism mass properties, ray picking.
- Validation and a versioned one-request JSON stdio bridge.

**Important:** `FaceId` and `EdgeId` are construction-local IDs, not persistent face/edge names across edits. A valid mesh is not an exact STEP model.

## Next vertical slices (planned, not implemented)

1. **Topology contract**: shell, coedge/half-edge, wire/loop, face orientation, manifold validator, stable typed handles, transaction semantics and diagnostic witnesses. Explicit support for multiple shells and internal voids comes later.
2. **Geometric layer**: parametrized planes, lines, circles, cylinders, then splines and NURBS with explicit domains, tolerances and surface/trim relationships.
3. **Modeling operations**: intersections, split, imprint and robust Boolean operations with independently testable predicates and topology provenance.
4. **CAD editing**: fillet/chamfer, offset, draft, direct face moves, deterministic selection and geometry-preserving history.
5. **Integration and exchange**: differential fixtures against OpenCascade, meshing tolerances, incremental regeneration and STEP support through license-compatible implementations.

## Invariants and tests

- No NaN/Inf, negative dimensions, zero-length boundary edges, self intersections or non-manifold edge incidences.
- Every mesh triangle has a valid owning topological face; every mesh edge is shared exactly twice in opposite orientations.
- The current single-shell topology has Euler characteristic `V - E + F = 2`.
- Analytic mass/bounds used as a baseline for polygonal prisms; clipping and picking must preserve face ownership.
- Add differential tests and metamorphic tests (rotation/reflection/scale/translation) before generalizing the primitive set.

## Intellectual-property boundary

Do not commit proprietary Siemens, Parasolid or decompiler artifacts to this repository. Treat published API names only as contextual feature taxonomy. Derive mathematical operators independently from permissively licensed literature and reproducible tests, documenting any third-party source and license.
