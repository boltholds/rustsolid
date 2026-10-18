# RustSolid G7a — strict single-wall circular-notch Boolean

G7a expands the source-owned analytic curved Boolean support beyond
G4b2b3b's four independent corners and G6's two independent strips.
It constructs the exact **one connected closed shell** remaining when
an upright full-height cylinder intersects exactly the MIN-X wall of
a rectangular analytic box, with all other walls strictly outside the
source cylinder disk.

## New UV/topological arrangement

The exposed curved boundary is one **major circular arc spanning more
than pi**, not a small quarter arc and not faceted line segments. Its
endpoints lie on the same original MIN-X wall; the retained portions
of that original wall form two physically distinct analytic planar
faces. They preserve the source face ID in provenance but are
independent result faces with different physical edge loops.

Every actual cut edge has two independently source-certifiable UV
carriers, from the original box plane and cylinder surface. Every
original retained straight edge must be found in the G4b2b3a
source-subdivided BrepModel. Every new cut arc must match a globally
reconciled G1 endpoint pair and its original G4b2b2 source face-pair
carrier; no proximity welding or virtual periodic chart vertex.

Result topology (strict family):

- 1 closed genus-zero shell
- 12 topological vertices
- 18 physical edges and 36 opposite coedges
- 8 trimmed analytic faces (2 caps, 5 original planar side
  pieces, 1 inward reversed cylinder)
- Euler characteristic 2
- 8 strict signed source-membership witnesses.

The source-specific analytic mass uses a circular segment formula,
including the first moments and exact total boundary area. These
fields are **not an independent G4c-style volume oracle**. G7a sets
`closed_shell_proven=true` only for this exact family and always
`can_apply=false`. No document commit, persistent schema or G5c
application endpoint is enabled by G7a alone.

## Supported / Held

Supported: canonical AABB, upright +Y cylinder with the canonical
X/-Z chart, finite cutter entirely spanning the original box height,
cutter crossing only MIN-X face at two transverse interior points,
and guaranteed separation from the remaining five exterior walls.

Held: grazing/tangent contact, MAX-X instead of MIN-X, tilted or
eccentric chart, Z-wall contact, through-cylinder hole contained
entirely in the box, multi-cylinder intersections, surface overlaps,
unresolved G1/G4 cut carriers or uncertain material-side samples.

Tests:

```bash
cargo test --test analytic_one_wall_notch
cargo test --all-targets
```

Next: G7b independent source-volume bracket and B-rep oriented flux,
then a versioned G5 executor/persistent artifact before exposing this
family as an executable Boolean. Tilted cylinder and general UV
arrangement should follow once their independent curvature and region
classifiers prove closure and mass.

## G7b independent source mass/volume oracle

G7b uses two mathematically independent calculations rather than the
G7a source circular-segment antiderivative:

1. From the original cylinder radius and the wall clipping plane,
   construct both a polygonal chord approximation **inside** the
   retained right-hand disk segment and a tangent-intersection polygon
   **outside** it. Their actual shoelace areas bound source material
   volume from below and above after subtracting from the original box.
   Dyadic angular subdivision continues to the source-relative numeric
   budget. The G7a analytic mass formula is never used for this path.
2. Integrate the signed flux F=(x-cx,0,0), div F=1, across the actual
   B-rep outer X wall, both retained pieces of the X-min source wall
   and the reversed cylinder. Actual YZ polygon areas come from
   physical edge/coedge loops; the cylindrical contribution integrates
   -radius^2 cos²(theta) from the two true major-arc Edge3 carriers.

The independently bracketed source material, actual RESULT boundary
flux, candidate component mass and total retained+removed box volume
must agree. Boundary orientation and source revisions are revalidated
before issuing `volume_certified=true`.

As in G4c/G6b, the convex geometric enclosures are exact in real
arithmetic, with practical compensated floating point and explicit
f64 error margin; this is **not** directed-rounding interval proof
for arbitrary IEEE-754 NURBS. Numerically unresolved input is Held.

The G7b oracle itself is read-only and maintains `can_apply=false`.
A G5 document executor/persistent schema must be added and separately
verified before this new family becomes directly editable in Gefest.

Tests: `cargo test --test independent_one_wall_volume`.
