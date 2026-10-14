# RustSolid v0.38 — certified imprint mesh, analytic face regions and bounded volume-changing cut

This increment builds on v0.37's closed, transverse plane/cylinder imprint
and v0.36's paired analytic UV witnesses. It is written independently
from public mathematics and the pre-existing RustSolid source tree. It
does **not** use proprietary Parasolid/Siemens code or implement general
body/body Boolean operations.

## 1. Four-face, shared-edge imprint mesh

The original cylinder is an analytic solid with a bottom circle, upper circle
and a single lateral surface. v0.37 splits the lateral face into two using
one physical spatial ellipse and two opposite coedges with cylindrical UV
carriers. v0.38 adds `ClosedCylinderImprint::tessellate(n)`:

- one sample ring on the bottom circular edge;
- one sample ring on the single physical elliptical cut edge;
- one sample ring on the upper circular edge;
- two cap centers (not mistaken for physical B-rep vertices);
- **no duplicate UV seam vertices** on the 3D mesh.

The two triangles of each lower-wall segment are owned by `FaceId(2)`;
the two of each upper-wall segment by `FaceId(3)`. Cap triangles use
`FaceId(0)` and `FaceId(1)`. For `n` angular segments:

```text
vertices = 3n + 2
triangles = 6n
triangle owners = [n, n, 2n, 2n]
χ(mesh) = V - E + F = 2
```

`ImprintMeshCertificate` is checked against the actual vertices, physical
edge carriers, two opposite coedge UV charts, per-face ownership, exterior
material normals, undirected manifold incidence, signed volume and Euler
characteristic. The triangle mesh is never treated as the source's exact mass.

For radius `r`, height `h`, and uniform angular resolution `n`, an
independent signed-volume oracle is the **inscribed polygonal prism**:

```text
A_polygon = n r² sin(2π/n)/2
V_polygon = A_polygon h
```

The calculated facet volume must match this oracle within a stated finite
floating-point budget. It is below the analytical cylinder volume except
for numerical rounding. All four face areas must be positive.

The certificate is revision-bound: after Undo/Redo, an old mesh is not
silently reused.

## 2. Analytic whole-face region classification

`classify_imprinted_cylinder` and `ObliqueHalfspaceEditor::classify`
recognize exactly one closed transverse section with a nonzero axial
derivative. Write the cylinder as

```text
S(θ,v) = O + r X cosθ + r Y sinθ + A v
v_cut(θ) = m + c cosθ + s sinθ
```

The validated finite separation condition is

```text
0 < m - sqrt(c²+s²)
m + sqrt(c²+s²) < h
```

with a positive numerical margin. Since the cutting plane equation is
affine along every cylinder generator, the entire lower cylindrical
face (`0 ≤ v ≤ v_cut`) stays in the **negative** halfspace. The entire
upper face (`v_cut ≤ v ≤ h`) stays in the **positive** halfspace.
The original bottom and top cap faces have the matching classifications.

This is a mathematical whole-cell classification for one supported
analytic family. It is **not** a general face arrangement, polygon/NURBS
region classifier or region result from an arbitrary multi-body Boolean.

## 3. Bounded, volume-changing Boolean with an infinite plane

`ObliqueHalfspaceEditor` constructs a new, closed, source-certified
analytic solid after choosing `ObliqueHalfspaceSide::Negative` or
`Positive`. The two outputs are the below-plane and above-plane
portions of the original finite cylinder. The input is v0.37's closed
imprint; the halfspace is unbounded, so this operation does not yet
support Boolean between two arbitrary finite solids.

Each retained output has one original circular cap, one **new elliptical
planar cap** in the actual cutting plane, and one cylindrical side face.
The new curved cap and cylinder wall share a single true 3D elliptical
edge with opposite coedges and distinct, source-bound 2D pcurves.
The original unretained cylinder cap and half wall are removed from the
candidate. The result has:

```text
vertices = 2, edges = 3, coedges = 6, loops = 3, faces = 3
closed shells = 1, Euler characteristic = 2, orientable
```

`BrepOrigin::AnalyticObliqueClip` distinguishes the result from the
uncut `AnalyticCylinder` source. New planar support frames are built
with the *correct outward normal* for the retained halfspace and
validated against their coedge UV curves.

### Analytic mass and centroid

Let `R` be the cylinder radius, `h` its height, and `m,c,s` the
axial harmonic coefficients. Define `d²=c²+s²` and
`A=πR²`.

```text
V_negative = A m
V_positive = A (h-m)
```

The *local* centroids in the cylinder frame `(X,Y,A_axis)` are:

```text
lower radial = (R c / (4m), R s / (4m))
lower axial  = m/2 + d²/(8m)

upper radial = (-R c / (4(h-m)), -R s / (4(h-m)))
upper axial  = (h+m)/2 - d²/(8(h-m))
```

The new elliptic cap area is `π |cosine × sine|`. The cylinder-wall
area is `2πR m` or `2πR(h-m)`. Bounds are derived from exact maxima
of the endpoint circular/elliptical harmonic rings in each Cartesian
axis, not from the original whole-cylinder bounding box.

Analytical volume and first moments are independently tested by polar
quadrature of the affine height over the actual disk, not by calling the
production formula a second time. Lower and upper volumes sum to
the source volume, and their weighted centroid equals the full
cylinder centroid.

### Atomicity

The editor constructs a private B-rep candidate; the full
`BrepModel::validate()` checks incidence, coedge/pcurve consistency,
material-facing orientation and closed Euler topology. Before any
candidate becomes observable, a **separate 64-segment face-owned
watertight mesh** is generated and independently certified against
the polygonal prism oracle, outward face normals and source rings.

Only then does the editor swap the candidate, update the revision
exactly once and record a bounded inverse snapshot. Undo and Redo
restore validated complete models and increment revisions
monotonically. Invalid inputs/failed candidate validation never
partially mutate the original imprinted model.

## 4. Three-face certified cut display

`ObliqueHalfspaceEditor::tessellate(n)` returns one canonical rim and
one elliptical rim with cap centers; its `CertifiedObliqueClipDisplay`
includes a provenance-bound `ObliqueClipMeshCertificate`. For `n`
angular segments:

```text
vertices = 2n + 2
triangles = 4n
triangle owners = [n, n, 2n]
χ(mesh) = 2
V_polygon_negative = A_polygon m
V_polygon_positive = A_polygon (h-m)
```

The two clipped meshes sum to the polygonal prism volume with the
original cylinder's height. Every elliptical mesh edge is used by the
cylindrical side and the new planar cap in opposite directions.
The mesh validator rejects incorrect source-ring coordinates, triangle
owners, normal reversal, open incidence, altered volume certificates
and stale revisions after Undo/Redo.

## 5. Supported and unsupported modeling operations

Supported:

- source-certified **one** closed non-tangent, non-coincident oblique
  plane/cylinder imprint;
- analytic two-region classification for that fixed topology;
- retain precisely the lower or upper connected portion of that
  *one* source cylinder;
- new closed B-rep, analytic ellipse, mass/centroid/bounds,
  watertight face-owned display and atomic Undo/Redo;
- failure/rollback without exposing partial material.

Not supported:

- generic `TopologyDeltaPlan::can_apply` is **still false**;
- arbitrary two-body Boolean, open traces, interacting cylinders,
  NURBS, general inner loops, multi-shell or nonmanifold topology;
- ellipsoid/ellipse boundary inverse-parameter classification in the
  general `IntersectionGraphBuilder`;
- coincidence, tangency and numerically uncertain cuts outside the
  narrow already-accepted source predicate;
- durable feature naming or automatic propagation into Gefest's
  general feature/history graph;
- a CLI/JSON command for this new bounded editing operation.

Existing v0.35 connected three-cylinder Boolean remains a separate
source-certified feature. Neither its assumptions nor this v0.38
halfspace adapter are silently generalized.

## Verification commands

```sh
cargo test --test imprint_mesh
cargo test --test oblique_halfspace
cargo test --test oblique_clip_mesh
cargo test --all-targets
cargo check --release --all-targets
```

Independent numerical witnesses validate tested scales and families,
not all floating-point inputs. No formal interval proof or
production-grade machining certificate is claimed.

## Next gate

1. Add a source-owned `BrepBody`/Gefest adapter and a versioned JSON
   request for the bounded halfspace Boolean, preserving provenance
   and canonical state.
2. Add generalized classification of multiple finite face cells,
   including inner loops and multiple intersection arcs.
3. Add source-independent face-pair intersection families and numerical
   degeneracy processors; keep tangency/near-coincidence Held until proven.
4. Extend `TopologyDeltaPlan` to a general mutable closed-shell executor
   only after mixed-source mass, watertightness, source revalidation and
   rollback contracts cover it.

## English summary

RustSolid v0.38 introduces a four-face certified mesh for the existing
oblique cylinder imprint, whole-face analytical halfspace classification,
and two bounded volume-changing cuts (retain the material above or
below the transverse plane). Both new solids are closed analytical
B-reps with an elliptical cap, exact-in-form mass and centroid and
certified three-face meshes. Undo/Redo is atomic. General Boolean
between arbitrary finite B-rep solids remains unsupported.
