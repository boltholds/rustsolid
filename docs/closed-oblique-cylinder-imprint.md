# RustSolid v0.37 — closed oblique cylinder imprint and owned coedge carriers

This is a **bounded topological editing milestone**, based entirely on
independent, first-party analytic geometry. The operator imprints one
**closed, transverse, oblique plane/cylinder section** onto the lateral
surface of a source-certified finite right cylinder. The source material
is unchanged. It does **not** execute general Boolean union, subtraction,
intersection, multi-face trimming, or arbitrary NURBS operations.

## Why this executable slice is mathematically tractable

A transverse plane cutting the infinite support of a right circular
cylinder has a one-period harmonic parameterization in the cylinder
angular coordinate:

```text
S(u,v) = O + r X cos(u) + r Y sin(u) + A v
v(u) = m + c cos(u) + s sin(u)
C(u) = (O + A m) + (r X + A c) cos(u) + (r Y + A s) sin(u)
```

The executor requires the **entire** intersection curve to lie
strictly between the real lower and upper axial bounds of the finite
cylindrical side, with an explicit tolerance margin:

```text
epsilon < m - hypot(c,s) <= m + hypot(c,s) < height - epsilon
```

A near-cap, tangent, split, open or numerically unresolved crossing is
`Held`, never silently rounded to an executable imprint.

## Persistable analytic carriers

The generic catalogs now represent three related carrier forms:

- `Curve3::Harmonic(HarmonicCurve3)`: exact-in-form 3D ellipse
  with center, cosine and sine vectors, parameterized by original
  cylinder angle. The basis vectors are generally **not orthogonal**.
- `Curve2::Harmonic(HarmonicCurve2)`: the same analytic curve
  expressed in the UV coordinates of the *intersecting plane*.
- `Curve2::CylinderSection(CylinderSection2)`: unwrapped `u`
  and harmonic axial `v(u)`; it supports reversed coedge traversal
  without changing the physical ellipse or introducing a fake seam edge.

`IntersectionPcurve::to_curve2()` maps the v0.36 preview into a
self-contained `Curve2` with no proprietary memory/layout dependencies.
Both preview pcurves reconstruct the same 3D arc. A generated plane
pcurve remains a witness of the cutter plane, **not a coedge of the
source-cylinder result**.

## One actual physical B-rep cut edge and two coedges

`ClosedCylinderImprint::imprint` accepts a source-bound
`TopologyDeltaPlan` plus its verified `IntersectionGraph`,
original plane face and original cylinder B-rep. These are *intent
inputs*, not general authorization to `can_apply=true`.

It reconstructs the source cylinder from `CylindricalBrep` and
compares it against the actual graph operand before any write.
Only a single closed non-tangent periodic arc with one virtual graph
anchor and the original side face `FaceId(2)` can proceed.
The source cannot already have been imprinted.

The new B-rep has:

| Entity | Before | After |
| --- | ---: | ---: |
| Vertices | 2 | 3 |
| Physical edges | 3 | 5 |
| Coedges | 6 | 10 |
| Loops | 3 | 4 |
| Faces | 3 | 4 |
| Closed shells | 1 | 1 |
| Euler characteristic | 2 | 2 |

The original seam physical line is split at `v(0)`, so the chosen
ellipse anchor has a real topological vertex and both new side
patches share the same physical seam junction. The new ellipse itself
is one true `EdgeId(4)`, supported by the harmonic spatial carrier,
with exactly two opposite uses `CoedgeId(8)` and `CoedgeId(9)`
in two distinct cylinder-side faces. Both coedges own an analytic
`Curve2::CylinderSection`; one evaluates `u`, the other
`2π−u`. The original bottom and top cap face/edge carriers remain
the same.

The canonical `BrepModel::validate()` checks each UV trim against
its physical 3D edge carrier and face surface at interior samples
and endpoints, twin reciprocity, closed directed loops, seam UV
distinction, orientability, and connectedness. Independently checked
analytic source mass and bounds must be identical before and after.
This is a numerical witness and finite-domain structural validation,
not a formal continuous robustness theorem.

## Transaction semantics

The editor creates a **private candidate** before publishing any
topology. A failed validator never modifies `current`.
On success the model swaps atomically, the revision increases once,
and one bounded (128) undo snapshot is retained. Undo/Redo each restore
a validated complete B-rep snapshot and increase the revision
monotonically. After Undo, a new successful imprint discards the redo
branch. A second imprint without Undo returns
`ImprintHold::AlreadyImprinted`.

For rejected geometry a typed `ImprintResult::Held(reason)` is
returned. Forged or mutated graph/delta source witnesses fail source
validation rather than being accepted. The public
`TopologyDeltaPlan::can_apply` remains **false**, and the general
`BrepArena::require_executable_plan` guard is not relaxed.

## Scope restrictions / next gate

1. The imprinted B-rep's analytic mass and envelope are inherited from
   the **unchanged cylinder**. A general Boolean **cannot** reuse
   these mass values for a shape-changing edit without independent
   source integration.
2. The stored cutter-plane pcurve is independent preview evidence.
   Cutter-plane faces are not part of the imprinted cylinder shell.
   A real body/body Boolean must modify both operands' topological
   material boundaries and perform cell classification first.
3. Existing source display mesh remains geometrically watertight, but
   a new **face-owned** mesher is needed to render the split cylinder
   wall with distinct upper and lower face IDs and shared new ellipse
   vertices. Mesh face ownership is **not** certified for the imprint.
4. Only a pristine canonical analytic cylinder, one transverse closed
   branch and a finite plane face are supported. No arbitrary oblique
   *open* arc terminating on an original cap rim, multi-branch
   interpolation, simultaneous contact, void creation or persistent
   cross-session topological names are claimed.
5. `Curve3::Harmonic` boundary splitting and inverse parameter
   classification are explicitly unsupported in generic planners until
   a separate robust numerical inverse is provided. They fail closed.
6. The editor compares all original canonical source fields
   structurally before imprinting; it is not a general heterogeneous
   mutable BrepModel transaction framework.
7. The existing rectangle-only `TrimmedIntersectionEngine` now checks
   the **actual four straight UV sides**, not only the number of
   coedges and parameter bounds. It explicitly rejects a newly
   imprinted harmonic side domain until a generalized trim-domain
   classifier is implemented. Otherwise, a later query could silently
   mistake the elliptically trimmed face for a rectangular one.

Next engineering gates:
- source-owned face meshing of the two trimmed cylinder side patches;
- internal/outer loop classification for cutter-plane regions;
- generalized `Curve3` inverse coordinate and physical edge subdivision;
- apply two-sided region classification and topology delta to a
  bounded **volume-changing** Boolean, independently integrating mass;
- keep rollback/Undo/Redo and all near-degenerate errors fail-closed.

## Verification

`tests/closed_cylinder_imprint.rs` independently checks:
closed intersection graph provenance, both analytic pcurve carriers,
physical edge/coedge/face topology, source mass, UV-to-3D agreement,
Undo/Redo invariance and monotonic revision, forged plan rejection,
unsupported source refusal, and near-cap holds.

```sh
cargo test --test closed_cylinder_imprint
cargo test --test oblique_plane_cylinder_preview
cargo test --all-targets
cargo check --release --all-targets
```

Verification on 2026-10-14 in an isolated WSL copy of the RustSolid
working tree (source based on parent `06ea004`):

- `cargo test --test closed_cylinder_imprint` — **4 passed, 0 failed**.
- `cargo test --test oblique_plane_cylinder_preview` — **5 passed, 0 failed**
  in the preceding v0.36 preview suite.
- `cargo test --all-targets` — **passed**, including the new imprint and
  original v0.35 regression suites.
- `cargo check --release --all-targets` — **passed**.
- The new regression suite also rejects B-rep trim tampering and refuses
  to classify a harmonic cylinder-side trim as a rectangle.

Pre-existing compiler warnings remain. No production CAD parity is
asserted, and no proprietary DLL was executed.

No proprietary Siemens/Parasolid code, decompiler output or private
binary layouts are used by this implementation.

## English summary

RustSolid now has a bounded executable *imprint* of a closed oblique
plane/cylinder intersection on an analytic cylinder: one true spatial
elliptical edge, two opposite coedges with analytic cylinder UV
pcurves, a consistent split of the cylindrical wall, structural
validation and atomic Undo/Redo. It does not execute a general Boolean
or claim face-owned display fidelity for the newly split patches.
