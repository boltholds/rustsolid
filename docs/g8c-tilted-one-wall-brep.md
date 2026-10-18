# RustSolid G8c — source-certified tilted MIN-X analytic B-rep and actual-shell flux

G8c stacks on the G8a physical source-event graph and G8b independent
source-volume intervals. It constructs a *real closed candidate* B-rep
for a tilted cylinder intersecting ONLY the MIN-X wall of an AABB.

This is still a **read-only staged** geometric result. Undo/Redo,
generational handles, producer DAG, archive format, Gefest integration
and any executable/atomic document Boolean remain out of scope.

## Physical topology

Exactly one orientable, closed, genus-zero shell:

| Property | Expected |
|---|---:|
| Shells | 1 |
| Vertices | 12 |
| Edges | 18 |
| Coedges | 36 |
| Faces | 8 |
| Loops | 8 |
| Euler characteristic | 2 |
| Genus | 0 |

Two physical retained MIN-X planar faces maintain the original source
face provenance as two DISTINCT result faces. The cutter's one inward
cylindrical face is a true trimmed analytic cylinder support; no mesh
facets form the topology.

Four exact source-owned cut carriers are reused as B-rep edges:

1. Bottom Y-cap major harmonic ellipse arc.
2. Top Y-cap major harmonic ellipse arc.
3. Upper MIN-X/cylinder line or harmonic arc.
4. Lower MIN-X/cylinder line or harmonic arc.

Each physical result edge maps to exactly one original G1
GlobalArcRef, and both endpoints refer to its original reconciled
GlobalEventIds. Original unaffected AABB edges are checked against
their source physical edges and original incident source-face IDs.

Each cut edge's two analytic Curve2 coedges come from the matching
original source surfaces. A reversed coedge is reparameterized
analytically, not sampled. The result cylinder face has one connected
unwrapped UV loop: source-chart shifts by integral multiples of 2pi
reconcile the periodic angular coordinates without creating synthetic
physical vertices or losing source identity. UV bounds include
harmonic-curve extrema, not only their physical endpoints.

Independent radial audit and the canonical shared BrepModel validator
must accept physical opposite twins, shell connectedness and
face/material orientation. The result cylinder is
FaceSense::Reversed; original planar result faces remain Forward.
The candidate is immutable outside the builder and cannot be applied.

## Independent result-shell volume

G8c computes an ACTUAL RESULT boundary-flux witness, not a second
source-only analytic notch formula.

Choose F=(x-xmin,0,0), div(F)=1. The exterior Xmax planar face
contributes its actual oriented coedge-loop YZ area times
(xmax-xmin). All other original planar faces carry zero X flux
(including the two MIN-X fragments because x=xmin).

The result cylindrical face contributes an inward flux. With the
actual trimmed result cylinder support

    S(u,v)=O+R*(E1*cos(u)+E2*sin(u))+A*v,

and n_x(u)=E1.x*cos(u)+E2.x*sin(u), the outward cutter normal
is n. For the inward material face the flux density is

    -R*n_x(u)*(O.x-xmin +R*n_x(u)+A.x*v).

Integrating the positive cutter-normal density first in v yields

    G(u,v)=R*n_x(u)*[(O.x-xmin+R*n_x(u))*v + A.x*v^2/2].

Green's theorem and the RESULT face coedge orientation reduce the
cylinder face's inward flux to

    Phi_cylinder=-integral_(actual cylinder UV coedge loop) G(u,v) du.

No original-source section area, G8b analytic antiderivative,
faceted mesh or candidate-mass field is used in this integral.

For every actual CylinderSection2 coedge, v(u)=a+b*cos(u)+c*sin(u).
The integrand G(u,v(u)) is a trigonometric polynomial of order at most
three, integrated analytically via a bounded -3..3 Fourier series.
True cylinder generator coedges have constant angular u and
contribute zero to the boundary integral. All four coedges must be
verified, and their UV endpoint chart must close consistently.

Total flux is exterior Xmax flux plus inward cylindrical flux.
It is independently required to agree with G8b's source-only retained
volume interval within explicit f64 error slack. It must also agree
with the STAGED result BrepModel.volume, which is the G8b interval
midpoint. Negative/wrong cylinder orientation and invalid loop
topology fail closed.

The trigonometric integration is analytic in real arithmetic.
Practical f64 floating-point roundoff is bounded by an engineering
guard, not formally directed IEEE interval arithmetic.

## Important mass-property boundary

BrepModel requires a positive MassProperties structure. G8c stages
the G8b interval midpoint as volume. **Surface area and centroid are
placeholders from the input box**, pending a source-independent
surface/first-moment integration before any executable CAD Boolean.
They are not certified result properties and must not be used
downstream as production mass evidence.

G8c shell and flux witnesses explicitly expose
`can_apply=false`. Shell closure is certified only for the
specific strict G8a/G8b family; this does NOT imply arbitrary tilted
NURBS, multiwall, tangency, cylinder-cylinder or multicut Boolean.

## Verification

    cargo test --test g8c_tilted_one_wall_brep
    cargo test --all-targets

Tests cover opposite signed X tilts, YZ-only tilt with true source
generator lines, offcenter/translated sources, every original G1
cut arc and dual face UV/3D carrier, physical UV closure, corrupt
face sense/edge twin/trim/source metadata and upstream G7 regression.

Next: G8d revision-atomic G5-style document and versioned,
source-reconstructed schema for this new bounded family, **after**
certifying actual surface area and centroid and testing archive
reopening/Undo/Redo against independently rebuilt source topology.
