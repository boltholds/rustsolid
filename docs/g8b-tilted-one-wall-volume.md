# RustSolid G8b — source-only volume bounds for tilted one-wall Boolean

G8b adds an independent volume witness for G8a's bounded, strictly
transverse MIN-X notch family. It does not construct a result shell, use
the display mesh, integrate B-rep flux, or authorize a Boolean commit.

## Source-only horizontal ellipse

For the original cylinder frame origin O, orthogonal radial axes E1/E2,
axis A and radius R,

    S(u,v)=O+R E1 cos(u)+R E2 sin(u)+A v.

At horizontal Y=h its section is an affine ellipse:

    C(h)=O+A*(h-O_y)/A_y
    B1=R*(E1-A*E1_y/A_y)
    B2=R*(E2-A*E2_y/A_y)
    Q_h(u)=C(h)+B1*cos(u)+B2*sin(u).

The source-only Jacobian J=abs(B1.x*B2.z-B2.x*B1.z) equals R^2/A_y,
which is rechecked independently. The X amplitude
a=hypot(B1.x,B2.x) is height invariant. Thus the distance parameter

    d(h)=(C(h).x-xmin)/a

is affine in h. G8a already excludes tangencies, other walls and
cylinder cap contact, guaranteeing 0<d(h)<1 for the full box height.
The removed source footprint is the ellipse region to the right of the
original MIN-X wall. Its physical arc span is

    L(h)=2*acos(-d(h))  (pi<L<2*pi).

## Geometric 2D brackets

N equally spaced ellipse boundary points joined by true chords plus
the physical wall chord form an inscribed polygon. The polygon formed
by neighboring tangent intersections and the same end points
circumscribes the convex footprint. Their signed shoelace sums simplify
to

    I_N(h)=J/2*(N*sin(L/N)-sin(L))
    O_N(h)=I_N(h)+J*N*sin(L/(2N))^3/cos(L/(2N)).

Consequently I_N <= A_removed <= O_N in exact real arithmetic. This
geometric method does not read the G8a result candidate, source mass
antiderivative or mesh. The difference O_N-I_N is monotone in L,
so its largest value at the two box caps bounds the gap for all heights.

## Height integration and remainder

The true footprint area is

    A(h)=J*(pi/2+asin(d)+d*sqrt(1-d^2)).

With d affine, its fourth derivative is

    A''''(h)=-6*J*d*(d')^4/(1-d^2)^(5/2).

For d>0 its maximum absolute value M4 occurs at the largest endpoint
d. With H=Ymax-Ymin, even m panels and step=H/m, composite Simpson
has a rigorous real-arithmetic integral remainder

    |int A dh - Simpson(A)| <= H*step^4*M4/180.

Since Simpson weights are all positive, we may replace each A node
with I_N and O_N and retain a valid enclosure:

    lower_removed=Simpson(I_N)-remainder-f64_guard
    upper_removed=Simpson(O_N)+remainder+f64_guard
    lower_retained=box_volume-upper_removed
    upper_retained=box_volume-lower_removed.

Both panel counts are adaptively doubled under explicit caps.
Unresolved geometry and numerical budgets are Held, never rounded
silently into a confirmed edit.

The **geometric polygon bounds and analytic Simpson remainder** are
rigorous in exact real arithmetic. The separate practical f64 roundoff
guard is explicitly recorded and is not a formal outward-rounded
interval arithmetic guarantee for arbitrary CAD curves.

## Second source-only mass check

A separate integral uses the analytic area primitive

    P(d)=pi*d/2 + d*asin(d)
       +sqrt(1-d^2)-(1-d^2)^(3/2)/3.

For a near-constant d, a high-order Simpson check of the same analytic
area avoids subtractive cancellation. This independently calculated
source volume must agree with the polygon/height bracket within
declared margins. It does not replace a boundary flux integral through
a real resulting B-rep. That remains the job of G8c.

## Revision, mutation and testing

The oracle revalidates the complete G8a source, checks original box
mass, units and finite-cylinder volume, and reconstructs every saved
source-derived numeric witness. The cylinder source is checked by full
geometry comparison because it has no standalone revision.

The only positive gate is source_bracket_established. Both
result_shell_certified and can_apply remain false. There is no
new document, archive schema, API or executable general Boolean.

Tests:

    cargo test --test independent_tilted_volume
    cargo test --all-targets

G8b tests include positive and negative X tilt, YZ-only tilt, an
explicit independent vertex shoelace calculation, translated/offcenter
geometry, source revisions, witness tampering, numerical refusal, and
compatibility with the existing upright G7 operation.

G8c will reconstruct a true analytic closed result B-rep using the
four physical G1 cut arcs, certify actual trimmed 3D/UV coedges and
shell orientation, and independently integrate the real shell flux
against the G8b source-only retained-volume interval. Only then may
a bounded revision-atomic Boolean executor be introduced.
