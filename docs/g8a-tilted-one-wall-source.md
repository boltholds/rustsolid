# G8a — conservative tilted one-wall source geometry (read-only)

G7 commits an upright +Y cylinder with two planar MIN-X fragments and a
single inward cylindrical patch. G8a introduces a **tilted** finite
right circular cutter while remaining strictly read-only. There is no
result B-rep, face patch assembly, independent volume oracle, edit
transaction, archive schema or geometry.v1 endpoint in this change.

## Exact-in-form source section

For a source cylinder with frame origin O, radial basis E1/E2,
unit axis A and radius R,

    S(u,v) = O + R(E1 cos(u) + E2 sin(u)) + A v

the intersection of the **side surface** with a horizontal box plane
Y=h has an exact harmonic carrier

    v(u) = (h - O_y - R E1_y cos(u) - R E2_y sin(u)) / A_y
    Q_h(u) = C_h + B1 cos(u) + B2 sin(u)
    C_h = O + A * (h-O_y)/A_y
    B1 = R * (E1 - A * E1_y/A_y)
    B2 = R * (E2 - A * E2_y/A_y).

The Y components of B1 and B2 are zero. The section is an ellipse in
world XZ, evaluated directly from original input geometry. No
tessellation, root search or fitted curve becomes a topology carrier.

On MIN-X at X=xmin, the two transverse analytic events satisfy

    A_x cos(u-phi) = xmin - C_h.x
    A_x = hypot(B1.x, B2.x), phi = atan2(B2.x, B1.x)
    u_low/high = phi -/+ acos((xmin-C_h.x)/A_x).

G8a requires the center to be strictly right of xmin and the ellipse to
cross xmin with positive margin. Hence the retained major arc has
span strictly between pi and 2*pi. Both cap events are source-backed.

## Conservative whole-height family gate

- Original left operand: an independently certified canonical AABB.
- Right operand: a valid analytic finite cylinder with unit axis
  A_y greater than 0.90 and a resolvable nonzero XZ tilt.
- All eight box corners, projected onto the cylinder axis, lie strictly
  **between** the two cutter cap planes. By affine extremality this
  excludes both cutter end caps from the entire box.
- Horizontal ellipse X and Z coordinate envelopes are affine in h.
  Checking their strict margins at both Y levels proves that every
  intermediate section intersects **only** MIN-X, never MAX-X/Z walls.
- Tangencies, coincident contacts, uncertain numerical margins, other
  operations, extraneous source arcs and any G1 ambiguity are Held.
- One transverse analytic ellipse-major arc on each Y cap and two
  independent X-wall intersection arcs are required.
- Each of the four **physical** G1 endpoints is source-topology
  reconciled, with degree 2 in one connected four-arc cycle. Periodic
  chart anchors do not count as physical endpoints.
- All 3D and source-face UV analytic carriers are checked against both
  original surfaces. X-wall branches may be harmonic curves or true
  generators if the tilt lies in the YZ plane.

The carrier registry is the existing G1/G4b2b2 source-certified
GlobalAnalyticCutPlan. validate_against rechecks the original source
revisions, the cut registry, cross-sections and physical cycle. The
return explicitly sets can_apply=false, closed_shell_proven=false.

Because source certification does not imply result topology, G8a makes
**no claims** about coedges, shell Euler characteristic, solid volume,
permitted document commits, G7 schema-v3 compatibility or generic CAD
Boolean. G7's published upright contracts remain unchanged.

## Source tests

Run:

    cargo test --test g8_tilted_one_wall_source
    cargo test --all-targets

Tests cover X-tilt, opposite X-tilt, translated/off-center sources,
YZ-axis tilt (line generators on MIN-X), implicit cylinder equation
independent of G1, reconstructed dual-face UV/3D curves, a four-event
cycle, revision/tamper rejection, cylinder cap intrusion, other-wall
contact, an unsupported operation and the existing upright G7 path.

## Next

G8b should construct source-only swept volume *brackets* for the
Y-dependent clipped elliptical footprint, independently cross-check
the signed planar cap integrals, and reject unresolved numerical
cases. G8c should assemble a closed 8-face B-rep from four actual
G1 arcs and G8b-certified volume, preserving generational history
and archive compatibility only after the new result is independently
validated.
