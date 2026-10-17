# G6a: two-wall analytic Boolean (strict vertical cylinder)

A bounded new certified family beyond the existing G4b2b3b four-corner
fixture. The world-Y upright circular cutter spans the full box height.
Its disk intersects only the two opposite X walls transversely and lies
strictly inside both Z walls. The remaining material has **two disconnected
closed shells**, not four.

## Geometry and invariants

G1/G2/G4b2b3a are rebuilt before assembling two retained curved-UV
prismatic face patches. Source line edges on box caps and side walls are
validated against exact G4b2b3a source edge splits. Every new physical
cut edge is matched by original face-pair ID, analytic 3D curve samples,
and real G1 endpoint identities to one unique G4b2b2 source arc.
Both physical coedges have independently evaluable source-surface
Curve2 carriers; periodic virtual chart anchors are never welded.

Each component has 8 vertices, 12 physical edges, 24 opposite coedges,
6 bounded faces (2 horizontal planar caps, 3 planar side walls and one
reversed cylinder wall) and one closed genus-zero shell. Totals: two
shells, V=16, E=24, F=12, loops=12, coedges=48, total Euler=4,
total genus=0. Material orientation is checked by strict positive and
negative side samples of all 12 faces, and radial vertex fans have
independent manifold certification.

The **source-specific analytic area and mass are candidate metadata**,
not independent G4c proof. General NURBS curved UV arrangement, tilted
cylinders, boundary tangencies, multiple cutters and intersection/union
remain deliberately held. can_apply=false.

Tests: cargo test --test analytic_two_wall_partition.
## G6b independent volume certification

AnalyticTwoWallVolumeOracle verifies both components without consulting
the candidate mass formula for either independent computation:

- Source-space integration of the surviving material strip
  width(z)=b-sqrt(radius^2-x^2) in a rectangular cross-section.
  Since the integrand is convex, midpoint integration is a lower
  estimate and trapezoid integration is an upper estimate with adaptive
  subdivision and an explicit floating-point budget.
- Oriented RESULT B-rep boundary-flux integration using F=(0,0,z-cz)
  and div F=1. The source-owned outer Z wall contributes its physical
  XY polygon area; the reversed cylinder contributes the analytic
  integral of -radius^2*sin(theta)^2 over each true circular trim.
- The per-shell bracket, B-rep flux, candidate analytic mass and total
  source-box balance must agree. Unsupported/tangent/numerically
  unresolved shapes remain Held. G6b is not directed-rounding interval
  arithmetic for general NURBS surfaces.

G6a/B build previews both keep can_apply=false. A separately audited
G5 operation and versioned persistent document schema are required
before new Boolean material can be committed.

Tests: cargo test --test independent_two_wall_volume.

