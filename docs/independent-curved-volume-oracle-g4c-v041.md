# RustSolid G4c1 — independent four-shell volume evidence

Status: bounded-family volume certification of G4b2b3b, NOT an
arbitrary multi-surface or NURBS Boolean oracle. G5 remains disabled.

Supported operation: Difference of a canonical AABB and full-height
upright analytic cylinder whose transverse wall cuts leave four
disconnected corner solids. A fresh G4b2b3b candidate is mandatory.

## Two independent mathematical paths

1. Original-source material membership integral. Compute for each
   corner its local positive wall distances a,b from the cylinder
   center, radius r, and box height h. Integrate the surviving strip
   width w(u)=b-sqrt(r^2-u^2) from sqrt(r^2-b^2) to a.
   Since w''(u)>0 over the strictly supported domain, midpoint sums
   bound the area from below and trapezoid sums from above.
   Adaptive dyadic subdivision stops only when a verified gap fits a
   source-scale volume tolerance. Conservative f64 error slack is
   explicit. Ill-conditioned cases or a panel budget overrun are Held.
   Neither candidate exact_volume nor model.mass is read for this path.

2. Result B-rep surface-flux integral. Independently traverse every
   actual shell, oriented face loop, original x-wall polygon and
   inward-facing analytic cylindrical patch. Integrate the vector
   field F=(x-cx,0,0), divergence(F)=1, over the retained closed shell.
   The x-wall flux comes from its real 3D projected YZ polygon area.
   The cylinder contribution comes from two physical circular edge
   parameter domains and its inward face normal, giving
   -r^2*h*integral(cos^2(theta)dtheta). The antiderivative here is
   theta/2 + sin(2theta)/4, mathematically distinct from the source
   strip integration and from G4b2b3b's source mass formula.
   All other planar faces must have zero x-normal. Two real circular
   arcs and two generators must bound each cylindrical patch.

Each component's signed flux must agree with its source interval,
and its declared G4b2b3b component volume must agree with BOTH paths.
The total source interval must contain the B-rep model's stored mass
within policy. The box's independently computed volume is partitioned
into retained and removed-in-box volumes; removed volume must be less
than the full finite source cutter volume. The latter is not
subtracted directly because its caps extend beyond the original box.

## Numeric interpretation

The convexity brackets are mathematical enclosures in exact real
arithmetic. Their floating-point implementation uses explicit
conservative error margins and compensated summation, but is not a
formally directed-rounding interval arithmetic system. Do not
interpret it as a proof of arbitrary IEEE-754 rounding error bounds.
The acceptance error derives from source feature extent and the
declared model GeometryTolerance, never from a claimed result volume.
Numerically unresolved cases are Held without a partial result.

## API and gates

IndependentVolumeOraclePlan::build(first,second,&candidate,tolerance)
returns either Ready or Held. The candidate itself is fully
source-revalidated. Result includes a SourceVolumeBracket and a
ShellFluxWitness per shell, the total source interval, the total
B-rep flux, the full original box volume, and the intersection volume
bounds inside the box. validate_against recomputes source and boundary
evidence and rejects stale revisions, changed tolerance or witnesses.

volume_certified=true applies only to this strict case. can_apply=false
always. No document mutation, Undo/Redo or G5 transaction is implied.

## Tests

cargo test --test independent_curved_volume_oracle
cargo test --all-targets

Acceptance includes all four shell volumes, two independent methods,
rectangle with translated asymmetric cylinder, f64 error control,
source mass conservation, tamper rejection and revision/tolerance
staleness. G4b2b3b remains independently regression tested.
