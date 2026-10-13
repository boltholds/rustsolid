# Crossing through/blind Boolean: analytical boundary decomposition (next gated stage)

This document records the **read-only seven-region analytical decomposition** developed after v0.30. **RustSolid v0.31 now also implements a separate, restricted executable B-rep** with real coedges, UV trims, a closed genus-one shell and a certified derived display mesh. See [v0.31 mixed crossing B-rep](mixed-crossing-brep.md). The support-area planner itself stays read-only; its `can_apply()` remains false.

## Geometry

The source is a certified axis-aligned block. One +Y through-cylinder spans its full height, one top-entry blind-cylinder stops at a genuine interior floor. Both horizontal disks intersect at exactly two transverse points, checked using the v0.30 source-bound overlap planner.

The resulting physical boundary splits into seven analytical support regions:

| Region | Physical material boundary | Area formula |
|---|---|---|
| Top cap | box top minus union of both disks | A_top - (At + Ab - Ai) |
| Bottom cap | box bottom minus through disk | A_bottom - At |
| Outer sides | four original box walls | A_sides |
| Lower through wall | full through cylinder from box bottom to blind floor | 2π rt (H-d) |
| Upper through wall | only the arc of through circle **outside** blind footprint, from floor to top | rt * angle_exposed_through * d |
| Upper blind wall | arc of blind circle **outside** through footprint, from floor to top | rb * angle_exposed_blind * d |
| Blind floor | blind disk minus through disk overlap, at y=top-d | Ab - Ai |

This support-region area model is **not an arbitrary faceted approximation**. It is an analytical partition by actual plane and cylinder supports. The formula is independently testable without materializing a canonical topological graph. The existing v0.30 overlap volume remains:

`V = V0 - At * H - (Ab - Ai) * d`.

The v0.30 exposed-arc witnesses delimit a single two-arc top aperture. The separate lower through-circle is one physical opening. The independently expected closed-shell invariants are **genus 1** and Euler characteristic **0**; the v0.31 dedicated B-rep constructor now checks those values together with coedge/UV integrity.

## Subsequent generalized topology milestones

One restricted **one-through/one-blind transverse pair** now has source-bound carriers, valid coedges, a two-arc top aperture, bottom through aperture, blind-floor crescent, correct face senses and source-atomic transactions in a dedicated `MixedCrossingBrep`. The original `MixedThroughBlindBrep` supports multiple **disjoint** components and still deliberately rejects through/blind contacts.

The full v0.31 restricted constructor verifies a connected orientable shell, reciprocal twins, loop cycles, sampled UV→3D consistency, Euler 0/genus 1, independent source mass and watertight tessellation. Extending it to arbitrary through/blind networks, general curved Boolean or persistent topology naming remains **unsupported**. The generic `TopologyDeltaPlan::can_apply()` and `MixedCrossingSurfacePlan::can_apply()` stay false because their general execution paths have not been implemented.

This uses standard independent computational geometry only; no reverse-engineered proprietary code or ABI information enters RustSolid.
