# Crossing through/blind Boolean: analytical boundary decomposition (next gated stage)

The first stage after RustSolid v0.30 is a **real analytical decomposition of the boundary support regions**, but it deliberately does **not** yet allocate coedges, faces or commit a B-rep. Calling the present source-bound output a completed solid would be incorrect.

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

The v0.30 exposed-arc witnesses delimit a single two-arc top aperture. The separate lower through-circle is one physical opening. The expected resulting connected closed shell has **genus 1**, Euler characteristic **0**. These values are predictions to be independently checked when the true B-rep is assembled, not proof that a B-rep exists today.

## Next mandatory topology construction

Before enabling this case in `MixedThroughBlindBrep`, construct source-bound 3D circle/line carriers, per-face UV p-curves, shared physical intersection edges and coedges, the upper two-arc aperture, the bottom circle, and the blind floor crescent. The through-cylinder support must be split at the blind-floor height. The cylindrical side patches and floor must have exact material-facing orientation.

Then enforce a single connected orientable manifold shell, reciprocal twins, next/previous loop cycles, physically consistent UV→3D mappings, Euler 0/genus 1, independent analytical volume and surface-area checks, watertight tessellation and source-lineage-safe Undo/Redo. Until all gates pass, `MixedCrossingSurfacePlan::can_apply()` is always false.

This uses standard independent computational geometry only; no reverse-engineered proprietary code or ABI information enters RustSolid.
