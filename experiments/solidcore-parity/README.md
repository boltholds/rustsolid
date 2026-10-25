# SolidCore / NVIDIA SMLib parity experiment

**Status:** isolated, experimental Rust adapter for comparing primitive construction and
independent B-rep measurements. This does **not** make SolidCore a RustSolid Boolean
backend and does not change Gefest's production geometry path.

This experiment builds the *same source-level box specification* independently in
RustSolid and [SolidCore](https://github.com/tangency-labs/solidcore), a Rust port of
the NVIDIA USD BRep SMLib/SM_API components. Both kernels own their B-reps; neither
kernel is used to convert or repair the other's result.

## Run

From the RustSolid repository root:

\`\`\`bash
cargo test --manifest-path experiments/solidcore-parity/Cargo.toml --release
cargo run --manifest-path experiments/solidcore-parity/Cargo.toml --release
\`\`\`

The SolidCore git dependency is pinned to
\`ba7a0d9b8f4032441f656cec24e17677fa2668ec\` (the publicly recorded
2026-10-25 19:22 UTC progress snapshot). Do not replace it with a floating branch
during qualification.

## What is verified

For three deterministic source parameter sets, the executable:

1. Builds an independent box in each Rust kernel and validates RustSolid;
2. Requires SolidCore's own \`is_manifold_solid\` to succeed;
3. Checks independently obtained face/edge/vertex counts (6/12/8);
4. Computes each kernel's volume/area and compares both with analytic formulae;
5. Compares separate bounding-box queries with analytic boundaries;
6. Returns a failing exit code instead of silently accepting inconsistent results.

The independent analytical formulae are the primary geometric oracle; agreement
between two kernels alone is **not** a proof of correctness. This exercise has not
been independently compiled or run in the author's environment. CI must verify it.

## Integration policy

- SolidCore is a **reference/experimental provider** first, never the owner of
  RustSolid's canonical B-rep, source stamps, topology lineage or document history.
- Do not reinterpret \`solidcore::Brep\` arena IDs as RustSolid canonical IDs.
  If import is built later, use a validated explicit mapping for physical edges,
  face orientation, coedges, pcurves, seam multiplicity and tolerance semantics.
- Only a fully validated RustSolid \`BrepDraft\` with independent source and
  boundary checks could be considered for a future \`CertifiedCandidate\`.
  A passing SolidCore Boolean does not by itself authorize commit.
- Do not use upstream SMLib or this port as the sole oracle for general
  surface/surface intersection. Its [independent SSI corpus](https://tangency-labs.github.io/solidcore/ssi.html)
  records incorrect and incomplete results in both C++ and Rust.
- Beyond boxes, add cylinder, torus and NURBS fixtures, then Boolean, fillet,
  offset and healing **only after** each component passes source-based and
  metamorphic qualification.
- The standalone Cargo package is excluded from the root crate's normal build:
  no extra SolidCore download or new runtime dependency for RustSolid/Gefest.

## Licensing and provenance

[SolidCore](https://github.com/tangency-labs/solidcore) is Apache-2.0, based on
[NVIDIA USD BRep](https://github.com/NVIDIA-Omniverse/usd-brep)
(\`source/SMLib\`, \`source/SM_API\`; Apache-2.0). See its LICENSE,
NOTICE and ATTRIBUTION.md. RustSolid is MPL-2.0. Keep copyright notices and
license attribution when redistributing the dependency or derived files.
No NVIDIA implementation code is copied into this experiment.

## Qualification TODO

- CI green at the pinned SolidCore commit and fixed RustSolid base commit;
- Seeded parameter stress / near-degenerate cases with explicitly reconciled
  tolerance policies;
- Topology graph comparison (including physical coedges, pcurve correspondence,
  oriented normal consistency, seams and genus), not ID equality;
- Revision-safe conversion into a temporary canonical BrepDraft (if needed);
- Independent property and material-classification oracles for Boolean cases;
- Compare against OCCT and Truck; Parasolid only through a licensed SDK.
