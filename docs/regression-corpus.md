# RustSolid self-contained regression corpus (v0.11)

This is a first-party, **source-controlled** mathematical regression suite. It does
not treat OpenCascade, Parasolid, any third-party CAD output, or the kernel's own
`validate()` result as the correctness oracle.

## Two complementary suites

1. **Curated regressions** (`tests/fixtures/regressions.v1.json`) pin a small,
   reviewed set of inputs, known expected accept/reject outcomes, independently
   calculated analytic mass/bounds, and metamorphic relations. A previously
   broken thin-face construction is pinned with its exact failing width.
   Every release and PR builds these models again.
2. **Boundary probes** (`accuracy.rs`) perturb focused families with a seeded
   generator. The seed is an input for *discovery*, not a persistent pass/fail
   snapshot. Interesting failures should be promoted to curated JSON fixtures,
   including their minimal reproducer.

A valid model must pass `validate`, be convertible to the common analytic B-rep,
match analytic mass/bounds, and agree with its curve/UV trim carriers. A model
that must be rejected is considered correct only when it is **rejected**.

### Independent expectations

Known mathematical formulas, specified in the fixtures or evaluated independently
in the probe runner:

- Cuboid `volume = w*h*d`, `area = 2*(w*h + w*d + h*d)`.
- Prism `volume = polygon_area*h`; polygon centroid from anchored moments.
- Cylinder `volume = π*r²*h`, `area = 2π*r*(r+h)`.
- Regular `n`-gon inscribed in a cylinder has volume
  `n*r²*sin(2π/n)*h/2`, explicitly separated from **analytic** volume.

Check both acceptance and precision; a successful kernel operation is never
sufficient evidence of correctness.

### Metamorphic relations

These do not require another CAD kernel:

- Translation preserves volume and area and translates centroid/bounding box.
- Uniform scaling by `k` multiplies volume by `k³`, area by `k²`, positions by `k`.
- Reversing profile winding preserves geometric mass and topology counts.
- Increasing cylinder display facets reduces the *expected* discretization
  error, while analytic geometry and mass stay unchanged.

Each relation is restricted to supported model types; invalid fixture definitions
are rejected **before** running any geometry operation. Numeric comparison uses
relative mass bounds and coordinate-scale-aware positional tolerances.

## Reproduce and promote a failure

```bash
cargo test --all-targets
cargo test --test regression_corpus
printf '%s' '{"operation":"regression_corpus"}' | cargo run --quiet --bin gefest-geometry
printf '%s' '{"operation":"probe_kernel","seed":13}' | cargo run --quiet --bin gefest-geometry
```

The report contains `fixture`, individual checks, failure reasons, and a
`minimized_probe` where a deterministic, bounded profile-vertex simplifier can
preserve the failure category. This simplifier is conservative: it does not yet
minimize arbitrary 3D features or keep failed metamorphic relations while
removing features. Do not mistake it for a general delta debugger.

If a seeded probe fails: capture its `case.input`, promote it under a unique ID
in `tests/fixtures/regressions.v1.json`, record the independent expected result,
and add the smallest useful metamorphic relations. Never make a known failure
pass by removing its fixture or merely enlarging tolerances.

## Limits

This corpus is designed to grow with new operators; it is not a proof of
production-grade robustness. It does not currently evaluate freeform NURBS,
Boolean solids, fillet behavior, mesh self-intersection, property conservation
under general rotation, or arbitrary topology-edit histories. Extend the
oracle and fixtures alongside each new operation. Keep fixture IDs stable;
version the schema for incompatible changes.
