# Parasolid official capability ↔ PK API ↔ RustSolid audit (2026-10-11)

This is a **versioned evidence map**, not a port of Parasolid's implementation.
Machine-readable, 14-row ledger: [parasolid-official-capabilities.v1.json](parasolid-official-capabilities.v1.json).
Current RustSolid comparison baseline: `v0.14.0`, main as inspected on 2026-10-11.

## Evidence levels and a critical version boundary

**Level A — confirmed public claim:** [Siemens official Parasolid product page](https://www.siemens.com/en-us/products/plm-components/parasolid/), [Siemens developer overview](https://developer.siemens.com/plm-components/overview.html), and dated Siemens technical release announcements.
The published product description says **900+ functions** and covers solid, sheet, surface/direct and mixed modeling. That is a vendor's product-level claim, **not** a public, function-by-function specification of PK signatures and preconditions.

**Level B — previously reported export-name observations:** our prior architectural note
[unified-brep-parasolid-comparison.md](unified-brep-parasolid-comparison.md)
identifies names from a user-supplied decompilation reported as Parasolid **38.01.207**.
The decompiled binary and vendor SDK reference manual are **not available for reinspection in this audit**.
These names therefore serve as *leads pending verification*, not as newly verified symbols, internal call graphs, official ABI declarations, or proof of underlying numerical algorithms. A PK export count need not equal the marketed number of documented SDK functions.

**Level C — independent RustSolid evidence:** modules and tests actually present in `main`.
`BrepModel`, `GeometryStore`, `SurfaceIntersectionEngine`,
`TrimmedIntersectionEngine`, `CommandHistory` and the pinned regression corpus.

**Version guard:** Siemens describes general-body and tracking improvements for **Parasolid 38.0** (August 4, 2025), while the new mesh↔B-rep mapping, mesh Euler utilities, configurable edge curve reuse and enhanced parallel enquiries were announced for **39.0** (August 25, 2026). Do **not** infer that features announced in 39.0 existed in 38.01.207 simply because they appear in the newer product description.

## Concrete crosswalk

| Officially declared capability and source | Prior PK API lead (NOT SDK-verified here) | RustSolid v0.14 and missing contract |
|---|---|---|
| [General body operations, v38.0](https://blogs.sw.siemens.com/plm-components/whats-new-in-parasolid-v-38-0/) | `PK_FIN_ask_next_in_loop`, `PK_FIN_ask_next_of_edge`, `PK_SESSION_ask_general_topology` | Directed loops/twins only for a closed manifold shell; need radial cycles, sheets, regions and local manifold checks |
| [B-rep and trimmed B-spline data, 2018](https://blogs.sw.siemens.com/plm-components/parasolid-data-exchange/) | `PK_FACE_ask_oriented_surf`, `PK_EDGE_ask_curve`, `PK_FIN_find_surf_parameters` | Analytic planes/cylinders and UV p-curves; need general B-spline carriers, multiple trims and robust reparameterization |
| [Tolerant modeling and import, 2018](https://blogs.sw.siemens.com/plm-components/parasolid-data-exchange/) | `PK_EDGE_ask_precision` | Tolerances exist but are not propagated and certified during topology-changing operations |
| [Boolean region tracking, v32.0](https://blogs.sw.siemens.com/plm-components/parasolid-v32-0-release-highlights/) | `PK_BODY_boolean` | No Boolean; require a typed `IntersectionGraph → RegionClassification → TopologyDelta` pipeline with lineage witnesses |
| [General geometry interrogation/sectioning, 2018](https://blogs.sw.siemens.com/plm-components/parasolid-data-exchange/) | `PK_SURF_intersect_surf` | Plane/plane and plane/cylinder support curves with restricted UV trimming; need broader pair coverage, completeness and topology-safe uncertainty |
| [Reported extensions and face-pattern origin tracking, v38.0](https://blogs.sw.siemens.com/plm-components/whats-new-in-parasolid-v-38-0/) | Specific function signatures not independently confirmed | Add optional reports for carrier reuse/new creation, trim extensions and source-to-derived face associations |
| [Classic/facet mixed models, v33.0](https://blogs.sw.siemens.com/plm-components/parasolid-v33-0-release-highlights/) | Previously noted `PK_MESH_*` / `PK_MFACET_*` family | Current rendering mesh is separate; no mixed, editable geometric region |
| [B-rep↔mesh associations and mesh Euler edits, v39.0](https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/) | **No v38.01.207 equivalence asserted** | Need persistent two-way mapping between fins/facets and analytic entities, local retessellation |
| [Edge curve creation/reuse policy, v39.0](https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/) | **No v38.01.207 equivalence asserted** | Need explicit edit policy and handle/name provenance, not implicit duplication |
| [Concurrent enquiries and partition transmissions, v39.0](https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/) | Earlier notes mention `PK_MARK_*`, `PK_PARTITION_*`, without rechecking concrete signatures | Single-body command history; future immutable inquiry snapshots, partitions, isolation and durable deltas |
| [Repair/geometry clean-up, 2018](https://blogs.sw.siemens.com/plm-components/parasolid-data-exchange/) | No exact observed repair symbol verified | Structured diagnostics and valid/invalid probes, but not tolerance-governed repair with proof of preservation |
| [Regression-based quality, 2020](https://blogs.sw.siemens.com/plm-components/parasolid-delivers-30-years-of-quality-in-industrial-use/) | No binary mapping required | Our first-party golden, metamorphic and boundary corpus is a starting point, not production-scale coverage |

The complete ledger lists a falsifiable mathematical acceptance test for each capability.
All newly proposed tests must remain **independent of Parasolid and OpenCascade as correctness oracles**.

## The architectural patterns that matter for Gefest

1. **Separate entity navigation, geometric carriers and parameter-space trims.** The public topology/geometry API shape is consistent with maintaining distinct, independently verified contracts. It does not reveal Siemens' private algorithms.
2. **Do not conflate precise analytic geometry with rendered facets.** Newer Siemens support for bidirectional mapping demonstrates a useful API boundary for selection, tessellation, mesh edits and source-face lineage.
3. **Model operation effects as a typed delta with provenance.** Origin tracking, retained/deleted/created region reports, and policy for reusing versus replacing carriers matter for direct modeling and stable selection.
4. **Treat tolerance as local provenance, not a universal epsilon.** Detect and *report* ambiguity, propagate admissible local error bounds through constructive edits, and distinguish exact mathematical predicates from inexact constructed coordinates.
5. **Generalize validity by region and topology kind.** Euler chi=2 for a single genus-zero shell is not a universal criterion for general-body regions.
6. **Separate read-only enquiry snapshots from mutating modeling sessions.** This supports concurrent agents, viewers and a future geometric query API without leaking partially changed topology.

## How to maintain this map

For each RustSolid feature or Parasolid release:
- Identify the **exact official page, publication date, and feature/release version**.
- If function-level comparison is needed, consult a legally accessible, version-matched Parasolid API manual or re-examine the matching exported symbol in an authorized research environment. Record whether observation is **verified**, **inferred**, or only present in **previous notes**; do not invent signatures or algorithms.
- Translate the capability into a **RustSolid-owned contract** with inputs, output classification, tolerance/uncertainty, stable identifiers, provenance, and independently specified regression expectations.
- Update the [machine-readable ledger](parasolid-official-capabilities.v1.json). CI validates its source provenance and explicit version labels but cannot certify vendor claims, proprietary behavior, or unmounted binaries.
- Keep the independent production implementation in Rust; don't introduce decompiled Siemens code, private layouts, or substantial copied implementations.

## Immediate consequences for RustSolid's next feature

Before a Boolean, design an `IntersectionGraph` that records every intersection vertex/branch, which source faces and carriers generated it, whether it is transverse/tangent/overlap/ambiguous, and the tolerances of all supporting observations. Generate a typed `TopologyDelta` only after classification, manifold checks and replayable provenance have passed first-party invariants. This combines the verified v0.14 UV trimming work with independently testable origin tracking and local-tolerance contracts.

The 39.0 mesh lookup/reuse/concurrency innovations are targets for future releases, **not** evidence about internals of the reported 38.01.207 binary.
