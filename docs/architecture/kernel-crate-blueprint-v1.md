# RFC 0002 — Complete crate architecture before further RustSolid extraction

**Status:** TARGET DESIGN ONLY; implementation NOT complete. **Baseline:** `53d8a8e3385e4b67076536465706e2bb34c3cc0f` (G8c / PR #41).

The current implementation branch includes 7 partially extracted crates and a transitional root `rustsolid` package. This RFC defines 20 eventual workspace packages (19 library crates and the `gefest-geometry` app) and maps every baseline 89 `src/**/*.rs` module. See `kernel-crate-plan-v1.json` for machine-readable ownership, dependency limits and mandatory splits.

## Why design precedes extraction

Moving files by name produced false dependencies. Existing `brep_arena.rs` includes Boolean dispatch, a sparse registry, document revisions, and Undo/Redo; `unified_brep.rs` combines a canonical B-rep view and special-case source constructors; volume oracles import G4/G7/G8 candidate plans; and mesh builders sometimes call editors. Whole-file extraction would create cycles and undermine independent certification. Physical extraction is frozen until these ownership boundaries and contracts are reviewed.

**Public API breakage is allowed.** Update Gefest and dependent callers instead of preserving the monolith via compatibility facades. Preserve geometrical invariants, source provenance and mathematically independent G7/G8 tests, not the old package paths.

## Proposed workspace (design maximum direct dependency edges)

The dependency lists are **allowed boundaries**, not requirements to eagerly add every listed dependency to every crate. Actual manifests must use the smallest subset. A package is instantiated only when it has executable independent contracts.

| Crate | State | Responsibilities | Allowed dependencies |
| --- | --- | --- | --- |
| `rustsolid-foundation` | partial | Coordinates, tolerance configuration and base errors | none |
| `rustsolid-predicates` | partial | Robust orientation and filtered numeric decisions | `rustsolid-foundation` |
| `rustsolid-geometry` | partial | Analytic curves/surfaces and UV chart domains | `rustsolid-foundation` |
| `rustsolid-uv` | partial | Source-neutral analytic UV arrangement and polygon-only geo adapter | `rustsolid-geometry`, `rustsolid-predicates` |
| `rustsolid-intersect` | partial | Untrimmed curve/surface and surface/surface intersections | `rustsolid-foundation`, `rustsolid-geometry`, `rustsolid-predicates` |
| `rustsolid-topology` | partial | Canonical B-rep graph, geometry attachments, scoped IDs and invariants | `rustsolid-foundation`, `rustsolid-geometry`, `rustsolid-predicates` |
| `rustsolid-source` | partial | Immutable source geometry parameters, SourceId and typed BooleanIntent; never constructs B-rep | `rustsolid-foundation`, `rustsolid-geometry` |
| `rustsolid-construct` | planned | Source primitive B-rep builders (box/cylinder/extrude, bores); source parameters live in rustsolid-source | `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-predicates`, `rustsolid-source` |
| `rustsolid-imprint` | planned | Trimmed source-face intersections, G1 graph and provenance | `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-intersect`, `rustsolid-uv`, `rustsolid-predicates`, `rustsolid-source` |
| `rustsolid-classification` | planned | G2/G3 cell and whole-solid material classification | `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-imprint`, `rustsolid-uv`, `rustsolid-predicates`, `rustsolid-source` |
| `rustsolid-validation` | planned | Independent source-volume brackets, result flux and invariants; no dependency on constructors | `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-uv`, `rustsolid-predicates`, `rustsolid-source` |
| `rustsolid-assembly` | planned | G4 source-bound topology changes, patch stitching and staged candidate B-rep | `rustsolid-construct`, `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-imprint`, `rustsolid-uv`, `rustsolid-classification`, `rustsolid-validation`, `rustsolid-source` |
| `rustsolid-boolean` | planned | Certified staged Boolean execution, no document commit | `rustsolid-construct`, `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-imprint`, `rustsolid-classification`, `rustsolid-assembly`, `rustsolid-validation`, `rustsolid-source` |
| `rustsolid-modeling` | planned | Euler operations, structural/carrier edits and feature intents | `rustsolid-construct`, `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-imprint`, `rustsolid-validation`, `rustsolid-boolean`, `rustsolid-source` |
| `rustsolid-meshing` | planned | Derived face-owned display mesh and facet certificates | `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-uv`, `rustsolid-validation` |
| `rustsolid-document` | planned | Atomic commit, persistent naming and Undo/Redo | `rustsolid-topology`, `rustsolid-construct`, `rustsolid-validation`, `rustsolid-boolean`, `rustsolid-modeling`, `rustsolid-source` |
| `rustsolid-archive` | planned | Versioned feature archive and reopen by recomputation | `rustsolid-document`, `rustsolid-geometry`, `rustsolid-source` |
| `rustsolid-query` | planned | Bounded read-only geometry query language | `rustsolid-topology`, `rustsolid-geometry` |
| `rustsolid-qualification` | planned | Independent G4–G8 fixtures, metamorphic tests and benchmarks | `rustsolid-archive`, `rustsolid-query`, `rustsolid-document`, `rustsolid-modeling`, `rustsolid-boolean`, `rustsolid-assembly`, `rustsolid-classification`, `rustsolid-imprint`, `rustsolid-meshing`, `rustsolid-validation`, `rustsolid-construct`, `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-intersect`, `rustsolid-uv`, `rustsolid-predicates`, `rustsolid-source` |
| `gefest-geometry` | planned | Gefest JSON/stdio CLI and adapters | `rustsolid-archive`, `rustsolid-query`, `rustsolid-document`, `rustsolid-modeling`, `rustsolid-boolean`, `rustsolid-meshing`, `rustsolid-construct`, `rustsolid-topology`, `rustsolid-geometry`, `rustsolid-intersect`, `rustsolid-classification`, `rustsolid-source` |

Dependency direction follows: numerical and analytic leaf crates → canonical model/source construction → imprint and UV → whole-solid classification → candidate assembly and independent validation → staged Boolean → modeling and transactions → archive/CLI/qualification. `validation` is a peer of `assembly`, never its client. The transitional root is not a permanent package in this graph.

## Types and authority: who owns the truth?

**`rustsolid-foundation`:** Primitive length/coordinate/tolerance semantics only. The current single `GeometryError` is a transition artifact: CAD-operation-specific failures move to owner crates as typed `Held`, `Unsupported` and invariant failures. Model units are defined in geometry; world axes and UV `(u,v)` must never be silently swapped.

**`rustsolid-geometry`:** Exact-in-form `Curve2`, `Curve3`, `Surface3`, harmonic and circular carriers, surface evaluation and chart domains. No reference to `Solid`, `BrepModel` or `FaceId`. `rustsolid-predicates` owns adaptive exact-sign tests separately from engineering-tolerance classification; analytic SSI belongs to `rustsolid-intersect`.

**`rustsolid-topology`:** Canonical `BrepModel` (geometry catalogs, coedge-specific PCurve, face/loop/coedge/edge/shell incidence, source IDs and orientation) plus process-scoped generational handles. Its topological identity is distinct from document feature identity. No `Mesh` stored as canonical body geometry; display triangulation is derived and may be cached externally. Canonical geometry projections from the old `GeometryStore` adapter migrate here; constructors migrate to `rustsolid-construct`.

**`rustsolid-source`:** Immutable `SourceSpec`, source parameters, `SourceId`, and typed `BooleanIntent`/operation descriptions. This leaf package depends on geometry and foundation only: it cannot construct B-rep, query a document, calculate the mass of a candidate, or run Boolean.

**`rustsolid-construct`:** Given immutable sources, construct canonical box/cylinder/extrusion/bore operands using `rustsolid-topology`. It owns no primary monolith file yet for the shared source types: they are carved out of existing mixed files. It must not depend on `rustsolid-validation`.

**`rustsolid-uv`:** Source-neutral analytic chart arrangement, including periodic seam handling, circular holes and exact Curve2 ownership. The *existing* `geo` adapter is for polygonal operations only: no automatic replacement of a circle/ellipse by faceted topology. Face-pair graphs and physical source ownership belong in `rustsolid-imprint`, which adapts its geometry into UV input contracts.

**`rustsolid-imprint`:** Trimmed 3D/2D intersections and G1 face-pair/global event graphs with source-certified coedge pcurves, vertex/edge reconciliation and typed ambiguity. An intersection event is merged only when source topology supplies a verifiable common physical object, never by proximity alone. No mutation or Boolean region selection.

**`rustsolid-classification`:** G2/G3 source-bound UV cell membership, in/out/boundary and global material classification. It is read-only and independent of the Boolean executor. Sampling alone does not prove uniform cell membership; ambiguous boundary events remain Held.

**`rustsolid-assembly`:** G4 topology-change intent, physical edge splits, face patch selection, stitching and preflight candidate shell construction. It never publishes document IDs or claims a successful Boolean solely because the candidate B-rep validates.

**`rustsolid-validation`:** Depends on immutable `rustsolid-source` data, **never on the `rustsolid-construct` builder package**. Independently checks source-only volume brackets against flux over the *actual* oriented B-rep boundary and reports topology and source-provenance violations. G8b source integration must not read the G8c candidate, and G8c flux must not reuse source-volume calculations or an assembler's declared mass. Strictly forbidden dependencies: `construct`, `assembly`, `boolean`, `classification`, `imprint`, `meshing`, `modeling`, `document`, `archive`.

**`rustsolid-boolean`:** Orchestrates staged G1–G4 and independent certification, returning an immutable certified candidate or typed failure. It does not own Undo/Redo, project schemas, rendering or a public commit API. Restricted through/blind, multi-cut, oblique, one-wall and two-wall executors are internal algorithms, not separate libraries.

**`rustsolid-modeling` and `rustsolid-document`:** Modeling owns Euler, local face/edge, carrier and structural edit transactions; `document` owns authoritative revisions, persistent feature lineage, bounded Undo/Redo and the mutation barrier. Generation-bound process handles must never be used as serialized feature identity.

**`rustsolid-meshing`:** Shared analytic physical edge samples, triangulation and per-triangle topological face ownership. Faceting must never determine canonical volume, shape or edges. It may depend on validation, but no kernel constructor or Boolean stage may call mesh builders.

**`rustsolid-archive`:** Versioned source/feature DAG, immutable parameter snapshots and independent evidence. Reopen replays document operations and re-certifies the model; serialized B-rep meshes or transient handles cannot authorize accepting a shape.

**`rustsolid-query`, `rustsolid-qualification`, `gefest-geometry`:** Read-only bounded model querying, offline/regression/stress/benchmark harness, and the production Gefest JSON client protocol respectively. Runtime production crates may not depend on `qualification`; diagnostic CLI commands become a separate opt-in executable.

## Shared type ownership (breaking-contract migration)

Naming a crate is insufficient if types still originate in a different layer. This table determines the target owner of the key public types. These are **design assignments**, not claims that code was already moved.

| Type or invariant | Target owner |
| --- | --- |
| `Point3 / Point2 (generic Cartesian)` | `rustsolid-foundation` |
| `GeometryTolerance` | `rustsolid-foundation` |
| `Orientation / PredicateClassification` | `rustsolid-predicates` |
| `Plane3 / PlaneIntersection` | `rustsolid-predicates` |
| `ModelUnits / ParameterRange / Point2Param` | `rustsolid-geometry` |
| `Curve2 / Curve3 / Surface3` | `rustsolid-geometry` |
| `BoundingBox (Aabb3) / GeometricMoments` | `rustsolid-geometry` |
| `Curve2Id / Curve3Id / SurfaceId (catalog handles)` | `rustsolid-topology` |
| `VertexId / EdgeId / CoedgeId / FaceId / LoopId / ShellId` | `rustsolid-topology` |
| `FaceSense / BrepModel / GeometryStore` | `rustsolid-topology` |
| `TopologyIssue / RuntimeGenerationalHandle` | `rustsolid-topology` |
| `SourceId / SourceSpec / PrimitiveParameters` | `rustsolid-source` |
| `BooleanOperation / BooleanIntent` | `rustsolid-source` |
| `SourceFaceSemanticRole (e.g. TopCap)` | `rustsolid-source` |
| `IntersectionCurve / IntersectionCertificate` | `rustsolid-intersect` |
| `TrimmedIntersection / GraphArc / ReconciledGlobalEvent` | `rustsolid-imprint` |
| `PhysicalCutProvenance / DualSourcePcurve` | `rustsolid-imprint` |
| `AnalyticUvCell / UvBoundaryCycle` | `rustsolid-uv` |
| `CertifiedCellLocation / WholeSolidCellWitness` | `rustsolid-classification` |
| `StagedTopologyDelta / FacePatch / GlobalStitchPlan` | `rustsolid-assembly` |
| `IndependentSourceVolumeBracket / ActualBoundaryFlux` | `rustsolid-validation` |
| `ShellValidationCertificate / InvariantViolation` | `rustsolid-validation` |
| `CertifiedCandidate / BooleanHeld` | `rustsolid-boolean` |
| `EditIntent / EulerCommand` | `rustsolid-modeling` |
| `FaceOwnedMesh / MeshCertificate` | `rustsolid-meshing` |
| `FeatureId / PersistentTopologyLineage` | `rustsolid-document` |
| `DocumentRevision / UndoRedoTransition` | `rustsolid-document` |
| `ProjectArchive / ReopenSchema` | `rustsolid-archive` |
| `QueryLimits / QueryError` | `rustsolid-query` |
| `RegressionFixture / StressScenario` | `rustsolid-qualification` |
| `GeometryRequest / GeometryResponse (CLI)` | `gefest-geometry` |

Mandatory intentional API breakage: `Curve2Id/Curve3Id/SurfaceId` are currently defined alongside pure geometry, but they identify entries in a *B-rep-owned catalog* and must move to topology. Likewise `FaceRole::TopCap/BottomCap/Side` is a source-semantic label, not a universal face-topology property; preserve `FaceSense` in topology and record semantics through source provenance. `Mesh` and its face-owned triangle IDs are derived in meshing, never embedded in canonical B-rep. A persisted `FeatureId` is not a process-local `TopoId` or generational handle.

Sources and certificate states must cross crate boundaries as immutable, versioned input/result values. Error ownership stays local to each stage; do not perpetuate the current all-purpose `GeometryError::InvalidEdit` for classification ambiguity or Boolean rejection.

## Proposed cross-crate protocol (not implemented Rust signatures)

```text
Source:          SourceSpec + BooleanIntent -> ImmutableSourceOperands
Construct:       ImmutableSourceOperands -> ValidatedSourceBrep
Intersections:   (Surface3, Surface3) -> UntrimmedBranches | Held
Imprint G1:      (BrepModel, BrepModel, Branches) -> ProvenantFacePairGraph
UV G2:          AnalyticChartTraces -> AnalyticUvCells | Held
Classify G3:    (AnalyticUvCells, ImmutableSources) -> CertifiedMaterialCells | Held
Assembly G4:    (ProvenantGraph, ClassifiedCells) -> ProvisionalBrep + Lineage
Validate G4c:   (ImmutableSourceParams, ActualProvisionalBrep) -> IndependentCertificate
Boolean:        (BooleanIntent, ValidatedSourceBrep) -> CertifiedCandidate | Held | Unsupported
Document G5:    (ExpectedRevision, CertifiedCandidate) -> AtomicCommit | Conflict
Meshing:        (CanonicalBrepModel, DisplayTolerance) -> FaceOwnedMesh
Archive:        (VersionedFeatureSources) -> ReplayedAndRecertifiedDocument | Rejected
```

Each `CertifiedCandidate` must carry the accepted original source lineage, complete physical coedge ownership, topology validation, the source-volume interval and the independent boundary-flux check. Reject-without-mutation is mandatory. The current numeric `f64` error guards are not advertised as rigorous interval-arithmetic proofs for general NURBS.

## Baseline source ownership ledger

Exactly one **primary target** is assigned to every existing Rust source module; multiple code destinations must be declared separately below. `rustsolid-source` intentionally has **zero primary legacy modules** because its DTOs are extracted from multiple modules (`boolean.rs`, `through_hole.rs`, `blind_hole.rs`, etc.), not migrated wholesale. `retire` files are transitional shims to delete after migration.

| Target package | Files | Modules (`src/` relative) |
| --- | ---: | --- |
| `rustsolid-qualification` | 2 | `accuracy.rs`, `regression.rs` |
| `rustsolid-assembly` | 11 | `analytic_bore_candidate.rs`, `analytic_face_patch_partition.rs`, `analytic_one_wall_notch.rs`, `analytic_source_edge_split.rs`, `analytic_tilted_one_wall_brep.rs`, `analytic_two_wall_partition.rs`, `canonical_box_candidate.rs`, `face_patch_assembly.rs`, `global_edge_stitch.rs`, `mixed_crossing_surfaces.rs`, `topology_delta_plan.rs` |
| `rustsolid-document` | 5 | `atomic_boolean_document.rs`, `brep_arena.rs`, `g5_two_wall_document.rs`, `g7_one_wall_document.rs`, `history.rs` |
| `gefest-geometry` | 11 | `bin/composed_crossing_cli/mod.rs`, `bin/connected_triple_cli/mod.rs`, `bin/g5_project_cli.rs`, `bin/g6_project_cli.rs`, `bin/g7_project_cli.rs`, `bin/gefest-geometry.rs`, `bin/mixed_crossing_cli/mod.rs`, `bin/mixed_footprint_cli/mod.rs`, `bin/mixed_regions_cli/mod.rs`, `bin/multi_cut_regions_cli/mod.rs`, `bin/oblique_cylinder_cli.rs` |
| `rustsolid-construct` | 4 | `blind_hole.rs`, `curved_brep.rs`, `multi_through_hole.rs`, `through_hole.rs` |
| `rustsolid-classification` | 4 | `boolean_region.rs`, `mixed_footprint.rs`, `multi_face_region_preview.rs`, `whole_solid_region_classification.rs` |
| `rustsolid-boolean` | 9 | `boolean.rs`, `composed_crossing_brep.rs`, `composite_blind.rs`, `compound_blind.rs`, `connected_triple_brep.rs`, `intersecting_blind.rs`, `mixed_crossing_brep.rs`, `mixed_regions.rs`, `oblique_halfspace.rs` |
| `rustsolid-topology` | 5 | `brep.rs`, `geometry.rs`, `identity.rs`, `predicates.rs`, `unified_brep.rs` |
| `rustsolid-modeling` | 7 | `carrier_editor.rs`, `closed_cylinder_imprint.rs`, `edit.rs`, `euler.rs`, `journal.rs`, `operator_contract.rs`, `structural_editor.rs` |
| `rustsolid-imprint` | 8 | `coedge_reconciliation.rs`, `connected_cut_transitions.rs`, `global_analytic_cut_registry.rs`, `global_vertex_reconciliation.rs`, `intersection_graph.rs`, `multi_face_intersection_graph.rs`, `tilted_one_wall_source.rs`, `trimmed_intersections.rs` |
| `rustsolid-meshing` | 9 | `composed_crossing_mesh.rs`, `compound_blind_mesh.rs`, `connected_triple_mesh.rs`, `face_region_mesher.rs`, `imprint_mesh.rs`, `intersecting_blind_mesh.rs`, `mixed_crossing_mesh.rs`, `mixed_regions_mesh.rs`, `oblique_clip_mesh.rs` |
| `rustsolid-archive` | 2 | `g5_project_archive.rs`, `g7_project_archive.rs` |
| `rustsolid-validation` | 6 | `independent_curved_volume_oracle.rs`, `independent_one_wall_volume.rs`, `independent_tilted_boundary_flux.rs`, `independent_tilted_volume.rs`, `independent_two_wall_volume.rs`, `radial_shell_audit.rs` |
| `retire` | 2 | `intersections.rs`, `lib.rs` |
| `rustsolid-uv` | 3 | `multi_cut_regions.rs`, `multi_trace_uv_arrangement.rs`, `planar_uv_arrangement.rs` |
| `rustsolid-query` | 1 | `query.rs` |

## Mandatory splits of mixed-layer files

**Do not move these files wholesale.** Define the lower-level types first. In particular source-only volume oracles must depend on `rustsolid-source` and actual B-rep carrier geometry, not constructors or candidate assembly plan objects.

| Baseline file | Required decomposition |
| --- | --- |
| `src/accuracy.rs` | move probes to non-production qualification, not CLI release dependency |
| `src/blind_hole.rs` | BlindHoleParameters -> rustsolid-source; source B-rep builder -> construct; transaction -> document/modeling; display -> meshing. |
| `src/boolean_region.rs` | material membership in classification; BooleanIntent lives below bool executor |
| `src/boolean.rs` | BooleanOperation and immutable operands -> rustsolid-source; staged CSG executor -> rustsolid-boolean. |
| `src/brep_arena.rs` | registry into topology; command executor/history into document; modeling commands into modeling |
| `src/closed_cylinder_imprint.rs` | read-only source curves into imprint; mutation into modeling |
| `src/compound_blind.rs` | Immutable blind region source specifications -> rustsolid-source; bounded B-rep Boolean -> rustsolid-boolean; transaction wrapper -> document. |
| `src/connected_triple_brep.rs` | Source parameter struct -> rustsolid-source; bounded multi-cut candidate builder -> rustsolid-boolean; edit history -> document. |
| `src/curved_brep.rs` | Cylinder source description -> rustsolid-source; canonical capped-cylinder constructor -> rustsolid-construct; display-only mesh -> rustsolid-meshing. |
| `src/face_region_mesher.rs` | accept BrepModel alone, not BrepStructuralEditor |
| `src/g5_two_wall_document.rs` | document state vs persisted archive schema/replay split |
| `src/geometry.rs` | remaining Solid-to-GeometryStore and carrier bindings into topology; pure curves already in geometry |
| `src/identity.rs` | runtime generational handles into topology; persistent feature lineage into document |
| `src/independent_curved_volume_oracle.rs` | replace G4 candidate input with immutable original source parameters + actual output B-rep |
| `src/independent_one_wall_volume.rs` | replace G7a plan input with immutable source parameters + actual B-rep flux |
| `src/independent_tilted_boundary_flux.rs` | G8c uses result B-rep only, never Boolean construction or G8b integration code |
| `src/independent_tilted_volume.rs` | G8b uses source geometry alone, not higher-stage plan type |
| `src/independent_two_wall_volume.rs` | replace G6 plan input with source primitives + actual B-rep |
| `src/intersections.rs` | delete already obsolete root re-export after callers migrate |
| `src/journal.rs` | topology edit inverse log into modeling; remove mesh edit coupling; committed history into document |
| `src/lib.rs` | delete legacy root exports at final cut; no permanent compatibility facade |
| `src/mixed_regions.rs` | Through/blind source parameters -> source; Boolean shell builder -> boolean; preview display -> meshing; transactions -> document. |
| `src/multi_trace_uv_arrangement.rs` | analytic chart/cell construction in uv; face-pair binding in imprint |
| `src/oblique_halfspace.rs` | Halfspace intent and side type -> source; CSG executor -> boolean; remove display dependency. |
| `src/operator_contract.rs` | Operation intent types -> source (when shared); executable Euler/feature commands -> modeling, coverage registry -> qualification. |
| `src/planar_uv_arrangement.rs` | analytic UV operations source-neutral in uv; source/topology adapter into imprint |
| `src/predicates.rs` | Solid-aware face support plane into topology; pure predicates already extracted |
| `src/radial_shell_audit.rs` | read-only manifold preflight in validation, generic topology inputs |
| `src/regression.rs` | move goldens, generators and shrinkers to qualification |
| `src/through_hole.rs` | HoleParameters -> rustsolid-source; immutable B-rep builder -> construct; transaction -> document/modeling; display -> meshing. |
| `src/unified_brep.rs` | canonical model into topology; concrete shape converter into construct; remove embedded display |
| `src/whole_solid_region_classification.rs` | G3 classification cannot import BooleanKernel |

## External library review (scope, not new dependencies yet)

| Library | Integration boundary | Adoption gate |
| --- | --- | --- |
| `robust` | predicates, already connected | exact-sign and tolerance semantics |
| `geo` | polygon-only UV adapter | polygon validity and analytic-source error budget; never silently faceted B-rep |
| `earcutr` | meshing only | shared physical edges, per-face ownership, consistent orientation |
| `rstar` | optional future spatial broad phase | measured acceleration and deterministic candidate ordering before making `rustsolid-spatial` |
| Truck (`truck-*`) | qualification/reference backend first | check licensing and run common G7/G8 regressions before replacing internal geometric types |
| STEP/IGES stack | future `rustsolid-step` adapter | schema/format import requirements and independent topology/shape healing corpus |

No empty packages for future NURBS, constraints, spatial indexing or STEP IO without executable requirements; extend existing geometry/intersect contracts first.

## Binding four cross-crate contracts (RFC 0003)

The full SourceSpec, BrepModel, TopologyLineage and CertifiedCandidate design is in [RFC 0003](kernel-contracts-v1.md); the [JSON contract inventory](kernel-contracts-v1.json) and [synthetic fail-closed scenarios](kernel-contract-scenarios-v1.json) are checked by `tests/validate_kernel_contracts.py`. **The tests inspect contract descriptions, not actual CAD implementation**. New physical extraction remains blocked until source versioning, mesh-free topology, independent G8 validation, persistent feature naming and commit semantics have reviewable implementation tests.

## CI and migration rules

1. Architecture gate: `python3 tests/validate_crate_blueprint.py` checks declared target packages, acyclic design, mandatory file splits, extracted manifest presence and coverage of *all current* `src/**/*.rs` paths. It does not imply unimplemented packages compile.
2. For an actual new package: update `tests/check_workspace_boundaries.py` (the real `cargo metadata` DAG), allocate one independent owner, provide positive/negative tests, and then add its matrix CI job.
3. After each extraction: `cargo test -p <new-package> --all-targets`, `cargo test --workspace --all-targets`, G7/G8 source-vs-boundary-flux assertions, shell/genus/coedge consistency, archived Undo/Redo/replay and Gefest JSON smoke tests.
4. Preserve *geometrical behavior*, not legacy Rust paths; never add a catch-all compatibility crate. Stage the algorithm migration separately from replacement of the algorithm by a third-party library.
5. Benchmark per-crate compile times, wall time of full workspace and incremental rebuild, peak memory, and test runtime. Crate splitting primarily controls ownership and test impact; speed claims require measurements.
6. The canonical model must be mesh-free before cutting Boolean/modeling/document. The source oracle must be independent before moving G4c/G7b/G8b/G8c into `validation`. If either invariant fails, **stop that extraction**.

## Sequencing

**Design gate (this RFC)** → immutable `source` contract + canonical topology representation + primitive construction → G1/imprint and source-neutral UV → whole-solid classification (G2/G3) → G4/assembly and independent source vs boundary validation → staged Boolean → modeling + mesh in parallel → document atomic history → archive/query/qualification/CLI adaptation → remove transitional root.

**Deferred, not allocated:** `rustsolid-spatial` after measured needs, `rustsolid-step` after import scope, `rustsolid-constraints` for sketch solving, `rustsolid-truck-adapter` for experimental comparison.

**Review checkpoint:** This is the full package ownership plan for current G4–G8 scope, not a declaration that all 19 packages exist. A new feature affecting cross-crate type ownership must update this RFC and the machine ledger *before* code is extracted.
