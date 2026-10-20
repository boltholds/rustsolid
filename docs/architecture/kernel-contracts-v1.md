# RFC 0003 — RustSolid G4–G8 cross-crate contracts v1

**Status:** PARTIAL IMPLEMENTATION. SourceSpec/SourceSnapshot implemented as a new crate in the stacked implementation PR; the remaining contracts are not yet production-integrated.  
**Parent architecture:** [RFC 0002](kernel-crate-blueprint-v1.md) and [source ownership ledger](kernel-crate-plan-v1.json).  
**Base implementation investigated:** G5, G7 and G8c on PR #41, head before this RFC \`45d9fee\`.  
**Machine-readable agreement:** [kernel-contracts-v1.json](kernel-contracts-v1.json).  
**Contract scenarios:** [kernel-contract-scenarios-v1.json](kernel-contract-scenarios-v1.json).

The public monolithic RustSolid API is **not** a compatibility target. Gefest consumers should migrate to new explicit types. However, existing G4–G8 geometric invariants, exact-in-form curve carriers, source ownership and independent volume/flux evidence remain acceptance requirements. This RFC settles four cross-crate contracts before more code moves.

## 0. Existing implementation evidence and mismatches

- **Current \`BrepModel\`** (\`src/unified_brep.rs\`) already stores vertices/edges/coedges/loops/faces/shells, curves3/curves2/surfaces, per-coedge geometry, tolerance, units and \`FaceSense\`; it also holds **derived mass, bbox and producer-specific \`BrepOrigin\`**. Its validator keys several orientation requirements off \`BrepOrigin\` and face indices, which is not a general B-rep invariant.
- **Current G7** (\`src/g7_one_wall_document.rs\`) has \`prepare\` → revision-, document-token- and source-stamp-bound \`commit\`. It checks candidate and an independent source/flux witness before publication and maintains Undo/Redo. Preserve the behavior, not separate G5/G6/G7 document classes.
- **Current G8c** (\`src/analytic_tilted_one_wall_brep.rs\`) carries actual physical G1 cut-edge ↔ source endpoint/face mappings, with \`can_apply = false\`. Neither a closed shell nor a staged volume alone authorizes commit.
- **Current G8b/G8c oracles** (\`src/independent_tilted_volume.rs\`, \`src/independent_tilted_boundary_flux.rs\`) import source/shell *plan types*. These references must be removed before independent validation can become a lower-level crate. G8c's flux integration must read the **actual candidate BrepModel**; the source interval is calculated from the original source parameters by a distinct path.
- **Current \`TopologyName\` and \`TopologyHandle\`** (\`src/identity.rs\`) combine feature naming with process-local generations. New persistent feature/lineage names move to \`document\`; local slot identities remain in \`topology\`.

## 1. Contract A — \`SourceSpec\`: authoritative immutable inputs

**Owner:** \`rustsolid-source\` (planned), which depends only on foundation and geometry. This crate **must not** import \`BrepModel\`, \`construct\`, \`boolean\`, document history or validation.

Design sketch (the declarations are deliberately **not compiled Rust** yet):

\`\`\`rust
pub struct SourceId(Uuid128);
pub struct SourceRevision(u64);

pub struct SourceSpec {
    id: SourceId,
    revision: SourceRevision,
    units: ModelUnits,
    primitive: SourceGeometry,
}

pub enum SourceGeometry {
    AxisAlignedBlock { min: Point3, max: Point3 },
    FiniteCylinder { frame: Frame3, radius: PositiveLength,
                     axial: OrderedLengthInterval },
    Extrusion { sketch: SourceSketch, frame: Frame3,
                distance: PositiveLength },
    // ImportedBrepRef is reserved and must return Unsupported until implemented.
}

pub struct SourceSnapshot {
    sources: NonEmptyMap<SourceId, SourceSpec>,
    snapshot_stamp: SourceStamp,
}

pub struct SourceRef { id: SourceId, revision: SourceRevision }

pub enum BooleanKind { Union, Intersection, Difference }
pub struct BooleanIntent {
    kind: BooleanKind,
    left: SourceRef,
    right: SourceRef,
    tolerance: GeometryTolerance,
}
\`\`\`

\`SourceSpec\` is immutable after validation; changing a parameter creates a new \`SourceRevision\`. Source IDs are stable across edits. \`SourceRef\` must resolve to an existing **exact** version in the same \`SourceSnapshot\`, left/right references. A Boolean may intentionally use the **same** versioned source on both sides (e.g. A−A); such identities must be handled explicitly rather than rejected solely because both operand IDs match. Units must be equal or converted through an explicit, recorded unit transform; implicit scale conversions are forbidden. Finite/positive dimensions, right-handed frames and legal tolerance budgets are required. General imported B-rep sources remain a separate future scope.

\`SourceStamp\` is a canonical, domain-separated digest of validated source IDs, revisions, numerical parameters, unit policy and BooleanIntent. Canonicalization rules are a **prerequisite** to relying on a digest: reject NaN/Infinity, normalize negative zero, assign a format version, encode numeric fields explicitly, and sort map entries by stable ID. A hash is **not** a geometric equivalence test. Until canonical digest v1 has been implemented and verified by golden vectors, commit must compare immutable source snapshots/revisions directly and revalidate; it must not rely on formatted \`Debug\` strings. The stamp is not an authentication signature.

**Constructor boundary:** \`rustsolid-construct\` takes a \`SourceSpec\` and produces an immutable validated source \`BrepModel\`. The separate \`rustsolid-validation\` crate calculates source-volume envelopes from \`SourceSpec\` **without invoking construct**. This is the central fix for the G8b dependency cycle.

## 2. Contract B — \`BrepModel\`: canonical analytic geometry and topology

**Owner:** \`rustsolid-topology\`. No Boolean, source, document, meshing or validation dependency. Pure geometry catalog values are from \`rustsolid-geometry\`; local catalog handles (\`Curve3Id\`, \`Curve2Id\`, \`SurfaceId\`) are owned by topology.

\`\`\`rust
pub struct BrepModel {
    identity: LocalBrepIncarnation,
    units: ModelUnits,
    tolerance: GeometryTolerance,
    vertices: DenseArena<Vertex>,
    edges: DenseArena<Edge>,
    coedges: DenseArena<Coedge>,
    loops: DenseArena<Loop>,
    faces: DenseArena<Face>,
    shells: DenseArena<Shell>,
    carriers: GeometryCatalog, // Curve3, Curve2 per-coedge, Surface3
}

pub struct Edge { vertices: [VertexId; 2], curve: Curve3Id,
                  interval: ParameterRange, uses: [CoedgeId; 2] }

pub struct Coedge { edge: EdgeId, face: FaceId, loop_id: LoopId,
                    pcurve: Curve2Id, uv_interval: ParameterRange,
                    edge_to_uv_parameter: EdgeToPcurveParameterMap,
                    reversed: bool, next: CoedgeId,
                    prev: CoedgeId, twin: CoedgeId }

pub struct Face { surface: SurfaceId, loops: NonEmptyVec<LoopId>,
                  orientation: FaceSense }

pub struct Shell { faces: NonEmptyVec<FaceId>, sense: ShellSense }

pub enum BuildBrep { Draft(BrepDraft), Checked(BrepModel) }
\`\`\`

Only \`BrepDraft\` is mutable during G4 assembly. \`BrepModel\` is published after supported manifold/geometry checks. The exact representation supports a closed circular edge with equal start/end vertices, opposite coedge uses on the *same face* for a periodic seam, interior UV loops, and multiple oriented shells/components. No universal genus=0, one-shell, or six-face invariant. Multiple shells require their material nesting/orientation to be proven by validation; otherwise return Held. Non-manifold bodies are outside the initial committed solid domain.

**Cross-domain contract:** For every directed coedge the 3D carrier \`Curve3(t)\` and its attached source-face \`Surface3(Curve2(s))\` describe the **same physical geometry** (with explicitly tracked parameter mapping and allowed error). Arc domains, periodic unwrap, reversed use and source cut provenance are all explicit. Two coedges are physical twins only when the shared edge and both source-face pcurves are source-consistent. An internal UV seam is not an extra physical cut. Near coincidence alone does not establish a physical vertex identity.

**Keep outside canonical \`BrepModel\`:** display \`Mesh\`, mass/centroid, visual AABB cache, \`BrepOrigin\` producer tag, feature history, source command, runtime document token, proof of Boolean semantic correctness. These are derived views, external caches or operation-specific metadata, not core shape ownership. Bounding-box or mass computation may be exposed as separate pure queries and independently validated, but must not silently become trusted stored data. \`FaceSense\` is geometric material orientation; \`FaceRole::TopCap\` is source semantics and moves to source provenance.

A closed oriented manifold and pointwise pcurve checks do **not** prove Boolean material correctness. Those are separate G3/G4c/validation obligations. Local topological entity IDs are revision-/incarnation-scoped; they are not saved as persistent CAD feature identifiers.

## 3. Contract C — \`TopologyLineage\`: persistent naming and remapping

**Owner:** \`rustsolid-document\`. It consumes physical source-provenance evidence from G1/G4 and result-local topology IDs, but \`assembly\`/\`boolean\` may not depend on document. Therefore assembly exports a **\`LineageEvidence\`** value containing only source references and local result references. Document converts it into a persistent \`TopologyLineage\` during a successful atomic commit.

\`\`\`rust
pub struct FeatureId(Uuid128);
pub struct PersistentEntityRef { feature: FeatureId, name: SemanticTopoName }
pub struct LocalTopoRef { incarnation: LocalBrepIncarnation, entity: TopoEntityId }

pub struct LineageEvidence {
    source_stamp: SourceStamp,
    result_incarnation: LocalBrepIncarnation,
    links: Vec<PhysicalDerivation>, // source G1 arcs, faces and original uses
}

pub enum LineageRelation {
    Preserved { from: PersistentEntityRef, to: PersistentEntityRef },
    SplitFrom { parent: PersistentEntityRef, children: NonEmptyVec<PersistentEntityRef> },
    MergedFrom { parents: NonEmptyVec<PersistentEntityRef>, child: PersistentEntityRef },
    GeneratedFrom { sources: NonEmptyVec<SourcePhysicalRef>, child: PersistentEntityRef },
    Deleted { former: PersistentEntityRef },
    Unresolved { former: PersistentEntityRef, candidates: Vec<PersistentEntityRef>,
                 reason: AmbiguityReason },
}
pub struct TopologyLineage {
    from_revision: DocumentRevision, to_revision: DocumentRevision,
    relations: Vec<LineageRelation>,
}
\`\`\`

One resulting physical topological vertex has one result-local ID; all relevant original face and G1 event sources must point to it through proven geometry/topology relationships. A physical split may be 1→many, merging may be many→1, and new cut edges often originate from **both operands**. This is **not** a bijection. Proximity-only welding and silent guessing of persistent names are forbidden.

On an ambiguous remap the engine reports \`Unresolved\`. The document may still commit an otherwise valid standalone shape **only if no live dependent feature, selection constraint or user-authorized reference requires the ambiguous name**. If a dependent feature would silently rebind, commit is Held and source state is unchanged. Retired entity names are tombstoned (with provenance), not recycled as other physical entities. Undo/Redo generates a new revision and new local handles; it may restore a semantic name only after a fresh lineage/replay check.

**Naming guarantee:** deterministic source features and proven lineage can retain *semantic* references; RustSolid must not promise universal persistent naming across arbitrary topology changes.

## 4. Contract D — \`CertifiedCandidate\`: immutable proof-bound staged result

**Owner:** \`rustsolid-boolean\`. \`CertifiedCandidate\` cannot be directly constructed or deserialized by clients. A producer may construct a \`StagedCandidate\`, but only \`rustsolid-validation\` can create the accepted validation evidence. \`boolean::prepare\` combines both through an internal gated constructor.

\`\`\`rust
pub enum PrepareOutcome {
    Ready(CertifiedCandidate),
    Held(HoldReason),       // geometric/numeric ambiguity, insufficient evidence
    Unsupported(UnsupportedReason),
}

pub enum CandidateGeometry {
    Empty(DisjointnessOrRemovalProof),
    Solids(BrepModel),      // closed, oriented, supported analytic B-rep
}

pub struct CertifiedCandidate {             // fields private; no Deserialize
    geometry: CandidateGeometry,
    source_stamp: SourceStamp,
    source_versions: NonEmptyVec<SourceRef>,
    intended_operation: BooleanIntent,
    validated_evidence: AcceptedValidation,
    lineage_evidence: LineageEvidence,
    preparation_anchor: SnapshotAnchor,
}

pub struct SnapshotAnchor {                  // lower-level, document supplied
    expected_revision: u64,
    document_incarnation: u128,
    immutable_source_stamp: SourceStamp,
}

pub struct AcceptedValidation {
    b_rep_invariants: ClosedManifoldEvidence,
    geometric_binding: AnalyticCarrierEvidence,
    complete_material_classification: MaterialWitness,
    source_face_coverage: ProvenanceCoverage,
    volume: IndependentVolumeComparison,
    assurance: AssuranceLevel,
}
\`\`\`

A candidate is only \`Ready\` when every *supported-domain* gate passes. Source parameters and operation stamps must match, every result face must be source-owned or proved generated, material-side classification must be complete, B-rep incidence/surface charts must be valid, and independent source-volume and **actual boundary-flux** intervals must agree within explicit error budgets. A shell can be topologically valid but be the *wrong body*; volume agreement alone does not certify Boolean correctness.

**Validation has two independent paths**, which share source/types but not candidate-building routines:

- **Source path:** compute the mathematically expected retained volume enclosure from immutable original \`SourceSpec\` inputs and Boolean operation (G8b family), with explicit domain, interval and rounding-assurance status. Do not consume the staged model or a G4/G8 candidate plan.
- **Actual-result path:** integrate the oriented boundary of the resulting \`BrepModel\` using its actual coedge 3D/UV carriers (G8c family). Do not consume the source path's intermediate integrals or trust a stored candidate mass.

The **comparison gate** may combine these independent reports; neither evaluator calls the other. Preserve an explicit difference between \`GuardedF64\` with conservative engineering guards and \`CertifiedInterval\` requiring a proven directed-rounding/error-bound implementation. Current G8 numerical witnesses are in the \`GuardedF64\` category; do not relabel them as formally certified intervals. A domain-scoped \`AcceptancePolicy\` decides which level is permitted (never by silent downgrade).

**Empty result:** \`CandidateGeometry::Empty\` is permitted only with a separately verified empty-material/complete-removal witness; no fake empty shell or fabricated zero-volume B-rep.

**Strong atomic commit:** \`rustsolid-document::commit(candidate)\` must check document incarnation, current revision, current source versions and current source stamp, consume the candidate once, and **revalidate against immutable current inputs before publication**. Generate lineage mappings, allocate registry/Undo space, and verify independent evidence **before** the publication barrier. A failure (including stale proposal, changed source, insufficient evidence, ambiguous depended-on topology) leaves canonical model, revisions and Undo/Redo state unchanged. Undo/Redo issue new revisions and invalidate stale local handles. Candidate bytes and numerical certificates are *not* trusted upon project reopen: archive persists immutable sources/feature DAG and selected evidence stamps, then replays and recertifies.

Avoid \`can_apply: bool\`, nullable domain states, publicly writable candidate fields and \`Debug\`-string equality as authority. Represent G7/G8 outcomes as exhaustive enum variants. The existing G8c shell-only result remains \`Held\` until a separate boundary-flux comparison and the chosen scoped acceptance policy succeed.

## 5. End-to-end protocol and failure behavior

\`\`\`text
Document source snapshot (document rev N, source stamps)
   └── source::SourceSnapshot [immutable; SourceSpec]
         ├── construct::build_sources(...) -> immutable validated input Breps
         ├── validation::source_volume(SourceSpec + BooleanIntent)
         └── boolean::prepare(intent, snapshot, anchor)
               ├── intersect::untrimmed SSI
               ├── imprint::trimmed G1 graph + physical provenance
               ├── uv::exact-in-form chart subdivision
               ├── classification::complete G2/G3 material cells
               ├── assembly::G4 provisional Brep + lineage evidence
               ├── topology::validate candidate analytic B-rep
               └── validation::actual_boundary_flux(candidate model)
                       + validation::compare_source_vs_result
                         └── PrepareOutcome::Ready(CertifiedCandidate)
Document::commit(expected revision, Ready candidate)
   ├── context/source stamps and no-dependent-ambiguity checks
   ├── source + B-rep + flux + semantic provenance revalidation
   ├── prepare persistent topology lineage & Undo/Redo checkpoint
   └── atomic publication (only infallible swaps after gate)
Reopen: archive source/feature replay -> recertify -> new local handles
\`\`\`

**Fail closed:** unresolved SSI/UV chart, non-manifold shell, incompatible units, source/3D/UV disagreement, unproved source cut provenance, absence of volume enclosure/flux, broken interval comparison, weaker-than-policy numerical assurance, changed source/revision, foreign document anchor or dependent-on-ambiguous name ⇒ \`Held\`/error **without document mutation**.

## 6. Contract tests and migration gates

\`tests/validate_kernel_contracts.py\` and \`kernel-contract-scenarios-v1.json\` verify the **architectural declaration**, reference-case state gates, owner crate DAG and rejection behavior using synthetic values. These checks **do not execute 3D intersections or prove numeric correctness**. Real G7/G8 regression testing remains mandatory before migrating code.

The minimum future implementation suite must include: G7 single-wall, G8 oblique cylinder notch, multiple closed shells, one-seam cylindrical face, circular edges with equal endpoints, hole loops, source-owned 3D/pcurve pair, invalid UV mismatch, almost tangent and overlapping curves, duplicate physical vertex source, source-version replay, cross-document candidate rejection, stale revision, incorrect face sense, independent volume/flux mismatch, ambiguous dependent lineage, closed-looking but wrong material face selection, and failed-commit unchanged revisions/Undo.

**Implementation order after review:** \`source\` types → canonical topology draft/check boundary → G1 physical provenance / UV source-neutral interface → independent validation adapters → \`boolean::prepare\` sealed candidate → \`document::commit\` and lineage → archive/replay and Gefest integration. No new crates are created by this RFC.

## 7. Explicit review decisions and deferred choices

Agreed for v1 design: immutable source versions; one canonical analytic B-rep; no mesh/mass in canonical graph; separate source-only and actual-result oracles; document-only persistent lineage; source-bound, private, one-shot certified candidate; fail-closed atomic commit; nullable state replacement with discriminated enum outcomes.

Deferred for implementation: exact binary canonical hashing layout/golden vectors, NURBS-level interval certification, general non-manifold topology, robust generalized shell nesting for imported solids, precise serialization across new archive versions, policy to accept weaker \`GuardedF64\` evidence outside pinned G7/G8 families, and legal/reference-backend compatibility. These deferrals **must not be silently resolved** by default \`f64\` epsilon or a re-export of old APIs.

## Implementation note — initial source slice

The new `rustsolid-source` package implements validated immutable AxisAlignedBlock and FiniteCylinder sources, explicit monotonic versions, functional snapshots and canonical *full byte* SourceStamp values. The stamp is collision-free by complete encoding, not a cryptographic digest; a future compact digest requires separate golden-vector verification. This first slice does not yet implement Extrusion, source-based builders, GeometryKind extensibility or document commit. Existing G7/G8 operations still run through the transitional model pending adoption.
