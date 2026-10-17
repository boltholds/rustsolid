# RustSolid v0.41: Gefest G5b source-certified persistent document bridge

Scope: one canonical AABB minus one upright full-height analytic cylinder,
where the Boolean leaves four independently closed corner solids.
This is a versioned and validated **RustSolid/geometry.v1 document archive
contract** for Gefest. UI integration with the separate Gefest front-end,
cross-process locks and durable database transactions are not claimed.

## Source-first archive

G5ProjectArchive (schema_version=1) writes only original parametric
source descriptors, immutable operation DAG, output topology names, an
independent G4c result/volume stamp per accepted history state, the
monotonic document revision and all previous/next Undo/Redo snapshots.
No render mesh, mutable B-rep graph, or process-local handle is stored.

DAG: source_a (source index 0), source_b (index 1), difference (inputs
source_a and source_b). Its output_bindings associate all 177 logical
TopologyName values with the difference operation. When the current
state is undone, bindings remain in the graph if an accepted result is
reachable through Redo or past history. Names use the owning feature_key
plus signed box corner and semantic topological role.

On open, every state is reconstructed independently:
- recreate exact bounded primitive inputs from validated dimensions;
- rebuild source G1/G2/G4b2b3b analytic intersections and shell;
- independently re-evaluate the G4c four-component source quadrature
  and B-rep boundary flux;
- compare all stored volume stamps against freshly rederived proof;
- verify timeline structural rules, source ownership, namespace,
  topology names and the canonical feature DAG;
- allocate a fresh document token and handle generations;
- restore the current document revision and both history branches.

The file is not a trusted instruction. Unsupported schema versions,
unexpected JSON fields, missing graph nodes/edges, stale or tampered
result stamps, invalid primitive dimensions and missing history states
are rejected before the restored document is exposed.

Bounded policy: schema 1 only; <= 64 saved Undo/Redo states;
<= 1 MiB JSON; unsupported original body types fail closed. The
current live G5BooleanDocument has no imposed 64-state limit, but its
G5b exporter refuses to persist histories beyond the safe bound.

## Serialization and file writes

G5ProjectArchive::from_document(&doc)
G5ProjectArchive::to_json()
G5ProjectArchive::from_json(&json)
G5ProjectArchive::reopen()
G5ProjectArchive::save_atomic(&path)
G5ProjectArchive::load_file(&path)

save_atomic revalidates the full document/graph/history before any
filesystem action; writes a temporary file in the same directory,
flushes the file, then renames to the final target. On a recoverable
pre-rename failure the previous target remains unchanged and the
temporary file is cleaned up. Atomic replacement and folder fsync
semantics are OS/filesystem dependent. It is not a distributed
transaction or a guarantee against power failure on all platforms.
load_file uses a bounded read and the same independent revalidation.

## geometry.v1 adapter

New CLI JSON operation g5_project accepts either
feature_key + left/right BrepInput with actions, or a previously
returned archive and actions. Actions: commit, undo, redo, reopen,
replace_source(first|second). A single call returns the updated
self-contained archive, graph, semantic_names, document_revision,
volume/topology and transition events. Future Gefest clients can
persist the returned archive with the save API, then load/reopen it in
a different process.

Sample first call:

{
  "operation":"g5_project",
  "feature_key":"assembly_a",
  "left":{"kind":"block","origin":[-1.5,1,-1.5],
          "width":3,"height":3,"depth":3},
  "right":{"kind":"analytic_cylinder","origin":[0,0,0],
           "radius":2,"height":5},
  "actions":[{"kind":"commit"}]
}

Next request: pass the returned data.archive into
{"operation":"g5_project","archive":...,
 "actions":[{"kind":"undo"},{"kind":"reopen"},{"kind":"redo"}]}

Any failed action returns the existing fail-closed geometry.v1 response
with ok=false, data=null, solid=null. No partial project state is
returned or written to a destination path on command failure.

## Tests

cargo test --test g5_project_archive
cargo test --test g5_project_cli_json
cargo test --all-targets

Coverage: source-certified real file roundtrip, existing
four-shell B-rep, graph linkage across Undo/Redo, fresh handles,
stable 177 names, asymmetric source edit, tamper rejection for
schema/graph/volume/history, invalid write destination preserving the
previous file and CLI separate-process archive handoff.

Outstanding: Gefest front-end integration, application feature graph
migration across kernel versions, storage-level concurrency/CAS,
arbitrary curved bodies, industrial crash-consistency/ACID and
unbounded history. PR remains draft until reviewed.
