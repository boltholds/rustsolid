# RustSolid v0.18 — dependency-aware generational B-rep arena

This is the **first editable arena boundary**, not a claim that the curved
analytic B-rep is already a fully writable or generalized Boolean kernel.
`BrepArena` owns a validated `CommandHistory` for the existing *planar polyhedral*
body, or an independently validated *read-only* analytic cylinder. Both can be
read through the shared `BrepModel`, using immutable revision-tagged
`ArenaSnapshot`s.

## Canonical entity identity

Each live topological entity has an `ArenaHandle { arena_token, slot,
generation }` and one deterministic `TopologyName`. The arena maintains
separate `by_name` and `by_entity` indexes over **sparse slots**. A deleted
entity retires its slot and increments the generation before reallocation,
preventing ABA when numerical IDs are reused. A surviving semantic name keeps
its arena handle even when a compacting inverse Euler operation reassigns its
internal `VertexId`, `FaceId`, or `EdgeId`. Foreign handles, forged generations,
and deleted handles are rejected. This generation is *process-local*, not a
persistent naming protocol or cross-regeneration shape matching.

## Strong atomic commands and replay

`BrepArena::begin_transaction()` stages type-checked, handle-based commands:

- `split_edge(edge, fraction)` — MEV, finite strictly interior fraction;
- `split_face(face, start_vertex, end_vertex)` — MEF;
- `kill_edge_vertex(vertex)` — KEV of a valence-two collinear split;
- `kill_edge_face(edge, removed_face)` — KEF of two coplanar neighbors.

The editor translates handles into **semantic topology names** and runs the
whole batch using the existing `CommandHistory`/`EditMemento` implementation.
**Neither the body nor sparse identities change while the command is staged.**
For this first implementation, the history and slot registry are cloned and
edited privately. Geometry validation and registry reconciliation complete
*before* the candidate is published. Any failure (including a later command
in a batch) drops the candidate and keeps exactly the previous model, revision,
undo/redo state and every arena handle. A successful command batch advances
one model revision and produces **one Undo/Redo step**. This is strong atomicity,
but copying the working history is not yet an O(changed entities) edit algorithm.

`ArenaCommit` returns the existing `EditReport` (including source entity
provenance), `created` and `retired` named generational bindings, and `remapped`
survivors whose compact dense indices changed. Failed and dropped
transactions produce no `ArenaCommit`. Undo/Redo retain or invalidate arena
handles according to semantic survival, and recreated elements get different
generations. Immutable snapshots remain readable as records of their original
revision but `check_snapshot` rejects them against newer arena revisions.

## Dependency-aware batches (v0.18)

The existing handle-only API is source compatible. New `*_ref` methods return
`ArenaStep`; `step.output(ArenaOutput::...)` produces a typed
`ArenaEntityRef` usable by strictly later commands of **the same**
transaction. Currently exposed outputs are `InsertedVertex` and
`RemainderEdge` from MEV, and `CreatedFace` and `DiagonalEdge` from
MEF. A produced reference is never a process-level `ArenaHandle`.

```rust
use rustsolid::{ArenaOutput, BrepArena, EdgeId, TopologyEntity};
// arena: &mut BrepArena, with a validated polyhedral body
let old_edge = arena.handle(TopologyEntity::Edge(EdgeId(0))).unwrap();
let mut tx = arena.begin_transaction()?;
let first = tx.split_edge_ref(old_edge, 0.25)?;
tx.split_edge_ref(first.output(ArenaOutput::RemainderEdge), 0.5)?;
let committed = tx.commit()?; // exactly one revision / Undo step
```

References are checked for transaction identity, operation order and output
role during staging. Any rejected stage permanently aborts that transaction,
even if the caller ignores the error. At commit, dependent batches execute
a private sequential Euler preflight and resolve created elements by the
**actual emitted `EditDelta` and `TopologyName`**, including after index
compaction. A retired output, invalid geometry, missing name, or divergent
preflight cancels the entire batch without changing the arena, revision,
handles or history. The compiled semantic-name commands execute as **one**
replayable `CommandHistory` batch, so Undo and Redo retain correct dependencies.

Dependent preflight presently performs extra work: one private Solid for
sequential dependency discovery, then the atomic history candidate. Pure
handle-only batches avoid that extra preflight. No claiming O(changed-elements)
geometry editing yet.

## Boundary with v0.16 TopologyDeltaPlan

`TopologyDeltaPlan` remains `can_apply=false` until region classification and
topology-general curve/subsurface editing exist. `BrepArena` explicitly denies
`require_executable_plan` even if clients tamper with `can_apply`. It cannot
silently turn the proposed analytic curves/periodic anchors into polygonal
edge splits. This v0.18 slice supports actual Euler subdivisions of **already
supported planar B-rep**, not arbitrary Boolean construction.

## First-party acceptance and integration

[tests/fixtures/arena_transactions.v1.json](../tests/fixtures/arena_transactions.v1.json)
pins nineteen deterministic scenarios, including compound MEV+MEF, late failure
rollback, transaction drop, ABA protection, compacting KEV/KEF, foreign handles,
read-only analytic geometry, forged handles, plan refusal, dependent edge/face
chains, same-batch inverse, cross-transaction outputs, wrong output roles, references across compaction,
retired outputs and same-batch MEF/KEF inverse.
[tests/brep_arena.rs](../tests/brep_arena.rs) executes them with analytically
known cuboid volume/area/centroid, Euler invariants, Undo/Redo, and immutable
snapshot checks. No other CAD kernel acts as an oracle.

`geometry.v1` exposes a stateless `arena_transaction` workflow using entity
selectors local to the request (not external serializations of process-local
handles). A batch is one atomic command group and Undo step. Integer `*_id`
selectors refer only to the pre-batch model. To target a created output in
the same batch, supply `*_ref` with `{"step":0,"output":"remainder_edge"}`
instead of the corresponding `*_id` selector. The step is zero-based, must
precede its use, and output roles are type-checked. For example:

```json
{"operation":"arena_transaction",
 "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
 "batches":[[
   {"kind":"split_edge","edge_id":0,"fraction":0.25},
   {"kind":"split_edge","edge_ref":{"step":0,"output":"remainder_edge"},"fraction":0.5}
 ]],"undo":1,"redo":1}
```

The result
includes the updated shared B-rep summary, committed source identifiers,
created/retired handles, and revision information. A failure returns an error
without a partial solid. The existing API operations remain unchanged.

## Next gates

1. Move canonical curved and planar geometry into a mutable shared sparse
   B-rep store with explicit curve/trim ownership, local tolerance propagation
   and persistent source-to-created lineage.
2. Implement formal Boolean region classification before authorizing v0.16
   `TopologyDeltaPlan` as an executable topology delta.
3. Optimize history duplication and dependent preflight into an arena-level
   transaction journal once stronger semantic/failure-injection regressions
   remain green.

### Rejected staging is an abort, not a partial batch

A staging error (including wrong entity kind, foreign/stale handle, or invalid
edge fraction) **poisons** the transaction. Subsequent staging and `commit()`
fail even when the caller ignores the earlier `Err`. Geometry, sparse slots,
generations, revision, Undo/Redo stack and the snapshot remain untouched.
See `staging_failure_poison_is_atomic_even_when_caller_ignores_error`.
