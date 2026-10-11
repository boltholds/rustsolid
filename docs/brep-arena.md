# RustSolid v0.17 — sparse generational B-rep identity arena

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

A batch can reference pre-existing handles. Referencing a vertex created
**earlier in the same batch** by an output placeholder is not yet available;
such workflows require a later dependency-aware command representation.

## Boundary with v0.16 TopologyDeltaPlan

`TopologyDeltaPlan` remains `can_apply=false` until region classification and
topology-general curve/subsurface editing exist. `BrepArena` explicitly denies
`require_executable_plan` even if clients tamper with `can_apply`. It cannot
silently turn the proposed analytic curves/periodic anchors into polygonal
edge splits. This v0.17 commit supports actual Euler subdivisions of **already
supported planar B-rep**, not arbitrary Boolean construction.

## First-party acceptance and integration

[tests/fixtures/arena_transactions.v1.json](../tests/fixtures/arena_transactions.v1.json)
pins ten deterministic scenarios for compound MEV+MEF, atomic rollback,
transaction drop, generational ABA prevention, compacting KEV/KEF identities,
foreign handles, read-only analytic geometry, forged handles and plan refusal.
[tests/brep_arena.rs](../tests/brep_arena.rs) executes them with analytically
known cuboid volume/area/centroid, Euler invariants, Undo/Redo, and immutable
snapshot checks. No other CAD kernel acts as an oracle.

`geometry.v1` exposes a stateless `arena_transaction` workflow using entity
selectors local to the request (not external serializations of process-local
handles). Each `batches` item is one atomic command and Undo step. The result
includes the updated shared B-rep summary, committed source identifiers,
created/retired handles, and revision information. A failure returns an error
without a partial solid. The existing API operations remain unchanged.

## Next gates

1. Introduce a dependency-aware edit batch using newly created vertex/edge
   placeholders without surrendering atomicity or immutable command replay.
2. Move canonical curved and planar geometry into a mutable shared sparse
   B-rep store with explicit curve/trim ownership, local tolerance propagation
   and persistent source-to-created lineage.
3. Implement formal Boolean region classification before authorizing v0.16
   `TopologyDeltaPlan` as an executable topology delta.
4. Optimize history duplication into an arena-level transaction journal once
   the first-party semantic and failure-injection regressions remain green.
