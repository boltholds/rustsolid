# RustSolid v0.41 G5a — source-bound atomic four-shell Boolean document

G5a is a **bounded in-process document executor** for the G4b2b3b
four-corner Box-minus-upright-Cylinder result independently certified
by G4c1. It does not replace BrepArena's existing polyhedral editor,
does not provide a general curved CSG transaction, and does not enable
disk persistence or application-level database atomicity.

## Owner and staged proposal

G5BooleanDocument owns two private BrepBody inputs, optional committed
analytic four-shell B-rep, a monotonic document revision, one process-local
document token, a semantic name/handle registry, and independent Undo/Redo
stacks. Its feature_key is an explicit stable document-owned identifier.

prepare() rebuilds G4b2b3b and G4c1 from the document's current input
geometry. An unsupported input returns Held without changing anything.
Successful preparation returns a G5PreparedBoolean containing the complete
originally certified proposal and volume evidence, document token,
expected revision, feature identity, and exact private in-process source
representation stamps. Source stamps are not cryptographic hashes,
serialization formats, or persistent cache keys.

commit(proposal) first checks document ownership, document revision,
absence of an already committed Boolean and both complete source stamps.
While holding the exclusive Rust &mut self owner, it performs a fresh
G4b2b3b full source-graph reconstruction and a fresh G4c independent
volume verification on the actual current inputs and the proposed result.
It then validates the B-rep topology, produces every semantic name and
generational handle, clones a coherent undo checkpoint, reserves enough
undo-stack space, and computes the next monotonic revision.

Only AFTER every fallible check passes does a short non-fallible
publication barrier swap the new document state, identity registry
and revision. A recoverable Error cannot publish a partial graph.
A concurrent application must serialize access through an exclusive
document owner or lock; database transactions and process crashes are
outside this guarantee.

## Identity

Each result entity receives a collision-checked logical TopologyName
derived from the document feature_key, Boolean operation, one of the
four signed corners, and its semantic role:

- vertex: bottom_corner, bottom_x_arc, bottom_z_arc, top equivalents;
- edge: planar cap lines, circular trim edges, vertical generators;
- coedge: each edge role plus the distinct owning face role;
- face, loop, shell, and body get equally scoped names.

All 177 result entities are covered (1 body, 24 vertices, 36 edges,
72 coedges, 20 loops, 20 faces, 4 shells). The names are independent
of the current dense array index and can serve as logical document
references for this bounded producer. This does not claim geometric
matching for arbitrary regenerated shapes.

G5TopologyHandle = {document_token, slot, generation} is process-local.
It resolves only while its exact incarnation is active; an Undo, Redo
or source edit retires old handles. Restored output gets new generations,
while logically identical result roles keep their TopologyName. This
prevents ABA where a numeric EdgeId or slot is reused for a different
material edge. Rebind via resolve_name() after regeneration.

## Document history

- A successful Boolean commit adds one complete prior-document undo
  checkpoint and clears Redo.
- replace_source() performs a separately atomic source mutation,
  retires the result and live handles, and increments revision even if
  the caller substitutes an identical source (stale preview rejection).
- Undo and Redo validate the saved source and complete G4c certificate
  before restoring a coherent state. Each transition increments the
  monotonic revision. No old generational handle is reactivated.
- Undo followed by a new commit discards the obsolete Redo branch.
- Empty history returns Ok(None) without changing anything.
- A malformed preview, foreign document token, stale revision, invalid
  source or geometry change returns Err/Held and preserves the document,
  result, history depths, identities and revision.

## Acceptance and limits

G5a permits mutation only after independent G4c volume_certified and
G4b2b3b closed_shell_proven pass. Both predecessor plans retain their
can_apply=false preview flags, since those generic plans do not publish
documents; only the explicitly bounded G5BooleanDocument owns this
narrower commit contract.

Supported: validated AABB minus one full-height upright cylinder
crossing all four walls and splitting the box into four disconnected,
genus-zero closed analytic components. No arbitrary multi-face curved
arrangement, multi-cylinder edits or NURBS kernels.

Not yet implemented: serialization with durable topological history,
distributed compare-and-swap, persistent database commits, replay across
different RustSolid versions, a generalized curved editable BrepArena,
or integration with the Gefest CAD document/feature graph.

Run: cargo test --test atomic_boolean_document
Then: cargo test --all-targets

Tests cover committed four-shell B-rep and all 177 semantic names,
successful Undo/Redo, stale handle invalidation, source edits and
revision conflicts (including identical-geometry ABA replacement),
tampered G4/G4c evidence, rejecting damaged inputs without changes,
foreign document tokens, Redo branching, and unsupported fail-closed
preview requests.
