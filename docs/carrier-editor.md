# RustSolid v0.22 — transactional canonical B-rep carrier editor

The first editable layer on the common `BrepModel` owns actual 3D edge
curves (`Curve3`), 2D coedge UV trims (`Curve2`) and face support surfaces
(`Surface3`), with source-validated atomic transactions. It works with
previously validated polyhedral, analytical-cylinder and cylindrical
through-hole B-reps (including the latter's genus-one inner loops).

This release enables **geometry-equivalent carrier and UV-chart editing**.
It does not change the physical boundary, create or retire topological
entities, or recalculate volume after deformation. That constraint is an
explicit correctness gate until generalized region construction, shape
regeneration and independent mass validation are available.

## Why geometry ownership needed to change

Prior `BrepModel` exported separate curve/trim/surface vectors as a
read-only interchange representation. A `CoedgeGeometry` stores an index
into the 2D curve catalog; `EdgeGeometry` stores a 3D curve index, and
`FaceGeometry` references a surface. Those numeric positions can change
when a shared carrier is replaced and unused catalog entries are removed.

`BrepCarrierEditor` now owns an isolated canonical graph and a validated
reference. `CarrierOwner` is a typed, topological reference:

- `CarrierOwner::Edge(EdgeId)` identifies the edge's 3D curve.
- `CarrierOwner::Coedge(CoedgeId)` identifies its own 2D UV trim.
- `CarrierOwner::Face(FaceId)` identifies the supporting surface.

Each has a process-local `CarrierHandle { editor_token, owner, generation }`.
Existing owner's handle remains stable when the internal `Curve3Id`,
`Curve2Id` or `SurfaceId` changes after compaction. Foreign editor tokens,
forged generations, and nonexistent owners are rejected.

Carrier replacement is **copy-on-write**, including if two coedges happen to
reference the same 2D trim curve. At commit, old unreferenced catalog entries
are discarded and the remaining dense IDs are consistently reindexed;
references to the same still-live original carrier remain shared. No
tombstone reclamation of topological owner slots is yet required because this
slice does not insert/delete owners; the generation registry establishes the
API boundary for that later phase.

## Supported edits

All edits are staged within `BrepCarrierTransaction`, using the handle of
the topological owner. A batch is a single revision and a single Undo step.

- `replace_edge(handle, Curve3, ParameterRange)`: replace one edge's 3D
  analytic carrier, with a new oriented bounded parameter interval.
- `replace_trim(handle, Curve2, ParameterRange)`: replace a single
  coedge's UV pcurve and its parameter interval without modifying the
  opposite coedge.
- `replace_surface(handle, Surface3)`: replace the analytic face surface
  within its current family. Multiple staged replacements may be required
  before the full candidate becomes internally consistent.
- `reparameterize_edge(handle, scale, offset)`: preserve the 3D edge path
  under `t = scale * u + offset`. Lines support positive affine mappings;
  circles support an angular phase shift with `scale=1`.
- `reparameterize_trim(handle, scale, offset)`: preserve one 2D trim
  trajectory. Lines support positive affine mappings; the current
  `Circle2` representation supports full-turn parameter shifts only.
- `rebase_face_uv(handle, du, dv)`: move a planar face's UV frame origin
  in-plane, or a cylinder's UV frame origin axially. All of the face's
  coedge pcurves and parameter bounds change together, preserving their
  physical 3D paths. Periodic angular `u` rebasing is explicitly refused
  until a full seam-chart transformation is available.

These operations work on existing circular and linear carriers and
on the cylindrical side of the v0.21 through-hole. `SphereSurface`
reparameterization is not implemented.

## Atomicity and numerical proof

The editor clones the canonical `BrepModel` into an isolated candidate. No
public references expose mutable access to either the current model or the
candidate's geometry arrays. On any stage failure the transaction becomes
poisoned: ignoring `Err` cannot cause its earlier operations to commit.

Successful candidate publication requires:

1. Complete copy-on-write catalog compaction and reassignment of every
   referenced carrier ID.
2. `BrepModel::validate()`, including coedge orientation, UV-to-3D
   agreement, manifold incidence, support bindings, periodic seam semantics
   and the genus-one inner-loop checks.
3. Confirmation that the same physical surface support/orientation survives
   each face edit. A cylinder must retain radius and axis; a plane must
   retain its oriented supporting plane.
4. Re-evaluation of **nine independent fractions** of every 3D edge and
   corresponding source-surface UV trim, comparing the existing and proposed
   world-space paths within bounded model-space numerical resolution.
5. Unchanged topology counts, vertex locations, bounding box, model units,
   tolerance, and analytic mass properties — compared both against the
   prior revision and the original canonical source.

A rejected commit drops its entire candidate. The original geometry,
carrier IDs, sparse owner generations, revision and Undo/Redo branches do
not change. Undo and Redo restore validated full-model Mementos with
**monotonically increasing runtime revisions**; stale snapshot views are
rejected. `MAX_EDITS=256` and Undo history is bounded at 128 steps.

The geometric equivalence comparison samples known supported analytic
line/circle/plane/cylinder families. It is a tolerance-aware correctness
check for the restricted edit contract, **not** a symbolic proof of shape
equivalence for arbitrary splines or a general-purpose CAD deformation
guarantee. Physical shape-changing edits must wait for an independent
mass/boundary recomputation contract.

## Rust API

```rust
use rustsolid::{BrepBody, BrepCarrierEditor, CarrierOwner,
    EdgeId, FaceId, Point3, Solid};
let solid=Solid::block(Point3{x:0.0,y:0.0,z:0.0},2.0,3.0,4.0)?;
let mut editor=BrepCarrierEditor::from_body(&BrepBody::from(solid))?;
let edge=editor.handle(CarrierOwner::Edge(EdgeId(0))).unwrap();
let face=editor.handle(CarrierOwner::Face(FaceId(0))).unwrap();
let mut transaction=editor.begin_transaction();
transaction.reparameterize_edge(edge,2.0,0.25)?;
transaction.rebase_face_uv(face,0.5,-0.25)?;
transaction.commit()?;
editor.validate()?;
assert_eq!(editor.model().mass.volume,24.0);
```

`BrepArena::detached_carrier_editor()` creates an independent editor from
an immutable snapshot of any existing arena. Its edits **do not** mutate the
source arena's Euler or parametric history; the detached editor explicitly
owns the resulting common B-rep geometry. General in-place unification of
Euler commands, curved parameter updates and carrier edits is the next gate.

## geometry.v1: carrier_transaction

```json
{
  "operation": "carrier_transaction",
  "source": {"kind":"analytic_through_hole","origin":[-3,-1,-3],
             "width":6,"height":4,"depth":6,
             "center":[0,0],"radius":1},
  "batches": [[
    {"kind":"reparameterize_edge","edge_id":12,"scale":1,"offset":0.5},
    {"kind":"rebase_face_uv","face_id":6,"du":0,"dv":0.75}
  ]],
  "undo": 1,
  "redo": 1
}
```

Each batch is one atomic revision/Undo step. `data.kind` equals
`carrier_transaction`. The response includes the revision, topology/genus,
analytic mass, new curve catalog counts and owner/carrier rebindings.
`shape_preserving:true` and `topology_mutated:false` are explicit.
No partial solid is emitted on failure. External JSON selectors are scoped
to that request's initially constructed source; process-local handles must
not be serialized between requests.

The CLI supports reparameterizing edge and coedge carriers and rebasing
face UV frames. Direct curve/trim/surface substitution requires the Rust
library's typed API and full commit validation.

## Acceptance and next gates

`tests/carrier_editor.rs` pins nine original regressions plus shared trim
copy-on-write, explicit replacement, repeated-owner budget exhaustion and
detached arena integration. CI runs both this integration test and an
independent JSON smoke/failure gate, in addition to the entire prior suite.

Next, a fully **structural** heterogeneous B-rep editor must support:
new and retired 3D curves, 2D trims, topological vertices/edges/coedges,
face loops and shell components; update of volume and surface moments;
all-pairs crossing/overlap region construction; and failure-injected
Boolean commits with persistent lineage. Only then can the existing
`TopologyDeltaPlan` become an executable Boolean transaction.

**Unchanged limitation:** `TopologyDeltaPlan::can_apply == false`.
Neither the v0.21 specialized through-hole builder nor this editor makes
general curved Boolean operations available.
