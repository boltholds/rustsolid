# geometry.v1 physical-edge picking witnesses

The bounded oblique-cylinder clip CLI publishes native B-rep edge witnesses
under `solid.analytic_brep.pick_edges` alongside its existing
`solid.triangle_face_ids` output.

Each `pick_edges` entry contains:

- `edge_id`: source `BrepModel::edges[].id`, never a display mesh edge index.
- `curve_kind`: `line`, `circle`, or `ellipse`, from the native `Curve3` carrier.
- `face_ids`: the two actual adjacent coedge face IDs. Periodic seams can
  legitimately have the same face ID twice.
- `start_vertex_id`, `end_vertex_id`: native topological vertex IDs.
- `points`: source `Curve3::evaluate()` samples over the bounded native
  `EdgeGeometry.domain`. Lines use 2 endpoints, curved edges use between
  13 and 129 samples, including both ends.

The sampling is for viewport hit testing only. It does not replace analytic
geometry or change the B-rep. In particular, no triangulation interior
diagonal is ever granted an Edge ID. Invalid/unbounded native curve domains
or missing geometry fail the response; no invented boundary is provided.

`triangle_face_ids[i]` continues to map indexed display triangle `i` to
physical `FaceId`. A Gefest viewport can use the indexed mesh raycast
`faceIndex` as a triangle array index to retrieve a source face ID. The
edge selector must raycast the explicit native edge witness polylines.

The CLI is stateless: its `history.revision` describes operations within
**one request**. It is not a durable project revision. Consumers must also
bind their selection to the enclosing project/preview generation and drop
picks after rebuilding or replacing the B-rep.

Validation:

```bash
cargo test --test oblique_cylinder_cli
```

The test verifies all 3 physical edges of each halfspace cut, bounded and
finite samples, source coedge ownership and closed analytical rings, as well
as the existing undo/redo and Held/no-partial-geometry behavior.
