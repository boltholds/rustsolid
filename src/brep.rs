//! Explicit oriented B-rep incidence. A coedge is one directed use of an edge
//! by a face loop; `next`/`prev` form its boundary cycle, and `twin` connects
//! the two oppositely directed uses on this closed orientable shell.
use crate::{Edge, EdgeId, Face, FaceId, GeometryError, Point3, Solid, VertexId};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CoedgeId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LoopId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ShellId(pub u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopRole { Outer, Inner }

#[derive(Debug, Clone)]
pub struct Coedge {
    pub id: CoedgeId,
    pub edge: EdgeId,
    pub face: FaceId,
    pub loop_id: LoopId,
    /// False: edge.start -> edge.end; true: edge.end -> edge.start.
    pub reversed: bool,
    pub next: CoedgeId,
    pub prev: CoedgeId,
    pub twin: CoedgeId,
}

impl Coedge {
    pub fn start_vertex(&self, edge: &Edge) -> VertexId {
        if self.reversed { edge.end } else { edge.start }
    }
    pub fn end_vertex(&self, edge: &Edge) -> VertexId {
        if self.reversed { edge.start } else { edge.end }
    }
}

#[derive(Debug, Clone)]
pub struct Loop {
    pub id: LoopId,
    pub face: FaceId,
    pub role: LoopRole,
    pub first_coedge: CoedgeId,
}

#[derive(Debug, Clone)]
pub struct Shell {
    pub id: ShellId,
    pub faces: Vec<FaceId>,
    pub closed: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TopologyEntity {
    Body,
    Vertex(VertexId),
    Edge(EdgeId),
    Coedge(CoedgeId),
    Loop(LoopId),
    Face(FaceId),
    Shell(ShellId),
}

/// First failing invariant, with a stable machine-readable code and a witness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyIssue {
    pub code: &'static str,
    pub entity: TopologyEntity,
    pub detail: String,
}

impl fmt::Display for TopologyIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at {:?}: {}", self.code, self.entity, self.detail)
    }
}

fn issue(code: &'static str, entity: TopologyEntity, detail: impl Into<String>) -> TopologyIssue {
    TopologyIssue { code, entity, detail: detail.into() }
}

/// Build deterministic edge/coedge incidence from each face's oriented legacy
/// vertex boundary. Only one closed outer loop per face is supported for now.
pub(crate) fn from_faces(faces: &mut [Face]) -> Result<(Vec<Edge>, Vec<Coedge>, Vec<Loop>, Vec<Shell>), GeometryError> {
    let mut uses: BTreeMap<(u32, u32), Vec<CoedgeId>> = BTreeMap::new();
    let mut coedges = Vec::new();
    let mut loops = Vec::new();
    for face in faces.iter_mut() {
        if face.boundary.len() < 3 {
            return Err(GeometryError::InvalidTopology(format!("face {:?} has too few vertices", face.id)));
        }
        let loop_id = LoopId(loops.len() as u32);
        let first = CoedgeId(coedges.len() as u32);
        loops.push(Loop { id: loop_id, face: face.id, role: LoopRole::Outer, first_coedge: first });
        face.loops = vec![loop_id];
        let n = face.boundary.len();
        for i in 0..n {
            let start = face.boundary[i];
            let end = face.boundary[(i + 1) % n];
            if start == end { return Err(GeometryError::InvalidTopology(format!("zero-length topological edge in face {:?}", face.id))); }
            let key = (start.0.min(end.0), start.0.max(end.0));
            let id = CoedgeId(coedges.len() as u32);
            coedges.push(Coedge {
                id, edge: EdgeId(u32::MAX), face: face.id, loop_id,
                reversed: start.0 > end.0,
                next: CoedgeId(first.0 + ((i + 1) % n) as u32),
                prev: CoedgeId(first.0 + ((i + n - 1) % n) as u32),
                twin: CoedgeId(u32::MAX),
            });
            uses.entry(key).or_default().push(id);
        }
    }
    let mut edges = Vec::new();
    for ((start, end), incidences) in uses {
        if incidences.len() != 2 {
            return Err(GeometryError::InvalidTopology(format!("edge ({start},{end}) has {} uses instead of two", incidences.len())));
        }
        let (a, b) = (incidences[0].0 as usize, incidences[1].0 as usize);
        if coedges[a].reversed == coedges[b].reversed || coedges[a].face == coedges[b].face {
            return Err(GeometryError::InvalidTopology(format!("edge ({start},{end}) is not shared by opposite face boundaries")));
        }
        let eid = EdgeId(edges.len() as u32);
        edges.push(Edge { id: eid, start: VertexId(start), end: VertexId(end), coedges: [incidences[0], incidences[1]] });
        coedges[a].edge = eid;
        coedges[b].edge = eid;
        coedges[a].twin = incidences[1];
        coedges[b].twin = incidences[0];
    }
    let shell = Shell { id: ShellId(0), faces: faces.iter().map(|f| f.id).collect(), closed: true };
    Ok((edges, coedges, loops, vec![shell]))
}

fn orient2(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

fn on_segment2(a: (f64, f64), p: (f64, f64), b: (f64, f64), length: f64, area: f64) -> bool {
    orient2(a, p, b).abs() <= area
        && p.0 >= a.0.min(b.0) - length && p.0 <= a.0.max(b.0) + length
        && p.1 >= a.1.min(b.1) - length && p.1 <= a.1.max(b.1) + length
}

fn intersect2(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64), length: f64, area: f64) -> bool {
    let (p, q, r, s) = (orient2(a,b,c), orient2(a,b,d), orient2(c,d,a), orient2(c,d,b));
    ((p > area && q < -area) || (p < -area && q > area))
        && ((r > area && s < -area) || (r < -area && s > area))
        || on_segment2(a,c,b,length,area) || on_segment2(a,d,b,length,area)
        || on_segment2(c,a,d,length,area) || on_segment2(c,b,d,length,area)
}

fn face_normal(body: &Solid, ring: &[VertexId]) -> Result<(Point3, f64), TopologyIssue> {
    let anchor = body.vertices[ring[0].0 as usize].position;
    let mut sum = Point3 { x: 0.0, y: 0.0, z: 0.0 };
    let mut min = anchor;
    let mut max = anchor;
    for i in 0..ring.len() {
        let a = body.vertices[ring[i].0 as usize].position;
        let b = body.vertices[ring[(i + 1) % ring.len()].0 as usize].position;
        sum = sum.add(a.sub(anchor).cross(b.sub(anchor)));
        min.x = min.x.min(a.x); min.y = min.y.min(a.y); min.z = min.z.min(a.z);
        max.x = max.x.max(a.x); max.y = max.y.max(a.y); max.z = max.z.max(a.z);
    }
    let extent = max.sub(min).length_squared().sqrt();
    let mag = sum.length_squared().sqrt();
    if !extent.is_finite() || !mag.is_finite() || mag <= 2.0 * body.tolerance.area_at(extent) {
        return Err(issue("face.degenerate", TopologyEntity::Face(FaceId(0)), "oriented face area is below tolerance"));
    }
    let normal = sum.scale(1.0 / mag);
    for &vertex in ring {
        let p = body.vertices[vertex.0 as usize].position;
        if p.sub(anchor).dot(normal).abs() > body.tolerance.length_at(extent) {
            return Err(issue("face.nonplanar", TopologyEntity::Vertex(vertex), "vertex is outside face plane tolerance"));
        }
    }
    // Planar self-crossing and self-touching boundaries are not valid trim loops.
    // Projection uses anchor-relative coordinates to avoid precision loss after
    // translating a small face far away from the world origin.
    let drop_axis = if normal.x.abs() >= normal.y.abs() && normal.x.abs() >= normal.z.abs() { 0 }
        else if normal.y.abs() >= normal.z.abs() { 1 } else { 2 };
    let pts: Vec<(f64,f64)> = ring.iter().map(|&id| {
        let p = body.vertices[id.0 as usize].position.sub(anchor);
        match drop_axis { 0 => (p.y,p.z), 1 => (p.x,p.z), _ => (p.x,p.y) }
    }).collect();
    let n = pts.len();
    for i in 0..n {
        for j in (i+1)..n {
            if j == i + 1 || (i == 0 && j == n - 1) { continue; }
            if intersect2(pts[i], pts[(i+1)%n], pts[j], pts[(j+1)%n],
                body.tolerance.length_at(extent), body.tolerance.area_at(extent)) {
                return Err(issue("face.self_intersection", TopologyEntity::Face(FaceId(0)), "face boundary self-crosses or self-touches"));
            }
        }
    }
    Ok((normal, extent))
}

impl Solid {
    /// Validate the canonical B-rep graph, returning a typed witness on error.
    /// This v0.2 validator deliberately accepts only one connected, orientable,
    /// closed genus-zero shell with one outer loop per face.
    pub fn check_topology(&self) -> Result<(), TopologyIssue> {
        if self.vertices.is_empty() || self.edges.is_empty() || self.faces.is_empty() || self.shells.len() != 1 {
            return Err(issue("body.invalid_counts", TopologyEntity::Body, "expected nonempty body with exactly one shell"));
        }
        if self.vertices.iter().enumerate().any(|(i, v)| v.id.0 as usize != i || !v.position.is_finite()) {
            return Err(issue("vertex.invalid", TopologyEntity::Body, "invalid vertex ID or non-finite position"));
        }
        for (i, edge) in self.edges.iter().enumerate() {
            if edge.id.0 as usize != i || edge.start.0 >= edge.end.0 || edge.end.0 as usize >= self.vertices.len() {
                return Err(issue("edge.invalid", TopologyEntity::Edge(edge.id), "edge has invalid ID or endpoint index"));
            }
        }
        for (i, face) in self.faces.iter().enumerate() {
            if face.id.0 as usize != i || face.loops.len() != 1 || face.boundary.len() < 3 {
                return Err(issue("face.invalid", TopologyEntity::Face(face.id), "only one outer loop with >=3 vertices is supported"));
            }
            if face.boundary.iter().any(|v| v.0 as usize >= self.vertices.len()) {
                return Err(issue("face.vertex_out_of_bounds", TopologyEntity::Face(face.id), "face boundary references missing vertex"));
            }
        }
        let mut usage = vec![0usize; self.coedges.len()];
        let mut owned_loops = vec![0usize; self.loops.len()];
        for (idx, face) in self.faces.iter().enumerate() {
            let lid = face.loops[0];
            if lid.0 as usize >= self.loops.len() {
                return Err(issue("face.loop_out_of_bounds", TopologyEntity::Face(face.id), "unknown loop"));
            }
            owned_loops[lid.0 as usize] += 1;
            if self.loops[lid.0 as usize].face != face.id {
                return Err(issue("face.loop_owner", TopologyEntity::Face(face.id), "loop belongs to a different face"));
            }
            if idx != face.id.0 as usize {
                return Err(issue("face.invalid_id", TopologyEntity::Face(face.id), "face index mismatch"));
            }
        }
        for (i, loop_data) in self.loops.iter().enumerate() {
            if loop_data.id.0 as usize != i || loop_data.face.0 as usize >= self.faces.len() || owned_loops[i] != 1 {
                return Err(issue("loop.invalid_owner", TopologyEntity::Loop(loop_data.id), "orphaned, duplicated or misindexed loop"));
            }
            if loop_data.role != LoopRole::Outer {
                return Err(issue("loop.unsupported_inner", TopologyEntity::Loop(loop_data.id), "inner loops are reserved for a later slice"));
            }
            if loop_data.first_coedge.0 as usize >= self.coedges.len() {
                return Err(issue("loop.invalid_start", TopologyEntity::Loop(loop_data.id), "invalid starting coedge"));
            }
            let mut ring = Vec::new();
            let mut visited = BTreeSet::new();
            let mut current = loop_data.first_coedge;
            loop {
                let i = current.0 as usize;
                if i >= self.coedges.len() || !visited.insert(i) {
                    return Err(issue("loop.broken_cycle", TopologyEntity::Loop(loop_data.id), "cycle escapes range or revisits coedge"));
                }
                let c = &self.coedges[i];
                if c.id != current || c.loop_id != loop_data.id || c.face != loop_data.face || c.edge.0 as usize >= self.edges.len() {
                    return Err(issue("coedge.incorrect_owner", TopologyEntity::Coedge(current), "invalid ID, loop, face or edge link"));
                }
                let next = c.next.0 as usize;
                let prev = c.prev.0 as usize;
                if next >= self.coedges.len() || prev >= self.coedges.len() || self.coedges[next].prev != current || self.coedges[prev].next != current {
                    return Err(issue("coedge.broken_links", TopologyEntity::Coedge(current), "next/prev reciprocity failed"));
                }
                let edge = &self.edges[c.edge.0 as usize];
                let next_c = &self.coedges[next];
                if next_c.edge.0 as usize >= self.edges.len() || c.end_vertex(edge) != next_c.start_vertex(&self.edges[next_c.edge.0 as usize]) {
                    return Err(issue("loop.open_boundary", TopologyEntity::Coedge(current), "consecutive coedges do not share vertex"));
                }
                ring.push(c.start_vertex(edge));
                usage[i] += 1;
                current = c.next;
                if current == loop_data.first_coedge { break; }
                if visited.len() > self.coedges.len() {
                    return Err(issue("loop.broken_cycle", TopologyEntity::Loop(loop_data.id), "cycle does not terminate"));
                }
            }
            if ring.len() < 3 || ring.iter().copied().collect::<BTreeSet<_>>().len() != ring.len() {
                return Err(issue("loop.degenerate", TopologyEntity::Loop(loop_data.id), "loop has <3 or repeated boundary vertices"));
            }
            if ring != self.faces[loop_data.face.0 as usize].boundary {
                return Err(issue("face.boundary_mismatch", TopologyEntity::Face(loop_data.face), "legacy boundary does not match canonical coedge cycle"));
            }
            // Face plane and area checks operate on local coordinate deltas.
            if let Err(mut error) = face_normal(self, &ring) {
                if error.code == "face.degenerate" || error.code == "face.self_intersection" { error.entity = TopologyEntity::Face(loop_data.face); }
                return Err(error);
            }
        }
        if usage.iter().any(|&n| n != 1) {
            return Err(issue("coedge.orphaned_or_reused", TopologyEntity::Body, "each coedge must appear in exactly one loop"));
        }
        let mut edge_use_counts = vec![0usize; self.edges.len()];
        for (i, coedge) in self.coedges.iter().enumerate() {
            if coedge.id.0 as usize != i || coedge.edge.0 as usize >= self.edges.len() || coedge.twin.0 as usize >= self.coedges.len() {
                return Err(issue("coedge.invalid", TopologyEntity::Coedge(coedge.id), "out of bounds coedge ID, edge or twin"));
            }
            edge_use_counts[coedge.edge.0 as usize] += 1;
            let twin = &self.coedges[coedge.twin.0 as usize];
            if twin.twin != coedge.id || twin.edge != coedge.edge || twin.reversed == coedge.reversed || twin.face == coedge.face {
                return Err(issue("coedge.inconsistent_twin", TopologyEntity::Coedge(coedge.id), "twin must reverse the same edge in another face"));
            }
        }
        for edge in &self.edges {
            let [a, b] = edge.coedges;
            if a == b || a.0 as usize >= self.coedges.len() || b.0 as usize >= self.coedges.len() || edge_use_counts[edge.id.0 as usize] != 2 {
                return Err(issue("edge.non_manifold", TopologyEntity::Edge(edge.id), "expected exactly two distinct coedges"));
            }
            let ca = &self.coedges[a.0 as usize];
            let cb = &self.coedges[b.0 as usize];
            if ca.edge != edge.id || cb.edge != edge.id || ca.twin != b || cb.twin != a {
                return Err(issue("edge.incidence_mismatch", TopologyEntity::Edge(edge.id), "edge's coedge pair disagrees with twin links"));
            }
            let pa = self.vertices[edge.start.0 as usize].position;
            let pb = self.vertices[edge.end.0 as usize].position;
            let dx = pb.sub(pa);
            let extent = self.bbox.max.sub(self.bbox.min).length_squared().sqrt();
            if !dx.length_squared().is_finite() || dx.length_squared().sqrt() <= self.tolerance.length_at(extent) {
                return Err(issue("edge.degenerate", TopologyEntity::Edge(edge.id), "edge length below model tolerance"));
            }
        }
        let shell = &self.shells[0];
        if shell.id != ShellId(0) || !shell.closed || shell.faces.len() != self.faces.len() {
            return Err(issue("shell.invalid", TopologyEntity::Shell(shell.id), "expected one closed shell with all faces"));
        }
        let mut face_count = vec![0usize; self.faces.len()];
        for id in &shell.faces {
            if id.0 as usize >= self.faces.len() {
                return Err(issue("shell.face_out_of_bounds", TopologyEntity::Shell(shell.id), "missing face"));
            }
            face_count[id.0 as usize] += 1;
        }
        if face_count.iter().any(|&n| n != 1) {
            return Err(issue("shell.face_ownership", TopologyEntity::Shell(shell.id), "each face must occur exactly once"));
        }
        // Every face must be reachable by crossing shared edges.
        let mut connected = BTreeSet::new();
        let mut pending = vec![shell.faces[0]];
        while let Some(fid) = pending.pop() {
            if !connected.insert(fid) { continue; }
            let lid = self.faces[fid.0 as usize].loops[0];
            let first = self.loops[lid.0 as usize].first_coedge;
            let mut c = first;
            loop {
                let coedge = &self.coedges[c.0 as usize];
                pending.push(self.coedges[coedge.twin.0 as usize].face);
                c = coedge.next;
                if c == first { break; }
            }
        }
        if connected.len() != self.faces.len() {
            return Err(issue("shell.disconnected", TopologyEntity::Shell(shell.id), "shell has disconnected face groups"));
        }
        // Edge-manifold alone does not rule out a pinch at a vertex. Verify
        // that all outgoing coedges at each vertex form a single fan.
        for vertex in &self.vertices {
            let outgoing: BTreeSet<_> = self.coedges.iter().filter(|c| {
                c.start_vertex(&self.edges[c.edge.0 as usize]) == vertex.id
            }).map(|c| c.id).collect();
            if outgoing.is_empty() {
                return Err(issue("vertex.orphaned", TopologyEntity::Vertex(vertex.id), "unused vertex"));
            }
            let first = *outgoing.iter().next().unwrap();
            let mut fan = BTreeSet::new();
            let mut current = first;
            loop {
                if !fan.insert(current) {
                    return Err(issue("vertex.broken_fan", TopologyEntity::Vertex(vertex.id), "vertex adjacency does not form one cycle"));
                }
                let previous = self.coedges[current.0 as usize].prev;
                current = self.coedges[previous.0 as usize].twin;
                if current == first { break; }
                if !outgoing.contains(&current) || fan.len() > outgoing.len() {
                    return Err(issue("vertex.broken_fan", TopologyEntity::Vertex(vertex.id), "vertex fan escapes its incident coedges"));
                }
            }
            if fan.len() != outgoing.len() {
                return Err(issue("vertex.multiple_fans", TopologyEntity::Vertex(vertex.id), "pinched vertex with disjoint face fans"));
            }
        }
        if self.euler_characteristic() != 2 {
            return Err(issue("shell.euler_characteristic", TopologyEntity::Shell(shell.id), "single genus-zero shell must have V-E+F=2"));
        }
        Ok(())
    }
}
