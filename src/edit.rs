//! Atomic edits of the oriented polyhedral B-rep.
//!
//! Edits apply to the borrowed body behind an exclusive transaction.
//! Operation-local inverse records restore touched objects on drop or error;
//! no full Solid clone is created. Original IDs survive local splits.
use crate::{
    Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, GeometryError, Loop, LoopId,
    LoopRole, Point3, Solid, TopologyEntity, Vertex, VertexId, NameChange, JournalStats,
};
use crate::journal::UndoFrame;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeSplit {
    pub original_edge: EdgeId,
    pub inserted_vertex: VertexId,
    pub created_edge: EdgeId,
    pub created_coedges: [CoedgeId; 2],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceSplit {
    pub original_face: FaceId,
    pub created_face: FaceId,
    pub created_loop: LoopId,
    pub diagonal_edge: EdgeId,
    pub created_coedges: [CoedgeId; 2],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditDelta {
    SplitEdge(EdgeSplit),
    SplitFace(FaceSplit),
    KillEdgeVertex(crate::KillEdgeVertex),
    KillEdgeFace(crate::KillEdgeFace),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EditReport {
    pub revision_before: u64,
    pub revision_after: u64,
    /// Append-only provenance for created objects and retained parents.
    pub changes: Vec<EditDelta>,
    pub named_changes: Vec<NameChange>,
    pub journal: JournalStats,
}

/// Exclusive edit session backed by local inverse journal frames.
pub struct EditTransaction<'a> {
    target: &'a mut Solid,
    frames: Vec<UndoFrame>,
    changes: Vec<EditDelta>,
    named_changes: Vec<NameChange>,
    failed: bool,
    committed: bool,
}

#[derive(Debug)]
pub(crate) enum EditMemento {
    /// Incremental inverse frames for append-only MEV/MEF commands.
    Delta { frames: Vec<UndoFrame>, body_token: u64 },
    /// Whole-body fallback for KEV/KEF operations which compact indices.
    Snapshot(Box<Solid>),
}
impl EditMemento {
    pub(crate) fn snapshot(before: Solid) -> Self { Self::Snapshot(Box::new(before)) }
    pub(crate) fn check_target(&self, solid: &Solid) -> Result<(), GeometryError> {
        let current = solid.topology_handle(TopologyEntity::Body)
            .ok_or_else(|| edit_error("missing body identity for Memento"))?;
        let origin = match self {
            Self::Delta { body_token, .. } => *body_token,
            Self::Snapshot(body) => body.topology_handle(TopologyEntity::Body)
                .ok_or_else(|| edit_error("snapshot body identity missing"))?.body_token,
        };
        if origin != current.body_token {
            return Err(edit_error("Memento belongs to another body incarnation"));
        }
        Ok(())
    }
    pub(crate) fn restore(self, solid: &mut Solid) -> Result<(), GeometryError> {
        self.check_target(solid)?;
        solid.validate()?;
        let revision = solid.revision.checked_add(1)
            .ok_or_else(|| edit_error("revision overflow"))?;
        match self {
            Self::Delta { frames, .. } => {
                for frame in frames.into_iter().rev() { frame.restore(solid); }
                solid.validate()?;
                solid.revision = revision;
            }
            Self::Snapshot(mut before) => {
                before.revision = revision;
                before.validate()?;
                *solid = *before;
            }
        }
        Ok(())
    }
}

impl Solid {
    pub fn begin_edit(&mut self) -> Result<EditTransaction<'_>, GeometryError> {
        self.validate()?;
        Ok(EditTransaction { target: self, frames: Vec::new(), changes: Vec::new(),
            named_changes: Vec::new(), failed: false, committed: false })
    }
    pub fn edit_atomic(
        &mut self,
        f: impl FnOnce(&mut EditTransaction<'_>) -> Result<(), GeometryError>,
    ) -> Result<EditReport, GeometryError> {
        let mut transaction = self.begin_edit()?;
        f(&mut transaction)?;
        transaction.commit()
    }

    /// Composite-style *ownership* traversal. Shared edges and vertices are
    /// referenced by coedges, not owned by each face (they are graph nodes).
    /// Only call on a validated body.
    pub fn topology_children(&self, entity: TopologyEntity) -> Vec<TopologyEntity> {
        match entity {
            TopologyEntity::Body => self.shells.iter().map(|s| TopologyEntity::Shell(s.id)).collect(),
            TopologyEntity::Shell(id) => self.shells.get(id.0 as usize).filter(|s| s.id == id)
                .map(|s| s.faces.iter().copied().map(TopologyEntity::Face).collect()).unwrap_or_default(),
            TopologyEntity::Face(id) => self.faces.get(id.0 as usize).filter(|f| f.id == id)
                .map(|f| f.loops.iter().copied().map(TopologyEntity::Loop).collect()).unwrap_or_default(),
            TopologyEntity::Loop(id) => {
                let Some(first) = self.loops.get(id.0 as usize).filter(|l| l.id == id)
                    .map(|l| l.first_coedge) else { return Vec::new(); };
                let mut ids = Vec::new();
                let mut current = first;
                for _ in 0..self.coedges.len() {
                    let Some(coedge) = self.coedges.get(current.0 as usize) else { return Vec::new(); };
                    if coedge.id != current || coedge.loop_id != id { return Vec::new(); }
                    ids.push(TopologyEntity::Coedge(current));
                    current = coedge.next;
                    if current == first { return ids; }
                }
                Vec::new()
            }
            // A coedge's edge, twin and vertices are references, not owned children.
            TopologyEntity::Coedge(_) | TopologyEntity::Edge(_) | TopologyEntity::Vertex(_) => Vec::new(),
        }
    }
}

impl EditTransaction<'_> {
    /// Returns the state behind the exclusive transaction borrow.
    pub fn preview(&self) -> &Solid { &*self.target }

    pub fn journal_stats(&self) -> JournalStats {
        let mut stats = JournalStats::default();
        for frame in &self.frames {
            let item = frame.stats();
            stats.frames += item.frames;
            stats.topology_snapshots += item.topology_snapshots;
            stats.triangle_snapshots += item.triangle_snapshots;
        }
        stats
    }

    pub fn split_edge(&mut self, edge: EdgeId, fraction: f64) -> Result<EdgeSplit, GeometryError> {
        if self.failed { return Err(edit_error("transaction already aborted")); }
        let frame = match UndoFrame::for_edge(self.target, edge) {
            Ok(frame) => frame,
            Err(err) => { self.failed = true; return Err(err); }
        };
        self.frames.push(frame);
        let result = split_edge_local(self.target, edge, fraction).and_then(|change| {
            let names = self.target.identities.register_delta(EditDelta::SplitEdge(change))?;
            self.target.validate()?;
            Ok((change, names))
        });
        match result {
            Ok((change,names)) => {
                self.changes.push(EditDelta::SplitEdge(change));
                self.named_changes.extend(names);
                Ok(change)
            }
            Err(err) => { self.failed = true; Err(err) }
        }
    }

    pub fn split_face(&mut self, face: FaceId, start: VertexId, end: VertexId) -> Result<FaceSplit, GeometryError> {
        if self.failed { return Err(edit_error("transaction already aborted")); }
        let frame = match UndoFrame::for_face(self.target, face) {
            Ok(frame) => frame,
            Err(err) => { self.failed = true; return Err(err); }
        };
        self.frames.push(frame);
        let result = split_face_local(self.target, face, start, end).and_then(|change| {
            let names = self.target.identities.register_delta(EditDelta::SplitFace(change))?;
            self.target.validate()?;
            Ok((change, names))
        });
        match result {
            Ok((change,names)) => {
                self.changes.push(EditDelta::SplitFace(change));
                self.named_changes.extend(names);
                Ok(change)
            }
            Err(err) => { self.failed = true; Err(err) }
        }
    }

    pub fn commit(self) -> Result<EditReport, GeometryError> {
        self.commit_inner(false).map(|(report, _)| report)
    }

    pub(crate) fn commit_recorded(self) -> Result<(EditReport, EditMemento), GeometryError> {
        let (report, memento) = self.commit_inner(true)?;
        let memento = memento.ok_or_else(|| edit_error("Memento recording failed"))?;
        Ok((report, memento))
    }

    fn commit_inner(mut self, record_inverse: bool)
        -> Result<(EditReport, Option<EditMemento>), GeometryError>
    {
        if self.failed { return Err(edit_error("cannot commit an aborted transaction")); }
        self.target.validate()?;
        let before = self.target.revision;
        let after = if self.changes.is_empty() { before }
            else { before.checked_add(1).ok_or_else(|| edit_error("revision overflow"))? };
        let journal = self.journal_stats();
        let token = if record_inverse {
            Some(self.target.topology_handle(TopologyEntity::Body)
                .ok_or_else(|| edit_error("body has no identity"))?.body_token)
        } else { None };
        self.target.revision = after;
        let memento = token.map(|body_token| EditMemento::Delta {
            frames: std::mem::take(&mut self.frames), body_token,
        });
        self.committed = true;
        Ok((EditReport { revision_before: before, revision_after: after,
            changes: std::mem::take(&mut self.changes),
            named_changes: std::mem::take(&mut self.named_changes), journal }, memento))
    }

}
impl Drop for EditTransaction<'_> {
    fn drop(&mut self) {
        if !self.committed {
            while let Some(frame) = self.frames.pop() { frame.restore(self.target); }
        }
    }
}

fn edit_error(message: impl Into<String>) -> GeometryError { GeometryError::InvalidEdit(message.into()) }

fn checked_id(count: usize) -> Result<u32, GeometryError> {
    u32::try_from(count).map_err(|_| edit_error("topology handle space exhausted"))
}

fn vertex_distance(a: Point3, b: Point3) -> f64 { a.sub(b).length_squared().sqrt() }

fn split_edge_local(body: &mut Solid, id: EdgeId, fraction: f64) -> Result<EdgeSplit, GeometryError> {
    let edge = body.edges.get(id.0 as usize).filter(|e| e.id == id)
        .ok_or_else(|| edit_error(format!("unknown edge {id:?}")))?.clone();
    if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 {
        return Err(edit_error("edge fraction must be finite and strictly between zero and one"));
    }
    let start = body.vertices[edge.start.0 as usize].position;
    let end = body.vertices[edge.end.0 as usize].position;
    let length = vertex_distance(start, end);
    let inserted = start.add(end.sub(start).scale(fraction));
    if !inserted.is_finite() || vertex_distance(start, inserted) <= body.tolerance.length_at(length)
        || vertex_distance(inserted, end) <= body.tolerance.length_at(length) {
        return Err(edit_error("split point is indistinguishable from an edge endpoint"));
    }
    let new_vertex = VertexId(checked_id(body.vertices.len())?);
    let new_edge = EdgeId(checked_id(body.edges.len())?);
    let new_coedge0 = CoedgeId(checked_id(body.coedges.len())?);
    let new_coedge1 = CoedgeId(checked_id(body.coedges.len() + 1)?);

    // Find the two triangles adjacent to the original *topological* edge and
    // subdivide each without changing winding or face ownership.
    // Modify only incident triangles. Inverse journal snapshots their slots.
    let mut split_count = 0;
    let original_count = body.mesh.triangles.len();
    for index in 0..original_count {
        let tri = body.mesh.triangles[index];
        let face = body.mesh.triangle_faces[index];
        let mut directed = None;
        for [u, v, w] in [[tri[0],tri[1],tri[2]], [tri[1],tri[2],tri[0]], [tri[2],tri[0],tri[1]]] {
            if (u == edge.start.0 && v == edge.end.0) || (u == edge.end.0 && v == edge.start.0) {
                directed = Some([u,v,w]);
                break;
            }
        }
        if let Some([u,v,w]) = directed {
            split_count += 1;
            body.mesh.triangles[index] = [u, new_vertex.0, w];
            body.mesh.triangles.push([new_vertex.0, v, w]);
            body.mesh.triangle_faces.push(face);
        }
    }
    if split_count != 2 { return Err(edit_error("edge must belong to exactly two mesh boundary triangles")); }
    body.vertices.push(Vertex { id: new_vertex, position: inserted });
    body.mesh.vertices.push(inserted);

    // The old edge ID stays on the segment adjacent to its old low-ID vertex.
    // The new edge is canonically (old_high_ID -> new_vertex_ID).
    let (positive, negative) = {
        let a = edge.coedges[0]; let b = edge.coedges[1];
        if body.coedges[a.0 as usize].reversed { (b, a) } else { (a, b) }
    };
    let p = body.coedges[positive.0 as usize].clone();
    let n = body.coedges[negative.0 as usize].clone();
    if p.edge != id || n.edge != id || p.twin != n.id || n.twin != p.id {
        return Err(edit_error("edge incidence is inconsistent before split"));
    }
    body.edges[id.0 as usize].end = new_vertex;
    body.edges[id.0 as usize].coedges = [positive, new_coedge1];
    body.edges.push(Edge {
        id: new_edge, start: edge.end, end: new_vertex, coedges: [new_coedge0, negative],
    });

    // Old coedges retain their IDs; each new coedge is inserted after its old
    // neighbour in the same directed loop. Their twin pairs are reconnected.
    body.coedges[positive.0 as usize].edge = id;
    body.coedges[positive.0 as usize].reversed = false;
    body.coedges[positive.0 as usize].next = new_coedge0;
    body.coedges[positive.0 as usize].twin = new_coedge1;
    body.coedges[negative.0 as usize].edge = new_edge;
    body.coedges[negative.0 as usize].reversed = false;
    body.coedges[negative.0 as usize].next = new_coedge1;
    body.coedges[negative.0 as usize].twin = new_coedge0;
    body.coedges[p.next.0 as usize].prev = new_coedge0;
    body.coedges[n.next.0 as usize].prev = new_coedge1;
    body.coedges.push(Coedge {
        id: new_coedge0, edge: new_edge, face: p.face, loop_id: p.loop_id,
        reversed: true, next: p.next, prev: positive, twin: negative,
    });
    body.coedges.push(Coedge {
        id: new_coedge1, edge: id, face: n.face, loop_id: n.loop_id,
        reversed: true, next: n.next, prev: negative, twin: positive,
    });
    for c in [&p, &n] {
        let ring = &mut body.faces[c.face.0 as usize].boundary;
        let original_start = if c.reversed { edge.end } else { edge.start };
        let original_end = if c.reversed { edge.start } else { edge.end };
        let where_at = (0..ring.len()).find(|&i| ring[i] == original_start && ring[(i+1)%ring.len()] == original_end)
            .ok_or_else(|| edit_error("cannot find split edge in owning face boundary"))?;
        ring.insert(where_at + 1, new_vertex);
    }
    Ok(EdgeSplit { original_edge: id, inserted_vertex: new_vertex,
        created_edge: new_edge, created_coedges: [new_coedge0, new_coedge1] })
}

#[derive(Clone, Copy)]
struct P2 { u: f64, v: f64 }
fn cross(a: P2, b: P2, c: P2) -> f64 {
    (b.u-a.u)*(c.v-a.v) - (b.v-a.v)*(c.u-a.u)
}
fn on_segment(a: P2, b: P2, p: P2, length_eps: f64, area_eps: f64) -> bool {
    cross(a,b,p).abs() <= area_eps && p.u >= a.u.min(b.u)-length_eps
        && p.u <= a.u.max(b.u)+length_eps && p.v >= a.v.min(b.v)-length_eps
        && p.v <= a.v.max(b.v)+length_eps
}
fn segments_intersect(a: P2, b: P2, c: P2, d: P2, length_eps:f64, area_eps:f64) -> bool {
    let (c1,c2,c3,c4)=(cross(a,b,c),cross(a,b,d),cross(c,d,a),cross(c,d,b));
    ((c1 > area_eps && c2 < -area_eps)||(c1 < -area_eps && c2 > area_eps))
        && ((c3 > area_eps && c4 < -area_eps)||(c3 < -area_eps && c4 > area_eps))
        || on_segment(a,b,c,length_eps,area_eps) || on_segment(a,b,d,length_eps,area_eps)
        || on_segment(c,d,a,length_eps,area_eps) || on_segment(c,d,b,length_eps,area_eps)
}

fn project_ring(body: &Solid, ring: &[VertexId]) -> Result<Vec<P2>, GeometryError> {
    if ring.len() < 3 { return Err(edit_error("face has less than three boundary vertices")); }
    let anchor=body.vertices[ring[0].0 as usize].position;
    let mut normal=Point3{x:0.0,y:0.0,z:0.0};
    for i in 0..ring.len() {
        let a=body.vertices[ring[i].0 as usize].position.sub(anchor);
        let b=body.vertices[ring[(i+1)%ring.len()].0 as usize].position.sub(anchor);
        normal=normal.add(a.cross(b));
    }
    let drop = if normal.x.abs() >= normal.y.abs() && normal.x.abs() >= normal.z.abs() {0}
        else if normal.y.abs() >= normal.z.abs() {1} else {2};
    Ok(ring.iter().map(|id| {
        let p=body.vertices[id.0 as usize].position.sub(anchor);
        match drop {0=>P2{u:p.y,v:p.z},1=>P2{u:p.x,v:p.z},_=>P2{u:p.x,v:p.y}}
    }).collect())
}
fn polygon_extent(points:&[P2]) -> f64 {
    let (mut min_u, mut max_u, mut min_v, mut max_v) =
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for p in points { min_u=min_u.min(p.u); max_u=max_u.max(p.u);
        min_v=min_v.min(p.v); max_v=max_v.max(p.v); }
    (max_u-min_u).max(max_v-min_v)
}
fn signed_area(points:&[P2]) -> f64 {
    let p0=points[0];
    (0..points.len()).map(|i| cross(p0,points[i],points[(i+1)%points.len()])).sum::<f64>()*0.5
}
fn lerp(a:P2,b:P2,t:f64)->P2{P2{u:a.u+(b.u-a.u)*t,v:a.v+(b.v-a.v)*t}}
fn strictly_inside(p:P2,points:&[P2],linear:f64,area:f64)->bool{
    let mut inside=false;
    for i in 0..points.len() {
        let a=points[i];let b=points[(i+1)%points.len()];
        if on_segment(a,b,p,linear,area) { return false; }
        if (a.v>p.v)!=(b.v>p.v) && p.u<a.u+(p.v-a.v)*(b.u-a.u)/(b.v-a.v) {inside=!inside;}
    }
    inside
}
fn diagonal_is_inside(body:&Solid,points:&[P2],i:usize,j:usize)->bool{
    let extent=polygon_extent(points);
    let linear=body.tolerance.length_at(extent);
    let area=body.tolerance.area_at(extent);
    let a=points[i];let b=points[j];
    if ((a.u-b.u).powi(2)+(a.v-b.v).powi(2)).sqrt() <= linear {return false;}
    for k in 0..points.len() {
        let next=(k+1)%points.len();
        if k==i || next==i || k==j || next==j {continue;}
        if segments_intersect(a,b,points[k],points[next],linear,area) {return false;}
    }
    for t in [0.0001, 0.5, 0.9999] {
        if !strictly_inside(lerp(a,b,t),points,linear,area) {return false;}
    }
    true
}
pub(crate) fn triangulate_ring(body:&Solid,ring:&[VertexId])->Result<Vec<[u32;3]>,GeometryError>{
    let pts=project_ring(body,ring)?;
    let area=signed_area(&pts);
    let eps=body.tolerance.area_at(polygon_extent(&pts));
    if !area.is_finite() || area.abs() <= eps { return Err(edit_error("split face has degenerate area")); }
    let sign=area.signum();
    let mut remain: Vec<usize> =(0..ring.len()).collect();
    let mut triangles=Vec::with_capacity(ring.len()-2);
    while remain.len()>3 {
        let mut found=None;
        for k in 0..remain.len() {
            let ai=remain[(k+remain.len()-1)%remain.len()];
            let bi=remain[k];
            let ci=remain[(k+1)%remain.len()];
            let (a,b,c)=(pts[ai],pts[bi],pts[ci]);
            if cross(a,b,c)*sign <= eps {continue;}
            // Collinear vertices on the *ear boundary* are allowed so long as
            // all triangles remain nondegenerate and no vertex is swallowed.
            if remain.iter().any(|&t|t!=ai && t!=bi && t!=ci &&
                cross(a,b,pts[t])*sign>=-eps && cross(b,c,pts[t])*sign>=-eps
                    && cross(c,a,pts[t])*sign>=-eps) {continue;}
            found=Some((k,ai,bi,ci));break;
        }
        let Some((k,ai,bi,ci))=found else {return Err(edit_error("face triangulation failed"));};
        triangles.push([ring[ai].0,ring[bi].0,ring[ci].0]);
        remain.remove(k);
    }
    let [a,b,c]:[usize;3]=remain.try_into().map_err(|_|edit_error("invalid final face triangle"))?;
    if cross(pts[a],pts[b],pts[c])*sign <= eps {return Err(edit_error("degenerate final face triangle"));}
    triangles.push([ring[a].0,ring[b].0,ring[c].0]);
    Ok(triangles)
}

fn split_face_local(body:&mut Solid,face:FaceId,start:VertexId,end:VertexId)->Result<FaceSplit,GeometryError>{
    let original=body.faces.get(face.0 as usize).filter(|f|f.id==face)
        .ok_or_else(||edit_error(format!("unknown face {face:?}")))?.clone();
    let n=original.boundary.len();
    let i=original.boundary.iter().position(|&x|x==start)
        .ok_or_else(||edit_error("start vertex is not on the face boundary"))?;
    let j=original.boundary.iter().position(|&x|x==end)
        .ok_or_else(||edit_error("end vertex is not on the face boundary"))?;
    if i==j || (i+1)%n==j || (j+1)%n==i {
        return Err(edit_error("split requires distinct nonadjacent face vertices"));
    }
    let (i,j)=if i<j{(i,j)}else{(j,i)};
    if !diagonal_is_inside(body,&project_ring(body,&original.boundary)?,i,j) {
        return Err(edit_error("face diagonal lies outside its interior or touches the boundary"));
    }
    let old_ring=original.boundary[i..=j].to_vec();
    let new_ring=original.boundary[j..].iter().chain(original.boundary[..=i].iter()).copied().collect::<Vec<_>>();
    // Certify both planar face regions *before* changing topology.
    let old_tris=triangulate_ring(body,&old_ring)?;
    let new_tris=triangulate_ring(body,&new_ring)?;
    let new_face=FaceId(checked_id(body.faces.len())?);
    let new_loop=LoopId(checked_id(body.loops.len())?);
    let new_edge=EdgeId(checked_id(body.edges.len())?);
    let coedge_old=CoedgeId(checked_id(body.coedges.len())?);
    let coedge_new=CoedgeId(checked_id(body.coedges.len()+1)?);
    let existing_loop=original.loops[0];
    let first=body.loops[existing_loop.0 as usize].first_coedge;
    let mut ordered=Vec::with_capacity(n);
    let mut current=first;
    for _ in 0..n {
        ordered.push(current);
        current=body.coedges[current.0 as usize].next;
    }
    // For each face boundary vertex, retrieve its outgoing coedge. Both arrays
    // are cyclic; the loop's first coedge need not be at boundary index zero.
    let mut at=vec![CoedgeId(u32::MAX);n];
    for id in ordered {
        let coedge=&body.coedges[id.0 as usize];
        let directed_start=coedge.start_vertex(&body.edges[coedge.edge.0 as usize]);
        let k=original.boundary.iter().position(|&v|v==directed_start)
            .ok_or_else(||edit_error("coedge is missing from face boundary"))?;
        at[k]=id;
    }
    if at.iter().any(|id|id.0==u32::MAX) {return Err(edit_error("face coedge lookup incomplete"));}
    let first_old=at[i];let last_old=at[j-1];
    let first_new=at[j];let last_new=at[(i+n-1)%n];
    for k in 0..n {
        if k>=j || k<i {let c=&mut body.coedges[at[k].0 as usize];c.face=new_face;c.loop_id=new_loop;}
    }
    body.coedges[last_old.0 as usize].next=coedge_old;
    body.coedges[first_old.0 as usize].prev=coedge_old;
    body.coedges[last_new.0 as usize].next=coedge_new;
    body.coedges[first_new.0 as usize].prev=coedge_new;
    let u=original.boundary[i];let v=original.boundary[j];
    let (low,high) = if u<v {(u,v)}else{(v,u)};
    let old_reverse= v != low;
    body.edges.push(Edge {id:new_edge,start:low,end:high,coedges:[coedge_old,coedge_new]});
    body.coedges.push(Coedge {id:coedge_old,edge:new_edge,face,loop_id:existing_loop,
        reversed:old_reverse,next:first_old,prev:last_old,twin:coedge_new});
    body.coedges.push(Coedge {id:coedge_new,edge:new_edge,face:new_face,loop_id:new_loop,
        reversed:!old_reverse,next:first_new,prev:last_new,twin:coedge_old});
    body.loops[existing_loop.0 as usize].first_coedge=first_old;
    body.loops.push(Loop{id:new_loop,face:new_face,role:LoopRole::Outer,first_coedge:first_new});
    body.faces[face.0 as usize].boundary=old_ring;
    body.faces.push(Face{id:new_face,role:original.role,boundary:new_ring,loops:vec![new_loop]});
    body.shells[0].faces.push(new_face);

    // The two resulting polygons have as many total triangles as the original.
    // Keep unaffected triangle indices stable.
    let old_indices = body.mesh.triangle_faces.iter().enumerate()
        .filter(|(_, &owner)| owner == face)
        .map(|(i, _)| i).collect::<Vec<_>>();
    let replacements = old_tris.into_iter().map(|tri| (tri,face))
        .chain(new_tris.into_iter().map(|tri| (tri,new_face))).collect::<Vec<_>>();
    if old_indices.len() != replacements.len() {
        return Err(edit_error("face split triangulation changed total triangle count"));
    }
    for (index, (tri,owner)) in old_indices.into_iter().zip(replacements) {
        body.mesh.triangles[index] = tri;
        body.mesh.triangle_faces[index] = owner;
    }
    Ok(FaceSplit {original_face:face,created_face:new_face,created_loop:new_loop,
        diagonal_edge:new_edge,created_coedges:[coedge_old,coedge_new]})
}
