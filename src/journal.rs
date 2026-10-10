//! Operation-local inverse journal. Captures *only* modified existing objects,
//! appended vector lengths and affected mesh triangles. A failed or dropped
//! transaction replays frames backwards; it never clones an entire `Solid`.
use crate::{Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, GeometryError, Loop,
    LoopId, Shell, Solid};
use std::collections::BTreeSet;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct JournalStats {
    pub frames: usize,
    /// Number of existing topology objects copied into the inverse journal.
    pub topology_snapshots: usize,
    /// Number of triangle + owner pairs copied; independent of total mesh size.
    pub triangle_snapshots: usize,
}

#[derive(Debug, Clone, Copy)]
struct Lengths {
    vertices: usize, edges: usize, coedges: usize, loops: usize,
    faces: usize, shells: usize, mesh_vertices: usize, triangles: usize,
    triangle_faces: usize,
}
impl Lengths {
    fn of(body: &Solid) -> Self {
        Self { vertices: body.vertices.len(), edges: body.edges.len(),
            coedges: body.coedges.len(), loops: body.loops.len(),
            faces: body.faces.len(), shells: body.shells.len(),
            mesh_vertices: body.mesh.vertices.len(), triangles: body.mesh.triangles.len(),
            triangle_faces: body.mesh.triangle_faces.len() }
    }
}

#[derive(Debug)]
pub(crate) struct UndoFrame {
    len: Lengths,
    identity_checkpoint: (usize, u64),
    edges: Vec<(usize, Edge)>,
    coedges: Vec<(usize, Coedge)>,
    loops: Vec<(usize, Loop)>,
    faces: Vec<(usize, Face)>,
    shells: Vec<(usize, Shell)>,
    triangles: Vec<(usize, [u32;3], FaceId)>,
}
impl UndoFrame {
    fn empty(body: &Solid) -> Self {
        Self { len: Lengths::of(body), identity_checkpoint: body.identities.checkpoint(),
            edges: Vec::new(), coedges: Vec::new(), loops: Vec::new(),
            faces: Vec::new(), shells: Vec::new(), triangles: Vec::new() }
    }

    pub(crate) fn for_edge(body: &Solid, id: EdgeId) -> Result<Self, GeometryError> {
        let edge = body.edges.get(id.0 as usize).filter(|e| e.id == id)
            .ok_or_else(|| GeometryError::InvalidEdit(format!("unknown edge {id:?}")))?;
        let mut out = Self::empty(body);
        out.edges.push((id.0 as usize, edge.clone()));
        let mut touched_coedges = BTreeSet::new();
        let mut touched_faces = BTreeSet::new();
        for id in edge.coedges {
            let coedge = &body.coedges[id.0 as usize];
            touched_coedges.insert(id);
            touched_coedges.insert(coedge.next);
            touched_faces.insert(coedge.face);
        }
        for id in touched_coedges {
            out.coedges.push((id.0 as usize, body.coedges[id.0 as usize].clone()));
        }
        for id in touched_faces {
            out.faces.push((id.0 as usize, body.faces[id.0 as usize].clone()));
        }
        out.triangles = body.mesh.triangles.iter().zip(&body.mesh.triangle_faces)
            .enumerate().filter(|(_, (tri, _))| tri.contains(&edge.start.0) && tri.contains(&edge.end.0))
            .map(|(index, (tri, owner))| (index, *tri, *owner)).collect();
        Ok(out)
    }

    pub(crate) fn for_face(body: &Solid, id: FaceId) -> Result<Self, GeometryError> {
        let face = body.faces.get(id.0 as usize).filter(|f| f.id == id)
            .ok_or_else(|| GeometryError::InvalidEdit(format!("unknown face {id:?}")))?;
        let mut out = Self::empty(body);
        out.faces.push((id.0 as usize, face.clone()));
        for loop_id in &face.loops {
            let loop_data = &body.loops[loop_id.0 as usize];
            out.loops.push((loop_id.0 as usize, loop_data.clone()));
            let first = loop_data.first_coedge;
            let mut current = first;
            for _ in 0..body.coedges.len() {
                let coedge = &body.coedges[current.0 as usize];
                out.coedges.push((current.0 as usize, coedge.clone()));
                current = coedge.next;
                if current == first { break; }
            }
        }
        for shell in &body.shells {
            if shell.faces.contains(&id) { out.shells.push((shell.id.0 as usize, shell.clone())); }
        }
        out.triangles = body.mesh.triangles.iter().zip(&body.mesh.triangle_faces)
            .enumerate().filter(|(_, (_, owner))| **owner == id)
            .map(|(index, (tri, owner))| (index, *tri, *owner)).collect();
        Ok(out)
    }

    pub(crate) fn stats(&self) -> JournalStats {
        JournalStats { frames: 1,
            topology_snapshots: self.edges.len() + self.coedges.len() + self.loops.len()
                + self.faces.len() + self.shells.len(),
            triangle_snapshots: self.triangles.len() }
    }

    pub(crate) fn restore(self, body: &mut Solid) {
        // First discard appended data; then replace precisely the touched
        // preexisting objects and mesh entries in their original slots.
        body.vertices.truncate(self.len.vertices);
        body.edges.truncate(self.len.edges);
        body.coedges.truncate(self.len.coedges);
        body.loops.truncate(self.len.loops);
        body.faces.truncate(self.len.faces);
        body.shells.truncate(self.len.shells);
        body.mesh.vertices.truncate(self.len.mesh_vertices);
        body.mesh.triangles.truncate(self.len.triangles);
        body.mesh.triangle_faces.truncate(self.len.triangle_faces);
        for (i, before) in self.edges { body.edges[i] = before; }
        for (i, before) in self.coedges { body.coedges[i] = before; }
        for (i, before) in self.loops { body.loops[i] = before; }
        for (i, before) in self.faces { body.faces[i] = before; }
        for (i, before) in self.shells { body.shells[i] = before; }
        for (i, tri, owner) in self.triangles {
            body.mesh.triangles[i] = tri;
            body.mesh.triangle_faces[i] = owner;
        }
        body.identities.rewind(self.identity_checkpoint);
    }
}
