//! Oriented B-rep incidence data, stable within a single model revision.
//! Pure topology data cannot depend on Booleans, meshing, or document storage.
use rustsolid_foundation::{GeometryError, Point3};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VertexId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FaceId(pub u32);

#[derive(Debug, Clone)]
pub struct Vertex {
    pub id: VertexId,
    pub position: Point3,
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub id: EdgeId,
    pub start: VertexId,
    pub end: VertexId,
    /// The two opposing directed usages of this edge on the closed shell.
    pub coedges: [CoedgeId; 2],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceRole { BottomCap, TopCap, Side(u32) }

#[derive(Debug, Clone)]
pub struct Face {
    pub id: FaceId,
    pub role: FaceRole,
    /// Boundary vertices ordered counterclockwise as viewed from the outside.
    pub boundary: Vec<VertexId>,
    /// Canonical face loops. `boundary` remains a checked compatibility view.
    pub loops: Vec<LoopId>,
}

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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reversed_coedge_uses_opposite_vertex() {
        let e = Edge { id: EdgeId(2), start: VertexId(3), end: VertexId(4),
            coedges: [CoedgeId(9), CoedgeId(10)] };
        let c = Coedge { id: CoedgeId(10), edge: e.id, face: FaceId(0),
            loop_id: LoopId(0), reversed: true, next: CoedgeId(10),
            prev: CoedgeId(10), twin: CoedgeId(9) };
        assert_eq!(c.start_vertex(&e), VertexId(4));
        assert_eq!(c.end_vertex(&e), VertexId(3));
        let owner = TopologyEntity::Coedge(c.id);
        assert_eq!(owner, TopologyEntity::Coedge(CoedgeId(10)));
    }
}

/// Assembles the closed oriented manifold incidence on a scratch copy.
/// On any error, the input face records and their loops remain unchanged.
pub fn from_faces(faces: &mut [Face]) -> Result<(Vec<Edge>, Vec<Coedge>, Vec<Loop>, Vec<Shell>), GeometryError> {
    let mut scratch = faces.to_vec();
    let result = build_on_scratch(&mut scratch)?;
    faces.clone_from_slice(&scratch);
    Ok(result)
}

/// Build deterministic edge/coedge incidence from each face's oriented legacy
/// vertex boundary. Only one closed outer loop per face is supported for now.
fn build_on_scratch(faces: &mut [Face]) -> Result<(Vec<Edge>, Vec<Coedge>, Vec<Loop>, Vec<Shell>), GeometryError> {
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


#[cfg(test)]
mod incidence_tests {
    use super::*;
    fn tetrahedron() -> Vec<Face> {
        let rings = [[0,2,1], [0,1,3], [1,2,3], [2,0,3]];
        rings.into_iter().enumerate().map(|(i,ring)| Face {
            id: FaceId(i as u32),
            role: FaceRole::Side(i as u32),
            boundary: ring.into_iter().map(VertexId).collect(),
            loops: vec![],
        }).collect()
    }
    #[test]
    fn tetrahedron_has_symmetric_twin_and_cycle_incidence() {
        let mut faces=tetrahedron();
        let (edges, coedges, loops, shells)=from_faces(&mut faces).unwrap();
        assert_eq!((edges.len(), coedges.len(), loops.len(), shells.len()), (6, 12, 4, 1));
        assert!(shells[0].closed);
        assert!(faces.iter().all(|face| face.loops.len()==1));
        for edge in &edges {
            assert_ne!(edge.coedges[0],edge.coedges[1]);
        }
        for c in &coedges {
            assert_eq!(coedges[c.next.0 as usize].prev,c.id);
            assert_eq!(coedges[c.prev.0 as usize].next,c.id);
            assert_eq!(coedges[c.twin.0 as usize].twin,c.id);
            assert_ne!(coedges[c.twin.0 as usize].reversed,c.reversed);
        }
    }
    #[test]
    fn failed_assembly_never_partially_edits_original_faces() {
        let mut faces=tetrahedron();
        faces[1].boundary.reverse();
        let before=faces.iter().map(|f|(f.id,f.boundary.clone(),f.loops.clone())).collect::<Vec<_>>();
        assert!(from_faces(&mut faces).is_err());
        let after=faces.iter().map(|f|(f.id,f.boundary.clone(),f.loops.clone())).collect::<Vec<_>>();
        assert_eq!(after,before);
    }
}

// Clean analytic B-rep v2. The old incidence types above remain only for the
// transitional G4-G8 monolithic caller while its clients are migrated.
pub mod canonical;
pub use canonical::{BrepDraft, BrepModel, BrepIssue, GeometryCatalog,
    BrepVertex, BrepEdge, BrepCoedge, BrepLoop, BrepFace, BrepShell,
    EdgeKind, FaceSense, ShellSense, AffineParameterMap, LocalBrepIncarnation,
    Curve2Id, Curve3Id, SurfaceId};
