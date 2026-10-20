//! Oriented B-rep incidence data, stable within a single model revision.
//! Pure topology data cannot depend on Booleans, meshing, or document storage.
use rustsolid_foundation::Point3;
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
