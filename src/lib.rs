//! RustSolid: an independent geometric modeling kernel for Gefest CAD.
//!
//! This module implements closed *planar polyhedral* boundary representations.
//! It is **not** a full analytic/NURBS B-rep or a Parasolid-compatible kernel.
//! The source is original and does not incorporate decompiled proprietary code.

use std::collections::BTreeMap;
use std::fmt;

mod brep;
mod tolerance;
mod edit;
mod identity;
mod journal;
mod history;
pub use brep::{Coedge, CoedgeId, Loop, LoopId, LoopRole, Shell, ShellId, TopologyEntity, TopologyIssue};
pub use tolerance::GeometryTolerance;
pub use edit::{EdgeSplit, FaceSplit, EditDelta, EditReport, EditTransaction};
pub use identity::{TopologyName, TopologyHandle, HandleError, NameChange, TopologyIdentity};
pub use journal::JournalStats;
pub use history::{EulerCommand, CommandBatch, CommandHistory, HistoryDirection, HistoryTransition};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeometryError {
    InvalidProfile(&'static str),
    InvalidDimension(&'static str),
    InvalidTopology(String),
    InvalidTolerance(&'static str),
    InvalidEdit(String),
    TriangulationFailed,
}

impl fmt::Display for GeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidProfile(reason) => write!(f, "invalid extrusion profile: {reason}"),
            Self::InvalidDimension(reason) => write!(f, "invalid dimension: {reason}"),
            Self::InvalidTopology(reason) => write!(f, "invalid topology: {reason}"),
            Self::InvalidTolerance(reason) => write!(f, "invalid tolerance: {reason}"),
            Self::InvalidEdit(reason) => write!(f, "invalid edit: {reason}"),
            Self::TriangulationFailed => write!(f, "polygon triangulation failed"),
        }
    }
}

impl std::error::Error for GeometryError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    pub x: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub fn array(self) -> [f64; 3] { [self.x, self.y, self.z] }
    pub fn is_finite(self) -> bool { self.x.is_finite() && self.y.is_finite() && self.z.is_finite() }
    pub(crate) fn add(self, v: Point3) -> Self {
        Self { x: self.x + v.x, y: self.y + v.y, z: self.z + v.z }
    }
    pub(crate) fn sub(self, v: Point3) -> Self {
        Self { x: self.x - v.x, y: self.y - v.y, z: self.z - v.z }
    }
    pub(crate) fn scale(self, s: f64) -> Self {
        Self { x: self.x * s, y: self.y * s, z: self.z * s }
    }
    pub(crate) fn dot(self, v: Point3) -> f64 { self.x * v.x + self.y * v.y + self.z * v.z }
    pub(crate) fn cross(self, v: Point3) -> Self {
        Self {
            x: self.y * v.z - self.z * v.y,
            y: self.z * v.x - self.x * v.z,
            z: self.x * v.y - self.y * v.x,
        }
    }
    pub(crate) fn length_squared(self) -> f64 { self.dot(self) }
}

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

#[derive(Debug, Clone)]
pub struct Mesh {
    pub vertices: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
    /// A triangle-to-B-rep-face ownership map for face selection.
    pub triangle_faces: Vec<FaceId>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min: Point3,
    pub max: Point3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    pub volume: f64,
    pub surface_area: f64,
    pub centroid: Point3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayHit {
    pub face_id: FaceId,
    pub distance: f64,
    pub point: Point3,
    pub normal: Point3,
}

#[derive(Debug, Clone)]
pub struct Solid {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub coedges: Vec<Coedge>,
    pub loops: Vec<Loop>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub(crate) identities: TopologyIdentity,
    pub tolerance: GeometryTolerance,
    /// Local edit revision (not a persistent feature identifier).
    pub revision: u64,
    pub mesh: Mesh,
    pub bbox: BoundingBox,
    pub mass: MassProperties,
}

fn cross2(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.z - a.z) - (b.z - a.z) * (c.x - a.x)
}
fn dist2(a: Point2, b: Point2) -> f64 {
    (a.x - b.x).powi(2) + (a.z - b.z).powi(2)
}
fn within(a: f64, p: f64, b: f64, eps: f64) -> bool {
    p >= a.min(b) - eps && p <= a.max(b) + eps
}
fn segments_intersect(a: Point2, b: Point2, c: Point2, d: Point2, eps: f64, area_eps: f64) -> bool {
    let ab_c = cross2(a, b, c);
    let ab_d = cross2(a, b, d);
    let cd_a = cross2(c, d, a);
    let cd_b = cross2(c, d, b);
    if ((ab_c > area_eps && ab_d < -area_eps) || (ab_c < -area_eps && ab_d > area_eps))
        && ((cd_a > area_eps && cd_b < -area_eps) || (cd_a < -area_eps && cd_b > area_eps)) {
        return true;
    }
    (ab_c.abs() <= area_eps && within(a.x, c.x, b.x, eps) && within(a.z, c.z, b.z, eps))
        || (ab_d.abs() <= area_eps && within(a.x, d.x, b.x, eps) && within(a.z, d.z, b.z, eps))
        || (cd_a.abs() <= area_eps && within(c.x, a.x, d.x, eps) && within(c.z, a.z, d.z, eps))
        || (cd_b.abs() <= area_eps && within(c.x, b.x, d.x, eps) && within(c.z, b.z, d.z, eps))
}

fn signed_area(points: &[Point2]) -> f64 {
    let first = points[0];
    let mut twice = 0.0;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        twice += (a.x - first.x) * (b.z - first.z) - (b.x - first.x) * (a.z - first.z);
    }
    0.5 * twice
}

fn normalized_profile(profile: &[Point2], tolerance: GeometryTolerance) -> Result<(Vec<Point2>, f64), GeometryError> {
    if profile.len() < 3 { return Err(GeometryError::InvalidProfile("need at least three vertices")); }
    if profile.len() > (u32::MAX as usize) / 6 { return Err(GeometryError::InvalidProfile("too many vertices")); }
    if !profile.iter().all(|p| p.x.is_finite() && p.z.is_finite()) {
        return Err(GeometryError::InvalidProfile("coordinates must be finite"));
    }
    let minx = profile.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let maxx = profile.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let minz = profile.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
    let maxz = profile.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);
    let extent = (maxx - minx).max(maxz - minz);
    if !extent.is_finite() || extent <= 0.0 {
        return Err(GeometryError::InvalidProfile("profile has zero extent"));
    }
    let eps = tolerance.length_at(extent);
    let area_eps = tolerance.area_at(extent);
    if !eps.is_finite() || !area_eps.is_finite() || extent <= eps {
        return Err(GeometryError::InvalidProfile("profile extent is below tolerance"));
    }
    let mut points = profile.to_vec();
    // Explicit closing vertex is accepted, but zero-length interior edges are rejected.
    if dist2(points[0], *points.last().unwrap()) <= eps * eps { points.pop(); }
    if points.len() < 3 { return Err(GeometryError::InvalidProfile("need three unique vertices")); }
    for i in 0..points.len() {
        if dist2(points[i], points[(i + 1) % points.len()]) <= eps * eps {
            return Err(GeometryError::InvalidProfile("duplicate consecutive vertices"));
        }
    }
    // Reject crossings (and touching non-neighbouring edges) before triangulation.
    let n = points.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if j == i + 1 || (i == 0 && j == n - 1) { continue; }
            if segments_intersect(points[i], points[(i + 1) % n],
                                  points[j], points[(j + 1) % n], eps, area_eps) {
                return Err(GeometryError::InvalidProfile("self-intersecting or self-touching boundary"));
            }
        }
    }
    // A straight-through vertex adds no geometric information and can create zero-area cap triangles.
    loop {
        if points.len() <= 3 { break; }
        let mut removed = false;
        for i in 0..points.len() {
            let prev = points[(i + points.len() - 1) % points.len()];
            let next = points[(i + 1) % points.len()];
            if cross2(prev, points[i], next).abs() <= area_eps
                && within(prev.x, points[i].x, next.x, eps)
                && within(prev.z, points[i].z, next.z, eps) {
                points.remove(i);
                removed = true;
                break;
            }
        }
        if !removed { break; }
    }
    let mut area = signed_area(&points);
    if !area.is_finite() || area.abs() <= area_eps {
        return Err(GeometryError::InvalidProfile("profile area is too small"));
    }
    if area < 0.0 { points.reverse(); area = -area; }
    Ok((points, area))
}

fn point_in_triangle(p: Point2, a: Point2, b: Point2, c: Point2, eps: f64) -> bool {
    cross2(a, b, p) >= -eps && cross2(b, c, p) >= -eps && cross2(c, a, p) >= -eps
}

fn triangulate(points: &[Point2], tolerance: GeometryTolerance) -> Result<Vec<[u32; 3]>, GeometryError> {
    let minx = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let maxx = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let minz = points.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
    let maxz = points.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);
    let area_eps = tolerance.area_at((maxx - minx).max(maxz - minz));
    let mut ring: Vec<u32> = (0..points.len() as u32).collect();
    let mut out = Vec::with_capacity(points.len() - 2);
    while ring.len() > 3 {
        let mut found = false;
        for i in 0..ring.len() {
            let ai = ring[(i + ring.len() - 1) % ring.len()];
            let bi = ring[i];
            let ci = ring[(i + 1) % ring.len()];
            let (a, b, c) = (points[ai as usize], points[bi as usize], points[ci as usize]);
            if cross2(a, b, c) <= area_eps { continue; }
            if ring.iter().any(|&p| p != ai && p != bi && p != ci
                && point_in_triangle(points[p as usize], a, b, c, area_eps)) { continue; }
            out.push([ai, bi, ci]);
            ring.remove(i);
            found = true;
            break;
        }
        if !found { return Err(GeometryError::TriangulationFailed); }
    }
    let (a, b, c) = (ring[0], ring[1], ring[2]);
    if cross2(points[a as usize], points[b as usize], points[c as usize]) <= area_eps {
        return Err(GeometryError::TriangulationFailed);
    }
    out.push([a, b, c]);
    Ok(out)
}

impl Solid {
    /// Extrudes a simple XZ-plane polygon toward +Y, matching Gefest's viewport convention.
    /// The input must not contain holes; the output is a closed, oriented planar B-rep.
    pub fn extrude_xz(profile: &[Point2], height: f64) -> Result<Self, GeometryError> {
        Self::extrude_xz_with_tolerance(profile, height, GeometryTolerance::default())
    }

    /// Build a prism with explicit, validated model-space tolerance semantics.
    pub fn extrude_xz_with_tolerance(profile: &[Point2], height: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !height.is_finite() || height <= 0.0 {
            return Err(GeometryError::InvalidDimension("extrusion height must be finite and positive"));
        }
        let (points, area) = normalized_profile(profile, tolerance)?;
        let cap_tris = triangulate(&points, tolerance)?;
        let extent = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
            - points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let extent_z = points.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max)
            - points.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
        if height <= tolerance.length_at(extent.max(extent_z).max(height)) {
            return Err(GeometryError::InvalidDimension("extrusion height is below tolerance"));
        }
        let n = points.len();
        let mut vertices = Vec::with_capacity(2 * n);
        for (i, point) in points.iter().enumerate() {
            vertices.push(Vertex { id: VertexId(i as u32), position: Point3 {x: point.x, y: 0.0, z: point.z} });
        }
        for (i, point) in points.iter().enumerate() {
            vertices.push(Vertex { id: VertexId((n + i) as u32), position: Point3 {x: point.x, y: height, z: point.z} });
        }
        let bottom = FaceId(0);
        let top = FaceId(1);
        let mut faces = vec![
            Face { id: bottom, role: FaceRole::BottomCap,
                boundary: (0..n as u32).map(VertexId).collect(), loops: Vec::new() },
            Face { id: top, role: FaceRole::TopCap,
                boundary: (0..n as u32).rev().map(|i| VertexId(n as u32 + i)).collect(), loops: Vec::new() },
        ];
        for i in 0..n {
            let j = (i + 1) % n;
            faces.push(Face {
                id: FaceId((i + 2) as u32), role: FaceRole::Side(i as u32),
                boundary: vec![VertexId(i as u32), VertexId((n + i) as u32),
                               VertexId((n + j) as u32), VertexId(j as u32)],
                loops: Vec::new(),
            });
        }
        let (edges, coedges, loops, shells) = brep::from_faces(&mut faces)?;
        let mut triangles = Vec::with_capacity(4 * n - 4);
        let mut triangle_faces = Vec::with_capacity(4 * n - 4);
        for [a, b, c] in cap_tris {
            triangles.push([a, b, c]); triangle_faces.push(bottom);
            triangles.push([n as u32 + c, n as u32 + b, n as u32 + a]);
            triangle_faces.push(top);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            let a = i as u32; let b = (n + i) as u32;
            let c = (n + j) as u32; let d = j as u32;
            triangles.push([a, b, c]); triangle_faces.push(FaceId((i + 2) as u32));
            triangles.push([a, c, d]); triangle_faces.push(FaceId((i + 2) as u32));
        }
        let perimeter: f64 = (0..n).map(|i| dist2(points[i], points[(i + 1) % n]).sqrt()).sum();
        // Compute moments relative to the first vertex. Absolute-coordinate shoelace
        // moments suffer catastrophic cancellation for small parts far from origin.
        let anchor = points[0];
        let (mut moment_x, mut moment_z) = (0.0, 0.0);
        for i in 0..n {
            let p = Point2 { x: points[i].x - anchor.x, z: points[i].z - anchor.z };
            let q = Point2 { x: points[(i + 1) % n].x - anchor.x,
                             z: points[(i + 1) % n].z - anchor.z };
            let cross = p.x * q.z - q.x * p.z;
            moment_x += (p.x + q.x) * cross;
            moment_z += (p.z + q.z) * cross;
        }
        let centroid = Point3 {
            x: anchor.x + moment_x / (6.0 * area),
            y: height / 2.0,
            z: anchor.z + moment_z / (6.0 * area),
        };
        let bbox = BoundingBox {
            min: Point3 { x: points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
                          y: 0.0, z: points.iter().map(|p| p.z).fold(f64::INFINITY, f64::min) },
            max: Point3 { x: points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
                          y: height, z: points.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max) },
        };
        let mesh = Mesh { vertices: vertices.iter().map(|v| v.position).collect(), triangles, triangle_faces };
        let mut solid = Solid { vertices, edges, coedges, loops, faces, shells, identities: TopologyIdentity::empty(), tolerance, revision: 0, mesh, bbox,
                            mass: MassProperties { volume: area * height,
                                surface_area: 2.0 * area + perimeter * height, centroid } };
        let mut identity = solid.identities.clone();
        identity.initialize(&solid)?;
        solid.identities = identity;
        solid.validate()?;
        Ok(solid)
    }

    pub fn block(origin: Point3, width: f64, height: f64, depth: f64) -> Result<Self, GeometryError> {
        Self::block_with_tolerance(origin, width, height, depth, GeometryTolerance::default())
    }

    pub fn block_with_tolerance(origin: Point3, width: f64, height: f64, depth: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !origin.is_finite() { return Err(GeometryError::InvalidDimension("origin must be finite")); }
        if !width.is_finite() || width <= 0.0 || !depth.is_finite() || depth <= 0.0 {
            return Err(GeometryError::InvalidDimension("block width and depth must be finite and positive"));
        }
        let profile = [Point2{x: 0.0,z: 0.0}, Point2{x: width,z: 0.0},
                       Point2{x: width,z: depth}, Point2{x: 0.0,z: depth}];
        Self::extrude_xz_with_tolerance(&profile, height, tolerance)?.translated(origin)
    }

    pub fn translated(mut self, delta: Point3) -> Result<Self, GeometryError> {
        if !delta.is_finite() { return Err(GeometryError::InvalidDimension("translation must be finite")); }
        for v in &mut self.vertices { v.position = v.position.add(delta); }
        for v in &mut self.mesh.vertices { *v = v.add(delta); }
        self.mass.centroid = self.mass.centroid.add(delta);
        self.bbox.min = self.bbox.min.add(delta);
        self.bbox.max = self.bbox.max.add(delta);
        if !self.vertices.iter().all(|v| v.position.is_finite()) {
            return Err(GeometryError::InvalidDimension("translation overflow"));
        }
        self.validate()?;
        Ok(self)
    }

    pub fn euler_characteristic(&self) -> isize {
        self.vertices.len() as isize - self.edges.len() as isize + self.faces.len() as isize
    }

    /// Checks B-rep graph invariants, geometric face tolerances, and mesh ownership.
    pub fn validate(&self) -> Result<(), GeometryError> {
        self.tolerance.validate()?;
        self.check_topology().map_err(|issue| GeometryError::InvalidTopology(issue.to_string()))?;
        self.identities.validate(self)?;
        if self.mesh.vertices.len() != self.vertices.len()
            || self.mesh.triangles.len() != self.mesh.triangle_faces.len() {
            return Err(GeometryError::InvalidTopology("mesh indices or ownership invalid".into()));
        }
        let extent = self.bbox.max.sub(self.bbox.min).length_squared().sqrt();
        let linear = self.tolerance.length_at(extent);
        if !extent.is_finite() || extent <= linear {
            return Err(GeometryError::InvalidTopology("body extent is below tolerance".into()));
        }
        for (vertex, source) in self.mesh.vertices.iter().zip(&self.vertices) {
            if !vertex.is_finite() || vertex.sub(source.position).length_squared().sqrt() > linear {
                return Err(GeometryError::InvalidTopology("mesh and B-rep vertex positions disagree".into()));
            }
        }
        let mut mesh_edge_uses: BTreeMap<(u32,u32), Vec<(bool,FaceId)>> = BTreeMap::new();
        // The outward signed volume of the triangle boundary must agree with
        // the analytic prism volume. Coordinates are anchored for stability.
        let anchor = self.vertices[0].position;
        let mut six_volume = 0.0;
        let mut computed_surface_area = 0.0;
        for (tri, owner) in self.mesh.triangles.iter().zip(&self.mesh.triangle_faces) {
            if owner.0 as usize >= self.faces.len() || tri.iter().any(|&i| i as usize >= self.vertices.len()) {
                return Err(GeometryError::InvalidTopology("mesh index outside bounds".into()));
            }
            let [a,b,c] = *tri;
            if a==b || b==c || c==a {
                return Err(GeometryError::InvalidTopology("degenerate mesh triangle".into()));
            }
            let v0 = self.mesh.vertices[a as usize];
            let v1 = self.mesh.vertices[b as usize];
            let v2 = self.mesh.vertices[c as usize];
            let cross = v1.sub(v0).cross(v2.sub(v0));
            let area2 = cross.length_squared().sqrt();
            if !area2.is_finite() || area2 <= 0.0 {
                return Err(GeometryError::InvalidTopology("zero-area or non-finite mesh triangle".into()));
            }
            let face = &self.faces[owner.0 as usize];
            // v0.2 tessellation uses only boundary vertices; a later curved-face
            // mesher will need a separate face-domain certification contract.
            if !tri.iter().all(|v| face.boundary.contains(&VertexId(*v))) {
                return Err(GeometryError::InvalidTopology("triangle references vertices outside owning face".into()));
            }
            let face_anchor = self.vertices[face.boundary[0].0 as usize].position;
            let mut face_sum = Point3 { x: 0.0, y: 0.0, z: 0.0 };
            for i in 0..face.boundary.len() {
                let x = self.vertices[face.boundary[i].0 as usize].position.sub(face_anchor);
                let y = self.vertices[face.boundary[(i + 1) % face.boundary.len()].0 as usize].position.sub(face_anchor);
                face_sum = face_sum.add(x.cross(y));
            }
            let face_norm = face_sum.length_squared().sqrt();
            if !face_norm.is_finite() || face_norm <= 0.0 || cross.dot(face_sum) <= 0.0 {
                return Err(GeometryError::InvalidTopology("triangle winding disagrees with owning face".into()));
            }
            // Compare unit directions by sine of the deviation angle.
            let sine = cross.cross(face_sum).length_squared().sqrt() / (area2 * face_norm);
            if !sine.is_finite() || sine > self.tolerance.angular.sin() {
                return Err(GeometryError::InvalidTopology("triangle normal deviates from owning face".into()));
            }
            six_volume += v0.sub(anchor).dot(v1.sub(anchor).cross(v2.sub(anchor)));
            computed_surface_area += 0.5 * area2;
            for (u,v) in [(a,b),(b,c),(c,a)] {
                mesh_edge_uses.entry((u.min(v),u.max(v))).or_default().push((u<v,*owner));
            }
        }
        if mesh_edge_uses.values().any(|uses| uses.len()!=2 || uses[0].0==uses[1].0) {
            return Err(GeometryError::InvalidTopology("mesh not closed or inconsistently wound".into()));
        }
        // A topological edge must be represented by exactly the pair of
        // opposing, correctly owned triangle edges in this vertex-sharing
        // polygonal tessellation. Internal triangulation diagonals are exempt.
        for edge in &self.edges {
            let pair=mesh_edge_uses.get(&(edge.start.0,edge.end.0))
                .ok_or_else(|| GeometryError::InvalidTopology("B-rep edge missing from tessellation".into()))?;
            let mut owners=[self.coedges[edge.coedges[0].0 as usize].face,
                self.coedges[edge.coedges[1].0 as usize].face];
            owners.sort();
            let mut mesh_owners=[pair[0].1,pair[1].1];
            mesh_owners.sort();
            if owners != mesh_owners {
                return Err(GeometryError::InvalidTopology("B-rep and mesh edge face ownership disagree".into()));
            }
        }
        let area_budget = self.mass.surface_area.abs()*1e-9
            + self.tolerance.length_at(extent) * extent * 10.0;
        if !computed_surface_area.is_finite() ||
            (computed_surface_area-self.mass.surface_area).abs()>area_budget {
            return Err(GeometryError::InvalidTopology("tessellation surface area disagrees with analytic solid".into()));
        }
        if !self.mass.volume.is_finite() || self.mass.volume <= 0.0
            || !self.mass.surface_area.is_finite() || self.mass.surface_area <= 0.0
            || !self.mass.centroid.is_finite() || !self.bbox.min.is_finite() || !self.bbox.max.is_finite() {
            return Err(GeometryError::InvalidTopology("invalid mass properties".into()));
        }
        // Keep an explicit numerical budget for accumulated signed tetrahedra.
        let computed_volume = six_volume / 6.0;
        let error_budget = self.mass.volume.abs() * 1.0e-9
            + self.tolerance.length_at(extent) * extent * extent * 10.0;
        if !computed_volume.is_finite() || computed_volume <= 0.0 ||
            (computed_volume - self.mass.volume).abs() > error_budget {
            return Err(GeometryError::InvalidTopology("mesh signed volume disagrees with analytic solid".into()));
        }
        Ok(())
    }

    /// Ray-picks the nearest triangle, returning its B-rep face ID.
    pub fn ray_cast(&self, origin: Point3, direction: Point3) -> Option<RayHit> {
        if !origin.is_finite() || !direction.is_finite() || direction.length_squared() <= 0.0 { return None; }
        // Normalize so RayHit.distance is a physical distance even for non-unit input rays.
        let direction=direction.scale(1.0/direction.length_squared().sqrt());
        let mut best: Option<RayHit> = None;
        for (&tri, &face_id) in self.mesh.triangles.iter().zip(&self.mesh.triangle_faces) {
            let a=self.mesh.vertices[tri[0] as usize];
            let b=self.mesh.vertices[tri[1] as usize];
            let c=self.mesh.vertices[tri[2] as usize];
            let e1=b.sub(a); let e2=c.sub(a);
            let p=direction.cross(e2); let determinant=e1.dot(p);
            if determinant.abs() <= 1e-14 { continue; }
            let inv_det=1.0/determinant;
            let t=origin.sub(a);
            let u=t.dot(p)*inv_det;
            if !(-1e-12..=1.0+1e-12).contains(&u) { continue; }
            let q=t.cross(e1);
            let v=direction.dot(q)*inv_det;
            if v < -1e-12 || u+v > 1.0+1e-12 { continue; }
            let distance=e2.dot(q)*inv_det;
            if distance < 0.0 || best.is_some_and(|hit| distance>=hit.distance) { continue; }
            let normal=e1.cross(e2); let norm=normal.length_squared().sqrt();
            if norm==0.0 { continue; }
            best=Some(RayHit { face_id, distance, point:origin.add(direction.scale(distance)),
                               normal:normal.scale(1.0/norm) });
        }
        best
    }
}
