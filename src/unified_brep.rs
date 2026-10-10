//! Shared geometry-backed B-rep graph for planar polyhedra and analytic solids.
//!
//! This is a validated, read-only *canonical interchange view*, not a polygon
//! conversion of curved geometry. The source kernels still own edit operations.
//! Face/edge/coedge reference one shared surface/curve/p-curve catalog.
//! It supports closed circle edges and a same-face pair of seam coedges, unlike
//! the older polyhedral-only validator. Higher radial valence and multi-shell
//! non-manifold models are deliberately not claimed yet.
use crate::{
    BoundingBox, Coedge, CoedgeGeometry, CoedgeId, Curve2, Curve3,
    CylindricalBrep, Edge, EdgeGeometry, EdgeId, FaceGeometry, FaceId,
    GeometryError, GeometryTolerance, Loop, LoopId, MassProperties,
    Mesh, ModelUnits, ParameterRange, Point3, Shell, Solid, Surface3, Vertex,
};
use std::collections::BTreeSet;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn distance(a: Point3, b: Point3) -> f64 {
    let d = a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}
fn parameter(range: ParameterRange, fraction: f64) -> Result<f64, GeometryError> {
    match range {
        ParameterRange::Bounded { start, end } if start.is_finite() && end.is_finite() && end > start => {
            Ok(start + (end - start) * fraction)
        }
        _ => Err(invalid("B-rep topological geometry requires finite bounded parameter intervals")),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepOrigin { Polyhedral, AnalyticCylinder }

/// Generic face: surface identity and trim loops, with no polygonal vertex-ring
/// requirement. Closed circle loops can consist of a single directed coedge.
#[derive(Debug, Clone)]
pub struct BrepFace {
    pub id: FaceId,
    pub loops: Vec<LoopId>,
    pub geometry: FaceGeometry,
}

#[derive(Debug, Clone)]
pub struct BrepGeometry {
    pub curves3: Vec<Curve3>,
    pub curves2: Vec<Curve2>,
    pub surfaces: Vec<Surface3>,
    pub edges: Vec<EdgeGeometry>,
    pub coedges: Vec<CoedgeGeometry>,
}

/// The same topology + geometry graph is built from either source. Curved
/// surfaces and circle edges remain analytic; the display mesh is not stored
/// or used to define B-rep topology.
#[derive(Debug, Clone)]
pub struct BrepModel {
    pub origin: BrepOrigin,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub coedges: Vec<Coedge>,
    pub loops: Vec<Loop>,
    pub faces: Vec<BrepFace>,
    pub shells: Vec<Shell>,
    pub geometry: BrepGeometry,
    pub mass: MassProperties,
    pub bbox: BoundingBox,
    pub tolerance: GeometryTolerance,
    pub units: ModelUnits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepSummary {
    pub vertices: usize,
    pub edges: usize,
    pub coedges: usize,
    pub loops: usize,
    pub faces: usize,
    pub shells: usize,
    pub plane_faces: usize,
    pub cylindrical_faces: usize,
    pub spherical_faces: usize,
    pub straight_edges: usize,
    pub circular_edges: usize,
    pub closed_edges: usize,
    pub seam_edges: usize,
    pub euler_characteristic: isize,
    pub genus: usize,
}

/// Extensible adapter for additional analytic or faceted B-rep producers.
/// It preserves the independent mathematical geometry in the shared model.
pub trait BrepSource {
    fn brep_model(&self) -> Result<BrepModel, GeometryError>;
    fn display_mesh(&self, cylinder_segments: usize) -> Result<Mesh, GeometryError>;
}
impl BrepSource for Solid {
    fn brep_model(&self) -> Result<BrepModel, GeometryError> {
        BrepModel::from_polyhedral(self, ModelUnits::default())
    }
    fn display_mesh(&self, _cylinder_segments: usize) -> Result<Mesh, GeometryError> {
        self.validate()?;
        Ok(self.mesh.clone())
    }
}
impl BrepSource for CylindricalBrep {
    fn brep_model(&self) -> Result<BrepModel, GeometryError> {
        BrepModel::from_cylinder(self)
    }
    fn display_mesh(&self, cylinder_segments: usize) -> Result<Mesh, GeometryError> {
        self.tessellate(cylinder_segments)
    }
}

/// Explicitly separates the origin-specific editable B-rep models from their
/// shared query and traversal interface. This enum is not the final mutable
/// heterogeneous kernel; edits still belong to the source representation.
#[derive(Debug, Clone)]
pub enum BrepBody {
    Polyhedral(Solid),
    Cylinder(CylindricalBrep),
}
impl From<Solid> for BrepBody {
    fn from(body: Solid) -> Self { Self::Polyhedral(body) }
}
impl From<CylindricalBrep> for BrepBody {
    fn from(body: CylindricalBrep) -> Self { Self::Cylinder(body) }
}
impl BrepBody {
    pub fn validate(&self) -> Result<(), GeometryError> {
        match self {
            Self::Polyhedral(body) => body.validate(),
            Self::Cylinder(body) => body.validate(),
        }
    }
    pub fn shared(&self) -> Result<BrepModel, GeometryError> {
        self.source().brep_model()
    }
    pub fn tessellate(&self, cylinder_segments: usize) -> Result<Mesh, GeometryError> {
        self.source().display_mesh(cylinder_segments)
    }
    pub fn source(&self) -> &dyn BrepSource {
        match self {
            Self::Polyhedral(body) => body,
            Self::Cylinder(body) => body,
        }
    }
}

impl BrepModel {
    pub fn from_polyhedral(body: &Solid, units: ModelUnits) -> Result<Self, GeometryError> {
        body.validate()?;
        let store = body.geometry_store(units)?;
        let edges = body.edges.iter().map(|e| {
            store.edges.get(&e.id).copied()
                .ok_or_else(|| invalid(format!("missing analytic binding for edge {:?}", e.id)))
        }).collect::<Result<Vec<_>,_>>()?;
        let coedges = body.coedges.iter().map(|c| {
            store.coedges.get(&c.id).copied()
                .ok_or_else(|| invalid(format!("missing UV binding for coedge {:?}", c.id)))
        }).collect::<Result<Vec<_>,_>>()?;
        let faces = body.faces.iter().map(|f| {
            let geometry = store.faces.get(&f.id).cloned()
                .ok_or_else(|| invalid(format!("missing surface binding for face {:?}", f.id)))?;
            Ok(BrepFace { id: f.id, loops: f.loops.clone(), geometry })
        }).collect::<Result<Vec<_>,GeometryError>>()?;
        let result = Self {
            origin: BrepOrigin::Polyhedral,
            vertices: body.vertices.clone(), edges: body.edges.clone(),
            coedges: body.coedges.clone(), loops: body.loops.clone(),
            faces, shells: body.shells.clone(),
            geometry: BrepGeometry {
                curves3: store.curves3, curves2: store.curves2,
                surfaces: store.surfaces, edges, coedges,
            },
            mass: body.mass, bbox: body.bbox, tolerance: body.tolerance, units,
        };
        result.validate()?;
        Ok(result)
    }

    pub fn from_cylinder(body: &CylindricalBrep) -> Result<Self, GeometryError> {
        body.validate()?;
        let faces = body.faces.iter().map(|f| BrepFace {
            id: f.id, loops: f.loops.clone(), geometry: f.geometry.clone(),
        }).collect();
        let result = Self {
            origin: BrepOrigin::AnalyticCylinder,
            vertices: body.vertices.clone(), edges: body.edges.clone(),
            coedges: body.coedges.clone(), loops: body.loops.clone(),
            faces, shells: body.shells.clone(),
            geometry: BrepGeometry {
                curves3: body.curves3.clone(), curves2: body.curves2.clone(),
                surfaces: body.surfaces.clone(), edges: body.edge_geometry.clone(),
                coedges: body.coedge_geometry.clone(),
            },
            mass: body.mass, bbox: body.bbox, tolerance: body.tolerance, units: body.units,
        };
        result.validate()?;
        Ok(result)
    }

    pub fn face_surface(&self, id: FaceId) -> Option<Surface3> {
        let face = self.faces.get(id.0 as usize).filter(|f| f.id == id)?;
        self.geometry.surfaces.get(face.geometry.surface.0 as usize).copied()
    }
    pub fn edge_curve(&self, id: EdgeId) -> Option<Curve3> {
        let edge = self.edges.get(id.0 as usize).filter(|e| e.id == id)?;
        let geom = self.geometry.edges.get(edge.id.0 as usize)?;
        self.geometry.curves3.get(geom.curve.0 as usize).copied()
    }
    pub fn coedge_pcurve(&self, id: CoedgeId) -> Option<Curve2> {
        let coedge = self.coedges.get(id.0 as usize).filter(|c| c.id == id)?;
        let geom = self.geometry.coedges.get(coedge.id.0 as usize)?;
        self.geometry.curves2.get(geom.pcurve.0 as usize).copied()
    }
    pub fn next_in_loop(&self, id: CoedgeId) -> Option<CoedgeId> {
        self.coedges.get(id.0 as usize).filter(|c| c.id == id).map(|c| c.next)
    }
    /// Two-use manifold radial traversal. A general radial cycle is reserved
    /// for a future storage model with arbitrary coedges per edge.
    pub fn next_of_edge(&self, id: CoedgeId) -> Option<CoedgeId> {
        self.coedges.get(id.0 as usize).filter(|c| c.id == id).map(|c| c.twin)
    }
    pub fn coedges_of_edge(&self, id: EdgeId) -> Option<[CoedgeId; 2]> {
        self.edges.get(id.0 as usize).filter(|e| e.id == id).map(|e| e.coedges)
    }
    pub fn loops_of_face(&self, id: FaceId) -> Option<&[LoopId]> {
        self.faces.get(id.0 as usize).filter(|f| f.id == id).map(|f| f.loops.as_slice())
    }
    pub fn euler_characteristic(&self) -> isize {
        self.vertices.len() as isize - self.edges.len() as isize + self.faces.len() as isize
    }
    pub fn summary(&self) -> BrepSummary {
        let mut out = BrepSummary {
            vertices: self.vertices.len(), edges: self.edges.len(),
            coedges: self.coedges.len(), loops: self.loops.len(),
            faces: self.faces.len(), shells: self.shells.len(),
            plane_faces: 0, cylindrical_faces: 0, spherical_faces: 0,
            straight_edges: 0, circular_edges: 0, closed_edges: 0,
            seam_edges: 0, euler_characteristic: self.euler_characteristic(), genus: 0,
        };
        for face in &self.faces {
            match self.face_surface(face.id) {
                Some(Surface3::Plane(_)) => out.plane_faces += 1,
                Some(Surface3::Cylinder(_)) => out.cylindrical_faces += 1,
                Some(Surface3::Sphere(_)) => out.spherical_faces += 1,
                None => {},
            }
        }
        for edge in &self.edges {
            if edge.start == edge.end { out.closed_edges += 1; }
            match self.edge_curve(edge.id) {
                Some(Curve3::Line(_)) => out.straight_edges += 1,
                Some(Curve3::Circle(_)) => out.circular_edges += 1,
                None => {},
            }
            if (edge.coedges[0].0 as usize) < self.coedges.len()
                && (edge.coedges[1].0 as usize) < self.coedges.len()
                && self.coedges[edge.coedges[0].0 as usize].face == self.coedges[edge.coedges[1].0 as usize].face {
                out.seam_edges += 1;
            }
        }
        let deficit = 2 - out.euler_characteristic;
        if deficit >= 0 && deficit % 2 == 0 { out.genus = (deficit / 2) as usize; }
        out
    }

    /// One shared validator for both input B-rep representations. In addition
    /// to topology graph connectivity, verifies 3D edge and UV-trim carriers
    /// at interior points. This checks sample consistency, not an exact proof
    /// for arbitrary curves, surface self-intersections or Boolean validity.
    pub fn validate(&self) -> Result<(), GeometryError> {
        self.tolerance.validate()?;
        if self.vertices.is_empty() || self.edges.is_empty() || self.faces.is_empty() || self.shells.len() != 1 {
            return Err(invalid("expected one nonempty closed B-rep shell"));
        }
        if self.geometry.edges.len() != self.edges.len()
            || self.geometry.coedges.len() != self.coedges.len()
            || self.faces.iter().any(|f| f.geometry.surface.0 as usize >= self.geometry.surfaces.len()) {
            return Err(invalid("incomplete geometric carrier references"));
        }
        for (i,v) in self.vertices.iter().enumerate() {
            if v.id.0 as usize != i || !v.position.is_finite() {
                return Err(invalid("invalid B-rep vertex identity or coordinates"));
            }
        }
        for (i,face) in self.faces.iter().enumerate() {
            if face.id.0 as usize != i || face.loops.is_empty()
                || face.geometry.outer_loop != face.loops[0] {
                return Err(invalid("invalid B-rep face identity or outer-loop owner"));
            }
        }
        let shell = &self.shells[0];
        if shell.id.0 != 0 || !shell.closed || shell.faces.len() != self.faces.len() {
            return Err(invalid("expected one closed shell owning every face"));
        }
        let mut face_owners = vec![0usize; self.faces.len()];
        for id in &shell.faces {
            let Some(value) = face_owners.get_mut(id.0 as usize) else {
                return Err(invalid("shell contains invalid face"));
            };
            *value += 1;
        }
        if face_owners.iter().any(|&n| n != 1) {
            return Err(invalid("shell face ownership is duplicated or incomplete"));
        }
        let mut loop_owners = vec![0usize; self.loops.len()];
        for face in &self.faces {
            for &id in &face.loops {
                let Some(lp) = self.loops.get(id.0 as usize) else {
                    return Err(invalid("face loop out of bounds"));
                };
                if lp.face != face.id { return Err(invalid("face loop owner mismatch")); }
                loop_owners[id.0 as usize] += 1;
            }
        }
        if loop_owners.iter().any(|&n| n != 1) { return Err(invalid("loop has zero or repeated owners")); }
        let mut counts = vec![0usize; self.coedges.len()];
        let mut face_adjacency = vec![BTreeSet::<FaceId>::new(); self.faces.len()];
        for (i,lp) in self.loops.iter().enumerate() {
            if lp.id.0 as usize != i || lp.first_coedge.0 as usize >= self.coedges.len() {
                return Err(invalid("invalid loop identity or starting coedge"));
            }
            let first = lp.first_coedge;
            let mut current = first;
            let mut seen = BTreeSet::new();
            loop {
                let index = current.0 as usize;
                let Some(coedge) = self.coedges.get(index) else { return Err(invalid("loop references invalid coedge")); };
                if !seen.insert(current) { return Err(invalid("coedge loop repeats before closure")); }
                if coedge.id != current || coedge.loop_id != lp.id || coedge.face != lp.face {
                    return Err(invalid("coedge belongs to wrong face or loop"));
                }
                if coedge.edge.0 as usize >= self.edges.len() { return Err(invalid("coedge edge out of bounds")); }
                if coedge.next.0 as usize >= self.coedges.len() || coedge.prev.0 as usize >= self.coedges.len() {
                    return Err(invalid("coedge next or prev out of bounds"));
                }
                if self.coedges[coedge.next.0 as usize].prev != current
                    || self.coedges[coedge.prev.0 as usize].next != current {
                    return Err(invalid("loop next/prev reciprocity violated"));
                }
                let edge = &self.edges[coedge.edge.0 as usize];
                let next_ref = &self.coedges[coedge.next.0 as usize];
                let next_edge = self.edges.get(next_ref.edge.0 as usize)
                    .ok_or_else(|| invalid("next coedge references a missing edge"))?;
                if coedge.end_vertex(edge) != self.coedges[coedge.next.0 as usize].start_vertex(next_edge) {
                    return Err(invalid("consecutive coedges do not share a vertex"));
                }
                counts[index] += 1;
                current = coedge.next;
                if current == first { break; }
                if seen.len() > self.coedges.len() { return Err(invalid("unbounded coedge loop")); }
            }
        }
        if counts.iter().any(|&n| n != 1) { return Err(invalid("orphaned or reused coedge")); }
        for (i,edge) in self.edges.iter().enumerate() {
            if edge.id.0 as usize != i || edge.start.0 as usize >= self.vertices.len()
                || edge.end.0 as usize >= self.vertices.len() {
                return Err(invalid("edge identity or endpoints are invalid"));
            }
            let [a,b] = edge.coedges;
            if a == b || a.0 as usize >= self.coedges.len() || b.0 as usize >= self.coedges.len() {
                return Err(invalid("edge must have two distinct directed uses"));
            }
            let a = &self.coedges[a.0 as usize];
            let b = &self.coedges[b.0 as usize];
            if a.edge != edge.id || b.edge != edge.id || a.twin != b.id
                || b.twin != a.id || a.reversed == b.reversed {
                return Err(invalid("invalid opposite coedge incidence"));
            }
            if edge.start == edge.end && !matches!(self.edge_curve(edge.id), Some(Curve3::Circle(_))) {
                return Err(invalid("closed edge needs a closed circular geometric carrier"));
            }
            if a.face == b.face {
                let face = self.faces.get(a.face.0 as usize).ok_or_else(|| invalid("seam face out of bounds"))?;
                let bounds = face.geometry.bounds;
                if !bounds.u_periodic && !bounds.v_periodic {
                    return Err(invalid("same-face edge uses require periodic surface domain"));
                }
                let trim_a = self.geometry.coedges[a.id.0 as usize];
                let trim_b = self.geometry.coedges[b.id.0 as usize];
                let ua = self.geometry.curves2.get(trim_a.pcurve.0 as usize)
                    .ok_or_else(|| invalid("seam pcurve missing"))?.evaluate(parameter(trim_a.domain,0.5)?)?;
                let ub = self.geometry.curves2.get(trim_b.pcurve.0 as usize)
                    .ok_or_else(|| invalid("seam second pcurve missing"))?.evaluate(parameter(trim_b.domain,0.5)?)?;
                let d = if bounds.u_periodic {(ua.u-ub.u).abs()} else {(ua.v-ub.v).abs()};
                if !d.is_finite() || d <= self.tolerance.angular {
                    return Err(invalid("periodic seam must have distinct UV representations"));
                }
            } else {
                face_adjacency[a.face.0 as usize].insert(b.face);
                face_adjacency[b.face.0 as usize].insert(a.face);
            }
        }
        let mut connected = BTreeSet::new();
        let mut pending = vec![FaceId(0)];
        while let Some(face) = pending.pop() {
            if !connected.insert(face) { continue; }
            for &other in &face_adjacency[face.0 as usize] { pending.push(other); }
        }
        if connected.len() != self.faces.len() { return Err(invalid("closed shell is not face-connected")); }
        let chi = self.euler_characteristic();
        if chi > 2 || (2-chi) % 2 != 0 {
            return Err(invalid("orientable connected shell has invalid Euler characteristic"));
        }
        let mut edge_usage_counts = vec![0usize; self.edges.len()];
        for (i,coedge) in self.coedges.iter().enumerate() {
            if let Some(usage) = edge_usage_counts.get_mut(coedge.edge.0 as usize) {
                *usage += 1;
            } else {
                return Err(invalid("coedge references missing parent edge"));
            }
            if coedge.id.0 as usize != i { return Err(invalid("misindexed coedge")); }
            let carrier = self.geometry.edges.get(coedge.edge.0 as usize)
                .ok_or_else(|| invalid("coedge missing edge carrier"))?;
            let trim_geom = self.geometry.coedges.get(i)
                .ok_or_else(|| invalid("coedge missing UV carrier"))?;
            if trim_geom.face != coedge.face || trim_geom.loop_id != coedge.loop_id {
                return Err(invalid("UV trim owner mismatch"));
            }
            let face = self.faces.get(coedge.face.0 as usize)
                .ok_or_else(|| invalid("UV trim face out of range"))?;
            let surface = self.geometry.surfaces.get(face.geometry.surface.0 as usize)
                .ok_or_else(|| invalid("face missing geometric surface"))?;
            let curve = self.geometry.curves3.get(carrier.curve.0 as usize)
                .ok_or_else(|| invalid("edge missing 3D curve"))?;
            let trim = self.geometry.curves2.get(trim_geom.pcurve.0 as usize)
                .ok_or_else(|| invalid("coedge missing 2D curve"))?;
            let edge = &self.edges[coedge.edge.0 as usize];
            // Permit inevitable floating-point reconstruction noise at very
            // large world coordinates, without using that slack to validate
            // topology or permit sub-tolerance physical features.
            let coordinate_scale = [self.vertices[edge.start.0 as usize].position,
                self.vertices[edge.end.0 as usize].position].into_iter()
                .flat_map(|p| [p.x.abs(),p.y.abs(),p.z.abs()])
                .fold(1.0_f64,f64::max);
            let eps = carrier.tolerance.length_at(self.bbox_extent())
                .max(16.0*f64::EPSILON*coordinate_scale);
            for fraction in [0.0,0.25,0.5,0.75,1.0] {
                let uv = trim.evaluate(parameter(trim_geom.domain,fraction)?)?;
                let on_surface = surface.evaluate(uv.u,uv.v)?.point;
                let oriented = if coedge.reversed {1.0-fraction} else {fraction};
                let on_edge = curve.evaluate(parameter(carrier.domain,oriented)?)?;
                if distance(on_surface,on_edge) > eps {
                    return Err(invalid(format!("coedge {:?} UV/3D representations disagree",coedge.id)));
                }
                if fraction == 0.0 || fraction == 1.0 {
                    let vertex = if fraction == 0.0 { coedge.start_vertex(edge) }
                        else { coedge.end_vertex(edge) };
                    let at_vertex = self.vertices[vertex.0 as usize].position;
                    if distance(on_surface,at_vertex) > eps {
                        return Err(invalid("coedge boundary misses topological vertex"));
                    }
                }
            }
        }
        if edge_usage_counts.iter().any(|&count| count != 2) {
            return Err(invalid("each manifold edge must have exactly two coedges"));
        }
        if !self.mass.volume.is_finite() || self.mass.volume <= 0.0
            || !self.mass.surface_area.is_finite() || self.mass.surface_area <= 0.0
            || !self.mass.centroid.is_finite() || !self.bbox.min.is_finite()
            || !self.bbox.max.is_finite() {
            return Err(invalid("invalid B-rep bounds or mass properties"));
        }
        Ok(())
    }
    fn bbox_extent(&self) -> f64 {
        let v = self.bbox.max.sub(self.bbox.min);
        v.x.hypot(v.y).hypot(v.z)
    }
}

