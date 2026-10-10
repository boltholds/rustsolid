//! Independent analytic geometry and parameter-domain layer for RustSolid.
//!
//! The existing polyhedral B-rep is projected into a revision-bound geometry
//! view. Analytic circles/cylinders/spheres are independently evaluable, but
//! curved faces and p-curves are not yet editable topology in v0.7.
use crate::{CoedgeId, EdgeId, FaceId, GeometryError, GeometryTolerance, LoopId, Point3, Solid, TopologyEntity, VertexId};
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, TAU};

fn invalid(reason: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(reason.into())
}
fn invalid_dimension(reason: &'static str) -> GeometryError {
    GeometryError::InvalidDimension(reason)
}
fn length(p: Point3) -> f64 { p.x.hypot(p.y).hypot(p.z) }
fn unit(v: Point3, epsilon: f64) -> Result<Point3, GeometryError> {
    let magnitude = length(v);
    if !magnitude.is_finite() || magnitude <= epsilon {
        return Err(invalid_dimension("geometry direction below tolerance"));
    }
    Ok(v.scale(1.0 / magnitude))
}
fn distance(a: Point3, b: Point3) -> f64 { length(a.sub(b)) }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit { Millimeter, Meter, Inch }
impl LengthUnit {
    pub fn meters_per_unit(self) -> f64 {
        match self { Self::Millimeter => 0.001, Self::Meter => 1.0, Self::Inch => 0.0254 }
    }
    pub fn convert(self, amount: f64, target: LengthUnit) -> Result<f64, GeometryError> {
        let converted = amount * (self.meters_per_unit() / target.meters_per_unit());
        if !converted.is_finite() { return Err(invalid_dimension("length unit conversion overflow")); }
        Ok(converted)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AngleUnit { Radian }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModelUnits { pub length: LengthUnit, pub angle: AngleUnit }
impl Default for ModelUnits {
    fn default() -> Self { Self { length: LengthUnit::Millimeter, angle: AngleUnit::Radian } }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParameterRange {
    Unbounded,
    Bounded { start: f64, end: f64 },
}
impl ParameterRange {
    pub fn bounded(start: f64, end: f64) -> Result<Self, GeometryError> {
        if !start.is_finite() || !end.is_finite() || start >= end {
            return Err(invalid_dimension("parameter interval must have finite ordered endpoints"));
        }
        Ok(Self::Bounded { start, end })
    }
    pub fn includes(self, parameter: f64, slack: f64) -> bool {
        if !parameter.is_finite() || !slack.is_finite() || slack < 0.0 { return false; }
        match self {
            Self::Unbounded => true,
            Self::Bounded { start, end } => parameter >= start - slack && parameter <= end + slack,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParameterDomain2 {
    pub u: ParameterRange,
    pub v: ParameterRange,
    pub u_periodic: bool,
    pub v_periodic: bool,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2Param { pub u: f64, pub v: f64 }
impl Point2Param {
    pub fn new(u: f64, v: f64) -> Result<Self, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("UV parameters must be finite")); }
        Ok(Self { u, v })
    }
}

/// Right-handed orthonormal frame. The stored X/Y directions are explicit;
/// `normal` is derived, not independently guessed from a surface's mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame3 {
    pub origin: Point3,
    pub x: Point3,
    pub y: Point3,
    pub normal: Point3,
}
impl Frame3 {
    pub fn from_axes(origin: Point3, x_axis: Point3, y_axis: Point3,
                     tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !origin.is_finite() || !x_axis.is_finite() || !y_axis.is_finite() {
            return Err(invalid_dimension("frame coordinates must be finite"));
        }
        let x = unit(x_axis, tolerance.absolute_length)?;
        let n = unit(x.cross(y_axis), tolerance.absolute_length * length(y_axis).max(1.0))?;
        let y = unit(n.cross(x), tolerance.absolute_length)?;
        if !x.is_finite() || !y.is_finite() || !n.is_finite() {
            return Err(invalid_dimension("invalid frame basis"));
        }
        Ok(Self { origin, x, y, normal: n })
    }
    pub fn position(self, u: f64, v: f64) -> Result<Point3, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("surface coordinates must be finite")); }
        let p = self.origin.add(self.x.scale(u)).add(self.y.scale(v));
        if !p.is_finite() { return Err(invalid_dimension("surface evaluation overflow")); }
        Ok(p)
    }
    pub fn project(self, point: Point3) -> Result<Point2Param, GeometryError> {
        if !point.is_finite() { return Err(invalid_dimension("projected point must be finite")); }
        let delta = point.sub(self.origin);
        Point2Param::new(delta.dot(self.x), delta.dot(self.y))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line3 { pub origin: Point3, pub direction: Point3 }
impl Line3 {
    pub fn between(a: Point3, b: Point3, tolerance: GeometryTolerance)
        -> Result<(Self, ParameterRange), GeometryError>
    {
        tolerance.validate()?;
        if !a.is_finite() || !b.is_finite() { return Err(invalid_dimension("line coordinates must be finite")); }
        let delta = b.sub(a);
        let span = length(delta);
        if !span.is_finite() || span <= tolerance.length_at(span) {
            return Err(invalid_dimension("line length below tolerance"));
        }
        Ok((Self { origin: a, direction: delta.scale(1.0 / span) }, ParameterRange::bounded(0.0, span)?))
    }
    pub fn evaluate(self, t: f64) -> Result<Point3, GeometryError> {
        if !t.is_finite() { return Err(invalid_dimension("line parameter must be finite")); }
        let p = self.origin.add(self.direction.scale(t));
        if !p.is_finite() { return Err(invalid_dimension("line evaluation overflow")); }
        Ok(p)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle3 { pub frame: Frame3, pub radius: f64 }
impl Circle3 {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("circle radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
    pub fn evaluate(self, angle: f64) -> Result<Point3, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("circle angle must be finite")); }
        let p = self.frame.origin.add(self.frame.x.scale(self.radius * angle.cos()))
            .add(self.frame.y.scale(self.radius * angle.sin()));
        if !p.is_finite() { return Err(invalid_dimension("circle evaluation overflow")); }
        Ok(p)
    }
    pub fn derivative(self, angle: f64) -> Result<Point3, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("circle angle must be finite")); }
        Ok(self.frame.x.scale(-self.radius * angle.sin())
            .add(self.frame.y.scale(self.radius * angle.cos())))
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve3 { Line(Line3), Circle(Circle3) }
impl Curve3 {
    pub fn evaluate(self, t: f64) -> Result<Point3, GeometryError> {
        match self { Self::Line(line) => line.evaluate(t), Self::Circle(circle) => circle.evaluate(t) }
    }
    pub fn derivative(self, t: f64) -> Result<Point3, GeometryError> {
        match self {
            Self::Line(line) => Ok(line.direction),
            Self::Circle(circle) => circle.derivative(t),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line2 { pub origin: Point2Param, pub delta: Point2Param }
impl Line2 {
    pub fn evaluate(self, t: f64) -> Result<Point2Param, GeometryError> {
        Point2Param::new(self.origin.u + t * self.delta.u, self.origin.v + t * self.delta.v)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2 { pub center: Point2Param, pub radius: f64 }
impl Circle2 {
    pub fn new(center: Point2Param, radius: f64, tol: GeometryTolerance) -> Result<Self, GeometryError> {
        tol.validate()?;
        if !radius.is_finite() || radius <= tol.absolute_length {
            return Err(invalid_dimension("parametric circle radius below tolerance"));
        }
        Ok(Self {center, radius})
    }
    pub fn evaluate(self, angle: f64) -> Result<Point2Param, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("parametric circle angle must be finite")); }
        Point2Param::new(self.center.u + self.radius * angle.cos(), self.center.v + self.radius * angle.sin())
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve2 { Line(Line2), Circle(Circle2) }
impl Curve2 {
    pub fn evaluate(self, t: f64) -> Result<Point2Param, GeometryError> {
        match self { Self::Line(l) => l.evaluate(t), Self::Circle(c) => c.evaluate(t) }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneSurface { pub frame: Frame3 }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSurface { pub frame: Frame3, pub radius: f64 }
impl CylinderSurface {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("cylinder radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereSurface { pub frame: Frame3, pub radius: f64 }
impl SphereSurface {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("sphere radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface3 {
    Plane(PlaneSurface),
    Cylinder(CylinderSurface),
    Sphere(SphereSurface),
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceDerivatives {
    pub point: Point3, pub du: Point3, pub dv: Point3,
}
impl Surface3 {
    pub fn domain(self) -> ParameterDomain2 {
        let periodic = ParameterRange::Bounded { start: 0.0, end: TAU };
        match self {
            Self::Plane(_) => ParameterDomain2 { u: ParameterRange::Unbounded, v: ParameterRange::Unbounded, u_periodic:false, v_periodic:false },
            Self::Cylinder(_) => ParameterDomain2 { u: periodic, v: ParameterRange::Unbounded, u_periodic:true, v_periodic:false },
            Self::Sphere(_) => ParameterDomain2 { u: periodic, v: ParameterRange::Bounded { start:-FRAC_PI_2, end:FRAC_PI_2 }, u_periodic:true, v_periodic:false },
        }
    }
    pub fn evaluate(self, u: f64, v: f64) -> Result<SurfaceDerivatives, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("surface parameters must be finite")); }
        let result = match self {
            Self::Plane(surface) => SurfaceDerivatives {
                point: surface.frame.position(u, v)?, du: surface.frame.x, dv: surface.frame.y,
            },
            Self::Cylinder(surface) => {
                let (su, cu) = u.sin_cos();
                let radial = surface.frame.x.scale(cu).add(surface.frame.y.scale(su));
                SurfaceDerivatives {
                    point: surface.frame.origin.add(radial.scale(surface.radius)).add(surface.frame.normal.scale(v)),
                    du: surface.frame.x.scale(-surface.radius * su).add(surface.frame.y.scale(surface.radius * cu)),
                    dv: surface.frame.normal,
                }
            }
            Self::Sphere(surface) => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let radial = surface.frame.x.scale(cu).add(surface.frame.y.scale(su));
                let tangent = surface.frame.x.scale(-su).add(surface.frame.y.scale(cu));
                SurfaceDerivatives {
                    point: surface.frame.origin.add(radial.scale(surface.radius * cv))
                        .add(surface.frame.normal.scale(surface.radius * sv)),
                    du: tangent.scale(surface.radius * cv),
                    dv: radial.scale(-surface.radius * sv)
                        .add(surface.frame.normal.scale(surface.radius * cv)),
                }
            }
        };
        if !result.point.is_finite() || !result.du.is_finite() || !result.dv.is_finite() {
            return Err(invalid_dimension("surface evaluation overflow"));
        }
        Ok(result)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve3Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve2Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SurfaceId(pub u32);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeGeometry {
    pub curve: Curve3Id,
    pub domain: ParameterRange,
    pub tolerance: GeometryTolerance,
}
#[derive(Debug, Clone, PartialEq)]
pub struct FaceGeometry {
    pub surface: SurfaceId,
    /// Projected UV bounds, distinct from the full unbounded analytic plane.
    pub bounds: ParameterDomain2,
    pub outer_loop: LoopId,
    pub tolerance: GeometryTolerance,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoedgeGeometry {
    pub pcurve: Curve2Id,
    pub domain: ParameterRange,
    pub face: FaceId,
    pub loop_id: LoopId,
}

/// Revision-bound B-rep geometry view. The topological identity stays in
/// `Solid`, while analytic carriers, parameter intervals, and p-curves are
/// explicit, independently evaluable entities in this store.
#[derive(Debug, Clone)]
pub struct GeometryStore {
    pub revision: u64,
    pub units: ModelUnits,
    pub curves3: Vec<Curve3>,
    pub curves2: Vec<Curve2>,
    pub surfaces: Vec<Surface3>,
    pub edges: BTreeMap<EdgeId, EdgeGeometry>,
    pub faces: BTreeMap<FaceId, FaceGeometry>,
    pub coedges: BTreeMap<CoedgeId, CoedgeGeometry>,
    body_token: u64,
    vertex_count: usize,
    model_tolerance: GeometryTolerance,
    vertex_tolerances: BTreeMap<VertexId, GeometryTolerance>,
    edge_tolerances: BTreeMap<EdgeId, GeometryTolerance>,
}
impl GeometryStore {
    /// Adapt validated polyhedral B-rep topology to exact planar/linear
    /// carriers, retaining a complete coedge-to-2D-trim mapping.
    pub fn from_solid(solid: &Solid, units: ModelUnits) -> Result<Self, GeometryError> {
        solid.validate()?;
        let token = solid.topology_handle(TopologyEntity::Body)
            .ok_or_else(|| invalid("body lacks identity"))?.body_token;
        let mut store = Self {
            revision: solid.revision, units,
            curves3: Vec::new(), curves2: Vec::new(), surfaces: Vec::new(),
            edges: BTreeMap::new(), faces: BTreeMap::new(), coedges: BTreeMap::new(),
            body_token: token, vertex_count: solid.vertices.len(), model_tolerance: solid.tolerance,
            vertex_tolerances: BTreeMap::new(), edge_tolerances: BTreeMap::new(),
        };
        for edge in &solid.edges {
            let a = solid.vertices[edge.start.0 as usize].position;
            let b = solid.vertices[edge.end.0 as usize].position;
            let (line, domain) = Line3::between(a, b, solid.tolerance)?;
            let curve = Curve3Id(u32::try_from(store.curves3.len())
                .map_err(|_| invalid("too many geometry curves"))?);
            store.curves3.push(Curve3::Line(line));
            store.edges.insert(edge.id, EdgeGeometry { curve, domain, tolerance: solid.tolerance });
        }
        for face in &solid.faces {
            let points: Vec<Point3> = face.boundary.iter()
                .map(|id| solid.vertices[id.0 as usize].position).collect();
            let anchor = points[0];
            let mut sum = Point3 { x:0.0, y:0.0, z:0.0 };
            for i in 0..points.len() {
                sum = sum.add(points[i].sub(anchor).cross(points[(i+1)%points.len()].sub(anchor)));
            }
            let mut axis = None;
            for &p in &points[1..] {
                let diff=p.sub(anchor);
                if length(diff)>solid.tolerance.absolute_length { axis=Some(diff); break; }
            }
            let x = axis.ok_or_else(|| invalid("face has no independent axis"))?;
            let frame = Frame3::from_axes(anchor, x, sum.cross(x), solid.tolerance)?;
            let surface=SurfaceId(u32::try_from(store.surfaces.len())
                .map_err(|_| invalid("too many face surfaces"))?);
            store.surfaces.push(Surface3::Plane(PlaneSurface { frame }));
            let mut uv=Vec::with_capacity(points.len());
            for &p in &points {
                let coords=frame.project(p)?;
                let back=frame.position(coords.u,coords.v)?;
                let extent=solid.bbox.max.sub(solid.bbox.min);
                if distance(back,p)>solid.tolerance.length_at(length(extent)) {
                    return Err(invalid(format!("face {:?} does not project onto supporting plane",face.id)));
                }
                uv.push(coords);
            }
            let umin=uv.iter().map(|p|p.u).fold(f64::INFINITY,f64::min);
            let umax=uv.iter().map(|p|p.u).fold(f64::NEG_INFINITY,f64::max);
            let vmin=uv.iter().map(|p|p.v).fold(f64::INFINITY,f64::min);
            let vmax=uv.iter().map(|p|p.v).fold(f64::NEG_INFINITY,f64::max);
            let bounds=ParameterDomain2 {
                u:ParameterRange::bounded(umin,umax)?,
                v:ParameterRange::bounded(vmin,vmax)?,
                u_periodic:false, v_periodic:false,
            };
            let outer_loop=*face.loops.first().ok_or_else(||invalid("face has no loop"))?;
            store.faces.insert(face.id, FaceGeometry { surface, bounds, outer_loop, tolerance:solid.tolerance });
        }
        for c in &solid.coedges {
            let e=&solid.edges[c.edge.0 as usize];
            let start=solid.vertices[c.start_vertex(e).0 as usize].position;
            let end=solid.vertices[c.end_vertex(e).0 as usize].position;
            let owner=store.faces.get(&c.face).ok_or_else(||invalid("coedge face binding missing"))?;
            let Surface3::Plane(plane) = store.surfaces[owner.surface.0 as usize] else {
                return Err(invalid("v0.7 only projects planar polygonal faces"));
            };
            let a=plane.frame.project(start)?;
            let b=plane.frame.project(end)?;
            let id=Curve2Id(u32::try_from(store.curves2.len()).map_err(|_|invalid("too many pcurves"))?);
            store.curves2.push(Curve2::Line(Line2 { origin:a,
                delta:Point2Param {u:b.u-a.u,v:b.v-a.v} }));
            store.coedges.insert(c.id, CoedgeGeometry {
                pcurve:id, domain:ParameterRange::Bounded {start:0.0,end:1.0},
                face:c.face, loop_id:c.loop_id,
            });
        }
        store.validate_bindings(solid)?;
        Ok(store)
    }

    /// Catches stale geometry views after a CAD transaction or another body.
    pub fn check_current(&self, solid: &Solid) -> Result<(), GeometryError> {
        let token=solid.topology_handle(TopologyEntity::Body)
            .ok_or_else(||invalid("body has no identity"))?.body_token;
        if self.body_token!=token || self.revision!=solid.revision {
            return Err(invalid("geometry view is stale or belongs to a different body"));
        }
        Ok(())
    }
    pub fn set_edge_tolerance(&mut self, edge: EdgeId, tol:GeometryTolerance)
        -> Result<(),GeometryError>
    {
        tol.validate()?;
        if !self.edges.contains_key(&edge){ return Err(invalid("unknown edge tolerance target")); }
        self.edge_tolerances.insert(edge,tol);
        self.edges.get_mut(&edge).expect("checked").tolerance=tol;
        Ok(())
    }
    pub fn set_vertex_tolerance(&mut self, vertex:VertexId, tol:GeometryTolerance)
        -> Result<(),GeometryError>
    {
        tol.validate()?;
        if vertex.0 as usize >= self.vertex_count {
            return Err(invalid("unknown vertex tolerance target"));
        }
        self.vertex_tolerances.insert(vertex,tol);
        Ok(())
    }
    pub fn vertex_tolerance(&self, vertex:VertexId) -> GeometryTolerance {
        self.vertex_tolerances.get(&vertex).copied().unwrap_or(self.model_tolerance)
    }
    pub fn edge_tolerance(&self, edge:EdgeId) -> GeometryTolerance {
        self.edge_tolerances.get(&edge).copied().unwrap_or(self.model_tolerance)
    }
    pub fn face_surface(&self, id:FaceId) -> Option<Surface3> {
        let binding=self.faces.get(&id)?;
        self.surfaces.get(binding.surface.0 as usize).copied()
    }
    pub fn edge_curve(&self, id:EdgeId) -> Option<Curve3> {
        let binding=self.edges.get(&id)?;
        self.curves3.get(binding.curve.0 as usize).copied()
    }
    pub fn trim_curve(&self, id:CoedgeId) -> Option<Curve2> {
        let binding=self.coedges.get(&id)?;
        self.curves2.get(binding.pcurve.0 as usize).copied()
    }
    /// Validate geometry/topology associations and scoped tolerance rules.
    pub fn validate_bindings(&self, solid:&Solid)->Result<(),GeometryError>{
        self.check_current(solid)?;
        if self.edges.len()!=solid.edges.len() || self.faces.len()!=solid.faces.len()
            || self.coedges.len()!=solid.coedges.len(){ return Err(invalid("incomplete analytic B-rep bindings")); }
        for (&id,_) in &self.vertex_tolerances {
            if id.0 as usize >= solid.vertices.len(){return Err(invalid("vertex tolerance references missing vertex"));}
        }
        for e in &solid.edges {
            let carrier=self.edge_curve(e.id).ok_or_else(||invalid("missing edge curve"))?;
            let binding=self.edges.get(&e.id).ok_or_else(||invalid("missing edge binding"))?;
            let a=solid.vertices[e.start.0 as usize].position;
            let b=solid.vertices[e.end.0 as usize].position;
            let edge_len=distance(a,b);
            let mut threshold=self.edge_tolerance(e.id).length_at(edge_len);
            threshold=threshold.max(self.vertex_tolerance(e.start).length_at(edge_len));
            threshold=threshold.max(self.vertex_tolerance(e.end).length_at(edge_len));
            if !edge_len.is_finite() || edge_len<=threshold {
                return Err(invalid(format!("edge {:?} is below scoped tolerance",e.id)));
            }
            let (lo,hi)=match binding.domain {
                ParameterRange::Bounded{start,end}=>(start,end),
                _=>return Err(invalid("topological edge must be parametrically bounded")),
            };
            let t=self.edge_tolerance(e.id).length_at(edge_len);
            if distance(carrier.evaluate(lo)?,a)>t || distance(carrier.evaluate(hi)?,b)>t {
                return Err(invalid(format!("edge {:?} curve misses its end vertices",e.id)));
            }
        }
        for coedge in &solid.coedges {
            let trim=self.trim_curve(coedge.id).ok_or_else(||invalid("missing coedge pcurve"))?;
            let binding=self.coedges.get(&coedge.id).ok_or_else(||invalid("missing trim binding"))?;
            if binding.face!=coedge.face || binding.loop_id!=coedge.loop_id {
                return Err(invalid("trim owner disagrees with B-rep"));
            }
            let surface=self.face_surface(coedge.face).ok_or_else(||invalid("missing supporting surface"))?;
            let edge=&solid.edges[coedge.edge.0 as usize];
            let start=solid.vertices[coedge.start_vertex(edge).0 as usize].position;
            let end=solid.vertices[coedge.end_vertex(edge).0 as usize].position;
            let ParameterRange::Bounded {start:t0,end:t1}=binding.domain else {
                return Err(invalid("coedge trim must be parametrically bounded"));
            };
            for (t,vertex) in [(t0,start),(t1,end)] {
                let uv=trim.evaluate(t)?;
                let back=surface.evaluate(uv.u,uv.v)?.point;
                let eps=self.model_tolerance.length_at(distance(start,end));
                if distance(back,vertex)>eps {
                    return Err(invalid(format!("coedge {:?} pcurve does not map to face boundary",coedge.id)));
                }
            }
        }
        Ok(())
    }
}

impl Solid {
    /// Return a geometry view of the currently validated B-rep; it is
    /// intentionally not cached through subsequent edit revisions.
    pub fn geometry_store(&self, units: ModelUnits) -> Result<GeometryStore, GeometryError> {
        GeometryStore::from_solid(self, units)
    }
}
