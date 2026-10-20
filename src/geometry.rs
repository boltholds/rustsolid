//! B-rep-to-analytic-carrier projection; the geometry types live in rustsolid-geometry.
//! These adapters depend on topology and must never be moved into the pure geometry crate.
use crate::{CoedgeId, EdgeId, FaceId, GeometryError, GeometryTolerance, LoopId, Point3, Solid, TopologyEntity, VertexId};
use std::collections::BTreeMap;
pub use rustsolid_geometry::{
    LengthUnit, AngleUnit, ModelUnits, ParameterRange, ParameterDomain2, Point2Param,
    Frame3, Line3, Circle3, HarmonicCurve3, Curve3, Line2, Circle2,
    HarmonicCurve2, CylinderSection2, Curve2,
    PlaneSurface, CylinderSurface, SphereSurface, Surface3, SurfaceDerivatives, OrientedCircle2,
    Curve3Id, Curve2Id, SurfaceId,
};

fn invalid(reason: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(reason.into())
}
fn length(p: Point3) -> f64 { p.x.hypot(p.y).hypot(p.z) }
fn distance(a: Point3, b: Point3) -> f64 { length(a.sub(b)) }

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
            // Use normalized directions before constructing the second axis.
            // Previously sum.cross(x) multiplied a tiny face-area vector by
            // a tiny first edge, underflowing the angular frame tolerance
            // even though the B-rep face itself was geometrically valid.
            let x_length=length(x);
            let normal_length=length(sum);
            if !x_length.is_finite() || !normal_length.is_finite()
                || x_length<=solid.tolerance.absolute_length || normal_length<=0.0 {
                return Err(invalid("face frame cannot be resolved"));
            }
            let unit_x=x.scale(1.0/x_length);
            let unit_normal=sum.scale(1.0/normal_length);
            let frame=Frame3::from_axes(anchor,unit_x,unit_normal.cross(unit_x),solid.tolerance)?;
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
