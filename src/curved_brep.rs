//! Exact analytic cylinder B-rep with periodic UV seam and a separate facet mesh.
//!
//! A circular edge is a single closed topological edge (start == end),
//! with two opposing coedges. The cylinder's seam has TWO coedges on the
//! SAME lateral face but on distinct periodic UV boundaries (0 and 2*pi).
//! The face geometry is analytic and independent of facet resolution.
//! This first builder supports one complete right circular cylinder only.
use crate::{
    BoundingBox, Circle2, Circle3, Coedge, CoedgeGeometry, CoedgeId, Curve2, Curve2Id,
    Curve3, Curve3Id, CylinderSurface, Edge, EdgeGeometry, EdgeId, FaceGeometry,
    FaceId, Frame3, GeometryError, GeometryTolerance, LengthUnit, Line2,
    Line3, Loop, LoopId, LoopRole, MassProperties, Mesh, ModelUnits, ParameterDomain2,
    ParameterRange, PlaneSurface, Point2Param, Point3, Shell, ShellId, Surface3,
    SurfaceId, Vertex, VertexId,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{PI, TAU};

fn failure(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn dimension(message: &'static str) -> GeometryError { GeometryError::InvalidDimension(message) }
fn point_distance(a: Point3, b: Point3) -> f64 { a.sub(b).length_squared().sqrt() }
fn point_close(a: Point3, b: Point3, eps: f64) -> bool {
    point_distance(a, b) <= eps
}
fn p2(u: f64, v: f64) -> Point2Param { Point2Param { u, v } }
fn trim_line(a: Point2Param, b: Point2Param) -> Curve2 {
    Curve2::Line(Line2 { origin: a, delta: p2(b.u - a.u, b.v - a.v) })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvedFaceKind { BottomCap, TopCap, CylinderSide }

/// `Face` for analytic periodic geometry. Unlike polyhedral `Face`, it does
/// not pretend a circular trim can be represented by a polygon of vertices.
#[derive(Debug, Clone)]
pub struct CurvedFace {
    pub id: FaceId,
    pub kind: CurvedFaceKind,
    pub geometry: FaceGeometry,
    pub loops: Vec<LoopId>,
}

/// Oriented analytic B-rep of a complete right circular cylinder. Public
/// topology is validated before use; the render mesh is derived on demand.
#[derive(Debug, Clone)]
pub struct CylindricalBrep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub coedges: Vec<Coedge>,
    pub loops: Vec<Loop>,
    pub faces: Vec<CurvedFace>,
    pub shells: Vec<Shell>,
    pub curves3: Vec<Curve3>,
    pub curves2: Vec<Curve2>,
    pub surfaces: Vec<Surface3>,
    pub edge_geometry: Vec<EdgeGeometry>,
    pub coedge_geometry: Vec<CoedgeGeometry>,
    pub frame: Frame3,
    pub radius: f64,
    pub height: f64,
    pub tolerance: GeometryTolerance,
    pub units: ModelUnits,
    pub mass: MassProperties,
    pub bbox: BoundingBox,
}

impl CylindricalBrep {
    /// A +Y extrusion from an XZ sketch. The analytic frame is right handed:
    /// frame.x=+X, frame.y=-Z, normal=+Y.
    pub fn upright(origin: Point3, radius: f64, height: f64,
                   tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        let frame = Frame3::from_axes(
            origin,
            Point3 { x: 1.0, y: 0.0, z: 0.0 },
            Point3 { x: 0.0, y: 0.0, z: -1.0 },
            tolerance,
        )?;
        Self::new(frame, radius, height, tolerance)
    }

    /// Build exact cylinder surfaces, cap circles, a shared seam line and six
    /// *oriented* trim curves. The cylinder seam occurs twice on face 2.
    pub fn new(frame: Frame3, radius: f64, height: f64,
               tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !frame.origin.is_finite() || !frame.x.is_finite() ||
            !frame.y.is_finite() || !frame.normal.is_finite() {
            return Err(dimension("non-finite cylinder frame"));
        }
        if !radius.is_finite() || !height.is_finite() || radius <= 0.0 || height <= 0.0 {
            return Err(dimension("cylinder radius and height must be finite and positive"));
        }
        if radius <= tolerance.length_at(radius.max(height)) ||
           height <= tolerance.length_at(radius.max(height)) {
            return Err(dimension("cylinder dimensions below model tolerance"));
        }
        // The frame is public; reject invalid frames rather than relying on
        // an unchecked assumption that it came from Frame3::from_axes().
        for axis in [frame.x, frame.y, frame.normal] {
            let norm = axis.length_squared().sqrt();
            if !norm.is_finite() || (norm - 1.0).abs() > tolerance.angular {
                return Err(dimension("cylinder frame axes must be unit length"));
            }
        }
        if frame.x.dot(frame.y).abs() > tolerance.angular ||
           frame.x.dot(frame.normal).abs() > tolerance.angular ||
           frame.y.dot(frame.normal).abs() > tolerance.angular ||
           frame.x.cross(frame.y).sub(frame.normal).length_squared().sqrt() > tolerance.angular {
            return Err(dimension("cylinder frame must be right-handed and orthogonal"));
        }
        let bottom_origin = frame.origin;
        let top_origin = frame.origin.add(frame.normal.scale(height));
        let b = bottom_origin.add(frame.x.scale(radius));
        let t = top_origin.add(frame.x.scale(radius));
        if !b.is_finite() || !t.is_finite() { return Err(dimension("cylinder coordinates overflow")); }

        // The bottom cap's frame is reversed so a positive UV circle is
        // outward oriented (-N); the top cap has +N orientation.
        let bottom_frame = Frame3::from_axes(bottom_origin, frame.x,
                                             frame.y.scale(-1.0), tolerance)?;
        let top_frame = Frame3 { origin: top_origin, ..frame };
        let side = Surface3::Cylinder(CylinderSurface::new(frame, radius, tolerance)?);
        let surfaces = vec![Surface3::Plane(PlaneSurface {frame:bottom_frame}),
                            Surface3::Plane(PlaneSurface {frame:top_frame}), side];
        let curves3 = vec![
            Curve3::Circle(Circle3::new(frame, radius, tolerance)?),
            Curve3::Circle(Circle3::new(top_frame, radius, tolerance)?),
            Curve3::Line(Line3 {origin:b, direction:frame.normal}),
        ];
        // Each coedge owns its own UV representation. The TWO seam trims
        // are distinct lines (u=TAU and u=0), though they share 3D edge 2.
        let curves2 = vec![
            trim_line(p2(0.0,0.0), p2(TAU,0.0)),
            Curve2::Circle(Circle2::new(p2(0.0,0.0),radius,tolerance)?),
            trim_line(p2(TAU,height), p2(0.0,height)),
            Curve2::Circle(Circle2::new(p2(0.0,0.0),radius,tolerance)?),
            trim_line(p2(TAU,0.0), p2(TAU,height)),
            trim_line(p2(0.0,height), p2(0.0,0.0)),
        ];
        let vertices = vec![
            Vertex{id:VertexId(0),position:b},
            Vertex{id:VertexId(1),position:t},
        ];
        let edges = vec![
            Edge{id:EdgeId(0),start:VertexId(0),end:VertexId(0),coedges:[CoedgeId(0),CoedgeId(1)]},
            Edge{id:EdgeId(1),start:VertexId(1),end:VertexId(1),coedges:[CoedgeId(2),CoedgeId(3)]},
            Edge{id:EdgeId(2),start:VertexId(0),end:VertexId(1),coedges:[CoedgeId(4),CoedgeId(5)]},
        ];
        let coedges = vec![
            // Side's rectangle in UV: bottom U+, seam right V+, top U-, seam left V-.
            Coedge{id:CoedgeId(0),edge:EdgeId(0),face:FaceId(2),loop_id:LoopId(2),reversed:false,
                next:CoedgeId(4),prev:CoedgeId(5),twin:CoedgeId(1)},
            Coedge{id:CoedgeId(1),edge:EdgeId(0),face:FaceId(0),loop_id:LoopId(0),reversed:true,
                next:CoedgeId(1),prev:CoedgeId(1),twin:CoedgeId(0)},
            Coedge{id:CoedgeId(2),edge:EdgeId(1),face:FaceId(2),loop_id:LoopId(2),reversed:true,
                next:CoedgeId(5),prev:CoedgeId(4),twin:CoedgeId(3)},
            Coedge{id:CoedgeId(3),edge:EdgeId(1),face:FaceId(1),loop_id:LoopId(1),reversed:false,
                next:CoedgeId(3),prev:CoedgeId(3),twin:CoedgeId(2)},
            Coedge{id:CoedgeId(4),edge:EdgeId(2),face:FaceId(2),loop_id:LoopId(2),reversed:false,
                next:CoedgeId(2),prev:CoedgeId(0),twin:CoedgeId(5)},
            Coedge{id:CoedgeId(5),edge:EdgeId(2),face:FaceId(2),loop_id:LoopId(2),reversed:true,
                next:CoedgeId(0),prev:CoedgeId(2),twin:CoedgeId(4)},
        ];
        let loops = vec![
            Loop {id:LoopId(0),face:FaceId(0),role:LoopRole::Outer,first_coedge:CoedgeId(1)},
            Loop {id:LoopId(1),face:FaceId(1),role:LoopRole::Outer,first_coedge:CoedgeId(3)},
            Loop {id:LoopId(2),face:FaceId(2),role:LoopRole::Outer,first_coedge:CoedgeId(0)},
        ];
        let full_period=ParameterRange::Bounded {start:0.0,end:TAU};
        let cap_extent=ParameterRange::Bounded {start:-radius,end:radius};
        let faces=vec![
            CurvedFace{id:FaceId(0),kind:CurvedFaceKind::BottomCap,loops:vec![LoopId(0)],
                geometry:FaceGeometry {surface:SurfaceId(0),outer_loop:LoopId(0),tolerance,
                    bounds:ParameterDomain2{u:cap_extent,v:cap_extent,u_periodic:false,v_periodic:false}}},
            CurvedFace{id:FaceId(1),kind:CurvedFaceKind::TopCap,loops:vec![LoopId(1)],
                geometry:FaceGeometry {surface:SurfaceId(1),outer_loop:LoopId(1),tolerance,
                    bounds:ParameterDomain2{u:cap_extent,v:cap_extent,u_periodic:false,v_periodic:false}}},
            CurvedFace{id:FaceId(2),kind:CurvedFaceKind::CylinderSide,loops:vec![LoopId(2)],
                geometry:FaceGeometry {surface:SurfaceId(2),outer_loop:LoopId(2),tolerance,
                    bounds:ParameterDomain2{u:full_period,
                        v:ParameterRange::Bounded {start:0.0,end:height},
                        u_periodic:true,v_periodic:false}}},
        ];
        let edges_geom=vec![
            EdgeGeometry {curve:Curve3Id(0),domain:full_period,tolerance},
            EdgeGeometry {curve:Curve3Id(1),domain:full_period,tolerance},
            EdgeGeometry {curve:Curve3Id(2),domain:ParameterRange::Bounded{start:0.0,end:height},tolerance},
        ];
        let coedge_geometry=(0..6).map(|i|CoedgeGeometry {
            pcurve:Curve2Id(i),
            domain:if i==1 || i==3 {full_period}else{ParameterRange::Bounded{start:0.0,end:1.0}},
            face:coedges[i].face,loop_id:coedges[i].loop_id,
        }).collect();
        let shells=vec![Shell{id:ShellId(0),faces:vec![FaceId(0),FaceId(1),FaceId(2)],closed:true}];
        let area=2.0*PI*radius*radius+TAU*radius*height;
        let volume=PI*radius*radius*height;
        let center=frame.origin.add(frame.normal.scale(height/2.0));
        let bounds=Self::bounds(frame,radius,height);
        let body=Self {vertices,edges,coedges,loops,faces,shells,curves3,curves2,surfaces,
            edge_geometry:edges_geom,coedge_geometry,frame,radius,height,tolerance,
            units:ModelUnits {length:LengthUnit::Millimeter,..ModelUnits::default()},
            mass:MassProperties{volume,surface_area:area,centroid:center},bbox:bounds};
        body.validate()?;
        Ok(body)
    }

    fn bounds(frame:Frame3,radius:f64,height:f64)->BoundingBox {
        let start=frame.origin;
        let end=frame.origin.add(frame.normal.scale(height));
        let amp=|x:f64,y:f64|radius*x.hypot(y);
        BoundingBox{
            min:Point3{x:start.x.min(end.x)-amp(frame.x.x,frame.y.x),
                y:start.y.min(end.y)-amp(frame.x.y,frame.y.y),
                z:start.z.min(end.z)-amp(frame.x.z,frame.y.z)},
            max:Point3{x:start.x.max(end.x)+amp(frame.x.x,frame.y.x),
                y:start.y.max(end.y)+amp(frame.x.y,frame.y.y),
                z:start.z.max(end.z)+amp(frame.x.z,frame.y.z)},
        }
    }

    pub fn euler_characteristic(&self) -> isize {
        self.vertices.len() as isize - self.edges.len() as isize + self.faces.len() as isize
    }
    pub fn side_surface(&self) -> Surface3 { self.surfaces[2] }
    pub fn seam_edge(&self) -> EdgeId { EdgeId(2) }
    pub fn face_surface(&self, face:FaceId) -> Option<Surface3> {
        self.faces.get(face.0 as usize).and_then(|f|self.surfaces.get(f.geometry.surface.0 as usize)).copied()
    }
    pub fn edge_curve(&self, edge:EdgeId) -> Option<Curve3> {
        self.edge_geometry.get(edge.0 as usize).and_then(|e|self.curves3.get(e.curve.0 as usize)).copied()
    }
    pub fn trim_curve(&self, coedge:CoedgeId) -> Option<Curve2> {
        self.coedge_geometry.get(coedge.0 as usize)
            .and_then(|g|self.curves2.get(g.pcurve.0 as usize)).copied()
    }
    /// Validate orientation, radial incidence, UV loop closure and consistency
    /// of 2D trims -> 3D surface -> 3D edge for samples including interiors.
    /// This is deliberately a constrained periodic-cylinder validator, not a
    /// generic NURBS intersection or multi-shell B-rep validator.
    pub fn validate(&self) -> Result<(), GeometryError> {
        self.tolerance.validate()?;
        if self.vertices.len()!=2 || self.edges.len()!=3 || self.coedges.len()!=6 ||
           self.loops.len()!=3 || self.faces.len()!=3 || self.shells.len()!=1 ||
           self.curves3.len()!=3 || self.curves2.len()!=6 || self.surfaces.len()!=3 ||
           self.edge_geometry.len()!=3 || self.coedge_geometry.len()!=6 {
            return Err(failure("cylinder topology and geometry entity counts disagree"));
        }
        if self.radius<=self.tolerance.length_at(self.radius.max(self.height)) ||
           self.height<=self.tolerance.length_at(self.radius.max(self.height)) ||
           !self.radius.is_finite() || !self.height.is_finite() {
            return Err(failure("invalid cylinder dimensions"));
        }
        if !self.shells[0].closed || self.shells[0].id!=ShellId(0) ||
           self.shells[0].faces!=vec![FaceId(0),FaceId(1),FaceId(2)] {
            return Err(failure("inconsistent cylinder shell face ownership"));
        }
        if self.euler_characteristic()!=2 { return Err(failure("invalid cylinder Euler characteristic")); }
        // Geometry evaluation at very large coordinates can lose sub-ULP
        // information. This bound is for *checking evaluated points*, not
        // for accepting sub-resolution features or classifying topology.
        let coordinate_scale=[self.frame.origin.x.abs(),self.frame.origin.y.abs(),
            self.frame.origin.z.abs()].into_iter().fold(1.0f64,f64::max);
        let eps=self.tolerance.length_at((2.0*self.radius).max(self.height))
            .max(16.0*f64::EPSILON*coordinate_scale);
        for (i,v) in self.vertices.iter().enumerate() {
            if v.id.0 as usize!=i || !v.position.is_finite() {
                return Err(failure("invalid cylinder vertex identity or coordinates"));
            }
        }
        for (i,e) in self.edges.iter().enumerate() {
            if e.id.0 as usize!=i || e.start.0 as usize>=self.vertices.len() ||
               e.end.0 as usize>=self.vertices.len() || e.coedges[0]==e.coedges[1] {
                return Err(failure("invalid cylinder edge incidence"));
            }
            let [first,second]=e.coedges;
            if first.0 as usize>=self.coedges.len() || second.0 as usize>=self.coedges.len(){
                return Err(failure("cylinder edge has out-of-bounds coedge"));
            }
            let a=&self.coedges[first.0 as usize];let b=&self.coedges[second.0 as usize];
            if self.edge_geometry[i].curve.0 as usize >= self.curves3.len() {
                return Err(failure("cylinder edge has invalid curve reference"));
            }
            if a.edge!=e.id || b.edge!=e.id || a.twin!=b.id || b.twin!=a.id || a.reversed==b.reversed {
                return Err(failure("cylinder edge coedge twins are inconsistent"));
            }
            if i<2 && e.start!=e.end {return Err(failure("circle edge must be closed"));}
            if i==2 && (e.start!=VertexId(0) || e.end!=VertexId(1)) {
                return Err(failure("cylinder seam edge endpoints invalid"));
            }
            // The periodic seam is intentionally shared by the SAME face.
            if i==2 && (a.face!=FaceId(2) || b.face!=FaceId(2)) {
                return Err(failure("periodic seam must be used twice by lateral face"));
            }
        }
        let expected_loops=[1usize,1usize,4usize];
        let mut usages=vec![0usize;self.coedges.len()];
        for (i,lp) in self.loops.iter().enumerate() {
            if lp.id.0 as usize!=i || lp.face.0 as usize>=self.faces.len() ||
               lp.role!=LoopRole::Outer || lp.first_coedge.0 as usize>=self.coedges.len() ||
               self.faces[lp.face.0 as usize].loops!=vec![lp.id] ||
               self.faces[lp.face.0 as usize].geometry.outer_loop!=lp.id {
                return Err(failure("cylinder face loop ownership mismatch"));
            }
            let start=lp.first_coedge;let mut current=start;let mut seen=BTreeSet::new();
            loop {
                let ix=current.0 as usize;
                if ix>=self.coedges.len() || !seen.insert(ix) {
                    return Err(failure("broken or repeated cylinder coedge loop"));
                }
                let c=&self.coedges[ix];
                if c.id!=current || c.face!=lp.face || c.loop_id!=lp.id ||
                   c.edge.0 as usize>=self.edges.len() || c.next.0 as usize>=self.coedges.len() ||
                   c.prev.0 as usize>=self.coedges.len() {
                    return Err(failure("invalid cylinder coedge owner or pointer"));
                }
                if self.coedges[c.next.0 as usize].prev!=current ||
                   self.coedges[c.prev.0 as usize].next!=current {
                    return Err(failure("cylinder loop next/prev mismatch"));
                }
                let next=&self.coedges[c.next.0 as usize];
                let edge=&self.edges[c.edge.0 as usize];
                let next_edge=&self.edges[next.edge.0 as usize];
                if c.end_vertex(edge)!=next.start_vertex(next_edge){
                    return Err(failure("cylinder loop has disconnected vertices"));
                }
                usages[ix]+=1;
                current=c.next;
                if current==start{break;}
                if seen.len()>self.coedges.len(){return Err(failure("unbounded cylinder loop"));}
            }
            if seen.len()!=expected_loops[i] {return Err(failure("cylinder loop has unexpected coedge count"));}
        }
        if usages.iter().any(|&count|count!=1) {return Err(failure("cylinder coedge orphan or duplicate"));}
        for (i,face) in self.faces.iter().enumerate() {
            if face.id.0 as usize!=i || face.geometry.surface.0 as usize>=self.surfaces.len() ||
               face.kind != [CurvedFaceKind::BottomCap,CurvedFaceKind::TopCap,CurvedFaceKind::CylinderSide][i] {
                return Err(failure("invalid cylindrical face geometry binding"));
            }
        }
        // Full parameter domain as a global trim. Every edge's *interior*
        // must match its parent curve AND the surface evaluation of its pcurve.
        for (i,c) in self.coedges.iter().enumerate() {
            let binding=&self.coedge_geometry[i];
            if binding.face!=c.face || binding.loop_id!=c.loop_id ||
                binding.pcurve.0 as usize>=self.curves2.len() ||
                c.edge.0 as usize>=self.edge_geometry.len() {
                return Err(failure("invalid cylinder trim geometry binding"));
            }
            let e=&self.edges[c.edge.0 as usize];
            let carrier=&self.edge_geometry[c.edge.0 as usize];
            let ParameterRange::Bounded{start:t0,end:t1}=carrier.domain else {
                return Err(failure("unbounded geometric edge carrier"));
            };
            let ParameterRange::Bounded{start:s0,end:s1}=binding.domain else {
                return Err(failure("unbounded cylinder pcurve"));
            };
            let curve=self.curves3[carrier.curve.0 as usize];
            let trim=self.curves2[binding.pcurve.0 as usize];
            if carrier.curve.0 as usize>=self.curves3.len() ||
               c.face.0 as usize>=self.faces.len() ||
               self.faces[c.face.0 as usize].geometry.surface.0 as usize>=self.surfaces.len() {
                return Err(failure("cylinder coedge has out-of-bounds geometry carrier"));
            }
            let surface=self.surfaces[self.faces[c.face.0 as usize].geometry.surface.0 as usize];
            let mut previous=None;
            for fraction in [0.0,0.25,0.5,0.75,1.0] {
                let parameter=s0+(s1-s0)*fraction;
                let uv=trim.evaluate(parameter)?;
                let point=surface.evaluate(uv.u,uv.v)?.point;
                let forward=if c.reversed {1.0-fraction}else{fraction};
                let on_edge=curve.evaluate(t0+(t1-t0)*forward)?;
                if !point_close(point,on_edge,eps){
                    return Err(failure(format!("coedge {} pcurve departs from 3D edge carrier",i)));
                }
                if fraction==0.0 {
                    if !point_close(point,self.vertices[c.start_vertex(e).0 as usize].position,eps) {
                        return Err(failure("cylinder trim start vertex mismatch"));
                    }
                } else if fraction==1.0 {
                    if !point_close(point,self.vertices[c.end_vertex(e).0 as usize].position,eps) {
                        return Err(failure("cylinder trim end vertex mismatch"));
                    }
                }
                if let Some(prev)=previous {
                    // Detect accidentally constant non-degenerate trim curves.
                    if fraction!=1.0 && point_distance(prev,point)<=self.tolerance.absolute_length {
                        return Err(failure("cylinder trim is degenerate"));
                    }
                }
                previous=Some(point);
            }
        }
        // For the side: S(0,v)==S(TAU,v), but the UV bounds remain distinct.
        let side=self.side_surface();
        for v in [0.0,self.height*0.5,self.height] {
            let a=side.evaluate(0.0,v)?.point;
            let b=side.evaluate(TAU,v)?.point;
            if !point_close(a,b,eps){return Err(failure("cylinder periodic seam does not close"));}
        }
        // Outward normals: bottom -N, top +N, side radial (du x dv).
        for (f,uv,expected) in [
            (0,p2(0.0,0.0),self.frame.normal.scale(-1.0)),
            (1,p2(0.0,0.0),self.frame.normal),
            (2,p2(0.0,self.height*0.5),self.frame.x),
        ] {
            let surf=self.face_surface(FaceId(f)).ok_or_else(||failure("missing cylinder surface"))?;
            let deriv=surf.evaluate(uv.u,uv.v)?;
            let normal=deriv.du.cross(deriv.dv);
            let norm=normal.length_squared().sqrt();
            if norm<=0.0 || normal.scale(1.0/norm).dot(expected) < 1.0-self.tolerance.angular {
                return Err(failure("cylinder face has incorrect orientation"));
            }
        }
        let expected_volume=PI*self.radius*self.radius*self.height;
        let expected_area=2.0*PI*self.radius*self.radius+TAU*self.radius*self.height;
        if !self.mass.volume.is_finite() || !self.mass.surface_area.is_finite() ||
           (self.mass.volume-expected_volume).abs()>expected_volume*1.0e-12 ||
           (self.mass.surface_area-expected_area).abs()>expected_area*1.0e-12 {
            return Err(failure("cylinder analytic mass does not match carrier surfaces"));
        }
        let center=self.frame.origin.add(self.frame.normal.scale(self.height/2.0));
        if !point_close(self.mass.centroid,center,eps) {return Err(failure("cylinder centroid inconsistent"));}
        let expected_bbox=Self::bounds(self.frame,self.radius,self.height);
        if !point_close(self.bbox.min,expected_bbox.min,eps) ||
           !point_close(self.bbox.max,expected_bbox.max,eps) {
            return Err(failure("cylinder bounding box inconsistent"));
        }
        Ok(())
    }

    /// Faceting is only for display. None of these triangles replaces the
    /// analytic circle/cylinder geometry, and each triangle retains FaceId.
    pub fn tessellate(&self, segments:usize) -> Result<Mesh, GeometryError> {
        self.validate()?;
        if !(8..=16384).contains(&segments) {
            return Err(dimension("cylinder tessellation segments must be in 8..=16384"));
        }
        let n=segments;
        let mut vertices=Vec::with_capacity(n*2+2);
        for v in [0.0,self.height] {
            for i in 0..n {
                let angle=TAU*i as f64/n as f64;
                vertices.push(self.side_surface().evaluate(angle,v)?.point);
            }
        }
        let bottom_center=vertices.len() as u32;
        vertices.push(self.frame.origin);
        let top_center=vertices.len() as u32;
        vertices.push(self.frame.origin.add(self.frame.normal.scale(self.height)));
        let mut triangles=Vec::with_capacity(4*n);
        let mut triangle_faces=Vec::with_capacity(4*n);
        for i in 0..n {
            let j=(i+1)%n;
            let a=i as u32;
            let b=j as u32;
            let c=(n+i) as u32;
            let d=(n+j) as u32;
            // Bottom is opposite the increasing-angle circle; top agrees.
            triangles.push([bottom_center,b,a]);
            triangle_faces.push(FaceId(0));
            triangles.push([top_center,c,d]);
            triangle_faces.push(FaceId(1));
            triangles.push([a,b,d]);
            triangle_faces.push(FaceId(2));
            triangles.push([a,d,c]);
            triangle_faces.push(FaceId(2));
        }
        let mesh=Mesh {vertices,triangles,triangle_faces};
        self.check_mesh(&mesh)?;
        Ok(mesh)
    }

    pub fn check_mesh(&self,mesh:&Mesh) -> Result<(), GeometryError> {
        if mesh.triangles.len()!=mesh.triangle_faces.len() ||
           mesh.vertices.iter().any(|v|!v.is_finite()) {
            return Err(failure("analytic cylinder mesh has malformed arrays"));
        }
        let mut incidents:BTreeMap<(u32,u32),Vec<bool>>=BTreeMap::new();
        let mut total_six_volume=0.0;
        let anchor=self.frame.origin;
        let mut areas=[0.0;3];
        for (tri,face) in mesh.triangles.iter().zip(&mesh.triangle_faces) {
            if face.0 > 2 || tri.iter().any(|v|*v as usize>=mesh.vertices.len()) ||
               tri[0]==tri[1] || tri[1]==tri[2] || tri[2]==tri[0] {
                return Err(failure("analytic cylinder mesh has invalid triangle"));
            }
            let a=mesh.vertices[tri[0] as usize];
            let b=mesh.vertices[tri[1] as usize];
            let c=mesh.vertices[tri[2] as usize];
            let normal=b.sub(a).cross(c.sub(a));
            let twice=normal.length_squared().sqrt();
            if !twice.is_finite() || twice<=0.0{return Err(failure("degenerate cylinder facet"));}
            areas[face.0 as usize]+=0.5*twice;
            total_six_volume+=a.sub(anchor).dot(b.sub(anchor).cross(c.sub(anchor)));
            for (u,v) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])] {
                incidents.entry((u.min(v),u.max(v))).or_default().push(u<v);
            }
        }
        if incidents.values().any(|uses|uses.len()!=2 || uses[0]==uses[1]) {
            return Err(failure("cylinder display mesh is not closed or consistently wound"));
        }
        let approx_volume=total_six_volume/6.0;
        // Inscribed facets *underestimate* the true analytical volume.
        if !approx_volume.is_finite() || approx_volume<=0.0 ||
           approx_volume>self.mass.volume*(1.0+1.0e-9) {
            return Err(failure("cylinder mesh volume is inconsistent with analytic solid"));
        }
        if areas.iter().any(|x|!x.is_finite() || *x<=0.0) {
            return Err(failure("cylinder mesh must cover all three analytic faces"));
        }
        Ok(())
    }
}
