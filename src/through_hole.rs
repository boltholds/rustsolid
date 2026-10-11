//! First executable *curved* Boolean result: an analytic cylindrical
//! through-hole in a six-face axis-aligned box.
//!
//! Canonical topology retains two true planar inner loops and a single
//! periodic cylindrical side with a seam. No cylinder is converted to facets
//! for B-rep ownership, curves, trims or mass. Meshes are display-only.
//! Parameter edits use an isolated candidate and a bounded undo history.
use crate::{
    BoundingBox, BrepFace, BrepModel, BrepOrigin, BrepSource, Coedge,
    CoedgeGeometry, CoedgeId, Circle2, Circle3, Curve2, Curve2Id,
    Curve3, Curve3Id, CylinderSurface, Edge, EdgeGeometry, EdgeId,
    FaceGeometry, FaceId, Frame3, GeometryError, GeometryTolerance,
    Line2, Line3, Loop, LoopId, LoopRole, MassProperties, Mesh,
    ModelUnits, ParameterDomain2, ParameterRange, PlaneSurface,
    Point2Param, Point3, ShellId, Solid, Surface3, SurfaceId, Vertex,
    VertexId,
};
use std::f64::consts::{PI,TAU};

fn invalid(msg:impl Into<String>)->GeometryError {GeometryError::InvalidEdit(msg.into())}
fn p2(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn uv_line(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:p2(b.u-a.u,b.v-a.v)})
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct HoleParameters {
    pub center_x:f64,
    pub center_z:f64,
    pub radius:f64,
}
#[derive(Debug,Clone)]
pub struct ThroughHoleBrep {
    template:Solid,
    parameters:HoleParameters,
    model:BrepModel,
    revision:u64,
    undo:Vec<HoleParameters>,
    redo:Vec<HoleParameters>,
}
impl ThroughHoleBrep {
    /// Arguments must represent an axis-aligned box and a +Y cylinder,
    /// with identical bottom/top planes. Non-through/coplanar/tangent cuts
    /// remain unsupported rather than being rounded to this primitive.
    pub(crate) fn from_operands(box_:&Solid,cylinder:&crate::CylindricalBrep)
        ->Result<Self,GeometryError>
    {
        box_.validate()?;cylinder.validate()?;
        let n=cylinder.frame.normal;
        let tol=box_.tolerance;
        let h=box_.bbox.max.y-box_.bbox.min.y;
        let eps=tol.length_at(h.max(cylinder.radius*2.0))
            .max(cylinder.tolerance.length_at(h.max(cylinder.radius*2.0)));
        if n.x.abs()>tol.angular || (n.y-1.0).abs()>tol.angular ||
            n.z.abs()>tol.angular {
            return Err(invalid("through-hole requires a +Y cylinder axis"));
        }
        if (cylinder.frame.origin.y-box_.bbox.min.y).abs()>eps ||
            (cylinder.height-h).abs()>eps {
            return Err(invalid("through-hole requires a cylinder spanning exactly the block height"));
        }
        let parameters=HoleParameters{center_x:cylinder.frame.origin.x,
            center_z:cylinder.frame.origin.z,radius:cylinder.radius};
        let model=Self::build(box_,parameters)?;
        Ok(Self{template:box_.clone(),parameters,model,revision:0,
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn from_box(box_:&Solid,parameters:HoleParameters)->Result<Self,GeometryError>{
        // Check box is exactly an axis-aligned six-face cuboid, not an
        // arbitrary solid reinterpreted as its AABB.
        crate::boolean::certify_box(box_).map_err(|e|invalid(e.to_string()))?;
        let model=Self::build(box_,parameters)?;
        Ok(Self{template:box_.clone(),parameters,model,revision:0,
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn parameters(&self)->HoleParameters{self.parameters}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn validate(&self)->Result<(),GeometryError>{
        let expected=Self::build(&self.template,self.parameters)?;
        if self.model.summary()!=expected.summary() ||
            self.model.geometry.curves3!=expected.geometry.curves3 ||
            self.model.geometry.curves2!=expected.geometry.curves2 ||
            self.model.geometry.surfaces!=expected.geometry.surfaces ||
            self.model.vertices.iter().map(|x|x.position).collect::<Vec<_>>()!=
                expected.vertices.iter().map(|x|x.position).collect::<Vec<_>>() {
            return Err(invalid("through-hole parameter/carrier topology disagrees with canonical source"));
        }
        self.model.validate()?;
        Ok(())
    }
    pub fn begin_edit(&mut self)->ThroughHoleTransaction<'_>{
        let staged=self.parameters;
        ThroughHoleTransaction{owner:self,staged,failed:false}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{
        let Some(params)=self.undo.last().copied() else{return Ok(false)};
        let model=Self::build(&self.template,params)?;
        let next=self.revision.checked_add(1).ok_or_else(||invalid("hole revision overflow"))?;
        self.undo.pop();
        self.redo.push(self.parameters);
        self.parameters=params;
        self.model=model;
        self.revision=next;
        Ok(true)
    }
    pub fn redo(&mut self)->Result<bool,GeometryError>{
        let Some(params)=self.redo.last().copied() else{return Ok(false)};
        let model=Self::build(&self.template,params)?;
        let next=self.revision.checked_add(1).ok_or_else(||invalid("hole revision overflow"))?;
        self.redo.pop();
        self.undo.push(self.parameters);
        self.parameters=params;
        self.model=model;
        self.revision=next;
        Ok(true)
    }

    fn build(block:&Solid,parameters:HoleParameters)->Result<BrepModel,GeometryError>{
        let tol=block.tolerance;
        tol.validate()?;
        if !parameters.center_x.is_finite() || !parameters.center_z.is_finite() ||
            !parameters.radius.is_finite() || parameters.radius<=0.0 {
            return Err(invalid("hole center/radius must be finite and radius positive"));
        }
        let bounds=block.bbox;
        let width=bounds.max.x-bounds.min.x;
        let height=bounds.max.y-bounds.min.y;
        let depth=bounds.max.z-bounds.min.z;
        let eps=tol.length_at(width.max(height).max(depth));
        if parameters.radius <= eps ||
            parameters.center_x-parameters.radius <= bounds.min.x+eps ||
            parameters.center_x+parameters.radius >= bounds.max.x-eps ||
            parameters.center_z-parameters.radius <= bounds.min.z+eps ||
            parameters.center_z+parameters.radius >= bounds.max.z-eps {
            return Err(invalid("cylindrical hole must be strictly inside the four block sides"));
        }
        let mut model=BrepModel::from_polyhedral(block,ModelUnits::default())?;
        if model.vertices.len()!=8 || model.edges.len()!=12 ||
           model.coedges.len()!=24 || model.faces.len()!=6 ||
           model.loops.len()!=6 || model.shells.len()!=1 {
            return Err(invalid("through-hole requires the canonical twelve-edge block"));
        }
        let axis_origin=Point3{x:parameters.center_x,y:bounds.min.y,z:parameters.center_z};
        let frame=Frame3::from_axes(axis_origin,
            Point3{x:1.0,y:0.0,z:0.0},Point3{x:0.0,y:0.0,z:-1.0},tol)?;
        let top_origin=Point3{y:bounds.max.y,..axis_origin};
        let top_frame=Frame3{origin:top_origin,..frame};
        let bottom_seam=Point3{x:parameters.center_x+parameters.radius,
            y:bounds.min.y,z:parameters.center_z};
        let top_seam=Point3{y:bounds.max.y,..bottom_seam};
        model.vertices.push(Vertex{id:VertexId(8),position:bottom_seam});
        model.vertices.push(Vertex{id:VertexId(9),position:top_seam});
        let geom=&mut model.geometry;
        let lower_curve=Curve3Id(geom.curves3.len() as u32);
        geom.curves3.push(Curve3::Circle(Circle3::new(frame,parameters.radius,tol)?));
        let upper_curve=Curve3Id(geom.curves3.len() as u32);
        geom.curves3.push(Curve3::Circle(Circle3::new(top_frame,parameters.radius,tol)?));
        let seam_curve=Curve3Id(geom.curves3.len() as u32);
        geom.curves3.push(Curve3::Line(Line3{origin:bottom_seam,direction:frame.normal}));
        model.edges.extend([
            Edge{id:EdgeId(12),start:VertexId(8),end:VertexId(8),
                coedges:[CoedgeId(26),CoedgeId(24)]},
            Edge{id:EdgeId(13),start:VertexId(9),end:VertexId(9),
                coedges:[CoedgeId(27),CoedgeId(25)]},
            Edge{id:EdgeId(14),start:VertexId(8),end:VertexId(9),
                coedges:[CoedgeId(28),CoedgeId(29)]},
        ]);
        geom.edges.extend([
            EdgeGeometry{curve:lower_curve,
                domain:ParameterRange::Bounded{start:0.0,end:TAU},tolerance:tol},
            EdgeGeometry{curve:upper_curve,
                domain:ParameterRange::Bounded{start:0.0,end:TAU},tolerance:tol},
            EdgeGeometry{curve:seam_curve,
                domain:ParameterRange::Bounded{start:0.0,end:height},tolerance:tol},
        ]);
        let cylinder_surface=SurfaceId(geom.surfaces.len() as u32);
        geom.surfaces.push(Surface3::Cylinder(CylinderSurface::new(
            frame,parameters.radius,tol)?));
        // The original block bottom and top local UV frames are opposite.
        // Independent circle2 loops on those planes must follow the 3D edge
        // with reversed orientation at the bottom and forward at the top.
        let lower_plane=match geom.surfaces[model.faces[0].geometry.surface.0 as usize]{
            Surface3::Plane(s)=>s,
            _=>return Err(invalid("block bottom is not a plane")),
        };
        let upper_plane=match geom.surfaces[model.faces[1].geometry.surface.0 as usize]{
            Surface3::Plane(s)=>s,
            _=>return Err(invalid("block top is not a plane")),
        };
        let lower_p=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(Curve2::Circle(Circle2::new(
            lower_plane.frame.project(axis_origin)?,parameters.radius,tol)?));
        let upper_p=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(Curve2::Circle(Circle2::new(
            upper_plane.frame.project(top_origin)?,parameters.radius,tol)?));
        // Side trims copy the exact periodic schema used by CylindricalBrep.
        let bottom_p=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(uv_line(p2(0.0,0.0),p2(TAU,0.0)));
        let top_p=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(uv_line(p2(TAU,height),p2(0.0,height)));
        let seam_a=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(uv_line(p2(TAU,0.0),p2(TAU,height)));
        let seam_b=Curve2Id(geom.curves2.len() as u32);
        geom.curves2.push(uv_line(p2(0.0,height),p2(0.0,0.0)));
        model.coedges.extend([
            Coedge{id:CoedgeId(24),edge:EdgeId(12),face:FaceId(0),
                loop_id:LoopId(6),reversed:true,next:CoedgeId(24),prev:CoedgeId(24),
                twin:CoedgeId(26)},
            Coedge{id:CoedgeId(25),edge:EdgeId(13),face:FaceId(1),
                loop_id:LoopId(7),reversed:false,next:CoedgeId(25),prev:CoedgeId(25),
                twin:CoedgeId(27)},
            Coedge{id:CoedgeId(26),edge:EdgeId(12),face:FaceId(6),
                loop_id:LoopId(8),reversed:false,next:CoedgeId(28),prev:CoedgeId(29),
                twin:CoedgeId(24)},
            Coedge{id:CoedgeId(27),edge:EdgeId(13),face:FaceId(6),
                loop_id:LoopId(8),reversed:true,next:CoedgeId(29),prev:CoedgeId(28),
                twin:CoedgeId(25)},
            Coedge{id:CoedgeId(28),edge:EdgeId(14),face:FaceId(6),
                loop_id:LoopId(8),reversed:false,next:CoedgeId(27),prev:CoedgeId(26),
                twin:CoedgeId(29)},
            Coedge{id:CoedgeId(29),edge:EdgeId(14),face:FaceId(6),
                loop_id:LoopId(8),reversed:true,next:CoedgeId(26),prev:CoedgeId(27),
                twin:CoedgeId(28)},
        ]);
        geom.coedges.extend([
            CoedgeGeometry{pcurve:lower_p,domain:ParameterRange::Bounded{start:0.0,end:TAU},
                face:FaceId(0),loop_id:LoopId(6)},
            CoedgeGeometry{pcurve:upper_p,domain:ParameterRange::Bounded{start:0.0,end:TAU},
                face:FaceId(1),loop_id:LoopId(7)},
            CoedgeGeometry{pcurve:bottom_p,domain:ParameterRange::Bounded{start:0.0,end:1.0},
                face:FaceId(6),loop_id:LoopId(8)},
            CoedgeGeometry{pcurve:top_p,domain:ParameterRange::Bounded{start:0.0,end:1.0},
                face:FaceId(6),loop_id:LoopId(8)},
            CoedgeGeometry{pcurve:seam_a,domain:ParameterRange::Bounded{start:0.0,end:1.0},
                face:FaceId(6),loop_id:LoopId(8)},
            CoedgeGeometry{pcurve:seam_b,domain:ParameterRange::Bounded{start:0.0,end:1.0},
                face:FaceId(6),loop_id:LoopId(8)},
        ]);
        model.loops.extend([
            Loop{id:LoopId(6),face:FaceId(0),role:LoopRole::Inner,first_coedge:CoedgeId(24)},
            Loop{id:LoopId(7),face:FaceId(1),role:LoopRole::Inner,first_coedge:CoedgeId(25)},
            Loop{id:LoopId(8),face:FaceId(6),role:LoopRole::Outer,first_coedge:CoedgeId(26)},
        ]);
        model.faces[0].loops.push(LoopId(6));
        model.faces[1].loops.push(LoopId(7));
        model.faces.push(BrepFace{
            id:FaceId(6),loops:vec![LoopId(8)],
            geometry:FaceGeometry{
                surface:cylinder_surface,outer_loop:LoopId(8),tolerance:tol,
                bounds:ParameterDomain2{
                    u:ParameterRange::Bounded{start:0.0,end:TAU},
                    v:ParameterRange::Bounded{start:0.0,end:height},
                    u_periodic:true,v_periodic:false,
                },
            },
        });
        model.shells[0].faces.push(FaceId(6));
        let removed=PI*parameters.radius*parameters.radius*height;
        let remaining=block.mass.volume-removed;
        if !remaining.is_finite() || remaining<=0.0 {
            return Err(invalid("hole exceeds block volume"));
        }
        let initial=block.mass.centroid;
        model.mass=MassProperties{
            volume:remaining,
            surface_area:block.mass.surface_area-2.0*PI*parameters.radius.powi(2)
                +TAU*parameters.radius*height,
            centroid:Point3 {
                x:initial.x+(initial.x-parameters.center_x)*removed/remaining,
                y:initial.y,
                z:initial.z+(initial.z-parameters.center_z)*removed/remaining,
            },
        };
        model.origin=BrepOrigin::AnalyticThroughHole;
        if model.euler_characteristic()!=0 {
            return Err(invalid("through-hole should have genus-one Euler characteristic"));
        }
        model.validate()?;
        Ok(model)
    }

    /// A facet mesh for CAD display only; NEVER used as Boolean topology.
    pub fn tessellate(&self,segments:usize)->Result<Mesh,GeometryError>{
        if !(12..=4096).contains(&segments){
            return Err(invalid("through-hole mesh needs 12..=4096 circle segments"));
        }
        let b=self.template.bbox;
        let c=self.parameters;
        let mut angles=(0..segments).map(|i|TAU*i as f64/segments as f64).collect::<Vec<_>>();
        for x in [b.min.x,b.max.x] {
            for z in [b.min.z,b.max.z] {
                angles.push((c.center_z-z).atan2(x-c.center_x).rem_euclid(TAU));
            }
        }
        angles.sort_by(f64::total_cmp);
        angles.dedup_by(|a,b|(*a-*b).abs()<1e-12);
        let n=angles.len();
        let mut mesh=Mesh{vertices:Vec::with_capacity(n*4),
            triangles:Vec::with_capacity(n*12),triangle_faces:Vec::with_capacity(n*12)};
        for &u in &angles {
            let (su,cu)=u.sin_cos();
            let rx=if cu>1e-12 {(b.max.x-c.center_x)/cu}
                else if cu< -1e-12 {(b.min.x-c.center_x)/cu}
                else {f64::INFINITY};
            let rz=if su>1e-12 {(c.center_z-b.min.z)/su}
                else if su< -1e-12 {(c.center_z-b.max.z)/su}
                else {f64::INFINITY};
            let t=rx.min(rz);
            let outer=(c.center_x+t*cu,c.center_z-t*su);
            let inner=(c.center_x+c.radius*cu,c.center_z-c.radius*su);
            mesh.vertices.extend([
                Point3{x:outer.0,y:b.min.y,z:outer.1},
                Point3{x:inner.0,y:b.min.y,z:inner.1},
                Point3{x:outer.0,y:b.max.y,z:outer.1},
                Point3{x:inner.0,y:b.max.y,z:inner.1},
            ]);
        }
        for i in 0..n {
            let j=(i+1)%n;
            let a=(4*i) as u32;let b=(4*j) as u32;
            let quads=[
                ([a,b,b+1,a+1],FaceId(0),false),
                ([a+2,b+2,b+3,a+3],FaceId(1),true),
                ([a+1,b+1,b+3,a+3],FaceId(6),false),
            ];
            for (q,face,flip) in quads {
                let mut tri=[[q[0],q[1],q[2]],[q[0],q[2],q[3]]];
                if flip{for t in &mut tri{t.swap(1,2);}}
                mesh.triangles.extend(tri);
                mesh.triangle_faces.extend([face,face]);
            }
            // Outer wall quads belong to one of the four source box sides.
            let aa=mesh.vertices[a as usize];
            let bb=mesh.vertices[b as usize];
            let eps=self.template.tolerance.length_at(
                (self.template.bbox.max.x-self.template.bbox.min.x)
                .max(self.template.bbox.max.z-self.template.bbox.min.z))*64.0;
            let source=self.template.faces.iter().find(|face|{
                face.id.0>=2 && face.boundary.iter().all(|id|{
                    let v=self.template.vertices[id.0 as usize].position;
                    ((v.x-aa.x).abs()<eps && (v.x-bb.x).abs()<eps) ||
                    ((v.z-aa.z).abs()<eps && (v.z-bb.z).abs()<eps)
                })
            }).ok_or_else(||invalid("display outer wall lost source-face ownership"))?.id;
            // For increasing cylinder angle the outward box wall normal is
            // opposite the radial axis: determine winding geometrically.
            let quad=[a,b,b+2,a+2];
            let v0=mesh.vertices[quad[0] as usize];
            let v1=mesh.vertices[quad[1] as usize];
            let v2=mesh.vertices[quad[2] as usize];
            let candidate=v1.sub(v0).cross(v2.sub(v0));
            let center=Point3{x:c.center_x,y:v0.y,z:c.center_z};
            let positive=candidate.dot(v0.sub(center))>0.0;
            let mut tris=[[quad[0],quad[1],quad[2]],[quad[0],quad[2],quad[3]]];
            if !positive{for t in &mut tris{t.swap(1,2);}}
            mesh.triangles.extend(tris);
            mesh.triangle_faces.extend([source,source]);
        }
        Ok(mesh)
    }
}
impl BrepSource for ThroughHoleBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.model.validate()?;
        Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        self.tessellate(segments)
    }
}

pub struct ThroughHoleTransaction<'a>{
    owner:&'a mut ThroughHoleBrep,
    staged:HoleParameters,
    failed:bool,
}
impl ThroughHoleTransaction<'_> {
    pub fn set_radius(&mut self,radius:f64)->Result<(),GeometryError>{
        if self.failed{return Err(invalid("hole transaction already aborted"));}
        if !radius.is_finite()||radius<=0.0{
            self.failed=true;
            return Err(invalid("hole radius must be finite and positive"));
        }
        self.staged.radius=radius;Ok(())
    }
    pub fn set_center(&mut self,x:f64,z:f64)->Result<(),GeometryError>{
        if self.failed{return Err(invalid("hole transaction already aborted"));}
        if !x.is_finite()||!z.is_finite(){
            self.failed=true;return Err(invalid("hole center must be finite"));
        }
        self.staged.center_x=x;self.staged.center_z=z;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.failed {return Err(invalid("cannot commit aborted hole transaction"));}
        let next=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("hole revision overflow"))?;
        // A new complete canonical candidate must validate before changing
        // either model, revision, history or parameters.
        let candidate=ThroughHoleBrep::build(&self.owner.template,self.staged)?;
        if self.owner.undo.len()>=128 {self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters);
        self.owner.redo.clear();
        self.owner.model=candidate;
        self.owner.parameters=self.staged;
        self.owner.revision=next;
        Ok(next)
    }
}
