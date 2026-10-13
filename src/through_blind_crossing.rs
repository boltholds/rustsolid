//! Executable narrow curved Boolean: one upright through cylinder crosses one
//! finite-depth top-entry blind cylinder inside a certified box.
//!
//! The upper aperture follows the *union* of the two source disks. At the
//! blind floor, only the crescent (blind \ through) remains. A separate full
//! lower through-cylinder face connects that floor level to the box bottom.
//! No NURBS or unrestricted curved Boolean execution is implied.

use std::f64::consts::{PI,TAU};
use crate::{
    BlindHoleParameters, BrepFace, BrepModel, BrepOrigin, Circle2,
    Circle3, ClockwiseCircle2, Coedge, CoedgeGeometry, CoedgeId, Curve2,
    Curve2Id, Curve3, Curve3Id, CylinderSurface, Edge, EdgeGeometry, EdgeId,
    FaceGeometry, FaceId, FaceSense, Frame3, GeometryError, HoleParameters,
    IntersectingBlindBrep, Line2, Line3, Loop, LoopId, LoopRole,
    MassProperties, ParameterDomain2, ParameterRange, PlaneSurface,
    Point2Param, Point3, Solid, Surface3, SurfaceId, Vertex, VertexId,
    MixedFootprintContact, MixedFootprintOverlapPlan, MixedFootprintOverlapPlanner,
};

fn invalid(message: impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
fn uv(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn uv_line(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:uv(b.u-a.u,b.v-a.v)})
}
fn bounded(a:f64,b:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(a,b)
}
fn push_uv(m:&mut BrepModel,curve:Curve2)->Result<Curve2Id,GeometryError>{
    let id=Curve2Id(u32::try_from(m.geometry.curves2.len())
        .map_err(|_|invalid("UV curve catalog exhausted"))?);
    m.geometry.curves2.push(curve);Ok(id)
}
fn push_3d(m:&mut BrepModel,curve:Curve3)->Result<Curve3Id,GeometryError>{
    let id=Curve3Id(u32::try_from(m.geometry.curves3.len())
        .map_err(|_|invalid("3D curve catalog exhausted"))?);
    m.geometry.curves3.push(curve);Ok(id)
}
fn angle(v:Point3,frame:Frame3)->f64{
    v.dot(frame.y).atan2(v.dot(frame.x)).rem_euclid(TAU)
}

const MAX_HISTORY:usize=128;
#[derive(Debug,Clone)]
pub struct ThroughBlindCrossingBrep {
    source:Solid,
    through:HoleParameters,
    blind:BlindHoleParameters,
    overlap:MixedFootprintOverlapPlan,
    model:BrepModel,
    revision:u64,
    undo:Vec<(HoleParameters,BlindHoleParameters)>,
    redo:Vec<(HoleParameters,BlindHoleParameters)>,
}
impl ThroughBlindCrossingBrep {
    pub fn from_box(source:&Solid,through:HoleParameters,blind:BlindHoleParameters)
        ->Result<Self,GeometryError>
    {
        let (overlap,model)=Self::rebuild(source,through,blind)?;
        Ok(Self{source:source.clone(),through,blind,overlap,model,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::ThroughBlindCrossingDisplay,GeometryError> {
        crate::ThroughBlindCrossingMesher::tessellate(self,segments)
    }
    pub fn overlap(&self)->&MixedFootprintOverlapPlan{&self.overlap}
    pub fn source_box(&self)->&Solid{&self.source}
    pub fn parameters(&self)->(HoleParameters,BlindHoleParameters){(self.through,self.blind)}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}

    fn rebuild(source:&Solid,through:HoleParameters,blind:BlindHoleParameters)
        ->Result<(MixedFootprintOverlapPlan,BrepModel),GeometryError>
    {
        let overlap=MixedFootprintOverlapPlanner::plan(source,through,blind)?;
        if overlap.contact!=MixedFootprintContact::TwoTransverseCrossings ||
            overlap.overlap_area.is_none() || overlap.exposed_arcs.is_none() {
            return Err(invalid("executable crossing requires two certified transverse circle crossings"));
        }
        let equivalent_through=BlindHoleParameters{
            center_x:through.center_x,center_z:through.center_z,
            radius:through.radius,depth:blind.depth,
        };
        // Reuse an existing *source-proven* upper union opening and its
        // cylinder patch carriers: only the finite-depth lower topology is
        // newly constructed below that certified arrangement.
        let donor=IntersectingBlindBrep::from_box(source,[equivalent_through,blind])?;
        let region=donor.region();
        let arc=region.arc_a;
        let mut m=donor.model().clone();
        let bounds=source.bbox;
        let floor_y=bounds.max.y-blind.depth;
        let lower_height=floor_y-bounds.min.y;
        let eps=source.tolerance.length_at(
            (bounds.max.x-bounds.min.x).max(bounds.max.y-bounds.min.y)
                .max(bounds.max.z-bounds.min.z));
        if lower_height<=eps {return Err(invalid("lower through-cylinder wall has zero height"));}
        let upper=match m.face_surface(FaceId(1)){
            Some(Surface3::Plane(p))=>p,
            _=>return Err(invalid("missing top planar source support")),
        };
        let _donor_floor=match m.face_surface(FaceId(8)){
            Some(Surface3::Plane(p))=>p,
            _=>return Err(invalid("missing existing analytic blind floor support")),
        };
        // Unlike the donor's disk-union floor, the B\T crescent has a
        // clockwise exterior B arc and a counterclockwise interior T arc.
        // Use a +Y support with Forward material sense and mixed-orientation
        // analytic pcurves, so trim winding agrees with actual material.
        let floor_frame=Frame3{
            origin:Point3{y:floor_y,..upper.frame.origin},
            ..upper.frame
        };
        let crescent_floor_surface=SurfaceId(u32::try_from(m.geometry.surfaces.len())
            .map_err(|_|invalid("floor surface catalog exhausted"))?);
        m.geometry.surfaces.push(Surface3::Plane(PlaneSurface{frame:floor_frame}));
        m.faces[8].geometry.surface=crescent_floor_surface;
        m.faces[8].geometry.bounds=m.faces[1].geometry.bounds;
        m.faces[8].sense=FaceSense::Forward;
        let bottom=match m.face_surface(FaceId(0)){
            Some(Surface3::Plane(p))=>p,
            _=>return Err(invalid("missing bottom planar source support")),
        };
        if m.vertices.len()!=12||m.edges.len()!=18||m.coedges.len()!=36 ||
            m.faces.len()!=9||m.loops.len()!=10 {
            return Err(invalid("unexpected donor crossing topology counts"));
        }
        let floor_seam=m.vertices[8].position;
        let base_seam=Point3{y:bounds.min.y,..floor_seam};
        m.vertices.push(Vertex{id:VertexId(12),position:base_seam});
        let lower_origin=Point3{
            x:through.center_x,y:bounds.min.y,z:through.center_z,
        };
        let lower_frame=Frame3{origin:lower_origin,..upper.frame};
        let floor_curve=m.geometry.edges[12].curve;
        let full_circle=push_3d(&mut m,Curve3::Circle(Circle3::new(
            lower_frame,through.radius,source.tolerance)?))?;
        let seam_curve=push_3d(&mut m,Curve3::Line(Line3{
            origin:base_seam,direction:upper.frame.normal,
        }))?;
        let lower_surface=SurfaceId(u32::try_from(m.geometry.surfaces.len())
            .map_err(|_|invalid("surface catalog exhausted"))?);
        m.geometry.surfaces.push(Surface3::Cylinder(CylinderSurface::new(
            lower_frame,through.radius,source.tolerance)?));
        let span_end=arc.parameter_start+TAU;
        let inside=bounded(arc.parameter_end,span_end)?;
        let bottom_circle=bounded(arc.parameter_start,span_end)?;
        let lower_interval=bounded(0.0,lower_height)?;
        m.edges.extend([
            // The floor crescent's inner edge is the physically hidden arc
            // of the through circle, not an arc copied from the exterior.
            Edge{id:EdgeId(18),start:VertexId(9),end:VertexId(8),
                coedges:[CoedgeId(36),CoedgeId(37)]},
            Edge{id:EdgeId(19),start:VertexId(12),end:VertexId(12),
                coedges:[CoedgeId(38),CoedgeId(39)]},
            Edge{id:EdgeId(20),start:VertexId(12),end:VertexId(8),
                coedges:[CoedgeId(40),CoedgeId(41)]},
        ]);
        m.geometry.edges.extend([
            EdgeGeometry{curve:floor_curve,domain:inside,tolerance:source.tolerance},
            EdgeGeometry{curve:full_circle,domain:bottom_circle,tolerance:source.tolerance},
            EdgeGeometry{curve:seam_curve,domain:lower_interval,tolerance:source.tolerance},
        ]);
        // Prior floor coedge 24 becomes a lower-through-wall use of the
        // exterior arc at the blind depth; its twin 28 is the upper wall.
        m.coedges[24].face=FaceId(9);
        m.coedges[24].loop_id=LoopId(11);
        m.coedges[24].next=CoedgeId(40);
        m.coedges[24].prev=CoedgeId(37);
        m.coedges[26].next=CoedgeId(36);
        m.coedges[26].prev=CoedgeId(36);
        m.loops[7].first_coedge=CoedgeId(26);
        m.geometry.coedges[24].face=FaceId(9);
        m.geometry.coedges[24].loop_id=LoopId(11);
        m.geometry.coedges[24].pcurve=push_uv(&mut m,uv_line(
            uv(arc.parameter_end,lower_height),
            uv(arc.parameter_start,lower_height)))?;
        m.geometry.coedges[24].domain=bounded(0.0,1.0)?;

        // The true blind floor is B\\T, with the existing clockwise B
        // exterior arc and a new counterclockwise hidden-T arc. Its
        // clockwise UV carrier is represented analytically rather than
        // approximated by a polyline.
        let floor_center=floor_frame.project(Point3{
            x:through.center_x,y:floor_y,z:through.center_z,
        })?;
        let crescent_pcurve=push_uv(&mut m,Curve2::Circle(
            Circle2::new(floor_center,through.radius,source.tolerance)?))?;
        let b_arc=region.arc_b;
        let b_center=floor_frame.project(Point3{
            x:blind.center_x,y:floor_y,z:blind.center_z,
        })?;
        let b_reverse_uv=push_uv(&mut m,Curve2::ClockwiseCircle(
            ClockwiseCircle2::new(b_center,blind.radius,
                b_arc.parameter_start+b_arc.parameter_end,source.tolerance)?))?;
        m.geometry.coedges[26].pcurve=b_reverse_uv;
        let lower_inside_uv=push_uv(&mut m,uv_line(
            uv(span_end,lower_height),uv(arc.parameter_end,lower_height)))?;
        let bottom_center=bottom.frame.project(lower_origin)?;
        let bottom_phase=angle(upper.frame.x,bottom.frame);
        let bottom_uv=push_uv(&mut m,Curve2::Circle(Circle2::new(
            bottom_center,through.radius,source.tolerance)?))?;
        let lower_circle_uv=push_uv(&mut m,uv_line(
            uv(arc.parameter_start,0.0),uv(span_end,0.0)))?;
        let seam_down_uv=push_uv(&mut m,uv_line(
            uv(span_end,lower_height),uv(span_end,0.0)))?;
        let seam_up_uv=push_uv(&mut m,uv_line(
            uv(arc.parameter_start,0.0),uv(arc.parameter_start,lower_height)))?;
        let extra=[
            // Floor: B exterior (26) -> through interior (36) -> B.
            (18,8,7,false,26,26,37),
            // Lower cylindrical wall: ... 41 -> 37 -> 24 ...
            (18,9,11,true,24,41,36),
            (19,0,10,true,38,38,39),
            (19,9,11,false,41,40,38),
            (20,9,11,true,39,24,41),
            (20,9,11,false,37,39,40),
        ];
        for (i,(edge,face,loop_id,reversed,next,prev,twin)) in extra.into_iter().enumerate() {
            m.coedges.push(Coedge{
                id:CoedgeId(36+i as u32),edge:EdgeId(edge),
                face:FaceId(face),loop_id:LoopId(loop_id),reversed,
                next:CoedgeId(next),prev:CoedgeId(prev),twin:CoedgeId(twin),
            });
        }
        m.geometry.coedges.extend([
            CoedgeGeometry{pcurve:crescent_pcurve,domain:inside,
                face:FaceId(8),loop_id:LoopId(7)},
            CoedgeGeometry{pcurve:lower_inside_uv,domain:bounded(0.0,1.0)?,
                face:FaceId(9),loop_id:LoopId(11)},
            CoedgeGeometry{pcurve:bottom_uv,
                domain:bounded(bottom_phase-span_end,
                    bottom_phase-arc.parameter_start)?,
                face:FaceId(0),loop_id:LoopId(10)},
            CoedgeGeometry{pcurve:lower_circle_uv,domain:bounded(0.0,1.0)?,
                face:FaceId(9),loop_id:LoopId(11)},
            CoedgeGeometry{pcurve:seam_down_uv,domain:bounded(0.0,1.0)?,
                face:FaceId(9),loop_id:LoopId(11)},
            CoedgeGeometry{pcurve:seam_up_uv,domain:bounded(0.0,1.0)?,
                face:FaceId(9),loop_id:LoopId(11)},
        ]);
        m.loops.extend([
            Loop{id:LoopId(10),face:FaceId(0),role:LoopRole::Inner,
                first_coedge:CoedgeId(38)},
            Loop{id:LoopId(11),face:FaceId(9),role:LoopRole::Outer,
                first_coedge:CoedgeId(24)},
        ]);
        m.faces[0].loops.push(LoopId(10));
        m.faces.push(BrepFace{
            id:FaceId(9),loops:vec![LoopId(11)],sense:FaceSense::Reversed,
            geometry:FaceGeometry{
                surface:lower_surface,outer_loop:LoopId(11),
                tolerance:source.tolerance,
                bounds:ParameterDomain2{
                    u:bottom_circle,v:lower_interval,
                    u_periodic:true,v_periodic:false,
                },
            },
        });
        m.shells[0].faces.push(FaceId(9));

        // Independent source-space mass: material removed by the through
        // cutter for H, plus the crescent blind-only area for depth d.
        let area_t=PI*through.radius.powi(2);
        let area_b=PI*blind.radius.powi(2);
        let area_i=overlap.overlap_area.unwrap();
        let crescent_area=area_b-area_i;
        let height=bounds.max.y-bounds.min.y;
        let volume=overlap.predicted_remaining_volume.ok_or_else(||
            invalid("crossed geometry has no verified remaining volume"))?;
        let cylinder_wall=TAU*through.radius*lower_height;
        let exposed_wall=(region.arc_a.angular_span()*through.radius+
            region.arc_b.angular_span()*blind.radius)*blind.depth;
        // Top removes area(T ∪ B), bottom removes area(T), floor restores
        // area(B\\T); net cap area removal = 2 area(T).
        let surface=source.mass.surface_area-2.0*area_t+
            cylinder_wall+exposed_wall;
        if !volume.is_finite()||volume<=0.0||
            !surface.is_finite()||surface<=0.0 ||
            crescent_area<=eps*eps {
            return Err(invalid("through/blind crossing mass is invalid"));
        }
        let old=source.mass.centroid;
        let crescent_moment_x=region.relative_area_moment[0]-
            area_t*(through.center_x-old.x);
        let crescent_moment_z=region.relative_area_moment[1]-
            area_t*(through.center_z-old.z);
        let bx=area_t*height*(through.center_x-old.x)+
            blind.depth*crescent_moment_x;
        let bz=area_t*height*(through.center_z-old.z)+
            blind.depth*crescent_moment_z;
        let cy=(bounds.max.y-blind.depth*0.5)-old.y;
        m.mass=MassProperties{
            volume,surface_area:surface,
            centroid:Point3{
                x:old.x-bx/volume,
                y:old.y-blind.depth*crescent_area*cy/volume,
                z:old.z-bz/volume,
            },
        };
        m.origin=BrepOrigin::AnalyticThroughBlindCrossing;
        m.validate()?;
        let s=m.summary();
        if s.vertices!=13||s.edges!=21||s.coedges!=42||s.loops!=12||
            s.faces!=10||s.shells!=1||s.genus!=1||
            s.euler_characteristic!=0||m.faces[1].loops.len()!=2||
            m.faces[0].loops.len()!=2||m.faces[8].loops.len()!=1 {
            return Err(invalid("through/blind crossing topology/genus certificate failed"));
        }
        let expected=area_t*height+(area_b-area_i)*blind.depth;
        if (source.mass.volume-volume-expected).abs()>
            source.tolerance.area_at(height.max(through.radius*2.0))*32.0 {
            return Err(invalid("through/blind source-volume identity failed"));
        }
        Ok((overlap,m))
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (plan,model)=Self::rebuild(&self.source,self.through,self.blind)?;
        if plan!=self.overlap || model.summary()!=self.model.summary() ||
            model.mass!=self.model.mass || model.origin!=self.model.origin ||
            model.vertices.iter().map(|v|v.position).collect::<Vec<_>>()!=
                self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>() ||
            model.geometry.curves2!=self.model.geometry.curves2 ||
            model.geometry.curves3!=self.model.geometry.curves3 ||
            model.geometry.surfaces!=self.model.geometry.surfaces {
            return Err(invalid("through/blind crossing disagrees with certified rebuild"));
        }
        Ok(())
    }
    pub fn begin_edit(&mut self)->ThroughBlindCrossingTransaction<'_>{
        let through=self.through;
        let blind=self.blind;
        ThroughBlindCrossingTransaction{
            owner:self,through,blind,poisoned:false,
        }
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let next=if redo{self.redo.last().copied()}else{self.undo.last().copied()};
        let Some((through,blind))=next else{return Ok(false)};
        let (overlap,model)=Self::rebuild(&self.source,through,blind)?;
        let revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("crossing revision overflow"))?;
        if redo{self.redo.pop();self.undo.push((self.through,self.blind));}
        else{self.undo.pop();self.redo.push((self.through,self.blind));}
        self.through=through;self.blind=blind;self.overlap=overlap;
        self.model=model;self.revision=revision;
        Ok(true)
    }
}
impl crate::BrepSource for ThroughBlindCrossingBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<crate::Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct ThroughBlindCrossingTransaction<'a>{
    owner:&'a mut ThroughBlindCrossingBrep,
    through:HoleParameters,
    blind:BlindHoleParameters,
    poisoned:bool,
}
impl ThroughBlindCrossingTransaction<'_>{
    pub fn replace_through(&mut self,next:HoleParameters)->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("crossing transaction poisoned"));}
        if !next.center_x.is_finite()||!next.center_z.is_finite()||
            !next.radius.is_finite()||next.radius<=0.0{
            self.poisoned=true;return Err(invalid("invalid staged through cylinder"));
        }
        self.through=next;Ok(())
    }
    pub fn replace_blind(&mut self,next:BlindHoleParameters)->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("crossing transaction poisoned"));}
        if !next.center_x.is_finite()||!next.center_z.is_finite()||
            !next.radius.is_finite()||next.radius<=0.0||
            !next.depth.is_finite()||next.depth<=0.0{
            self.poisoned=true;return Err(invalid("invalid staged blind cylinder"));
        }
        self.blind=next;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned crossing transaction"));}
        let (overlap,model)=ThroughBlindCrossingBrep::rebuild(
            &self.owner.source,self.through,self.blind)?;
        let revision=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("crossing revision overflow"))?;
        if self.owner.undo.len()>=MAX_HISTORY{self.owner.undo.remove(0);}
        self.owner.undo.push((self.owner.through,self.owner.blind));
        self.owner.redo.clear();
        self.owner.through=self.through;self.owner.blind=self.blind;
        self.owner.overlap=overlap;self.owner.model=model;
        self.owner.revision=revision;
        Ok(revision)
    }
}
