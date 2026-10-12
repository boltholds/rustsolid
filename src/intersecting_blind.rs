//! First *executable* crossing-trim arrangement: the union of two partially
//! overlapping, strictly vertical, equal-depth blind cylinders cut from a
//! certified axis-aligned box. This is NOT a general curved Boolean engine.
//!
//! The two cylinder/cylinder intersection branches are vertical lines. They
//! split four 3D circle edges into finite exposed arcs, two transverse seam
//! edges, a single two-arc opening loop and an independent two-arc floor.
//! Both cylindrical patches are real trimmed analytical surfaces; no polygon
//! approximation is used for their shared B-rep topology or analytical mass.
use crate::{
    BlindHoleParameters, BrepFace, BrepGeometry, BrepModel, BrepOrigin,
    Circle2, Circle3, Coedge, CoedgeGeometry, CoedgeId, Curve2, Curve2Id,
    Curve3, Curve3Id, CylinderSurface, Edge, EdgeGeometry, EdgeId,
    FaceGeometry, FaceId, FaceSense, Frame3, GeometryError, Line2, Line3,
    Loop, LoopId, LoopRole, MassProperties, ModelUnits, ParameterDomain2,
    ParameterRange, PlaneSurface, Point2Param, Point3, Solid, Surface3,
    SurfaceId, Vertex, VertexId,
};
use std::f64::consts::{PI,TAU};

fn invalid(msg:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(msg.into())}
fn norm(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn coord(p:Point3)->(f64,f64){(p.x,p.z)}
fn uv(a:f64,b:f64)->Point2Param{Point2Param{u:a,v:b}}
fn uv_line(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:uv(b.u-a.u,b.v-a.v)})
}
fn angle(v:Point3,frame:Frame3)->f64{
    v.dot(frame.y).atan2(v.dot(frame.x)).rem_euclid(TAU)
}
fn range_for(start:f64,end:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(start,end)
}
fn fresh_index(len:usize)->Result<u32,GeometryError>{
    u32::try_from(len).map_err(|_|invalid("analytical arrangement carrier catalog exhausted"))
}
fn add_curve3(model:&mut BrepModel,curve:Curve3)->Result<Curve3Id,GeometryError>{
    let id=Curve3Id(fresh_index(model.geometry.curves3.len())?);
    model.geometry.curves3.push(curve);Ok(id)
}
fn add_curve2(model:&mut BrepModel,curve:Curve2)->Result<Curve2Id,GeometryError>{
    let id=Curve2Id(fresh_index(model.geometry.curves2.len())?);
    model.geometry.curves2.push(curve);Ok(id)
}
fn add_surface(model:&mut BrepModel,surface:Surface3)->Result<SurfaceId,GeometryError>{
    let id=SurfaceId(fresh_index(model.geometry.surfaces.len())?);
    model.geometry.surfaces.push(surface);Ok(id)
}
fn point_on_circle(p:BlindHoleParameters,frame:Frame3,u:f64,y:f64)->Point3{
    Point3{x:p.center_x,y,z:p.center_z}
        .add(frame.x.scale(p.radius*u.cos()))
        .add(frame.y.scale(p.radius*u.sin()))
}
fn circular_segment_area(radius:f64,chord_x:f64,chord_height:f64)->f64 {
    let q=(chord_x/radius).clamp(-1.0,1.0);
    radius*radius*q.acos()-chord_x*chord_height
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ExposedCircularArc {
    pub cutter_index:usize,
    pub parameter_start:f64,
    pub parameter_end:f64,
    pub start_crossing:usize,
    pub end_crossing:usize,
    pub radius:f64,
}
impl ExposedCircularArc {
    pub fn angular_span(self)->f64{self.parameter_end-self.parameter_start}
    pub fn arc_length(self)->f64{self.angular_span()*self.radius}
}
/// Verified planar set arrangement, with exactly two *transverse* circle
/// intersections. Each exposed arc is a true parametrically bounded circle,
/// not a tessellated segment. The union of two disk interiors is one region.
#[derive(Debug,Clone)]
pub struct CrossedDiskRegion {
    pub arc_a:ExposedCircularArc,
    pub arc_b:ExposedCircularArc,
    pub crossings_xz:[[f64;2];2],
    pub disk_union_area:f64,
    pub exposed_perimeter:f64,
    /// ∫_union (X - boxCentroid) dA, projected onto world X/Z.
    pub relative_area_moment:[f64;2],
    pub cut_depth:f64,
}
impl CrossedDiskRegion {
    pub fn from_box(box_:&Solid,holes:[BlindHoleParameters;2])
        ->Result<Self,GeometryError>
    {
        crate::boolean::certify_box(box_).map_err(|e|invalid(e.to_string()))?;
        let b=box_.bbox;
        let span=norm(b.max.sub(b.min));
        let eps=box_.tolerance.length_at(span);
        for h in holes {
            // Reuse the existing independently source-validated top opening
            // and true circular blind floor contract.
            let donor=crate::BlindHoleBrep::from_box(box_,h)?;
            donor.opening_trace()?;
        }
        let [a,c]=holes;
        let depth_budget=16.0*f64::EPSILON*
            [a.depth.abs(),c.depth.abs(),b.max.y.abs()].into_iter()
            .fold(1.0_f64,f64::max);
        if (a.depth-c.depth).abs()>depth_budget {
            return Err(invalid("crossing blind cuts require mathematically equal floor depths"));
        }
        let depth=a.depth;
        let dx=c.center_x-a.center_x;
        let dz=c.center_z-a.center_z;
        let d=dx.hypot(dz);
        if !d.is_finite() || d<=eps || d<=(a.radius-c.radius).abs()+eps ||
            d>=a.radius+c.radius-eps {
            return Err(invalid("two blind disks must partially overlap at two transverse crossings, without tangency, containment or coincidence"));
        }
        let axis=Point3{x:dx/d,y:0.0,z:dz/d};
        let tangent=Point3{x:-axis.z,y:0.0,z:axis.x};
        let x_a=(d*d+a.radius*a.radius-c.radius*c.radius)/(2.0*d);
        let x_b=d-x_a;
        if x_a.abs()>=a.radius || x_b.abs()>=c.radius {
            return Err(invalid("circle chord parameters exceed real disk boundary"));
        }
        let chord=(a.radius*a.radius-x_a*x_a).max(0.0).sqrt();
        let chord_b=(c.radius*c.radius-x_b*x_b).max(0.0).sqrt();
        if chord<=eps*32.0 || (chord-chord_b).abs()>eps*4.0 {
            return Err(invalid("circle crossing is nearly tangent or numerically unresolved"));
        }
        let center=Point3{x:a.center_x+x_a*axis.x,
            y:b.max.y,z:a.center_z+x_a*axis.z};
        let intersections=[
            center.add(tangent.scale(chord)),
            center.sub(tangent.scale(chord)),
        ];
        for p in intersections {
            for h in holes {
                let err=(p.x-h.center_x).hypot(p.z-h.center_z)-h.radius;
                if err.abs()>eps*4.0 {
                    return Err(invalid("circle-crossing witness leaves an analytical source carrier"));
                }
            }
        }
        let top=match crate::BrepModel::from_polyhedral(box_,ModelUnits::default())?
            .face_surface(FaceId(1)){
            Some(Surface3::Plane(surface))=>surface,
            _=>return Err(invalid("top cut source must be a real planar face"))
        };
        let upper=top.frame;
        if (upper.normal.y-1.0).abs()>16.0*f64::EPSILON ||
            upper.x.y.abs()>16.0*f64::EPSILON ||
            upper.y.y.abs()>16.0*f64::EPSILON {
            return Err(invalid("only the true upward top plane supports crossing blind cuts"));
        }
        let arcs=[a,c].into_iter().enumerate().map(|(index,h)|{
            let other=holes[1-index];
            let angles=intersections.map(|point|{
                angle(point.sub(Point3{x:h.center_x,y:point.y,z:h.center_z}),upper)
            });
            // Select the forward physical arc whose interior does NOT lie
            // in the opposite disk; the hidden arc is never emitted as B-rep.
            let candidate=|start:usize,end:usize|->ExposedCircularArc{
                let u0=angles[start];
                let u1=u0+(angles[end]-u0).rem_euclid(TAU);
                ExposedCircularArc{cutter_index:index,parameter_start:u0,
                    parameter_end:u1,start_crossing:start,end_crossing:end,
                    radius:h.radius}
            };
            let c0=candidate(0,1);let c1=candidate(1,0);
            let exposed=|arc:ExposedCircularArc|{
                let mid=(arc.parameter_start+arc.parameter_end)/2.0;
                let sample=point_on_circle(h,upper,mid,b.max.y);
                (sample.x-other.center_x).hypot(sample.z-other.center_z)
            };
            let good0=exposed(c0)>other.radius+eps;
            let good1=exposed(c1)>other.radius+eps;
            if good0==good1 {
                Err(invalid("unable to classify the two exposed circle arc interiors"))
            }else{Ok(if good0{c0}else{c1})}
        }).collect::<Result<Vec<_>,_>>()?;
        let arc_a=arcs[0];let arc_b=arcs[1];
        if arc_a.start_crossing!=arc_b.end_crossing ||
           arc_a.end_crossing!=arc_b.start_crossing {
            return Err(invalid("crossing arcs do not close into a single oriented top trim loop"));
        }
        for (arc,h,other) in [(arc_a,a,c),(arc_b,c,a)] {
            if arc.angular_span()<=box_.tolerance.angular*16.0 ||
               arc.angular_span()>=TAU-box_.tolerance.angular*16.0 {
                return Err(invalid("exposed analytic arc degenerates at the circle crossing"));
            }
            for fraction in [0.125,0.25,0.5,0.75,0.875] {
                let u=arc.parameter_start+arc.angular_span()*fraction;
                let p=point_on_circle(h,upper,u,b.max.y);
                if (p.x-other.center_x).hypot(p.z-other.center_z)<=other.radius+eps {
                    return Err(invalid("visible circle arc enters hidden overlapping disk region"));
                }
            }
        }
        let seg_a=circular_segment_area(a.radius,x_a,chord);
        let seg_b=circular_segment_area(c.radius,x_b,chord);
        let union=PI*(a.radius*a.radius+c.radius*c.radius)-seg_a-seg_b;
        let moment_a=(2.0/3.0)*chord.powi(3);
        let moment_b=(2.0/3.0)*chord.powi(3);
        let bx=box_.mass.centroid.x;let bz=box_.mass.centroid.z;
        // Disk 1 + disk 2 minus their overlapping lens. The lens is the
        // disjoint union of two chord-delimited circular segments.
        let relative=[
            (PI*a.radius*a.radius-seg_a)*(a.center_x-bx)+
                (PI*c.radius*c.radius-seg_b)*(c.center_x-bx)-
                (moment_a-moment_b)*axis.x,
            (PI*a.radius*a.radius-seg_a)*(a.center_z-bz)+
                (PI*c.radius*c.radius-seg_b)*(c.center_z-bz)-
                (moment_a-moment_b)*axis.z,
        ];
        let perimeter=arc_a.arc_length()+arc_b.arc_length();
        let box_area=(b.max.x-b.min.x)*(b.max.z-b.min.z);
        if !union.is_finite() || union<=eps*eps ||
            union>=box_area-eps*eps || !perimeter.is_finite() ||
            perimeter<=0.0 || relative.iter().any(|x|!x.is_finite()) {
            return Err(invalid("closed crossing-disk arrangement has invalid mass moments"));
        }
        Ok(Self{arc_a,arc_b,
            crossings_xz:intersections.map(|p|[p.x,p.z]),
            disk_union_area:union,exposed_perimeter:perimeter,
            relative_area_moment:relative,cut_depth:depth})
    }
}

#[derive(Debug,Clone)]
pub struct IntersectingBlindBrep {
    source:Solid,
    parameters:[BlindHoleParameters;2],
    region:CrossedDiskRegion,
    model:BrepModel,
    revision:u64,
    undo:Vec<[BlindHoleParameters;2]>,
    redo:Vec<[BlindHoleParameters;2]>,
}
impl IntersectingBlindBrep {
    pub fn from_box(source:&Solid,parameters:[BlindHoleParameters;2])
        ->Result<Self,GeometryError>
    {
        let region=CrossedDiskRegion::from_box(source,parameters)?;
        let model=Self::assemble(source,parameters,&region)?;
        Ok(Self{source:source.clone(),parameters,region,model,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn parameters(&self)->[BlindHoleParameters;2]{self.parameters}
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn region(&self)->&CrossedDiskRegion{&self.region}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let region=CrossedDiskRegion::from_box(&self.source,self.parameters)?;
        let independent=Self::assemble(&self.source,self.parameters,&region)?;
        if self.model.summary()!=independent.summary() ||
           self.model.mass!=independent.mass ||
           self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>()!=
                independent.vertices.iter().map(|v|v.position).collect::<Vec<_>>() ||
           self.model.geometry.curves3!=independent.geometry.curves3 ||
           self.model.geometry.curves2!=independent.geometry.curves2 ||
           self.model.geometry.surfaces!=independent.geometry.surfaces ||
           self.model.faces.iter().map(|f|f.sense).collect::<Vec<_>>()!=
                independent.faces.iter().map(|f|f.sense).collect::<Vec<_>>() {
            return Err(invalid("crossed-circle B-rep disagrees with independently reconstructed source arrangement"));
        }
        Ok(())
    }
    pub fn begin_edit(&mut self)->IntersectingBlindTransaction<'_>{
        let staged=self.parameters;
        IntersectingBlindTransaction{owner:self,staged,poisoned:false}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let target=if redo{&self.redo}else{&self.undo};
        let Some(next)=target.last().copied()else{return Ok(false)};
        let region=CrossedDiskRegion::from_box(&self.source,next)?;
        let model=Self::assemble(&self.source,next,&region)?;
        let next_revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("crossed-circle edit revision overflow"))?;
        if redo{self.redo.pop();self.undo.push(self.parameters);}
        else{self.undo.pop();self.redo.push(self.parameters);}
        self.parameters=next;
        self.region=region;
        self.model=model;
        self.revision=next_revision;
        Ok(true)
    }
    fn assemble(source:&Solid,holes:[BlindHoleParameters;2],region:&CrossedDiskRegion)
        ->Result<BrepModel,GeometryError>
    {
        let mut m=BrepModel::from_polyhedral(source,ModelUnits::default())?;
        if m.vertices.len()!=8 || m.edges.len()!=12 || m.coedges.len()!=24 ||
            m.loops.len()!=6 || m.faces.len()!=6 || m.shells.len()!=1 {
            return Err(invalid("crossed-circle assembly requires six canonical box faces"));
        }
        let top=match m.face_surface(FaceId(1)){
            Some(Surface3::Plane(p))=>p,
            _=>return Err(invalid("source top face lacks an analytic plane")),
        };
        let bottom=match m.face_surface(FaceId(0)){
            Some(Surface3::Plane(p))=>p,
            _=>return Err(invalid("source bottom face lacks an analytic plane")),
        };
        if (top.frame.normal.y-1.0).abs()>16.0*f64::EPSILON ||
            (bottom.frame.normal.y+1.0).abs()>16.0*f64::EPSILON {
            return Err(invalid("crossed-circle builder requires top +Y and bottom -Y oriented cap supports"));
        }
        let top_y=source.bbox.max.y;
        let floor_y=top_y-region.cut_depth;
        let bottom_frame=Frame3{origin:Point3{y:floor_y,..bottom.frame.origin},..bottom.frame};
        let floor_sid=add_surface(&mut m,Surface3::Plane(PlaneSurface{frame:bottom_frame}))?;
        let arcs=[region.arc_a,region.arc_b];
        let junctions=region.crossings_xz;
        let lo=arcs[0].start_crossing;
        let hi=arcs[0].end_crossing;
        for (index,crossing) in [lo,hi].into_iter().enumerate(){
            let point=junctions[crossing];
            m.vertices.push(Vertex{id:VertexId(8+index as u32),
                position:Point3{x:point[0],y:floor_y,z:point[1]}});
        }
        for (index,crossing) in [lo,hi].into_iter().enumerate(){
            let point=junctions[crossing];
            m.vertices.push(Vertex{id:VertexId(10+index as u32),
                position:Point3{x:point[0],y:top_y,z:point[1]}});
        }
        let mut side_surfaces=[SurfaceId(0);2];
        let mut circle3=[Curve3Id(0);4]; // bottom A, top A, bottom B, top B
        let mut top_uv=[Curve2Id(0);2];
        let mut floor_uv=[Curve2Id(0);2];
        let mut side_uv=[Curve2Id(0);8]; // bottom,top,low,high × two sides
        let mut circle_domains=[ParameterRange::Unbounded;2];
        for (i,arc) in arcs.iter().enumerate(){
            let p=holes[i];
            let base=Frame3{origin:Point3{x:p.center_x,y:floor_y,z:p.center_z},
                ..top.frame};
            let at_top=Frame3{origin:Point3{y:top_y,..base.origin},..base};
            side_surfaces[i]=add_surface(&mut m,Surface3::Cylinder(
                CylinderSurface::new(base,p.radius,source.tolerance)?))?;
            circle3[2*i]=add_curve3(&mut m,Curve3::Circle(Circle3::new(
                base,p.radius,source.tolerance)?))?;
            circle3[2*i+1]=add_curve3(&mut m,Curve3::Circle(Circle3::new(
                at_top,p.radius,source.tolerance)?))?;
            circle_domains[i]=range_for(arc.parameter_start,arc.parameter_end)?;
            top_uv[i]=add_curve2(&mut m,Curve2::Circle(Circle2::new(
                top.frame.project(at_top.origin)?,p.radius,source.tolerance)?))?;
            let phase=angle(base.x,bottom_frame);
            let floor_curve=Circle2::new(bottom_frame.project(base.origin)?,
                p.radius,source.tolerance)?;
            floor_uv[i]=add_curve2(&mut m,Curve2::Circle(floor_curve))?;
            // The floor face has -Y parameterization but +Y material
            // normal. Its directed coedge uses the reverse of the 3D edge.
            // Phase phi = phase-u reconciles the opposite planar UV frames.
            let _=phase;
            let side_start=arc.parameter_start;
            let side_end=arc.parameter_end;
            side_uv[4*i]=add_curve2(&mut m,uv_line(
                uv(side_start,0.0),uv(side_end,0.0)))?;
            side_uv[4*i+1]=add_curve2(&mut m,uv_line(
                uv(side_end,region.cut_depth),uv(side_start,region.cut_depth)))?;
            side_uv[4*i+2]=add_curve2(&mut m,uv_line(
                uv(side_start,region.cut_depth),uv(side_start,0.0)))?;
            side_uv[4*i+3]=add_curve2(&mut m,uv_line(
                uv(side_end,0.0),uv(side_end,region.cut_depth)))?;
        }
        let vertical_low=add_curve3(&mut m,Curve3::Line(Line3 {
            origin:m.vertices[8].position,direction:Point3{x:0.0,y:1.0,z:0.0},
        }))?;
        let vertical_high=add_curve3(&mut m,Curve3::Line(Line3 {
            origin:m.vertices[9].position,direction:Point3{x:0.0,y:1.0,z:0.0},
        }))?;
        m.edges.extend([
            Edge{id:EdgeId(12),start:VertexId(8),end:VertexId(9),
                coedges:[CoedgeId(24),CoedgeId(28)]},
            Edge{id:EdgeId(13),start:VertexId(10),end:VertexId(11),
                coedges:[CoedgeId(25),CoedgeId(29)]},
            Edge{id:EdgeId(14),start:VertexId(9),end:VertexId(8),
                coedges:[CoedgeId(26),CoedgeId(30)]},
            Edge{id:EdgeId(15),start:VertexId(11),end:VertexId(10),
                coedges:[CoedgeId(27),CoedgeId(31)]},
            Edge{id:EdgeId(16),start:VertexId(8),end:VertexId(10),
                coedges:[CoedgeId(32),CoedgeId(33)]},
            Edge{id:EdgeId(17),start:VertexId(9),end:VertexId(11),
                coedges:[CoedgeId(34),CoedgeId(35)]},
        ]);
        for i in 0..2{
            for v in 0..2 {
                m.geometry.edges.push(EdgeGeometry {
                    curve:circle3[2*i+v],
                    domain:circle_domains[i],
                    tolerance:source.tolerance,
                });
            }
        }
        for curve in [vertical_low,vertical_high] {
            m.geometry.edges.push(EdgeGeometry{
                curve,domain:range_for(0.0,region.cut_depth)?,
                tolerance:source.tolerance,
            });
        }
        // Top cap: A forward -> B forward (one inner loop).
        // Floor: A reverse -> B reverse (one outer two-arc loop).
        // Each side has four edges: lower arc, vertical crossing,
        // reversed upper arc, and the opposite vertical crossing.
        let tuples=[
            (12,8,7,true,26,26,28),
            (13,1,6,false,27,27,29),
            (14,8,7,true,24,24,30),
            (15,1,6,false,25,25,31),
            (12,6,8,false,34,32,24),
            (13,6,8,true,32,34,25),
            (14,7,9,false,33,35,26),
            (15,7,9,true,35,33,27),
            (16,6,8,true,28,29,33),
            (16,7,9,false,31,30,32),
            (17,6,8,false,29,28,35),
            (17,7,9,true,30,31,34),
        ];
        for (i,(edge,face,lp,reversed,next,prev,twin)) in tuples.into_iter().enumerate(){
            m.coedges.push(Coedge{id:CoedgeId(24+i as u32),
                edge:EdgeId(edge),face:FaceId(face),loop_id:LoopId(lp),
                reversed,next:CoedgeId(next),prev:CoedgeId(prev),twin:CoedgeId(twin)});
        }
        // The exact B-rep 3D circle is parameterized in the top-plane
        // orientation. The bottom support normal is reversed; floor pcurve
        // angles therefore increase as the 3D arc parameter decreases.
        for (i,arc) in arcs.iter().enumerate(){
            let plane_phase=angle(top.frame.x,bottom_frame);
            let u_start=plane_phase-arc.parameter_end;
            let u_end=plane_phase-arc.parameter_start;
            m.geometry.coedges.push(CoedgeGeometry{
                pcurve:floor_uv[i],domain:range_for(u_start,u_end)?,
                face:FaceId(8),loop_id:LoopId(7)
            });
            m.geometry.coedges.push(CoedgeGeometry{
                pcurve:top_uv[i],domain:circle_domains[i],
                face:FaceId(1),loop_id:LoopId(6)
            });
        }
        // Coedge owner order already matches indices 24..27 exactly.
        for (id,curve,face,lp) in [
            (28,side_uv[0],6,8),(29,side_uv[1],6,8),
            (30,side_uv[4],7,9),(31,side_uv[5],7,9),
            (32,side_uv[2],6,8),(33,side_uv[7],7,9),
            (34,side_uv[3],6,8),(35,side_uv[6],7,9),
        ]{
            let _=id;
            m.geometry.coedges.push(CoedgeGeometry{
                pcurve:curve,domain:range_for(0.0,1.0)?,
                face:FaceId(face),loop_id:LoopId(lp),
            });
        }
        m.loops.extend([
            Loop{id:LoopId(6),face:FaceId(1),
                role:LoopRole::Inner,first_coedge:CoedgeId(25)},
            Loop{id:LoopId(7),face:FaceId(8),
                role:LoopRole::Outer,first_coedge:CoedgeId(24)},
            Loop{id:LoopId(8),face:FaceId(6),
                role:LoopRole::Outer,first_coedge:CoedgeId(28)},
            Loop{id:LoopId(9),face:FaceId(7),
                role:LoopRole::Outer,first_coedge:CoedgeId(30)},
        ]);
        m.faces[1].loops.push(LoopId(6));
        for i in 0..2{
            let side=FaceId(6+i as u32);
            let lp=LoopId(8+i as u32);
            m.faces.push(BrepFace{
                id:side,loops:vec![lp],sense:FaceSense::Reversed,
                geometry:FaceGeometry{
                    surface:side_surfaces[i],outer_loop:lp,tolerance:source.tolerance,
                    bounds:ParameterDomain2{
                        u:circle_domains[i],
                        v:range_for(0.0,region.cut_depth)?,
                        u_periodic:false,v_periodic:false,
                    },
                },
            });
        }
        let mut floor_geometry=m.faces[0].geometry.clone();
        floor_geometry.surface=floor_sid;
        floor_geometry.outer_loop=LoopId(7);
        m.faces.push(BrepFace{id:FaceId(8),loops:vec![LoopId(7)],
            sense:FaceSense::Reversed,geometry:floor_geometry});
        m.shells[0].faces.extend([FaceId(6),FaceId(7),FaceId(8)]);
        let volume_removed=region.disk_union_area*region.cut_depth;
        let volume=source.mass.volume-volume_removed;
        let surface_area=source.mass.surface_area+
            region.exposed_perimeter*region.cut_depth;
        if !volume.is_finite()||volume<=0.0||!surface_area.is_finite()||
           surface_area<=0.0 {
            return Err(invalid("crossed-circle committed mass is nonpositive"));
        }
        let cy=source.bbox.max.y-region.cut_depth*0.5;
        m.mass=MassProperties{
            volume,surface_area,
            centroid:Point3{
                x:source.mass.centroid.x-
                    region.relative_area_moment[0]*region.cut_depth/volume,
                y:source.mass.centroid.y-
                    volume_removed*(cy-source.mass.centroid.y)/volume,
                z:source.mass.centroid.z-
                    region.relative_area_moment[1]*region.cut_depth/volume,
            },
        };
        m.origin=BrepOrigin::AnalyticIntersectingBlindHoles;
        m.validate()?;
        let s=m.summary();
        if s.vertices!=12 || s.edges!=18 || s.coedges!=36 ||
            s.loops!=10 || s.faces!=9 || s.shells!=1 || s.genus!=0 ||
            s.euler_characteristic!=2 || s.circular_edges!=4 ||
            s.cylindrical_faces!=2 || m.faces[1].loops.len()!=2 ||
            m.loops.iter().filter(|l|l.role==LoopRole::Inner).count()!=1{
            return Err(invalid("overlapping-disk cut has incorrect analytic B-rep topology"));
        }
        Ok(m)
    }
}

pub struct IntersectingBlindTransaction<'a>{
    owner:&'a mut IntersectingBlindBrep,
    staged:[BlindHoleParameters;2],
    poisoned:bool,
}
impl IntersectingBlindTransaction<'_>{
    pub fn replace(&mut self,index:usize,parameters:BlindHoleParameters)
        ->Result<(),GeometryError>
    {
        if self.poisoned{return Err(invalid("crossed-circle transaction is poisoned"));}
        if index>=2 || !parameters.radius.is_finite()||
            !parameters.depth.is_finite()||!parameters.center_x.is_finite()||
            !parameters.center_z.is_finite()||parameters.radius<=0.0||
            parameters.depth<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid index or staged intersecting-circle parameters"));
        }
        self.staged[index]=parameters;
        Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned crossed-circle transaction"));}
        let region=CrossedDiskRegion::from_box(&self.owner.source,self.staged)?;
        let model=IntersectingBlindBrep::assemble(
            &self.owner.source,self.staged,&region)?;
        let next=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("crossed-circle revision overflow"))?;
        if self.owner.undo.len()>=128{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters);
        self.owner.redo.clear();
        self.owner.parameters=self.staged;
        self.owner.region=region;
        self.owner.model=model;
        self.owner.revision=next;
        Ok(next)
    }
}
