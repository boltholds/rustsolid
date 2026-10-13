//! A strictly scoped, executable analytic B-rep for ONE transversely
//! overlapping upright through-cylinder and top-entry blind cylinder.
//!
//! This creates real shared edges, coedges, UV p-curves and one closed genus-1
//! shell. It does not enable the general curved Boolean/TopologyDeltaPlan path.
//! The existing MixedThroughBlindBrep multi-hole executor remains fail-closed.
use crate::{
    BlindHoleParameters, BrepFace, BrepModel, BrepOrigin, BrepSource,
    Circle3, Coedge, CoedgeGeometry, CoedgeId, Curve2, Curve2Id, Curve3,
    Curve3Id, CylinderSurface, Edge, EdgeGeometry, EdgeId, FaceGeometry,
    FaceId, FaceSense, Frame3, GeometryError, HoleParameters, Line2, Line3,
    Loop, LoopId, LoopRole, MassProperties, Mesh, MixedCrossingSurfacePlan,
    MixedCrossingSurfacePlanner, ModelUnits, OrientedCircle2,
    ParameterDomain2, ParameterRange, PlaneSurface, Point2Param, Point3,
    Solid, Surface3, SurfaceId, Vertex, VertexId,
};
use std::f64::consts::{PI,TAU};

const MAX_HISTORY: usize=128;
const EMPTY_COEDGE:CoedgeId=CoedgeId(u32::MAX);

fn invalid(msg:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(msg.into())
}
pub(crate) fn bounded(a:f64,b:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(a,b)
}
fn next_id(len:usize)->Result<u32,GeometryError>{
    u32::try_from(len).map_err(|_|invalid("B-rep index overflow"))
}
fn p2(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn line2(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:p2(b.u-a.u,b.v-a.v)})
}
pub(crate) fn add_surface(model:&mut BrepModel,surface:Surface3)->Result<SurfaceId,GeometryError>{
    let id=SurfaceId(next_id(model.geometry.surfaces.len())?);
    model.geometry.surfaces.push(surface);Ok(id)
}
fn add_curve3(model:&mut BrepModel,curve:Curve3)->Result<Curve3Id,GeometryError>{
    let id=Curve3Id(next_id(model.geometry.curves3.len())?);
    model.geometry.curves3.push(curve);Ok(id)
}
fn add_curve2(model:&mut BrepModel,curve:Curve2)->Result<Curve2Id,GeometryError>{
    let id=Curve2Id(next_id(model.geometry.curves2.len())?);
    model.geometry.curves2.push(curve);Ok(id)
}
pub(crate) fn add_vertex(model:&mut BrepModel,point:Point3)->Result<VertexId,GeometryError>{
    let id=VertexId(next_id(model.vertices.len())?);
    model.vertices.push(Vertex{id,position:point});Ok(id)
}
fn add_edge(model:&mut BrepModel,start:VertexId,end:VertexId,
            carrier:Curve3,domain:ParameterRange)->Result<EdgeId,GeometryError>{
    let id=EdgeId(next_id(model.edges.len())?);
    let cid=add_curve3(model,carrier)?;
    model.edges.push(Edge{id,start,end,coedges:[EMPTY_COEDGE;2]});
    model.geometry.edges.push(EdgeGeometry{
        curve:cid,domain,tolerance:model.tolerance,
    });
    Ok(id)
}
pub(crate) fn add_circle(model:&mut BrepModel,frame:Frame3,radius:f64,y:f64,
              angles:(f64,f64),start:VertexId,end:VertexId)->Result<EdgeId,GeometryError>{
    let local=Frame3{origin:Point3{y,..frame.origin},..frame};
    let circle=Circle3::new(local,radius,model.tolerance)?;
    add_edge(model,start,end,Curve3::Circle(circle),bounded(angles.0,angles.1)?)
}
pub(crate) fn add_vertical(model:&mut BrepModel,start:VertexId,end:VertexId)->Result<EdgeId,GeometryError>{
    let a=model.vertices[start.0 as usize].position;
    let b=model.vertices[end.0 as usize].position;
    let (line,domain)=Line3::between(a,b,model.tolerance)?;
    add_edge(model,start,end,Curve3::Line(line),domain)
}

/// Unwrap the physical angular coordinate into the actual chart of its
/// owning face, preserving exact 2*pi seam identity where required.
fn chart_angle(raw:f64,face_domain:ParameterRange)->Result<f64,GeometryError>{
    if let ParameterRange::Bounded{start,end}=face_domain {
        let target=(start+end)/2.0;
        let candidate=raw+TAU*((target-raw)/TAU).round();
        if candidate<start-1e-7||candidate>end+1e-7 {
            return Err(invalid("circle angular coordinate is outside the source patch"));
        }
        Ok(candidate)
    }else{Err(invalid("cylinder patch needs finite angular bounds"))}
}
fn oriented_uv_circle(plane:PlaneSurface,circle:Circle3,
                      interval:ParameterRange,reversed:bool,
                      tol:crate::GeometryTolerance)->Result<(Curve2,ParameterRange),GeometryError>{
    let ParameterRange::Bounded{start,end}=interval else{
        return Err(invalid("circle needs finite angular parameters"));
    };
    let origin=plane.frame.project(circle.frame.origin)?;
    let axis=plane.frame.project(circle.evaluate(0.0)?)?;
    let phase=(axis.v-origin.v).atan2(axis.u-origin.u);
    let n=plane.frame.normal;
    let k=circle.frame.normal;
    let dot=n.x*k.x+n.y*k.y+n.z*k.z;
    if (dot.abs()-1.0).abs()>tol.angular {
        return Err(invalid("plane and circular carrier are not coplanar"));
    }
    let handed=if dot>=0.0{1.0}else{-1.0};
    let (phase,direction)=if reversed{
        (phase+handed*(start+end),-handed)
    }else{(phase,handed)};
    Ok((Curve2::OrientedCircle(
        OrientedCircle2::new(origin,circle.radius,phase,direction,tol)?
    ),interval))
}
fn coedge_trim(model:&mut BrepModel,face:FaceId,edge:EdgeId,
               reversed:bool)->Result<CoedgeGeometry,GeometryError>{
    let face_data=&model.faces[face.0 as usize];
    let support=model.geometry.surfaces[face_data.geometry.surface.0 as usize];
    let edge_data=model.geometry.edges[edge.0 as usize];
    let carrier=model.geometry.curves3[edge_data.curve.0 as usize];
    let src_edge=&model.edges[edge.0 as usize];
    let tol=model.tolerance;
    let (curve,domain)=match (support,carrier) {
        (Surface3::Plane(plane),Curve3::Circle(circle)) =>
            oriented_uv_circle(plane,circle,edge_data.domain,reversed,tol)?,
        (Surface3::Cylinder(cylinder),Curve3::Circle(circle))=>{
            let ParameterRange::Bounded{start,end}=edge_data.domain else{
                return Err(invalid("circular edge has no parameter range"));
            };
            let shift=TAU*(match face_data.geometry.bounds.u {
                ParameterRange::Bounded{start:u0,end:u1}=>((u0+u1)/2.0-(start+end)/2.0)/TAU,
                _=>return Err(invalid("cylinder face needs angular domain")),
            }).round();
            let (a,b)=if reversed{(end+shift,start+shift)}
                else{(start+shift,end+shift)};
            let v=(circle.frame.origin.y-cylinder.frame.origin.y)*
                cylinder.frame.normal.y;
            (line2(p2(a,v),p2(b,v)),bounded(0.0,1.0)?)
        }
        (Surface3::Cylinder(cylinder),Curve3::Line(_))=>{
            let start=model.vertices[src_edge.start.0 as usize].position;
            let end=model.vertices[src_edge.end.0 as usize].position;
            let (start,end)=if reversed{(end,start)}else{(start,end)};
            let resolve=|p:Point3|->Result<Point2Param,GeometryError>{
                let raw=(p.z-cylinder.frame.origin.z)
                    .atan2(p.x-cylinder.frame.origin.x);
                let u=chart_angle(raw,face_data.geometry.bounds.u)?;
                let v=(p.y-cylinder.frame.origin.y)*cylinder.frame.normal.y;
                Ok(p2(u,v))
            };
            (line2(resolve(start)?,resolve(end)?),bounded(0.0,1.0)?)
        }
        _=>return Err(invalid("unsupported new edge/face analytical trim binding")),
    };
    let pc=add_curve2(model,curve)?;
    Ok(CoedgeGeometry{pcurve:pc,domain,face,loop_id:LoopId(u32::MAX)})
}
pub(crate) fn append_loop(model:&mut BrepModel,face:FaceId,role:LoopRole,
               parts:&[(EdgeId,bool)])->Result<LoopId,GeometryError>{
    if parts.len()<2 {return Err(invalid("analytic loop needs at least two edges"));}
    let id=LoopId(next_id(model.loops.len())?);
    let first=next_id(model.coedges.len())?;
    let n=parts.len();
    for (i,&(edge,reversed)) in parts.iter().enumerate() {
        let cid=CoedgeId(first+i as u32);
        let bindings=&mut model.edges[edge.0 as usize].coedges;
        if bindings[0]==EMPTY_COEDGE {bindings[0]=cid;}
        else if bindings[1]==EMPTY_COEDGE {bindings[1]=cid;}
        else{return Err(invalid("analytical edge has more than two uses"));}
        let mut uv=coedge_trim(model,face,edge,reversed)?;
        uv.loop_id=id;
        model.geometry.coedges.push(uv);
        model.coedges.push(Coedge{
            id:cid,edge,face,loop_id:id,reversed,
            next:CoedgeId(first+(i+1).rem_euclid(n) as u32),
            prev:CoedgeId(first+(i+n-1).rem_euclid(n) as u32),
            twin:EMPTY_COEDGE,
        });
    }
    model.loops.push(Loop{id,face,role,first_coedge:CoedgeId(first)});
    let target=&mut model.faces[face.0 as usize];
    if target.loops.is_empty() {
        if role!=LoopRole::Outer {return Err(invalid("face starts with nonouter loop"));}
        target.geometry.outer_loop=id;
    }else if role!=LoopRole::Inner {
        return Err(invalid("only first face loop may be outer"));
    }
    target.loops.push(id);
    Ok(id)
}
pub(crate) fn twins(model:&mut BrepModel,from:usize)->Result<(),GeometryError>{
    for e in &model.edges[from..]{
        let [a,b]=e.coedges;
        if a==EMPTY_COEDGE||b==EMPTY_COEDGE||a==b {
            return Err(invalid("unpaired analytical topological edge"));
        }
        let ca=&model.coedges[a.0 as usize];
        let cb=&model.coedges[b.0 as usize];
        if ca.reversed==cb.reversed {return Err(invalid("edge has equal-direction coedge uses"));}
        model.coedges[a.0 as usize].twin=b;
        model.coedges[b.0 as usize].twin=a;
    }
    Ok(())
}

/// A dedicated source-owner for exactly one transverse crossing of a through
/// and a blind cylindrical cut. Keeps its own strong-atomic source history.
#[derive(Debug,Clone)]
pub struct MixedCrossingBrep {
    source:Solid,
    through:HoleParameters,
    blind:BlindHoleParameters,
    model:BrepModel,
    plan:MixedCrossingSurfacePlan,
    revision:u64,
    undo:Vec<(HoleParameters,BlindHoleParameters)>,
    redo:Vec<(HoleParameters,BlindHoleParameters)>,
}
impl MixedCrossingBrep {
    pub fn from_box(source:&Solid,t:HoleParameters,b:BlindHoleParameters)
        ->Result<Self,GeometryError>{
        let (model,plan)=Self::build(source,t,b)?;
        Ok(Self{source:source.clone(),through:t,blind:b,model,plan,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn plan(&self)->&MixedCrossingSurfacePlan{&self.plan}
    pub fn through(&self)->HoleParameters{self.through}
    pub fn blind(&self)->BlindHoleParameters{self.blind}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (expected,plan)=Self::build(&self.source,self.through,self.blind)?;
        // A matching entity count and identical carriers are not sufficient:
        // verify every physical endpoint, directed incidence and UV owner
        // against the independently rebuilt canonical source.
        let matches_source=
            self.model.origin==expected.origin &&
            self.model.summary()==expected.summary() &&
            self.model.bbox==expected.bbox &&
            self.model.tolerance==expected.tolerance &&
            self.model.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>()==
                expected.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>() &&
            self.model.edges.len()==expected.edges.len() &&
            self.model.edges.iter().zip(&expected.edges).all(|(a,b)|
                a.id==b.id && a.start==b.start && a.end==b.end &&
                a.coedges==b.coedges) &&
            self.model.coedges.len()==expected.coedges.len() &&
            self.model.coedges.iter().zip(&expected.coedges).all(|(a,b)|
                a.id==b.id && a.edge==b.edge && a.face==b.face &&
                a.loop_id==b.loop_id && a.reversed==b.reversed &&
                a.next==b.next && a.prev==b.prev && a.twin==b.twin) &&
            self.model.loops.len()==expected.loops.len() &&
            self.model.loops.iter().zip(&expected.loops).all(|(a,b)|
                a.id==b.id && a.face==b.face && a.role==b.role &&
                a.first_coedge==b.first_coedge) &&
            self.model.faces.len()==expected.faces.len() &&
            self.model.faces.iter().zip(&expected.faces).all(|(a,b)|
                a.id==b.id && a.loops==b.loops && a.sense==b.sense &&
                a.geometry==b.geometry) &&
            self.model.shells.len()==expected.shells.len() &&
            self.model.shells.iter().zip(&expected.shells).all(|(a,b)|
                a.id==b.id && a.faces==b.faces && a.closed==b.closed) &&
            self.model.geometry.curves3==expected.geometry.curves3 &&
            self.model.geometry.curves2==expected.geometry.curves2 &&
            self.model.geometry.surfaces==expected.geometry.surfaces &&
            self.model.geometry.edges==expected.geometry.edges &&
            self.model.geometry.coedges==expected.geometry.coedges &&
            self.model.mass==expected.mass &&
            self.plan==plan;
        if !matches_source {
            return Err(invalid("crossed through/blind topology or mass differs from source"));
        }
        Ok(())
    }
    pub fn begin_edit(&mut self)->MixedCrossingTransaction<'_>{
        let through=self.through;
        let blind=self.blind;
        MixedCrossingTransaction{owner:self,through,blind,poisoned:false}
    }
    pub(crate) fn source_box(&self)->&Solid{&self.source}
    pub fn tessellate(&self,segments:usize)->Result<crate::MixedCrossingDisplay,GeometryError>{
        crate::MixedCrossingMesher::tessellate(self,segments)
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let desired=if redo{&self.redo}else{&self.undo};
        let Some(&(t,b))=desired.last()else{return Ok(false)};
        let (model,plan)=Self::build(&self.source,t,b)?;
        let revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("mixed-crossing revision overflow"))?;
        let old=(self.through,self.blind);
        if redo{self.redo.pop();self.undo.push(old);}
        else{self.undo.pop();self.redo.push(old);}
        self.through=t;self.blind=b;self.model=model;self.plan=plan;
        self.revision=revision;
        Ok(true)
    }
    fn build(source:&Solid,t:HoleParameters,b:BlindHoleParameters)
        ->Result<(BrepModel,MixedCrossingSurfacePlan),GeometryError>{
        let plan=MixedCrossingSurfacePlanner::plan(source,t,b)?;
        plan.validate_against(source,t,b)?;
        let arc=plan.footprint.exposed_arcs
            .ok_or_else(||invalid("crossed source has no exposed arcs"))?;
        let points=plan.footprint.crossing_points
            .ok_or_else(||invalid("crossed source has no point witnesses"))?;
        let through_arc=arc[0];let blind_arc=arc[1];
        if through_arc.end_crossing!=blind_arc.start_crossing ||
            blind_arc.end_crossing!=through_arc.start_crossing {
            return Err(invalid("source arcs do not form one aperture"));
        }
        let mut model=BrepModel::from_polyhedral(source,ModelUnits::default())?;
        if model.vertices.len()!=8 || model.edges.len()!=12 ||
            model.coedges.len()!=24 || model.faces.len()!=6 ||
            model.loops.len()!=6 || model.shells.len()!=1 {
            return Err(invalid("crossing requires canonical six-face block"));
        }
        let tol=source.tolerance;
        let y0=source.bbox.min.y;
        let yf=plan.splitting_height_y;
        let yt=source.bbox.max.y;
        let x_axis=Point3{x:1.0,y:0.0,z:0.0};
        let z_axis=Point3{x:0.0,y:0.0,z:1.0};
        // The circle frames use world X/Z angle CCW (normal -Y); both
        // cylinders are parameterized DOWNWARD from the box's upper plane.
        let tf=Frame3::from_axes(Point3{x:t.center_x,y:yt,z:t.center_z},
            x_axis,z_axis,tol)?;
        let bf=Frame3::from_axes(Point3{x:b.center_x,y:yt,z:b.center_z},
            x_axis,z_axis,tol)?;
        let ff=Frame3::from_axes(Point3{x:0.0,y:yf,z:0.0},
            x_axis,z_axis,tol)?;
        let ts=add_surface(&mut model,Surface3::Cylinder(
            CylinderSurface::new(tf,t.radius,tol)?))?;
        let bs=add_surface(&mut model,Surface3::Cylinder(
            CylinderSurface::new(bf,b.radius,tol)?))?;
        let fs=add_surface(&mut model,Surface3::Plane(PlaneSurface{frame:ff}))?;
        let t0=through_arc.start_angle;
        let t1=t0+through_arc.sweep_radians;
        let t2=t0+TAU;
        let b0=blind_arc.start_angle;
        let b1=b0+blind_arc.sweep_radians;
        let bounds=|u0:f64,u1:f64,v0:f64,v1:f64|->Result<ParameterDomain2,GeometryError>{
            Ok(ParameterDomain2{u:bounded(u0,u1)?,v:bounded(v0,v1)?,
                u_periodic:false,v_periodic:false})
        };
        let h=yt-y0;
        let d=b.depth;
        let supports=[
            (ts,bounds(t0,t1,d,h)?),
            (ts,bounds(t1,t2,d,h)?),
            (ts,bounds(t0,t1,0.0,d)?),
            (bs,bounds(b0,b1,0.0,d)?),
            (fs,ParameterDomain2{u:ParameterRange::Unbounded,
                v:ParameterRange::Unbounded,u_periodic:false,v_periodic:false})
        ];
        for (sid,uv) in supports {
            let fid=FaceId(next_id(model.faces.len())?);
            model.faces.push(BrepFace {
                id:fid,loops:Vec::new(),sense:FaceSense::Reversed,
                geometry:FaceGeometry{
                    surface:sid,bounds:uv,outer_loop:LoopId(u32::MAX),
                    tolerance:tol
                }
            });
            model.shells[0].faces.push(fid);
        }
        // Logical crossing 0/1 are defined by through_arc traversal. A
        // through-paired bottom, floor, top set shares identical X/Z data.
        let crossings=[points[through_arc.start_crossing],points[through_arc.end_crossing]];
        let mut v=[VertexId(0);6];
        for (i,y) in [y0,yf,yt].into_iter().enumerate() {
            for (j,c) in crossings.into_iter().enumerate() {
                v[2*i+j]=add_vertex(&mut model,Point3{x:c.x,y,z:c.z})?;
            }
        }
        // Circle physical edges are all nonclosed intervals and retain the
        // exact analytical through/blind carrier; no tessellation.
        let mut e=[EdgeId(0);11];
        e[0]=add_circle(&mut model,tf,t.radius,y0,(t0,t1),v[0],v[1])?;
        e[1]=add_circle(&mut model,tf,t.radius,y0,(t1,t2),v[1],v[0])?;
        e[2]=add_circle(&mut model,tf,t.radius,yf,(t0,t1),v[2],v[3])?;
        e[3]=add_circle(&mut model,tf,t.radius,yf,(t1,t2),v[3],v[2])?;
        e[4]=add_circle(&mut model,bf,b.radius,yf,(b0,b1),v[3],v[2])?;
        e[5]=add_circle(&mut model,tf,t.radius,yt,(t0,t1),v[4],v[5])?;
        e[6]=add_circle(&mut model,bf,b.radius,yt,(b0,b1),v[5],v[4])?;
        e[7]=add_vertical(&mut model,v[0],v[2])?;
        e[8]=add_vertical(&mut model,v[1],v[3])?;
        e[9]=add_vertical(&mut model,v[2],v[4])?;
        e[10]=add_vertical(&mut model,v[3],v[5])?;
        // Bottom and top each have their ORIGINAL box outer loop plus one
        // genuine inner loop. Floor is the blind-minus-through crescent.
        append_loop(&mut model,FaceId(0),LoopRole::Inner,
            &[(e[0],true),(e[1],true)])?;
        append_loop(&mut model,FaceId(1),LoopRole::Inner,
            &[(e[5],false),(e[6],false)])?;
        append_loop(&mut model,FaceId(10),LoopRole::Outer,
            &[(e[4],true),(e[3],false)])?;
        // Every cylindrical patch is a real four-sided UV trim. Splitting
        // the bottom through-cylinder into TWO source-equivalent patches
        // avoids any zero-length or nonphysical periodic seam edge.
        append_loop(&mut model,FaceId(6),LoopRole::Outer,
            &[(e[2],true),(e[7],true),(e[0],false),(e[8],false)])?;
        append_loop(&mut model,FaceId(7),LoopRole::Outer,
            &[(e[3],true),(e[8],true),(e[1],false),(e[7],false)])?;
        append_loop(&mut model,FaceId(8),LoopRole::Outer,
            &[(e[5],true),(e[9],true),(e[2],false),(e[10],false)])?;
        append_loop(&mut model,FaceId(9),LoopRole::Outer,
            &[(e[6],true),(e[10],true),(e[4],false),(e[9],false)])?;
        twins(&mut model,12)?;
        // Exact-in-form mass is independent of the boundary curves and any
        // eventual facet mesh. The overlap first moment uses the two
        // circular segments' opposite first moments, which cancel.
        let center_distance=(b.center_x-t.center_x).hypot(b.center_z-t.center_z);
        let x1=(center_distance*center_distance+t.radius*t.radius-b.radius*b.radius)
            /(2.0*center_distance);
        let x2=center_distance-x1;
        let chord=(t.radius*t.radius-x1*x1).sqrt();
        let area_1=t.radius*t.radius*chord.atan2(x1)-x1*chord;
        let area_2=b.radius*b.radius*chord.atan2(x2)-x2*chord;
        let overlap=plan.footprint.overlap_area
            .ok_or_else(||invalid("missing proven disk overlap"))?;
        if (area_1+area_2-overlap).abs()>tol.area_at(t.radius+b.radius)*16.0 {
            return Err(invalid("disk lens moments disagree with source area"));
        }
        let at=PI*t.radius*t.radius;
        let ab=PI*b.radius*b.radius;
        let v_through=at*h;
        let v_blind=(ab-overlap)*d;
        let c0=source.mass.centroid;
        let first_moment_x=(t.center_x-c0.x)*v_through +
            ((b.center_x-c0.x)*ab -
             ((t.center_x-c0.x)*area_1+(b.center_x-c0.x)*area_2))*d;
        let first_moment_z=(t.center_z-c0.z)*v_through +
            ((b.center_z-c0.z)*ab -
             ((t.center_z-c0.z)*area_1+(b.center_z-c0.z)*area_2))*d;
        let yc=yt-d/2.0;
        let volume=plan.predicted_remaining_volume;
        if !volume.is_finite()||volume<=0.0{return Err(invalid("invalid crossed volume"));}
        model.mass=MassProperties {
            volume,surface_area:plan.predicted_surface_area,
            centroid:Point3 {
                x:c0.x-first_moment_x/volume,
                y:c0.y-v_blind*(yc-c0.y)/volume,
                z:c0.z-first_moment_z/volume,
            }
        };
        model.origin=BrepOrigin::AnalyticCrossingThroughBlind;
        model.validate()?;
        let summary=model.summary();
        if summary.vertices!=14 || summary.edges!=23 ||
           summary.coedges!=46 || summary.faces!=11 ||
           summary.loops!=13 || summary.shells!=1 ||
           summary.euler_characteristic!=0 || summary.genus!=1 ||
           summary.circular_edges!=7 || summary.cylindrical_faces!=4 ||
           model.faces[0].loops.len()!=2 ||
           model.faces[1].loops.len()!=2 ||
           model.loops.iter().filter(|l|l.role==LoopRole::Inner).count()!=2 {
            return Err(invalid("crossing analytic shell has incorrect topology"));
        }
        Ok((model,plan))
    }
}

impl BrepSource for MixedCrossingBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct MixedCrossingTransaction<'a>{
    owner:&'a mut MixedCrossingBrep,
    through:HoleParameters,
    blind:BlindHoleParameters,
    poisoned:bool,
}
impl MixedCrossingTransaction<'_>{
    pub fn replace_through(&mut self,value:HoleParameters)->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("crossing transaction is poisoned"));}
        if !value.center_x.is_finite()||!value.center_z.is_finite()||
            !value.radius.is_finite()||value.radius<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid staged through cut"));
        }
        self.through=value;Ok(())
    }
    pub fn replace_blind(&mut self,value:BlindHoleParameters)->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("crossing transaction is poisoned"));}
        if !value.center_x.is_finite()||!value.center_z.is_finite()||
            !value.radius.is_finite()||value.radius<=0.0 ||
            !value.depth.is_finite()||value.depth<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid staged blind cut"));
        }
        self.blind=value;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned crossing transaction"));}
        let (model,plan)=MixedCrossingBrep::build(
            &self.owner.source,self.through,self.blind)?;
        let next=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("crossing revision overflow"))?;
        if self.owner.undo.len()>=MAX_HISTORY{self.owner.undo.remove(0);}
        self.owner.undo.push((self.owner.through,self.owner.blind));
        self.owner.redo.clear();
        self.owner.through=self.through;
        self.owner.blind=self.blind;
        self.owner.model=model;
        self.owner.plan=plan;
        self.owner.revision=next;
        Ok(next)
    }
}
