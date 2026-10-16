//! G4b2b3b: exact, source-certified, four-component curved Boolean shell.
//!
//! Bounded strict family: an AABB minus a full-height upright circular
//! cylinder whose disk cuts ALL four walls and contains NO box corner.
//! This results in four independent, closed, genus-zero corner shells.
//! Source G1 vertices and both source UV cutting curves are independently
//! revalidated before introducing genuine new physical cut edges.
//! This is an auditable pre-G4c/G5 construction; it is never committed.

use std::collections::{BTreeMap,BTreeSet};
use crate::{
    AnalyticSourceEdgeSplitHold, BrepBody, BrepModel, BrepGeometry,
    BrepFace, BrepOrigin, BooleanError, BooleanOperation, Coedge,
    CoedgeGeometry, CoedgeId, Curve2, Curve2Id, Curve3, Curve3Id,
    Edge, EdgeGeometry, EdgeId, FaceGeometry, FaceId, FaceSense,
    Frame3, GeometryError, GeometryTolerance, GlobalAnalyticCutPlan,
    GlobalAnalyticEdgeSplitPlan, GlobalAnalyticEdgeSplitResult,
    GlobalArcRef, GlobalEventId, GlobalEventKind, GraphOperand,
    HarmonicCurve2, Line2, Line3, Loop, LoopId, LoopRole,
    MassProperties, ParameterDomain2, ParameterRange, Point2Param,
    Point3, Shell, ShellId, Surface3, SurfaceId, Vertex, VertexId,
    Circle3, BoundingBox,
};
use crate::radial_shell_audit::{
    certify_radial_shell,RadialHalfedge,RadialLoop,
};
use crate::{RadialTopologyCertificate,RadialTopologyHold};

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn bound(lo:f64,hi:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(lo,hi)
}
fn parts(range:ParameterRange)->Result<(f64,f64),GeometryError>{
    match range{
        ParameterRange::Bounded{start,end}=>Ok((start,end)),
        _=>Err(invalid("G4b2b3b requires bounded physical cut carriers")),
    }
}
fn uv(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2 {
        origin:a,delta:Point2Param {u:b.u-a.u,v:b.v-a.v},
    })
}
fn finite(v:f64)->bool{v.is_finite()}
fn unit_normal(v:Point3)->Result<Point3,GeometryError>{
    let norm=v.length_squared().sqrt();
    if !norm.is_finite()||norm<=0.0 {return Err(invalid("invalid face normal"));}
    Ok(v.scale(1.0/norm))
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticFacePatchHold {
    UpstreamSourceSplit(AnalyticSourceEdgeSplitHold),
    UnsupportedOperation,
    UnsupportedOperand,
    NonCanonicalBox,
    NonUprightCylinder,
    NoFourSeparateCorners,
    AmbiguousCutGeometry,
    MissingOriginalMaterialEdge,
    MissingSourceCutArc,
    DuplicateOrExtraCutArc,
    IncorrectGlobalEvent,
    InvalidRadialTopology(RadialTopologyHold),
    InvalidComponentTopology,
    DegenerateAnalyticMass,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct AnalyticPatchFaceSource {
    pub result_face: FaceId,
    pub shell: ShellId,
    pub source_operand: GraphOperand,
    pub source_face: FaceId,
    pub flipped_for_difference: bool,
}
#[derive(Debug,Clone,PartialEq)]
pub struct AnalyticCornerPatch {
    /// -1 is the lower world-axis box wall; +1 is the upper.
    pub x_side: i8,
    pub z_side: i8,
    pub shell: ShellId,
    /// Ordered [bottom cap, top cap, x wall, z wall, cylinder].
    pub faces: [FaceId;5],
    /// Ordered [bottom C-A, bottom C-B, bottom arc,
    /// top C-A, top C-B, top arc, corner vertical, A vertical, B vertical].
    pub edges: [EdgeId;9],
    /// Certified original G4b2b2 arcs for [bottom, top, A, B] new edges.
    pub cut_arcs: [GlobalArcRef;4],
    /// G1 event identities for BOTH endpoints of each real cut edge,
    /// ordered along that edge's canonical 3D carrier.
    pub global_endpoints: [[GlobalEventId;2];4],
    pub exact_area_xy: f64,
    pub exact_volume: f64,
}
#[derive(Debug,Clone)]
pub struct AnalyticFacePatchPartitionPlan {
    model: BrepModel,
    pub corners: Vec<AnalyticCornerPatch>,
    pub source_faces: Vec<AnalyticPatchFaceSource>,
    pub radial: RadialTopologyCertificate,
    pub edge_splits: GlobalAnalyticEdgeSplitPlan,
    /// Strict source-validated curved shell closure, not a general Boolean.
    pub closed_shell_proven: bool,
    /// Remains false until independent G4c and G5.
    pub can_apply: bool,
}
#[derive(Debug,Clone)]
pub enum AnalyticFacePatchPartitionResult {
    Ready(AnalyticFacePatchPartitionPlan),
    Held(AnalyticFacePatchHold),
}
impl AnalyticFacePatchPartitionPlan {
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        self.edge_splits.validate_against(first,second,tolerance)?;
        let AnalyticFacePatchPartitionResult::Ready(rebuilt)=
            Self::build(first,second,BooleanOperation::Difference,tolerance)? else{
                return Err(invalid("G4b2b3b changed source family or graph"));
            };
        if self.corners!=rebuilt.corners
            ||self.source_faces!=rebuilt.source_faces
            ||self.radial!=rebuilt.radial
            ||self.closed_shell_proven!=rebuilt.closed_shell_proven
            ||self.can_apply!=rebuilt.can_apply
            ||format!("{:?}",self.model)!=format!("{:?}",rebuilt.model)
        {
            return Err(invalid("G4b2b3b physical curved shell differs from source rebuild"));
        }
        self.model.validate()
    }

    pub fn build(
        first:&BrepBody,second:&BrepBody,
        operation:BooleanOperation,tolerance:GeometryTolerance,
    )->Result<AnalyticFacePatchPartitionResult,GeometryError>{
        tolerance.validate()?;
        if operation!=BooleanOperation::Difference {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::UnsupportedOperation));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cylinder))=(first,second)
            else{return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::UnsupportedOperand));};
        match crate::boolean::certify_box(block) {
            Ok(())=>{},
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::NonCanonicalBox)),
        }
        let split=match GlobalAnalyticEdgeSplitPlan::build(first,second,tolerance)?{
            GlobalAnalyticEdgeSplitResult::Ready(plan)=>plan,
            GlobalAnalyticEdgeSplitResult::Held(hold)=>
                return Ok(AnalyticFacePatchPartitionResult::Held(
                    AnalyticFacePatchHold::UpstreamSourceSplit(hold))),
        };
        // This bounded producer uses a fixed world-axis extrusion. A rotated
        // cylinder chart or tilted axis needs a separate UV recharting path.
        let frame=cylinder.frame;
        let orientation=distance(frame.x,p(1.,0.,0.))
            .max(distance(frame.y,p(0.,0.,-1.)))
            .max(distance(frame.normal,p(0.,1.,0.)));
        if orientation>tolerance.angular*16.0{
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::NonUprightCylinder));
        }
        let box_bounds=block.bbox;
        let extent=box_bounds.max.sub(box_bounds.min)
            .length_squared().sqrt().max(cylinder.radius).max(cylinder.height).max(1.);
        let eps=tolerance.length_at(extent)*64.0;
        let ymin=box_bounds.min.y;
        let ymax=box_bounds.max.y;
        if frame.origin.y>=ymin-eps
            ||frame.origin.y+cylinder.height<=ymax+eps
            ||!(frame.origin.x>box_bounds.min.x+eps
                &&frame.origin.x<box_bounds.max.x-eps
                &&frame.origin.z>box_bounds.min.z+eps
                &&frame.origin.z<box_bounds.max.z-eps)
        {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::NoFourSeparateCorners));
        }
        let source_box=first.shared()?;
        let source_cyl=second.shared()?;
        let source_cyl_side=source_cyl.faces.iter()
            .filter(|f|matches!(source_cyl.face_surface(f.id),Some(Surface3::Cylinder(_))))
            .map(|f|f.id).collect::<Vec<_>>();
        if source_cyl_side.len()!=1{
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::UnsupportedOperand));
        }
        let wall=source_cyl_side[0];
        let original_cut=&split.registry;
        let mut builder=PatchBuilder::new(box_bounds,tolerance,source_box.units);
        let mut corners=Vec::new();
        let mut sources=Vec::new();
        let mut used_cut_arcs=BTreeSet::<GlobalArcRef>::new();
        let mut original_edge_uses=BTreeSet::<EdgeId>::new();
        let mut event_vertices=BTreeMap::<GlobalEventId,VertexId>::new();
        let mut mass_area=0.0;
        let mut area3=0.0;
        let mut mass_x=0.0;
        let mut mass_y=0.0;
        let mut mass_z=0.0;
        for sx in [-1_i8,1_i8] {
            for sz in [-1_i8,1_i8] {
                let xlim=if sx<0{box_bounds.min.x}else{box_bounds.max.x};
                let zlim=if sz<0{box_bounds.min.z}else{box_bounds.max.z};
                let a=(xlim-frame.origin.x)*f64::from(sx);
                let b=(zlim-frame.origin.z)*f64::from(sz);
                let r=cylinder.radius;
                if a<=eps||b<=eps||r<=a+eps||r<=b+eps
                    ||a.hypot(b)<=r+eps {
                    return Ok(AnalyticFacePatchPartitionResult::Held(
                        AnalyticFacePatchHold::NoFourSeparateCorners));
                }
                let v0=(r*r-a*a).sqrt();
                let u0=(r*r-b*b).sqrt();
                if !v0.is_finite()||!u0.is_finite()
                    ||v0<=eps||u0<=eps||v0>=b-eps||u0>=a-eps
                    ||a-u0<=eps||b-v0<=eps {
                    return Ok(AnalyticFacePatchPartitionResult::Held(
                        AnalyticFacePatchHold::NoFourSeparateCorners));
                }
                let height=ymax-ymin;
                // Exact source-coordinate area/first moments of the corner
                // outside a circle; this is bounded-family metadata, NOT G4c.
                let f=|u:f64|->f64{
                    let root=(r*r-u*u).max(0.).sqrt();
                    0.5*(u*root+r*r*(u/r).clamp(-1.,1.).asin())
                };
                let area=b*(a-u0)-(f(a)-f(u0));
                let mu=0.5*b*(a*a-u0*u0)-(b*b*b-v0*v0*v0)/3.0;
                let mv=0.5*((b*b-r*r)*(a-u0)+(a*a*a-u0*u0*u0)/3.0);
                if ![area,mu,mv].iter().all(|x|finite(*x))
                    ||area<=eps*eps*16.0 ||mu<=0.0||mv<=0.0 {
                    return Ok(AnalyticFacePatchPartitionResult::Held(
                        AnalyticFacePatchHold::DegenerateAnalyticMass));
                }
                let xc=frame.origin.x;
                let zc=frame.origin.z;
                let c=p(xlim,ymin,zlim);
                let a0=p(xlim,ymin,zc+f64::from(sz)*v0);
                let b0=p(xc+f64::from(sx)*u0,ymin,zlim);
                let ct=p(xlim,ymax,zlim);
                let at=p(a0.x,ymax,a0.z);
                let bt=p(b0.x,ymax,b0.z);
                let ta=(-(a0.z-zc)).atan2(a0.x-xc);
                let tb=(-(b0.z-zc)).atan2(b0.x-xc);
                if (ta-tb).abs()>=std::f64::consts::PI
                    ||(ta-tb).abs()<=tolerance.angular*32.0 {
                    return Ok(AnalyticFacePatchPartitionResult::Held(
                        AnalyticFacePatchHold::AmbiguousCutGeometry));
                }
                let tmin=ta.min(tb);
                let tmax=ta.max(tb);
                let mid=0.5*(ta+tb);
                let radial=p(mid.cos(),0.,-mid.sin());
                let side_area=r*(tmax-tmin)*height;
                let surface_area=2.0*area+height*((a-u0)+(b-v0))+side_area;
                let bottom_face=find_source_plane(
                    &source_box,1,ymin,-1.,eps)
                    .ok_or_else(||invalid("G4b2b3b original bottom face missing"))?;
                let top_face=find_source_plane(
                    &source_box,1,ymax,1.,eps)
                    .ok_or_else(||invalid("G4b2b3b original top face missing"))?;
                let x_face=find_source_plane(
                    &source_box,0,xlim,f64::from(sx),eps)
                    .ok_or_else(||invalid("G4b2b3b original x-wall face missing"))?;
                let z_face=find_source_plane(
                    &source_box,2,zlim,f64::from(sz),eps)
                    .ok_or_else(||invalid("G4b2b3b original z-wall face missing"))?;
                let specs=[
                    FaceSpec{local:vec![0,1,2],source:bottom_face,operand:GraphOperand::First,
                        normal:p(0.,-1.,0.)},
                    FaceSpec{local:vec![3,4,5],source:top_face,operand:GraphOperand::First,
                        normal:p(0.,1.,0.)},
                    FaceSpec{local:vec![0,1,4,3],source:x_face,operand:GraphOperand::First,
                        normal:p(f64::from(sx),0.,0.)},
                    FaceSpec{local:vec![0,2,5,3],source:z_face,operand:GraphOperand::First,
                        normal:p(0.,0.,f64::from(sz))},
                    FaceSpec{local:vec![1,2,5,4],source:wall,operand:GraphOperand::Second,
                        normal:radial.scale(-1.)},
                ];
                let coords=[c,a0,b0,ct,at,bt];
                let shell_id=ShellId(u32::try_from(builder.model.shells.len())
                    .map_err(|_|invalid("too many corner shells"))?);
                let out=builder.add_corner(&coords,ta,tb,frame,r,&specs,
                    &source_box,&source_cyl)?;
                // Five original material edges MUST correspond to real,
                // already source-subdivided edges, not fabricated XYZ joins.
                let expected_faces=[
                    [bottom_face,x_face],[bottom_face,z_face],
                    [top_face,x_face],[top_face,z_face],
                    [x_face,z_face],
                ];
                for (index,faces) in [0usize,1,3,4,6].into_iter().zip(expected_faces) {
                    let new_edge=out.edges[index];
                    let edge=&builder.model.edges[new_edge.0 as usize];
                    if !certify_original_material_edge(
                        &split.split_sources[0],edge,
                        &builder.model.vertices,faces,eps,
                        &mut original_edge_uses)
                    {
                        return Ok(AnalyticFacePatchPartitionResult::Held(
                            AnalyticFacePatchHold::MissingOriginalMaterialEdge));
                    }
                }
                let cut_specs=[
                    (out.edges[2],bottom_face),
                    (out.edges[5],top_face),
                    (out.edges[7],x_face),
                    (out.edges[8],z_face),
                ];
                let mut arc_ids=Vec::new();
                let mut events=Vec::new();
                for (edge_id,face_id) in cut_specs {
                    let edge=&builder.model.edges[edge_id.0 as usize];
                    let curve=builder.model.edge_curve(edge_id)
                        .ok_or_else(||invalid("G4b2b3b missing cut Curve3"))?;
                    let domain=builder.model.geometry.edges[edge_id.0 as usize].domain;
                    let endpoints=[
                        builder.model.vertices[edge.start.0 as usize].position,
                        builder.model.vertices[edge.end.0 as usize].position,
                    ];
                    let certified=certify_new_cut(
                        original_cut,face_id,wall,curve,domain,endpoints,
                        [edge.start,edge.end],eps,&mut used_cut_arcs,
                        &mut event_vertices);
                    let Some((arc,ids))=certified? else{
                        return Ok(AnalyticFacePatchPartitionResult::Held(
                            AnalyticFacePatchHold::MissingSourceCutArc));
                    };
                    arc_ids.push(arc);
                    events.push(ids);
                }
                for (face_index,&fid) in out.faces.iter().enumerate(){
                    sources.push(AnalyticPatchFaceSource {
                        result_face:fid,shell:shell_id,
                        source_operand:specs[face_index].operand,
                        source_face:specs[face_index].source,
                        flipped_for_difference:face_index==4,
                    });
                }
                let component_volume=area*height;
                mass_area+=component_volume;
                area3+=surface_area;
                mass_x+=component_volume*(xc+f64::from(sx)*mu/area);
                mass_z+=component_volume*(zc+f64::from(sz)*mv/area);
                mass_y+=component_volume*(0.5*(ymin+ymax));
                corners.push(AnalyticCornerPatch {
                    x_side:sx,z_side:sz,shell:shell_id,faces:out.faces,
                    edges:out.edges,cut_arcs:arc_ids.try_into()
                        .map_err(|_|invalid("missing four cut arcs"))?,
                    global_endpoints:events.try_into()
                        .map_err(|_|invalid("missing four G1 endpoint pairs"))?,
                    exact_area_xy:area,exact_volume:component_volume,
                });
            }
        }
        // The G1/G2 registry was built from *all* face pairs. A strict
        // four-corner partition must consume every transverse source arc
        // exactly once; ignoring an extra arc would lose material topology.
        if used_cut_arcs.len()!=original_cut.arcs.len()
            ||event_vertices.len()!=original_cut.physical_event_count {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::DuplicateOrExtraCutArc));
        }
        if corners.len()!=4 || !mass_area.is_finite() || mass_area<=0.0 {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::InvalidComponentTopology));
        }
        builder.model.mass=MassProperties {
            volume:mass_area,surface_area:area3,
            centroid:p(mass_x/mass_area,mass_y/mass_area,mass_z/mass_area),
        };
        let built=builder.finish()?;
        if built.shells.len()!=4 || built.faces.len()!=20
            || built.vertices.len()!=24 || built.edges.len()!=36
            || built.coedges.len()!=72 || built.euler_characteristic()!=8 {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::InvalidComponentTopology));
        }
        built.validate()?;
        let radial=radial_certificate(&built,eps)
            .map_err(|e|invalid(format!("invalid G4b2b3b radial proof: {e:?}")))?;
        if radial.shell_components!=4 ||radial.total_genus!=0
            ||radial.euler_characteristic!=8
            ||radial.vertices!=24 ||radial.edges!=36
            ||radial.material_face_cells!=20 {
            return Ok(AnalyticFacePatchPartitionResult::Held(
                AnalyticFacePatchHold::InvalidComponentTopology));
        }
        Ok(AnalyticFacePatchPartitionResult::Ready(Self {
            model:built,corners,source_faces:sources,radial,
            edge_splits:split,closed_shell_proven:true,can_apply:false,
        }))
    }
}

fn find_source_plane(
    model:&BrepModel,axis:u8,coordinate:f64,normal_sign:f64,eps:f64,
)->Option<FaceId>{
    let mut matches=Vec::new();
    for face in &model.faces{
        let Some(Surface3::Plane(plane))=model.face_surface(face.id) else{continue};
        let normal=model.oriented_face_normal(face.id,0.,0.).ok()?;
        let (origin,sign)=match axis {
            0=>(plane.frame.origin.x,normal.x),
            1=>(plane.frame.origin.y,normal.y),
            _=>(plane.frame.origin.z,normal.z),
        };
        if (origin-coordinate).abs()<=eps&&(sign-normal_sign).abs()<1e-8 {
            matches.push(face.id);
        }
    }
    if matches.len()==1{Some(matches[0])}else{None}
}
fn certify_original_material_edge(
    source:&BrepModel,result:&Edge,vertices:&[Vertex],
    faces:[FaceId;2],eps:f64,used:&mut BTreeSet<EdgeId>,
)->bool{
    let a=vertices[result.start.0 as usize].position;
    let b=vertices[result.end.0 as usize].position;
    let mut matching=Vec::new();
    for edge in &source.edges {
        let other=[
            source.coedges[edge.coedges[0].0 as usize].face,
            source.coedges[edge.coedges[1].0 as usize].face,
        ];
        if !((other==faces)||(other==[faces[1],faces[0]])){continue;}
        let c=source.vertices[edge.start.0 as usize].position;
        let d=source.vertices[edge.end.0 as usize].position;
        if distance(a,c)<=eps&&distance(b,d)<=eps
            ||distance(a,d)<=eps&&distance(b,c)<=eps {
            matching.push(edge.id);
        }
    }
    if matching.len()!=1 ||used.contains(&matching[0]){return false;}
    used.insert(matching[0]);
    true
}
fn certify_new_cut(
    registry:&GlobalAnalyticCutPlan,first_face:FaceId,second_face:FaceId,
    curve:Curve3,domain:ParameterRange,
    endpoints:[Point3;2],vertices:[VertexId;2],eps:f64,
    used:&mut BTreeSet<GlobalArcRef>,
    event_vertices:&mut BTreeMap<GlobalEventId,VertexId>,
)->Result<Option<(GlobalArcRef,[GlobalEventId;2])>,GeometryError>{
    let (lo,hi)=parts(domain)?;
    let mut found=Vec::new();
    for arc in &registry.arcs{
        if arc.first_face!=first_face||arc.second_face!=second_face
            ||used.contains(&arc.source) {continue;}
        let Some(ids)=arc.physical_endpoints else{continue;};
        let (start,end)=parts(arc.parameter_domain)?;
        let a=arc.carrier.evaluate(start)?;
        let b=arc.carrier.evaluate(end)?;
        let reversed=if distance(endpoints[0],a)<=eps
            &&distance(endpoints[1],b)<=eps {false}
            else if distance(endpoints[0],b)<=eps
                &&distance(endpoints[1],a)<=eps {true}
            else {continue;};
        let mut curve_ok=true;
        for fraction in [0.,0.125,0.25,0.5,0.75,0.875,1.]{
            let t=lo+(hi-lo)*fraction;
            let ref_frac=if reversed{1.-fraction}else{fraction};
            let ref_t=start+(end-start)*ref_frac;
            if distance(curve.evaluate(t)?,arc.carrier.evaluate(ref_t)?)>eps{
                curve_ok=false;break;
            }
        }
        if curve_ok {found.push((arc.source,if reversed {[ids[1],ids[0]]}else{ids}));}
    }
    if found.len()!=1{return Ok(None);}
    let (source,events)=found[0];
    for (id,vertex,point) in [
        (events[0],vertices[0],endpoints[0]),
        (events[1],vertices[1],endpoints[1]),
    ]{
        let event=registry.global_nodes.events.iter()
            .find(|e|e.id==id && e.kind==GlobalEventKind::BoundaryEndpoint)
            .ok_or_else(||invalid("G4b2b3b cut edge missing true G1 event"))?;
        if distance(point,event.point)>eps {
            return Ok(None);
        }
        if let Some(previous)=event_vertices.get(&id){
            if *previous!=vertex{return Ok(None);}
        }else{event_vertices.insert(id,vertex);}
    }
    used.insert(source);
    Ok(Some((source,events)))
}

struct FaceSpec {
    local:Vec<usize>,
    source:FaceId,
    operand:GraphOperand,
    normal:Point3,
}
struct AddedCorner{faces:[FaceId;5],edges:[EdgeId;9]}
struct PatchBuilder{model:BrepModel,tolerance:GeometryTolerance}
impl PatchBuilder{
    fn new(bbox:BoundingBox,tolerance:GeometryTolerance,units:crate::ModelUnits)->Self{
        Self {model:BrepModel{
            origin:BrepOrigin::AnalyticPartialBoxCylinderDifference,
            vertices:Vec::new(),edges:Vec::new(),coedges:Vec::new(),
            loops:Vec::new(),faces:Vec::new(),shells:Vec::new(),
            geometry:BrepGeometry {
                curves3:Vec::new(),curves2:Vec::new(),
                surfaces:Vec::new(),edges:Vec::new(),coedges:Vec::new(),
            },
            mass:MassProperties {volume:1.,surface_area:1.,centroid:p(0.,0.,0.)},
            bbox,tolerance,units,
        },tolerance}
    }
    fn add_edge(
        &mut self,start:VertexId,end:VertexId,
        curve:Option<(Circle3,ParameterRange)>,
    )->Result<EdgeId,GeometryError>{
        let a=self.model.vertices[start.0 as usize].position;
        let b=self.model.vertices[end.0 as usize].position;
        let (carrier,domain)=match curve{
            Some((circle,domain))=>{
                let (lo,hi)=parts(domain)?;
                if distance(circle.evaluate(lo)?,a)>self.tolerance.absolute_length*32.
                    ||distance(circle.evaluate(hi)?,b)>self.tolerance.absolute_length*32. {
                    return Err(invalid("G4b2b3b circle endpoints mismatch"));
                }
                (Curve3::Circle(circle),domain)
            }
            None=>{
                let (line,dom)=Line3::between(a,b,self.tolerance)?;
                (Curve3::Line(line),dom)
            }
        };
        let id=EdgeId(self.model.edges.len() as u32);
        let curve_id=Curve3Id(self.model.geometry.curves3.len() as u32);
        self.model.geometry.curves3.push(carrier);
        self.model.geometry.edges.push(EdgeGeometry {
            curve:curve_id,domain,tolerance:self.tolerance,
        });
        self.model.edges.push(Edge {
            id,start,end,coedges:[CoedgeId(u32::MAX);2],
        });
        Ok(id)
    }
    fn add_corner(
        &mut self,coords:&[Point3;6],ta:f64,tb:f64,
        frame:Frame3,r:f64,specs:&[FaceSpec;5],
        first:&BrepModel,second:&BrepModel,
    )->Result<AddedCorner,GeometryError>{
        let vbase=self.model.vertices.len();
        for &point in coords {
            let id=VertexId(self.model.vertices.len() as u32);
            self.model.vertices.push(Vertex {id,position:point});
        }
        let v=|index:usize|VertexId((vbase+index) as u32);
        let circle_endpoints=if ta<tb { [1usize,2usize] } else { [2usize,1usize] };
        let mk_circle=|y:f64|{
            Circle3 {frame:Frame3{origin:p(frame.origin.x,y,frame.origin.z),..frame},
                radius:r}
        };
        let circle_domain=bound(ta.min(tb),ta.max(tb))?;
        let edge_ids=[
            self.add_edge(v(0),v(1),None)?,
            self.add_edge(v(0),v(2),None)?,
            self.add_edge(v(circle_endpoints[0]),v(circle_endpoints[1]),
                Some((mk_circle(coords[0].y),circle_domain)))?,
            self.add_edge(v(3),v(4),None)?,
            self.add_edge(v(3),v(5),None)?,
            self.add_edge(v(circle_endpoints[0]+3),v(circle_endpoints[1]+3),
                Some((mk_circle(coords[3].y),circle_domain)))?,
            self.add_edge(v(0),v(3),None)?,
            self.add_edge(v(1),v(4),None)?,
            self.add_edge(v(2),v(5),None)?,
        ];
        let edge_pairs=[
            (0,1),(0,2),(circle_endpoints[0],circle_endpoints[1]),
            (3,4),(3,5),(circle_endpoints[0]+3,circle_endpoints[1]+3),
            (0,3),(1,4),(2,5),
        ];
        let shell=ShellId(self.model.shells.len() as u32);
        let mut face_ids=Vec::new();
        for (face_index,spec) in specs.iter().enumerate() {
            let source=if spec.operand==GraphOperand::First{first}else{second};
            let surface=source.face_surface(spec.source)
                .ok_or_else(||invalid("G4b2b3b missing source supporting surface"))?;
            if (face_index==4)!=matches!(surface,Surface3::Cylinder(_)){
                return Err(invalid("G4b2b3b original supporting surface mislabeled"));
            }
            let fid=FaceId(self.model.faces.len() as u32);
            let lid=LoopId(self.model.loops.len() as u32);
            let sid=SurfaceId(self.model.geometry.surfaces.len() as u32);
            self.model.geometry.surfaces.push(surface);
            let mut local=spec.local.clone();
            let q0=coords[local[0]];
            let q1=coords[local[1]];
            let q2=coords[local[2]];
            let candidate_normal=q1.sub(q0).cross(q2.sub(q1));
            let sign=candidate_normal.dot(spec.normal);
            if !finite(sign)||sign.abs()<=self.tolerance.absolute_length.powi(2) {
                return Err(invalid("G4b2b3b degenerate oriented patch"));
            }
            if sign<0.0 {local.reverse();}
            let mut bounds_points=Vec::new();
            for &index in &local{
                let point=coords[index];
                bounds_points.push(match surface{
                    Surface3::Plane(plane)=>plane.frame.project(point)?,
                    Surface3::Cylinder(_)=>Point2Param {
                        u:if index==1||index==4{ta}else{tb},
                        v:point.y-frame.origin.y,
                    },
                    _=>return Err(invalid("G4b2b3b unsupported surface")),
                });
            }
            let umin=bounds_points.iter().map(|p|p.u).fold(f64::INFINITY,f64::min);
            let umax=bounds_points.iter().map(|p|p.u).fold(f64::NEG_INFINITY,f64::max);
            let vmin=bounds_points.iter().map(|p|p.v).fold(f64::INFINITY,f64::min);
            let vmax=bounds_points.iter().map(|p|p.v).fold(f64::NEG_INFINITY,f64::max);
            let epsuv=self.tolerance.absolute_length*2.;
            let face_bounds=ParameterDomain2{
                u:bound(umin-epsuv,umax+epsuv)?,
                v:bound(vmin-epsuv,vmax+epsuv)?,
                u_periodic:false,v_periodic:false,
            };
            self.model.faces.push(BrepFace {
                id:fid,loops:vec![lid],
                sense:if face_index==4{FaceSense::Reversed}else{FaceSense::Forward},
                geometry:FaceGeometry {
                    surface:sid,bounds:face_bounds,outer_loop:lid,
                    tolerance:self.tolerance,
                },
            });
            let first_coedge=CoedgeId(self.model.coedges.len() as u32);
            self.model.loops.push(Loop {
                id:lid,face:fid,role:LoopRole::Outer,first_coedge,
            });
            let n=local.len();
            for index in 0..n{
                let start=local[index];
                let end=local[(index+1)%n];
                let (edge_index,reversed)=edge_pairs.iter().enumerate()
                    .find_map(|(e,&(a,b))|{
                        if a==start && b==end {Some((e,false))}
                        else if a==end && b==start {Some((e,true))}
                        else {None}
                    }).ok_or_else(||invalid("G4b2b3b patch boundary lacks physical edge"))?;
                let eid=edge_ids[edge_index];
                let id=CoedgeId(self.model.coedges.len() as u32);
                let first=first_coedge.0;
                self.model.coedges.push(Coedge{
                    id,edge:eid,face:fid,loop_id:lid,reversed,
                    next:CoedgeId(first+((index+1)%n) as u32),
                    prev:CoedgeId(first+((index+n-1)%n) as u32),
                    twin:CoedgeId(u32::MAX),
                });
                let from=coords[start];
                let to=coords[end];
                let (pcurve,domain)=match surface {
                    Surface3::Plane(plane)=>{
                        if edge_index==2||edge_index==5 {
                            // A rotated and possibly reversed true analytic
                            // planar pcurve of the source cylinder circle.
                            let center=plane.frame.project(
                                p(frame.origin.x,from.y,frame.origin.z))?;
                            let cosine=Point2Param{
                                u:r*frame.x.dot(plane.frame.x),
                                v:r*frame.x.dot(plane.frame.y),
                            };
                            let sine=Point2Param{
                                u:r*frame.y.dot(plane.frame.x),
                                v:r*frame.y.dot(plane.frame.y),
                            };
                            let phase=if start==1||start==4{ta}else{tb};
                            let direction=if (start==1||start==4) == (ta<tb) {1.}else{-1.};
                            let (s,c)=phase.sin_cos();
                            let rotated_cos=Point2Param{
                                u:cosine.u*c+sine.u*s,
                                v:cosine.v*c+sine.v*s,
                            };
                            let rotated_sin=Point2Param{
                                u:direction*(-cosine.u*s+sine.u*c),
                                v:direction*(-cosine.v*s+sine.v*c),
                            };
                            (Curve2::Harmonic(HarmonicCurve2{
                                center,cosine:rotated_cos,sine:rotated_sin,
                            }),bound(0.,(ta-tb).abs())?)
                        }else{
                            (uv(plane.frame.project(from)?,plane.frame.project(to)?),
                                bound(0.,1.)?)
                        }
                    },
                    Surface3::Cylinder(_)=>{
                        let a=Point2Param {
                            u:if start==1||start==4{ta}else{tb},
                            v:from.y-frame.origin.y,
                        };
                        let b=Point2Param {
                            u:if end==1||end==4{ta}else{tb},
                            v:to.y-frame.origin.y,
                        };
                        (uv(a,b),bound(0.,1.)?)
                    },
                    _=>return Err(invalid("unsupported patch pcurve surface")),
                };
                let pcurve_id=Curve2Id(self.model.geometry.curves2.len() as u32);
                self.model.geometry.curves2.push(pcurve);
                self.model.geometry.coedges.push(CoedgeGeometry {
                    pcurve:pcurve_id,domain,face:fid,loop_id:lid,
                });
            }
            face_ids.push(fid);
        }
        for &edge in &edge_ids {
            let uses=self.model.coedges.iter()
                .filter(|co|co.edge==edge).map(|co|co.id).collect::<Vec<_>>();
            if uses.len()!=2{
                return Err(invalid("G4b2b3b physical edge has nonmanifold incidence"));
            }
            let a=uses[0];
            let b=uses[1];
            if self.model.coedges[a.0 as usize].reversed
                ==self.model.coedges[b.0 as usize].reversed {
                return Err(invalid("G4b2b3b physical cut edge has same-direction uses"));
            }
            self.model.coedges[a.0 as usize].twin=b;
            self.model.coedges[b.0 as usize].twin=a;
            self.model.edges[edge.0 as usize].coedges=[a,b];
        }
        self.model.shells.push(Shell {
            id:shell,faces:face_ids.clone(),closed:true,
        });
        Ok(AddedCorner {
            faces:face_ids.try_into().map_err(|_|invalid("missing corner faces"))?,
            edges:edge_ids,
        })
    }
    fn finish(self)->Result<BrepModel,GeometryError>{
        self.model.validate()?;
        Ok(self.model)
    }
}
fn radial_certificate(
    body:&BrepModel,eps:f64,
)->Result<RadialTopologyCertificate,RadialTopologyHold>{
    let uses=body.coedges.iter().map(|co|{
        let edge=&body.edges[co.edge.0 as usize];
        let (a,b)=if co.reversed{(edge.end,edge.start)}
            else{(edge.start,edge.end)};
        RadialHalfedge {
            start:body.vertices[a.0 as usize].position,
            end:body.vertices[b.0 as usize].position,
            face:co.face.0 as usize,
        }
    }).collect::<Vec<_>>();
    let mut loops=Vec::new();
    for lp in &body.loops{
        let mut ids=Vec::new();
        let mut at=lp.first_coedge;
        loop {
            ids.push(at.0 as usize);
            at=body.coedges[at.0 as usize].next;
            if at==lp.first_coedge{break;}
        }
        loops.push(RadialLoop {
            face:lp.face.0 as usize,uses:ids,
            inner:lp.role==LoopRole::Inner,
        });
    }
    let pairs=body.edges.iter().map(|edge|
        [edge.coedges[0].0 as usize,edge.coedges[1].0 as usize]
    ).collect::<Vec<_>>();
    certify_radial_shell(&uses,&loops,&pairs,body.faces.len(),eps)
}
