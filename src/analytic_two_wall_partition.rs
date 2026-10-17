//! G6a: source-audited analytic two-wall Box minus Cylinder difference.
//! Retained planar UV regions have circular trims, physical arc coedges,
//! and two independent closed genus-zero shells, not display-only facets.

use std::collections::{BTreeMap,BTreeSet};
use crate::{
    BrepBody,BrepModel,BrepOrigin,BrepFace,BooleanOperation,BooleanError,
    Coedge,CoedgeGeometry,CoedgeId,Curve2,Curve2Id,Curve3,Curve3Id,
    Edge,EdgeGeometry,EdgeId,FaceGeometry,FaceId,FaceSense,Frame3,
    GeometryError,GeometryTolerance,GlobalAnalyticEdgeSplitPlan,
    GlobalAnalyticEdgeSplitResult,AnalyticSourceEdgeSplitHold,
    GlobalAnalyticCutPlan,GlobalArcRef,GlobalEventId,GraphOperand,
    HarmonicCurve2,Line2,Line3,Loop,LoopId,LoopRole,
    MassProperties,ParameterDomain2,ParameterRange,Point2Param,
    Point3,Shell,ShellId,Surface3,SurfaceId,Vertex,VertexId,Circle3,
    RadialTopologyCertificate,
};
use crate::analytic_face_patch_partition::{
    PatchBuilder,FaceSpec,find_source_plane,certify_new_cut,
    certify_original_material_edge,radial_certificate,
};
fn invalid(message:impl Into<String>)->GeometryError {
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
    match range {
        ParameterRange::Bounded{start,end}=>Ok((start,end)),
        _=>Err(invalid("G6a requires bounded curve domains")),
    }
}
fn uv(a:Point2Param,b:Point2Param)->Curve2 {
    Curve2::Line(Line2{origin:a,delta:Point2Param {u:b.u-a.u,v:b.v-a.v}})
}
fn finite(v:f64)->bool{v.is_finite()}

#[derive(Debug,Clone)]
struct AddedTwoWall{faces:[FaceId;6],edges:[EdgeId;12]}
impl PatchBuilder {
    fn add_two_wall(
        &mut self,coords:&[Point3;8],ta:f64,tb:f64,
        frame:Frame3,r:f64,specs:&[FaceSpec;6],
        first:&BrepModel,second:&BrepModel,
    )->Result<AddedTwoWall,GeometryError>{
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
            self.add_edge(v(circle_endpoints[0]),v(circle_endpoints[1]),
                Some((mk_circle(coords[0].y),circle_domain)))?,
            self.add_edge(v(2),v(3),None)?,
            self.add_edge(v(3),v(0),None)?,
            self.add_edge(v(4),v(5),None)?,
            self.add_edge(v(circle_endpoints[0]+4),v(circle_endpoints[1]+4),
                Some((mk_circle(coords[4].y),circle_domain)))?,
            self.add_edge(v(6),v(7),None)?,
            self.add_edge(v(7),v(4),None)?,
            self.add_edge(v(0),v(4),None)?,
            self.add_edge(v(1),v(5),None)?,
            self.add_edge(v(2),v(6),None)?,
            self.add_edge(v(3),v(7),None)?,
        ];
        let edge_pairs=[
            (0,1),(circle_endpoints[0],circle_endpoints[1]),(2,3),(3,0),
            (4,5),(circle_endpoints[0]+4,circle_endpoints[1]+4),(6,7),(7,4),
            (0,4),(1,5),(2,6),(3,7),
        ];
        let shell=ShellId(self.model.shells.len() as u32);
        let mut face_ids=Vec::new();
        for (face_index,spec) in specs.iter().enumerate() {
            let source=if spec.operand==GraphOperand::First{first}else{second};
            let surface=source.face_surface(spec.source)
                .ok_or_else(||invalid("G4b2b3b missing source supporting surface"))?;
            if (face_index==5)!=matches!(surface,Surface3::Cylinder(_)){
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
                        u:if index==1||index==5{ta}else{tb},
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
                sense:if face_index==5{FaceSense::Reversed}else{FaceSense::Forward},
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
                        if edge_index==1||edge_index==5 {
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
                            let phase=if start==1||start==5{ta}else{tb};
                            let direction=if (start==1||start==5) == (ta<tb) {1.}else{-1.};
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
                            u:if start==1||start==5{ta}else{tb},
                            v:from.y-frame.origin.y,
                        };
                        let b=Point2Param {
                            u:if end==1||end==5{ta}else{tb},
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
        Ok(AddedTwoWall {
            faces:face_ids.try_into().map_err(|_|invalid("missing corner faces"))?,
            edges:edge_ids,
        })
    }

}


#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticTwoWallHold {
    UnsupportedOperation,
    UnsupportedOperand,
    NonCanonicalBox,
    UpstreamSplit(AnalyticSourceEdgeSplitHold),
    NotTwoWallTransverse,
    AmbiguousAnalyticArc,
    MissingSourceWall,
    MissingOriginalMaterialEdge,
    MissingCertifiedCut,
    UnconsumedGlobalArc,
    MaterialSideNotCertified,
    InconsistentShellTopology,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct TwoWallFaceSource {
    pub result_face:FaceId,
    pub shell:ShellId,
    pub operand:GraphOperand,
    pub source_face:FaceId,
    pub reversed_for_difference:bool,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TwoWallMaterialSide {
    pub face:FaceId,
    pub sample:Point3,
    pub retained:Point3,
    pub removed:Point3,
    pub outward:Point3,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TwoWallComponent {
    pub z_side:i8,
    pub shell:ShellId,
    /// bottom,top,left-x,right-x,outer-z,cylinder
    pub faces:[FaceId;6],
    /// bottom [left,arc,right,z]; top same; 4 vertical generators
    pub edges:[EdgeId;12],
    /// [bottom circle, top circle, left vertical cut, right vertical cut]
    pub global_arcs:[GlobalArcRef;4],
    pub global_endpoints:[[GlobalEventId;2];4],
    pub analytic_area:f64,
    pub analytic_volume:f64,
}
#[derive(Debug,Clone)]
pub struct AnalyticTwoWallPartitionPlan {
    model:BrepModel,
    pub components:Vec<TwoWallComponent>,
    pub source_faces:Vec<TwoWallFaceSource>,
    pub material_sides:Vec<TwoWallMaterialSide>,
    pub split_sources:GlobalAnalyticEdgeSplitPlan,
    pub radial:RadialTopologyCertificate,
    pub closed_shell_proven:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone)]
pub enum AnalyticTwoWallPartitionResult {
    Ready(AnalyticTwoWallPartitionPlan),
    Held(AnalyticTwoWallHold),
}
fn signed_material(
    point:Point3,bounds:crate::BoundingBox,cx:f64,cz:f64,
    r:f64,eps:f64,
)->Option<bool>{
    if !point.is_finite(){return None;}
    let margins=[
        point.x-bounds.min.x,bounds.max.x-point.x,
        point.y-bounds.min.y,bounds.max.y-point.y,
        point.z-bounds.min.z,bounds.max.z-point.z,
    ];
    let radial=(point.x-cx).hypot(point.z-cz);
    if margins.iter().any(|x|*x < -eps)||radial<r-eps{
        return Some(false);
    }
    if margins.iter().all(|x|*x>eps)&&radial>r+eps{
        return Some(true);
    }
    None
}
impl AnalyticTwoWallPartitionPlan {
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        self.split_sources.validate_against(first,second,tolerance)?;
        let AnalyticTwoWallPartitionResult::Ready(fresh)=Self::build(
            first,second,BooleanOperation::Difference,tolerance)? else{
            return Err(invalid("G6a original sources no longer yield two connected patches"));
        };
        if self.components!=fresh.components
            ||self.source_faces!=fresh.source_faces
            ||self.material_sides!=fresh.material_sides
            ||self.radial!=fresh.radial
            ||self.closed_shell_proven!=fresh.closed_shell_proven
            ||self.can_apply!=fresh.can_apply
            ||format!("{:?}",self.model)!=format!("{:?}",fresh.model)
        {
            return Err(invalid("G6a result differs from original source rebuilding"));
        }
        self.model.validate()
    }

    pub fn build(
        first:&BrepBody,second:&BrepBody,
        operation:BooleanOperation,tolerance:GeometryTolerance,
    )->Result<AnalyticTwoWallPartitionResult,GeometryError>{
        tolerance.validate()?;
        first.validate()?;second.validate()?;
        if operation!=BooleanOperation::Difference {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::UnsupportedOperation));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cylinder))=
            (first,second) else{
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::UnsupportedOperand));
        };
        match crate::boolean::certify_box(block) {
            Ok(())=>{},
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::NonCanonicalBox)),
        }
        let frame=cylinder.frame;
        let scale=block.bbox.max.sub(block.bbox.min).length_squared().sqrt()
            .max(cylinder.radius).max(cylinder.height).max(1.);
        let eps=tolerance.length_at(scale)*64.;
        if distance(frame.x,p(1.,0.,0.))>tolerance.angular*16.
            ||distance(frame.y,p(0.,0.,-1.))>tolerance.angular*16.
            ||distance(frame.normal,p(0.,1.,0.))>tolerance.angular*16.
        {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::NotTwoWallTransverse));
        }
        let bounds=block.bbox;
        let (xmin,xmax,zmin,zmax)=(
            bounds.min.x,bounds.max.x,bounds.min.z,bounds.max.z
        );
        let cx=frame.origin.x;let cz=frame.origin.z;
        let r=cylinder.radius;let h=bounds.max.y-bounds.min.y;
        let a_left=cx-xmin;let a_right=xmax-cx;
        let b_neg=cz-zmin;let b_pos=zmax-cz;
        if ![a_left,a_right,b_neg,b_pos,h,r].iter().all(|x|finite(*x))
            ||a_left<=eps||a_right<=eps||b_neg<=eps||b_pos<=eps
            ||r<=a_left+eps||r<=a_right+eps
            ||r>=b_neg-eps||r>=b_pos-eps
            ||frame.origin.y>=bounds.min.y-eps
            ||frame.origin.y+cylinder.height<=bounds.max.y+eps
        {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::NotTwoWallTransverse));
        }
        let split=match GlobalAnalyticEdgeSplitPlan::build(first,second,tolerance)?{
            GlobalAnalyticEdgeSplitResult::Ready(x)=>x,
            GlobalAnalyticEdgeSplitResult::Held(why)=>
                return Ok(AnalyticTwoWallPartitionResult::Held(
                    AnalyticTwoWallHold::UpstreamSplit(why))),
        };
        let original=&split.registry;
        let source_a=first.shared()?;
        let source_b=second.shared()?;
        let cylinder_faces=source_b.faces.iter().filter(|f|
            matches!(source_b.face_surface(f.id),Some(Surface3::Cylinder(_)))
        ).map(|f|f.id).collect::<Vec<_>>();
        if cylinder_faces.len()!=1{
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::UnsupportedOperand));
        }
        let cutter_face=cylinder_faces[0];
        let bottom=find_source_plane(&source_a,1,bounds.min.y,-1.,eps);
        let top=find_source_plane(&source_a,1,bounds.max.y,1.,eps);
        let xleft=find_source_plane(&source_a,0,xmin,-1.,eps);
        let xright=find_source_plane(&source_a,0,xmax,1.,eps);
        let zneg=find_source_plane(&source_a,2,zmin,-1.,eps);
        let zpos=find_source_plane(&source_a,2,zmax,1.,eps);
        let (Some(bottom),Some(top),Some(xleft),Some(xright),
            Some(zneg),Some(zpos))=(bottom,top,xleft,xright,zneg,zpos)
        else{
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::MissingSourceWall));
        };
        let mut builder=PatchBuilder::new(bounds,tolerance,source_a.units);
        builder.model.origin=BrepOrigin::AnalyticTwoWallCylinderDifference;
        let mut components=Vec::new();
        let mut sources=Vec::new();
        let mut material_sides=Vec::new();
        let mut used_edges=BTreeSet::<EdgeId>::new();
        let mut used_arcs=BTreeSet::<GlobalArcRef>::new();
        let mut events=BTreeMap::<GlobalEventId,VertexId>::new();
        let mut volume=0.;
        let mut area_surface=0.;
        let mut moment_x=0.;
        let mut moment_z=0.;
        for sz in [-1_i8,1_i8]{
            let limit=if sz<0{zmin}else{zmax};
            let b=(limit-cz).abs();
            let zl=((r-a_left)*(r+a_left)).sqrt();
            let zr=((r-a_right)*(r+a_right)).sqrt();
            let delta_x=xmax-xmin;
            // Analytic source-area metadata only; an independent source
            // quadrature and flux is required by G6b before G5 commits.
            let f=|u:f64| {
                let root=((r-u)*(r+u)).max(0.0).sqrt();
                0.5*(u*root+r*r*(u/r).clamp(-1.,1.).asin())
            };
            let ua=xmin-cx;let ub=xmax-cx;
            let area=b*delta_x-(f(ub)-f(ua));
            let m_u=0.5*b*(ub*ub-ua*ua)
                +(((r*r-ub*ub).max(0.)).powf(1.5)
                    -((r*r-ua*ua).max(0.)).powf(1.5))/3.;
            let m_v=0.5*((b*b-r*r)*delta_x
                +(ub*ub*ub-ua*ua*ua)/3.);
            if !area.is_finite()||area<=eps*eps*16.
                ||!m_u.is_finite()||!m_v.is_finite()||m_v<=0. {
                return Ok(AnalyticTwoWallPartitionResult::Held(
                    AnalyticTwoWallHold::NotTwoWallTransverse));
            }
            let pl=p(xmin,bounds.min.y,limit);
            let al=p(xmin,bounds.min.y,cz+f64::from(sz)*zl);
            let ar=p(xmax,bounds.min.y,cz+f64::from(sz)*zr);
            let pr=p(xmax,bounds.min.y,limit);
            let pt=p(pl.x,bounds.max.y,pl.z);
            let atl=p(al.x,bounds.max.y,al.z);
            let atr=p(ar.x,bounds.max.y,ar.z);
            let ptr=p(pr.x,bounds.max.y,pr.z);
            let ta=(-(al.z-cz)).atan2(al.x-cx);
            let tb=(-(ar.z-cz)).atan2(ar.x-cx);
            let span=(ta-tb).abs();
            if span<=tolerance.angular*32.
                ||span>=std::f64::consts::PI-tolerance.angular*32. {
                return Ok(AnalyticTwoWallPartitionResult::Held(
                    AnalyticTwoWallHold::AmbiguousAnalyticArc));
            }
            let middle=0.5*(ta+tb);
            let normal=p(-middle.cos(),0.,middle.sin());
            let zwall=if sz<0{zneg}else{zpos};
            let specs=[
                FaceSpec{local:vec![0,1,2,3],source:bottom,
                    operand:GraphOperand::First,normal:p(0.,-1.,0.)},
                FaceSpec{local:vec![4,5,6,7],source:top,
                    operand:GraphOperand::First,normal:p(0.,1.,0.)},
                FaceSpec{local:vec![0,1,5,4],source:xleft,
                    operand:GraphOperand::First,normal:p(-1.,0.,0.)},
                FaceSpec{local:vec![2,3,7,6],source:xright,
                    operand:GraphOperand::First,normal:p(1.,0.,0.)},
                FaceSpec{local:vec![0,3,7,4],source:zwall,
                    operand:GraphOperand::First,normal:p(0.,0.,f64::from(sz))},
                FaceSpec{local:vec![1,2,6,5],source:cutter_face,
                    operand:GraphOperand::Second,normal},
            ];
            let coords=[pl,al,ar,pr,pt,atl,atr,ptr];
            let shell=ShellId(builder.model.shells.len() as u32);
            let result=builder.add_two_wall(
                &coords,ta,tb,frame,r,&specs,&source_a,&source_b)?;
            let original_spec=[
                (0usize,[bottom,xleft]),(2,[bottom,xright]),(3,[bottom,zwall]),
                (4,[top,xleft]),(6,[top,xright]),(7,[top,zwall]),
                (8,[xleft,zwall]),(11,[xright,zwall]),
            ];
            for (index,faces) in original_spec{
                let id=result.edges[index];
                let e=&builder.model.edges[id.0 as usize];
                if !certify_original_material_edge(
                    &split.split_sources[0],e,
                    &builder.model.vertices,faces,eps,&mut used_edges)
                {
                    return Ok(AnalyticTwoWallPartitionResult::Held(
                        AnalyticTwoWallHold::MissingOriginalMaterialEdge));
                }
            }
            let cuts=[
                (1usize,bottom),(5,top),(9,xleft),(10,xright),
            ];
            let mut source_arcs=Vec::new();
            let mut source_events=Vec::new();
            for (index,first_face) in cuts{
                let eid=result.edges[index];
                let e=&builder.model.edges[eid.0 as usize];
                let physical=builder.model.edge_curve(eid)
                    .ok_or_else(||invalid("G6a missing analytic cut carrier"))?;
                let domain=builder.model.geometry.edges[eid.0 as usize].domain;
                let xy=[
                    builder.model.vertices[e.start.0 as usize].position,
                    builder.model.vertices[e.end.0 as usize].position,
                ];
                let matched=certify_new_cut(
                    original,first_face,cutter_face,physical,domain,
                    xy,[e.start,e.end],eps,&mut used_arcs,&mut events)?;
                let Some((arc,ids))=matched else{
                    return Ok(AnalyticTwoWallPartitionResult::Held(
                        AnalyticTwoWallHold::MissingCertifiedCut));
                };
                source_arcs.push(arc);
                source_events.push(ids);
            }
            for (i,&fid) in result.faces.iter().enumerate(){
                sources.push(TwoWallFaceSource{
                    result_face:fid,shell,
                    operand:specs[i].operand,source_face:specs[i].source,
                    reversed_for_difference:i==5,
                });
            }
            let step=eps*32.;
            if b-r<=step*8. {
                return Ok(AnalyticTwoWallPartitionResult::Held(
                    AnalyticTwoWallHold::MaterialSideNotCertified));
            }
            let zmid=cz+f64::from(sz)*(r+b)/2.;
            let xc_mid=(xmin+xmax)/2.;
            let ym=(bounds.min.y+bounds.max.y)/2.;
            let samples=[
                p(xc_mid,bounds.min.y,zmid),p(xc_mid,bounds.max.y,zmid),
                p(xmin,ym,cz+f64::from(sz)*(zl+b)/2.),
                p(xmax,ym,cz+f64::from(sz)*(zr+b)/2.),
                p(xc_mid,ym,limit),
                p(cx+r*middle.cos(),ym,cz-r*middle.sin()),
            ];
            let normals=[
                p(0.,-1.,0.),p(0.,1.,0.),p(-1.,0.,0.),
                p(1.,0.,0.),p(0.,0.,f64::from(sz)),normal,
            ];
            for (i,sample) in samples.into_iter().enumerate(){
                let n=normals[i];
                let retained=sample.sub(n.scale(step));
                let removed=sample.add(n.scale(step));
                if signed_material(retained,bounds,cx,cz,r,eps)!=Some(true)
                    ||signed_material(removed,bounds,cx,cz,r,eps)!=Some(false){
                    return Ok(AnalyticTwoWallPartitionResult::Held(
                        AnalyticTwoWallHold::MaterialSideNotCertified));
                }
                material_sides.push(TwoWallMaterialSide{
                    face:result.faces[i],sample,retained,removed,outward:n,
                });
            }
            let v=area*h;
            volume+=v;
            moment_x+=v*(cx+m_u/area);
            moment_z+=v*(cz+f64::from(sz)*m_v/area);
            area_surface+=2.*area+h*((b-zl)+(b-zr)+delta_x+r*span);
            components.push(TwoWallComponent{
                z_side:sz,shell,faces:result.faces,edges:result.edges,
                global_arcs:source_arcs.try_into().map_err(|_|
                    invalid("G6a must account for four true cut edges"))?,
                global_endpoints:source_events.try_into().map_err(|_|
                    invalid("G6a must account for eight endpoints"))?,
                analytic_area:area,analytic_volume:v,
            });
        }
        if components.len()!=2||used_arcs.len()!=original.arcs.len()
            ||events.len()!=original.physical_event_count {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::UnconsumedGlobalArc));
        }
        builder.model.mass=MassProperties{
            volume,surface_area:area_surface,
            centroid:p(moment_x/volume,(bounds.min.y+bounds.max.y)/2.,
                moment_z/volume),
        };
        let model=builder.finish()?;
        if model.shells.len()!=2 ||model.faces.len()!=12
            ||model.loops.len()!=12||model.edges.len()!=24
            ||model.vertices.len()!=16||model.coedges.len()!=48
            ||model.euler_characteristic()!=4 {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::InconsistentShellTopology));
        }
        let radial=radial_certificate(&model,eps)
            .map_err(|e|invalid(format!("G6a radial certification: {e:?}")))?;
        if radial.shell_components!=2||radial.total_genus!=0
            ||radial.euler_characteristic!=4||radial.vertices!=16
            ||radial.edges!=24||radial.material_face_cells!=12 {
            return Ok(AnalyticTwoWallPartitionResult::Held(
                AnalyticTwoWallHold::InconsistentShellTopology));
        }
        Ok(AnalyticTwoWallPartitionResult::Ready(Self{
            model,components,source_faces:sources,material_sides,
            split_sources:split,radial,closed_shell_proven:true,can_apply:false,
        }))
    }
}
