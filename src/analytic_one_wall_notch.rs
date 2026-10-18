//! G7a: bounded analytic one-wall cylinder notch partition.
//!
//! The exact source circle crosses ONLY one of the box's X walls.
//! Its retained major arc is a true Circle3/Curve2 carrier (> pi),
//! not a polygonal approximation; one shell contains two distinct
//! retained planar patches of the very same original X-wall face.

use std::collections::{BTreeMap,BTreeSet};
use crate::{
    BrepBody,BrepModel,BrepOrigin,BrepFace,BooleanOperation,BooleanError,
    Coedge,CoedgeGeometry,CoedgeId,Curve2,Curve2Id,Curve3,Curve3Id,
    Edge,EdgeGeometry,EdgeId,FaceGeometry,FaceId,FaceSense,Frame3,
    GeometryError,GeometryTolerance,GlobalAnalyticEdgeSplitPlan,
    GlobalAnalyticEdgeSplitResult,AnalyticSourceEdgeSplitHold,
    GlobalAnalyticCutPlan,GlobalArcRef,GlobalEventId,GraphOperand,
    HarmonicCurve2,Line2,Line3,Loop,LoopId,LoopRole,MassProperties,
    ParameterDomain2,ParameterRange,Point2Param,Point3,Shell,ShellId,
    Surface3,SurfaceId,Vertex,VertexId,Circle3,RadialTopologyCertificate,
};
use crate::analytic_face_patch_partition::{
    PatchBuilder,FaceSpec,find_source_plane,certify_new_cut,
    certify_original_material_edge,radial_certificate,
};
fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn bound(lo:f64,hi:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(lo,hi)
}
fn finite(v:f64)->bool{v.is_finite()}
fn uv(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:Point2Param{u:b.u-a.u,v:b.v-a.v}})
}
#[derive(Debug,Clone)]
struct AddedOneWall{
    faces:[FaceId;8],edges:[EdgeId;18],
}
impl PatchBuilder {
    fn add_one_wall(
        &mut self,coords:&[Point3;12],ta:f64,tb:f64,
        frame:Frame3,r:f64,specs:&[FaceSpec;8],
        first:&BrepModel,second:&BrepModel,
    )->Result<AddedOneWall,GeometryError>{
        let vbase=self.model.vertices.len();
        for &point in coords {
            let id=VertexId(self.model.vertices.len() as u32);
            self.model.vertices.push(Vertex {id,position:point});
        }
        let v=|index:usize|VertexId((vbase+index) as u32);
        let circle_endpoints=if ta<tb { [4usize,5usize] } else { [5usize,4usize] };
        let mk_circle=|y:f64|{
            Circle3 {frame:Frame3{origin:p(frame.origin.x,y,frame.origin.z),..frame},
                radius:r}
        };
        let circle_domain=bound(ta.min(tb),ta.max(tb))?;
        let edge_ids=[
            self.add_edge(v(0),v(1),None)?,
            self.add_edge(v(1),v(2),None)?,
            self.add_edge(v(2),v(3),None)?,
            self.add_edge(v(3),v(4),None)?,
            self.add_edge(v(circle_endpoints[0]),v(circle_endpoints[1]),
                Some((mk_circle(coords[0].y),circle_domain)))?,
            self.add_edge(v(5),v(0),None)?,
            self.add_edge(v(6),v(7),None)?,
            self.add_edge(v(7),v(8),None)?,
            self.add_edge(v(8),v(9),None)?,
            self.add_edge(v(9),v(10),None)?,
            self.add_edge(v(circle_endpoints[0]+6),v(circle_endpoints[1]+6),
                Some((mk_circle(coords[6].y),circle_domain)))?,
            self.add_edge(v(11),v(6),None)?,
            self.add_edge(v(0),v(6),None)?,
            self.add_edge(v(1),v(7),None)?,
            self.add_edge(v(2),v(8),None)?,
            self.add_edge(v(3),v(9),None)?,
            self.add_edge(v(4),v(10),None)?,
            self.add_edge(v(5),v(11),None)?,
        ];
        let edge_pairs=[
            (0,1),(1,2),(2,3),(3,4),
            (circle_endpoints[0],circle_endpoints[1]),(5,0),
            (6,7),(7,8),(8,9),(9,10),
            (circle_endpoints[0]+6,circle_endpoints[1]+6),(11,6),
            (0,6),(1,7),(2,8),(3,9),(4,10),(5,11),
        ];
        let shell=ShellId(self.model.shells.len() as u32);
        let mut face_ids=Vec::new();
        for (face_index,spec) in specs.iter().enumerate() {
            let source=if spec.operand==GraphOperand::First{first}else{second};
            let surface=source.face_surface(spec.source)
                .ok_or_else(||invalid("G4b2b3b missing source supporting surface"))?;
            if (face_index==7)!=matches!(surface,Surface3::Cylinder(_)){
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
                        u:if index==4||index==10{ta}else{tb},
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
                sense:if face_index==7{FaceSense::Reversed}else{FaceSense::Forward},
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
                        if edge_index==4||edge_index==10 {
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
                            let phase=if start==4||start==10{ta}else{tb};
                            let direction=if (start==4||start==10) == (ta<tb) {1.}else{-1.};
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
                            u:if start==4||start==10{ta}else{tb},
                            v:from.y-frame.origin.y,
                        };
                        let b=Point2Param {
                            u:if end==4||end==10{ta}else{tb},
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
        Ok(AddedOneWall {
            faces:face_ids.try_into().map_err(|_|invalid("missing corner faces"))?,
            edges:edge_ids,
        })
    }
}


#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticOneWallHold{
    UnsupportedOperation,
    UnsupportedOperand,
    NonCanonicalBox,
    UpstreamSourceSplit(AnalyticSourceEdgeSplitHold),
    NotSingleWallTransverse,
    AmbiguousMajorArc,
    MissingSourceFace,
    MissingOriginalEdge,
    MissingSourceCut,
    UnusedSourceArc,
    MaterialSideUncertain,
    InvalidShellTopology,
    UncertainAnalyticMass,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct OneWallPatchSource{
    pub result_face:FaceId,pub shell:ShellId,
    pub source_operand:GraphOperand,pub source_face:FaceId,
    pub reversed_for_difference:bool,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct OneWallMaterialWitness{
    pub face:FaceId,pub surface_point:Point3,
    pub retained_point:Point3,pub removed_point:Point3,
    pub outward:Point3,
}
#[derive(Debug,Clone,PartialEq)]
pub struct OneWallComponent{
    pub shell:ShellId,
    /// bottom cap, top cap, zmin, xmax, zmax, xleft upper/lower, cylinder
    pub faces:[FaceId;8],
    /// bottom: 6 edges; top: 6 edges; vertical: 6 edges.
    pub edges:[EdgeId;18],
    /// [bottom circle, top circle, upper generator, lower generator]
    pub new_cut_arcs:[GlobalArcRef;4],
    pub new_cut_events:[[GlobalEventId;2];4],
    pub major_arc_parameter_span:f64,
    pub exact_area_xz:f64,
    pub exact_volume:f64,
}
#[derive(Debug,Clone)]
pub struct AnalyticOneWallPlan {
    model:BrepModel,
    pub component:OneWallComponent,
    pub source_faces:Vec<OneWallPatchSource>,
    pub material_witnesses:Vec<OneWallMaterialWitness>,
    pub edge_splits:GlobalAnalyticEdgeSplitPlan,
    pub radial:RadialTopologyCertificate,
    pub closed_shell_proven:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone)]
pub enum AnalyticOneWallResult {
    Ready(AnalyticOneWallPlan),
    Held(AnalyticOneWallHold),
}
fn strict_material(
    sample:Point3,bounds:crate::BoundingBox,cx:f64,cz:f64,
    r:f64,eps:f64,
)->Option<bool>{
    if !sample.is_finite(){return None;}
    let margins=[
        sample.x-bounds.min.x,bounds.max.x-sample.x,
        sample.y-bounds.min.y,bounds.max.y-sample.y,
        sample.z-bounds.min.z,bounds.max.z-sample.z,
    ];
    let radial=(sample.x-cx).hypot(sample.z-cz);
    if margins.iter().any(|m|*m < -eps)||radial<r-eps{
        return Some(false);
    }
    if margins.iter().all(|m|*m>eps)&&radial>r+eps{
        return Some(true);
    }
    None
}
impl AnalyticOneWallPlan {
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        self.edge_splits.validate_against(first,second,tol)?;
        let AnalyticOneWallResult::Ready(rebuilt)=
            Self::build(first,second,BooleanOperation::Difference,tol)? else{
                return Err(invalid("G7a source graph no longer forms one-wall notch"));
            };
        if self.component!=rebuilt.component
            ||self.source_faces!=rebuilt.source_faces
            ||self.material_witnesses!=rebuilt.material_witnesses
            ||self.radial!=rebuilt.radial
            ||self.closed_shell_proven!=rebuilt.closed_shell_proven
            ||self.can_apply!=rebuilt.can_apply
            ||format!("{:?}",self.model)!=format!("{:?}",rebuilt.model){
            return Err(invalid("G7a physical shell differs from original sources"));
        }
        self.model.validate()
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        operation:BooleanOperation,tol:GeometryTolerance,
    )->Result<AnalyticOneWallResult,GeometryError>{
        tol.validate()?;first.validate()?;second.validate()?;
        if operation!=BooleanOperation::Difference{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::UnsupportedOperation));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cyl))=
            (first,second) else{
                return Ok(AnalyticOneWallResult::Held(
                    AnalyticOneWallHold::UnsupportedOperand));
        };
        match crate::boolean::certify_box(block){
            Ok(())=>{},
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::NonCanonicalBox)),
        }
        let frame=cyl.frame;
        let bounds=block.bbox;
        let scale=bounds.max.sub(bounds.min).length_squared().sqrt()
            .max(cyl.radius).max(cyl.height).max(1.);
        let eps=tol.length_at(scale)*64.;
        if distance(frame.x,p(1.,0.,0.))>tol.angular*16.
            ||distance(frame.y,p(0.,0.,-1.))>tol.angular*16.
            ||distance(frame.normal,p(0.,1.,0.))>tol.angular*16. {
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::NotSingleWallTransverse));
        }
        let (xmin,xmax,zmin,zmax)=(
            bounds.min.x,bounds.max.x,bounds.min.z,bounds.max.z
        );
        let cx=frame.origin.x;let cz=frame.origin.z;let r=cyl.radius;
        let ymin=bounds.min.y;let ymax=bounds.max.y;
        let height=ymax-ymin;
        let d=cx-xmin;
        if ![d,r,height,xmax,zmin,zmax,cx,cz].iter().all(|v|finite(*v))
            ||d<=eps ||r<=d+eps
            ||xmax-cx<=r+eps
            ||cz-zmin<=r+eps
            ||zmax-cz<=r+eps
            ||frame.origin.y>=ymin-eps
            ||frame.origin.y+cyl.height<=ymax+eps{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::NotSingleWallTransverse));
        }
        let zcut=((r-d)*(r+d)).sqrt();
        if !finite(zcut)||zcut<=eps*16.
            ||zmax-cz-zcut<=eps*32.
            ||cz-zmin-zcut<=eps*32.{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::AmbiguousMajorArc));
        }
        let ta=(-(zcut)).atan2(-d);
        let tb=zcut.atan2(-d);
        let span=tb-ta;
        if !(span>std::f64::consts::PI+tol.angular*32.
            &&span<std::f64::consts::TAU-tol.angular*32.){
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::AmbiguousMajorArc));
        }
        let split=match GlobalAnalyticEdgeSplitPlan::build(first,second,tol)?{
            GlobalAnalyticEdgeSplitResult::Ready(x)=>x,
            GlobalAnalyticEdgeSplitResult::Held(reason)=>
                return Ok(AnalyticOneWallResult::Held(
                    AnalyticOneWallHold::UpstreamSourceSplit(reason))),
        };
        let a=first.shared()?;
        let b=second.shared()?;
        let cyl_faces=b.faces.iter()
            .filter(|f|matches!(b.face_surface(f.id),Some(Surface3::Cylinder(_))))
            .map(|f|f.id).collect::<Vec<_>>();
        if cyl_faces.len()!=1{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::UnsupportedOperand));
        }
        let wall_cyl=cyl_faces[0];
        let bottom=find_source_plane(&a,1,ymin,-1.,eps);
        let top=find_source_plane(&a,1,ymax,1.,eps);
        let zlow=find_source_plane(&a,2,zmin,-1.,eps);
        let zhigh=find_source_plane(&a,2,zmax,1.,eps);
        let xright=find_source_plane(&a,0,xmax,1.,eps);
        let xleft=find_source_plane(&a,0,xmin,-1.,eps);
        let (Some(bottom),Some(top),Some(zlow),Some(zhigh),
            Some(xright),Some(xleft))=
            (bottom,top,zlow,zhigh,xright,xleft) else{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::MissingSourceFace));
        };
        let area_box=(xmax-xmin)*(zmax-zmin);
        let cap=r*r*(d/r).acos()-d*zcut;
        let inside_disk=std::f64::consts::PI*r*r-cap;
        let area=area_box-inside_disk;
        let moment_u=2.*zcut*zcut*zcut/3.;
        let centroid_x=((xmin+xmax)/2.*area_box
            -(cx*inside_disk+moment_u))/area;
        let centroid_z=((zmin+zmax)/2.*area_box-cz*inside_disk)/area;
        let side_perimeter=2.*((xmax-xmin)+(zmax-zmin))
            -2.*zcut+r*span;
        if ![cap,inside_disk,area,centroid_x,centroid_z,side_perimeter]
            .iter().all(|v|finite(*v))
            ||cap<=eps*eps
            ||inside_disk<=0.||area<=eps*eps*32.
            ||side_perimeter<=0.{
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::UncertainAnalyticMass));
        }
        let mut builder=PatchBuilder::new(bounds,tol,a.units);
        builder.model.origin=BrepOrigin::AnalyticOneWallCylinderDifference;
        let lower=p(xmin,ymin,zmin);
        let lower_right=p(xmax,ymin,zmin);
        let upper_right=p(xmax,ymin,zmax);
        let upper=p(xmin,ymin,zmax);
        let upper_cut=p(xmin,ymin,cz+zcut);
        let lower_cut=p(xmin,ymin,cz-zcut);
        let coords=[
            lower,lower_right,upper_right,upper,upper_cut,lower_cut,
            p(lower.x,ymax,lower.z),p(lower_right.x,ymax,lower_right.z),
            p(upper_right.x,ymax,upper_right.z),p(upper.x,ymax,upper.z),
            p(upper_cut.x,ymax,upper_cut.z),p(lower_cut.x,ymax,lower_cut.z),
        ];
        let specs=[
            FaceSpec{local:vec![0,1,2,3,4,5],source:bottom,
                operand:GraphOperand::First,normal:p(0.,-1.,0.)},
            FaceSpec{local:vec![6,7,8,9,10,11],source:top,
                operand:GraphOperand::First,normal:p(0.,1.,0.)},
            FaceSpec{local:vec![0,1,7,6],source:zlow,
                operand:GraphOperand::First,normal:p(0.,0.,-1.)},
            FaceSpec{local:vec![1,2,8,7],source:xright,
                operand:GraphOperand::First,normal:p(1.,0.,0.)},
            FaceSpec{local:vec![2,3,9,8],source:zhigh,
                operand:GraphOperand::First,normal:p(0.,0.,1.)},
            FaceSpec{local:vec![3,4,10,9],source:xleft,
                operand:GraphOperand::First,normal:p(-1.,0.,0.)},
            FaceSpec{local:vec![5,0,6,11],source:xleft,
                operand:GraphOperand::First,normal:p(-1.,0.,0.)},
            FaceSpec{local:vec![4,5,11,10],source:wall_cyl,
                operand:GraphOperand::Second,normal:p(-1.,0.,0.)},
        ];
        let output=builder.add_one_wall(
            &coords,ta,tb,frame,r,&specs,&a,&b)?;
        let orig=[
            (0usize,[bottom,zlow]),(1,[bottom,xright]),
            (2,[bottom,zhigh]),(3,[bottom,xleft]),(5,[bottom,xleft]),
            (6,[top,zlow]),(7,[top,xright]),(8,[top,zhigh]),
            (9,[top,xleft]),(11,[top,xleft]),
            (12,[zlow,xleft]),(13,[zlow,xright]),
            (14,[zhigh,xright]),(15,[zhigh,xleft]),
        ];
        let mut used_original=BTreeSet::<EdgeId>::new();
        for (index,faces) in orig{
            let id=output.edges[index];
            let edge=&builder.model.edges[id.0 as usize];
            if !certify_original_material_edge(
                &split.split_sources[0],edge,&builder.model.vertices,
                faces,eps,&mut used_original,
            ){
                return Ok(AnalyticOneWallResult::Held(
                    AnalyticOneWallHold::MissingOriginalEdge));
            }
        }
        let mut arcs=BTreeSet::<GlobalArcRef>::new();
        let mut endpoints=BTreeMap::<GlobalEventId,VertexId>::new();
        let mut arc_refs=Vec::new();
        let mut events=Vec::new();
        for (edge_index,first_face) in [
            (4usize,bottom),(10,top),(16,xleft),(17,xleft)
        ]{
            let eid=output.edges[edge_index];
            let edge=&builder.model.edges[eid.0 as usize];
            let carrier=builder.model.edge_curve(eid)
                .ok_or_else(||invalid("G7a missing source-certified cut curve"))?;
            let domain=builder.model.geometry.edges[eid.0 as usize].domain;
            let ends=[
                builder.model.vertices[edge.start.0 as usize].position,
                builder.model.vertices[edge.end.0 as usize].position,
            ];
            let Some((reference,ids))=certify_new_cut(
                &split.registry,first_face,wall_cyl,
                carrier,domain,ends,[edge.start,edge.end],eps,
                &mut arcs,&mut endpoints,
            )? else{
                return Ok(AnalyticOneWallResult::Held(
                    AnalyticOneWallHold::MissingSourceCut));
            };
            arc_refs.push(reference);events.push(ids);
        }
        if arcs.len()!=split.registry.arcs.len()
            ||endpoints.len()!=split.registry.physical_event_count {
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::UnusedSourceArc));
        }
        let centre_y=(ymin+ymax)/2.;
        let farx=(cx+r+xmax)/2.;
        let side_upper=(cz+zcut+zmax)/2.;
        let side_lower=(zmin+cz-zcut)/2.;
        let sample=[
            p(farx,ymin,cz),p(farx,ymax,cz),
            p(cx,centre_y,zmin),p(xmax,centre_y,cz),
            p(cx,centre_y,zmax),p(xmin,centre_y,side_upper),
            p(xmin,centre_y,side_lower),p(cx+r,centre_y,cz),
        ];
        let normals=[
            p(0.,-1.,0.),p(0.,1.,0.),p(0.,0.,-1.),
            p(1.,0.,0.),p(0.,0.,1.),p(-1.,0.,0.),
            p(-1.,0.,0.),p(-1.,0.,0.),
        ];
        let step=eps*32.;
        if (xmax-cx-r)<=step*8.
            ||(zmax-cz-zcut)<=step*8.
            ||(cz-zmin-zcut)<=step*8. {
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::MaterialSideUncertain));
        }
        let mut witnesses=Vec::new();
        let mut sources=Vec::new();
        for (i,&fid) in output.faces.iter().enumerate(){
            let point=sample[i];
            let n=normals[i];
            let retained=point.sub(n.scale(step));
            let removed=point.add(n.scale(step));
            if strict_material(retained,bounds,cx,cz,r,eps)!=Some(true)
                ||strict_material(removed,bounds,cx,cz,r,eps)!=Some(false){
                return Ok(AnalyticOneWallResult::Held(
                    AnalyticOneWallHold::MaterialSideUncertain));
            }
            witnesses.push(OneWallMaterialWitness{
                face:fid,surface_point:point,
                retained_point:retained,removed_point:removed,outward:n,
            });
            sources.push(OneWallPatchSource{
                result_face:fid,shell:ShellId(0),
                source_operand:specs[i].operand,
                source_face:specs[i].source,
                reversed_for_difference:i==7,
            });
        }
        builder.model.mass=MassProperties{
            volume:area*height,
            surface_area:2.*area+side_perimeter*height,
            centroid:p(centroid_x,centre_y,centroid_z),
        };
        let result=builder.finish()?;
        let radial=radial_certificate(&result,eps)
            .map_err(|e|invalid(format!("G7a radial manifold validation: {e:?}")))?;
        if result.shells.len()!=1||result.faces.len()!=8
            ||result.loops.len()!=8||result.vertices.len()!=12
            ||result.edges.len()!=18||result.coedges.len()!=36
            ||result.euler_characteristic()!=2
            ||radial.shell_components!=1||radial.total_genus!=0
            ||radial.euler_characteristic!=2 {
            return Ok(AnalyticOneWallResult::Held(
                AnalyticOneWallHold::InvalidShellTopology));
        }
        Ok(AnalyticOneWallResult::Ready(Self{
            model:result,component:OneWallComponent{
                shell:ShellId(0),faces:output.faces,edges:output.edges,
                new_cut_arcs:arc_refs.try_into().map_err(|_|
                    invalid("G7a does not have four source cut carriers"))?,
                new_cut_events:events.try_into().map_err(|_|
                    invalid("G7a missing physical G1 cut endpoints"))?,
                major_arc_parameter_span:span,
                exact_area_xz:area,exact_volume:area*height,
            },
            source_faces:sources,material_witnesses:witnesses,
            edge_splits:split,radial,closed_shell_proven:true,
            can_apply:false,
        }))
    }
}
