//! Canonical analytic B-rep v2; no Mesh, MassProperties, feature history or producer tags.
//! Local incidence and sampled 3D/UV bindings are checked here. Boolean material
//! correctness, global self-intersections and independent volume belong to validation.

use crate::{CoedgeId, EdgeId, FaceId, LoopId, LoopRole, ShellId, VertexId};
use rustsolid_foundation::{GeometryError, GeometryTolerance, Point3};
use rustsolid_geometry::{Curve2, Curve3, ModelUnits, ParameterRange, Surface3};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::f64::consts::TAU;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve2Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve3Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SurfaceId(pub u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalBrepIncarnation(u64);
impl LocalBrepIncarnation {
    pub fn new(value: u64) -> Result<Self, BrepIssue> {
        if value == 0 { return Err(BrepIssue::Invalid("zero incarnation")); }
        Ok(Self(value))
    }
    pub fn get(self)->u64 {self.0}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeKind { Physical, PeriodicSeam }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceSense { Forward, Reversed }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellSense { Outward, Cavity }

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineParameterMap { pub scale:f64, pub offset:f64 }
impl AffineParameterMap {
    pub fn map(self, physical_edge_t:f64)->Result<f64,BrepIssue>{
        if !self.scale.is_finite()||self.scale==0.0||!self.offset.is_finite(){
            return Err(BrepIssue::Invalid("invalid coedge parameter map"));
        }
        let value=self.scale*physical_edge_t+self.offset;
        if !value.is_finite(){return Err(BrepIssue::Invalid("parameter mapping overflow"));}
        Ok(value)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepVertex { pub id:VertexId, pub point:Point3 }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepEdge {
    pub id:EdgeId,
    pub vertices:[VertexId;2],
    pub curve:Curve3Id,
    pub domain:ParameterRange,
    pub uses:[CoedgeId;2],
    pub kind:EdgeKind,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepCoedge {
    pub id:CoedgeId, pub edge:EdgeId, pub face:FaceId, pub loop_id:LoopId,
    pub pcurve:Curve2Id,
    /// Maps the PHYSICAL 3D curve parameter to this face's Curve2 parameter.
    /// Loop reversal is stored independently in reversed.
    pub edge_to_pcurve:AffineParameterMap,
    pub reversed:bool, pub next:CoedgeId, pub prev:CoedgeId, pub twin:CoedgeId,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepLoop {
    pub id:LoopId, pub face:FaceId, pub role:LoopRole, pub first:CoedgeId,
}
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFace {
    pub id:FaceId, pub surface:SurfaceId,
    pub loops:Vec<LoopId>, pub sense:FaceSense,
}
#[derive(Debug, Clone, PartialEq)]
pub struct BrepShell {
    pub id:ShellId, pub faces:Vec<FaceId>, pub sense:ShellSense,
}
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeometryCatalog {
    pub curves3:Vec<Curve3>,
    pub curves2:Vec<Curve2>,
    pub surfaces:Vec<Surface3>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BrepIssue {
    Invalid(&'static str),
    Geometry(GeometryError),
    CarrierMismatch{coedge:CoedgeId, distance:f64, tolerance:f64},
}
impl fmt::Display for BrepIssue {
    fn fmt(&self,f:&mut fmt::Formatter<'_>)->fmt::Result{write!(f,"{self:?}")}
}
impl std::error::Error for BrepIssue {}
impl From<GeometryError> for BrepIssue {
    fn from(v:GeometryError)->Self{Self::Geometry(v)}
}

/// The ONLY mutable modeling input. Invalid drafts cannot be converted to BrepModel.
#[derive(Debug, Clone)]
pub struct BrepDraft {
    pub incarnation:LocalBrepIncarnation,
    pub units:ModelUnits,
    pub tolerance:GeometryTolerance,
    pub vertices:Vec<BrepVertex>,
    pub edges:Vec<BrepEdge>,
    pub coedges:Vec<BrepCoedge>,
    pub loops:Vec<BrepLoop>,
    pub faces:Vec<BrepFace>,
    pub shells:Vec<BrepShell>,
    pub carriers:GeometryCatalog,
}
impl BrepDraft {
    pub fn check(self)->Result<BrepModel,BrepIssue>{
        let model=BrepModel{
            incarnation:self.incarnation, units:self.units, tolerance:self.tolerance,
            vertices:self.vertices,edges:self.edges,coedges:self.coedges,
            loops:self.loops,faces:self.faces,shells:self.shells,carriers:self.carriers,
        };
        model.validate()?;
        Ok(model)
    }
}

/// Immutable, incidence-checked analytic shape. This is NOT Boolean certification.
#[derive(Debug, Clone)]
pub struct BrepModel {
    incarnation:LocalBrepIncarnation,
    units:ModelUnits,
    tolerance:GeometryTolerance,
    vertices:Vec<BrepVertex>,
    edges:Vec<BrepEdge>,
    coedges:Vec<BrepCoedge>,
    loops:Vec<BrepLoop>,
    faces:Vec<BrepFace>,
    shells:Vec<BrepShell>,
    carriers:GeometryCatalog,
}
impl BrepModel {
    pub fn incarnation(&self)->LocalBrepIncarnation{self.incarnation}
    pub fn units(&self)->ModelUnits{self.units}
    pub fn tolerance(&self)->GeometryTolerance{self.tolerance}
    pub fn vertices(&self)->&[BrepVertex]{&self.vertices}
    pub fn edges(&self)->&[BrepEdge]{&self.edges}
    pub fn coedges(&self)->&[BrepCoedge]{&self.coedges}
    pub fn loops(&self)->&[BrepLoop]{&self.loops}
    pub fn faces(&self)->&[BrepFace]{&self.faces}
    pub fn shells(&self)->&[BrepShell]{&self.shells}
    pub fn carriers(&self)->&GeometryCatalog{&self.carriers}
    pub fn validate(&self)->Result<(),BrepIssue>{
        self.tolerance.validate()?;
        if self.vertices.is_empty()||self.edges.is_empty()||self.coedges.is_empty()
            ||self.faces.is_empty()||self.loops.is_empty()||self.shells.is_empty(){
            return Err(BrepIssue::Invalid("empty canonical B-rep graph"));
        }
        // This compact format deliberately requires dense, incidence-local IDs;
        // stable CAD feature identities are maintained in the document layer.
        for (i,v) in self.vertices.iter().enumerate(){
            if v.id.0 as usize!=i ||!v.point.is_finite(){
                return Err(BrepIssue::Invalid("invalid vertex identity or coordinates"));
            }
        }
        for (i,e) in self.edges.iter().enumerate(){
            if e.id.0 as usize!=i {return Err(BrepIssue::Invalid("edge IDs not dense"));}
        }
        for (i,c) in self.coedges.iter().enumerate(){
            if c.id.0 as usize!=i {return Err(BrepIssue::Invalid("coedge IDs not dense"));}
        }
        for (i,l) in self.loops.iter().enumerate(){
            if l.id.0 as usize!=i {return Err(BrepIssue::Invalid("loop IDs not dense"));}
        }
        for (i,f) in self.faces.iter().enumerate(){
            if f.id.0 as usize!=i {return Err(BrepIssue::Invalid("face IDs not dense"));}
        }
        for (i,s) in self.shells.iter().enumerate(){
            if s.id.0 as usize!=i {return Err(BrepIssue::Invalid("shell IDs not dense"));}
        }

        let mut face_shell=vec![None;self.faces.len()];
        for shell in &self.shells {
            if shell.faces.is_empty(){return Err(BrepIssue::Invalid("empty shell"));}
            for face_id in &shell.faces {
                let slot=face_shell.get_mut(face_id.0 as usize)
                    .ok_or(BrepIssue::Invalid("shell references missing face"))?;
                if slot.replace(shell.id).is_some(){
                    return Err(BrepIssue::Invalid("face belongs to multiple shells"));
                }
            }
        }
        if face_shell.iter().any(Option::is_none) {
            return Err(BrepIssue::Invalid("face unassigned to shell"));
        }

        let mut loop_owner=vec![0usize;self.loops.len()];
        for face in &self.faces {
            if face.loops.is_empty() || face.surface.0 as usize>=self.carriers.surfaces.len(){
                return Err(BrepIssue::Invalid("face missing loop or analytic surface"));
            }
            let mut outer=0usize;
            for loop_id in &face.loops {
                let slot=loop_owner.get_mut(loop_id.0 as usize)
                    .ok_or(BrepIssue::Invalid("face references nonexistent loop"))?;
                let loop_item=&self.loops[loop_id.0 as usize];
                if loop_item.face!=face.id {
                    return Err(BrepIssue::Invalid("loop attached to wrong face"));
                }
                if loop_item.role==LoopRole::Outer {outer+=1;}
                *slot+=1;
            }
            if outer!=1{return Err(BrepIssue::Invalid("face requires one outer loop"));}
        }
        if loop_owner.iter().any(|&n|n!=1){
            return Err(BrepIssue::Invalid("loop not owned by exactly one face"));
        }

        let mut coedge_uses=vec![0usize;self.coedges.len()];
        for loop_item in &self.loops {
            let mut current=loop_item.first;
            // A valid cycle cannot have more uses than the entire coedge catalog.
            let mut visited=BTreeSet::new();
            loop {
                let c=self.coedges.get(current.0 as usize)
                    .ok_or(BrepIssue::Invalid("loop coedge out of bounds"))?;
                if c.loop_id!=loop_item.id || c.face!=loop_item.face {
                    return Err(BrepIssue::Invalid("coedge owner mismatch"));
                }
                if !visited.insert(current){
                    return Err(BrepIssue::Invalid("loop has repeated coedge before closure"));
                }
                coedge_uses[current.0 as usize]+=1;
                let next=self.coedges.get(c.next.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge.next out of bounds"))?;
                let prev=self.coedges.get(c.prev.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge.prev out of bounds"))?;
                if next.prev!=c.id || prev.next!=c.id {
                    return Err(BrepIssue::Invalid("next/prev incidence inconsistent"));
                }
                let current_edge=self.edges.get(c.edge.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge refers to missing edge"))?;
                let next_edge=self.edges.get(next.edge.0 as usize)
                    .ok_or(BrepIssue::Invalid("next coedge refers to missing edge"))?;
                let finish=if c.reversed {current_edge.vertices[0]} else {current_edge.vertices[1]};
                let next_start=if next.reversed {next_edge.vertices[1]} else {next_edge.vertices[0]};
                if finish!=next_start {return Err(BrepIssue::Invalid("loop is not vertex-contiguous"));}
                current=c.next;
                if current==loop_item.first {break;}
                if visited.len()>=self.coedges.len(){
                    return Err(BrepIssue::Invalid("loop cycle exceeds coedge catalog"));
                }
            }
        }
        if coedge_uses.iter().any(|&n|n!=1){
            return Err(BrepIssue::Invalid("coedge missing or shared by multiple loops"));
        }

        // A shell is a connected face component. Edge incidence cannot connect
        // two independent components; a separate vertex-only touch is also rejected.
        let mut adjacent:Vec<BTreeSet<FaceId>>=vec![BTreeSet::new();self.faces.len()];
        let mut incident_faces:Vec<BTreeSet<FaceId>>=vec![BTreeSet::new();self.vertices.len()];
        let mut vertex_links:Vec<Vec<(FaceId,FaceId)>>=vec![Vec::new();self.vertices.len()];
        for edge in &self.edges {
            let [start,end]=edge.vertices;
            let a=self.vertices.get(start.0 as usize)
                .ok_or(BrepIssue::Invalid("edge start vertex out of bounds"))?;
            let b=self.vertices.get(end.0 as usize)
                .ok_or(BrepIssue::Invalid("edge end vertex out of bounds"))?;
            if edge.uses[0]==edge.uses[1]{
                return Err(BrepIssue::Invalid("edge must have two distinct coedges"));
            }
            let u=self.coedges.get(edge.uses[0].0 as usize)
                .ok_or(BrepIssue::Invalid("edge.use[0] out of bounds"))?;
            let v=self.coedges.get(edge.uses[1].0 as usize)
                .ok_or(BrepIssue::Invalid("edge.use[1] out of bounds"))?;
            if u.edge!=edge.id || v.edge!=edge.id || u.twin!=v.id || v.twin!=u.id
                || u.reversed==v.reversed {
                return Err(BrepIssue::Invalid("edge uses/twins are not opposite"));
            }
            match edge.kind {
                EdgeKind::PeriodicSeam if u.face!=v.face =>
                    return Err(BrepIssue::Invalid("periodic seam must be same-face")),
                EdgeKind::Physical if u.face==v.face =>
                    return Err(BrepIssue::Invalid("same-face seam cannot be physical cut")),
                _=>{},
            }
            if face_shell[u.face.0 as usize]!=face_shell[v.face.0 as usize] {
                return Err(BrepIssue::Invalid("edge crosses independent shells"));
            }
            adjacent[u.face.0 as usize].insert(v.face);
            adjacent[v.face.0 as usize].insert(u.face);
            for vertex in [start,end] {
                incident_faces[vertex.0 as usize].insert(u.face);
                incident_faces[vertex.0 as usize].insert(v.face);
                vertex_links[vertex.0 as usize].push((u.face,v.face));
            }

            let curve=*self.carriers.curves3.get(edge.curve.0 as usize)
                .ok_or(BrepIssue::Invalid("edge has missing analytic curve"))?;
            let (t0,t1)=match edge.domain {
                ParameterRange::Bounded {start,end}
                    if start.is_finite() && end.is_finite() && start<end => (start,end),
                _=>return Err(BrepIssue::Invalid("physical edge requires finite domain"))
            };
            if start==end {
                if !matches!(curve,Curve3::Circle(_) | Curve3::Harmonic(_)) ||
                    ((t1-t0)-TAU).abs()>self.tolerance.angular.max(1e-9) {
                    return Err(BrepIssue::Invalid("closed edge must have a full analytic revolution"));
                }
            }
            let scale=a.point.sub(b.point);
            let extent=scale.x.hypot(scale.y).hypot(scale.z)
                ; // origin-independent local extent
            let budget=self.tolerance.length_at(extent.max(1.0))*32.0;
            for (point,t) in [(a.point,t0),(b.point,t1)] {
                let model_point=curve.evaluate(t)?;
                let delta=model_point.sub(point);
                let residual=delta.x.hypot(delta.y).hypot(delta.z);
                if !residual.is_finite() || residual>budget {
                    return Err(BrepIssue::Invalid("curve endpoint does not match physical vertex"));
                }
            }
            for c in [u,v]{
                let chart=*self.carriers.curves2.get(c.pcurve.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge missing planar pcurve"))?;
                let face=self.faces.get(c.face.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge missing owning face"))?;
                let surface=*self.carriers.surfaces.get(face.surface.0 as usize)
                    .ok_or(BrepIssue::Invalid("coedge face missing support surface"))?;
                for alpha in [0.0,0.25,0.5,0.75,1.0]{
                    let t=t0+(t1-t0)*alpha;
                    let xyz=curve.evaluate(t)?;
                    let param=c.edge_to_pcurve.map(t)?;
                    let uv=chart.evaluate(param)?;
                    let back=surface.evaluate(uv.u,uv.v)?.point;
                    let d=xyz.sub(back);
                    let residual=d.x.hypot(d.y).hypot(d.z);
                    if !residual.is_finite() || residual>budget {
                        return Err(BrepIssue::CarrierMismatch{
                            coedge:c.id,distance:residual,tolerance:budget,
                        });
                    }
                }
            }
        }

        for (i,faces) in incident_faces.iter().enumerate(){
            if faces.is_empty(){return Err(BrepIssue::Invalid("orphan B-rep vertex"));}
            let first=*faces.iter().next().expect("nonempty");
            let mut reached=BTreeSet::from([first]);
            let mut stack=vec![first];
            while let Some(face)=stack.pop(){
                for &(a,b) in &vertex_links[i]{
                    let neighbor=if face==a {Some(b)}else if face==b {Some(a)}else{None};
                    if let Some(id)=neighbor {
                        if reached.insert(id){stack.push(id);}
                    }
                }
            }
            if &reached!=faces {
                return Err(BrepIssue::Invalid("nonmanifold/disconnected vertex fan"));
            }
        }
        for shell in &self.shells {
            let face_set:BTreeSet<_>=shell.faces.iter().copied().collect();
            let mut reached=BTreeSet::from([shell.faces[0]]);
            let mut stack=vec![shell.faces[0]];
            while let Some(face)=stack.pop(){
                for &n in &adjacent[face.0 as usize]{
                    if face_set.contains(&n) && reached.insert(n){stack.push(n);}
                }
            }
            if reached!=face_set {
                return Err(BrepIssue::Invalid("shell contains disconnected face components"));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsolid_geometry::{Frame3, Line2, Line3, PlaneSurface, Point2Param};
    fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
    fn tetrahedron()->BrepDraft {
        let tol=GeometryTolerance::default();
        let positions=[p(0.,0.,0.),p(1.,0.,0.),p(0.,1.,0.),p(0.,0.,1.)];
        let vertices=positions.iter().enumerate().map(|(i,&point)|
            BrepVertex{id:VertexId(i as u32),point}).collect::<Vec<_>>();
        let triangles=[[0u32,2,1],[0,1,3],[1,2,3],[2,0,3]];
        let mut catalog=GeometryCatalog::default();
        let mut faces=Vec::new();
        let mut loops=Vec::new();
        let mut coedges=Vec::new();
        let mut incidences:BTreeMap<(u32,u32),Vec<CoedgeId>>=BTreeMap::new();
        for (i,ring) in triangles.iter().enumerate(){
            let origin=positions[ring[0] as usize];
            let x=positions[ring[1] as usize].sub(origin);
            let y=positions[ring[2] as usize].sub(origin);
            let frame=Frame3::from_axes(origin,x,y,tol).unwrap();
            let sid=SurfaceId(catalog.surfaces.len() as u32);
            catalog.surfaces.push(Surface3::Plane(PlaneSurface{frame}));
            let lid=LoopId(i as u32);
            let first=CoedgeId(coedges.len() as u32);
            loops.push(BrepLoop{id:lid,face:FaceId(i as u32),role:LoopRole::Outer,first});
            faces.push(BrepFace{id:FaceId(i as u32),surface:sid,loops:vec![lid],sense:FaceSense::Forward});
            for j in 0..3 {
                let a=ring[j];let b=ring[(j+1)%3];
                let key=(a.min(b),a.max(b));
                let id=CoedgeId(coedges.len() as u32);
                coedges.push(BrepCoedge{
                    id,edge:EdgeId(u32::MAX),face:FaceId(i as u32),loop_id:lid,
                    pcurve:Curve2Id(u32::MAX),
                    edge_to_pcurve:AffineParameterMap{scale:1.0,offset:0.0},
                    reversed:a>b,
                    next:CoedgeId(first.0+((j+1)%3) as u32),
                    prev:CoedgeId(first.0+((j+2)%3) as u32),
                    twin:CoedgeId(u32::MAX),
                });
                incidences.entry(key).or_default().push(id);
            }
        }
        let mut edges=Vec::new();
        for ((start,end),uses) in incidences {
            assert_eq!(uses.len(),2);
            let id=EdgeId(edges.len() as u32);
            let (line,domain)=Line3::between(positions[start as usize],positions[end as usize],tol).unwrap();
            let span=match domain {ParameterRange::Bounded{end,..}=>end,_=>unreachable!()};
            catalog.curves3.push(Curve3::Line(line));
            edges.push(BrepEdge{
                id,vertices:[VertexId(start),VertexId(end)],
                curve:Curve3Id(catalog.curves3.len() as u32-1),
                domain,uses:[uses[0],uses[1]],kind:EdgeKind::Physical,
            });
            for (&cid,&other) in uses.iter().zip(uses.iter().rev()){
                let c=&mut coedges[cid.0 as usize];
                c.edge=id;c.twin=other;
                let f=&faces[c.face.0 as usize];
                let Surface3::Plane(s)=catalog.surfaces[f.surface.0 as usize] else {unreachable!()};
                let uv0=s.frame.project(positions[start as usize]).unwrap();
                let uv1=s.frame.project(positions[end as usize]).unwrap();
                let uv_delta=Point2Param{u:uv1.u-uv0.u,v:uv1.v-uv0.v};
                c.pcurve=Curve2Id(catalog.curves2.len() as u32);
                c.edge_to_pcurve=AffineParameterMap{scale:1.0/span,offset:0.0};
                catalog.curves2.push(Curve2::Line(Line2{origin:uv0,delta:uv_delta}));
            }
        }
        BrepDraft{
            incarnation:LocalBrepIncarnation::new(1).unwrap(),
            units:ModelUnits::default(),tolerance:tol,vertices,edges,coedges,loops,faces,
            shells:vec![BrepShell{id:ShellId(0),faces:(0..4).map(FaceId).collect(),sense:ShellSense::Outward}],
            carriers:catalog,
        }
    }
    #[test]
    fn tetrahedral_analytic_brep_has_no_mesh_or_mass_state() {
        let model=tetrahedron().check().unwrap();
        assert_eq!((model.vertices().len(),model.edges().len(),model.coedges().len(),
            model.loops().len(),model.faces().len(),model.shells().len()),(4,6,12,4,4,1));
        assert!(model.validate().is_ok());
        assert!(model.carriers().curves3.iter().all(|c|matches!(c,Curve3::Line(_))));
    }
    #[test]
    fn incomplete_uv_catalog_is_rejected() {
        let mut d=tetrahedron();
        d.coedges[0].pcurve=Curve2Id(999);
        assert!(d.check().is_err());
    }
    #[test]
    fn falsified_uv_carrier_is_rejected() {
        let mut d=tetrahedron();
        d.carriers.curves2[0]=Curve2::Line(Line2{
            origin:Point2Param{u:98.,v:0.},
            delta:Point2Param{u:0.,v:0.},
        });
        assert!(matches!(d.check(),Err(BrepIssue::CarrierMismatch{..})));
    }
    #[test]
    fn incomplete_twins_and_nonmanifold_edges_are_rejected() {
        let mut d=tetrahedron();
        d.coedges[0].twin=CoedgeId(0);
        assert!(d.check().is_err());
        let mut d=tetrahedron();
        d.edges[0].uses=[CoedgeId(0),CoedgeId(0)];
        assert!(d.check().is_err());
    }
    #[test]
    fn orphan_or_disconnected_faces_are_rejected() {
        let mut d=tetrahedron();
        d.shells[0].faces.pop();
        assert!(d.check().is_err());
        let mut d=tetrahedron();
        d.coedges[0].next=CoedgeId(11);
        assert!(d.check().is_err());
    }
    #[test]
    fn invalid_same_face_cut_is_not_treated_as_a_periodic_seam() {
        let mut d=tetrahedron();
        let e=d.edges[0].uses[1];
        d.coedges[e.0 as usize].face=d.coedges[d.edges[0].uses[0].0 as usize].face;
        assert!(d.check().is_err());
    }
}
