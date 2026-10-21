//! G7/G8 direct event-side assembly of the new canonical BrepDraft.
//! Events are recorded when their physical source geometry is constructed;
//! this module never converts a finished legacy BrepModel.
//! Retire the old builders and this dual-write scaffold after G4 consolidation.
use crate::{GeometryError,GeometryTolerance,Point3,ModelUnits,FaceSense,Curve2,Curve3,
    ParameterRange,Surface3,VertexId,EdgeId,CoedgeId,LoopId,LoopRole,FaceId,ShellId};
use rustsolid_topology as topo;
use std::sync::atomic::{AtomicU64,Ordering};

static NEXT_CANONICAL:AtomicU64=AtomicU64::new(1);
fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn extent(a:Point3,b:Point3)->f64{
    let d=a.sub(b); d.x.hypot(d.y).hypot(d.z)
}
fn bound(r:ParameterRange)->Result<(f64,f64),GeometryError>{
    match r{
        ParameterRange::Bounded{start,end} if start.is_finite()&&end.is_finite()&&start<end=>Ok((start,end)),
        _=>Err(invalid("canonical direct assembly requires finite ordered carrier domains")),
    }
}
#[derive(Debug)]
pub(crate) struct CanonicalDraftEmitter {
    draft:topo::BrepDraft,
}
impl CanonicalDraftEmitter {
    pub(crate) fn new(tolerance:GeometryTolerance,units:ModelUnits)->Result<Self,GeometryError>{
        let incarnation=topo::LocalBrepIncarnation::new(
            NEXT_CANONICAL.fetch_add(1,Ordering::Relaxed)
        ).map_err(|e|invalid(format!("{e}")))?;
        Ok(Self{draft:topo::BrepDraft{
            incarnation,units,tolerance,
            vertices:Vec::new(),edges:Vec::new(),coedges:Vec::new(),
            loops:Vec::new(),faces:Vec::new(),shells:Vec::new(),
            carriers:topo::GeometryCatalog::default(),
        }})
    }
    pub(crate) fn vertex(&mut self,id:VertexId,p:Point3)->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.vertices.len()||!p.is_finite(){
            return Err(invalid("direct B-rep vertex sequence inconsistent"));
        }
        self.draft.vertices.push(topo::BrepVertex{id,point:p});
        Ok(())
    }
    pub(crate) fn edge(
        &mut self,id:EdgeId,start:VertexId,end:VertexId,curve:Curve3,domain:ParameterRange,
    )->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.edges.len() {return Err(invalid("direct edge ID mismatch"));}
        bound(domain)?;
        let id3=topo::Curve3Id(self.draft.carriers.curves3.len() as u32);
        self.draft.carriers.curves3.push(curve);
        self.draft.edges.push(topo::BrepEdge{
            id,vertices:[start,end],curve:id3,domain,kind:topo::EdgeKind::Physical,
            uses:[CoedgeId(u32::MAX);2],
        });
        Ok(())
    }
    pub(crate) fn face(
        &mut self,id:FaceId,surface_id:u32,surface:Surface3,sense:FaceSense,loops:Vec<LoopId>
    )->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.faces.len()||
           surface_id as usize!=self.draft.carriers.surfaces.len(){
            return Err(invalid("direct face or source surface ID mismatch"));
        }
        self.draft.carriers.surfaces.push(surface);
        self.draft.faces.push(topo::BrepFace{
            id,surface:topo::SurfaceId(surface_id),loops,
            sense:match sense{
                FaceSense::Forward=>topo::FaceSense::Forward,
                FaceSense::Reversed=>topo::FaceSense::Reversed,
            },
        });
        Ok(())
    }
    pub(crate) fn loop_(
        &mut self,id:LoopId,face:FaceId,role:LoopRole,first:CoedgeId
    )->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.loops.len(){return Err(invalid("direct loop ID mismatch"));}
        self.draft.loops.push(topo::BrepLoop{id,face,role,first});
        Ok(())
    }
    pub(crate) fn coedge(
        &mut self,id:CoedgeId,edge:EdgeId,face:FaceId,loop_id:LoopId,
        reversed:bool,next:CoedgeId,prev:CoedgeId,curve:Curve2,trim_domain:ParameterRange,
    )->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.coedges.len(){
            return Err(invalid("direct coedge ID mismatch"));
        }
        let physical=self.draft.edges.get(edge.0 as usize)
            .ok_or_else(||invalid("missing physical edge in direct builder"))?;
        let (t0,t1)=bound(physical.domain)?;
        let (u0,u1)=bound(trim_domain)?;
        let physical_curve=*self.draft.carriers.curves3.get(physical.curve.0 as usize)
            .ok_or_else(||invalid("direct edge missing physical curve"))?;
        let support=self.draft.faces.get(face.0 as usize)
            .and_then(|f|self.draft.carriers.surfaces.get(f.surface.0 as usize))
            .copied().ok_or_else(||invalid("missing analytic source-face support"))?;
        let slope=(u1-u0)/(t1-t0);
        let maps=[
            topo::AffineParameterMap{scale:slope,offset:u0-slope*t0},
            topo::AffineParameterMap{scale:-slope,offset:u1+slope*t0},
        ];
        let mut best=None;
        for candidate in maps {
            let mut peak=0.0_f64;
            for fraction in [0.0,0.25,0.5,0.75,1.0] {
                let t=t0+(t1-t0)*fraction;
                let a=physical_curve.evaluate(t)?;
                let uv=curve.evaluate(candidate.map(t).map_err(|e|invalid(format!("{e}")))?)?;
                let b=support.evaluate(uv.u,uv.v)?.point;
                let d=extent(a,b);
                if !d.is_finite(){return Err(invalid("nonfinite direct curve/UV residual"));}
                peak=peak.max(d);
            }
            if best.map_or(true,|(_,previous)|peak<previous){best=Some((candidate,peak));}
        }
        let (map,error)=best.ok_or_else(||invalid("cannot bind source UV curve"))?;
        let a=self.draft.vertices.get(physical.vertices[0].0 as usize)
            .ok_or_else(||invalid("missing first physical vertex"))?.point;
        let b=self.draft.vertices.get(physical.vertices[1].0 as usize)
            .ok_or_else(||invalid("missing second physical vertex"))?.point;
        if error>32.*self.draft.tolerance.length_at(extent(a,b).max(1.)){
            return Err(invalid(format!("physical source curve and UV coedge diverge by {error}")));
        }
        let uv_id=topo::Curve2Id(self.draft.carriers.curves2.len() as u32);
        self.draft.carriers.curves2.push(curve);
        self.draft.coedges.push(topo::BrepCoedge{
            id,edge,face,loop_id,pcurve:uv_id,edge_to_pcurve:map,
            reversed,next,prev,twin:CoedgeId(u32::MAX),
        });
        Ok(())
    }
    pub(crate) fn twins(
        &mut self,edge:EdgeId,left:CoedgeId,right:CoedgeId,
    )->Result<(),GeometryError>{
        let e=self.draft.edges.get_mut(edge.0 as usize)
            .ok_or_else(||invalid("linking missing canonical edge"))?;
        if e.id!=edge{return Err(invalid("incorrect physical edge in twin pairing"));}
        e.uses=[left,right];
        let a=self.draft.coedges.get_mut(left.0 as usize)
            .ok_or_else(||invalid("linking missing first coedge"))?;
        a.twin=right;
        let b=self.draft.coedges.get_mut(right.0 as usize)
            .ok_or_else(||invalid("linking missing second coedge"))?;
        b.twin=left;
        Ok(())
    }
    pub(crate) fn updated_trim(&mut self,id:CoedgeId,curve:Curve2)->Result<(),GeometryError>{
        let item=self.draft.coedges.get(id.0 as usize)
            .ok_or_else(||invalid("missing reparameterized G8 UV coedge"))?;
        let slot=self.draft.carriers.curves2.get_mut(item.pcurve.0 as usize)
            .ok_or_else(||invalid("missing reparameterized G8 UV curve"))?;
        *slot=curve;
        Ok(())
    }
    pub(crate) fn shell(&mut self,id:ShellId,faces:Vec<FaceId>)->Result<(),GeometryError>{
        if id.0 as usize!=self.draft.shells.len(){return Err(invalid("direct shell ID mismatch"));}
        self.draft.shells.push(topo::BrepShell{id,faces,sense:topo::ShellSense::Outward});
        Ok(())
    }
    pub(crate) fn finish(self)->Result<topo::BrepModel,GeometryError>{
        self.draft.check().map_err(|e|invalid(format!("direct G7/G8 canonical B-rep: {e}")))
    }
}
