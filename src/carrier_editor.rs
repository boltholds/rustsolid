//! Strong-atomic, source-validated edits of canonical B-rep geometric carriers.
//!
//! Carrier handles are keyed to their *owning topological edge/coedge/face*,
//! never to the dense shared curve catalog position. Replacement uses
//! copy-on-write; commit removes unreferenced carriers and remaps dense IDs.
//! Supported operations reparameterize existing geometry without changing the
//! physical solid or its analytically computed mass. Generalized topology
//! creation, Booleans, and shape-deforming edits remain separate milestones.
use crate::{
    BrepBody, BrepModel, CoedgeId, Curve2, Curve2Id, Curve3, Curve3Id,
    EdgeId, FaceId, Frame3, GeometryError, HistoryDirection, Line2, Line3,
    ParameterRange, PlaneSurface, Point2Param, Point3, Surface3, SurfaceId,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_CARRIER_EDITOR: AtomicU64 = AtomicU64::new(1);
static NEXT_CARRIER_GENERATION: AtomicU64 = AtomicU64::new(1);
const MAX_EDITS: usize = 256;
const HISTORY_LIMIT: usize = 128;
const SAMPLES: [f64; 9] = [0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0];

fn invalid(message:impl Into<String>)->GeometryError {
    GeometryError::InvalidEdit(message.into())
}
fn length(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn distance(a:Point3,b:Point3)->f64{length(a.sub(b))}
fn dot(a:Point3,b:Point3)->f64{a.dot(b)}
fn bounded(range:ParameterRange)->Result<(f64,f64),GeometryError>{
    match range {
        ParameterRange::Bounded{start,end} if start.is_finite() &&
            end.is_finite() && start<end => Ok((start,end)),
        _=>Err(invalid("carrier requires a finite increasing parameter domain")),
    }
}
fn lerp(range:ParameterRange,fraction:f64)->Result<f64,GeometryError>{
    let (start,end)=bounded(range)?;
    Ok(start+(end-start)*fraction)
}
fn transform_range(domain:ParameterRange,scale:f64,offset:f64)
    ->Result<ParameterRange,GeometryError>
{
    if !scale.is_finite() || scale<=0.0 || !offset.is_finite(){
        return Err(invalid("affine parameter scale must be positive and finite"));
    }
    let (a,b)=bounded(domain)?;
    ParameterRange::bounded((a-offset)/scale,(b-offset)/scale)
}
fn shift_range(range:ParameterRange,offset:f64)
    ->Result<ParameterRange,GeometryError>
{
    match range {
        ParameterRange::Unbounded=>Ok(range),
        ParameterRange::Bounded{start,end}=>ParameterRange::bounded(
            start-offset,end-offset),
    }
}
fn carrier_budget(model:&BrepModel)->f64 {
    let d=model.bbox.max.sub(model.bbox.min);
    let extent=length(d);
    let coord=[model.bbox.min,model.bbox.max].into_iter()
        .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
        .fold(1.0_f64,f64::max);
    model.tolerance.length_at(extent).max(32.0*f64::EPSILON*coord)
}
fn same_support(old:Surface3,new:Surface3,tol:f64,angular:f64)
    ->Result<(),GeometryError>
{
    match (old,new) {
        (Surface3::Plane(a),Surface3::Plane(b))=>{
            let n=a.frame.normal;
            if length(n.sub(b.frame.normal))>angular ||
                (b.frame.origin.sub(a.frame.origin)).dot(n).abs()>tol {
                return Err(invalid("edited plane changes its physical support or orientation"));
            }
        }
        (Surface3::Cylinder(a),Surface3::Cylinder(b))=>{
            let axis=a.frame.normal;
            let displacement=b.frame.origin.sub(a.frame.origin);
            let off_axis=displacement.sub(axis.scale(displacement.dot(axis)));
            if length(axis.sub(b.frame.normal))>angular ||
                length(off_axis)>tol || (a.radius-b.radius).abs()>tol {
                return Err(invalid("edited cylinder changes its radius, axis or orientation"));
            }
        }
        (Surface3::Sphere(a),Surface3::Sphere(b))=>{
            if a!=b {return Err(invalid("sphere parameter edits are not yet supported"));}
        }
        _=>return Err(invalid("changing the surface type would change physical topology")),
    }
    Ok(())
}

/// Scope of the carrier's semantic identity. The underlying catalog ID may
/// move after compaction; this topology-owner key remains the same.
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash)]
pub enum CarrierOwner {
    Edge(EdgeId),
    Coedge(CoedgeId),
    Face(FaceId),
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash)]
pub struct CarrierHandle {
    pub editor_token:u64,
    pub owner:CarrierOwner,
    pub generation:u64,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum CarrierId {
    Curve3(Curve3Id),
    Curve2(Curve2Id),
    Surface(SurfaceId),
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct CarrierBinding {
    pub handle:CarrierHandle,
    pub id:CarrierId,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum CarrierHandleError {ForeignEditor,StaleCarrier}
impl std::fmt::Display for CarrierHandleError {
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{
        match self {
            Self::ForeignEditor=>write!(f,"carrier handle belongs to another editor"),
            Self::StaleCarrier=>write!(f,"carrier handle has an invalid incarnation"),
        }
    }
}
impl std::error::Error for CarrierHandleError {}

#[derive(Debug,Clone)]
pub struct CarrierSnapshot {
    pub editor_token:u64,
    pub revision:u64,
    pub model:BrepModel,
}
#[derive(Debug,Clone)]
pub struct CarrierCommit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub changed:Vec<CarrierBinding>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct CarrierTransition {
    pub direction:HistoryDirection,
    pub revision_before:u64,
    pub revision_after:u64,
}

/// An independently owned, mutable canonical B-rep geometry graph. A caller
/// can start from any already validated polyhedral or analytic body, without
/// giving ownership of its original topology or parametric feature history.
pub struct BrepCarrierEditor {
    token:u64,
    revision:u64,
    baseline:BrepModel,
    current:BrepModel,
    generations:BTreeMap<CarrierOwner,u64>,
    undo:Vec<BrepModel>,
    redo:Vec<BrepModel>,
}
impl BrepCarrierEditor {
    pub fn from_body(body:&BrepBody)->Result<Self,GeometryError>{
        Self::new(body.shared()?)
    }
    pub fn new(model:BrepModel)->Result<Self,GeometryError>{
        model.validate()?;
        let token=NEXT_CARRIER_EDITOR.fetch_add(1,Ordering::Relaxed);
        if token==0 || token==u64::MAX {return Err(invalid("carrier editor token exhausted"));}
        let mut generations=BTreeMap::new();
        let owners=model.edges.iter().map(|e|CarrierOwner::Edge(e.id))
            .chain(model.coedges.iter().map(|c|CarrierOwner::Coedge(c.id)))
            .chain(model.faces.iter().map(|f|CarrierOwner::Face(f.id)));
        for owner in owners {
            let generation=NEXT_CARRIER_GENERATION.fetch_add(1,Ordering::Relaxed);
            if generation==0 || generation==u64::MAX {
                return Err(invalid("carrier generation counter exhausted"));
            }
            generations.insert(owner,generation);
        }
        Ok(Self{token,revision:0,baseline:model.clone(),current:model,
            generations,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn editor_token(&self)->u64{self.token}
    pub fn revision(&self)->u64{self.revision}
    pub fn model(&self)->&BrepModel{&self.current}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn handle(&self,owner:CarrierOwner)->Option<CarrierHandle>{
        self.generations.get(&owner).map(|&generation|
            CarrierHandle{editor_token:self.token,owner,generation})
    }
    pub fn resolve(&self,handle:CarrierHandle)->Result<CarrierBinding,CarrierHandleError>{
        if handle.editor_token!=self.token {
            return Err(CarrierHandleError::ForeignEditor);
        }
        if self.generations.get(&handle.owner)!=Some(&handle.generation){
            return Err(CarrierHandleError::StaleCarrier);
        }
        let id=match handle.owner {
            CarrierOwner::Edge(e)=>self.current.geometry.edges.get(e.0 as usize)
                .map(|geom|CarrierId::Curve3(geom.curve)),
            CarrierOwner::Coedge(c)=>self.current.geometry.coedges.get(c.0 as usize)
                .map(|geom|CarrierId::Curve2(geom.pcurve)),
            CarrierOwner::Face(f)=>self.current.faces.get(f.0 as usize)
                .map(|face|CarrierId::Surface(face.geometry.surface)),
        }.ok_or(CarrierHandleError::StaleCarrier)?;
        Ok(CarrierBinding{handle,id})
    }
    pub fn snapshot(&self)->CarrierSnapshot{
        CarrierSnapshot{editor_token:self.token,revision:self.revision,
            model:self.current.clone()}
    }
    pub fn check_snapshot(&self,snapshot:&CarrierSnapshot)->Result<(),GeometryError>{
        if snapshot.editor_token!=self.token{return Err(invalid("foreign carrier editor snapshot"));}
        if snapshot.revision!=self.revision{return Err(invalid("stale carrier editor snapshot"));}
        Ok(())
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.current.validate()?;
        certify_equivalent(&self.baseline,&self.current)
    }
    pub fn begin_transaction(&mut self)->BrepCarrierTransaction<'_>{
        let candidate=self.current.clone();
        BrepCarrierTransaction{editor:self,staged:candidate,failed:false,
            staged_count:0,changed:BTreeSet::new()}
    }
    pub fn undo(&mut self)->Result<Option<CarrierTransition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<CarrierTransition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)->Result<Option<CarrierTransition>,GeometryError>{
        let source=if redo {&self.redo}else{&self.undo};
        let Some(candidate)=source.last() else{return Ok(None)};
        candidate.validate()?;
        certify_equivalent(&self.baseline,candidate)?;
        let before=self.revision;
        let after=before.checked_add(1).ok_or_else(||invalid("carrier revision overflow"))?;
        let previous=if redo{self.redo.pop().unwrap()}else{self.undo.pop().unwrap()};
        let original=std::mem::replace(&mut self.current,previous);
        if redo {self.undo.push(original);}else{self.redo.push(original);}
        self.revision=after;
        Ok(Some(CarrierTransition{direction:if redo{HistoryDirection::Redo}
            else{HistoryDirection::Undo},revision_before:before,revision_after:after}))
    }
}

/// A poisoned stage can never be committed, even if its Err was ignored.
/// Every stage changes only the isolated canonical candidate.
pub struct BrepCarrierTransaction<'a>{
    editor:&'a mut BrepCarrierEditor,
    staged:BrepModel,
    failed:bool,
    staged_count:usize,
    changed:BTreeSet<CarrierOwner>,
}
impl BrepCarrierTransaction<'_> {
    fn verify(&mut self,handle:CarrierHandle,expected:&'static str)
        ->Result<usize,GeometryError>
    {
        if self.failed{return Err(invalid("carrier transaction already failed"));}
        let result=self.editor.resolve(handle).map_err(|e|invalid(e.to_string()))
            .and_then(|binding|{
                let (actual,index)=match binding.handle.owner{
                    CarrierOwner::Edge(id)=>("edge",id.0 as usize),
                    CarrierOwner::Coedge(id)=>("coedge",id.0 as usize),
                    CarrierOwner::Face(id)=>("face",id.0 as usize),
                };
                if actual!=expected {return Err(invalid(format!(
                    "expected {expected} carrier, found {actual}")));}
                Ok(index)
            });
        if result.is_err(){self.failed=true;}
        result
    }
    fn stage<F>(&mut self,owner:CarrierOwner,change:F)->Result<(),GeometryError>
    where F:FnOnce(&mut BrepModel)->Result<(),GeometryError> {
        if self.failed{return Err(invalid("carrier transaction already failed"));}
        if self.staged_count>=MAX_EDITS {
            self.failed=true;return Err(invalid("carrier batch exceeds edit budget"));
        }
        let result=change(&mut self.staged);
        if result.is_err(){self.failed=true;}else{
            self.staged_count+=1;
            self.changed.insert(owner);
        }
        result
    }

    /// Replace an edge's owned 3D carrier and bounded parameter interval.
    /// Copy-on-write keeps other edges untouched if they shared a catalog id.
    pub fn replace_edge(&mut self,handle:CarrierHandle,curve:Curve3,
        domain:ParameterRange)->Result<(),GeometryError>
    {
        let index=self.verify(handle,"edge")?;
        self.stage(handle.owner,move |model|{
            bounded(domain)?;
            let old=model.edge_curve(EdgeId(index as u32))
                .ok_or_else(||invalid("missing original edge carrier"))?;
            if !matches!((old,curve),(Curve3::Line(_),Curve3::Line(_)) |
                (Curve3::Circle(_),Curve3::Circle(_))) {
                return Err(invalid("changing the edge carrier family is not supported"));
            }
            curve.evaluate(lerp(domain,0.5)?)?;
            let id=Curve3Id(u32::try_from(model.geometry.curves3.len())
                .map_err(|_|invalid("3D curve catalog exhausted"))?);
            model.geometry.curves3.push(curve);
            model.geometry.edges[index].curve=id;
            model.geometry.edges[index].domain=domain;
            Ok(())
        })
    }
    /// Replace one coedge's 2D trim and parameter interval. The opposite
    /// coedge has its OWN UV binding, including on periodic seams.
    pub fn replace_trim(&mut self,handle:CarrierHandle,curve:Curve2,
        domain:ParameterRange)->Result<(),GeometryError>
    {
        let index=self.verify(handle,"coedge")?;
        self.stage(handle.owner,move |model|{
            bounded(domain)?;
            let old=model.coedge_pcurve(CoedgeId(index as u32))
                .ok_or_else(||invalid("missing coedge trim"))?;
            if !matches!((old,curve),(Curve2::Line(_),Curve2::Line(_)) |
                (Curve2::Circle(_),Curve2::Circle(_))) {
                return Err(invalid("changing the coedge trim family is not supported"));
            }
            curve.evaluate(lerp(domain,0.5)?)?;
            let id=Curve2Id(u32::try_from(model.geometry.curves2.len())
                .map_err(|_|invalid("2D curve catalog exhausted"))?);
            model.geometry.curves2.push(curve);
            model.geometry.coedges[index].pcurve=id;
            model.geometry.coedges[index].domain=domain;
            Ok(())
        })
    }
    /// Replace an existing analytic support surface. A complete candidate
    /// must still represent the SAME physical plane or cylinder at commit.
    pub fn replace_surface(&mut self,handle:CarrierHandle,surface:Surface3)
        ->Result<(),GeometryError>
    {
        let index=self.verify(handle,"face")?;
        self.stage(handle.owner,move |model|{
            let original=model.face_surface(FaceId(index as u32))
                .ok_or_else(||invalid("missing original face surface"))?;
            if std::mem::discriminant(&original)!=std::mem::discriminant(&surface) {
                return Err(invalid("changing the supporting surface family is unsupported"));
            }
            surface.evaluate(0.0,0.0)?;
            let id=SurfaceId(u32::try_from(model.geometry.surfaces.len())
                .map_err(|_|invalid("surface catalog exhausted"))?);
            model.geometry.surfaces.push(surface);
            model.faces[index].geometry.surface=id;
            Ok(())
        })
    }

    /// Reparameterize a 3D edge's line by affine t=scale*u+offset; for
    /// circles, only an angular phase offset with scale=1 is supported.
    pub fn reparameterize_edge(&mut self,handle:CarrierHandle,
        scale:f64,offset:f64)->Result<(),GeometryError>
    {
        let index=self.verify(handle,"edge")?;
        if !scale.is_finite() || scale<=0.0 || !offset.is_finite() {
            self.failed=true;return Err(invalid("invalid edge parameter mapping"));
        }
        let geom=self.staged.geometry.edges[index];
        let curve=self.staged.edge_curve(EdgeId(index as u32))
            .ok_or_else(||invalid("missing edge carrier"))?;
        let new=match curve {
            Curve3::Line(line)=>Curve3::Line(Line3{
                origin:line.origin.add(line.direction.scale(offset)),
                direction:line.direction.scale(scale),
            }),
            Curve3::Circle(circle)=>{
                if scale!=1.0 {
                    self.failed=true;
                    return Err(invalid("circle reparameterization only supports positive phase shift"));
                }
                let (sin,cos)=offset.sin_cos();
                let x=circle.frame.x.scale(cos).add(circle.frame.y.scale(sin));
                let y=circle.frame.y.scale(cos).sub(circle.frame.x.scale(sin));
                Curve3::Circle(crate::Circle3{
                    frame:Frame3{x,y,..circle.frame},radius:circle.radius})
            }
        };
        let domain=match transform_range(geom.domain,scale,offset){
            Ok(value)=>value,Err(e)=>{self.failed=true;return Err(e);}
        };
        self.replace_edge(handle,new,domain)
    }
    /// For 2D lines apply the same affine remapping. Native Circle2 does not
    /// store a phase; it may be shifted only by full 2*pi revolutions.
    pub fn reparameterize_trim(&mut self,handle:CarrierHandle,
        scale:f64,offset:f64)->Result<(),GeometryError>
    {
        let index=self.verify(handle,"coedge")?;
        if !scale.is_finite() || scale<=0.0 || !offset.is_finite() {
            self.failed=true;return Err(invalid("invalid trim parameter mapping"));
        }
        let geom=self.staged.geometry.coedges[index];
        let curve=self.staged.coedge_pcurve(CoedgeId(index as u32))
            .ok_or_else(||invalid("missing trim carrier"))?;
        let new=match curve {
            Curve2::Line(line)=>Curve2::Line(Line2 {
                origin:Point2Param{u:line.origin.u+line.delta.u*offset,
                    v:line.origin.v+line.delta.v*offset},
                delta:Point2Param{u:line.delta.u*scale,v:line.delta.v*scale},
            }),
            Curve2::Circle(circle)=>{
                if scale!=1.0 || (offset-(offset/TAU).round()*TAU).abs()>1e-12 {
                    self.failed=true;
                    return Err(invalid("Circle2 accepts only full-turn trim parameter shifts"));
                }
                Curve2::Circle(circle)
            }
        };
        let domain=match transform_range(geom.domain,scale,offset){
            Ok(value)=>value,Err(e)=>{self.failed=true;return Err(e);}
        };
        self.replace_trim(handle,new,domain)
    }

    /// Move the origin of one face's UV *chart*, not its physical support.
    /// Rebase all owned coedge trims and face bounds together, as a single
    /// atomic stage. The cylinder periodic u chart currently stays fixed.
    pub fn rebase_face_uv(&mut self,handle:CarrierHandle,
        du:f64,dv:f64)->Result<(),GeometryError>
    {
        let index=self.verify(handle,"face")?;
        self.stage(handle.owner,move |model|{
            if !du.is_finite() || !dv.is_finite(){
                return Err(invalid("UV offset must be finite"));
            }
            let old=model.face_surface(FaceId(index as u32))
                .ok_or_else(||invalid("missing face surface"))?;
            let (surface,delta)=match old {
                Surface3::Plane(p)=>{
                    let new_origin=p.frame.position(du,dv)?;
                    (Surface3::Plane(PlaneSurface{
                        frame:Frame3{origin:new_origin,..p.frame}}),
                    Point2Param{u:du,v:dv})
                }
                Surface3::Cylinder(c)=>{
                    if du!=0.0 {return Err(invalid(
                        "periodic cylinder angular UV rebasing requires a seam chart editor"));}
                    let new_origin=c.frame.origin.add(c.frame.normal.scale(dv));
                    (Surface3::Cylinder(crate::CylinderSurface{
                        frame:Frame3{origin:new_origin,..c.frame},radius:c.radius}),
                     Point2Param{u:0.0,v:dv})
                }
                Surface3::Sphere(_)=>return Err(invalid(
                    "sphere face UV rebasing is not supported")),
            };
            let sid=SurfaceId(u32::try_from(model.geometry.surfaces.len())
                .map_err(|_|invalid("surface catalog exhausted"))?);
            model.geometry.surfaces.push(surface);
            model.faces[index].geometry.surface=sid;
            let bounds=&mut model.faces[index].geometry.bounds;
            bounds.u=shift_range(bounds.u,delta.u)?;
            bounds.v=shift_range(bounds.v,delta.v)?;
            let coedge_indices=model.coedges.iter().enumerate()
                .filter(|(_,co)|co.face.0 as usize==index)
                .map(|(i,_)|i).collect::<Vec<_>>();
            for coedge_index in coedge_indices {
                let old_id=model.geometry.coedges[coedge_index].pcurve;
                let old=*model.geometry.curves2.get(old_id.0 as usize)
                    .ok_or_else(||invalid("missing owned UV trim"))?;
                let updated=match old {
                    Curve2::Line(l)=>Curve2::Line(Line2{
                        origin:Point2Param{u:l.origin.u-delta.u,v:l.origin.v-delta.v},
                        delta:l.delta,
                    }),
                    Curve2::Circle(c)=>Curve2::Circle(crate::Circle2{
                        center:Point2Param{u:c.center.u-delta.u,
                            v:c.center.v-delta.v},radius:c.radius,
                    }),
                };
                let new_id=Curve2Id(u32::try_from(model.geometry.curves2.len())
                    .map_err(|_|invalid("UV curve catalog exhausted"))?);
                model.geometry.curves2.push(updated);
                model.geometry.coedges[coedge_index].pcurve=new_id;
            }
            Ok(())
        })?;
        // The rebasing also changed all the face's coedge trim bindings;
        // track those changes in the commit summary.
        for coedge in self.staged.coedges.iter().filter(|c|
            c.face.0 as usize==index) {
            self.changed.insert(CarrierOwner::Coedge(coedge.id));
        }
        Ok(())
    }
    pub fn staged_len(&self)->usize{self.staged_count}
    pub fn commit(mut self)->Result<CarrierCommit,GeometryError>{
        if self.failed {return Err(invalid("cannot commit poisoned carrier transaction"));}
        if self.changed.is_empty(){return Err(invalid("empty carrier transaction"));}
        compact_referenced(&mut self.staged)?;
        self.staged.validate()?;
        // Shape and mass must remain invariant. In a later feature-changing
        // editor this gate must be replaced by full independent mass proofs.
        certify_equivalent(&self.editor.current,&self.staged)?;
        certify_equivalent(&self.editor.baseline,&self.staged)?;
        let before=self.editor.revision;
        let after=before.checked_add(1).ok_or_else(||invalid("carrier revision overflow"))?;
        let changes=self.changed.into_iter().map(|owner|{
            let handle=self.editor.handle(owner).expect("validated carrier owner");
            let id=match owner{
                CarrierOwner::Edge(e)=>CarrierId::Curve3(self.staged.geometry.edges[e.0 as usize].curve),
                CarrierOwner::Coedge(c)=>CarrierId::Curve2(self.staged.geometry.coedges[c.0 as usize].pcurve),
                CarrierOwner::Face(f)=>CarrierId::Surface(self.staged.faces[f.0 as usize].geometry.surface),
            };
            CarrierBinding{handle,id}
        }).collect();
        let original=std::mem::replace(&mut self.editor.current,self.staged);
        if self.editor.undo.len()>=HISTORY_LIMIT {self.editor.undo.remove(0);}
        self.editor.undo.push(original);
        self.editor.redo.clear();
        self.editor.revision=after;
        Ok(CarrierCommit{revision_before:before,revision_after:after,changed:changes})
    }
}

fn compact_referenced(model:&mut BrepModel)->Result<(),GeometryError>{
    let old3=std::mem::take(&mut model.geometry.curves3);
    let mut seen3=BTreeMap::<u32,Curve3Id>::new();
    let mut new3=Vec::new();
    for geom in &mut model.geometry.edges {
        let original=geom.curve.0;
        geom.curve=if let Some(&id)=seen3.get(&original){id}else{
            let curve=*old3.get(original as usize)
                .ok_or_else(||invalid("edge references missing curve3"))?;
            let id=Curve3Id(u32::try_from(new3.len())
                .map_err(|_|invalid("curve3 compaction overflow"))?);
            new3.push(curve);seen3.insert(original,id);id
        };
    }
    model.geometry.curves3=new3;
    let old2=std::mem::take(&mut model.geometry.curves2);
    let mut seen2=BTreeMap::<u32,Curve2Id>::new();
    let mut new2=Vec::new();
    for geom in &mut model.geometry.coedges {
        let original=geom.pcurve.0;
        geom.pcurve=if let Some(&id)=seen2.get(&original){id}else{
            let curve=*old2.get(original as usize)
                .ok_or_else(||invalid("coedge references missing curve2"))?;
            let id=Curve2Id(u32::try_from(new2.len())
                .map_err(|_|invalid("curve2 compaction overflow"))?);
            new2.push(curve);seen2.insert(original,id);id
        };
    }
    model.geometry.curves2=new2;
    let old_s=std::mem::take(&mut model.geometry.surfaces);
    let mut seen_s=BTreeMap::<u32,SurfaceId>::new();
    let mut new_s=Vec::new();
    for face in &mut model.faces {
        let original=face.geometry.surface.0;
        face.geometry.surface=if let Some(&id)=seen_s.get(&original){id}else{
            let surface=*old_s.get(original as usize)
                .ok_or_else(||invalid("face references missing surface"))?;
            let id=SurfaceId(u32::try_from(new_s.len())
                .map_err(|_|invalid("surface compaction overflow"))?);
            new_s.push(surface);seen_s.insert(original,id);id
        };
    }
    model.geometry.surfaces=new_s;
    Ok(())
}
fn certify_equivalent(before:&BrepModel,after:&BrepModel)->Result<(),GeometryError>{
    if before.summary()!=after.summary() ||
        before.mass!=after.mass || before.bbox!=after.bbox ||
        before.tolerance!=after.tolerance || before.units!=after.units ||
        before.vertices.len()!=after.vertices.len() ||
        before.vertices.iter().zip(&after.vertices).any(|(a,b)|
            a.id!=b.id || a.position!=b.position) {
        return Err(invalid("carrier edit changed canonical topology, bounds or mass"));
    }
    let eps=carrier_budget(before);
    let angular=before.tolerance.angular;
    for i in 0..before.faces.len(){
        let f=FaceId(i as u32);
        same_support(before.face_surface(f).ok_or_else(||invalid("missing old support"))?,
            after.face_surface(f).ok_or_else(||invalid("missing new support"))?,
            eps,angular)?;
    }
    for i in 0..before.edges.len(){
        let a=before.geometry.edges[i];
        let b=after.geometry.edges[i];
        let old=before.edge_curve(EdgeId(i as u32)).ok_or_else(||invalid("old edge absent"))?;
        let new=after.edge_curve(EdgeId(i as u32)).ok_or_else(||invalid("new edge absent"))?;
        for f in SAMPLES {
            let p=old.evaluate(lerp(a.domain,f)?)?;
            let q=new.evaluate(lerp(b.domain,f)?)?;
            if distance(p,q)>eps {
                return Err(invalid(format!("edge {i} carrier changed physical geometry")));
            }
        }
    }
    for i in 0..before.coedges.len(){
        let co=&before.coedges[i];
        let a=before.geometry.coedges[i];
        let b=after.geometry.coedges[i];
        let left=before.coedge_pcurve(co.id).ok_or_else(||invalid("missing old pcurve"))?;
        let right=after.coedge_pcurve(co.id).ok_or_else(||invalid("missing new pcurve"))?;
        let face=co.face;
        let old_surface=before.face_surface(face).ok_or_else(||invalid("old face absent"))?;
        let new_surface=after.face_surface(face).ok_or_else(||invalid("new face absent"))?;
        for f in SAMPLES {
            let uv_left=left.evaluate(lerp(a.domain,f)?)?;
            let uv_right=right.evaluate(lerp(b.domain,f)?)?;
            let p=old_surface.evaluate(uv_left.u,uv_left.v)?.point;
            let q=new_surface.evaluate(uv_right.u,uv_right.v)?.point;
            if distance(p,q)>eps {
                return Err(invalid(format!("coedge {i} trim changed physical geometry")));
            }
        }
    }
    Ok(())
}
