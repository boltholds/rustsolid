//! Strong-atomic structural editing for independently removable analytic
//! through-holes and ordinary curve/UV reparameterizations in ONE batch.
//!
//! Bounded, exact mathematical scope: a certified axis-aligned box with up
//! to 32 non-overlapping complete +Y cylindrical through-holes. Structural
//! rebuilds create/retire true vertices, edges, coedges, inner cap loops and
//! cylindrical faces; volume, surface and centroid are rederived from source
//! primitives. A common BrepCarrierEditor then applies shape-preserving
//! carrier edits before the candidate and history become observable.
//!
//! Hole handles are process-scoped and generational, distinct from dense
//! topology IDs that may change when an earlier hole is removed. Retired or
//! abandoned staged handles NEVER resurrect through Undo, Redo or slot reuse.
use crate::{
    BrepBody, BrepCarrierEditor, BrepModel, CarrierOwner, CoedgeId, EdgeId,
    FaceId, GeometryError, HistoryDirection, HoleParameters, Solid,
};
use crate::multi_through_hole::{build_multiple,check_hole_set,topology_for,HoleTopology};
use std::sync::atomic::{AtomicU64,Ordering};

static NEXT_EDITOR:AtomicU64=AtomicU64::new(1);
static NEXT_GENERATION:AtomicU64=AtomicU64::new(1);
const MAX_EDITS:usize=256;
const MAX_CARRIER_HISTORY:usize=2048;
const UNDO_LIMIT:usize=128;

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn fresh_generation()->Result<u64,GeometryError>{
    let value=NEXT_GENERATION.fetch_add(1,Ordering::Relaxed);
    if value==0 || value==u64::MAX{return Err(invalid("structural hole generation exhausted"));}
    Ok(value)
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash)]
pub struct StructuralHoleHandle {
    pub editor_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash,serde::Deserialize)]
#[serde(rename_all="snake_case")]
pub enum HoleEdgeRole {BottomCircle,TopCircle,Seam}
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash,serde::Deserialize)]
#[serde(rename_all="snake_case")]
pub enum HoleTrimRole {
    BottomCap,TopCap,BottomSide,TopSide,SeamPositive,SeamNegative,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum StructuralEdgeRef {
    Base(EdgeId),
    Hole{hole:StructuralHoleHandle,role:HoleEdgeRole},
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum StructuralTrimRef {
    Base(CoedgeId),
    Hole{hole:StructuralHoleHandle,role:HoleTrimRole},
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum StructuralFaceRef {
    Base(FaceId),
    HoleSide(StructuralHoleHandle),
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct StructuralHoleTopology {
    pub side_face:FaceId,
    pub bottom_edge:EdgeId,
    pub top_edge:EdgeId,
    pub seam_edge:EdgeId,
    pub bottom_loop:crate::LoopId,
    pub top_loop:crate::LoopId,
    pub side_loop:crate::LoopId,
}
impl From<HoleTopology> for StructuralHoleTopology {
    fn from(t:HoleTopology)->Self{
        Self{side_face:t.side_face,bottom_edge:t.bottom_edge,
            top_edge:t.top_edge,seam_edge:t.seam_edge,
            bottom_loop:t.bottom_loop,top_loop:t.top_loop,side_loop:t.side_loop}
    }
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct StructuralHoleBinding {
    pub handle:StructuralHoleHandle,
    pub params:HoleParameters,
    pub topology:StructuralHoleTopology,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct StructuralRemap {
    pub handle:StructuralHoleHandle,
    pub before:StructuralHoleTopology,
    pub after:StructuralHoleTopology,
}
#[derive(Debug,Clone)]
pub struct StructuralCommit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub command_count:usize,
    pub created:Vec<StructuralHoleBinding>,
    pub retired:Vec<StructuralHoleBinding>,
    pub remapped:Vec<StructuralRemap>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct StructuralTransition {
    pub direction:HistoryDirection,
    pub revision_before:u64,
    pub revision_after:u64,
}
#[derive(Debug,Clone)]
pub struct StructuralSnapshot {
    pub editor_token:u64,
    pub revision:u64,
    pub model:BrepModel,
    pub holes:Vec<StructuralHoleBinding>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum StructuralHandleError {ForeignEditor,StaleHole}
impl std::fmt::Display for StructuralHandleError {
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{
        match self {
            Self::ForeignEditor=>write!(f,"hole handle belongs to another editor"),
            Self::StaleHole=>write!(f,"hole handle refers to a retired incarnation"),
        }
    }
}
impl std::error::Error for StructuralHandleError {}

#[derive(Debug,Clone,Copy)]
struct HoleSlot {
    generation:u64,
    current:Option<HoleParameters>,
}
#[derive(Debug,Clone,Copy)]
enum PrivateEdge {
    Base(EdgeId),
    Hole{slot:u32,role:HoleEdgeRole},
}
#[derive(Debug,Clone,Copy)]
enum PrivateTrim {
    Base(CoedgeId),
    Hole{slot:u32,role:HoleTrimRole},
}
#[derive(Debug,Clone,Copy)]
enum PrivateFace {
    Base(FaceId),
    HoleSide(u32),
}
#[derive(Debug,Clone,Copy)]
enum CarrierIntent {
    Edge{target:PrivateEdge,scale:f64,offset:f64},
    Trim{target:PrivateTrim,scale:f64,offset:f64},
    Face{target:PrivateFace,du:f64,dv:f64},
}
impl CarrierIntent {
    fn references(&self,slot:u32)->bool{
        match self {
            Self::Edge{target:PrivateEdge::Hole{slot:s,..},..} |
            Self::Trim{target:PrivateTrim::Hole{slot:s,..},..} |
            Self::Face{target:PrivateFace::HoleSide(s),..} => *s==slot,
            _=>false,
        }
    }
}
#[derive(Debug,Clone,Copy)]
enum CompiledCarrierIntent {
    Edge{handle:crate::CarrierHandle,scale:f64,offset:f64},
    Trim{handle:crate::CarrierHandle,scale:f64,offset:f64},
    Face{handle:crate::CarrierHandle,du:f64,dv:f64},
}
#[derive(Clone)]
struct State {
    slots:Vec<HoleSlot>,
    intents:Vec<CarrierIntent>,
    model:BrepModel,
}
fn active_holes(slots:&[HoleSlot])->Vec<HoleParameters>{
    slots.iter().filter_map(|s|s.current).collect()
}
fn ordinal(slots:&[HoleSlot],slot:u32)->Result<usize,GeometryError>{
    let index=slot as usize;
    if slots.get(index).and_then(|s|s.current).is_none(){
        return Err(invalid("carrier operation refers to a retired hole"));
    }
    Ok(slots[..index].iter().filter(|s|s.current.is_some()).count())
}
fn edge_owner(slots:&[HoleSlot],target:PrivateEdge)
    ->Result<EdgeId,GeometryError>
{
    Ok(match target{
        PrivateEdge::Base(id) if id.0<12=>id,
        PrivateEdge::Base(_)=>return Err(invalid("box source edge selector outside 0..12")),
        PrivateEdge::Hole{slot,role}=>{
            let t=topology_for(ordinal(slots,slot)?)?;
            match role {
                HoleEdgeRole::BottomCircle=>t.bottom_edge,
                HoleEdgeRole::TopCircle=>t.top_edge,
                HoleEdgeRole::Seam=>t.seam_edge,
            }
        }
    })
}
fn trim_owner(slots:&[HoleSlot],target:PrivateTrim)
    ->Result<CoedgeId,GeometryError>
{
    Ok(match target{
        PrivateTrim::Base(id) if id.0<24=>id,
        PrivateTrim::Base(_)=>return Err(invalid("box source coedge selector outside 0..24")),
        PrivateTrim::Hole{slot,role}=>{
            let t=topology_for(ordinal(slots,slot)?)?;
            match role{
                HoleTrimRole::BottomCap=>t.bottom_trim,
                HoleTrimRole::TopCap=>t.top_trim,
                HoleTrimRole::BottomSide=>t.bottom_side_trim,
                HoleTrimRole::TopSide=>t.top_side_trim,
                HoleTrimRole::SeamPositive=>t.seam_positive,
                HoleTrimRole::SeamNegative=>t.seam_negative,
            }
        }
    })
}
fn face_owner(slots:&[HoleSlot],target:PrivateFace)
    ->Result<FaceId,GeometryError>
{
    Ok(match target{
        PrivateFace::Base(id) if id.0<6=>id,
        PrivateFace::Base(_)=>return Err(invalid("box source face selector outside 0..6")),
        PrivateFace::HoleSide(slot)=>topology_for(ordinal(slots,slot)?)?.side_face,
    })
}
fn prepare(block:&Solid,slots:&[HoleSlot],intents:&[CarrierIntent])
    ->Result<BrepModel,GeometryError>
{
    let holes=active_holes(slots);
    let model=build_multiple(block,&holes)?;
    if intents.is_empty(){return Ok(model);}
    if intents.len()>MAX_CARRIER_HISTORY{
        return Err(invalid("structural carrier replay history exceeds 2048 intents"));
    }
    let mut carrier=BrepCarrierEditor::new(model)?;
    // Compile each intent against the current dense graph *before* taking
    // the exclusive carrier transaction borrow. Handles refer to owners,
    // so catalog compaction between bounded chunks does not change meaning.
    for chunk in intents.chunks(MAX_EDITS) {
        let compiled=chunk.iter().map(|intent|match *intent {
            CarrierIntent::Edge{target,scale,offset}=>{
                let id=edge_owner(slots,target)?;
                let handle=carrier.handle(CarrierOwner::Edge(id))
                    .ok_or_else(||invalid("missing compiled edge owner"))?;
                Ok(CompiledCarrierIntent::Edge{handle,scale,offset})
            }
            CarrierIntent::Trim{target,scale,offset}=>{
                let id=trim_owner(slots,target)?;
                let handle=carrier.handle(CarrierOwner::Coedge(id))
                    .ok_or_else(||invalid("missing compiled UV trim owner"))?;
                Ok(CompiledCarrierIntent::Trim{handle,scale,offset})
            }
            CarrierIntent::Face{target,du,dv}=>{
                let id=face_owner(slots,target)?;
                let handle=carrier.handle(CarrierOwner::Face(id))
                    .ok_or_else(||invalid("missing compiled face owner"))?;
                Ok(CompiledCarrierIntent::Face{handle,du,dv})
            }
        }).collect::<Result<Vec<_>,GeometryError>>()?;
        let mut tx=carrier.begin_transaction();
        for cmd in compiled {
            match cmd {
                CompiledCarrierIntent::Edge{handle,scale,offset}=>{
                    tx.reparameterize_edge(handle,scale,offset)?;
                }
                CompiledCarrierIntent::Trim{handle,scale,offset}=>{
                    tx.reparameterize_trim(handle,scale,offset)?;
                }
                CompiledCarrierIntent::Face{handle,du,dv}=>{
                    tx.rebase_face_uv(handle,du,dv)?;
                }
            }
        }
        tx.commit()?;
    }
    carrier.validate()?;
    Ok(carrier.model().clone())
}
fn bindings(token:u64,slots:&[HoleSlot])->Vec<StructuralHoleBinding>{
    let mut ordinal=0usize;
    let mut result=Vec::new();
    for (i,s) in slots.iter().enumerate(){
        if let Some(params)=s.current{
            let t=topology_for(ordinal).expect("validated bounded hole topology");
            result.push(StructuralHoleBinding{
                handle:StructuralHoleHandle{editor_token:token,slot:i as u32,
                    generation:s.generation},
                params,topology:t.into(),
            });
            ordinal+=1;
        }
    }
    result
}
/// Transactional structural owner. Source block and geometry history belong
/// to this editor; neither the caller's source nor a donor arena is modified.
pub struct BrepStructuralEditor {
    token:u64,
    revision:u64,
    template:Solid,
    state:State,
    undo:Vec<State>,
    redo:Vec<State>,
}
impl BrepStructuralEditor {
    pub fn new(block:Solid)->Result<Self,GeometryError>{
        crate::boolean::certify_box(&block).map_err(|e|invalid(e.to_string()))?;
        let token=NEXT_EDITOR.fetch_add(1,Ordering::Relaxed);
        if token==0||token==u64::MAX{return Err(invalid("structural editor token exhausted"));}
        let model=build_multiple(&block,&[])?;
        Ok(Self{token,revision:0,template:block,
            state:State{slots:Vec::new(),intents:Vec::new(),model},
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn from_through_hole(hole:&crate::ThroughHoleBrep)->Result<Self,GeometryError>{
        hole.validate()?;
        let mut out=Self::new(hole.template_solid().clone())?;
        let generation=fresh_generation()?;
        let params=hole.parameters();
        let slots=vec![HoleSlot{generation,current:Some(params)}];
        let model=prepare(&out.template,&slots,&[])?;
        out.state=State{slots,intents:Vec::new(),model};
        Ok(out)
    }
    pub fn from_body(body:&BrepBody)->Result<Self,GeometryError>{
        match body{
            BrepBody::Polyhedral(s)=>Self::new(s.clone()),
            BrepBody::ThroughHole(h)=>Self::from_through_hole(h),
            BrepBody::Cylinder(_) | BrepBody::BlindHole(_) | BrepBody::IntersectingBlind(_)=>Err(invalid("only box and full through-hole are supported structural sources")),
        }
    }
    pub fn revision(&self)->u64{self.revision}
    pub fn editor_token(&self)->u64{self.token}
    pub fn model(&self)->&BrepModel{&self.state.model}
    /// Geometry ownership stays with this editor; mesh is display-only.
    pub(crate) fn source_box(&self)->&Solid{&self.template}
    /// Independently certify planar cap loops and emit a watertight display mesh.
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::CertifiedDisplayMesh,GeometryError>{
        crate::FaceRegionMesher::tessellate(self,segments)
    }
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn hole_count(&self)->usize{active_holes(&self.state.slots).len()}
    pub fn bindings(&self)->Vec<StructuralHoleBinding>{
        bindings(self.token,&self.state.slots)
    }
    pub fn handle(&self,slot:u32)->Option<StructuralHoleHandle>{
        self.state.slots.get(slot as usize).filter(|s|s.current.is_some())
            .map(|s|StructuralHoleHandle{editor_token:self.token,slot,
                generation:s.generation})
    }
    pub fn resolve(&self,handle:StructuralHoleHandle)
        ->Result<StructuralHoleBinding,StructuralHandleError>
    {
        if handle.editor_token!=self.token{return Err(StructuralHandleError::ForeignEditor);}
        let slot=self.state.slots.get(handle.slot as usize)
            .filter(|s|s.current.is_some() && s.generation==handle.generation)
            .ok_or(StructuralHandleError::StaleHole)?;
        let t=topology_for(ordinal(&self.state.slots,handle.slot)
            .map_err(|_|StructuralHandleError::StaleHole)?)
            .map_err(|_|StructuralHandleError::StaleHole)?;
        Ok(StructuralHoleBinding{handle,params:slot.current.unwrap(),topology:t.into()})
    }
    pub fn snapshot(&self)->StructuralSnapshot{
        StructuralSnapshot{editor_token:self.token,revision:self.revision,
            model:self.state.model.clone(),holes:self.bindings()}
    }
    pub fn check_snapshot(&self,s:&StructuralSnapshot)->Result<(),GeometryError>{
        if s.editor_token!=self.token{return Err(invalid("foreign structural snapshot"));}
        if s.revision!=self.revision{return Err(invalid("stale structural snapshot revision"));}
        Ok(())
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.state.model.validate()?;
        let canonical=prepare(&self.template,&self.state.slots,&self.state.intents)?;
        if canonical.summary()!=self.state.model.summary() ||
           canonical.mass!=self.state.model.mass ||
           canonical.bbox!=self.state.model.bbox ||
           canonical.geometry.curves3!=self.state.model.geometry.curves3 ||
           canonical.geometry.curves2!=self.state.model.geometry.curves2 ||
           canonical.geometry.surfaces!=self.state.model.geometry.surfaces {
            return Err(invalid("structural state diverges from independently rebuilt geometry"));
        }
        Ok(())
    }
    pub fn begin_transaction(&mut self)->BrepStructuralTransaction<'_>{
        let slots=self.state.slots.clone();
        let intents=self.state.intents.clone();
        BrepStructuralTransaction{owner:self,slots,intents,failed:false,
            staged_commands:0}
    }
    pub fn undo(&mut self)->Result<Option<StructuralTransition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<StructuralTransition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)->Result<Option<StructuralTransition>,GeometryError>{
        let next=if redo{&self.redo}else{&self.undo};
        let Some(target)=next.last() else{return Ok(None)};
        target.model.validate()?;
        // Historical generation values are NEVER restored for an element
        // resurrected by Undo or Redo. Surviving feature handles stay stable.
        let mut candidate=target.clone();
        for (i,slot) in candidate.slots.iter_mut().enumerate(){
            let current=self.state.slots.get(i);
            if slot.current.is_some(){
                // Occupying the same sparse slot does not imply semantic
                // survival. Remove+reinsert can replace hole A with B in
                // one commit. Preserve generation only for identical live
                // incarnations; otherwise mint a new generation on Undo/Redo.
                slot.generation=if current.is_some_and(|s|
                    s.current.is_some() && s.generation==slot.generation
                ){
                    current.unwrap().generation
                }else{fresh_generation()?};
            }else if let Some(existing)=current{
                slot.generation=existing.generation;
            }
        }
        let before=self.revision;
        let after=before.checked_add(1).ok_or_else(||invalid("structural revision overflow"))?;
        let prior=std::mem::replace(&mut self.state,candidate);
        if redo{self.redo.pop();self.undo.push(prior);}else{
            self.undo.pop();self.redo.push(prior);
        }
        self.revision=after;
        Ok(Some(StructuralTransition{
            direction:if redo{HistoryDirection::Redo}else{HistoryDirection::Undo},
            revision_before:before,revision_after:after,
        }))
    }
}

/// Staged structural + COW analytic edits. A failed stage poisons this entire
/// transaction, even if its Err is ignored. Nothing is published until the
/// *whole* topology, analytic mass and carrier graph have validated.
pub struct BrepStructuralTransaction<'a>{
    owner:&'a mut BrepStructuralEditor,
    slots:Vec<HoleSlot>,
    intents:Vec<CarrierIntent>,
    failed:bool,
    staged_commands:usize,
}
impl BrepStructuralTransaction<'_> {
    fn guard(&mut self)->Result<(),GeometryError>{
        if self.failed{return Err(invalid("structural transaction already aborted"));}
        if self.staged_commands>=MAX_EDITS{
            self.failed=true;
            return Err(invalid("structural batch exceeds 256 operations"));
        }
        Ok(())
    }
    fn accept(&mut self){
        self.staged_commands+=1;
    }
    fn checked(&mut self,handle:StructuralHoleHandle)
        ->Result<usize,GeometryError>
    {
        if self.failed{return Err(invalid("structural transaction already aborted"));}
        let result=if handle.editor_token!=self.owner.token {
            Err(invalid("hole handle belongs to another editor"))
        }else{
            self.slots.get(handle.slot as usize)
                .filter(|s|s.current.is_some() && s.generation==handle.generation)
                .map(|_|handle.slot as usize)
                .ok_or_else(||invalid("hole handle is stale or retired"))
        };
        if result.is_err(){self.failed=true;}
        result
    }
    fn validate_set(&mut self)->Result<(),GeometryError>{
        let r=check_hole_set(&self.owner.template,&active_holes(&self.slots));
        if r.is_err(){self.failed=true;}
        r
    }
    pub fn handle(&self,slot:u32)->Option<StructuralHoleHandle>{
        if self.failed{return None;}
        self.slots.get(slot as usize).filter(|s|s.current.is_some())
            .map(|s|StructuralHoleHandle{editor_token:self.owner.token,
                slot,generation:s.generation})
    }
    pub fn add_through_hole(&mut self,parameters:HoleParameters)
        ->Result<StructuralHoleHandle,GeometryError>
    {
        self.guard()?;
        let index=self.slots.iter().position(|s|s.current.is_none())
            .unwrap_or(self.slots.len());
        let generation=match fresh_generation(){
            Ok(g)=>g,Err(e)=>{self.failed=true;return Err(e);}
        };
        if index==self.slots.len(){
            self.slots.push(HoleSlot{generation,current:Some(parameters)});
        }else{
            self.slots[index]=HoleSlot{generation,current:Some(parameters)};
        }
        self.validate_set()?;
        self.accept();
        Ok(StructuralHoleHandle{editor_token:self.owner.token,
            slot:index as u32,generation})
    }
    pub fn remove_through_hole(&mut self,handle:StructuralHoleHandle)
        ->Result<(),GeometryError>
    {
        self.guard()?;
        let index=self.checked(handle)?;
        self.slots[index].current=None;
        self.intents.retain(|op|!op.references(handle.slot));
        self.accept();
        Ok(())
    }
    pub fn resize_through_hole(&mut self,handle:StructuralHoleHandle,
        parameters:HoleParameters)->Result<(),GeometryError>
    {
        self.guard()?;
        let index=self.checked(handle)?;
        self.slots[index].current=Some(parameters);
        self.validate_set()?;
        self.accept();
        Ok(())
    }
    fn private_edge(&mut self,r:StructuralEdgeRef)->Result<PrivateEdge,GeometryError>{
        match r {
            StructuralEdgeRef::Base(id) if id.0<12=>Ok(PrivateEdge::Base(id)),
            StructuralEdgeRef::Base(_)=>{self.failed=true;Err(invalid("non-base edge selector"))},
            StructuralEdgeRef::Hole{hole,role}=>{
                self.checked(hole)?;
                Ok(PrivateEdge::Hole{slot:hole.slot,role})
            }
        }
    }
    fn private_trim(&mut self,r:StructuralTrimRef)->Result<PrivateTrim,GeometryError>{
        match r {
            StructuralTrimRef::Base(id) if id.0<24=>Ok(PrivateTrim::Base(id)),
            StructuralTrimRef::Base(_)=>{self.failed=true;Err(invalid("non-base UV coedge selector"))},
            StructuralTrimRef::Hole{hole,role}=>{
                self.checked(hole)?;
                Ok(PrivateTrim::Hole{slot:hole.slot,role})
            }
        }
    }
    fn private_face(&mut self,r:StructuralFaceRef)->Result<PrivateFace,GeometryError>{
        match r {
            StructuralFaceRef::Base(id) if id.0<6=>Ok(PrivateFace::Base(id)),
            StructuralFaceRef::Base(_)=>{self.failed=true;Err(invalid("non-base face selector"))},
            StructuralFaceRef::HoleSide(hole)=>{
                self.checked(hole)?;Ok(PrivateFace::HoleSide(hole.slot))
            }
        }
    }
    fn check_params(&mut self,values:&[f64])->Result<(),GeometryError>{
        if values.iter().any(|v|!v.is_finite()){
            self.failed=true;return Err(invalid("nonfinite structural carrier edit"));
        }
        Ok(())
    }
    pub fn reparameterize_edge(&mut self,target:StructuralEdgeRef,
        scale:f64,offset:f64)->Result<(),GeometryError>
    {
        self.guard()?;self.check_params(&[scale,offset])?;
        if scale<=0.0{self.failed=true;return Err(invalid("invalid edge parameter scale"));}
        let target=self.private_edge(target)?;
        self.intents.push(CarrierIntent::Edge{target,scale,offset});
        self.accept();
        Ok(())
    }
    pub fn reparameterize_trim(&mut self,target:StructuralTrimRef,
        scale:f64,offset:f64)->Result<(),GeometryError>
    {
        self.guard()?;self.check_params(&[scale,offset])?;
        if scale<=0.0{self.failed=true;return Err(invalid("invalid trim parameter scale"));}
        let target=self.private_trim(target)?;
        self.intents.push(CarrierIntent::Trim{target,scale,offset});
        self.accept();
        Ok(())
    }
    pub fn rebase_face_uv(&mut self,target:StructuralFaceRef,
        du:f64,dv:f64)->Result<(),GeometryError>
    {
        self.guard()?;self.check_params(&[du,dv])?;
        let target=self.private_face(target)?;
        self.intents.push(CarrierIntent::Face{target,du,dv});
        self.accept();
        Ok(())
    }
    pub fn staged_len(&self)->usize{self.staged_commands}
    pub fn commit(self)->Result<StructuralCommit,GeometryError>{
        if self.failed {return Err(invalid("cannot commit poisoned structural transaction"));}
        if self.staged_commands==0{return Err(invalid("empty structural transaction"));}
        if self.intents.len()>MAX_CARRIER_HISTORY{
            return Err(invalid("structural carrier replay intent budget exceeded"));
        }
        let model=prepare(&self.owner.template,&self.slots,&self.intents)?;
        let before=self.owner.revision;
        let after=before.checked_add(1).ok_or_else(||invalid("structural revision overflow"))?;
        let old=bindings(self.owner.token,&self.owner.state.slots);
        let current=bindings(self.owner.token,&self.slots);
        let mut created=Vec::new();
        let mut retired=Vec::new();
        let mut remapped=Vec::new();
        for now in &current {
            match old.iter().find(|b|b.handle.slot==now.handle.slot) {
                Some(prior) if prior.handle==now.handle=>{
                    if prior.topology!=now.topology{
                        remapped.push(StructuralRemap{handle:now.handle,
                            before:prior.topology,after:now.topology});
                    }
                }
                _=>created.push(*now),
            }
        }
        for prior in &old {
            if !current.iter().any(|b|b.handle==prior.handle) {
                retired.push(*prior);
            }
        }
        let next=State{slots:self.slots,intents:self.intents,model};
        let previous=std::mem::replace(&mut self.owner.state,next);
        if self.owner.undo.len()>=UNDO_LIMIT{self.owner.undo.remove(0);}
        self.owner.undo.push(previous);
        self.owner.redo.clear();
        self.owner.revision=after;
        Ok(StructuralCommit{revision_before:before,revision_after:after,
            command_count:self.staged_commands,created,retired,remapped})
    }
}
