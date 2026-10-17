//! G5a: in-memory, revision-atomic four-shell Boolean document.
//!
//! Stages only source-validated G4b2b3b/G4c candidates for a strict AABB
//! minus full-height upright cylinder. Result identity uses named semantic
//! corner/role paths, never an index assumed persistent outside this model.
//! Scoped generational handles are invalidated on Undo, Redo and regeneration.
//!
//! All fallible validation, identity allocation and history capacity checks
//! precede publication. No partial observable document state on recoverable
//! errors. This is NOT a disk/database transaction, and general G5 curved
//! B-rep editing is still unsupported.

use crate::{
    AnalyticFacePatchPartitionPlan,AnalyticFacePatchPartitionResult,
    AnalyticFacePatchHold,IndependentVolumeOraclePlan,
    IndependentVolumeOracleResult,IndependentVolumeHold,
    BooleanOperation,BrepBody,BrepModel,GeometryError,GeometryTolerance,
    HistoryDirection,TopologyEntity,TopologyName,VertexId,EdgeId,CoedgeId,
    FaceId,ShellId,LoopId,GraphOperand,
};
use std::collections::{BTreeMap,BTreeSet};
use std::sync::atomic::{AtomicU64,Ordering};

static NEXT_G5_INCAR:AtomicU64=AtomicU64::new(1);

fn invalid(s:impl Into<String>)->GeometryError {
    GeometryError::InvalidEdit(s.into())
}
fn fresh()->Result<u64,GeometryError>{
    let id=NEXT_G5_INCAR.fetch_add(1,Ordering::Relaxed);
    if id==0||id==u64::MAX {
        return Err(invalid("G5 identity generation exhausted"));
    }
    Ok(id)
}
fn operand_stamp(body:&BrepBody)->String{
    // Exact private, in-process source comparison. NOT a serialized or
    // cryptographic file digest. Includes full cylinder geometry even though
    // that source family has no independent monotonically changing revision.
    format!("{body:?}")
}
fn valid_feature(key:&str)->bool{
    !key.is_empty()&&key.len()<=96&&key.bytes()
        .all(|b|b.is_ascii_alphanumeric()||b==b'-'||b==b'_')
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G5Hold{
    ExistingCommittedResult,
    UnsupportedCandidate(AnalyticFacePatchHold),
    UncertifiedVolume(IndependentVolumeHold),
}
#[derive(Debug,Clone)]
pub enum G5PrepareResult {
    Ready(G5PreparedBoolean),Held(G5Hold),
}
/// An offline proposal. Even when fields are altered by the caller, no
/// commit can bypass source reconstruction or independent G4c checks.
#[derive(Debug,Clone)]
pub struct G5PreparedBoolean {
    pub candidate:AnalyticFacePatchPartitionPlan,
    pub volume:IndependentVolumeOraclePlan,
    document_token:u64,
    expected_revision:u64,
    feature_key:String,
    input_stamps:[String;2],
}
impl G5PreparedBoolean{
    pub fn expected_revision(&self)->u64{self.expected_revision}
    pub fn document_token(&self)->u64{self.document_token}
    pub fn candidate(&self)->&AnalyticFacePatchPartitionPlan{&self.candidate}
    pub fn volume(&self)->&IndependentVolumeOraclePlan{&self.volume}
}
#[derive(Debug,Clone)]
struct AcceptedBoolean {
    candidate:AnalyticFacePatchPartitionPlan,
    volume:IndependentVolumeOraclePlan,
}
#[derive(Debug,Clone)]
struct DocumentState {
    inputs:[BrepBody;2],
    accepted:Option<AcceptedBoolean>,
}
impl DocumentState {
    fn validate(&self,tolerance:GeometryTolerance)->Result<(),GeometryError>{
        self.inputs[0].validate()?;
        self.inputs[1].validate()?;
        if let Some(done)=&self.accepted{
            done.candidate.validate_against(&self.inputs[0],&self.inputs[1],tolerance)?;
            done.volume.validate_against(&self.inputs[0],&self.inputs[1],
                &done.candidate,tolerance)?;
            if !done.volume.volume_certified||done.volume.can_apply
                ||!done.candidate.closed_shell_proven||done.candidate.can_apply {
                return Err(invalid("G5 state includes an unverified application"));
            }
        }
        Ok(())
    }
    fn model(&self)->Option<&BrepModel>{
        self.accepted.as_ref().map(|accepted|accepted.candidate.brep())
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash)]
pub struct G5TopologyHandle {
    pub document_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug,Clone,PartialEq,Eq)]
pub enum G5TopologyHandleError { ForeignDocument,StaleHandle }
#[derive(Debug,Clone,PartialEq,Eq)]
pub struct G5TopologyBinding {
    pub handle:G5TopologyHandle,
    /// Stable logical feature name, valid across model regeneration, Undo and
    /// Redo; unlike a handle, it can be saved in the owning document.
    pub name:TopologyName,
    /// Result-local numeric index; this MAY change after regeneration.
    pub entity:TopologyEntity,
}
#[derive(Debug,Clone)]
struct NamedRegistry {
    token:u64,
    by_name:BTreeMap<TopologyName,G5TopologyBinding>,
    by_entity:BTreeMap<TopologyEntity,TopologyName>,
    by_slot:Vec<G5TopologyBinding>,
}
impl NamedRegistry {
    fn empty(token:u64)->Self{
        Self{token,by_name:BTreeMap::new(),by_entity:BTreeMap::new(),
            by_slot:Vec::new()}
    }
    // Pure derivation. Validating a document must never consume process-wide
    // generation IDs or create unused handles.
    fn semantic_names(
        candidate:Option<&AnalyticFacePatchPartitionPlan>,
        feature_key:&str,
    )->Result<BTreeMap<TopologyEntity,TopologyName>,GeometryError>{
        let Some(candidate)=candidate else{return Ok(BTreeMap::new())};
        let m=candidate.brep();
        m.validate()?;
        if !candidate.closed_shell_proven||candidate.can_apply
            ||candidate.corners.len()!=4||m.shells.len()!=4 {
            return Err(invalid("G5 names require four closed source-certified shells"));
        }
        let mut semantic=BTreeMap::<TopologyEntity,TopologyName>::new();
        let mut unique=BTreeSet::<TopologyName>::new();
        let mut add=|entity:TopologyEntity,name:String|->Result<(),GeometryError>{
            let key=TopologyName(name);
            if semantic.contains_key(&entity)||!unique.insert(key.clone()){
                return Err(invalid("G5 duplicate semantic name or material entity"));
            }
            semantic.insert(entity,key);
            Ok(())
        };
        add(TopologyEntity::Body,format!("{feature_key}/difference/body"))?;
        const VERTEX_ROLES:[&str;6]=[
            "bottom_corner","bottom_x_arc","bottom_z_arc",
            "top_corner","top_x_arc","top_z_arc",
        ];
        const EDGE_ROLES:[&str;9]=[
            "bottom_x","bottom_z","bottom_circle",
            "top_x","top_z","top_circle",
            "corner_generator","x_generator","z_generator",
        ];
        const FACE_ROLES:[&str;5]=[
            "bottom","top","x_wall","z_wall","cylinder_wall",
        ];
        let mut labels=BTreeSet::new();
        for corner in &candidate.corners{
            if ![-1_i8,1_i8].contains(&corner.x_side)
                ||![-1_i8,1_i8].contains(&corner.z_side)
                ||!labels.insert((corner.x_side,corner.z_side)){
                return Err(invalid("G5 invalid semantic corner labels"));
            }
            let sx=if corner.x_side<0 {"neg"}else{"pos"};
            let sz=if corner.z_side<0 {"neg"}else{"pos"};
            let prefix=format!("{feature_key}/difference/corner/x_{sx}/z_{sz}");
            let sh=m.shells.get(corner.shell.0 as usize)
                .filter(|s|s.id==corner.shell)
                .ok_or_else(||invalid("G5 invalid shell material identity"))?;
            if sh.faces.len()!=5||sh.faces!=corner.faces{
                return Err(invalid("G5 material shell ownership changed"));
            }
            add(TopologyEntity::Shell(corner.shell),format!("{prefix}/shell"))?;
            let edge_for=|index:usize|->Result<&crate::Edge,GeometryError>{
                let id=corner.edges[index];
                m.edges.get(id.0 as usize).filter(|e|e.id==id)
                    .ok_or_else(||invalid("G5 corner edge ID invalid"))
            };
            let va=[edge_for(0)?.start,edge_for(0)?.end,edge_for(1)?.end,
                    edge_for(3)?.start,edge_for(3)?.end,edge_for(4)?.end];
            if edge_for(0)?.start!=edge_for(1)?.start
                ||edge_for(3)?.start!=edge_for(4)?.start {
                return Err(invalid("G5 corner material vertices do not have semantic roles"));
            }
            for (index,&id) in va.iter().enumerate(){
                if m.vertices.get(id.0 as usize).is_none(){
                    return Err(invalid("G5 missing source-certified corner vertex"));
                }
                add(TopologyEntity::Vertex(id),format!(
                    "{prefix}/vertex/{}",VERTEX_ROLES[index]))?;
            }
            for (index,&fid) in corner.faces.iter().enumerate(){
                let face=m.faces.get(fid.0 as usize)
                    .filter(|f|f.id==fid&&f.loops.len()==1)
                    .ok_or_else(||invalid("G5 missing face with one analytic loop"))?;
                let lid=face.loops[0];
                if m.loops.get(lid.0 as usize)
                    .is_none_or(|lp|lp.id!=lid||lp.face!=fid) {
                    return Err(invalid("G5 face loop owner changed"));
                }
                add(TopologyEntity::Face(fid),format!(
                    "{prefix}/face/{}",FACE_ROLES[index]))?;
                add(TopologyEntity::Loop(lid),format!(
                    "{prefix}/loop/{}",FACE_ROLES[index]))?;
            }
            for (index,&eid) in corner.edges.iter().enumerate(){
                let e=edge_for(index)?;
                add(TopologyEntity::Edge(eid),format!(
                    "{prefix}/edge/{}",EDGE_ROLES[index]))?;
                for &coid in &e.coedges{
                    let co=m.coedges.get(coid.0 as usize)
                        .filter(|co|co.id==coid&&co.edge==eid)
                        .ok_or_else(||invalid("G5 missing incident coedge"))?;
                    let face_role=corner.faces.iter()
                        .position(|id|*id==co.face)
                        .ok_or_else(||invalid("G5 coedge owner outside shell"))?;
                    add(TopologyEntity::Coedge(coid),format!(
                        "{prefix}/edge/{}/use/{}",
                        EDGE_ROLES[index],FACE_ROLES[face_role]
                    ))?;
                }
            }
        }
        if labels.len()!=4 {
            return Err(invalid("G5 result omitted one or more material corners"));
        }
        let expected=1+m.vertices.len()+m.edges.len()+m.coedges.len()
            +m.loops.len()+m.faces.len()+m.shells.len();
        if semantic.len()!=expected {
            return Err(invalid("G5 source topology has unnamed entities"));
        }
        Ok(semantic)
    }
    fn from_candidate(
        candidate:Option<&AnalyticFacePatchPartitionPlan>,
        feature_key:&str,token:u64,
    )->Result<Self,GeometryError>{
        let semantic=Self::semantic_names(candidate,feature_key)?;
        let Some(candidate)=candidate else{return Ok(Self::empty(token))};
        let m=candidate.brep();
        let mut out=Self::empty(token);
        out.by_slot.try_reserve(semantic.len())
            .map_err(|_|invalid("G5 unable to allocate identity registry"))?;
        for (entity,name) in semantic {
            let slot=u32::try_from(out.by_slot.len())
                .map_err(|_|invalid("G5 material identity capacity exceeded"))?;
            let binding=G5TopologyBinding {
                handle:G5TopologyHandle{
                    document_token:token,slot,generation:fresh()?,
                },
                name:name.clone(),entity,
            };
            out.by_entity.insert(entity,name.clone());
            out.by_name.insert(name,binding.clone());
            out.by_slot.push(binding);
        }
        out.validate(m)?;
        Ok(out)
    }
    fn validate(&self,model:&BrepModel)->Result<(),GeometryError>{
        if self.by_name.len()!=self.by_entity.len()
            ||self.by_name.len()!=self.by_slot.len(){
            return Err(invalid("G5 identity maps have inconsistent cardinality"));
        }
        for (i,binding) in self.by_slot.iter().enumerate(){
            if binding.handle.document_token!=self.token
                ||binding.handle.slot as usize!=i
                ||binding.handle.generation==0
                ||self.by_name.get(&binding.name)!=Some(binding)
                ||self.by_entity.get(&binding.entity)!=Some(&binding.name){
                return Err(invalid("G5 semantic identity registry inconsistency"));
            }
        }
        let size=1+model.vertices.len()+model.edges.len()+model.coedges.len()
            +model.loops.len()+model.faces.len()+model.shells.len();
        if self.by_slot.len()!=size {
            return Err(invalid("G5 registry does not cover result topology"));
        }
        Ok(())
    }
    fn resolve(&self,handle:G5TopologyHandle)
        ->Result<G5TopologyBinding,G5TopologyHandleError>{
        if handle.document_token!=self.token {
            return Err(G5TopologyHandleError::ForeignDocument);
        }
        let Some(current)=self.by_slot.get(handle.slot as usize)else{
            return Err(G5TopologyHandleError::StaleHandle);
        };
        if current.handle!=handle{
            return Err(G5TopologyHandleError::StaleHandle);
        }
        Ok(current.clone())
    }
    fn by_name(&self,name:&TopologyName)->Option<G5TopologyBinding>{
        self.by_name.get(name).cloned()
    }
}

#[derive(Debug,Clone)]
pub struct G5DocumentSnapshot {
    pub document_token:u64,
    pub revision:u64,
    pub feature_key:String,
    pub has_result:bool,
    pub result_volume:Option<f64>,
    pub names:Vec<G5TopologyBinding>,
    pub undo_depth:usize,
    pub redo_depth:usize,
}
#[derive(Debug,Clone)]
pub struct G5Commit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub committed_volume:f64,
    pub shell_components:usize,
    pub bound_entities:usize,
    pub changed_names:Vec<G5TopologyBinding>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct G5Transition {
    pub direction:HistoryDirection,
    pub revision_before:u64,
    pub revision_after:u64,
    pub has_result:bool,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G5DocumentOperand{First,Second}
impl G5DocumentOperand {
    fn index(self)->usize{match self{Self::First=>0,Self::Second=>1}}
}

/// A single in-process owner; &mut self acts as an exclusive commit lock.
/// The feature_key is caller-owned document identity; TopologyName is stable
/// across undo/redo. Handles remain process-local and generation-scoped.
pub struct G5BooleanDocument {
    token:u64,
    feature_key:String,
    revision:u64,
    tolerance:GeometryTolerance,
    state:DocumentState,
    registry:NamedRegistry,
    undo:Vec<DocumentState>,
    redo:Vec<DocumentState>,
}
impl G5BooleanDocument {
    pub fn new(
        feature_key:&str,left:BrepBody,right:BrepBody,
        tolerance:GeometryTolerance,
    )->Result<Self,GeometryError>{
        if !valid_feature(feature_key) {
            return Err(invalid("G5 document feature key must be ASCII alphanumeric, - or _"));
        }
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let token=fresh()?;
        Ok(Self{
            token,feature_key:feature_key.to_owned(),revision:0,tolerance,
            state:DocumentState{inputs:[left,right],accepted:None},
            registry:NamedRegistry::empty(token),
            undo:Vec::new(),redo:Vec::new(),
        })
    }
    pub fn document_token(&self)->u64{self.token}
    pub fn revision(&self)->u64{self.revision}
    pub fn feature_key(&self)->&str{&self.feature_key}
    pub fn model(&self)->Option<&BrepModel>{self.state.model()}
    pub fn source(&self,operand:G5DocumentOperand)->&BrepBody{
        &self.state.inputs[operand.index()]
    }
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn identity_count(&self)->usize{self.registry.by_slot.len()}
    pub fn bindings(&self)->Vec<G5TopologyBinding>{self.registry.by_slot.clone()}
    pub fn resolve_name(&self,name:&TopologyName)->Option<G5TopologyBinding>{
        self.registry.by_name(name)
    }
    pub fn resolve_handle(&self,handle:G5TopologyHandle)
        ->Result<G5TopologyBinding,G5TopologyHandleError>{
        self.registry.resolve(handle)
    }
    pub fn snapshot(&self)->G5DocumentSnapshot{
        G5DocumentSnapshot{
            document_token:self.token,revision:self.revision,
            feature_key:self.feature_key.clone(),
            has_result:self.model().is_some(),
            result_volume:self.model().map(|m|m.mass.volume),
            names:self.bindings(),
            undo_depth:self.undo_depth(),redo_depth:self.redo_depth(),
        }
    }
    pub fn check_snapshot(&self,snapshot:&G5DocumentSnapshot)
        ->Result<(),GeometryError>{
        if snapshot.document_token!=self.token {
            return Err(invalid("foreign G5 document snapshot"));
        }
        if snapshot.revision!=self.revision {
            return Err(invalid("stale G5 document snapshot revision"));
        }
        // Snapshot is diagnostic, not authority for mutation.
        Ok(())
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.tolerance.validate()?;
        self.state.validate(self.tolerance)?;
        let expected=NamedRegistry::semantic_names(
            self.state.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key,
        )?;
        if self.registry.by_entity!=expected {
            return Err(invalid("G5 live semantic names differ from the proven B-rep"));
        }
        if let Some(m)=self.model(){
            self.registry.validate(m)?;
        }else if !self.registry.by_name.is_empty()
            ||!self.registry.by_entity.is_empty()
            ||!self.registry.by_slot.is_empty(){
            return Err(invalid("G5 result handles still live without a model"));
        }
        Ok(())
    }
    pub fn prepare(&self)->Result<G5PrepareResult,GeometryError>{
        if self.state.accepted.is_some(){
            return Ok(G5PrepareResult::Held(G5Hold::ExistingCommittedResult));
        }
        self.state.inputs[0].validate()?;
        self.state.inputs[1].validate()?;
        let candidate=match AnalyticFacePatchPartitionPlan::build(
            &self.state.inputs[0],&self.state.inputs[1],
            BooleanOperation::Difference,self.tolerance)?
        {
            AnalyticFacePatchPartitionResult::Ready(v)=>v,
            AnalyticFacePatchPartitionResult::Held(why)=>
                return Ok(G5PrepareResult::Held(G5Hold::UnsupportedCandidate(why))),
        };
        let volume=match IndependentVolumeOraclePlan::build(
            &self.state.inputs[0],&self.state.inputs[1],&candidate,self.tolerance)?
        {
            IndependentVolumeOracleResult::Ready(v)=>v,
            IndependentVolumeOracleResult::Held(why)=>
                return Ok(G5PrepareResult::Held(G5Hold::UncertifiedVolume(why))),
        };
        if !candidate.closed_shell_proven||candidate.can_apply
            ||!volume.volume_certified||volume.can_apply {
            return Err(invalid("G5 preparation received an invalid G4c volume gate"));
        }
        Ok(G5PrepareResult::Ready(G5PreparedBoolean {
            candidate,volume,document_token:self.token,
            expected_revision:self.revision,
            feature_key:self.feature_key.clone(),
            input_stamps:[
                operand_stamp(&self.state.inputs[0]),
                operand_stamp(&self.state.inputs[1]),
            ],
        }))
    }
    pub fn commit(&mut self,proposal:G5PreparedBoolean)
        ->Result<G5Commit,GeometryError>{
        if proposal.document_token!=self.token
            ||proposal.feature_key!=self.feature_key {
            return Err(invalid("G5 proposal belongs to a foreign document"));
        }
        if proposal.expected_revision!=self.revision
            ||self.state.accepted.is_some(){
            return Err(invalid("G5 stale proposed revision or Boolean already committed"));
        }
        let latest=[
            operand_stamp(&self.state.inputs[0]),
            operand_stamp(&self.state.inputs[1]),
        ];
        if proposal.input_stamps!=latest{
            return Err(invalid("G5 original source geometry changed since preparation"));
        }
        // This check is immediately before the only owner-side publication.
        // Re-derive all G1/G2/G3/G4/G4c evidence from the CURRENT sources.
        proposal.candidate.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],self.tolerance)?;
        proposal.volume.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],
            &proposal.candidate,self.tolerance)?;
        if !proposal.volume.volume_certified||proposal.volume.can_apply
            ||!proposal.candidate.closed_shell_proven||proposal.candidate.can_apply{
            return Err(invalid("G5 refuses any unverified analytical Boolean"));
        }
        let result=AcceptedBoolean{
            candidate:proposal.candidate,volume:proposal.volume,
        };
        let volume=result.candidate.brep().mass.volume;
        let shell_components=result.candidate.brep().shells.len();
        let names=NamedRegistry::from_candidate(
            Some(&result.candidate),&self.feature_key,self.token)?;
        let new_state=DocumentState {
            inputs:self.state.inputs.clone(),accepted:Some(result),
        };
        new_state.validate(self.tolerance)?;
        let after=self.revision.checked_add(1)
            .ok_or_else(||invalid("G5 document revision overflow"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G5 cannot reserve atomic undo checkpoint"))?;
        let report=G5Commit {
            revision_before:self.revision,revision_after:after,
            committed_volume:volume,shell_components,
            bound_entities:names.by_slot.len(),
            changed_names:names.by_slot.clone(),
        };
        // COMMIT BARRIER: no fallible operations after this point.
        let old=std::mem::replace(&mut self.state,new_state);
        self.undo.push(old);
        self.redo.clear();
        self.registry=names;
        self.revision=after;
        Ok(report)
    }
    /// Any in-document source edit invalidates the active Boolean result.
    /// This is a separate atomic operation with its own undo checkpoint.
    /// Even replacing with geometrically identical input increments revision,
    /// preventing ABA-style reuse of a proposal prepared earlier.
    pub fn replace_source(
        &mut self,operand:G5DocumentOperand,new_source:BrepBody,
    )->Result<u64,GeometryError>{
        new_source.validate()?;
        let mut proposed=self.state.clone();
        proposed.inputs[operand.index()]=new_source;
        proposed.accepted=None;
        proposed.validate(self.tolerance)?;
        let after=self.revision.checked_add(1)
            .ok_or_else(||invalid("G5 document revision overflow"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G5 unable to reserve source edit checkpoint"))?;
        let old=std::mem::replace(&mut self.state,proposed);
        self.undo.push(old);
        self.redo.clear();
        self.registry=NamedRegistry::empty(self.token);
        self.revision=after;
        Ok(after)
    }
    pub fn undo(&mut self)->Result<Option<G5Transition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<G5Transition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)
        ->Result<Option<G5Transition>,GeometryError>{
        let stack=if redo{&self.redo}else{&self.undo};
        let Some(target)=stack.last()else{return Ok(None)};
        target.validate(self.tolerance)?;
        let registry=NamedRegistry::from_candidate(
            target.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key,self.token)?;
        let next=self.revision.checked_add(1)
            .ok_or_else(||invalid("G5 history revision overflow"))?;
        // Reserve the opposite stack BEFORE changing anything.
        if redo {
            self.undo.try_reserve(1)
                .map_err(|_|invalid("G5 Undo history capacity exhausted"))?;
        }else{
            self.redo.try_reserve(1)
                .map_err(|_|invalid("G5 Redo history capacity exhausted"))?;
        }
        let direction=if redo{HistoryDirection::Redo}else{HistoryDirection::Undo};
        let before=self.revision;
        // COMMIT BARRIER: undo/redo source, result and registry publish in
        // one exclusive-owner update, no validation after mutation.
        let restored=if redo{
            self.redo.pop().unwrap()
        }else{
            self.undo.pop().unwrap()
        };
        let previous=std::mem::replace(&mut self.state,restored);
        if redo{self.undo.push(previous);}else{self.redo.push(previous);}
        self.registry=registry;
        self.revision=next;
        Ok(Some(G5Transition{
            direction,revision_before:before,revision_after:next,
            has_result:self.state.accepted.is_some(),
        }))
    }
}
