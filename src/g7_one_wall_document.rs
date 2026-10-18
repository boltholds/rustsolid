//! G7c: an exclusively-owned source-certified one-wall Boolean document.
//!
//! This is an intentionally bounded G5 transaction: it commits only a
//! G7a closed-shell candidate with a fresh G7b source/flux volume proof.
//! Generic G7a/G7b plans retain can_apply=false; the document owns the only
//! mutation barrier. G5 v1/four-corner and G6 v2/two-strip APIs are unchanged.

use std::collections::{BTreeMap,BTreeSet};
use std::sync::atomic::{AtomicU64,Ordering};
use crate::{
    AnalyticOneWallPlan,AnalyticOneWallResult,AnalyticOneWallHold,
    AnalyticOneWallVolumeOracle,AnalyticOneWallVolumeResult,
    OneWallVolumeHold,BooleanOperation,BrepBody,BrepModel,
    GeometryError,GeometryTolerance,HistoryDirection,TopologyName,
    TopologyEntity,VertexId,
};

static NEXT_G7_ID:AtomicU64=AtomicU64::new(1);
fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
fn fresh()->Result<u64,GeometryError>{
    let id=NEXT_G7_ID.fetch_add(1,Ordering::Relaxed);
    if id==0||id==u64::MAX {
        return Err(invalid("G7 process-local identity capacity exhausted"));
    }
    Ok(id)
}
fn valid_key(key:&str)->bool {
    !key.is_empty() && key.len()<=96 && key.bytes().all(|x|
        x.is_ascii_alphanumeric()||x==b'_'||x==b'-')
}
fn stamp(source:&BrepBody)->String {
    // Private in-memory source stamp, NOT a cryptographic/durable digest.
    format!("{source:?}")
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G7Hold {
    UnsupportedGeometry(AnalyticOneWallHold),
    UncertifiedVolume(OneWallVolumeHold),
    AlreadyCommitted,
}
#[derive(Debug,Clone)]
pub enum G7PrepareResult {
    Ready(G7PreparedBoolean),
    Held(G7Hold),
}
#[derive(Debug,Clone)]
pub struct G7PreparedBoolean {
    pub candidate:AnalyticOneWallPlan,
    pub volume:AnalyticOneWallVolumeOracle,
    token:u64,
    expected_revision:u64,
    feature_key:String,
    source_stamps:[String;2],
}
impl G7PreparedBoolean{
    pub fn expected_revision(&self)->u64{self.expected_revision}
}

#[derive(Debug,Clone)]
pub(crate) struct G7Accepted {
    pub candidate:AnalyticOneWallPlan,
    pub volume:AnalyticOneWallVolumeOracle,
}
#[derive(Debug,Clone)]
pub(crate) struct G7State {
    pub inputs:[BrepBody;2],
    pub accepted:Option<G7Accepted>,
}
impl G7State{
    pub(crate) fn model(&self)->Option<&BrepModel>{
        self.accepted.as_ref().map(|x|x.candidate.brep())
    }
    pub(crate) fn validate(&self,tolerance:GeometryTolerance)
        ->Result<(),GeometryError>{
        self.inputs[0].validate()?;
        self.inputs[1].validate()?;
        if let Some(ref accepted)=self.accepted{
            accepted.candidate.validate_against(
                &self.inputs[0],&self.inputs[1],tolerance)?;
            accepted.volume.validate_against(
                &self.inputs[0],&self.inputs[1],
                &accepted.candidate,tolerance)?;
            if !accepted.candidate.closed_shell_proven
                ||!accepted.volume.volume_certified
                ||accepted.candidate.can_apply||accepted.volume.can_apply{
                return Err(invalid("G7 requires separately audited source/flux certificates"));
            }
        }
        Ok(())
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash)]
pub struct G7Handle {
    pub document_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G7HandleError {ForeignDocument,StaleHandle}
#[derive(Debug,Clone,PartialEq,Eq)]
pub struct G7TopologyBinding {
    pub handle:G7Handle,
    pub name:TopologyName,
    /// This numeric ID is local to one B-rep incarnation.
    pub entity:TopologyEntity,
}
#[derive(Debug,Clone)]
pub(crate) struct G7IdentityRegistry{
    token:u64,
    pub by_entity:BTreeMap<TopologyEntity,TopologyName>,
    by_name:BTreeMap<TopologyName,G7TopologyBinding>,
    entries:Vec<G7TopologyBinding>,
}
impl G7IdentityRegistry{
    pub(crate) fn empty(token:u64)->Self{
        Self{token,by_entity:BTreeMap::new(),
            by_name:BTreeMap::new(),entries:Vec::new()}
    }
    /// Pure logical name derivation. No state mutation or handle allocation.
    pub(crate) fn names(
        candidate:Option<&AnalyticOneWallPlan>,feature_key:&str,
    )->Result<BTreeMap<TopologyEntity,TopologyName>,GeometryError>{
        let Some(plan)=candidate else{return Ok(BTreeMap::new())};
        let m=plan.brep();
        m.validate()?;
        if m.shells.len()!=1||plan.can_apply||!plan.closed_shell_proven
            ||plan.source_faces.len()!=8||plan.component.faces.len()!=8{
            return Err(invalid("G7 semantic identities require a complete one-wall shell"));
        }
        let mut names=BTreeMap::<TopologyEntity,TopologyName>::new();
        let mut unique=BTreeSet::<TopologyName>::new();
        let mut add=|entity:TopologyEntity,key:String|->Result<(),GeometryError>{
            let name=TopologyName(key);
            if names.insert(entity,name.clone()).is_some()
                ||!unique.insert(name){
                return Err(invalid("G7 duplicate logical topology role"));
            }
            Ok(())
        };
        let prefix=format!("{feature_key}/difference/notch");
        add(TopologyEntity::Body,format!("{prefix}/body"))?;
        add(TopologyEntity::Shell(plan.component.shell),
            format!("{prefix}/shell"))?;
        const V_ROLES:[&str;6]=[
            "corner_zmin_xmin","corner_zmin_xmax","corner_zmax_xmax",
            "corner_zmax_xmin","wall_upper_cut","wall_lower_cut",
        ];
        const E_ROLES:[&str;18]=[
            "bottom_zmin","bottom_xmax","bottom_zmax","bottom_xmin_upper",
            "bottom_major_circle","bottom_xmin_lower",
            "top_zmin","top_xmax","top_zmax","top_xmin_upper",
            "top_major_circle","top_xmin_lower",
            "corner_zmin_xmin_generator","corner_zmin_xmax_generator",
            "corner_zmax_xmax_generator","corner_zmax_xmin_generator",
            "upper_cut_generator","lower_cut_generator",
        ];
        const F_ROLES:[&str;8]=[
            "bottom","top","outer_zmin","outer_xmax","outer_zmax",
            "original_xmin_upper","original_xmin_lower","inner_cylinder",
        ];
        // Six physical generators pair each bottom vertex with exactly
        // one upper vertex, independent of any numeric vertex ordering.
        for i in 0..6{
            let edge=plan.component.edges[12+i];
            let generator=m.edges.get(edge.0 as usize)
                .filter(|e|e.id==edge)
                .ok_or_else(||invalid("G7 missing original vertical generator"))?;
            for (vertex,level) in [(generator.start,"bottom"),(generator.end,"top")]{
                let actual=m.vertices.get(vertex.0 as usize)
                    .filter(|v|v.id==vertex)
                    .ok_or_else(||invalid("G7 missing semantic vertex"))?;
                if !actual.position.is_finite(){
                    return Err(invalid("G7 vertex has invalid coordinates"));
                }
                add(TopologyEntity::Vertex(vertex),
                    format!("{prefix}/vertex/{level}/{}",V_ROLES[i]))?;
            }
        }
        for (i,&fid) in plan.component.faces.iter().enumerate(){
            let face=m.faces.get(fid.0 as usize)
                .filter(|f|f.id==fid&&f.loops.len()==1)
                .ok_or_else(||invalid("G7 face lacks its trimmed physical loop"))?;
            let loop_id=face.loops[0];
            if m.loops.get(loop_id.0 as usize).is_none_or(|x|x.face!=fid){
                return Err(invalid("G7 face loop ownership changed"));
            }
            add(TopologyEntity::Face(fid),format!("{prefix}/face/{}",F_ROLES[i]))?;
            add(TopologyEntity::Loop(loop_id),format!("{prefix}/loop/{}",F_ROLES[i]))?;
        }
        for (i,&eid) in plan.component.edges.iter().enumerate(){
            let edge=m.edges.get(eid.0 as usize)
                .filter(|x|x.id==eid)
                .ok_or_else(||invalid("G7 physical edge missing"))?;
            add(TopologyEntity::Edge(eid),format!("{prefix}/edge/{}",E_ROLES[i]))?;
            for &coid in &edge.coedges{
                let co=m.coedges.get(coid.0 as usize)
                    .filter(|x|x.id==coid&&x.edge==eid)
                    .ok_or_else(||invalid("G7 physical coedge missing"))?;
                let role=plan.component.faces.iter().position(|x|*x==co.face)
                    .ok_or_else(||invalid("G7 coedge belongs outside the result shell"))?;
                add(TopologyEntity::Coedge(coid),format!(
                    "{prefix}/edge/{}/use/{}",E_ROLES[i],F_ROLES[role]))?;
            }
        }
        let expected=1+m.shells.len()+m.vertices.len()+m.edges.len()
            +m.coedges.len()+m.faces.len()+m.loops.len();
        if names.len()!=expected || expected!=84{
            return Err(invalid("G7 logical names fail to cover 84 entities"));
        }
        Ok(names)
    }
    pub(crate) fn from_candidate(
        candidate:Option<&AnalyticOneWallPlan>,key:&str,token:u64,
    )->Result<Self,GeometryError>{
        let names=Self::names(candidate,key)?;
        let mut out=Self::empty(token);
        out.entries.try_reserve(names.len())
            .map_err(|_|invalid("G7 identity allocation failed"))?;
        for (entity,name) in names{
            let slot=u32::try_from(out.entries.len())
                .map_err(|_|invalid("G7 logical entity count overflow"))?;
            let binding=G7TopologyBinding {
                handle:G7Handle{
                    document_token:token,slot,generation:fresh()?,
                },name:name.clone(),entity,
            };
            out.by_entity.insert(entity,name.clone());
            out.by_name.insert(name,binding.clone());
            out.entries.push(binding);
        }
        Ok(out)
    }
    fn resolve(&self,h:G7Handle)->Result<G7TopologyBinding,G7HandleError>{
        if h.document_token!=self.token{
            return Err(G7HandleError::ForeignDocument);
        }
        self.entries.get(h.slot as usize)
            .filter(|b|b.handle==h)
            .cloned()
            .ok_or(G7HandleError::StaleHandle)
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G7Operand {First,Second}
impl G7Operand{
    fn index(self)->usize{match self{Self::First=>0,Self::Second=>1}}
}
#[derive(Debug,Clone)]
pub struct G7Commit {
    pub revision_before:u64,
    pub revision_after:u64,
    pub result_volume:f64,
    pub shell_components:usize,
    pub topology_names:Vec<G7TopologyBinding>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct G7Transition {
    pub direction:HistoryDirection,
    pub revision_before:u64,
    pub revision_after:u64,
    pub has_result:bool,
}
pub struct G7BooleanDocument{
    pub(crate) token:u64,
    pub(crate) feature_key:String,
    pub(crate) revision:u64,
    pub(crate) tolerance:GeometryTolerance,
    pub(crate) state:G7State,
    pub(crate) registry:G7IdentityRegistry,
    pub(crate) undo:Vec<G7State>,
    pub(crate) redo:Vec<G7State>,
}
impl G7BooleanDocument{
    pub fn new(
        feature_key:&str,first:BrepBody,second:BrepBody,
        tolerance:GeometryTolerance,
    )->Result<Self,GeometryError>{
        if !valid_key(feature_key){
            return Err(invalid("G7 document feature key is invalid"));
        }
        tolerance.validate()?;
        first.validate()?;second.validate()?;
        let token=fresh()?;
        Ok(Self{token,feature_key:feature_key.into(),revision:0,tolerance,
            state:G7State{inputs:[first,second],accepted:None},
            registry:G7IdentityRegistry::empty(token),
            undo:Vec::new(),redo:Vec::new(),
        })
    }
    pub fn feature_key(&self)->&str{&self.feature_key}
    pub fn revision(&self)->u64{self.revision}
    pub fn document_token(&self)->u64{self.token}
    pub fn model(&self)->Option<&BrepModel>{self.state.model()}
    pub fn source(&self,operand:G7Operand)->&BrepBody{
        &self.state.inputs[operand.index()]
    }
    pub fn bindings(&self)->Vec<G7TopologyBinding>{self.registry.entries.clone()}
    pub fn resolve_name(&self,name:&TopologyName)->Option<G7TopologyBinding>{
        self.registry.by_name.get(name).cloned()
    }
    pub fn resolve_handle(&self,h:G7Handle)->Result<G7TopologyBinding,G7HandleError>{
        self.registry.resolve(h)
    }
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn validate(&self)->Result<(),GeometryError>{
        self.state.validate(self.tolerance)?;
        let names=G7IdentityRegistry::names(
            self.state.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key)?;
        if names!=self.registry.by_entity
            ||self.registry.entries.len()!=names.len()
            ||self.registry.by_name.len()!=names.len(){
            return Err(invalid("G7 identity registry does not cover current B-rep"));
        }
        for (i,b) in self.registry.entries.iter().enumerate(){
            if b.handle.document_token!=self.token
                ||b.handle.slot as usize!=i
                ||b.handle.generation==0
                ||self.registry.by_name.get(&b.name)!=Some(b){
                return Err(invalid("G7 scoped identity no longer matches its handle"));
            }
        }
        Ok(())
    }
    pub fn prepare(&self)->Result<G7PrepareResult,GeometryError>{
        if self.state.accepted.is_some(){
            return Ok(G7PrepareResult::Held(G7Hold::AlreadyCommitted));
        }
        let candidate=match AnalyticOneWallPlan::build(
            &self.state.inputs[0],&self.state.inputs[1],
            BooleanOperation::Difference,self.tolerance)?
        {
            AnalyticOneWallResult::Ready(x)=>x,
            AnalyticOneWallResult::Held(why)=>
                return Ok(G7PrepareResult::Held(G7Hold::UnsupportedGeometry(why))),
        };
        let volume=match AnalyticOneWallVolumeOracle::build(
            &self.state.inputs[0],&self.state.inputs[1],&candidate,self.tolerance)?
        {
            AnalyticOneWallVolumeResult::Ready(x)=>x,
            AnalyticOneWallVolumeResult::Held(why)=>
                return Ok(G7PrepareResult::Held(G7Hold::UncertifiedVolume(why))),
        };
        if !candidate.closed_shell_proven||!volume.volume_certified
            ||candidate.can_apply||volume.can_apply{
            return Err(invalid("G7 could not establish all Boolean prerequisites"));
        }
        Ok(G7PrepareResult::Ready(G7PreparedBoolean{
            candidate,volume,token:self.token,
            expected_revision:self.revision,feature_key:self.feature_key.clone(),
            source_stamps:[stamp(&self.state.inputs[0]),stamp(&self.state.inputs[1])],
        }))
    }
    pub fn commit(&mut self,proposal:G7PreparedBoolean)
        ->Result<G7Commit,GeometryError>{
        if proposal.token!=self.token||proposal.feature_key!=self.feature_key{
            return Err(invalid("G7 Boolean proposal belongs to another document"));
        }
        if proposal.expected_revision!=self.revision||self.state.accepted.is_some(){
            return Err(invalid("G7 proposal stale or result already committed"));
        }
        if proposal.source_stamps !=
            [stamp(&self.state.inputs[0]),stamp(&self.state.inputs[1])]{
            return Err(invalid("G7 Boolean sources changed after staging"));
        }
        // Recheck G1/G2/G7a topological evidence and G7b independent volume
        // from the current sources while holding exclusive &mut self.
        proposal.candidate.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],self.tolerance)?;
        proposal.volume.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],
            &proposal.candidate,self.tolerance)?;
        if !proposal.candidate.closed_shell_proven||!proposal.volume.volume_certified
            ||proposal.candidate.can_apply||proposal.volume.can_apply{
            return Err(invalid("G7 source-owned Boolean has not been certified"));
        }
        let next=G7State{
            inputs:self.state.inputs.clone(),
            accepted:Some(G7Accepted{candidate:proposal.candidate,volume:proposal.volume}),
        };
        next.validate(self.tolerance)?;
        let identities=G7IdentityRegistry::from_candidate(
            next.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key,self.token)?;
        let next_revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("G7 document revision exhausted"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G7 cannot reserve atomic Undo checkpoint"))?;
        let result=next.model().ok_or_else(||invalid("G7 result missing"))?;
        let report=G7Commit{
            revision_before:self.revision,
            revision_after:next_revision,
            result_volume:result.mass.volume,
            shell_components:result.shells.len(),
            topology_names:identities.entries.clone(),
        };
        // Publication barrier: no recoverable error after this point.
        let previous=std::mem::replace(&mut self.state,next);
        self.undo.push(previous);
        self.redo.clear();
        self.registry=identities;
        self.revision=next_revision;
        Ok(report)
    }
    pub fn replace_source(
        &mut self,operand:G7Operand,source:BrepBody,
    )->Result<u64,GeometryError>{
        source.validate()?;
        let mut next=self.state.clone();
        next.inputs[operand.index()]=source;
        next.accepted=None;
        next.validate(self.tolerance)?;
        let revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("G7 source revision overflow"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G7 source history checkpoint unavailable"))?;
        let previous=std::mem::replace(&mut self.state,next);
        self.undo.push(previous);
        self.redo.clear();
        self.registry=G7IdentityRegistry::empty(self.token);
        self.revision=revision;
        Ok(revision)
    }
    pub fn undo(&mut self)->Result<Option<G7Transition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<G7Transition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)->Result<Option<G7Transition>,GeometryError>{
        let stack=if redo{&self.redo}else{&self.undo};
        let Some(target)=stack.last()else{return Ok(None)};
        target.validate(self.tolerance)?;
        let registry=G7IdentityRegistry::from_candidate(
            target.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key,self.token)?;
        let after=self.revision.checked_add(1)
            .ok_or_else(||invalid("G7 history revision overflow"))?;
        if redo{
            self.undo.try_reserve(1)
                .map_err(|_|invalid("G7 Redo history reservation failed"))?;
        }else{
            self.redo.try_reserve(1)
                .map_err(|_|invalid("G7 Undo history reservation failed"))?;
        }
        let before=self.revision;
        // Transition barrier: an unchanged document on every failed check.
        let target=if redo{self.redo.pop().unwrap()}else{self.undo.pop().unwrap()};
        let old=std::mem::replace(&mut self.state,target);
        if redo{self.undo.push(old);}else{self.redo.push(old);}
        self.registry=registry;
        self.revision=after;
        Ok(Some(G7Transition{
            direction:if redo{HistoryDirection::Redo}else{HistoryDirection::Undo},
            revision_before:before,revision_after:after,
            has_result:self.state.accepted.is_some(),
        }))
    }
}
pub(crate) fn g7_restore_token()->Result<u64,GeometryError>{fresh()}
pub(crate) fn g7_valid_key(key:&str)->bool{valid_key(key)}
