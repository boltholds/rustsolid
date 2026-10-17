//! G6c: G5-revision-safe atomic document for the new two-shell Boolean family.
//! Adds a second typed G5 geometry family without changing or weakening G5b
//! schema v1's existing four-corner source and identity invariants.
//! All fallible source/volume/identity/history checks precede publication.

use std::collections::{BTreeMap,BTreeSet};
use std::sync::atomic::{AtomicU64,Ordering};
use crate::{
    AnalyticTwoWallPartitionPlan,AnalyticTwoWallPartitionResult,
    AnalyticTwoWallHold,AnalyticTwoWallVolumeOracle,
    AnalyticTwoWallVolumeResult,TwoWallVolumeHold,
    BooleanOperation,BrepBody,BrepModel,GeometryError,GeometryTolerance,
    TopologyEntity,TopologyName,HistoryDirection,VertexId,
};
static NEXT_G6_DOCUMENT:AtomicU64=AtomicU64::new(1);
fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(s.into())
}
fn fresh()->Result<u64,GeometryError>{
    let id=NEXT_G6_DOCUMENT.fetch_add(1,Ordering::Relaxed);
    if id==0||id==u64::MAX {return Err(invalid("G6c process identity exhausted"));}
    Ok(id)
}
fn valid_key(key:&str)->bool {
    !key.is_empty()&&key.len()<=96&&key.bytes().all(|b|
        b.is_ascii_alphanumeric()||b==b'-'||b==b'_')
}
fn source_stamp(body:&BrepBody)->String{
    format!("{body:?}")
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G5TwoWallHold {
    UnsupportedGeometry(AnalyticTwoWallHold),
    VolumeUncertified(TwoWallVolumeHold),
    AlreadyCommitted,
}
#[derive(Debug,Clone)]
pub enum G5TwoWallPrepareResult {
    Ready(G5TwoWallPrepared),Held(G5TwoWallHold),
}
#[derive(Debug,Clone)]
pub struct G5TwoWallPrepared {
    pub candidate:AnalyticTwoWallPartitionPlan,
    pub volume:AnalyticTwoWallVolumeOracle,
    token:u64,
    expected_revision:u64,
    feature_key:String,
    input_stamps:[String;2],
}
impl G5TwoWallPrepared{
    pub fn expected_revision(&self)->u64{self.expected_revision}
}
#[derive(Debug,Clone)]
struct Accepted{
    candidate:AnalyticTwoWallPartitionPlan,
    oracle:AnalyticTwoWallVolumeOracle,
}
#[derive(Debug,Clone)]
struct State{
    inputs:[BrepBody;2],
    accepted:Option<Accepted>,
}
impl State{
    fn model(&self)->Option<&BrepModel>{
        self.accepted.as_ref().map(|x|x.candidate.brep())
    }
    fn validate(&self,tol:GeometryTolerance)->Result<(),GeometryError>{
        self.inputs[0].validate()?;self.inputs[1].validate()?;
        if let Some(result)=&self.accepted{
            result.candidate.validate_against(&self.inputs[0],&self.inputs[1],tol)?;
            result.oracle.validate_against(
                &self.inputs[0],&self.inputs[1],&result.candidate,tol)?;
            if !result.candidate.closed_shell_proven
                ||!result.oracle.volume_certified
                ||result.candidate.can_apply||result.oracle.can_apply{
                return Err(invalid("G6c result missing independent G6b certification"));
            }
        }
        Ok(())
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash)]
pub struct G5TwoWallHandle{
    pub document_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G5TwoWallHandleError{ForeignDocument,StaleHandle}
#[derive(Debug,Clone,PartialEq,Eq)]
pub struct G5TwoWallBinding{
    pub handle:G5TwoWallHandle,
    pub name:TopologyName,
    pub entity:TopologyEntity,
}
#[derive(Debug,Clone)]
struct Registry{
    token:u64,
    named:BTreeMap<TopologyName,G5TwoWallBinding>,
    by_entity:BTreeMap<TopologyEntity,TopologyName>,
    slots:Vec<G5TwoWallBinding>,
}
impl Registry{
    fn empty(token:u64)->Self{
        Self{token,named:BTreeMap::new(),by_entity:BTreeMap::new(),slots:Vec::new()}
    }
    fn names(candidate:Option<&AnalyticTwoWallPartitionPlan>,key:&str)
        ->Result<BTreeMap<TopologyEntity,TopologyName>,GeometryError>{
        let Some(candidate)=candidate else{return Ok(BTreeMap::new())};
        candidate.brep().validate()?;
        if candidate.components.len()!=2
            ||!candidate.closed_shell_proven||candidate.can_apply{
            return Err(invalid("G6c requires exactly two fully stitched shells"));
        }
        let m=candidate.brep();
        let mut out=BTreeMap::<TopologyEntity,TopologyName>::new();
        let mut seen=BTreeSet::new();
        let mut add=|entity:TopologyEntity,path:String|->Result<(),GeometryError>{
            let name=TopologyName(path);
            if out.contains_key(&entity)||!seen.insert(name.clone()){
                return Err(invalid("G6c duplicate semantic topology identity"));
            }
            out.insert(entity,name);
            Ok(())
        };
        add(TopologyEntity::Body,format!("{key}/difference/body"))?;
        let vroles=["bottom_left_corner","bottom_left_circle",
            "bottom_right_circle","bottom_right_corner",
            "top_left_corner","top_left_circle",
            "top_right_circle","top_right_corner"];
        let eroles=["bottom_left_line","bottom_circle","bottom_right_line",
            "bottom_z_wall","top_left_line","top_circle",
            "top_right_line","top_z_wall",
            "left_corner_generator","left_cut_generator",
            "right_cut_generator","right_corner_generator"];
        let froles=["bottom","top","left_x_wall","right_x_wall",
            "outer_z_wall","inner_cylinder_wall"];
        let mut parts=BTreeSet::new();
        for c in &candidate.components {
            if ![-1_i8,1_i8].contains(&c.z_side)||!parts.insert(c.z_side){
                return Err(invalid("G6c duplicate material component label"));
            }
            let prefix=format!(
                "{key}/difference/strip/z_{}",
                if c.z_side<0{"neg"}else{"pos"}
            );
            let shell=m.shells.get(c.shell.0 as usize)
                .filter(|s|s.id==c.shell&&s.faces==c.faces)
                .ok_or_else(||invalid("G6c shell not owned by result patch"))?;
            add(TopologyEntity::Shell(shell.id),format!("{prefix}/shell"))?;
            let edge_at=|i:usize|->Result<&crate::Edge,GeometryError>{
                let id=c.edges[i];
                m.edges.get(id.0 as usize).filter(|e|e.id==id)
                    .ok_or_else(||invalid("G6c missing candidate edge"))
            };
            let vertices=[
                edge_at(0)?.start,edge_at(0)?.end,
                edge_at(2)?.start,edge_at(2)?.end,
                edge_at(4)?.start,edge_at(4)?.end,
                edge_at(6)?.start,edge_at(6)?.end,
            ];
            for (i,vertex) in vertices.into_iter().enumerate(){
                if m.vertices.get(vertex.0 as usize).is_none(){
                    return Err(invalid("G6c semantic corner vertex missing"));
                }
                add(TopologyEntity::Vertex(vertex),format!(
                    "{prefix}/vertex/{}",vroles[i]))?;
            }
            for (i,&fid) in c.faces.iter().enumerate(){
                let face=m.faces.get(fid.0 as usize)
                    .filter(|f|f.id==fid&&f.loops.len()==1)
                    .ok_or_else(||invalid("G6c missing bounded source face"))?;
                let lid=face.loops[0];
                add(TopologyEntity::Face(fid),format!("{prefix}/face/{}",froles[i]))?;
                add(TopologyEntity::Loop(lid),format!("{prefix}/loop/{}",froles[i]))?;
            }
            for (i,&eid) in c.edges.iter().enumerate(){
                let e=edge_at(i)?;
                add(TopologyEntity::Edge(eid),format!("{prefix}/edge/{}",eroles[i]))?;
                for &cid in &e.coedges{
                    let co=m.coedges.get(cid.0 as usize)
                        .filter(|co|co.id==cid&&co.edge==eid)
                        .ok_or_else(||invalid("G6c missing physical coedge"))?;
                    let role=c.faces.iter().position(|f|*f==co.face)
                        .ok_or_else(||invalid("G6c coedge crosses material shells"))?;
                    add(TopologyEntity::Coedge(cid),format!(
                        "{prefix}/edge/{}/use/{}",eroles[i],froles[role]))?;
                }
            }
        }
        let expected=1+m.vertices.len()+m.edges.len()+m.coedges.len()
            +m.faces.len()+m.loops.len()+m.shells.len();
        if parts.len()!=2||out.len()!=expected||expected!=115{
            return Err(invalid("G6c incomplete semantic topology identity"));
        }
        Ok(out)
    }
    fn from_candidate(
        candidate:Option<&AnalyticTwoWallPartitionPlan>,key:&str,token:u64,
    )->Result<Self,GeometryError>{
        let names=Self::names(candidate,key)?;
        let mut out=Self::empty(token);
        out.slots.try_reserve(names.len())
            .map_err(|_|invalid("G6c unable to allocate identity slots"))?;
        for (entity,name) in names{
            let slot=u32::try_from(out.slots.len())
                .map_err(|_|invalid("G6c topology slots exhausted"))?;
            let binding=G5TwoWallBinding{
                handle:G5TwoWallHandle {
                    document_token:token,slot,generation:fresh()?,
                },
                name:name.clone(),entity,
            };
            out.named.insert(name.clone(),binding.clone());
            out.by_entity.insert(entity,name);
            out.slots.push(binding);
        }
        Ok(out)
    }
    fn resolve(&self,handle:G5TwoWallHandle)
        ->Result<G5TwoWallBinding,G5TwoWallHandleError>{
        if handle.document_token!=self.token{
            return Err(G5TwoWallHandleError::ForeignDocument);
        }
        self.slots.get(handle.slot as usize)
            .filter(|b|b.handle==handle).cloned()
            .ok_or(G5TwoWallHandleError::StaleHandle)
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G5TwoWallOperand{First,Second}
impl G5TwoWallOperand{
    fn index(self)->usize{match self{Self::First=>0,Self::Second=>1}}
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct G5TwoWallTransition{
    pub direction:HistoryDirection,
    pub revision_before:u64,pub revision_after:u64,
    pub has_result:bool,
}
#[derive(Debug,Clone)]
pub struct G5TwoWallCommit{
    pub revision_before:u64,pub revision_after:u64,
    pub volume:f64,pub shells:usize,
    pub identities:Vec<G5TwoWallBinding>,
}
pub struct G5TwoWallDocument{
    token:u64,feature_key:String,revision:u64,
    tolerance:GeometryTolerance,state:State,
    registry:Registry,undo:Vec<State>,redo:Vec<State>,
}
impl G5TwoWallDocument {
    pub fn new(
        key:&str,first:BrepBody,second:BrepBody,tol:GeometryTolerance,
    )->Result<Self,GeometryError>{
        if !valid_key(key){return Err(invalid("G6c invalid document feature key"));}
        tol.validate()?;first.validate()?;second.validate()?;
        let token=fresh()?;
        Ok(Self {
            token,feature_key:key.into(),revision:0,tolerance:tol,
            state:State{inputs:[first,second],accepted:None},
            registry:Registry::empty(token),undo:Vec::new(),redo:Vec::new(),
        })
    }
    pub fn feature_key(&self)->&str{&self.feature_key}
    pub fn revision(&self)->u64{self.revision}
    pub fn token(&self)->u64{self.token}
    pub fn model(&self)->Option<&BrepModel>{self.state.model()}
    pub fn source(&self,which:G5TwoWallOperand)->&BrepBody{
        &self.state.inputs[which.index()]
    }
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn bindings(&self)->Vec<G5TwoWallBinding>{self.registry.slots.clone()}
    pub fn resolve_name(&self,name:&TopologyName)->Option<G5TwoWallBinding>{
        self.registry.named.get(name).cloned()
    }
    pub fn resolve_handle(&self,handle:G5TwoWallHandle)
        ->Result<G5TwoWallBinding,G5TwoWallHandleError>{
        self.registry.resolve(handle)
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.state.validate(self.tolerance)?;
        let names=Registry::names(
            self.state.accepted.as_ref().map(|x|&x.candidate),
            &self.feature_key,
        )?;
        if names!=self.registry.by_entity
            ||names.len()!=self.registry.slots.len()
            ||self.registry.named.len()!=self.registry.slots.len(){
            return Err(invalid("G6c live identity registry inconsistent"));
        }
        for (i,b) in self.registry.slots.iter().enumerate(){
            if b.handle.document_token!=self.token
                ||b.handle.slot as usize!=i
                ||b.handle.generation==0
                ||self.registry.named.get(&b.name)!=Some(b){
                return Err(invalid("G6c corrupted scoped identity registry"));
            }
        }
        Ok(())
    }
    pub fn prepare(&self)->Result<G5TwoWallPrepareResult,GeometryError>{
        if self.state.accepted.is_some(){
            return Ok(G5TwoWallPrepareResult::Held(G5TwoWallHold::AlreadyCommitted));
        }
        let candidate=match AnalyticTwoWallPartitionPlan::build(
            &self.state.inputs[0],&self.state.inputs[1],
            BooleanOperation::Difference,self.tolerance)?
        {
            AnalyticTwoWallPartitionResult::Ready(x)=>x,
            AnalyticTwoWallPartitionResult::Held(reason)=>
                return Ok(G5TwoWallPrepareResult::Held(
                    G5TwoWallHold::UnsupportedGeometry(reason))),
        };
        let volume=match AnalyticTwoWallVolumeOracle::build(
            &self.state.inputs[0],&self.state.inputs[1],
            &candidate,self.tolerance)?
        {
            AnalyticTwoWallVolumeResult::Ready(x)=>x,
            AnalyticTwoWallVolumeResult::Held(reason)=>
                return Ok(G5TwoWallPrepareResult::Held(
                    G5TwoWallHold::VolumeUncertified(reason))),
        };
        if !candidate.closed_shell_proven||!volume.volume_certified
            ||candidate.can_apply||volume.can_apply{
            return Err(invalid("G6c cannot prepare unproven material"));
        }
        Ok(G5TwoWallPrepareResult::Ready(G5TwoWallPrepared{
            candidate,volume,token:self.token,
            expected_revision:self.revision,
            feature_key:self.feature_key.clone(),
            input_stamps:[
                source_stamp(&self.state.inputs[0]),
                source_stamp(&self.state.inputs[1]),
            ],
        }))
    }
    pub fn commit(&mut self,proposal:G5TwoWallPrepared)
        ->Result<G5TwoWallCommit,GeometryError>{
        if proposal.token!=self.token||proposal.feature_key!=self.feature_key{
            return Err(invalid("G6c proposal belongs to a different document"));
        }
        if proposal.expected_revision!=self.revision||self.state.accepted.is_some(){
            return Err(invalid("G6c stale material source revision"));
        }
        let fresh=[
            source_stamp(&self.state.inputs[0]),
            source_stamp(&self.state.inputs[1]),
        ];
        if fresh!=proposal.input_stamps {
            return Err(invalid("G6c original material changed before commit"));
        }
        // Exclusive &mut owner and full source revalidation before publication.
        proposal.candidate.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],self.tolerance)?;
        proposal.volume.validate_against(
            &self.state.inputs[0],&self.state.inputs[1],
            &proposal.candidate,self.tolerance)?;
        if !proposal.volume.volume_certified||!proposal.candidate.closed_shell_proven
            ||proposal.volume.can_apply||proposal.candidate.can_apply{
            return Err(invalid("G6c refuses an unverified Boolean result"));
        }
        let accepted=Accepted {candidate:proposal.candidate,oracle:proposal.volume};
        let next_state=State {
            inputs:self.state.inputs.clone(),accepted:Some(accepted),
        };
        next_state.validate(self.tolerance)?;
        let registry=Registry::from_candidate(
            next_state.accepted.as_ref().map(|x|&x.candidate),
            &self.feature_key,self.token)?;
        let after=self.revision.checked_add(1)
            .ok_or_else(||invalid("G6c document revision overflow"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G6c history reservation failed"))?;
        let result=next_state.model().unwrap();
        let report=G5TwoWallCommit{
            revision_before:self.revision,revision_after:after,
            volume:result.mass.volume,shells:result.shells.len(),
            identities:registry.slots.clone(),
        };
        // COMMIT BARRIER. No fallible operations after this point.
        let old=std::mem::replace(&mut self.state,next_state);
        self.undo.push(old);
        self.redo.clear();
        self.registry=registry;
        self.revision=after;
        Ok(report)
    }
    pub fn replace_source(
        &mut self,which:G5TwoWallOperand,source:BrepBody,
    )->Result<u64,GeometryError>{
        source.validate()?;
        let mut state=self.state.clone();
        state.inputs[which.index()]=source;
        state.accepted=None;
        state.validate(self.tolerance)?;
        let next=self.revision.checked_add(1)
            .ok_or_else(||invalid("G6c source revision overflow"))?;
        self.undo.try_reserve(1)
            .map_err(|_|invalid("G6c source undo checkpoint reserve failed"))?;
        let previous=std::mem::replace(&mut self.state,state);
        self.undo.push(previous);self.redo.clear();
        self.registry=Registry::empty(self.token);
        self.revision=next;
        Ok(next)
    }
    pub fn undo(&mut self)->Result<Option<G5TwoWallTransition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<G5TwoWallTransition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)
        ->Result<Option<G5TwoWallTransition>,GeometryError>{
        let stack=if redo{&self.redo}else{&self.undo};
        let Some(previous)=stack.last()else{return Ok(None)};
        previous.validate(self.tolerance)?;
        let registry=Registry::from_candidate(
            previous.accepted.as_ref().map(|x|&x.candidate),
            &self.feature_key,self.token)?;
        let after=self.revision.checked_add(1)
            .ok_or_else(||invalid("G6c history revision overflow"))?;
        if redo {
            self.undo.try_reserve(1)
                .map_err(|_|invalid("G6c cannot reserve redo checkpoint"))?;
        }else{
            self.redo.try_reserve(1)
                .map_err(|_|invalid("G6c cannot reserve undo checkpoint"))?;
        }
        let before=self.revision;
        // COMMIT BARRIER.
        let restored=if redo{self.redo.pop().unwrap()}else{self.undo.pop().unwrap()};
        let old=std::mem::replace(&mut self.state,restored);
        if redo{self.undo.push(old);}else{self.redo.push(old);}
        self.registry=registry;self.revision=after;
        Ok(Some(G5TwoWallTransition {
            direction:if redo{HistoryDirection::Redo}else{HistoryDirection::Undo},
            revision_before:before,revision_after:after,
            has_result:self.state.accepted.is_some(),
        }))
    }
}
