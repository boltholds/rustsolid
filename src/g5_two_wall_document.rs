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


// G6d persistent format is version 2 and explicitly NOT the original G5b
// four-corner archive schema 1. Both are source-rebuilt and fail-closed.
use serde::{Serialize,Deserialize};
use std::path::Path;
use std::fs::{self,File,OpenOptions};
use std::io::{Read,Write};
use crate::{G5ProjectSource,G5FeatureNode,G5FeatureGraph,G5GraphBinding};
const G6_ARCHIVE_SCHEMA:u32=2;
const G6_ARCHIVE_BYTES:usize=1024*1024;
const G6_HISTORY_CAP:usize=64;
static NEXT_G6_FILE:AtomicU64=AtomicU64::new(1);

#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5TwoWallVolumeStamp {
    pub mass_volume:f64,
    pub component_volumes:[f64;2],
    pub source_lower:f64,
    pub source_upper:f64,
    pub flux_volume:f64,
}
impl G5TwoWallVolumeStamp{
    fn from_accepted(a:&Accepted)->Result<Self,GeometryError>{
        if a.candidate.components.len()!=2 {
            return Err(invalid("G6d archive accepts precisely two components"));
        }
        Ok(Self{
            mass_volume:a.candidate.brep().mass.volume,
            component_volumes:[
                a.candidate.components[0].analytic_volume,
                a.candidate.components[1].analytic_volume,
            ],
            source_lower:a.oracle.source_total.lower,
            source_upper:a.oracle.source_total.upper,
            flux_volume:a.oracle.boundary_total,
        })
    }
    fn valid(&self)->bool{
        self.mass_volume.is_finite()&&self.mass_volume>0.
            &&self.component_volumes.iter().all(|x|x.is_finite()&&*x>0.)
            &&self.source_lower.is_finite()&&self.source_lower>0.
            &&self.source_upper.is_finite()&&self.source_upper>self.source_lower
            &&self.flux_volume.is_finite()&&self.flux_volume>0.
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5TwoWallArchiveState {
    pub sources:[G5ProjectSource;2],
    pub result:Option<G5TwoWallVolumeStamp>,
}
impl G5TwoWallArchiveState{
    fn from_state(s:&State)->Result<Self,GeometryError>{
        if let BrepBody::Polyhedral(first)=&s.inputs[0]{
            if first.revision!=0{
                return Err(invalid("G6d does not yet persist nonzero source-solid revisions"));
            }
        }
        Ok(Self{
            sources:[
                G5ProjectSource::from_body(&s.inputs[0])?,
                G5ProjectSource::from_body(&s.inputs[1])?,
            ],
            result:s.accepted.as_ref()
                .map(G5TwoWallVolumeStamp::from_accepted)
                .transpose()?,
        })
    }
    fn rebuild(&self,tol:GeometryTolerance)->Result<State,GeometryError>{
        let inputs=[self.sources[0].to_body()?,self.sources[1].to_body()?];
        let accepted=if let Some(record)=&self.result {
            if !record.valid(){
                return Err(invalid("G6d invalid stored mass certificate"));
            }
            let candidate=match AnalyticTwoWallPartitionPlan::build(
                &inputs[0],&inputs[1],BooleanOperation::Difference,tol)?
            {
                AnalyticTwoWallPartitionResult::Ready(x)=>x,
                AnalyticTwoWallPartitionResult::Held(_)=>
                    return Err(invalid("G6d archived two-wall geometry no longer certifies")),
            };
            let oracle=match AnalyticTwoWallVolumeOracle::build(
                &inputs[0],&inputs[1],&candidate,tol)?
            {
                AnalyticTwoWallVolumeResult::Ready(x)=>x,
                AnalyticTwoWallVolumeResult::Held(_)=>
                    return Err(invalid("G6d archived independent mass no longer certifies")),
            };
            let next=Accepted{candidate,oracle};
            let certified=G5TwoWallVolumeStamp::from_accepted(&next)?;
            if certified!=*record{
                return Err(invalid("G6d archived mass differs from reconstructed B-rep"));
            }
            Some(next)
        }else{None};
        let state=State{inputs,accepted};
        state.validate(tol)?;
        Ok(state)
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5TwoWallProjectArchive{
    pub schema_version:u32,
    pub feature_key:String,
    pub document_revision:u64,
    pub tolerance:GeometryTolerance,
    pub graph:G5FeatureGraph,
    pub current:G5TwoWallArchiveState,
    pub undo:Vec<G5TwoWallArchiveState>,
    pub redo:Vec<G5TwoWallArchiveState>,
}
impl G5TwoWallProjectArchive{
    fn graph_for(
        document:&G5TwoWallDocument,
    )->Result<G5FeatureGraph,GeometryError>{
        let accepted=document.state.accepted.iter()
            .chain(document.undo.iter().filter_map(|s|s.accepted.as_ref()))
            .chain(document.redo.iter().filter_map(|s|s.accepted.as_ref()))
            .next();
        let names=Registry::names(
            accepted.map(|a|&a.candidate),&document.feature_key)?;
        if !names.is_empty()&&names.len()!=115 {
            return Err(invalid("G6d partial semantic graph not supported"));
        }
        let mut bindings=names.into_values().map(|name|
            G5GraphBinding{topology_name:name,producing_node:"difference".into()}
        ).collect::<Vec<_>>();
        bindings.sort_by(|a,b|a.topology_name.cmp(&b.topology_name));
        Ok(G5FeatureGraph{
            nodes:vec![
                G5FeatureNode::Source{id:"source_a".into(),operand_index:0},
                G5FeatureNode::Source{id:"source_b".into(),operand_index:1},
                G5FeatureNode::Difference{
                    id:"difference".into(),left:"source_a".into(),
                    right:"source_b".into(),
                    output_scope:format!("{}/difference",document.feature_key),
                },
            ],
            output_bindings:bindings,
        })
    }
    pub fn from_document(doc:&G5TwoWallDocument)
        ->Result<Self,GeometryError>{
        doc.validate()?;
        if doc.undo.len()+doc.redo.len()>G6_HISTORY_CAP{
            return Err(invalid("G6d persistent history exceeds 64 states"));
        }
        let archive=Self {
            schema_version:G6_ARCHIVE_SCHEMA,
            feature_key:doc.feature_key.clone(),
            document_revision:doc.revision,
            tolerance:doc.tolerance,
            graph:Self::graph_for(doc)?,
            current:G5TwoWallArchiveState::from_state(&doc.state)?,
            undo:doc.undo.iter().map(G5TwoWallArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
            redo:doc.redo.iter().map(G5TwoWallArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
        };
        let _=archive.to_json()?;
        Ok(archive)
    }
    pub fn reopen(&self)->Result<G5TwoWallDocument,GeometryError>{
        if self.schema_version!=G6_ARCHIVE_SCHEMA
            ||!valid_key(&self.feature_key)
            ||self.document_revision==u64::MAX
            ||self.document_revision < self.undo.len() as u64
            ||self.undo.len()+self.redo.len()>G6_HISTORY_CAP{
            return Err(invalid("G6d archive version/history does not match policy"));
        }
        self.tolerance.validate()?;
        let current=self.current.rebuild(self.tolerance)?;
        let undo=self.undo.iter()
            .map(|s|s.rebuild(self.tolerance))
            .collect::<Result<Vec<_>,_>>()?;
        let redo=self.redo.iter()
            .map(|s|s.rebuild(self.tolerance))
            .collect::<Result<Vec<_>,_>>()?;
        let mut temporal=undo.iter().collect::<Vec<_>>();
        temporal.push(&current);
        temporal.extend(redo.iter().rev());
        for pair in temporal.windows(2){
            if pair[0].accepted.is_some()&&pair[1].accepted.is_some(){
                return Err(invalid("G6d archived graph skips required source edit"));
            }
        }
        let token=fresh()?;
        let registry=Registry::from_candidate(
            current.accepted.as_ref().map(|a|&a.candidate),
            &self.feature_key,token)?;
        let document=G5TwoWallDocument{
            token,feature_key:self.feature_key.clone(),
            revision:self.document_revision,tolerance:self.tolerance,
            state:current,registry,undo,redo,
        };
        document.validate()?;
        if Self::graph_for(&document)?!=self.graph{
            return Err(invalid("G6d archived producer graph/topology names differ"));
        }
        Ok(document)
    }
    pub fn to_json(&self)->Result<String,GeometryError>{
        if self.schema_version!=G6_ARCHIVE_SCHEMA{
            return Err(invalid("G6d only schema version 2 can be serialized"));
        }
        let data=serde_json::to_string_pretty(self)
            .map_err(|e|invalid(format!("G6d serialize: {e}")))?;
        if data.len()>G6_ARCHIVE_BYTES{
            return Err(invalid("G6d archive exceeds one MiB"));
        }
        Ok(data)
    }
    pub fn from_json(raw:&str)->Result<Self,GeometryError>{
        if raw.len()>G6_ARCHIVE_BYTES{
            return Err(invalid("G6d archive exceeds one MiB"));
        }
        let archive:Self=serde_json::from_str(raw)
            .map_err(|e|invalid(format!("G6d invalid project document: {e}")))?;
        // Never trust stored material or names.
        let _=archive.reopen()?;
        Ok(archive)
    }
    pub fn save_atomic(&self,path:impl AsRef<Path>)
        ->Result<(),GeometryError>{
        let body=self.to_json()?;
        let _=self.reopen()?;
        let dest=path.as_ref();
        let parent=dest.parent().filter(|p|!p.as_os_str().is_empty())
            .unwrap_or_else(||Path::new("."));
        if !parent.is_dir(){
            return Err(invalid("G6d missing destination directory"));
        }
        let name=dest.file_name().ok_or_else(||
            invalid("G6d missing destination name"))?;
        let ordinal=NEXT_G6_FILE.fetch_add(1,Ordering::Relaxed);
        if ordinal==0||ordinal==u64::MAX{
            return Err(invalid("G6d temp file counter exhausted"));
        }
        let candidate=parent.join(format!(
            ".{}.g6d-{}-{}.tmp",
            name.to_string_lossy(),std::process::id(),ordinal
        ));
        let result=(||->Result<(),GeometryError>{
            let mut f=OpenOptions::new()
                .create_new(true).write(true).open(&candidate)
                .map_err(|e|invalid(format!("G6d open temp: {e}")))?;
            f.write_all(body.as_bytes())
                .map_err(|e|invalid(format!("G6d write temp: {e}")))?;
            f.sync_all()
                .map_err(|e|invalid(format!("G6d flush temp: {e}")))?;
            drop(f);
            fs::rename(&candidate,dest)
                .map_err(|e|invalid(format!("G6d rename: {e}")))?;
            Ok(())
        })();
        if result.is_err(){let _=fs::remove_file(&candidate);}
        result?;
        let _=File::open(parent).and_then(|folder|folder.sync_all());
        Ok(())
    }
    pub fn load_file(path:impl AsRef<Path>)->Result<Self,GeometryError>{
        let mut f=File::open(path)
            .map_err(|e|invalid(format!("G6d open project: {e}")))?;
        let mut bounded=(&mut f).take((G6_ARCHIVE_BYTES+1) as u64);
        let mut bytes=Vec::new();
        bounded.read_to_end(&mut bytes)
            .map_err(|e|invalid(format!("G6d read project: {e}")))?;
        if bytes.len()>G6_ARCHIVE_BYTES{
            return Err(invalid("G6d project exceeds one MiB"));
        }
        let text=std::str::from_utf8(&bytes)
            .map_err(|_|invalid("G6d project is not UTF-8"))?;
        Self::from_json(text)
    }
}
