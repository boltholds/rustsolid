//! Versioned sparse topology registry and strong-atomic modeling transactions.
//!
//! This first arena layer owns the existing validated polyhedral CommandHistory
//! and exposes analytical cylinders as read-only B-rep snapshots. The registry
//! identifies entities by deterministic semantic TopologyName rather than dense
//! array position; surviving names keep handles across compaction, while retired
//! slots increase their generation before they can be reused. Neither a handle
//! nor an ID from this process is a durable serialized CAD feature identity.
//!
//! Edits execute on an isolated clone of CommandHistory plus registry, and
//! become observable only when both the modeling command and reconciliation
//! have passed. This intentionally favors correctness over edit throughput;
//! native persistent sparse editable curved geometry is a later milestone.
use crate::{
    BrepBody, BrepModel, CommandBatch, CommandHistory, CylindricalBrep,
    EditReport, EulerCommand, GeometryError, HistoryTransition, ModelUnits,
    Solid, TopologyDeltaPlan, TopologyEntity, TopologyName,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ARENA: AtomicU64 = AtomicU64::new(1);

fn invalid(msg:impl Into<String>)->GeometryError {GeometryError::InvalidEdit(msg.into())}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ArenaHandle {
    pub arena_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArenaHandleError {ForeignArena,StaleHandle}
impl fmt::Display for ArenaHandleError {
    fn fmt(&self,f:&mut fmt::Formatter<'_>)->fmt::Result {
        match self {
            Self::ForeignArena=>write!(f,"handle belongs to another BrepArena"),
            Self::StaleHandle=>write!(f,"handle refers to a retired arena incarnation"),
        }
    }
}
impl std::error::Error for ArenaHandleError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArenaBinding {
    pub handle:ArenaHandle,
    pub entity:TopologyEntity,
    pub name:TopologyName,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArenaRemap {
    pub handle:ArenaHandle,
    pub name:TopologyName,
    pub before:TopologyEntity,
    pub after:TopologyEntity,
}

#[derive(Debug, Clone)]
struct SparseSlot {
    generation:u64,
    current:Option<(TopologyEntity,TopologyName)>,
}
#[derive(Debug, Clone)]
struct Registry {
    token:u64,
    slots:Vec<SparseSlot>,
    free:Vec<u32>,
    by_name:BTreeMap<TopologyName,ArenaHandle>,
    by_entity:BTreeMap<TopologyEntity,ArenaHandle>,
}
impl Registry {
    fn empty()->Self {
        let token=NEXT_ARENA.fetch_add(1,Ordering::Relaxed);
        assert!(token!=0 && token!=u64::MAX,"arena token generator exhausted");
        Self{token,slots:Vec::new(),free:Vec::new(),by_name:BTreeMap::new(),by_entity:BTreeMap::new()}
    }
    fn binding(&self,handle:ArenaHandle)->Result<ArenaBinding,ArenaHandleError>{
        if handle.arena_token!=self.token{return Err(ArenaHandleError::ForeignArena);}
        let slot=self.slots.get(handle.slot as usize).ok_or(ArenaHandleError::StaleHandle)?;
        let (entity,name)=slot.current.as_ref().ok_or(ArenaHandleError::StaleHandle)?;
        if slot.generation!=handle.generation{return Err(ArenaHandleError::StaleHandle);}
        Ok(ArenaBinding{handle,entity:*entity,name:name.clone()})
    }
    fn handle(&self,entity:TopologyEntity)->Option<ArenaHandle>{self.by_entity.get(&entity).copied()}
    fn by_name(&self,name:&TopologyName)->Option<ArenaHandle>{self.by_name.get(name).copied()}
    fn live(&self)->Vec<ArenaBinding> {
        self.by_name.values().map(|&handle|self.binding(handle).expect("live map must be consistent")).collect()
    }
    fn reconcile(&mut self,named:&[(TopologyEntity,TopologyName)])->Result<(),GeometryError>{
        let mut identities=BTreeMap::new();
        let mut entities=BTreeSet::new();
        for (entity,name) in named {
            if identities.insert(name.clone(),*entity).is_some() || !entities.insert(*entity){
                return Err(invalid("duplicate topological name or entity in registry candidate"));
            }
        }
        // Retire removed identities first. The free list can then reuse slots
        // safely because retirement increases generation, eliminating ABA.
        let retired=self.by_name.keys().filter(|name|!identities.contains_key(*name))
            .cloned().collect::<Vec<_>>();
        for name in retired {
            let handle=self.by_name.remove(&name).expect("registered name");
            let slot=&mut self.slots[handle.slot as usize];
            slot.current=None;
            slot.generation=slot.generation.checked_add(1)
                .ok_or_else(||invalid("arena handle generation exhausted"))?;
            self.free.push(handle.slot);
        }
        self.by_entity.clear();
        for (name,entity) in identities {
            let handle=if let Some(&old)=self.by_name.get(&name) {
                let slot=&mut self.slots[old.slot as usize];
                slot.current=Some((entity,name.clone()));
                old
            }else {
                let index=if let Some(i)=self.free.pop(){i}else{
                    let next=u32::try_from(self.slots.len())
                        .map_err(|_|invalid("arena slot capacity exhausted"))?;
                    self.slots.push(SparseSlot{generation:1,current:None});next
                };
                let slot=&mut self.slots[index as usize];
                if slot.generation==0 {return Err(invalid("invalid zero arena generation"));}
                slot.current=Some((entity,name.clone()));
                let handle=ArenaHandle{arena_token:self.token,slot:index,generation:slot.generation};
                self.by_name.insert(name.clone(),handle);
                handle
            };
            if self.by_entity.insert(entity,handle).is_some(){
                return Err(invalid("two semantic names resolve to one topology entity"));
            }
        }
        if self.by_entity.len()!=named.len(){return Err(invalid("incomplete arena identity reconciliation"));}
        Ok(())
    }
}

enum ArenaKernel {
    Editable(CommandHistory),
    CurvedReadOnly(CylindricalBrep),
}
impl ArenaKernel {
    fn revision(&self)->u64{match self{Self::Editable(h)=>h.solid().revision,Self::CurvedReadOnly(_)=>0}}
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        match self {
            Self::Editable(h)=>BrepModel::from_polyhedral(h.solid(),ModelUnits::default()),
            Self::CurvedReadOnly(body)=>BrepModel::from_cylinder(body),
        }
    }
    fn named(&self)->Result<Vec<(TopologyEntity,TopologyName)>,GeometryError>{
        let model=self.brep_model()?;
        let mut entities=Vec::new();
        entities.push(TopologyEntity::Body);
        entities.extend(model.vertices.iter().map(|x|TopologyEntity::Vertex(x.id)));
        entities.extend(model.edges.iter().map(|x|TopologyEntity::Edge(x.id)));
        entities.extend(model.coedges.iter().map(|x|TopologyEntity::Coedge(x.id)));
        entities.extend(model.loops.iter().map(|x|TopologyEntity::Loop(x.id)));
        entities.extend(model.faces.iter().map(|x|TopologyEntity::Face(x.id)));
        entities.extend(model.shells.iter().map(|x|TopologyEntity::Shell(x.id)));
        let mut out=Vec::with_capacity(entities.len());
        for e in entities {
            let name=match self {
                Self::Editable(h)=>h.solid().topology_name(e)
                    .ok_or_else(||invalid(format!("missing semantic name for {e:?}")))?.clone(),
                Self::CurvedReadOnly(_)=>{
                    let (ty,idx)=match e {
                        TopologyEntity::Body=>("body",0),
                        TopologyEntity::Vertex(x)=>("vertex",x.0),
                        TopologyEntity::Edge(x)=>("edge",x.0),
                        TopologyEntity::Coedge(x)=>("coedge",x.0),
                        TopologyEntity::Loop(x)=>("loop",x.0),
                        TopologyEntity::Face(x)=>("face",x.0),
                        TopologyEntity::Shell(x)=>("shell",x.0),
                    };
                    TopologyName(format!("analytic_cylinder/{ty}/{idx}"))
                },
            };
            out.push((e,name));
        }
        Ok(out)
    }
}

/// An immutable revision-tagged geometric view. The arena never returns
/// mutable references to its canonical topology or internal sparse slots.
#[derive(Debug, Clone)]
pub struct ArenaSnapshot {
    pub arena_token:u64,
    pub revision:u64,
    pub model:BrepModel,
    pub handles:BTreeMap<TopologyEntity,ArenaHandle>,
}
impl ArenaSnapshot {
    pub fn handle(&self,entity:TopologyEntity)->Option<ArenaHandle>{self.handles.get(&entity).copied()}
}

/// Full committed topology plus sparse, generational semantic identities.
/// The polyhedral body is editable using the validated Command/Memento backend;
/// curved geometry may be queried, but not changed by this first arena slice.
pub struct BrepArena {
    storage:ArenaKernel,
    registry:Registry,
}
impl BrepArena {
    pub fn new(body:BrepBody)->Result<Self,GeometryError>{
        let storage=match body {
            BrepBody::Polyhedral(solid)=>ArenaKernel::Editable(CommandHistory::new(solid)?),
            BrepBody::Cylinder(cylinder)=>{
                cylinder.validate()?;ArenaKernel::CurvedReadOnly(cylinder)
            }
        };
        let mut registry=Registry::empty();
        registry.reconcile(&storage.named()?)?;
        Ok(Self{storage,registry})
    }
    pub fn from_solid(body:Solid)->Result<Self,GeometryError>{Self::new(BrepBody::Polyhedral(body))}
    pub fn from_cylinder(body:CylindricalBrep)->Result<Self,GeometryError>{Self::new(BrepBody::Cylinder(body))}
    pub fn arena_token(&self)->u64{self.registry.token}
    pub fn revision(&self)->u64{self.storage.revision()}
    pub fn is_editable(&self)->bool{matches!(&self.storage,ArenaKernel::Editable(_))}
    pub fn model(&self)->Result<BrepModel,GeometryError>{self.storage.brep_model()}
    pub fn snapshot(&self)->Result<ArenaSnapshot,GeometryError>{
        let model=self.model()?;
        Ok(ArenaSnapshot{arena_token:self.arena_token(),revision:self.revision(),model,
            handles:self.registry.by_entity.clone()})
    }
    pub fn check_snapshot(&self,s:&ArenaSnapshot)->Result<(),GeometryError>{
        if s.arena_token!=self.arena_token() {return Err(invalid("foreign arena snapshot"));}
        if s.revision!=self.revision(){return Err(invalid("stale arena snapshot revision"));}
        Ok(())
    }
    pub fn resolve(&self,handle:ArenaHandle)->Result<ArenaBinding,ArenaHandleError>{self.registry.binding(handle)}
    pub fn handle(&self,entity:TopologyEntity)->Option<ArenaHandle>{self.registry.handle(entity)}
    pub fn handle_name(&self,name:&TopologyName)->Option<ArenaHandle>{self.registry.by_name(name)}
    pub fn bindings(&self)->Vec<ArenaBinding>{self.registry.live()}
    pub fn live_slots(&self)->usize{self.registry.by_entity.len()}
    pub fn allocated_slots(&self)->usize{self.registry.slots.len()}
    pub fn undo_depth(&self)->usize{match &self.storage {ArenaKernel::Editable(h)=>h.undo_depth(),_=>0}}
    pub fn redo_depth(&self)->usize{match &self.storage {ArenaKernel::Editable(h)=>h.redo_depth(),_=>0}}
    pub fn validate(&self)->Result<(),GeometryError>{
        let named=self.storage.named()?;
        if named.len()!=self.live_slots(){return Err(invalid("arena live count mismatches B-rep"));}
        for (entity,name) in named {
            let handle=self.handle(entity).ok_or_else(||invalid("unbound B-rep entity"))?;
            let binding=self.resolve(handle).map_err(|e|invalid(e.to_string()))?;
            if binding.name!=name || binding.entity!=entity || self.handle_name(&name)!=Some(handle){
                return Err(invalid("arena identity inconsistent with B-rep topology"));
            }
        }
        Ok(())
    }
    pub fn begin_transaction(&mut self)->Result<ArenaTransaction<'_>,GeometryError>{
        if !self.is_editable(){return Err(invalid("analytical curved B-rep is read-only in v0.17"));}
        Ok(ArenaTransaction{arena:self,ops:Vec::new()})
    }
    /// The v0.16 intersection plan is a mathematical *candidate*, not a
    /// committed topology delta. It cannot bypass missing region classification,
    /// curved trim editor or the explicit can_apply=false guard.
    pub fn require_executable_plan(&self,plan:&TopologyDeltaPlan)->Result<(),GeometryError>{
        if !plan.can_apply || !plan.requirements.is_empty(){
            return Err(invalid("TopologyDeltaPlan is preview-only and has unmet prerequisites"));
        }
        Err(invalid("generalized curved/Boolean delta executor is not implemented"))
    }
    fn commit_ops(&mut self,ops:Vec<ArenaOperation>)->Result<ArenaCommit,GeometryError>{
        if ops.is_empty(){return Err(invalid("empty arena transaction"));}
        let before=self.registry.clone();
        let mut batch=Vec::with_capacity(ops.len());
        for op in ops {
            let named=|h:ArenaHandle,kind:&str|->Result<TopologyName,GeometryError>{
                let binding=self.resolve(h).map_err(|err|invalid(err.to_string()))?;
                let found=match binding.entity {
                    TopologyEntity::Edge(_)=>"edge",TopologyEntity::Vertex(_)=>"vertex",
                    TopologyEntity::Face(_)=>"face",_=>"other",
                };
                if found!=kind{return Err(invalid(format!("expected {kind}, found {found}")));}
                Ok(binding.name)
            };
            batch.push(match op {
                ArenaOperation::SplitEdge{edge,fraction}=>EulerCommand::MakeEdgeVertex{
                    edge:named(edge,"edge")?,fraction},
                ArenaOperation::SplitFace{face,start,end}=>EulerCommand::MakeEdgeFace{
                    face:named(face,"face")?,start:named(start,"vertex")?,end:named(end,"vertex")?},
                ArenaOperation::KillEdgeVertex{vertex}=>EulerCommand::KillEdgeVertex{
                    vertex:named(vertex,"vertex")?},
                ArenaOperation::KillEdgeFace{edge,removed_face}=>EulerCommand::KillEdgeFace{
                    edge:named(edge,"edge")?,removed_face:named(removed_face,"face")?},
            });
        }
        let ArenaKernel::Editable(current)=&self.storage else{return Err(invalid("arena is read-only"));};
        let mut candidate=current.clone();
        let report=candidate.execute(CommandBatch{commands:batch})?;
        // Reconcile *before* the swap: any error, including registry allocation
        // or identity mismatch, preserves BOTH the body and identity indexes.
        let storage=ArenaKernel::Editable(candidate);
        let mut registry=before.clone();
        registry.reconcile(&storage.named()?)?;
        let summary=commit_summary(&before,&registry,report);
        self.storage=storage;
        self.registry=registry;
        debug_assert!(self.validate().is_ok());
        Ok(summary)
    }
    pub fn undo(&mut self)->Result<Option<HistoryTransition>,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<Option<HistoryTransition>,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<Option<HistoryTransition>,GeometryError>{
        let ArenaKernel::Editable(current)=&self.storage else {return Ok(None);};
        let mut candidate=current.clone();
        let transition=if redo{candidate.redo()?}else{candidate.undo()?};
        let Some(value)=transition else{return Ok(None)};
        let next=ArenaKernel::Editable(candidate);
        let mut registry=self.registry.clone();
        registry.reconcile(&next.named()?)?;
        self.storage=next;
        self.registry=registry;
        debug_assert!(self.validate().is_ok());
        Ok(Some(value))
    }
}

fn commit_summary(before:&Registry,after:&Registry,report:EditReport)->ArenaCommit{
    let old=before.by_name.keys().cloned().collect::<BTreeSet<_>>();
    let new=after.by_name.keys().cloned().collect::<BTreeSet<_>>();
    let created=new.difference(&old).map(|n|after.binding(after.by_name(n).unwrap()).unwrap()).collect();
    let retired=old.difference(&new).map(|n|before.binding(before.by_name(n).unwrap()).unwrap()).collect();
    let remapped=old.intersection(&new).filter_map(|n|{
        let a=before.binding(before.by_name(n).unwrap()).unwrap();
        let b=after.binding(after.by_name(n).unwrap()).unwrap();
        (a.entity!=b.entity).then_some(ArenaRemap{handle:b.handle,name:n.clone(),before:a.entity,after:b.entity})
    }).collect();
    ArenaCommit{report,created,retired,remapped}
}
#[derive(Debug,Clone)]
pub struct ArenaCommit {
    pub report:EditReport,
    pub created:Vec<ArenaBinding>,
    pub retired:Vec<ArenaBinding>,
    pub remapped:Vec<ArenaRemap>,
}

#[derive(Debug,Clone)]
enum ArenaOperation {
    SplitEdge{edge:ArenaHandle,fraction:f64},
    SplitFace{face:ArenaHandle,start:ArenaHandle,end:ArenaHandle},
    KillEdgeVertex{vertex:ArenaHandle},
    KillEdgeFace{edge:ArenaHandle,removed_face:ArenaHandle},
}

/// Staged modeling edits. Dropping without commit has no effect. A committed
/// batch is one CommandHistory revision and one undoable Memento.
pub struct ArenaTransaction<'a> {
    arena:&'a mut BrepArena,
    ops:Vec<ArenaOperation>,
}
impl ArenaTransaction<'_> {
    fn checked(&self,handle:ArenaHandle,kind:&str)->Result<(),GeometryError>{
        let found=self.arena.resolve(handle).map_err(|e|invalid(e.to_string()))?.entity;
        let actual=match found {
            TopologyEntity::Vertex(_)=>"vertex",TopologyEntity::Edge(_)=>"edge",
            TopologyEntity::Face(_)=>"face",_=>"other",
        };
        if actual!=kind {return Err(invalid(format!("expected {kind}, found {actual}")));}
        Ok(())
    }
    pub fn split_edge(&mut self,edge:ArenaHandle,fraction:f64)->Result<(),GeometryError>{
        self.checked(edge,"edge")?;
        if !fraction.is_finite()||fraction<=0.0||fraction>=1.0 {
            return Err(invalid("split fraction must be strictly interior and finite"));
        }
        self.ops.push(ArenaOperation::SplitEdge{edge,fraction});Ok(())
    }
    pub fn split_face(&mut self,face:ArenaHandle,start:ArenaHandle,end:ArenaHandle)
        ->Result<(),GeometryError>{
        self.checked(face,"face")?;self.checked(start,"vertex")?;self.checked(end,"vertex")?;
        self.ops.push(ArenaOperation::SplitFace{face,start,end});Ok(())
    }
    pub fn kill_edge_vertex(&mut self,vertex:ArenaHandle)->Result<(),GeometryError>{
        self.checked(vertex,"vertex")?;
        self.ops.push(ArenaOperation::KillEdgeVertex{vertex});Ok(())
    }
    pub fn kill_edge_face(&mut self,edge:ArenaHandle,removed_face:ArenaHandle)
        ->Result<(),GeometryError>{
        self.checked(edge,"edge")?;self.checked(removed_face,"face")?;
        self.ops.push(ArenaOperation::KillEdgeFace{edge,removed_face});Ok(())
    }
    pub fn staged_len(&self)->usize{self.ops.len()}
    pub fn commit(self)->Result<ArenaCommit,GeometryError>{self.arena.commit_ops(self.ops)}
}
