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
    BrepCarrierEditor, BrepStructuralEditor,
    ThroughHoleBrep, HoleParameters, BlindHoleBrep, BlindHoleParameters,
    IntersectingBlindBrep, HistoryDirection,
    BooleanError, BooleanKernel, BooleanOperation, BooleanResult,
    EditDelta, EditReport, EulerCommand, GeometryError, HistoryTransition, ModelUnits,
    Solid, TopologyDeltaPlan, TopologyEntity, TopologyName,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ARENA: AtomicU64 = AtomicU64::new(1);
static NEXT_TRANSACTION: AtomicU64 = AtomicU64::new(1);

fn invalid(msg:impl Into<String>)->GeometryError {GeometryError::InvalidEdit(msg.into())}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ArenaHandle {
    pub arena_token:u64,
    pub slot:u32,
    pub generation:u64,
}
/// Identifies a staged operation, scoped to one transaction instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ArenaStep {
    transaction_token: u64,
    index: usize,
}
impl ArenaStep {
    pub fn index(self) -> usize { self.index }
    pub fn output(self, role: ArenaOutput) -> ArenaEntityRef {
        ArenaEntityRef::Created { step: self, role }
    }
}

/// Only outputs with a meaningful role in the current Euler operators can
/// become operands of later commands.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArenaOutput {
    InsertedVertex,
    RemainderEdge,
    CreatedFace,
    DiagonalEdge,
}

/// Existing handles are resolved by semantic name; created references are
/// transaction-local and have no persistent arena handle before commit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArenaEntityRef {
    Existing(ArenaHandle),
    Created { step: ArenaStep, role: ArenaOutput },
}
impl From<ArenaHandle> for ArenaEntityRef {
    fn from(value: ArenaHandle) -> Self { Self::Existing(value) }
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
    AnalyticHole(ThroughHoleBrep),
    AnalyticBlind(BlindHoleBrep),
    AnalyticCrossing(IntersectingBlindBrep),
}
impl ArenaKernel {
    fn revision(&self)->u64{match self{Self::Editable(h)=>h.solid().revision,Self::CurvedReadOnly(_)=>0,
        Self::AnalyticHole(h)=>h.revision(),Self::AnalyticBlind(h)=>h.revision(),
        Self::AnalyticCrossing(h)=>h.revision()}}
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        match self {
            Self::Editable(h)=>BrepModel::from_polyhedral(h.solid(),ModelUnits::default()),
            Self::CurvedReadOnly(body)=>BrepModel::from_cylinder(body),
            Self::AnalyticHole(body)=>Ok(body.model().clone()),
            Self::AnalyticBlind(body)=>Ok(body.model().clone()),
            Self::AnalyticCrossing(body)=>Ok(body.model().clone()),
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
                Self::CurvedReadOnly(_) | Self::AnalyticHole(_) | Self::AnalyticBlind(_) |
                    Self::AnalyticCrossing(_)=>{
                    let (ty,idx)=match e {
                        TopologyEntity::Body=>("body",0),
                        TopologyEntity::Vertex(x)=>("vertex",x.0),
                        TopologyEntity::Edge(x)=>("edge",x.0),
                        TopologyEntity::Coedge(x)=>("coedge",x.0),
                        TopologyEntity::Loop(x)=>("loop",x.0),
                        TopologyEntity::Face(x)=>("face",x.0),
                        TopologyEntity::Shell(x)=>("shell",x.0),
                    };
                    let prefix=match self {Self::AnalyticHole(_)=>"analytic_through_hole",
                        Self::AnalyticBlind(_)=>"analytic_blind_hole",
                        Self::AnalyticCrossing(_)=>"analytic_intersecting_blind_holes",
                        _=>"analytic_cylinder"};
                    TopologyName(format!("{prefix}/{ty}/{idx}"))
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
            },
            BrepBody::ThroughHole(body)=>{
                body.validate()?;ArenaKernel::AnalyticHole(body)
            },
            BrepBody::BlindHole(body)=>{
                body.validate()?;ArenaKernel::AnalyticBlind(body)
            },
            BrepBody::IntersectingBlind(body)=>{
                body.validate()?;ArenaKernel::AnalyticCrossing(body)
            }
        };
        let mut registry=Registry::empty();
        registry.reconcile(&storage.named()?)?;
        Ok(Self{storage,registry})
    }
    pub fn from_solid(body:Solid)->Result<Self,GeometryError>{Self::new(BrepBody::Polyhedral(body))}
    pub fn from_cylinder(body:CylindricalBrep)->Result<Self,GeometryError>{Self::new(BrepBody::Cylinder(body))}
    pub fn from_through_hole(body:ThroughHoleBrep)->Result<Self,GeometryError>{Self::new(BrepBody::ThroughHole(body))}
    pub fn from_blind_hole(body:BlindHoleBrep)->Result<Self,GeometryError>{Self::new(BrepBody::BlindHole(body))}
    pub fn from_intersecting_blind(body:IntersectingBlindBrep)->Result<Self,GeometryError>{
        Self::new(BrepBody::IntersectingBlind(body))
    }
    pub fn arena_token(&self)->u64{self.registry.token}
    pub fn revision(&self)->u64{self.storage.revision()}
    pub fn is_editable(&self)->bool{matches!(&self.storage,ArenaKernel::Editable(_))}
    pub fn model(&self)->Result<BrepModel,GeometryError>{self.storage.brep_model()}
    /// Create a detached, independently revisioned canonical geometry editor
    /// from the current validated arena state. The source arena and its Euler
    /// or parametric feature history are unchanged; importing physical-shape
    /// changes back into an existing arena is intentionally not implemented.
    pub fn detached_carrier_editor(&self)->Result<BrepCarrierEditor,GeometryError>{
        BrepCarrierEditor::new(self.model()?)
    }
    /// Start a separate structural + carrier editor from a validated arena.
    /// Its hole features can change genus and analytic mass, unlike the
    /// shape-preserving detached carrier editor. Source arena history and
    /// existing feature handles remain untouched.
    pub fn detached_structural_editor(&self)->Result<BrepStructuralEditor,GeometryError>{
        match &self.storage {
            ArenaKernel::Editable(history)=>
                BrepStructuralEditor::new(history.solid().clone()),
            ArenaKernel::AnalyticHole(hole)=>
                BrepStructuralEditor::from_through_hole(hole),
            ArenaKernel::CurvedReadOnly(_) | ArenaKernel::AnalyticBlind(_) |
                ArenaKernel::AnalyticCrossing(_)=>
                Err(invalid("only box and full-height through-hole are supported structural sources")),
        }
    }
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
    pub fn undo_depth(&self)->usize{match &self.storage {ArenaKernel::Editable(h)=>h.undo_depth(),
        ArenaKernel::AnalyticHole(h)=>h.undo_depth(),
        ArenaKernel::AnalyticBlind(h)=>h.undo_depth(),
        ArenaKernel::AnalyticCrossing(h)=>h.undo_depth(),_=>0}}
    pub fn redo_depth(&self)->usize{match &self.storage {ArenaKernel::Editable(h)=>h.redo_depth(),
        ArenaKernel::AnalyticHole(h)=>h.redo_depth(),
        ArenaKernel::AnalyticBlind(h)=>h.redo_depth(),
        ArenaKernel::AnalyticCrossing(h)=>h.redo_depth(),_=>0}}
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
    /// Non-destructive, exact planar box Boolean. The resulting Solid, when
    /// present, is an independent B-rep with a fresh identity incarnation;
    /// operands and their revision/Undo stacks stay completely unchanged.
    /// General polyhedral or analytic curved Booleans are not yet executable.
    pub fn boolean(&self,other:&BrepArena,op:BooleanOperation)
        ->Result<BooleanResult,BooleanError>
    {
        match (&self.storage,&other.storage) {
            (ArenaKernel::Editable(a),ArenaKernel::Editable(b)) =>
                BooleanKernel::apply_solids(a.solid(),b.solid(),op),
            (ArenaKernel::Editable(a),ArenaKernel::CurvedReadOnly(b)) =>
                BooleanKernel::apply(&BrepBody::Polyhedral(a.solid().clone()),
                    &BrepBody::Cylinder(b.clone()),op),
            _ => Err(BooleanError::UnsupportedOperand(
                "analytical curved B-rep Boolean requires mutable curved topology")),
        }
    }
    /// Strong-atomic parametric curved edit. A fresh validated B-rep and a
    /// reconciled generational registry publish together on success.
    pub fn update_hole(&mut self,parameters:HoleParameters)->Result<u64,GeometryError>{
        let ArenaKernel::AnalyticHole(current)=&self.storage else {
            return Err(invalid("arena does not own a parametrically editable hole"));
        };
        let mut candidate=current.clone();
        let mut tx=candidate.begin_edit();
        tx.set_center(parameters.center_x,parameters.center_z)?;
        tx.set_radius(parameters.radius)?;
        let revision=tx.commit()?;
        let mut registry=self.registry.clone();
        registry.reconcile(&ArenaKernel::AnalyticHole(candidate.clone()).named()?)?;
        self.storage=ArenaKernel::AnalyticHole(candidate);
        self.registry=registry;
        debug_assert!(self.validate().is_ok());
        Ok(revision)
    }
    /// Atomic replacement of a certified top-entry blind bore's true
    /// analytic depth/center/radius, with entire registry reconciliation
    /// completed before publishing the replacement geometry.
    pub fn update_blind_hole(&mut self,parameters:BlindHoleParameters)
        ->Result<u64,GeometryError>
    {
        let ArenaKernel::AnalyticBlind(current)=&self.storage else {
            return Err(invalid("arena does not own an editable blind bore"));
        };
        let mut candidate=current.clone();
        let mut tx=candidate.begin_edit();
        tx.set_center(parameters.center_x,parameters.center_z)?;
        tx.set_radius(parameters.radius)?;
        tx.set_depth(parameters.depth)?;
        let revision=tx.commit()?;
        let mut registry=self.registry.clone();
        registry.reconcile(&ArenaKernel::AnalyticBlind(candidate.clone()).named()?)?;
        self.storage=ArenaKernel::AnalyticBlind(candidate);
        self.registry=registry;
        debug_assert!(self.validate().is_ok());
        Ok(revision)
    }
    /// Rebuild an intersecting-circle analytic B-rep and reconcile its
    /// generational topology registry before publishing either object.
    pub fn update_intersecting_blind(&mut self,
        params:[crate::BlindHoleParameters;2])->Result<u64,GeometryError>
    {
        let ArenaKernel::AnalyticCrossing(current)=&self.storage else{
            return Err(invalid("arena does not own a crossed-blind analytic feature"));
        };
        let mut candidate=current.clone();
        let mut tx=candidate.begin_edit();
        tx.replace(0,params[0])?;
        tx.replace(1,params[1])?;
        let revision=tx.commit()?;
        let mut registry=self.registry.clone();
        registry.reconcile(&ArenaKernel::AnalyticCrossing(candidate.clone()).named()?)?;
        self.storage=ArenaKernel::AnalyticCrossing(candidate);
        self.registry=registry;
        debug_assert!(self.validate().is_ok());
        Ok(revision)
    }
    pub fn begin_transaction(&mut self)->Result<ArenaTransaction<'_>,GeometryError>{
        if !self.is_editable(){return Err(invalid("analytical curved B-rep is read-only in v0.17"));}
        let token = NEXT_TRANSACTION.fetch_add(1, Ordering::Relaxed);
        if token == 0 || token == u64::MAX {
            return Err(invalid("arena transaction token generator exhausted"));
        }
        Ok(ArenaTransaction{arena:self,token,ops:Vec::new(),failed:false})
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
    fn commit_ops(&mut self, ops: Vec<ArenaOperation>, transaction_token: u64)
        -> Result<ArenaCommit, GeometryError>
    {
        if ops.is_empty() || ops.len() > 256 {
            return Err(invalid("arena batch must have between 1 and 256 operations"));
        }
        let ArenaKernel::Editable(current) = &self.storage else {
            return Err(invalid("arena is read-only"));
        };
        let before = self.registry.clone();
        // Never expose partially compiled, edited or reconciled state.
        let mut candidate = current.clone();
        let needs_preview = ops.iter().any(ArenaOperation::uses_created);
        // Only dependent batches require a second private Solid. Use actual
        // EditDelta outputs; never infer newly minted numeric IDs or names.
        let mut preview = needs_preview.then(|| current.solid().clone());
        let mut produced = BTreeMap::<(usize, ArenaOutput), TopologyName>::new();
        let mut commands = Vec::with_capacity(ops.len());
        for (index, operation) in ops.iter().enumerate() {
            let command = compile_operation(
                &before, operation, transaction_token, index,
                preview.as_ref(), &produced,
            )?;
            if let Some(body) = preview.as_mut() {
                let (deltas, _) = command.apply_to_candidate(body)?;
                if deltas.len() != 1 {
                    return Err(invalid("Euler preview produced an unexpected delta count"));
                }
                for (role, entity) in emitted_entities(&deltas[0]) {
                    let name = body.topology_name(entity)
                        .ok_or_else(||invalid("new Euler output has no semantic name"))?.clone();
                    if produced.insert((index, role), name).is_some() {
                        return Err(invalid("duplicate Euler step output"));
                    }
                }
            }
            commands.push(command);
        }
        // One replayable batch -> one revision and a single undo/redo step.
        let report = candidate.execute(CommandBatch{commands})?;
        if let Some(preview) = preview {
            // Divergence in identity liveness between preflight and committed
            // batch must abort before publication.
            for name in produced.values() {
                let compiled_live = preview.resolve_topology_name(name).is_some();
                let replay_live = candidate.solid().resolve_topology_name(name).is_some();
                if compiled_live != replay_live {
                    return Err(invalid("dependent Euler preview diverged from atomic replay"));
                }
            }
        }
        let storage = ArenaKernel::Editable(candidate);
        let mut registry = before.clone();
        registry.reconcile(&storage.named()?)?;
        let summary = commit_summary(&before, &registry, report);
        self.storage = storage;
        self.registry = registry;
        debug_assert!(self.validate().is_ok());
        Ok(summary)
    }
    pub fn undo(&mut self)->Result<Option<HistoryTransition>,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<Option<HistoryTransition>,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<Option<HistoryTransition>,GeometryError>{
        if let ArenaKernel::AnalyticCrossing(current)=&self.storage {
            let mut candidate=current.clone();
            let before=candidate.revision();
            let changed=if redo {candidate.redo()?}else{candidate.undo()?};
            if !changed{return Ok(None);}
            let after=candidate.revision();
            let mut registry=self.registry.clone();
            registry.reconcile(&ArenaKernel::AnalyticCrossing(candidate.clone()).named()?)?;
            self.storage=ArenaKernel::AnalyticCrossing(candidate);
            self.registry=registry;
            debug_assert!(self.validate().is_ok());
            return Ok(Some(HistoryTransition{
                direction:if redo{HistoryDirection::Redo}else{HistoryDirection::Undo},
                revision_before:before,revision_after:after,command_count:1
            }));
        }
        if let ArenaKernel::AnalyticBlind(current)=&self.storage {
            let mut candidate=current.clone();
            let before=candidate.revision();
            let changed=if redo {candidate.redo()?}else{candidate.undo()?};
            if !changed {return Ok(None);}
            let after=candidate.revision();
            let mut registry=self.registry.clone();
            registry.reconcile(&ArenaKernel::AnalyticBlind(candidate.clone()).named()?)?;
            self.storage=ArenaKernel::AnalyticBlind(candidate);
            self.registry=registry;
            debug_assert!(self.validate().is_ok());
            return Ok(Some(HistoryTransition{direction:if redo{HistoryDirection::Redo}
                else{HistoryDirection::Undo},revision_before:before,
                revision_after:after,command_count:1}));
        }
        if let ArenaKernel::AnalyticHole(current)=&self.storage {
            let mut candidate=current.clone();
            let before=candidate.revision();
            let changed=if redo {candidate.redo()?}else{candidate.undo()?};
            if !changed {return Ok(None)};
            let after=candidate.revision();
            let mut registry=self.registry.clone();
            registry.reconcile(&ArenaKernel::AnalyticHole(candidate.clone()).named()?)?;
            self.storage=ArenaKernel::AnalyticHole(candidate);
            self.registry=registry;
            debug_assert!(self.validate().is_ok());
            return Ok(Some(HistoryTransition {direction:if redo{HistoryDirection::Redo}
                else{HistoryDirection::Undo},revision_before:before,revision_after:after,
                command_count:1}));
        }
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

/// Classifies topology without relying on the current dense numeric ID.
fn entity_kind(entity: TopologyEntity) -> &'static str {
    match entity {
        TopologyEntity::Vertex(_) => "vertex",
        TopologyEntity::Edge(_) => "edge",
        TopologyEntity::Face(_) => "face",
        _ => "other",
    }
}

fn resolve_operand(
    registry: &Registry,
    operand: ArenaEntityRef,
    expected_kind: &str,
    transaction_token: u64,
    before_step: usize,
    preview: Option<&Solid>,
    produced: &BTreeMap<(usize, ArenaOutput), TopologyName>,
) -> Result<TopologyName, GeometryError> {
    let name = match operand {
        ArenaEntityRef::Existing(handle) => {
            let binding = registry.binding(handle).map_err(|e|invalid(e.to_string()))?;
            if entity_kind(binding.entity) != expected_kind {
                return Err(invalid(format!("expected {expected_kind}, found {}",entity_kind(binding.entity))));
            }
            binding.name
        }
        ArenaEntityRef::Created { step, role } => {
            if step.transaction_token != transaction_token || step.index >= before_step {
                return Err(invalid("foreign or forward Euler output reference"));
            }
            produced.get(&(step.index, role)).cloned()
                .ok_or_else(|| invalid("Euler step output was not produced"))?
        }
    };
    if let Some(body) = preview {
        // Kills may compact IDs or retire an earlier output. A name must
        // still resolve in the state at the point of its use.
        let handle = body.resolve_topology_name(&name)
            .ok_or_else(||invalid(format!("Euler operand retired before use: {name}")))?;
        let entity = body.resolve_topology_handle(handle)
            .map_err(|error|invalid(error.to_string()))?;
        if entity_kind(entity) != expected_kind {
            return Err(invalid(format!("expected {expected_kind}, found {}",entity_kind(entity))));
        }
    }
    Ok(name)
}

fn compile_operation(
    registry: &Registry,
    operation: &ArenaOperation,
    transaction_token: u64,
    index: usize,
    preview: Option<&Solid>,
    produced: &BTreeMap<(usize, ArenaOutput), TopologyName>,
) -> Result<EulerCommand, GeometryError> {
    let name = |reference, expected| resolve_operand(
        registry, reference, expected, transaction_token, index, preview, produced);
    Ok(match *operation {
        ArenaOperation::SplitEdge { edge, fraction } => EulerCommand::MakeEdgeVertex {
            edge: name(edge, "edge")?, fraction,
        },
        ArenaOperation::SplitFace { face, start, end } => EulerCommand::MakeEdgeFace {
            face: name(face, "face")?,
            start: name(start, "vertex")?,
            end: name(end, "vertex")?,
        },
        ArenaOperation::KillEdgeVertex { vertex } => EulerCommand::KillEdgeVertex {
            vertex: name(vertex, "vertex")?,
        },
        ArenaOperation::KillEdgeFace { edge, removed_face } => EulerCommand::KillEdgeFace {
            edge: name(edge, "edge")?,
            removed_face: name(removed_face, "face")?,
        },
    })
}

fn emitted_entities(delta: &EditDelta) -> Vec<(ArenaOutput, TopologyEntity)> {
    match delta {
        EditDelta::SplitEdge(split) => vec![
            (ArenaOutput::InsertedVertex, TopologyEntity::Vertex(split.inserted_vertex)),
            (ArenaOutput::RemainderEdge, TopologyEntity::Edge(split.created_edge)),
        ],
        EditDelta::SplitFace(split) => vec![
            (ArenaOutput::CreatedFace, TopologyEntity::Face(split.created_face)),
            (ArenaOutput::DiagonalEdge, TopologyEntity::Edge(split.diagonal_edge)),
        ],
        EditDelta::KillEdgeVertex(_) | EditDelta::KillEdgeFace(_) => Vec::new(),
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

#[derive(Debug, Clone, Copy)]
enum ArenaOperation {
    SplitEdge { edge: ArenaEntityRef, fraction: f64 },
    SplitFace { face: ArenaEntityRef, start: ArenaEntityRef, end: ArenaEntityRef },
    KillEdgeVertex { vertex: ArenaEntityRef },
    KillEdgeFace { edge: ArenaEntityRef, removed_face: ArenaEntityRef },
}
impl ArenaOperation {
    fn output_kind(&self, output: ArenaOutput) -> Option<&'static str> {
        match (self, output) {
            (Self::SplitEdge { .. }, ArenaOutput::InsertedVertex) => Some("vertex"),
            (Self::SplitEdge { .. }, ArenaOutput::RemainderEdge) => Some("edge"),
            (Self::SplitFace { .. }, ArenaOutput::CreatedFace) => Some("face"),
            (Self::SplitFace { .. }, ArenaOutput::DiagonalEdge) => Some("edge"),
            _ => None,
        }
    }
    fn uses_created(&self) -> bool {
        let created = |r: ArenaEntityRef| matches!(r, ArenaEntityRef::Created { .. });
        match *self {
            Self::SplitEdge { edge, .. } => created(edge),
            Self::SplitFace { face, start, end } =>
                created(face) || created(start) || created(end),
            Self::KillEdgeVertex { vertex } => created(vertex),
            Self::KillEdgeFace { edge, removed_face } =>
                created(edge) || created(removed_face),
        }
    }
}

/// Staged modeling edits. A placeholder may refer only to an output from a
/// strictly earlier stage in the same transaction. Neither staging nor failed
/// dependency preflight changes the arena. A successful batch is one Memento.
pub struct ArenaTransaction<'a> {
    arena: &'a mut BrepArena,
    token: u64,
    ops: Vec<ArenaOperation>,
    // Even if the caller ignores Err, a failed stage makes commit impossible.
    failed: bool,
}
impl ArenaTransaction<'_> {
    fn checked(&mut self, reference: ArenaEntityRef, expected_kind: &str)
        -> Result<(), GeometryError>
    {
        if self.failed { return Err(invalid("arena transaction already failed")); }
        let result = match reference {
            ArenaEntityRef::Existing(handle) => self.arena.resolve(handle)
                .map_err(|error|invalid(error.to_string()))
                .and_then(|binding| {
                    let actual = entity_kind(binding.entity);
                    if actual == expected_kind { Ok(()) }
                    else { Err(invalid(format!("expected {expected_kind}, found {actual}"))) }
                }),
            ArenaEntityRef::Created { step, role } => {
                if step.transaction_token != self.token {
                    Err(invalid("output belongs to another arena transaction"))
                } else if step.index >= self.ops.len() {
                    Err(invalid("Euler output must come from an earlier stage"))
                } else {
                    match self.ops[step.index].output_kind(role) {
                        Some(actual) if actual == expected_kind => Ok(()),
                        Some(actual) => Err(invalid(format!(
                            "expected {expected_kind}, found produced {actual}"))),
                        None => Err(invalid("operation does not produce the requested Euler output")),
                    }
                }
            }
        };
        if result.is_err() { self.failed = true; }
        result
    }
    fn stage(&mut self, operation: ArenaOperation) -> Result<ArenaStep, GeometryError> {
        if self.failed { return Err(invalid("arena transaction already failed")); }
        if self.ops.len() >= 256 {
            self.failed = true;
            return Err(invalid("arena transaction cannot exceed 256 commands"));
        }
        let step = ArenaStep { transaction_token: self.token, index: self.ops.len() };
        self.ops.push(operation);
        Ok(step)
    }

    /// Backwards-compatible handle-only API; see split_edge_ref for outputs.
    pub fn split_edge(&mut self, edge: ArenaHandle, fraction: f64)
        -> Result<(), GeometryError>
    {
        self.split_edge_ref(edge, fraction).map(|_| ())
    }
    /// The returned step can supply InsertedVertex and RemainderEdge operands.
    pub fn split_edge_ref(&mut self, edge: impl Into<ArenaEntityRef>, fraction: f64)
        -> Result<ArenaStep, GeometryError>
    {
        let edge = edge.into();
        self.checked(edge, "edge")?;
        if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 {
            self.failed = true;
            return Err(invalid("split fraction must be strictly interior and finite"));
        }
        self.stage(ArenaOperation::SplitEdge { edge, fraction })
    }
    pub fn split_face(&mut self, face: ArenaHandle, start: ArenaHandle, end: ArenaHandle)
        -> Result<(), GeometryError>
    {
        self.split_face_ref(face, start, end).map(|_| ())
    }
    /// The returned step can supply CreatedFace and DiagonalEdge operands.
    pub fn split_face_ref(
        &mut self,
        face: impl Into<ArenaEntityRef>,
        start: impl Into<ArenaEntityRef>,
        end: impl Into<ArenaEntityRef>,
    ) -> Result<ArenaStep, GeometryError> {
        let (face, start, end) = (face.into(), start.into(), end.into());
        self.checked(face, "face")?;
        self.checked(start, "vertex")?;
        self.checked(end, "vertex")?;
        self.stage(ArenaOperation::SplitFace { face, start, end })
    }
    pub fn kill_edge_vertex(&mut self, vertex: ArenaHandle)
        -> Result<(), GeometryError>
    {
        self.kill_edge_vertex_ref(vertex).map(|_| ())
    }
    pub fn kill_edge_vertex_ref(&mut self, vertex: impl Into<ArenaEntityRef>)
        -> Result<ArenaStep, GeometryError>
    {
        let vertex = vertex.into();
        self.checked(vertex, "vertex")?;
        self.stage(ArenaOperation::KillEdgeVertex { vertex })
    }
    pub fn kill_edge_face(&mut self, edge: ArenaHandle, removed_face: ArenaHandle)
        -> Result<(), GeometryError>
    {
        self.kill_edge_face_ref(edge, removed_face).map(|_| ())
    }
    pub fn kill_edge_face_ref(
        &mut self, edge: impl Into<ArenaEntityRef>,
        removed_face: impl Into<ArenaEntityRef>,
    ) -> Result<ArenaStep, GeometryError> {
        let (edge, removed_face) = (edge.into(), removed_face.into());
        self.checked(edge, "edge")?;
        self.checked(removed_face, "face")?;
        self.stage(ArenaOperation::KillEdgeFace { edge, removed_face })
    }
    pub fn staged_len(&self) -> usize { self.ops.len() }
    pub fn commit(self) -> Result<ArenaCommit, GeometryError> {
        if self.failed {
            return Err(invalid("cannot commit an arena transaction after a rejected stage"));
        }
        self.arena.commit_ops(self.ops, self.token)
    }
}
