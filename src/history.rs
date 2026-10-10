//! Command + Memento orchestration for the supported Euler-preserving edits.
//!
//! Commands keep the ordered, replayable *intent* using feature-relative names.
//! An opaque Memento keeps local inverse frames from the original transaction.
//! Undo restores those frames, deleting the inserted edge/vertex or edge/face
//! pair; redo replays the command to mint new runtime handle generations.
//! KEV/KEF are constrained inverse topology operators: only collinear
//! valence-two vertices and adjacent coplanar faces are removable.

use crate::{
    EditDelta, EditReport, EditTransaction, GeometryError, JournalStats, NameChange,
    Solid, TopologyEntity, TopologyName,
};
use crate::edit::EditMemento;
use serde::{Deserialize, Serialize};

const DEFAULT_UNDO_LIMIT: usize = 128;
const MAX_COMMANDS_PER_BATCH: usize = 256;

fn command_error(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidEdit(message.into())
}

/// Reproducible modeling intent. Names, rather than process-local handles,
/// survive serialization and replay from the same feature-construction history.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum EulerCommand {
    /// Adds one vertex and one edge. Undo removes exactly those additions.
    SplitEdge { edge: TopologyName, fraction: f64 },
    /// Adds one edge and one face. Undo removes exactly those additions.
    SplitFace {
        face: TopologyName,
        start: TopologyName,
        end: TopologyName,
    },
    /// MEV-compatible edge subdivision. Alias of the historical SplitEdge.
    MakeEdgeVertex { edge: TopologyName, fraction: f64 },
    /// MEF-compatible face division. Alias of the historical SplitFace.
    MakeEdgeFace { face: TopologyName, start: TopologyName, end: TopologyName },
    /// KEV for a degree-two collinear vertex; no Undo-stack dependency.
    KillEdgeVertex { vertex: TopologyName },
    /// KEF for two incident coplanar polygonal faces.
    KillEdgeFace { edge: TopologyName, removed_face: TopologyName },
}

impl EulerCommand {
    fn resolve(tx: &EditTransaction<'_>, name: &TopologyName) -> Result<TopologyEntity, GeometryError> {
        let handle = tx.preview().resolve_topology_name(name)
            .ok_or_else(|| command_error(format!("unresolved topology name: {name}")))?;
        tx.preview().resolve_topology_handle(handle)
            .map_err(|error| command_error(error.to_string()))
    }

    fn resolve_body(solid:&Solid,name:&TopologyName)->Result<TopologyEntity,GeometryError>{
        let handle=solid.resolve_topology_name(name)
            .ok_or_else(||command_error(format!("unresolved topology name: {name}")))?;
        solid.resolve_topology_handle(handle).map_err(|e|command_error(e.to_string()))
    }

    fn is_kill(&self)->bool {
        matches!(self,Self::KillEdgeVertex{..}|Self::KillEdgeFace{..})
    }

    /// Called on a disposable candidate body when the batch requires
    /// compacting/renumbering Euler kills. Make commands still use local edits.
    fn apply_to_candidate(&self,solid:&mut Solid)
        ->Result<(Vec<EditDelta>,Vec<NameChange>),GeometryError>
    {
        match self {
            Self::KillEdgeVertex {vertex} => {
                let id=match Self::resolve_body(solid,vertex)? {
                    TopologyEntity::Vertex(id)=>id,
                    other=>return Err(command_error(format!("KEV expects vertex; got {other:?}"))),
                };
                let killed=solid.kill_edge_vertex(id)?;
                Ok((vec![EditDelta::KillEdgeVertex(killed)],Vec::new()))
            }
            Self::KillEdgeFace {edge,removed_face} => {
                let edge_id=match Self::resolve_body(solid,edge)? {
                    TopologyEntity::Edge(id)=>id,
                    other=>return Err(command_error(format!("KEF expects edge; got {other:?}"))),
                };
                let face_id=match Self::resolve_body(solid,removed_face)? {
                    TopologyEntity::Face(id)=>id,
                    other=>return Err(command_error(format!("KEF expects face; got {other:?}"))),
                };
                let killed=solid.kill_edge_face(edge_id,face_id)?;
                Ok((vec![EditDelta::KillEdgeFace(killed)],Vec::new()))
            }
            _=>{
                let mut tx=solid.begin_edit()?;
                self.apply(&mut tx)?;
                let report=tx.commit()?;
                Ok((report.changes,report.named_changes))
            }
        }
    }

    fn apply(&self, tx: &mut EditTransaction<'_>) -> Result<(), GeometryError> {
        match self {
            Self::SplitEdge { edge, fraction } | Self::MakeEdgeVertex { edge, fraction } => {
                let id = match Self::resolve(tx, edge)? {
                    TopologyEntity::Edge(id) => id,
                    other => return Err(command_error(format!("expected edge, found {other:?}"))),
                };
                tx.split_edge(id, *fraction)?;
            }
            Self::SplitFace { face, start, end } | Self::MakeEdgeFace { face, start, end } => {
                let f = match Self::resolve(tx, face)? {
                    TopologyEntity::Face(id) => id,
                    other => return Err(command_error(format!("expected face, found {other:?}"))),
                };
                let a = match Self::resolve(tx, start)? {
                    TopologyEntity::Vertex(id) => id,
                    other => return Err(command_error(format!("expected start vertex, found {other:?}"))),
                };
                let b = match Self::resolve(tx, end)? {
                    TopologyEntity::Vertex(id) => id,
                    other => return Err(command_error(format!("expected end vertex, found {other:?}"))),
                };
                tx.split_face(f, a, b)?;
            }
            Self::KillEdgeVertex { .. } | Self::KillEdgeFace { .. } => {
                return Err(command_error("kill operator requires snapshot-backed Euler transaction"));
            }
        }
        Ok(())
    }
}

/// A user-visible modeling step. Its commands commit as one revision and one
/// Memento, so Undo/Redo always respects batch boundaries.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CommandBatch {
    pub commands: Vec<EulerCommand>,
}

impl CommandBatch {
    pub fn single(command: EulerCommand) -> Self { Self { commands: vec![command] } }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistoryDirection { Undo, Redo }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HistoryTransition {
    pub direction: HistoryDirection,
    pub revision_before: u64,
    pub revision_after: u64,
    pub command_count: usize,
}

struct AppliedBatch {
    batch: CommandBatch,
    memento: EditMemento,
}

/// Owns the working body and enforces single-writer command history.
/// `archive` retains the complete applied command log when old inverse
/// Mementos are pruned to bound memory; archived steps cannot be undone.
pub struct CommandHistory {
    solid: Solid,
    archive: Vec<CommandBatch>,
    undo: Vec<AppliedBatch>,
    redo: Vec<CommandBatch>,
    undo_limit: usize,
}

impl CommandHistory {
    pub fn new(solid: Solid) -> Result<Self, GeometryError> {
        Self::with_undo_limit(solid, DEFAULT_UNDO_LIMIT)
    }

    pub fn with_undo_limit(solid: Solid, undo_limit: usize) -> Result<Self, GeometryError> {
        if undo_limit == 0 {
            return Err(command_error("undo limit must be positive"));
        }
        solid.validate()?;
        Ok(Self { solid, archive: Vec::new(), undo: Vec::new(),
            redo: Vec::new(), undo_limit })
    }

    pub fn solid(&self) -> &Solid { &self.solid }
    pub fn into_solid(self) -> Solid { self.solid }
    pub fn undo_depth(&self) -> usize { self.undo.len() }
    pub fn redo_depth(&self) -> usize { self.redo.len() }
    pub fn archived_depth(&self) -> usize { self.archive.len() }

    /// The applied command sequence, including older non-undoable commands.
    /// The redo branch is deliberately excluded after Undo.
    pub fn command_log(&self) -> Vec<CommandBatch> {
        self.archive.iter().cloned().chain(self.undo.iter().map(|item| item.batch.clone())).collect()
    }

    /// Deterministically reconstruct from a newly generated primitive with
    /// the same feature key and equivalent topology construction order.
    pub fn replay(solid: Solid, batches: &[CommandBatch]) -> Result<Self, GeometryError> {
        let mut history = Self::new(solid)?;
        for batch in batches { history.execute(batch.clone())?; }
        Ok(history)
    }

    /// Command executes and stores its inverse Memento on success only.
    /// Executing on an undone state discards the redo branch.
    pub fn execute(&mut self, batch: CommandBatch) -> Result<EditReport, GeometryError> {
        self.apply(batch, true)
    }

    fn apply(&mut self, batch: CommandBatch, clear_redo: bool) -> Result<EditReport, GeometryError> {
        if batch.commands.is_empty() || batch.commands.len() > MAX_COMMANDS_PER_BATCH {
            return Err(command_error("command batch must have between 1 and 256 commands"));
        }
        let (report,memento)=if batch.commands.iter().any(EulerCommand::is_kill) {
            // Kill operators change compact slot IDs and need a validated,
            // snapshot-backed candidate. Failed batches leave self untouched.
            let before=self.solid.revision;
            let after=before.checked_add(1).ok_or_else(||command_error("revision overflow"))?;
            let mut candidate=self.solid.clone();
            let mut changes=Vec::new();
            let mut named_changes=Vec::new();
            for command in &batch.commands {
                let (d,n)=command.apply_to_candidate(&mut candidate)?;
                changes.extend(d);named_changes.extend(n);
            }
            candidate.revision=after;
            candidate.validate()?;
            let previous=std::mem::replace(&mut self.solid,candidate);
            let stats=JournalStats{
                frames:1,
                topology_snapshots:previous.vertices.len()+previous.edges.len()
                    +previous.coedges.len()+previous.loops.len()+previous.faces.len()+previous.shells.len(),
                triangle_snapshots:previous.mesh.triangles.len(),
            };
            (EditReport{revision_before:before,revision_after:after,changes,
                named_changes,journal:stats},EditMemento::snapshot(previous))
        } else {
            let mut tx = self.solid.begin_edit()?;
            for command in &batch.commands { command.apply(&mut tx)?; }
            tx.commit_recorded()?
        };
        if clear_redo { self.redo.clear(); }
        self.undo.push(AppliedBatch { batch, memento });
        if self.undo.len() > self.undo_limit {
            let oldest = self.undo.remove(0);
            self.archive.push(oldest.batch);
        }
        Ok(report)
    }

    /// Restores an opaque delta-Memento, including its name index and geometry.
    /// An Undo advances the *runtime revision* monotonically; it does not
    /// restore the previous integer revision number.
    pub fn undo(&mut self) -> Result<Option<HistoryTransition>, GeometryError> {
        let Some(latest) = self.undo.last() else { return Ok(None); };
        latest.memento.check_target(&self.solid)?;
        self.solid.validate()?;
        let before = self.solid.revision;
        before.checked_add(1).ok_or_else(|| command_error("revision overflow"))?;
        let applied = self.undo.pop().expect("checked above");
        let count = applied.batch.commands.len();
        applied.memento.restore(&mut self.solid)?;
        self.redo.push(applied.batch);
        Ok(Some(HistoryTransition { direction: HistoryDirection::Undo,
            revision_before: before, revision_after: self.solid.revision, command_count: count }))
    }

    /// Re-executes the command, giving every recreated entity a fresh runtime
    /// generation. A stale handle from before Undo never resurrects.
    pub fn redo(&mut self) -> Result<Option<HistoryTransition>, GeometryError> {
        let Some(batch) = self.redo.last().cloned() else { return Ok(None); };
        let before = self.solid.revision;
        let count = batch.commands.len();
        self.apply(batch, false)?;
        self.redo.pop();
        Ok(Some(HistoryTransition { direction: HistoryDirection::Redo,
            revision_before: before, revision_after: self.solid.revision, command_count: count }))
    }
}

// These operations preserve Euler characteristic by adding paired elements:
// split_edge: ΔV=+1, ΔE=+1; split_face: ΔE=+1, ΔF=+1.
// Their inverse Memento removes those exact paired additions.
// KEV and KEF can also be applied independently to compatible earlier edits;
// arbitrary topological Euler pairs and trimmed-surface Booleans remain future work.
