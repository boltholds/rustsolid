//! G4b1: source-certified cross-face radial incidence for G4a patches.
//!
//! The two sides of a G2 cut on ONE source face are not physical twins.
//! This module pairs a selected G4a use with a different source face,
//! independently validates both input patches, and preserves analytic
//! source-curve ownership. Incomplete source-face coverage is reported rather
//! than inventing a closed Boolean B-rep. No persistent IDs are allocated.

use crate::{
    BrepBody, BrepModel, BooleanOperation, CoedgeId, EdgeId,
    FaceId, FacePatchAssemblyPlan, FacePairGraph, GeometryError,
    GeometryTolerance, GraphArc, GraphNodeKind, GraphOperand,
    IntersectionCurve, MultiUvBoundaryProvenance, MultiUvBoundarySpan,
    Point3, Surface3,
};
use std::collections::{BTreeMap, BTreeSet};

const MAX_STITCH_PATCHES: usize = 512;
const MAX_STITCH_USES: usize = 32_768;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn distance(a: Point3, b: Point3) -> f64 {
    let d = a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}
fn within(a: Point3, b: Point3, eps: f64) -> bool {
    a.is_finite() && b.is_finite() && distance(a,b) <= eps
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StitchFaceOwner {
    pub operand: GraphOperand,
    pub face: FaceId,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StitchUseRef {
    pub patch: usize,
    pub use_index: usize,
}
#[derive(Debug, Clone, Copy)]
pub struct StitchPatchInput<'a> {
    pub plan: &'a FacePatchAssemblyPlan,
    pub pairs: &'a [FacePairGraph],
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StitchedEdgeOrigin {
    /// Existing canonical edge, with two separate real source coedges.
    Source { operand: GraphOperand, edge: EdgeId },
    /// Intersection shared by a source A face and source B face.
    /// Local graph arc IDs here are canonical only within that face pair.
    Cut { first_face: FaceId, second_face: FaceId, arc: crate::GraphArcId },
}
#[derive(Debug, Clone, PartialEq)]
pub struct StitchedEdgeDraft {
    /// Provisional edge index only, not a persistent EdgeId.
    pub draft_id: u32,
    pub source: StitchedEdgeOrigin,
    /// Two different source faces. Directed uses must oppose in world XYZ.
    pub radial_uses: [StitchUseRef; 2],
    pub source_faces: [StitchFaceOwner; 2],
    pub start: Point3,
    pub end: Point3,
    /// Source pcurves retain distinct UV charts even for the same 3D edge.
    pub face_pcurves: [MultiUvBoundarySpan; 2],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpenStitchReason {
    MissingSourceCoedgeTwin,
    MissingOppositeCutPatch,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenStitchUse {
    pub use_ref: StitchUseRef,
    pub reason: OpenStitchReason,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlobalStitchHold {
    BudgetExceeded,
    DuplicateSourceFace,
    DamagedPatchLoop,
    NonManifoldIncidence,
    ContradictorySourceCoedge,
    CutCarrierMismatch,
    UncertifiedCutEndpoint,
    UncertifiedIntermediateCutSplit,
    SameDirectionUses,
}
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalEdgeStitchPlan {
    pub operation: BooleanOperation,
    pub edges: Vec<StitchedEdgeDraft>,
    pub open_uses: Vec<OpenStitchUse>,
    /// Source faces without a G4a plan are NEVER assumed discarded.
    pub missing_source_faces: Vec<StitchFaceOwner>,
    pub source_patch_count: usize,
    pub selected_coedge_uses: usize,
    pub verified_loops: usize,
    pub radially_paired: bool,
    /// Even a two-use radial graph is not a certified B-rep shell:
    /// vertex-fan continuity, all face domains and an independent volume
    /// proof must be established before G4 can claim a candidate.
    pub closed_shell_proven: bool,
    pub can_apply: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub enum GlobalEdgeStitchResult {
    Preview(GlobalEdgeStitchPlan),
    Held(GlobalStitchHold),
}

#[derive(Clone, Copy)]
struct OriginalUse {
    id: StitchUseRef,
    coedge: CoedgeId,
}
#[derive(Clone, Copy)]
struct CutUse {
    id: StitchUseRef,
    pair: (FaceId,FaceId),
    arc: GraphArc,
    /// Retains the two actual source-face trim/chart contexts.
    boundary_evidence: [bool; 2],
}

fn endpoint_witnesses(
    graph: &crate::IntersectionGraph,
    arc: GraphArc,
) -> Option<[bool; 2]> {
    if arc.closed || !matches!(arc.curve, IntersectionCurve::Line(_)) {
        return None;
    }
    let start = graph.nodes.get(arc.from.0 as usize)
        .filter(|node| node.id == arc.from && node.kind == GraphNodeKind::Endpoint)?;
    let end = graph.nodes.get(arc.to.0 as usize)
        .filter(|node| node.id == arc.to && node.kind == GraphNodeKind::Endpoint)?;
    // The endpoint must be physically bounded on EACH original face. A
    // chart seam or merely coincident XYZ cannot prove material vertices.
    let both = |node: &crate::GraphNode| {
        node.boundary_sources.iter().any(|w| w.operand == GraphOperand::First) &&
        node.boundary_sources.iter().any(|w| w.operand == GraphOperand::Second)
    };
    Some([both(start), both(end)])
}

fn arc_owns_entire_piece(
    arc: GraphArc, use_: &crate::FacePatchDraftUse, eps: f64,
) -> Result<bool, GeometryError> {
    // An internal G2 line-crossing can divide one 3D GraphArc into several
    // pieces. That requires global vertex reconciliation (G1b) before a
    // physical edge may be minted: refuse it in this bounded G4b slice.
    let first = arc.evaluate(0.0)?;
    let last = arc.evaluate(1.0)?;
    Ok((within(use_.start, first, eps) && within(use_.end, last, eps)) ||
       (within(use_.start, last, eps) && within(use_.end, first, eps)))
}

fn pair_geometry(
    a: &crate::FacePatchDraftUse,
    b: &crate::FacePatchDraftUse,
    eps: f64,
) -> bool {
    // Source-normal and Boolean orientation are reflected in these
    // directed endpoints. Never pair same-direction coedges.
    within(a.start, b.end, eps) && within(a.end, b.start, eps)
}

/// A G4a plan is a source-checked **single face**, with global pairs
/// evaluated afresh against both source operands in the validator.
impl GlobalEdgeStitchPlan {
    pub fn build(
        left: &BrepBody, right: &BrepBody,
        input: &[StitchPatchInput<'_>], tolerance: GeometryTolerance,
    ) -> Result<GlobalEdgeStitchResult, GeometryError> {
        tolerance.validate()?;
        left.validate()?; right.validate()?;
        let a=left.shared()?;
        let b=right.shared()?;
        if input.is_empty() || input.len()>MAX_STITCH_PATCHES {
            return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::BudgetExceeded));
        }
        let mut owners=BTreeSet::new();
        let mut operation=None;
        let mut source=BTreeMap::<(GraphOperand,EdgeId),Vec<OriginalUse>>::new();
        let mut cuts=BTreeMap::<(FaceId,FaceId,crate::GraphArcId),Vec<CutUse>>::new();
        let mut selected=0usize;
        let mut verified_loops=0usize;
        for (patch_index, entry) in input.iter().enumerate() {
            entry.plan.validate_against(left,right,entry.pairs,tolerance)?;
            let owner=StitchFaceOwner {
                operand:entry.plan.source_operand,face:entry.plan.source_face,
            };
            if !owners.insert(owner) {
                return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DuplicateSourceFace));
            }
            if operation.is_some_and(|op| op!=entry.plan.operation) {
                return Err(invalid("G4b source patches disagree on Boolean operation"));
            }
            operation=Some(entry.plan.operation);
            let model=if owner.operand==GraphOperand::First {&a}else{&b};
            let extent=model.bbox.max.sub(model.bbox.min).length_squared().sqrt();
            let eps=tolerance.length_at(extent.max(1.0))*16.0;
            if entry.plan.can_apply {
                return Err(invalid("G4b may not consume a prematurely committed G4a plan"));
            }
            let mut used=BTreeSet::new();
            for cycle in &entry.plan.retained_loops {
                if cycle.uses.is_empty() || !matches!(model.face_surface(owner.face),Some(Surface3::Plane(_))) {
                    return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DamagedPatchLoop));
                }
                for (i,&at) in cycle.uses.iter().enumerate() {
                    let Some(use_)=entry.plan.uses.get(at) else {
                        return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DamagedPatchLoop));
                    };
                    if !use_.kept || !used.insert(at) ||
                        use_.cell_index!=cycle.cell_index ||
                        use_.cycle_index!=cycle.cycle_index
                    {
                        return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DamagedPatchLoop));
                    }
                    let next=entry.plan.uses.get(cycle.uses[(i+1)%cycle.uses.len()])
                        .ok_or_else(||invalid("G4b next source coedge not present"))?;
                    if !within(use_.end,next.start,eps) {
                        return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DamagedPatchLoop));
                    }
                    let id=StitchUseRef{patch:patch_index,use_index:at};
                    match use_.provenance {
                        MultiUvBoundaryProvenance::Source{edge,coedge,..}=>{
                            source.entry((owner.operand,edge))
                                .or_default().push(OriginalUse{id,coedge});
                        },
                        MultiUvBoundaryProvenance::Cut{segment,source:arc_ref}=>{
                            let Some(cut)=entry.plan.cut_pairs.iter().find(|cut|
                                cut.segment==segment && cut.source_arc==arc_ref &&
                                cut.selected_use==at) else {
                                return Ok(GlobalEdgeStitchResult::Held(
                                    GlobalStitchHold::NonManifoldIncidence));
                            };
                            let Some(pair)=entry.pairs.iter().find(|pair|
                                pair.id==arc_ref.pair) else {
                                return Err(invalid("G4b missing certified source face pair"));
                            };
                            let Some(arc)=pair.graph.arcs.iter().find(|arc|
                                arc.id==arc_ref.local) else {
                                return Err(invalid("G4b source pair lost graph carrier"));
                            };
                            if cut.counterpart_operand==owner.operand ||
                                (owner.operand==GraphOperand::First &&
                                 (owner.face!=pair.first_face ||
                                  cut.counterpart_face!=pair.second_face)) ||
                                (owner.operand==GraphOperand::Second &&
                                 (owner.face!=pair.second_face ||
                                  cut.counterpart_face!=pair.first_face))
                            {
                                return Ok(GlobalEdgeStitchResult::Held(
                                    GlobalStitchHold::CutCarrierMismatch));
                            }
                            let witnesses=endpoint_witnesses(&pair.graph,*arc)
                                .unwrap_or([false,false]);
                            cuts.entry((pair.first_face,pair.second_face,arc.id))
                                .or_default().push(CutUse {
                                    id,pair:(pair.first_face,pair.second_face),
                                    arc:*arc,boundary_evidence:witnesses,
                                });
                        }
                    }
                    selected+=1;
                    if selected>MAX_STITCH_USES {
                        return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::BudgetExceeded));
                    }
                }
                verified_loops+=1;
            }
            if used.len()!=entry.plan.uses.iter().filter(|u|u.kept).count() {
                return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::DamagedPatchLoop));
            }
        }
        let operation=operation.ok_or_else(||invalid("G4b no operation"))?;
        let mut edges=Vec::new();
        let mut open=Vec::new();
        let use_at=|id:StitchUseRef|->&crate::FacePatchDraftUse{
            &input[id.patch].plan.uses[id.use_index]
        };
        let owner_at=|id:StitchUseRef|->StitchFaceOwner{
            let plan=input[id.patch].plan;
            StitchFaceOwner{operand:plan.source_operand,face:plan.source_face}
        };
        // Original physical source edges have stable operand/EdgeId
        // namespaces. A valid matching split requires the actual original
        // coedge's own topological TWIN and reversed XYZ endpoints.
        for ((operand,edge),mut uses) in source {
            uses.sort_by_key(|u|u.id);
            let model=if operand==GraphOperand::First {&a}else{&b};
            let Some(native_edge)=model.edges.get(edge.0 as usize)
                .filter(|e|e.id==edge) else {
                return Err(invalid("G4b source EdgeId not present"));
            };
            let extent=model.bbox.max.sub(model.bbox.min).length_squared().sqrt();
            let eps=tolerance.length_at(extent.max(1.0))*16.0;
            let mut matched=BTreeSet::new();
            for i in 0..uses.len() {
                if matched.contains(&i) {continue;}
                let current=uses[i];
                let co=model.coedges.get(current.coedge.0 as usize)
                    .filter(|c|c.id==current.coedge && c.edge==edge)
                    .ok_or_else(||invalid("G4b source coedge missing"))?;
                if !native_edge.coedges.contains(&co.id) {
                    return Ok(GlobalEdgeStitchResult::Held(
                        GlobalStitchHold::ContradictorySourceCoedge));
                }
                let mut candidates=Vec::new();
                for j in i+1..uses.len() {
                    let other=uses[j];
                    if matched.contains(&j) {continue;}
                    let reverse_co=model.coedges.get(other.coedge.0 as usize)
                        .filter(|c|c.id==other.coedge && c.edge==edge)
                        .ok_or_else(||invalid("G4b opposite source coedge missing"))?;
                    if reverse_co.id==co.twin && reverse_co.twin==co.id &&
                        owner_at(current.id)!=owner_at(other.id) &&
                        pair_geometry(use_at(current.id),use_at(other.id),eps)
                    {
                        candidates.push(j);
                    }
                }
                if candidates.len()>1 {
                    return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::NonManifoldIncidence));
                }
                if let Some(&j)=candidates.first() {
                    matched.insert(i);matched.insert(j);
                    let other=uses[j];
                    edges.push(StitchedEdgeDraft {
                        draft_id:edges.len() as u32,
                        source:StitchedEdgeOrigin::Source{operand,edge},
                        radial_uses:[current.id,other.id],
                        source_faces:[owner_at(current.id),owner_at(other.id)],
                        start:use_at(current.id).start,end:use_at(current.id).end,
                        face_pcurves:[use_at(current.id).uv_span,use_at(other.id).uv_span],
                    });
                } else {
                    // Two uses on the same source EdgeId in the same direction
                    // are not radial twins; fail, do not treat the second
                    // as a harmless unmatched source face.
                    if uses.iter().enumerate().any(|(j,other)| j!=i &&
                        other.coedge==co.twin &&
                        owner_at(current.id)!=owner_at(other.id) &&
                        within(use_at(current.id).start,use_at(other.id).start,eps) &&
                        within(use_at(current.id).end,use_at(other.id).end,eps)
                    ) {
                        return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::SameDirectionUses));
                    }
                    open.push(OpenStitchUse {
                        use_ref:current.id,reason:OpenStitchReason::MissingSourceCoedgeTwin,
                    });
                }
            }
        }
        // Two *different operand* patches must use the same original graph
        // face pair and analytic arc, with opposite 3D orientation.
        for ((first,second,arc_id),mut uses) in cuts {
            uses.sort_by_key(|u|u.id);
            if uses.len()>2 {
                return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::NonManifoldIncidence));
            }
            let arc=uses[0].arc;
            let mut eps=arc.tolerance.evaluation_budget;
            eps=eps.max(tolerance.absolute_length*16.0);
            for u in &uses {
                if u.pair!=(first,second) || u.arc!=arc {
                    return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::CutCarrierMismatch));
                }
                if !arc_owns_entire_piece(arc,use_at(u.id),eps)? {
                    return Ok(GlobalEdgeStitchResult::Held(
                        GlobalStitchHold::UncertifiedIntermediateCutSplit));
                }
                if !u.boundary_evidence.iter().all(|x|*x) {
                    return Ok(GlobalEdgeStitchResult::Held(
                        GlobalStitchHold::UncertifiedCutEndpoint));
                }
            }
            if uses.len()==1 {
                open.push(OpenStitchUse {
                    use_ref:uses[0].id,reason:OpenStitchReason::MissingOppositeCutPatch,
                });
                continue;
            }
            let x=uses[0];let y=uses[1];
            let owx=owner_at(x.id);let owy=owner_at(y.id);
            if owx.operand==owy.operand ||
                owx.face==owy.face && owx.operand==owy.operand ||
                !(owx.operand==GraphOperand::First && owy.operand==GraphOperand::Second ||
                  owx.operand==GraphOperand::Second && owy.operand==GraphOperand::First)
            {
                return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::NonManifoldIncidence));
            }
            if !pair_geometry(use_at(x.id),use_at(y.id),eps) {
                return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::SameDirectionUses));
            }
            edges.push(StitchedEdgeDraft {
                draft_id:edges.len() as u32,
                source:StitchedEdgeOrigin::Cut {
                    first_face:first,second_face:second,arc:arc_id,
                },
                radial_uses:[x.id,y.id],
                source_faces:[owx,owy],
                start:use_at(x.id).start,end:use_at(x.id).end,
                face_pcurves:[use_at(x.id).uv_span,use_at(y.id).uv_span],
            });
        }
        let mut missing=Vec::new();
        for face in &a.faces {
            let owner=StitchFaceOwner{operand:GraphOperand::First,face:face.id};
            if !owners.contains(&owner){missing.push(owner);}
        }
        for face in &b.faces {
            let owner=StitchFaceOwner{operand:GraphOperand::Second,face:face.id};
            if !owners.contains(&owner){missing.push(owner);}
        }
        // Every matched physical edge uses two distinct source faces.
        // The remaining single uses are reported, never omitted.
        if 2*edges.len()+open.len()!=selected {
            return Ok(GlobalEdgeStitchResult::Held(GlobalStitchHold::NonManifoldIncidence));
        }
        let radially_paired=open.is_empty();
        Ok(GlobalEdgeStitchResult::Preview(Self {
            operation,edges,open_uses:open,missing_source_faces:missing,
            source_patch_count:input.len(),selected_coedge_uses:selected,
            verified_loops,radially_paired,
            closed_shell_proven:false,can_apply:false,
        }))
    }

    pub fn validate_against(
        &self,left:&BrepBody,right:&BrepBody,input:&[StitchPatchInput<'_>],
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        let rebuilt=Self::build(left,right,input,tolerance)?;
        if rebuilt!=GlobalEdgeStitchResult::Preview(self.clone()) {
            return Err(invalid("G4b radial incidence plan differs from current source"));
        }
        Ok(())
    }
}
