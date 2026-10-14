//! Source-bound, all-or-nothing multi-face intersection and local region
//! preview, now with optional complete planar UV cells for supported holes
//! and physical edge/dual-coedge reconciliation witnesses.
//!
//! This is NOT a global two-body face-cell arrangement or an executable
//! generalized B-rep Boolean; can_apply remains false for every outcome.
use crate::{
    BooleanOperation, BooleanRegionClassifier, BrepBody, BrepModel,
    CoedgeReconciliationPlan, Curve2, EdgeReconciliationResult,
    FaceId, GeometryError, GeometryTolerance, GraphOperand, GraphRegionPlan,
    GraphRegionState, GraphState, IntersectionGraph,
    IntersectionGraphBuilder, UvFaceArrangement, UvFaceRegionPlan,
};
use crate::RegionHold;

const MAX_FACES_PER_OPERAND:usize=64;
const MAX_FACE_PAIRS:usize=512;
fn invalid(msg:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(msg.into())
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiFaceRegionHold {
    BudgetExceeded,
    UnsupportedInnerLoop,
    UnsupportedSurfacePair,
    AmbiguousIntersection,
    IsolatedContact,
    AmbiguousSideClassification(RegionHold),
    UnresolvedCoedgeReconciliation,
    UnsupportedFaceCellArrangement,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiFaceRegionState{
    Empty,
    ClassifiedPreview,
    Held(MultiFaceRegionHold),
}
#[derive(Debug,Clone,PartialEq)]
pub struct FacePairRegionWitness{
    pub left_face:FaceId,
    pub right_face:FaceId,
    /// Bounded analytic 3D intersections and endpoint provenance.
    pub graph:IntersectionGraph,
    /// Local adjacent-side samples, not yet whole-body region labels.
    pub regions:GraphRegionPlan,
    /// Physical shared curve, paired opposite 3D coedge directions, and
    /// independent owner UV pcurves. NO new persistent B-rep entity IDs.
    pub reconciled:CoedgeReconciliationPlan,
    /// Whole planar UV cell loops when one source side permits a single
    /// strictly transverse line across a convex outer loop/circular holes.
    pub first_uv_cells:Option<UvFaceRegionPlan>,
    pub second_uv_cells:Option<UvFaceRegionPlan>,
}
#[derive(Debug,Clone,PartialEq)]
pub struct MultiFaceRegionPlan{
    pub state:MultiFaceRegionState,
    pub operation:BooleanOperation,
    /// Explicitly selected face pairs; recomputation checks the same set.
    /// In plan() this is the Cartesian product of the original operands.
    pub selected_pairs:Vec<(FaceId,FaceId)>,
    pub pairs:Vec<FacePairRegionWitness>,
    pub total_analytic_arcs:usize,
    pub total_side_samples:usize,
    pub total_reconciled_edges:usize,
    pub total_uv_cells:usize,
    pub can_apply:bool,
}
fn supported_inner_loops(model:&BrepModel,face:FaceId)->bool{
    let Some(f)=model.faces.get(face.0 as usize) else{return false;};
    for &id in f.loops.iter().skip(1){
        let Some(lp)=model.loops.get(id.0 as usize) else{return false;};
        let co=lp.first_coedge;
        let Some(c)=model.coedges.get(co.0 as usize) else{return false;};
        if c.face!=face || c.loop_id!=id || c.next!=co ||
            !matches!(model.coedge_pcurve(co),Some(Curve2::Circle(_))){
            return false;
        }
    }
    true
}

impl MultiFaceRegionPlan{
    fn held(operation:BooleanOperation,selected_pairs:Vec<(FaceId,FaceId)>,
        why:MultiFaceRegionHold)->Self
    {
        Self{
            state:MultiFaceRegionState::Held(why),operation,selected_pairs,
            pairs:vec![],total_analytic_arcs:0,total_side_samples:0,
            total_reconciled_edges:0,total_uv_cells:0,can_apply:false,
        }
    }
    pub fn validate_against(&self,left:&BrepBody,right:&BrepBody,
        tolerance:GeometryTolerance)->Result<(),GeometryError>
    {
        let reconstructed=Self::plan_selected(left,right,self.operation,
            tolerance,&self.selected_pairs)?;
        if self!=&reconstructed {
            return Err(invalid("multi-face preview differs from the selected source faces"));
        }
        Ok(())
    }
    /// All original source face pairs, with an all-or-nothing Held result if
    /// *any* pair is uncertain or unsupported. Never leaks earlier successes.
    pub fn plan(left:&BrepBody,right:&BrepBody,operation:BooleanOperation,
        tolerance:GeometryTolerance)->Result<Self,GeometryError>
    {
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let a=left.shared()?;let b=right.shared()?;
        let mut pairs=Vec::new();
        for fa in &a.faces{
            for fb in &b.faces{
                if pairs.len()>MAX_FACE_PAIRS{
                    return Ok(Self::held(operation,vec![],
                        MultiFaceRegionHold::BudgetExceeded));
                }
                pairs.push((fa.id,fb.id));
            }
        }
        Self::plan_selected(left,right,operation,tolerance,&pairs)
    }
    /// Useful for a bounded, source-chosen set of truly intersecting face
    /// pairs: e.g. the top cap of a through-bored box and a crossing plane.
    /// This is still only a local preview, not proof of a complete Boolean.
    pub fn plan_selected(left:&BrepBody,right:&BrepBody,
        operation:BooleanOperation,tolerance:GeometryTolerance,
        selections:&[(FaceId,FaceId)],
    )->Result<Self,GeometryError>{
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let lm=left.shared()?;let rm=right.shared()?;
        let mut selected_pairs=selections.to_vec();
        selected_pairs.sort_by_key(|(a,b)|(a.0,b.0));
        selected_pairs.dedup();
        if lm.faces.len()>MAX_FACES_PER_OPERAND ||
            rm.faces.len()>MAX_FACES_PER_OPERAND ||
            selected_pairs.len()>MAX_FACE_PAIRS {
            return Ok(Self::held(operation,selected_pairs,
                MultiFaceRegionHold::BudgetExceeded));
        }
        for &(fa,fb) in &selected_pairs{
            if lm.faces.get(fa.0 as usize).map(|f|f.id)!=Some(fa) ||
               rm.faces.get(fb.0 as usize).map(|f|f.id)!=Some(fb){
                return Err(invalid("selected face IDs are not present in source bodies"));
            }
            if !supported_inner_loops(&lm,fa)||!supported_inner_loops(&rm,fb){
                return Ok(Self::held(operation,selected_pairs,
                    MultiFaceRegionHold::UnsupportedInnerLoop));
            }
        }
        let engine=IntersectionGraphBuilder::new(tolerance)?;
        let mut pairs=Vec::new();
        let mut arcs=0usize;
        let mut samples=0usize;
        let mut coedge_count=0usize;
        let mut uv_count=0usize;
        for &(fa,fb) in &selected_pairs {
            let graph=match engine.intersect_faces(&lm,fa,&rm,fb){
                Ok(g)=>g,
                Err(GeometryError::InvalidEdit(_)) |
                Err(GeometryError::InvalidTopology(_))=>{
                    return Ok(Self::held(operation,selected_pairs,
                        MultiFaceRegionHold::UnsupportedSurfacePair));
                }
                Err(error)=>return Err(error),
            };
            match graph.state {
                GraphState::Empty=>continue,
                GraphState::Ambiguous(_)=>{
                    return Ok(Self::held(operation,selected_pairs,
                        MultiFaceRegionHold::AmbiguousIntersection));
                }
                GraphState::Ready=>{}
            }
            if graph.arcs.is_empty(){
                return Ok(Self::held(operation,selected_pairs,
                    MultiFaceRegionHold::IsolatedContact));
            }
            let regions=BooleanRegionClassifier::classify_graph(
                &graph,left,fa,right,fb,operation)?;
            match regions.state {
                GraphRegionState::Held=>return Ok(Self::held(
                    operation,selected_pairs,
                    MultiFaceRegionHold::AmbiguousSideClassification(
                        regions.hold.unwrap_or(RegionHold::UnresolvedProbe)))),
                GraphRegionState::NoCut=>return Ok(Self::held(
                    operation,selected_pairs,MultiFaceRegionHold::AmbiguousIntersection)),
                GraphRegionState::ClassifiedPreview=>{}
            }
            regions.validate_against(&graph,left,fa,right,fb)?;
            let reconciled=match CoedgeReconciliationPlan::build(
                &lm,&rm,&graph,tolerance)?{
                EdgeReconciliationResult::Ready(data)=>data,
                EdgeReconciliationResult::Held(_)=>return Ok(Self::held(
                    operation,selected_pairs,
                    MultiFaceRegionHold::UnresolvedCoedgeReconciliation)),
            };
            let get_uv=|operand:GraphOperand|
                ->Result<Option<UvFaceRegionPlan>,GeometryError>{
                match UvFaceRegionPlan::build(&lm,&rm,&graph,operand,tolerance)?{
                    UvFaceArrangement::Ready(data)=>Ok(Some(data)),
                    UvFaceArrangement::Held(_)=>Ok(None),
                }
            };
            let require_left=lm.faces[fa.0 as usize].loops.len()>1;
            let require_right=rm.faces[fb.0 as usize].loops.len()>1;
            let first_uv_cells=get_uv(GraphOperand::First)?;
            let second_uv_cells=get_uv(GraphOperand::Second)?;
            if (require_left&&first_uv_cells.is_none())||
               (require_right&&second_uv_cells.is_none()){
                return Ok(Self::held(operation,selected_pairs,
                    MultiFaceRegionHold::UnsupportedFaceCellArrangement));
            }
            // Even a preview must prove that its 2D material cut has
            // precisely the same physical source curve as the two
            // face-owned coedge candidates, with opposite in-cell uses.
            for region in [&first_uv_cells,&second_uv_cells].into_iter().flatten(){
                region.certify_physical_twins(
                    &lm,&rm,&graph,&reconciled,tolerance)?;
            }
            arcs=arcs.checked_add(graph.arcs.len())
                .ok_or_else(||invalid("multi-face arc count overflow"))?;
            samples=samples.checked_add(regions.samples.len())
                .ok_or_else(||invalid("multi-face sample count overflow"))?;
            coedge_count=coedge_count.checked_add(reconciled.edges.len())
                .ok_or_else(||invalid("multi-face coedge budget overflow"))?;
            uv_count=uv_count.checked_add(first_uv_cells.as_ref()
                .map(|x|x.cells.len()).unwrap_or(0)+
                second_uv_cells.as_ref().map(|x|x.cells.len()).unwrap_or(0))
                .ok_or_else(||invalid("multi-face cell budget overflow"))?;
            pairs.push(FacePairRegionWitness{
                left_face:fa,right_face:fb,graph,regions,reconciled,
                first_uv_cells,second_uv_cells,
            });
        }
        Ok(Self{
            state:if pairs.is_empty(){MultiFaceRegionState::Empty}
                else{MultiFaceRegionState::ClassifiedPreview},
            operation,selected_pairs,pairs,
            total_analytic_arcs:arcs,total_side_samples:samples,
            total_reconciled_edges:coedge_count,total_uv_cells:uv_count,
            can_apply:false,
        })
    }
}
