//! Immutable multi-face intersection/region-preview batch. It is an
//! all-or-nothing *set of pair-local witnesses*, not a globally partitioned
//! face-cell arrangement and not authorization to execute a B-rep Boolean.
use crate::{
    BooleanOperation, BooleanRegionClassifier, BrepBody, FaceId, GeometryError,
    GeometryTolerance, GraphRegionPlan, GraphRegionState, GraphState,
    IntersectionGraph, IntersectionGraphBuilder,
};
use crate::RegionHold;

const MAX_FACES_PER_OPERAND:usize=64;
const MAX_FACE_PAIRS:usize=512;

fn invalid(msg:impl Into<String>)->GeometryError {
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
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiFaceRegionState {
    Empty,
    ClassifiedPreview,
    Held(MultiFaceRegionHold),
}
#[derive(Debug,Clone,PartialEq)]
pub struct FacePairRegionWitness {
    pub left_face:FaceId,
    pub right_face:FaceId,
    /// Direct source geometry and UV witnesses, never mesh guesses.
    pub graph:IntersectionGraph,
    /// Local three-sample-per-arc classification on both source faces.
    /// This is not a global material-region partition.
    pub regions:GraphRegionPlan,
}
#[derive(Debug,Clone,PartialEq)]
pub struct MultiFaceRegionPlan {
    pub state:MultiFaceRegionState,
    pub operation:BooleanOperation,
    pub pairs:Vec<FacePairRegionWitness>,
    pub total_analytic_arcs:usize,
    pub total_side_samples:usize,
    /// Always false: complete global face cells, topological reconciliation,
    /// and source mass still need separate mathematical certification.
    pub can_apply:bool,
}
impl MultiFaceRegionPlan {
    fn held(operation:BooleanOperation,why:MultiFaceRegionHold)->Self{
        Self{state:MultiFaceRegionState::Held(why),operation,
            pairs:vec![],total_analytic_arcs:0,total_side_samples:0,
            can_apply:false}
    }
    pub fn validate_against(&self,left:&BrepBody,right:&BrepBody,
        tolerance:GeometryTolerance)->Result<(),GeometryError>
    {
        let original=Self::plan(left,right,self.operation,tolerance)?;
        if self!=&original {return Err(invalid(
            "multi-face preview has changed, or source geometry/tolerances no longer match"));
        }
        Ok(())
    }
    pub fn plan(left:&BrepBody,right:&BrepBody,operation:BooleanOperation,
        tolerance:GeometryTolerance)->Result<Self,GeometryError>
    {
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let lm=left.shared()?;let rm=right.shared()?;
        let count=lm.faces.len().checked_mul(rm.faces.len())
            .ok_or_else(||invalid("face pair count overflow"))?;
        if lm.faces.len()>MAX_FACES_PER_OPERAND ||
            rm.faces.len()>MAX_FACES_PER_OPERAND ||
            count>MAX_FACE_PAIRS {
            return Ok(Self::held(operation,MultiFaceRegionHold::BudgetExceeded));
        }
        // One unsupported inner trim would invalidate any "complete" face
        // arrangement, even if earlier individual face pairs were solved.
        // Return no partial arcs/witnesses rather than hiding the hole.
        if lm.faces.iter().any(|f|f.loops.len()!=1) ||
            rm.faces.iter().any(|f|f.loops.len()!=1) {
            return Ok(Self::held(operation,MultiFaceRegionHold::UnsupportedInnerLoop));
        }

        let engine=IntersectionGraphBuilder::new(tolerance)?;
        let mut pairs=Vec::new();
        let mut arcs=0usize;
        let mut samples=0usize;
        for left_face in &lm.faces {
            for right_face in &rm.faces {
                let (fa,fb)=(left_face.id,right_face.id);
                let graph=match engine.intersect_faces(&lm,fa,&rm,fb){
                    Ok(g)=>g,
                    Err(GeometryError::InvalidEdit(_)) |
                    Err(GeometryError::InvalidTopology(_))=>{
                        return Ok(Self::held(operation,MultiFaceRegionHold::UnsupportedSurfacePair));
                    }
                    Err(e)=>return Err(e),
                };
                match graph.state {
                    GraphState::Empty=>continue,
                    GraphState::Ambiguous(_)=>return Ok(Self::held(
                        operation,MultiFaceRegionHold::AmbiguousIntersection)),
                    GraphState::Ready=>{}
                }
                if graph.arcs.is_empty(){
                    // Contacts do not create a Boolean cut by themselves.
                    return Ok(Self::held(operation,MultiFaceRegionHold::IsolatedContact));
                }
                let regions=BooleanRegionClassifier::classify_graph(
                    &graph,left,fa,right,fb,operation)?;
                match regions.state {
                    GraphRegionState::Held=>{
                        return Ok(Self::held(operation,
                            MultiFaceRegionHold::AmbiguousSideClassification(
                                regions.hold.unwrap_or(RegionHold::UnresolvedProbe))));
                    }
                    GraphRegionState::NoCut=>return Ok(Self::held(
                        operation,MultiFaceRegionHold::AmbiguousIntersection)),
                    GraphRegionState::ClassifiedPreview=>{}
                }
                regions.validate_against(&graph,left,fa,right,fb)?;
                arcs=arcs.checked_add(graph.arcs.len())
                    .ok_or_else(||invalid("multi-face arc count overflow"))?;
                samples=samples.checked_add(regions.samples.len())
                    .ok_or_else(||invalid("multi-face sample count overflow"))?;
                pairs.push(FacePairRegionWitness{
                    left_face:fa,right_face:fb,graph,regions,
                });
            }
        }
        Ok(Self{
            state:if pairs.is_empty(){MultiFaceRegionState::Empty}
                else{MultiFaceRegionState::ClassifiedPreview},
            operation,pairs,total_analytic_arcs:arcs,
            total_side_samples:samples,can_apply:false,
        })
    }
}
