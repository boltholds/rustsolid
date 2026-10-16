//! Read-only, operand-namespaced multi-face intersection collection.
//! Pair graphs remain independently certified. Cross-pair intersections are
//! NOT welded: global vertex equivalence requires a separate certified stage.
use crate::{
    BrepModel, FaceId, GeometryError, GeometryTolerance,
    GraphArcId, GraphNodeId, GraphState, IntersectionGraph,
    IntersectionGraphBuilder, Point3,
};

const MAX_GRAPH_PAIRS: usize = 512;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FacePairId(pub u32);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GlobalNodeRef {
    pub pair: FacePairId,
    pub local: GraphNodeId,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GlobalArcRef {
    pub pair: FacePairId,
    pub local: GraphArcId,
}
#[derive(Debug, Clone, PartialEq)]
pub struct FacePairGraph {
    pub id: FacePairId,
    pub first_face: FaceId,
    pub second_face: FaceId,
    pub graph: IntersectionGraph,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiFaceGraphState {
    Empty,
    Ready,
    Held,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MultiFaceIntersectionGraph {
    pub pairs: Vec<FacePairGraph>,
    pub state: MultiFaceGraphState,
    /// Original construction policy for deterministic source-bound checks.
    pub tolerance: GeometryTolerance,
}
impl MultiFaceIntersectionGraph {
    pub fn node(&self, reference: GlobalNodeRef) -> Option<Point3> {
        let pair=self.pairs.get(reference.pair.0 as usize)?;
        if pair.id!=reference.pair {return None;}
        let node=pair.graph.nodes.get(reference.local.0 as usize)?;
        (node.id==reference.local).then_some(node.point)
    }
    pub fn arc(&self, reference: GlobalArcRef) -> Option<&crate::GraphArc> {
        let pair=self.pairs.get(reference.pair.0 as usize)?;
        if pair.id!=reference.pair {return None;}
        pair.graph.arcs.get(reference.local.0 as usize)
            .filter(|arc|arc.id==reference.local)
    }
    pub fn validate_against(&self, first:&BrepModel, second:&BrepModel)
        -> Result<(),GeometryError>
    {
        first.validate()?;second.validate()?;
        self.tolerance.validate()?;
        let builder=IntersectionGraphBuilder::new(self.tolerance)?;
        let expected_pairs=first.faces.len().checked_mul(second.faces.len())
            .ok_or_else(||GeometryError::InvalidTopology("face pair count overflow".into()))?;
        if self.pairs.len()!=expected_pairs {
            return Err(GeometryError::InvalidTopology("incomplete multi-face graph".into()));
        }
        let mut held=false;
        let mut ready=false;
        for (index,pair) in self.pairs.iter().enumerate() {
            let i=index/second.faces.len();
            let j=index%second.faces.len();
            if pair.id.0 as usize!=index ||
                pair.first_face!=first.faces[i].id ||
                pair.second_face!=second.faces[j].id {
                return Err(GeometryError::InvalidTopology("multi-face pair provenance changed".into()));
            }
            pair.graph.validate_against(first,pair.first_face,second,pair.second_face)?;
            // The empty graph has no geometric witness of completeness.
            // Recompute, rather than trusting a user-supplied Empty state.
            let expected=builder.intersect_faces(first,pair.first_face,second,pair.second_face)?;
            if pair.graph!=expected {
                return Err(GeometryError::InvalidTopology(
                    "multi-face graph differs from recomputed source intersections".into()
                ));
            }
            match pair.graph.state {
                GraphState::Empty=>{},
                GraphState::Ready=>ready=true,
                GraphState::Ambiguous(_)=>held=true,
            }
        }
        let state=if held{MultiFaceGraphState::Held}
            else if ready{MultiFaceGraphState::Ready}
            else{MultiFaceGraphState::Empty};
        if self.state!=state {
            return Err(GeometryError::InvalidTopology("multi-face state mismatch".into()));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy)]
pub struct MultiFaceIntersectionGraphBuilder {
    pair_builder: IntersectionGraphBuilder,
}
impl MultiFaceIntersectionGraphBuilder {
    pub fn new(tolerance:GeometryTolerance)->Result<Self,GeometryError>{
        Ok(Self{pair_builder:IntersectionGraphBuilder::new(tolerance)?})
    }
    pub fn build(&self,first:&BrepModel,second:&BrepModel)
        ->Result<MultiFaceIntersectionGraph,GeometryError>
    {
        first.validate()?;second.validate()?;
        let count=first.faces.len().checked_mul(second.faces.len())
            .ok_or_else(||GeometryError::InvalidTopology("face pair count overflow".into()))?;
        if count>MAX_GRAPH_PAIRS {
            return Err(GeometryError::InvalidTopology(
                "bounded multi-face graph pair budget exceeded".into()));
        }
        let mut pairs=Vec::with_capacity(count);
        let mut held=false;
        let mut ready=false;
        for a in &first.faces {
            for b in &second.faces {
                let graph=self.pair_builder.intersect_faces(first,a.id,second,b.id)?;
                match graph.state {
                    GraphState::Empty=>{},
                    GraphState::Ready=>ready=true,
                    GraphState::Ambiguous(_)=>held=true,
                }
                pairs.push(FacePairGraph{
                    id:FacePairId(pairs.len() as u32),
                    first_face:a.id,second_face:b.id,graph,
                });
            }
        }
        let state=if held{MultiFaceGraphState::Held}
            else if ready{MultiFaceGraphState::Ready}
            else{MultiFaceGraphState::Empty};
        Ok(MultiFaceIntersectionGraph{pairs,state,tolerance:self.pair_builder.tolerance})
    }
}
