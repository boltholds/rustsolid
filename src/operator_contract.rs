//! Operator inventory used by the mandatory first-party mathematical quality gate.
//! A new EulerCommand variant must be added to the exhaustive operator_kind()
//! match; the independent JSON scenario inventory must then cover it.
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ModelingOperator {
    Block, ExtrudeXz, AnalyticCylinder,
    SplitEdge, SplitFace, MakeEdgeVertex, MakeEdgeFace,
    KillEdgeVertex, KillEdgeFace,
}

pub const ALL_MODELING_OPERATORS: [ModelingOperator; 9] = [
    ModelingOperator::Block, ModelingOperator::ExtrudeXz,
    ModelingOperator::AnalyticCylinder, ModelingOperator::SplitEdge,
    ModelingOperator::SplitFace, ModelingOperator::MakeEdgeVertex,
    ModelingOperator::MakeEdgeFace, ModelingOperator::KillEdgeVertex,
    ModelingOperator::KillEdgeFace,
];

impl ModelingOperator {
    pub const fn is_command(self) -> bool {
        !matches!(self, Self::Block | Self::ExtrudeXz | Self::AnalyticCylinder)
    }
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Block=>"block", Self::ExtrudeXz=>"extrude_xz",
            Self::AnalyticCylinder=>"analytic_cylinder",
            Self::SplitEdge=>"split_edge", Self::SplitFace=>"split_face",
            Self::MakeEdgeVertex=>"make_edge_vertex",
            Self::MakeEdgeFace=>"make_edge_face",
            Self::KillEdgeVertex=>"kill_edge_vertex",
            Self::KillEdgeFace=>"kill_edge_face",
        }
    }
}
