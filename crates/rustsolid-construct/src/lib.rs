//! Source-owned analytic construction of validated B-rep drafts.
//! A constructed shape is staged only and cannot be committed before
//! independent source-volume, actual-boundary flux, material and lineage proof.
mod g7_one_wall;
pub use g7_one_wall::{g7_one_wall_direct,G7BuildOutcome,G7BuildHold,G7StagedBrep};
