//! Independent source-only and result-boundary evidence; not Boolean acceptance.
pub mod source_volume;
pub use source_volume::{evaluate_tilted_one_wall, SourceVolumeEvidence, TiltedVolumeHold, VolumeInterval};

pub mod boundary_flux;
pub use boundary_flux::{evaluate_tilted_boundary_flux,BoundaryFluxEvidence,BoundaryFluxHold};
