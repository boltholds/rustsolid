//! Independent mathematical evidence. Validation never depends on constructors,
//! source imprint/UV plans, G4 assembly, Booleans, documents or display meshes.
mod g8_source;
pub use g8_source::{g8_source_volume, G8SourceHold, G8SourceVolumeEvidence,
    G8SourceVolumeOutcome, GuardedSourceVolume, TiltedSectionAreaBracket};

mod g8_flux;
pub use g8_flux::{g8_boundary_flux, G8BoundaryFlux, G8FluxHold, G8FluxOutcome};

mod g8_compare;
pub use g8_compare::{compare_g8_volume, G8ComparisonHold, G8ComparisonOutcome,
    G8GuardedVolumeAgreement};
