//! Independent mathematical evidence. Validation never depends on constructors,
//! source imprint/UV plans, G4 assembly, Booleans, documents or display meshes.
mod g8_source;
pub use g8_source::{g8_source_volume, G8SourceHold, G8SourceVolumeEvidence,
    G8SourceVolumeOutcome, GuardedSourceVolume, TiltedSectionAreaBracket};
