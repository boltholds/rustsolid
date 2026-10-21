//! Compare independently computed evidence. Matching volume is a necessary but
//! not sufficient condition for Boolean correctness. No CertifiedCandidate here.
use rustsolid_source::{BooleanIntent, SourceError, SourceSnapshot, SourceStamp};
use rustsolid_topology::BrepModel;
use super::{g8_boundary_flux, g8_source_volume, G8FluxHold, G8FluxOutcome,
    G8SourceHold, G8SourceVolumeOutcome};

#[derive(Debug, Clone, PartialEq)]
pub enum G8ComparisonHold {
    Source(G8SourceHold),
    Boundary(G8FluxHold),
    IncompatibleVolumes,
    NumericalBudgetExceeded,
}
#[derive(Debug, Clone, PartialEq)]
pub struct G8GuardedVolumeAgreement {
    /// Both reports are guarded f64 only; never claim directed interval rounding.
    pub source_interval:[f64;2],
    pub actual_boundary_volume:f64,
    pub numerical_slack:f64,
    pub model_incarnation:u64,
    pub source_stamp:SourceStamp,
    pub assurance_level: &'static str,
}
#[derive(Debug, Clone, PartialEq)]
pub enum G8ComparisonOutcome {
    VolumeCompatible(G8GuardedVolumeAgreement),
    Held(G8ComparisonHold),
}

/// Recompute BOTH reports from immutable inputs and the actual canonical B-rep.
/// The result deliberately has NO 'can_apply' or Boolean approval status.
pub fn compare_g8_volume(
    snapshot:&SourceSnapshot,
    intent:&BooleanIntent,
    candidate:&BrepModel,
)->Result<G8ComparisonOutcome,SourceError>{
    let source=match g8_source_volume(snapshot,intent)? {
        G8SourceVolumeOutcome::Ready(e)=>e,
        G8SourceVolumeOutcome::Held(hold)=>
            return Ok(G8ComparisonOutcome::Held(G8ComparisonHold::Source(hold))),
    };
    let actual=match g8_boundary_flux(candidate) {
        G8FluxOutcome::Ready(e)=>e,
        G8FluxOutcome::Held(hold)=>
            return Ok(G8ComparisonOutcome::Held(G8ComparisonHold::Boundary(hold))),
    };
    let scale=actual.length_scale;
    let tol=intent.tolerance();
    let budget=4.0*tol.length_at(scale)*scale*scale;
    let slack=budget+actual.numerical_guard+source.numerical_guard;
    if !slack.is_finite()||slack<=0.0 {
        return Ok(G8ComparisonOutcome::Held(G8ComparisonHold::NumericalBudgetExceeded));
    }
    if actual.retained_volume < source.retained.lower-slack ||
       actual.retained_volume > source.retained.upper+slack {
        return Ok(G8ComparisonOutcome::Held(G8ComparisonHold::IncompatibleVolumes));
    }
    Ok(G8ComparisonOutcome::VolumeCompatible(G8GuardedVolumeAgreement {
        source_interval:[source.retained.lower,source.retained.upper],
        actual_boundary_volume:actual.retained_volume,
        numerical_slack:slack,
        model_incarnation:actual.model_incarnation,
        source_stamp:source.source_stamp().clone(),
        assurance_level:"GuardedF64",
    }))
}
