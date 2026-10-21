//! G8 numerical consistency witness. Both evaluators are independent;
//! this is NOT Boolean material/lineage certification or a formal interval proof.
use rustsolid_foundation::GeometryTolerance;
use rustsolid_source::{BooleanIntent,SourceSnapshot,SourceStamp,SourceRef};
use crate::{SourceVolumeEvidence,BoundaryFluxEvidence};
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AssuranceLevel { GuardedF64 }
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AgreementHold {
    StaleSource, DifferentSource, InconsistentSourceBounds,
    InconsistentActualFlux, InvalidBudget,
}
#[derive(Debug,Clone)]
pub struct VolumeAgreement {
    stamp:SourceStamp,
    intent_stamp:SourceStamp,
    source_refs:[SourceRef;2],
    retained_interval:[f64;2],
    measured_flux:f64,
    budget:f64,
    assurance:AssuranceLevel,
}
impl VolumeAgreement {
    pub fn source_stamp(&self)->&SourceStamp{&self.stamp}
    pub fn intent_stamp(&self)->&SourceStamp{&self.intent_stamp}
    pub fn source_refs(&self)->[SourceRef;2]{self.source_refs}
    pub fn retained_interval(&self)->[f64;2]{self.retained_interval}
    pub fn measured_flux(&self)->f64{self.measured_flux}
    pub fn budget(&self)->f64{self.budget}
    pub fn assurance(&self)->AssuranceLevel{self.assurance}
}
/// This gate verifies only independent volumes. It is deliberately incapable
/// of certifying source-face ownership or material region classification.
pub fn compare_g8_numeric(
    snapshot:&SourceSnapshot, intent:&BooleanIntent,
    source:&SourceVolumeEvidence, flux:&BoundaryFluxEvidence,
)->Result<VolumeAgreement,AgreementHold>{
    snapshot.check_intent(intent).map_err(|_|AgreementHold::StaleSource)?;
    if source.source_stamp()!=snapshot.stamp()||
        source.intent_stamp()!=intent.intent_stamp()||
        source.source_refs()!=[intent.operands().0,intent.operands().1]{
        return Err(AgreementHold::DifferentSource);
    }
    let interval=source.interval();
    let [lo,hi]=[interval.lower,interval.upper];
    if !lo.is_finite()||!hi.is_finite()||lo<=0.||hi<=lo||
        !source.original_box_volume().is_finite()||
        source.original_box_volume()<=hi{
        return Err(AgreementHold::InconsistentSourceBounds);
    }
    let value=flux.result_volume();
    if !value.is_finite()||value<=0.{
        return Err(AgreementHold::InconsistentActualFlux);
    }
    let tol:GeometryTolerance=intent.tolerance();
    let scale=source.original_box_volume().cbrt().max(1.0);
    let budget=4.*tol.length_at(scale)*scale*scale
        +source.numerical_guard()+flux.numerical_guard();
    if !budget.is_finite()||budget<0.||budget>source.original_box_volume()*1e-4{
        return Err(AgreementHold::InvalidBudget);
    }
    if value<lo-budget||value>hi+budget{
        return Err(AgreementHold::InconsistentActualFlux);
    }
    Ok(VolumeAgreement{
        stamp:snapshot.stamp().clone(),intent_stamp:intent.intent_stamp().clone(),
        source_refs:[intent.operands().0,intent.operands().1],
        retained_interval:[lo,hi],measured_flux:value,budget,
        assurance:AssuranceLevel::GuardedF64,
    })
}
