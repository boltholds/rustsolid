//! Immutable, versioned source geometry. No B-rep, Boolean executor, or document dependency.
//! A SourceStamp is collision-free canonical bytes, NOT a cryptographic hash.
//! CAD geometry is never inferred from display tessellation.

use rustsolid_foundation::{GeometryError, GeometryTolerance, Point3};
use rustsolid_geometry::{Frame3, ModelUnits, ParameterRange, LengthUnit, AngleUnit};
use std::collections::BTreeMap;
use std::fmt;
use std::num::NonZeroU128;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourceError {
    InvalidId,
    InvalidRevision,
    InvalidGeometry(&'static str),
    InvalidFrame,
    DuplicateId(SourceId),
    MissingId(SourceId),
    StaleVersion { id: SourceId, expected: SourceRevision, actual: SourceRevision },
    RevisionNotSuccessor,
    UnitConversionRequired,
    InvalidUnitTransform,
    SnapshotMismatch,
    Geometry(GeometryError),
}
impl fmt::Display for SourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for SourceError {}
impl From<GeometryError> for SourceError {
    fn from(e: GeometryError) -> Self { Self::Geometry(e) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceId(NonZeroU128);
impl SourceId {
    pub fn new(value: u128) -> Result<Self, SourceError> {
        Ok(Self(NonZeroU128::new(value).ok_or(SourceError::InvalidId)?))
    }
    pub fn get(self) -> u128 { self.0.get() }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceRevision(u64);
impl SourceRevision {
    pub fn new(value: u64) -> Result<Self, SourceError> {
        if value == 0 { return Err(SourceError::InvalidRevision); }
        Ok(Self(value))
    }
    pub fn get(self) -> u64 { self.0 }
    pub fn next(self) -> Result<Self, SourceError> {
        Self::new(self.0.checked_add(1).ok_or(SourceError::InvalidRevision)?)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceRef {
    pub id: SourceId,
    pub revision: SourceRevision,
}

/// Geometry parameters describe original material, never a tessellated approximation.
/// Extrusions and imports are intentionally absent from this first G7/G8 source slice.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SourceGeometry {
    AxisAlignedBlock { min: Point3, max: Point3 },
    FiniteCylinder { frame: Frame3, radius: f64, axial: ParameterRange },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SourceSpec {
    id: SourceId,
    revision: SourceRevision,
    units: ModelUnits,
    geometry: SourceGeometry,
    validation_tolerance: GeometryTolerance,
}

impl SourceSpec {
    pub fn new(
        id: SourceId, revision: SourceRevision, units: ModelUnits,
        geometry: SourceGeometry, tolerance: GeometryTolerance,
    ) -> Result<Self, SourceError> {
        tolerance.validate()?;
        validate_geometry(&geometry, tolerance)?;
        Ok(Self { id, revision, units, geometry, validation_tolerance: tolerance })
    }
    pub fn id(&self) -> SourceId { self.id }
    pub fn revision(&self) -> SourceRevision { self.revision }
    pub fn units(&self) -> ModelUnits { self.units }
    pub fn geometry(&self) -> &SourceGeometry { &self.geometry }
    pub fn tolerance(&self) -> GeometryTolerance { self.validation_tolerance }
    pub fn reference(&self) -> SourceRef { SourceRef { id:self.id, revision:self.revision } }
    /// Returns a new source generation; previously issued specs and snapshots remain unchanged.
    pub fn successor(
        &self, new_geometry: SourceGeometry, tolerance: GeometryTolerance,
    ) -> Result<Self, SourceError> {
        Self::new(self.id, self.revision.next()?, self.units, new_geometry, tolerance)
    }
}

/// Validates geometrical *source parameters*, not resulting B-rep correctness.
fn validate_geometry(value: &SourceGeometry, tolerance: GeometryTolerance) -> Result<(), SourceError> {
    match value {
        SourceGeometry::AxisAlignedBlock { min, max } => {
            if !min.is_finite() || !max.is_finite() {
                return Err(SourceError::InvalidGeometry("box coordinates must be finite"));
            }
            for (lo, hi) in [(min.x,max.x),(min.y,max.y),(min.z,max.z)] {
                let extent = hi - lo;
                if !extent.is_finite() || extent <= tolerance.length_at(extent.abs()) {
                    return Err(SourceError::InvalidGeometry("box dimensions below local tolerance"));
                }
            }
        }
        SourceGeometry::FiniteCylinder { frame, radius, axial } => {
            validate_frame(frame, tolerance)?;
            let (ParameterRange::Bounded {start, end}) = *axial else {
                return Err(SourceError::InvalidGeometry("cylinder axial interval must be bounded"));
            };
            let span = end - start;
            if !start.is_finite() || !end.is_finite() || !span.is_finite()
                || span <= tolerance.length_at(span.abs())
                || !radius.is_finite() || *radius <= tolerance.absolute_length
            {
                return Err(SourceError::InvalidGeometry("invalid cylinder dimensions"));
            }
        }
    }
    Ok(())
}

/// Frame3 fields are publicly constructible in the geometry crate: verify the actual basis.
fn validate_frame(f: &Frame3, tol: GeometryTolerance) -> Result<(), SourceError> {
    if !f.origin.is_finite() || !f.x.is_finite() || !f.y.is_finite() || !f.normal.is_finite() {
        return Err(SourceError::InvalidFrame);
    }
    let norm = |v: Point3| v.x.hypot(v.y).hypot(v.z);
    let slack = tol.angular.min(1.0e-3) + 64.0 * f64::EPSILON;
    if [f.x, f.y, f.normal].into_iter().any(|v| (norm(v)-1.0).abs() > slack)
        || f.x.dot(f.y).abs() > slack
        || f.x.dot(f.normal).abs() > slack
        || f.y.dot(f.normal).abs() > slack
        || norm(f.x.cross(f.y).sub(f.normal)) > slack
    {
        return Err(SourceError::InvalidFrame);
    }
    Ok(())
}

/// Canonical binary encoding with explicit ordering, endian and float-zero rules.
/// This is a full byte string (not a digest) to avoid hash collision assumptions.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SourceStamp(Arc<[u8]>);
impl SourceStamp {
    pub fn as_bytes(&self) -> &[u8] { &self.0 }
}
fn number(bytes: &mut Vec<u8>, x: f64) {
    // All scalar values were validated before this helper is invoked.
    debug_assert!(x.is_finite());
    bytes.extend_from_slice(&(if x == 0.0 { 0.0 } else { x }).to_le_bytes());
}
fn point(bytes: &mut Vec<u8>, p: Point3) {
    number(bytes,p.x); number(bytes,p.y); number(bytes,p.z);
}
fn units(bytes: &mut Vec<u8>, u: ModelUnits) {
    bytes.push(match u.length { LengthUnit::Millimeter=>0, LengthUnit::Meter=>1, LengthUnit::Inch=>2 });
    bytes.push(match u.angle { AngleUnit::Radian=>0 });
}
fn spec_bytes(bytes: &mut Vec<u8>, spec: &SourceSpec) {
    bytes.extend_from_slice(&spec.id.get().to_be_bytes());
    bytes.extend_from_slice(&spec.revision.get().to_le_bytes());
    units(bytes,spec.units);
    number(bytes,spec.validation_tolerance.absolute_length);
    number(bytes,spec.validation_tolerance.relative_length);
    number(bytes,spec.validation_tolerance.angular);
    match &spec.geometry {
        SourceGeometry::AxisAlignedBlock {min,max} => {
            bytes.push(1);
            point(bytes,*min); point(bytes,*max);
        }
        SourceGeometry::FiniteCylinder {frame,radius,axial} => {
            bytes.push(2);
            point(bytes,frame.origin); point(bytes,frame.x);
            point(bytes,frame.y); point(bytes,frame.normal);
            number(bytes,*radius);
            let ParameterRange::Bounded {start,end} = *axial else { unreachable!("validated") };
            number(bytes,start); number(bytes,end);
        }
    }
}

#[derive(Clone, Debug)]
pub struct SourceSnapshot {
    sources: Arc<BTreeMap<SourceId, Arc<SourceSpec>>>,
    stamp: SourceStamp,
}
impl SourceSnapshot {
    pub fn empty() -> Self {
        Self::new(std::iter::empty::<SourceSpec>()).expect("empty snapshot is valid")
    }
    pub fn new(items: impl IntoIterator<Item=SourceSpec>) -> Result<Self, SourceError> {
        let mut sources = BTreeMap::new();
        for item in items {
            let id = item.id();
            if sources.insert(id, Arc::new(item)).is_some() {
                return Err(SourceError::DuplicateId(id));
            }
        }
        Self::from_map(sources)
    }
    fn from_map(sources: BTreeMap<SourceId, Arc<SourceSpec>>) -> Result<Self, SourceError> {
        let len = u32::try_from(sources.len())
            .map_err(|_| SourceError::InvalidGeometry("too many source specifications"))?;
        let mut bytes = b"rustsolid.source.snapshot.v1\0".to_vec();
        bytes.extend_from_slice(&len.to_le_bytes());
        for spec in sources.values() { spec_bytes(&mut bytes, spec); }
        Ok(Self { sources: Arc::new(sources), stamp: SourceStamp(bytes.into()) })
    }
    pub fn stamp(&self) -> &SourceStamp { &self.stamp }
    pub fn len(&self) -> usize { self.sources.len() }
    pub fn is_empty(&self) -> bool { self.sources.is_empty() }
    pub fn get(&self, id: SourceId) -> Option<&SourceSpec> {
        self.sources.get(&id).map(AsRef::as_ref)
    }
    pub fn resolve(&self, reference: SourceRef) -> Result<&SourceSpec, SourceError> {
        let current = self.get(reference.id).ok_or(SourceError::MissingId(reference.id))?;
        if current.revision != reference.revision {
            return Err(SourceError::StaleVersion {
                id:reference.id, expected:reference.revision, actual:current.revision,
            });
        }
        Ok(current)
    }
    pub fn with_added(&self, next: SourceSpec) -> Result<Self, SourceError> {
        if self.sources.contains_key(&next.id) { return Err(SourceError::DuplicateId(next.id)); }
        let mut sources = (*self.sources).clone();
        sources.insert(next.id,Arc::new(next));
        Self::from_map(sources)
    }
    /// Strong atomic functional edit: replaces exactly one source with its immediate successor.
    pub fn with_replacement(&self, next: SourceSpec) -> Result<Self, SourceError> {
        let current = self.get(next.id).ok_or(SourceError::MissingId(next.id))?;
        if current.revision.next()? != next.revision {
            return Err(SourceError::RevisionNotSuccessor);
        }
        let mut sources = (*self.sources).clone();
        sources.insert(next.id,Arc::new(next));
        Self::from_map(sources)
    }
    /// Constructs a version-bound operation without changing this snapshot.
    pub fn intent(
        &self, operation: BooleanKind, left: SourceRef, right: SourceRef,
        tolerance: GeometryTolerance, policy: UnitPolicy,
    ) -> Result<BooleanIntent, SourceError> {
        tolerance.validate()?;
        let a = self.resolve(left)?;
        let b = self.resolve(right)?;
        validate_policy(a.units,b.units,policy)?;
        let mut bytes = b"rustsolid.source.intent.v1\0".to_vec();
        bytes.extend_from_slice(self.stamp.as_bytes());
        bytes.push(match operation { BooleanKind::Union=>1, BooleanKind::Intersection=>2,
                                      BooleanKind::Difference=>3 });
        for r in [left,right] {
            bytes.extend_from_slice(&r.id.get().to_be_bytes());
            bytes.extend_from_slice(&r.revision.get().to_le_bytes());
        }
        number(&mut bytes,tolerance.absolute_length);
        number(&mut bytes,tolerance.relative_length);
        number(&mut bytes,tolerance.angular);
        match policy {
            UnitPolicy::RequireIdentical => bytes.push(0),
            UnitPolicy::RecordedRightToLeft{factor} => {
                bytes.push(1); number(&mut bytes,factor);
            }
        }
        Ok(BooleanIntent {
            operation,left,right,tolerance,unit_policy:policy,
            source_stamp:self.stamp.clone(), intent_stamp:SourceStamp(bytes.into()),
        })
    }
    pub fn check_intent(&self, intent: &BooleanIntent) -> Result<(), SourceError> {
        if self.stamp != intent.source_stamp { return Err(SourceError::SnapshotMismatch); }
        let expected = self.intent(intent.operation,intent.left,intent.right,
            intent.tolerance,intent.unit_policy)?;
        if expected.intent_stamp != intent.intent_stamp {
            return Err(SourceError::SnapshotMismatch);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanKind { Union, Intersection, Difference }
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitPolicy {
    RequireIdentical,
    /// Explicit right-operand numerical scale into the left-operand model units.
    /// The later constructor must actually apply this conversion.
    RecordedRightToLeft {factor:f64},
}
fn validate_policy(left:ModelUnits,right:ModelUnits,policy:UnitPolicy)->Result<(),SourceError>{
    match policy {
        UnitPolicy::RequireIdentical if left != right => Err(SourceError::UnitConversionRequired),
        UnitPolicy::RequireIdentical => Ok(()),
        UnitPolicy::RecordedRightToLeft {factor} => {
            let expected = right.length.meters_per_unit()/left.length.meters_per_unit();
            if !factor.is_finite() || factor<=0.0 || !expected.is_finite()
                || (factor-expected).abs()>32.0*f64::EPSILON*expected.max(1.0) {
                return Err(SourceError::InvalidUnitTransform);
            }
            Ok(())
        }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct BooleanIntent {
    operation:BooleanKind,
    left:SourceRef,
    right:SourceRef,
    tolerance:GeometryTolerance,
    unit_policy:UnitPolicy,
    source_stamp:SourceStamp,
    intent_stamp:SourceStamp,
}
impl BooleanIntent {
    pub fn operation(&self)->BooleanKind{self.operation}
    pub fn operands(&self)->(SourceRef,SourceRef){(self.left,self.right)}
    pub fn tolerance(&self)->GeometryTolerance{self.tolerance}
    pub fn unit_policy(&self)->UnitPolicy{self.unit_policy}
    pub fn source_stamp(&self)->&SourceStamp{&self.source_stamp}
    pub fn intent_stamp(&self)->&SourceStamp{&self.intent_stamp}
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
    fn id(i:u128)->SourceId{SourceId::new(i).unwrap()}
    fn rev(i:u64)->SourceRevision{SourceRevision::new(i).unwrap()}
    fn box_geom(x:f64)->SourceGeometry {
        SourceGeometry::AxisAlignedBlock {min:p(x,0.,0.),max:p(x+2.,3.,4.)}
    }
    fn sample(i:u128,r:u64,x:f64)->SourceSpec {
        SourceSpec::new(id(i),rev(r),ModelUnits::default(),box_geom(x),
            GeometryTolerance::default()).unwrap()
    }
    fn cylinder()->SourceSpec {
        let frame=Frame3::from_axes(p(2.,0.,1.),p(1.,0.,0.),p(0.,1.,0.),
            GeometryTolerance::default()).unwrap();
        SourceSpec::new(id(3),rev(1),ModelUnits::default(),
            SourceGeometry::FiniteCylinder{frame,radius:1.,axial:ParameterRange::bounded(-1.,3.).unwrap()},
            GeometryTolerance::default()).unwrap()
    }
    #[test]
    fn immutable_snapshots_version_and_insertion_order() {
        let old=SourceSnapshot::new([sample(5,1,0.),cylinder()]).unwrap();
        let same=SourceSnapshot::new([cylinder(),sample(5,1,0.)]).unwrap();
        assert_eq!(old.stamp(),same.stamp());
        let before=old.stamp().clone();
        let replacement=old.get(id(5)).unwrap().successor(box_geom(10.),
            GeometryTolerance::default()).unwrap();
        let next=old.with_replacement(replacement).unwrap();
        assert_eq!(old.stamp(),&before);
        assert_ne!(old.stamp(),next.stamp());
        assert_eq!(old.resolve(SourceRef{id:id(5),revision:rev(1)}).unwrap().revision(),rev(1));
        assert_eq!(next.resolve(SourceRef{id:id(5),revision:rev(1)}),
            Err(SourceError::StaleVersion{id:id(5),expected:rev(1),actual:rev(2)}));
        assert_eq!(next.get(id(5)).unwrap().revision(),rev(2));
    }
    #[test]
    fn reject_invalid_sources_and_updates() {
        assert_eq!(SourceId::new(0),Err(SourceError::InvalidId));
        assert_eq!(SourceRevision::new(0),Err(SourceError::InvalidRevision));
        assert!(SourceRevision::new(u64::MAX).unwrap().next().is_err());
        assert!(SourceSpec::new(id(1),rev(1),ModelUnits::default(),
            SourceGeometry::AxisAlignedBlock{min:p(0.,0.,0.),max:p(0.,1.,1.)},
            GeometryTolerance::default()).is_err());
        assert!(SourceSpec::new(id(1),rev(1),ModelUnits::default(),
            SourceGeometry::AxisAlignedBlock{min:p(f64::NAN,0.,0.),max:p(1.,1.,1.)},
            GeometryTolerance::default()).is_err());
        let a=sample(1,1,0.);
        let snap=SourceSnapshot::new([a.clone()]).unwrap();
        assert_eq!(snap.with_replacement(a).err(),Some(SourceError::RevisionNotSuccessor));
        assert!(snap.with_added(sample(1,2,0.)).is_err());
        assert!(SourceSnapshot::new([sample(2,1,0.),sample(2,1,0.)]).is_err());
        let mut bad=Frame3::from_axes(p(0.,0.,0.),p(1.,0.,0.),p(0.,1.,0.),
            GeometryTolerance::default()).unwrap();
        bad.normal=p(1.,0.,0.);
        assert!(SourceSpec::new(id(2),rev(1),ModelUnits::default(),
            SourceGeometry::FiniteCylinder{frame:bad,radius:1.,axial:ParameterRange::bounded(0.,2.).unwrap()},
            GeometryTolerance::default()).is_err());
    }
    #[test]
    fn binary_stamp_canonicalizes_signed_zero_and_source_order() {
        let a=sample(7,1,0.);
        let b=sample(7,1,-0.);
        assert_eq!(SourceSnapshot::new([a]).unwrap().stamp(),
                   SourceSnapshot::new([b]).unwrap().stamp());
    }
    #[test]
    fn intents_bind_exact_versions_units_and_boolean_kind() {
        let snap=SourceSnapshot::new([sample(1,1,0.),cylinder()]).unwrap();
        let a=snap.get(id(1)).unwrap().reference();
        let b=snap.get(id(3)).unwrap().reference();
        let tol=GeometryTolerance::default();
        let difference=snap.intent(BooleanKind::Difference,a,b,tol,UnitPolicy::RequireIdentical).unwrap();
        assert!(snap.check_intent(&difference).is_ok());
        let union=snap.intent(BooleanKind::Union,a,b,tol,UnitPolicy::RequireIdentical).unwrap();
        assert_ne!(difference.intent_stamp(),union.intent_stamp());
        let identity=snap.intent(BooleanKind::Difference,a,a,tol,UnitPolicy::RequireIdentical).unwrap();
        assert!(snap.check_intent(&identity).is_ok());
        let replaced=snap.with_replacement(snap.get(id(1)).unwrap().successor(box_geom(3.),tol).unwrap()).unwrap();
        assert_eq!(replaced.check_intent(&difference),Err(SourceError::SnapshotMismatch));
        assert!(replaced.intent(BooleanKind::Difference,a,b,tol,UnitPolicy::RequireIdentical).is_err());
    }
    #[test]
    fn unit_transform_must_be_recorded_and_correct() {
        let mm=sample(1,1,0.);
        let inch=SourceSpec::new(id(2),rev(1),
            ModelUnits{length:LengthUnit::Inch,angle:AngleUnit::Radian},
            box_geom(0.),GeometryTolerance::default()).unwrap();
        let snap=SourceSnapshot::new([mm.clone(),inch.clone()]).unwrap();
        let tol=GeometryTolerance::default();
        let a=mm.reference();let b=inch.reference();
        assert_eq!(snap.intent(BooleanKind::Union,a,b,tol,UnitPolicy::RequireIdentical).err(),
            Some(SourceError::UnitConversionRequired));
        assert!(snap.intent(BooleanKind::Union,a,b,tol,
            UnitPolicy::RecordedRightToLeft{factor:1.}).is_err());
        let intent=snap.intent(BooleanKind::Union,a,b,tol,
            UnitPolicy::RecordedRightToLeft{factor:25.4}).unwrap();
        assert!(snap.check_intent(&intent).is_ok());
    }
}
