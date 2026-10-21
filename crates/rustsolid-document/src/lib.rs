//! Strong-atomic, in-memory G8 document revision and candidate publication.
//! Persistent archive/replay and generalized feature identity are not in this slice.
use rustsolid_boolean::{CertifiedCandidate,PrepareOutcome,PreparationAnchor,
    FaceSourceWitness,prepare_g8_tilted};
use rustsolid_source::{SourceSnapshot,SourceSpec,SourceStamp,BooleanIntent,SourceError};
use rustsolid_topology::{BrepModel,FaceId};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64,Ordering};

static NEXT_DOC:AtomicU64=AtomicU64::new(1);

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum DocumentError {
    StaleRevision, ForeignDocument, StaleSource,
    InvalidCertificate, ExistingResult, InvalidSource,
    RevisionExhausted, HistoryAllocationFailed, EmptyHistory,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum Transition { Undo, Redo }
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash)]
pub struct FaceHandle {
    document_token:u128,revision:u64,face:FaceId,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum FaceHandleError {ForeignDocument,StaleHandle,MissingFace}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct ScopedFaceLineage {
    pub local_face:FaceId,
    pub source:rustsolid_source::SourceRef,
    pub role:rustsolid_boolean::FaceProvenance,
    /// Within a source role, ordinal follows deterministic B-rep face order.
    pub ordinal:u32,
}
/// G8 feature-scoped identity. This is intentionally not universal persistent
/// naming across arbitrary geometric topology changes.
#[derive(Debug,Clone,PartialEq,Eq)]
pub struct TopologyLineage {
    pub document_revision:u64,
    pub faces:Vec<ScopedFaceLineage>,
}
#[derive(Debug,Clone)]
struct DocumentState {
    sources:SourceSnapshot,
    committed:Option<Arc<CertifiedCandidate>>,
}
#[derive(Debug)]
pub struct ModelingDocument {
    token:u128,
    revision:u64,
    current:DocumentState,
    undo:Vec<DocumentState>,
    redo:Vec<DocumentState>,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CommitReceipt {
    pub before:u64,
    pub after:u64,
    pub certified_volume:f64,
    pub face_count:usize,
}
impl ModelingDocument {
    pub fn new(sources:SourceSnapshot)->Result<Self,DocumentError>{
        let token=NEXT_DOC.fetch_add(1,Ordering::Relaxed);
        if token==0||token==u64::MAX{return Err(DocumentError::RevisionExhausted);}
        Ok(Self{token:token as u128,revision:0,
            current:DocumentState{sources,committed:None},
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn token(&self)->u128{self.token}
    pub fn revision(&self)->u64{self.revision}
    pub fn snapshot(&self)->&SourceSnapshot{&self.current.sources}
    pub fn snapshot_stamp(&self)->&SourceStamp{self.current.sources.stamp()}
    pub fn model(&self)->Option<&BrepModel>{
        self.current.committed.as_ref().map(|c|c.model())
    }
    pub fn prepare_g8(
        &self,intent:&BooleanIntent,model:BrepModel,
    )->Result<PrepareOutcome,DocumentError>{
        if self.current.committed.is_some(){return Err(DocumentError::ExistingResult);}
        let anchor=PreparationAnchor::new(self.token,self.revision)
            .ok_or(DocumentError::RevisionExhausted)?;
        Ok(prepare_g8_tilted(&self.current.sources,intent,model,anchor))
    }
    fn next_revision(&self)->Result<u64,DocumentError>{
        self.revision.checked_add(1).ok_or(DocumentError::RevisionExhausted)
    }
    pub fn commit(&mut self,candidate:CertifiedCandidate)->Result<CommitReceipt,DocumentError>{
        let anchor=candidate.anchor();
        if anchor.document_token()!=self.token{return Err(DocumentError::ForeignDocument);}
        if anchor.expected_revision()!=self.revision{return Err(DocumentError::StaleRevision);}
        if self.current.committed.is_some(){return Err(DocumentError::ExistingResult);}
        if candidate.source_stamp()!=self.snapshot_stamp(){
            return Err(DocumentError::StaleSource);
        }
        // The whole model/source/flux/provenance must revalidate before mutation.
        candidate.recheck(&self.current.sources)
            .map_err(|_|DocumentError::InvalidCertificate)?;
        let after=self.next_revision()?;
        self.undo.try_reserve(1).map_err(|_|DocumentError::HistoryAllocationFailed)?;
        let receipt=CommitReceipt{
            before:self.revision,after,
            certified_volume:candidate.volume().measured_flux(),
            face_count:candidate.model().faces().len(),
        };
        let mut next=self.current.clone();
        next.committed=Some(Arc::new(candidate));
        // Publication barrier: no fallible operations after this point.
        let previous=std::mem::replace(&mut self.current,next);
        self.undo.push(previous);
        self.redo.clear();
        self.revision=after;
        Ok(receipt)
    }
    /// A source edit invalidates the accepted result and all runtime handles;
    /// the input source snapshot remains immutable to callers.
    pub fn replace_source(&mut self,next:SourceSpec)->Result<u64,DocumentError>{
        let next_sources=self.current.sources.with_replacement(next)
            .map_err(|_|DocumentError::InvalidSource)?;
        let revision=self.next_revision()?;
        self.undo.try_reserve(1).map_err(|_|DocumentError::HistoryAllocationFailed)?;
        let next_state=DocumentState{sources:next_sources,committed:None};
        let prev=std::mem::replace(&mut self.current,next_state);
        self.undo.push(prev);
        self.redo.clear();
        self.revision=revision;
        Ok(revision)
    }
    pub fn transition(&mut self,direction:Transition)->Result<u64,DocumentError>{
        let from=match direction{
            Transition::Undo=>&self.undo,
            Transition::Redo=>&self.redo,
        };
        let Some(next)=from.last()else{return Err(DocumentError::EmptyHistory);};
        if let Some(cert)=&next.committed{
            cert.recheck(&next.sources).map_err(|_|DocumentError::InvalidCertificate)?;
        }
        let after=self.next_revision()?;
        match direction {
            Transition::Undo=>self.redo.try_reserve(1),
            Transition::Redo=>self.undo.try_reserve(1),
        }.map_err(|_|DocumentError::HistoryAllocationFailed)?;
        // Infallible mutation barrier.
        let next=match direction{
            Transition::Undo=>self.undo.pop().expect("checked"),
            Transition::Redo=>self.redo.pop().expect("checked"),
        };
        let old=std::mem::replace(&mut self.current,next);
        match direction{Transition::Undo=>self.redo.push(old),Transition::Redo=>self.undo.push(old)}
        self.revision=after;
        Ok(after)
    }
    pub fn face_handle(&self,id:FaceId)->Result<FaceHandle,FaceHandleError>{
        let model=self.model().ok_or(FaceHandleError::MissingFace)?;
        if model.faces().get(id.0 as usize).map(|f|f.id)!=Some(id) {
            return Err(FaceHandleError::MissingFace);
        }
        Ok(FaceHandle{document_token:self.token,revision:self.revision,face:id})
    }
    pub fn resolve_face(&self,h:FaceHandle)->Result<FaceSourceWitness,FaceHandleError>{
        if h.document_token!=self.token{return Err(FaceHandleError::ForeignDocument);}
        if h.revision!=self.revision{return Err(FaceHandleError::StaleHandle);}
        let cand=self.current.committed.as_ref().ok_or(FaceHandleError::MissingFace)?;
        cand.face_lineage().iter().find(|w|w.result_face==h.face)
            .copied().ok_or(FaceHandleError::MissingFace)
    }
    pub fn lineage(&self)->Option<TopologyLineage>{
        let cand=self.current.committed.as_ref()?;
        let mut tally=[0u32;7];
        let mut entries=Vec::new();
        for f in cand.face_lineage() {
            let idx=match f.role{
                rustsolid_boolean::FaceProvenance::BoxXMin=>0,
                rustsolid_boolean::FaceProvenance::BoxXMax=>1,
                rustsolid_boolean::FaceProvenance::BoxYMin=>2,
                rustsolid_boolean::FaceProvenance::BoxYMax=>3,
                rustsolid_boolean::FaceProvenance::BoxZMin=>4,
                rustsolid_boolean::FaceProvenance::BoxZMax=>5,
                rustsolid_boolean::FaceProvenance::CylinderSide=>6,
            };
            let ordinal=tally[idx];
            tally[idx]+=1;
            entries.push(ScopedFaceLineage{
                local_face:f.result_face,source:f.source,role:f.role,ordinal
            });
        }
        Some(TopologyLineage{document_revision:self.revision,faces:entries})
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rustsolid_source::{SourceId,SourceRevision,SourceGeometry};
    use rustsolid_foundation::{Point3,GeometryTolerance};
    use rustsolid_geometry::ModelUnits;
    fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
    fn spec()->SourceSpec {
        SourceSpec::new(SourceId::new(1).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),SourceGeometry::AxisAlignedBlock{
                min:p(0.,0.,0.),max:p(3.,4.,5.)},GeometryTolerance::default()).unwrap()
    }
    #[test]
    fn source_replacement_undo_redo_are_atomic_and_revision_monotonic(){
        let initial=spec();
        let mut doc=ModelingDocument::new(SourceSnapshot::new([initial.clone()]).unwrap()).unwrap();
        let old=doc.snapshot_stamp().clone();
        let next=initial.successor(SourceGeometry::AxisAlignedBlock{
            min:p(0.,0.,0.),max:p(6.,4.,5.)},GeometryTolerance::default()).unwrap();
        assert_eq!(doc.replace_source(next).unwrap(),1);
        let changed=doc.snapshot_stamp().clone();
        assert_ne!(old,changed);
        assert!(doc.replace_source(initial.clone()).is_err());
        assert_eq!(doc.revision(),1);
        assert_eq!(doc.transition(Transition::Undo).unwrap(),2);
        assert_eq!(doc.snapshot_stamp(),&old);
        assert_eq!(doc.transition(Transition::Redo).unwrap(),3);
        assert_eq!(doc.snapshot_stamp(),&changed);
    }
}
