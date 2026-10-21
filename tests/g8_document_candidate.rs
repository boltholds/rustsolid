//! Real G8 candidate staging and atomic document acceptance regression.
use rustsolid::{BrepBody,BooleanOperation,CylindricalBrep,Frame3,GeometryTolerance,
    IndependentTiltedVolumeOracle,IndependentTiltedVolumeResult,
    Point3,Solid,TiltedOneWallSourcePlan,TiltedOneWallSourceResult,
    TiltedOneWallShellPlan,TiltedOneWallShellResult};
use rustsolid_source::{BooleanKind,BooleanIntent,SourceGeometry,SourceId,SourceRevision,
    SourceSnapshot,SourceSpec,UnitPolicy};
use rustsolid_topology::{BrepModel,FaceId};
use rustsolid_geometry::{ModelUnits,ParameterRange};
use rustsolid_boolean::PrepareOutcome;
use rustsolid_document::{ModelingDocument,DocumentError,Transition,FaceHandleError};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn fixture()->(SourceSnapshot,BooleanIntent,BrepModel){
    let angle=0.04_f64;
    let frame=Frame3::from_axes(p(0.,-2.,0.),p(angle.cos(),-angle.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let first=BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
    let second=BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap());
    let TiltedOneWallSourceResult::Ready(source)=TiltedOneWallSourcePlan::build(
        &first,&second,BooleanOperation::Difference,tol()).unwrap() else{panic!("G8a held")};
    let IndependentTiltedVolumeResult::Ready(volume)=IndependentTiltedVolumeOracle::build(
        &first,&second,&source,tol()).unwrap() else{panic!("G8b held")};
    let TiltedOneWallShellResult::Ready(result)=TiltedOneWallShellPlan::build(
        &first,&second,&source,&volume,tol()).unwrap() else{panic!("G8c held")};
    let model=result.checked_canonical().unwrap();
    let box_source=SourceSpec::new(SourceId::new(1).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),SourceGeometry::AxisAlignedBlock{
            min:p(-1.,1.,-3.),max:p(3.,4.,3.)},tol()).unwrap();
    let cyl_source=SourceSpec::new(SourceId::new(2).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),SourceGeometry::FiniteCylinder{frame,radius:1.5,
            axial:ParameterRange::bounded(0.,10.).unwrap()},tol()).unwrap();
    let snapshot=SourceSnapshot::new([box_source.clone(),cyl_source.clone()]).unwrap();
    let intent=snapshot.intent(BooleanKind::Difference,box_source.reference(),
        cyl_source.reference(),tol(),UnitPolicy::RequireIdentical).unwrap();
    (snapshot,intent,model)
}
#[test]
fn verified_g8_staging_commits_once_and_preserves_scoped_face_lineage(){
    let (snapshot,intent,model)=fixture();
    let mut doc=ModelingDocument::new(snapshot).unwrap();
    let PrepareOutcome::Ready(candidate)=doc.prepare_g8(&intent,model).unwrap()
        else{panic!("G8 source/topology/flux candidate held")};
    assert_eq!(candidate.face_lineage().len(),8);
    assert_eq!(candidate.volume().assurance(),rustsolid_validation::AssuranceLevel::GuardedF64);
    let receipt=doc.commit(candidate).unwrap();
    assert_eq!((receipt.before,receipt.after,receipt.face_count),(0,1,8));
    assert_eq!(doc.model().unwrap().faces().len(),8);
    let lineage=doc.lineage().unwrap();
    assert_eq!((lineage.document_revision,lineage.faces.len()),(1,8));
    let old=doc.face_handle(FaceId(7)).unwrap();
    assert_eq!(doc.resolve_face(old).unwrap().result_face,FaceId(7));
    assert_eq!(doc.transition(Transition::Undo).unwrap(),2);
    assert!(doc.model().is_none());
    assert_eq!(doc.transition(Transition::Redo).unwrap(),3);
    assert!(doc.model().is_some());
    assert_eq!(doc.resolve_face(old),Err(FaceHandleError::StaleHandle));
}
#[test]
fn foreign_and_stale_candidates_cannot_mutate_a_document(){
    let (snapshot,intent,model)=fixture();
    let first=ModelingDocument::new(snapshot.clone()).unwrap();
    let mut second=ModelingDocument::new(snapshot).unwrap();
    let PrepareOutcome::Ready(foreign)=first.prepare_g8(&intent,model.clone()).unwrap()
        else{panic!("G8 should be accepted")};
    assert_eq!(second.commit(foreign).unwrap_err(),DocumentError::ForeignDocument);
    assert_eq!(second.revision(),0);
    assert!(second.model().is_none());
    let PrepareOutcome::Ready(stale)=second.prepare_g8(&intent,model).unwrap()
        else{panic!("G8 should be accepted")};
    let original=second.snapshot().get(SourceId::new(1).unwrap()).unwrap();
    let updated=original.successor(SourceGeometry::AxisAlignedBlock{
        min:p(-1.,1.,-3.),max:p(4.,4.,3.)},tol()).unwrap();
    assert_eq!(second.replace_source(updated).unwrap(),1);
    assert_eq!(second.commit(stale).unwrap_err(),DocumentError::StaleRevision);
    assert_eq!(second.revision(),1);
    assert!(second.model().is_none());
}
#[test]
fn union_cannot_be_promoted_to_tilted_difference_candidate(){
    let (snapshot,intent,model)=fixture();
    let document=ModelingDocument::new(snapshot.clone()).unwrap();
    let (l,r)=intent.operands();
    let union=snapshot.intent(BooleanKind::Union,l,r,tol(),UnitPolicy::RequireIdentical).unwrap();
    assert!(matches!(document.prepare_g8(&union,model).unwrap(),PrepareOutcome::Unsupported));
    assert_eq!(document.revision(),0);
}
