//! Narrow, sealed and source-bound G8 tilted one-wall candidate.
//! This path is allowed only for the complete eight-face strict source family.
//! GuardedF64 is domain-scoped evidence; not a general NURBS Boolean theorem.
use rustsolid_foundation::{GeometryTolerance,Point3};
use rustsolid_geometry::{Surface3,Frame3};
use rustsolid_source::{SourceSnapshot,BooleanIntent,BooleanKind,SourceGeometry,SourceRef,SourceStamp};
use rustsolid_topology::{BrepModel,FaceId,FaceSense};
use rustsolid_validation::{evaluate_tilted_one_wall,evaluate_tilted_boundary_flux,
    compare_g8_numeric,VolumeAgreement,AssuranceLevel};

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum FaceProvenance {
    BoxXMin,BoxXMax,BoxYMin,BoxYMax,BoxZMin,BoxZMax,CylinderSide,
}
impl FaceProvenance {
    fn index(self)->usize{match self{
        Self::BoxXMin=>0,Self::BoxXMax=>1,Self::BoxYMin=>2,Self::BoxYMax=>3,
        Self::BoxZMin=>4,Self::BoxZMax=>5,Self::CylinderSide=>6
    }}
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct FaceSourceWitness {
    pub result_face:FaceId,
    pub source:SourceRef,
    pub role:FaceProvenance,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum HoldReason {
    InvalidSnapshot,InvalidSourceFamily,InvalidCandidateTopology,
    InvalidSourceSurface,UnresolvedMaterialSense,
    InvalidProvenance,IncompleteVolumeEvidence,
}
#[derive(Debug)]
pub enum PrepareOutcome {
    Ready(CertifiedCandidate),
    Held(HoldReason),
    Unsupported,
}
/// Runtime document anchor; forged anchors cannot bypass Document::commit
/// which always compares its private current incarnation and revision.
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct PreparationAnchor {
    document_token:u128,
    expected_revision:u64,
}
impl PreparationAnchor {
    pub fn new(document_token:u128,expected_revision:u64)->Option<Self>{
        if document_token==0{return None;}
        Some(Self{document_token,expected_revision})
    }
    pub fn document_token(self)->u128{self.document_token}
    pub fn expected_revision(self)->u64{self.expected_revision}
}
#[derive(Debug)]
pub struct CertifiedCandidate {
    model:BrepModel,
    source_stamp:SourceStamp,
    intent: BooleanIntent,
    anchor:PreparationAnchor,
    volume:VolumeAgreement,
    lineage:Vec<FaceSourceWitness>,
}
impl CertifiedCandidate {
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn source_stamp(&self)->&SourceStamp{&self.source_stamp}
    pub fn intent(&self)->&BooleanIntent{&self.intent}
    pub fn anchor(&self)->PreparationAnchor{self.anchor}
    pub fn volume(&self)->&VolumeAgreement{&self.volume}
    pub fn face_lineage(&self)->&[FaceSourceWitness]{&self.lineage}
    /// Re-run ALL source, analytic binding, actual flux and material support
    /// checks on the immutable current source snapshot.
    pub fn recheck(&self,snapshot:&SourceSnapshot)->Result<(),HoldReason>{
        let PrepareOutcome::Ready(fresh)=prepare_g8_tilted(
            snapshot,&self.intent,self.model.clone(),self.anchor)
        else{return Err(HoldReason::InvalidSnapshot)};
        if fresh.source_stamp!=self.source_stamp ||
            fresh.lineage!=self.lineage ||
            fresh.volume.retained_interval()!=self.volume.retained_interval() ||
            (fresh.volume.measured_flux()-self.volume.measured_flux()).abs()>1e-10 ||
            fresh.volume.assurance()!=self.volume.assurance() {
            return Err(HoldReason::IncompleteVolumeEvidence);
        }
        Ok(())
    }
}
fn dist(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn matches_frame(a:Frame3,b:Frame3,epsilon:f64)->bool{
    dist(a.origin,b.origin)<=epsilon &&
        dist(a.x,b.x)<=1e-9 &&dist(a.y,b.y)<=1e-9 &&
        dist(a.normal,b.normal)<=1e-9
}
fn source_faces(
    m:&BrepModel,box_min:Point3,box_max:Point3,
    cutter:Frame3,radius:f64,left:SourceRef,right:SourceRef,
    tol:GeometryTolerance,
)->Result<Vec<FaceSourceWitness>,HoldReason>{
    m.validate().map_err(|_|HoldReason::InvalidCandidateTopology)?;
    if m.faces().len()!=8||m.shells().len()!=1||m.vertices().len()!=12
        ||m.edges().len()!=18||m.coedges().len()!=36||m.loops().len()!=8{
        return Err(HoldReason::InvalidCandidateTopology);
    }
    if m.vertices().len() as isize -m.edges().len() as isize
        +m.faces().len() as isize != 2 {
        return Err(HoldReason::InvalidCandidateTopology);
    }
    let scale=dist(box_min,box_max).max(radius).max(1.);
    let eps=tol.length_at(scale)*128.;
    for v in m.vertices(){
        let p=v.point;
        if p.x<box_min.x-eps||p.x>box_max.x+eps||
            p.y<box_min.y-eps||p.y>box_max.y+eps||
            p.z<box_min.z-eps||p.z>box_max.z+eps{
            return Err(HoldReason::InvalidSourceSurface);
        }
    }
    let mut counts=[0usize;7];
    let mut mapped=Vec::new();
    for face in m.faces(){
        let support=m.carriers().surfaces[face.surface.0 as usize];
        let (role,owner)=match support{
            Surface3::Cylinder(actual)=>{
                if face.sense!=FaceSense::Reversed ||
                    (actual.radius-radius).abs()>eps ||
                    !matches_frame(actual.frame,cutter,eps) {
                    return Err(HoldReason::InvalidSourceSurface);
                }
                (FaceProvenance::CylinderSide,right)
            },
            Surface3::Plane(s)=>{
                let mut normal=s.frame.normal;
                if face.sense==FaceSense::Reversed{normal=normal.scale(-1.);}
                let p=s.frame.origin;
                let intended=[
                    (FaceProvenance::BoxXMin,p.x-box_min.x,Point3{x:-1.,y:0.,z:0.}),
                    (FaceProvenance::BoxXMax,p.x-box_max.x,Point3{x:1.,y:0.,z:0.}),
                    (FaceProvenance::BoxYMin,p.y-box_min.y,Point3{x:0.,y:-1.,z:0.}),
                    (FaceProvenance::BoxYMax,p.y-box_max.y,Point3{x:0.,y:1.,z:0.}),
                    (FaceProvenance::BoxZMin,p.z-box_min.z,Point3{x:0.,y:0.,z:-1.}),
                    (FaceProvenance::BoxZMax,p.z-box_max.z,Point3{x:0.,y:0.,z:1.}),
                ];
                let Some((role,_,_))=intended.into_iter()
                    .find(|(_,offset,direction)|
                        offset.abs()<=eps &&dist(normal,*direction)<=1e-8)
                else{return Err(HoldReason::UnresolvedMaterialSense)};
                (role,left)
            },
            Surface3::Sphere(_)=>return Err(HoldReason::InvalidSourceSurface),
        };
        counts[role.index()]+=1;
        mapped.push(FaceSourceWitness{result_face:face.id,source:owner,role});
    }
    if counts!=[2,1,1,1,1,1,1]{
        return Err(HoldReason::InvalidProvenance);
    }
    Ok(mapped)
}
/// Only complete G8 tilted one-wall results can become Ready here.
/// No caller-supplied numeric/material proofs are accepted.
pub fn prepare_g8_tilted(
    snapshot:&SourceSnapshot,intent:&BooleanIntent,
    model:BrepModel,anchor:PreparationAnchor,
)->PrepareOutcome{
    if snapshot.check_intent(intent).is_err(){
        return PrepareOutcome::Held(HoldReason::InvalidSnapshot);
    }
    if intent.operation()!=BooleanKind::Difference{
        return PrepareOutcome::Unsupported;
    }
    let (left,right)=intent.operands();
    let (Ok(a),Ok(b))=(snapshot.resolve(left),snapshot.resolve(right))
        else{return PrepareOutcome::Held(HoldReason::InvalidSnapshot)};
    let (SourceGeometry::AxisAlignedBlock{min,max},
         SourceGeometry::FiniteCylinder{frame,radius,..})=(a.geometry(),b.geometry())
        else{return PrepareOutcome::Unsupported};
    if a.units()!=b.units() ||model.units()!=a.units(){
        return PrepareOutcome::Held(HoldReason::InvalidSourceFamily);
    }
    let source=match evaluate_tilted_one_wall(snapshot,intent){
        Ok(s)=>s,Err(_)=>return PrepareOutcome::Held(HoldReason::InvalidSourceFamily)
    };
    let lineage=match source_faces(&model,*min,*max,*frame,*radius,left,right,intent.tolerance()){
        Ok(x)=>x,Err(why)=>return PrepareOutcome::Held(why)
    };
    let flux=match evaluate_tilted_boundary_flux(&model,min.x,intent.tolerance()){
        Ok(x)=>x,Err(_)=>return PrepareOutcome::Held(HoldReason::IncompleteVolumeEvidence)
    };
    let volume=match compare_g8_numeric(snapshot,intent,&source,&flux){
        Ok(v)=>v,Err(_)=>return PrepareOutcome::Held(HoldReason::IncompleteVolumeEvidence)
    };
    if volume.assurance()!=AssuranceLevel::GuardedF64{
        return PrepareOutcome::Held(HoldReason::IncompleteVolumeEvidence);
    }
    PrepareOutcome::Ready(CertifiedCandidate{
        model,source_stamp:snapshot.stamp().clone(),
        intent:intent.clone(),anchor,volume,lineage,
    })
}
