//! Direct G7 MIN-X analytic one-wall Boolean geometry construction.
//! Constructs topology::BrepDraft from immutable source parameters.
//! The resulting G7StagedBrep is TOPOLOGICALLY checked, but cannot be committed
//! as a CertifiedCandidate without independent volume and material/lineage proof.

use std::collections::BTreeMap;
use rustsolid_foundation::{GeometryTolerance,Point3};
use rustsolid_geometry::{
    Circle3,Curve2,Curve3,CylinderSection2,Frame3,HarmonicCurve2,
    Line2,Line3,ModelUnits,ParameterRange,PlaneSurface,Point2Param,
    Surface3,CylinderSurface,
};
use rustsolid_source::{
    BooleanIntent,BooleanKind,SourceError,SourceGeometry,
    SourceSnapshot,SourceStamp,UnitPolicy,
};
use rustsolid_topology::{
    AffineParameterMap,BrepCoedge,BrepDraft,BrepEdge,BrepFace,
    BrepLoop,BrepModel,BrepShell,BrepVertex,CoedgeId,Curve2Id,
    Curve3Id,EdgeId,EdgeKind,FaceId,FaceSense,GeometryCatalog,
    LocalBrepIncarnation,LoopId,LoopRole,ShellId,ShellSense,SurfaceId,VertexId,
};

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum G7BuildHold{
    UnsupportedOperation,
    UnsupportedSourceFamily,
    MixedUnits,
    CylinderCapsInMaterial,
    CylinderOrientationUnsupported,
    NotStrictOneWallNotch,
    AmbiguousCircularCut,
    TopologyAssemblyFailed,
}
#[derive(Debug,Clone)]
pub struct G7StagedBrep {
    model:BrepModel,
    source_stamp:SourceStamp,
    intent_stamp:SourceStamp,
    major_arc:[f64;2],
}
impl G7StagedBrep {
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn source_stamp(&self)->&SourceStamp{&self.source_stamp}
    pub fn intent_stamp(&self)->&SourceStamp{&self.intent_stamp}
    pub fn major_arc(&self)->[f64;2]{self.major_arc}
}
#[derive(Debug,Clone)]
pub enum G7BuildOutcome {
    Staged(G7StagedBrep),
    Held(G7BuildHold),
}
fn point(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn size(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn distance(a:Point3,b:Point3)->f64{size(a.sub(b))}

pub fn g7_one_wall_direct(
    snapshot:&SourceSnapshot,intent:&BooleanIntent,
)->Result<G7BuildOutcome,SourceError>{
    use G7BuildHold as Held;
    snapshot.check_intent(intent)?;
    if intent.operation()!=BooleanKind::Difference{
        return Ok(G7BuildOutcome::Held(Held::UnsupportedOperation));
    }
    let (left,right)=intent.operands();
    let a=snapshot.resolve(left)?;
    let b=snapshot.resolve(right)?;
    if a.units()!=b.units() || !matches!(intent.unit_policy(),UnitPolicy::RequireIdentical) {
        return Ok(G7BuildOutcome::Held(Held::MixedUnits));
    }
    let (
        SourceGeometry::AxisAlignedBlock{min,max},
        SourceGeometry::FiniteCylinder{frame,radius,axial},
    )=(a.geometry(),b.geometry()) else {
        return Ok(G7BuildOutcome::Held(Held::UnsupportedSourceFamily));
    };
    let (ParameterRange::Bounded{start,end})=*axial else{
        return Ok(G7BuildOutcome::Held(Held::UnsupportedSourceFamily));
    };
    let (min,max,frame,radius)=(*min,*max,*frame,*radius);
    let tol=intent.tolerance();
    let ext=max.sub(min);
    let scale=size(ext).max(radius).max(end-start).max(1.);
    let eps=tol.length_at(scale)*128.;
    if distance(frame.x,point(1.,0.,0.)).max(
        distance(frame.y,point(0.,0.,-1.))).max(
        distance(frame.normal,point(0.,1.,0.)))>tol.angular*16. {
        return Ok(G7BuildOutcome::Held(Held::CylinderOrientationUnsupported));
    }
    if min.y-frame.origin.y<=start+eps ||
        max.y-frame.origin.y>=end-eps {
        return Ok(G7BuildOutcome::Held(Held::CylinderCapsInMaterial));
    }
    let cx=frame.origin.x;let cz=frame.origin.z;
    let dx=cx-min.x;
    if dx<=eps ||dx>=radius-eps ||cx+radius>=max.x-eps ||
        cz-radius<=min.z+eps ||cz+radius>=max.z-eps {
        return Ok(G7BuildOutcome::Held(Held::NotStrictOneWallNotch));
    }
    let q=(radius*radius-dx*dx).sqrt();
    if !q.is_finite()||q<=eps || cz+q>=max.z-eps ||cz-q<=min.z+eps{
        return Ok(G7BuildOutcome::Held(Held::AmbiguousCircularCut));
    }
    let upper=(-q).atan2(-dx);
    let lower=q.atan2(-dx);
    let span=lower-upper;
    if !(span>std::f64::consts::PI+tol.angular &&
         span<std::f64::consts::TAU-tol.angular){
        return Ok(G7BuildOutcome::Held(Held::AmbiguousCircularCut));
    }
    let positions=[
        point(min.x,min.y,min.z),
        point(max.x,min.y,min.z),
        point(max.x,min.y,max.z),
        point(min.x,min.y,max.z),
        point(min.x,min.y,cz+q),
        point(min.x,min.y,cz-q),
        point(min.x,max.y,min.z),
        point(max.x,max.y,min.z),
        point(max.x,max.y,max.z),
        point(min.x,max.y,max.z),
        point(min.x,max.y,cz+q),
        point(min.x,max.y,cz-q),
    ];
    let vertices=positions.iter().enumerate().map(|(i,&p)|
        BrepVertex{id:VertexId(i as u32),point:p}).collect::<Vec<_>>();
    // The canonical 3D curve direction is stable irrespective of face use.
    let mut catalog=GeometryCatalog::default();
    let edges_def=[
        (0usize,1usize),(1,2),(2,3),(3,4),(4,5),(5,0),
        (6,7),(7,8),(8,9),(9,10),(10,11),(11,6),
        (0,6),(1,7),(2,8),(3,9),(4,10),(5,11),
    ];
    let mut edge_lookup=BTreeMap::new();
    let mut edges=Vec::new();
    for &(v0,v1) in &edges_def {
        let id=EdgeId(edges.len() as u32);
        let (curve,domain)=if (v0,v1)==(4,5)||(v0,v1)==(10,11) {
            let center=point(cx,positions[v0].y,cz);
            let f=Frame3{origin:center,..frame};
            (Curve3::Circle(Circle3{frame:f,radius}),
                ParameterRange::Bounded{start:upper,end:lower})
        }else{
            let Ok((line,range))=Line3::between(positions[v0],positions[v1],tol) else {
                return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
            };
            (Curve3::Line(line),range)
        };
        let cid=Curve3Id(catalog.curves3.len() as u32);
        catalog.curves3.push(curve);
        edges.push(BrepEdge{
            id,vertices:[VertexId(v0 as u32),VertexId(v1 as u32)],
            curve:cid,domain,uses:[CoedgeId(u32::MAX);2],kind:EdgeKind::Physical,
        });
        edge_lookup.insert((v0.min(v1),v0.max(v1)),id);
    }
    // Each polygonal face contributes physical coedges. The two caps include
    // one genuine circular major arc; the cylinder side has the reversed arc.
    let rings: [&[usize];8]=[
        &[0,1,2,3,4,5],
        &[6,11,10,9,8,7],
        &[1,0,6,7],
        &[2,1,7,8],
        &[3,2,8,9],
        &[4,3,9,10],
        &[0,5,11,6],
        &[5,4,10,11],
    ];
    let mut faces=Vec::new();let mut loops=Vec::new();let mut coedges=Vec::new();
    let mut edge_use_count=vec![0usize;edges.len()];
    for (face_index,ring) in rings.iter().enumerate(){
        let fid=FaceId(face_index as u32);
        let lid=LoopId(face_index as u32);
        let sid=SurfaceId(catalog.surfaces.len() as u32);
        let support=if face_index==7 {
            Surface3::Cylinder(CylinderSurface{frame,radius})
        }else{
            let first=positions[ring[0]];
            let x_axis=positions[ring[1]].sub(first);
            let y_axis=positions[ring[2]].sub(first);
            let Ok(plane_frame)=Frame3::from_axes(first,x_axis,y_axis,tol) else{
                return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
            };
            Surface3::Plane(PlaneSurface{frame:plane_frame})
        };
        catalog.surfaces.push(support);
        faces.push(BrepFace{
            id:fid,surface:sid,loops:vec![lid],
            sense:if face_index==7{FaceSense::Reversed}else{FaceSense::Forward},
        });
        let first=CoedgeId(coedges.len() as u32);
        loops.push(BrepLoop{id:lid,face:fid,role:LoopRole::Outer,first});
        for (i,&from) in ring.iter().enumerate(){
            let to=ring[(i+1)%ring.len()];
            let eid=*edge_lookup.get(&(from.min(to),from.max(to)))
                .expect("all source ring edges are present");
            let e=&edges[eid.0 as usize];
            let from_vertex=e.vertices[0].0 as usize;
            let to_vertex=e.vertices[1].0 as usize;
            let reversed=from!=from_vertex ||to!=to_vertex;
            let is_circle=matches!(catalog.curves3[e.curve.0 as usize],Curve3::Circle(_));
            let (pcurve,map)=match support {
                Surface3::Plane(plane)=>{
                    if is_circle{
                        let center=plane.frame.project(point(cx,positions[from].y,cz))
                            .map_err(SourceError::Geometry)?;
                        let cosine=Point2Param{
                            u:radius*frame.x.dot(plane.frame.x),
                            v:radius*frame.x.dot(plane.frame.y),
                        };
                        let sine=Point2Param{
                            u:radius*frame.y.dot(plane.frame.x),
                            v:radius*frame.y.dot(plane.frame.y),
                        };
                        (Curve2::Harmonic(HarmonicCurve2{center,cosine,sine}),
                            AffineParameterMap{scale:1.,offset:0.})
                    }else{
                        let u=plane.frame.project(positions[from_vertex])
                            .map_err(SourceError::Geometry)?;
                        let v=plane.frame.project(positions[to_vertex])
                            .map_err(SourceError::Geometry)?;
                        let (ParameterRange::Bounded{start,end})=e.domain else{
                            unreachable!("all source edges bounded");
                        };
                        let t=end-start;
                        (Curve2::Line(Line2{origin:u,delta:Point2Param{
                            u:v.u-u.u,v:v.v-u.v,
                        }}),AffineParameterMap{scale:1./t,offset:-start/t})
                    }
                }
                Surface3::Cylinder(_)=>{
                    if is_circle {
                        let axial=positions[from_vertex].y-frame.origin.y;
                        (Curve2::CylinderSection(CylinderSection2{
                            angle_offset:0.,direction:1.,
                            axial_center:axial,axial_cosine:0.,axial_sine:0.,
                        }),AffineParameterMap{scale:1.,offset:0.})
                    }else{
                        // The only straight cylinder-side curves are true
                        // generators at the physical circle-wall events.
                        if (from_vertex,to_vertex)!=(4,10) &&
                           (from_vertex,to_vertex)!=(5,11) {
                            return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
                        }
                        let angle=if from_vertex==4 {upper}else{lower};
                        let v0=positions[from_vertex].y-frame.origin.y;
                        let v1=positions[to_vertex].y-frame.origin.y;
                        let (ParameterRange::Bounded{start,end})=e.domain else{
                            unreachable!("bounded source generator");
                        };
                        let t=end-start;
                        (Curve2::Line(Line2{
                            origin:Point2Param{u:angle,v:v0},
                            delta:Point2Param{u:0.,v:v1-v0},
                        }),AffineParameterMap{scale:1./t,offset:-start/t})
                    }
                }
                Surface3::Sphere(_)=>unreachable!(),
            };
            let pcurve_id=Curve2Id(catalog.curves2.len() as u32);
            catalog.curves2.push(pcurve);
            let cid=CoedgeId(coedges.len() as u32);
            let first_index=first.0;
            coedges.push(BrepCoedge{
                id:cid,edge:eid,face:fid,loop_id:lid,pcurve:pcurve_id,
                edge_to_pcurve:map,reversed,
                next:CoedgeId(first_index+((i+1)%ring.len()) as u32),
                prev:CoedgeId(first_index+((i+ring.len()-1)%ring.len()) as u32),
                twin:CoedgeId(u32::MAX),
            });
            let count=&mut edge_use_count[eid.0 as usize];
            if *count>=2 {
                return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
            }
            edges[eid.0 as usize].uses[*count]=cid;
            *count+=1;
        }
    }
    if edge_use_count.iter().any(|&n|n!=2) {
        return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
    }
    for edge in &edges {
        let [first,second]=edge.uses;
        coedges[first.0 as usize].twin=second;
        coedges[second.0 as usize].twin=first;
    }
    let draft=BrepDraft{
        incarnation:LocalBrepIncarnation::new(1)
            .expect("nonzero local B-rep incarnation"),
        units:a.units(),tolerance:tol,vertices,edges,coedges,
        loops,faces,shells:vec![BrepShell{
            id:ShellId(0),faces:(0..8).map(FaceId).collect(),
            sense:ShellSense::Outward,
        }],carriers:catalog,
    };
    let model=match draft.check(){
        Ok(model)=>model,
        Err(issue)=>{
            #[cfg(test)]
            eprintln!("G7 direct canonical draft rejected: {issue:?}");
            return Ok(G7BuildOutcome::Held(Held::TopologyAssemblyFailed));
        }
    };
    Ok(G7BuildOutcome::Staged(G7StagedBrep{
        model,source_stamp:snapshot.stamp().clone(),
        intent_stamp:intent.intent_stamp().clone(),
        major_arc:[upper,lower],
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsolid_source::{SourceId,SourceRevision,SourceSpec};
    fn fixture(radius:f64)->(SourceSnapshot,BooleanIntent){
        let tol=GeometryTolerance::default();
        let a=SourceSpec::new(SourceId::new(1).unwrap(),
            SourceRevision::new(1).unwrap(),ModelUnits::default(),
            SourceGeometry::AxisAlignedBlock{
                min:point(-1.,1.,-3.),max:point(3.,4.,3.),
            },tol).unwrap();
        let frame=Frame3::from_axes(point(0.,0.,0.),point(1.,0.,0.),
            point(0.,0.,-1.),tol).unwrap();
        let b=SourceSpec::new(SourceId::new(2).unwrap(),
            SourceRevision::new(1).unwrap(),ModelUnits::default(),
            SourceGeometry::FiniteCylinder{
                frame,radius,axial:ParameterRange::bounded(0.,5.).unwrap(),
            },tol).unwrap();
        let snapshot=SourceSnapshot::new([a.clone(),b.clone()]).unwrap();
        let intent=snapshot.intent(BooleanKind::Difference,a.reference(),
            b.reference(),tol,UnitPolicy::RequireIdentical).unwrap();
        (snapshot,intent)
    }
    #[test]
    fn one_wall_analytic_direct_builder_has_expected_physical_counts(){
        let (s,intent)=fixture(1.5);
        let out=g7_one_wall_direct(&s,&intent).unwrap();
        let G7BuildOutcome::Staged(stage)=out else{
            panic!("source-only G7 direct builder unexpectedly held: {out:?}");
        };
        let m=stage.model();
        assert_eq!((m.vertices().len(),m.edges().len(),m.coedges().len(),
            m.faces().len(),m.loops().len()),(12,18,36,8,8));
        assert!(m.validate().is_ok());
        assert!(stage.major_arc()[1]-stage.major_arc()[0]>std::f64::consts::PI);
        assert_eq!(m.carriers().curves3.iter().filter(|c|matches!(c,Curve3::Circle(_))).count(),2);
        assert_eq!(stage.source_stamp(),s.stamp());
    }
    #[test]
    fn one_wall_tangency_is_rejected_without_a_brep(){
        let (s,i)=fixture(1.0);
        assert!(matches!(g7_one_wall_direct(&s,&i).unwrap(),
            G7BuildOutcome::Held(G7BuildHold::NotStrictOneWallNotch)));
    }
}
