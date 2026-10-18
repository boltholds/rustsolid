//! G8a: source-certified analytic preview for a tilted one-wall notch.
//!
//! This is deliberately NOT an executable Boolean or a closed result shell.
//! It derives finite horizontal ellipse sections from the original cylinder
//! frame, proves that every section crosses only MIN-X (with conservative
//! separation from other walls), and audits four physical G1 cut arcs.
//! Source-face UV and 3D carriers remain analytic and are independently
//! evaluated against both original operands.

use std::collections::{BTreeMap,BTreeSet};
use crate::{
    AnalyticCutRegistryHold, BooleanError, BooleanOperation, BrepBody,
    Curve2, Curve3, FaceId, GeometryError, GeometryTolerance,
    GlobalAnalyticCutPlan, GlobalAnalyticCutResult, GlobalArcRef,
    GlobalEventId, GlobalEventKind, Point3, Surface3,
};
use crate::analytic_face_patch_partition::find_source_plane;

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b); d.x.hypot(d.y).hypot(d.z)
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TiltedOneWallHold {
    UnsupportedOperation,
    UnsupportedOperand,
    NonCanonicalBox,
    NotPositiveYTransverseTilt,
    CutterCapsReachBox,
    NotIsolatedMinXWall,
    AmbiguousWallCrossing,
    MissingSourceFace,
    UpstreamCuts(AnalyticCutRegistryHold),
    UnexpectedSourceCuts,
    UnreconciledPhysicalCycle,
    SourceCarrierMismatch,
}

/// Horizontal section of the infinite cutter side at one original box cap.
/// P(u)=center+cosine*cos(u)+sine*sin(u), with zero Y amplitudes.
/// Intersections with MIN-X occur at angles [lo, hi] and the retained
/// x>=xmin ellipse arc has angular span hi-lo > PI.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TiltedOneWallSection {
    pub height_y:f64,
    pub center:Point3,
    pub cosine:Point3,
    pub sine:Point3,
    pub wall_angles:[f64;2],
    pub wall_points:[Point3;2],
    pub retained_major_span:f64,
}
impl TiltedOneWallSection {
    pub fn evaluate(&self,angle:f64)->Result<Point3,GeometryError>{
        if !angle.is_finite(){return Err(invalid("nonfinite G8 section parameter"));}
        let (s,c)=angle.sin_cos();
        let q=self.center.add(self.cosine.scale(c)).add(self.sine.scale(s));
        if !q.is_finite(){return Err(invalid("nonfinite G8 section evaluation"));}
        Ok(q)
    }
    fn build(
        origin:Point3, x_axis:Point3, y_axis:Point3, axis:Point3,
        radius:f64, level:f64, box_bounds:crate::BoundingBox,
        margin:f64, angular:f64,
    )->Option<Self>{
        // Axis has a positive nonzero Y component by the caller's contract.
        let inverse_y=1.0/axis.y;
        let shift=(level-origin.y)*inverse_y;
        let center=origin.add(axis.scale(shift));
        let cosine=x_axis.sub(axis.scale(x_axis.y*inverse_y)).scale(radius);
        let sine=y_axis.sub(axis.scale(y_axis.y*inverse_y)).scale(radius);
        if !center.is_finite()||!cosine.is_finite()||!sine.is_finite(){
            return None;
        }
        let ax=cosine.x.hypot(sine.x);
        let az=cosine.z.hypot(sine.z);
        let xmin=box_bounds.min.x;let xmax=box_bounds.max.x;
        let zmin=box_bounds.min.z;let zmax=box_bounds.max.z;
        if !ax.is_finite()||!az.is_finite()||ax<=margin ||
            center.x-xmin<=margin ||
            center.x-ax>=xmin-margin ||
            center.x+ax>=xmax-margin ||
            center.z-az<=zmin+margin ||
            center.z+az>=zmax-margin
        {return None;}
        // Analytic line/ellipse events; no root-search or sampled topology.
        let ratio=(xmin-center.x)/ax;
        if !ratio.is_finite()||ratio<=-1.0+angular*32.0||
            ratio>=-angular*32.0 {return None;}
        let phase=sine.x.atan2(cosine.x);
        let delta=ratio.acos();
        let lo=phase-delta;let hi=phase+delta;
        if hi-lo<=std::f64::consts::PI+angular*64.0 ||
            hi-lo>=std::f64::consts::TAU-angular*64.0{return None;}
        let mut section=Self{
            height_y:level,center,cosine,sine,
            wall_angles:[lo,hi],wall_points:[origin;2],
            retained_major_span:hi-lo,
        };
        let a=section.evaluate(lo).ok()?;
        let b=section.evaluate(hi).ok()?;
        if (a.x-xmin).abs()>margin/8.0 ||
            (b.x-xmin).abs()>margin/8.0 ||
            (a.y-level).abs()>margin/8.0 ||
            (b.y-level).abs()>margin/8.0 ||
            (a.z-b.z).abs()<=margin*2.0 ||
            [a.z,b.z].iter().any(|z|*z<=zmin+margin || *z>=zmax-margin){
            return None;
        }
        section.wall_points=[a,b];
        Some(section)
    }
}

#[derive(Debug,Clone,PartialEq)]
pub struct TiltedOneWallSourcePlan {
    pub bottom:TiltedOneWallSection,
    pub top:TiltedOneWallSection,
    pub box_bottom:FaceId,
    pub box_top:FaceId,
    pub box_min_x:FaceId,
    pub cutter_side:FaceId,
    /// Exactly four original face-pair G1 arcs; no generated B-rep edges.
    pub cut_arcs:[GlobalArcRef;4],
    /// Four source-topology-reconciled physical endpoint events.
    pub events:[GlobalEventId;4],
    pub source_cuts:GlobalAnalyticCutPlan,
    pub can_apply:bool,
    pub closed_shell_proven:bool,
}
#[derive(Debug,Clone)]
pub enum TiltedOneWallSourceResult {
    Ready(TiltedOneWallSourcePlan),
    Held(TiltedOneWallHold),
}

impl TiltedOneWallSourcePlan {
    pub fn validate_against(
        &self, first:&BrepBody, second:&BrepBody,
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        self.source_cuts.validate_against(first,second,tolerance)?;
        match Self::build(first,second,BooleanOperation::Difference,tolerance)? {
            TiltedOneWallSourceResult::Ready(ref rebuilt) if self==rebuilt=>Ok(()),
            _=>Err(invalid("G8a source geometry or analytic cut cycle changed")),
        }
    }

    /// Read-only. Strong conservative source restrictions intentionally
    /// exclude opposing-wall, cap-contact, tangent and uncertain cuts.
    pub fn build(
        first:&BrepBody,second:&BrepBody,operation:BooleanOperation,
        tolerance:GeometryTolerance,
    )->Result<TiltedOneWallSourceResult,GeometryError>{
        use TiltedOneWallHold as Hold;
        tolerance.validate()?; first.validate()?; second.validate()?;
        if operation!=BooleanOperation::Difference{
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnsupportedOperation));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cylinder))=
            (first,second) else {
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnsupportedOperand));
        };
        match crate::boolean::certify_box(block) {
            Ok(())=>{},
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(TiltedOneWallSourceResult::Held(Hold::NonCanonicalBox)),
        }
        let bounds=block.bbox;
        let frame=cylinder.frame;
        let axis=frame.normal;
        let tilt=axis.x.hypot(axis.z);
        if axis.y<=0.90 || tilt<=tolerance.angular*64.0 {
            return Ok(TiltedOneWallSourceResult::Held(Hold::NotPositiveYTransverseTilt));
        }
        let extent=bounds.max.sub(bounds.min);
        let scale=extent.length_squared().sqrt().max(cylinder.radius)
            .max(cylinder.height).max(1.0);
        let margin=tolerance.length_at(scale)*256.0;
        if !margin.is_finite()||margin<=0.0 ||
            !cylinder.radius.is_finite() || cylinder.radius<=margin*8.0 {
            return Ok(TiltedOneWallSourceResult::Held(Hold::AmbiguousWallCrossing));
        }
        // The extremal projections of a convex AABB onto a linear axis
        // occur at its 8 corners. Strictly containing ALL box corners in
        // the cutter's axial slab proves that cylinder end-cap planes
        // cannot touch or trim the box, without sampling.
        let mut minimum=f64::INFINITY;
        let mut maximum=f64::NEG_INFINITY;
        for x in [bounds.min.x,bounds.max.x] {
            for y in [bounds.min.y,bounds.max.y] {
                for z in [bounds.min.z,bounds.max.z] {
                    let t=p(x,y,z).sub(frame.origin).dot(axis);
                    minimum=minimum.min(t); maximum=maximum.max(t);
                }
            }
        }
        if !minimum.is_finite()||!maximum.is_finite()||
            minimum<=margin||maximum>=cylinder.height-margin{
            return Ok(TiltedOneWallSourceResult::Held(Hold::CutterCapsReachBox));
        }

        // All ellipse coefficients are independent of Y and their centers
        // are affine in Y. Strict coordinate envelopes at BOTH endpoints
        // therefore prove separation from Xmax and Zmin/Zmax at every Y.
        let make=|y:f64|TiltedOneWallSection::build(
            frame.origin,frame.x,frame.y,axis,cylinder.radius,
            y,bounds,margin,tolerance.angular);
        let (Some(bottom),Some(top))=(make(bounds.min.y),make(bounds.max.y))
            else{return Ok(TiltedOneWallSourceResult::Held(Hold::NotIsolatedMinXWall));};

        let a=first.shared()?; let b=second.shared()?;
        let ylo=find_source_plane(&a,1,bounds.min.y,-1.0,margin);
        let yhi=find_source_plane(&a,1,bounds.max.y,1.0,margin);
        let xlo=find_source_plane(&a,0,bounds.min.x,-1.0,margin);
        let (Some(ylo),Some(yhi),Some(xlo))=(ylo,yhi,xlo)else{
            return Ok(TiltedOneWallSourceResult::Held(Hold::MissingSourceFace));
        };
        let sides=b.faces.iter().filter(|f|
            matches!(b.face_surface(f.id),Some(Surface3::Cylinder(_)))
        ).map(|f|f.id).collect::<Vec<_>>();
        if sides.len()!=1 {
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnsupportedOperand));
        }
        let side=sides[0];
        let registry=match GlobalAnalyticCutPlan::build(first,second,tolerance)? {
            GlobalAnalyticCutResult::Ready(v)=>v,
            GlobalAnalyticCutResult::Held(r)=>
                return Ok(TiltedOneWallSourceResult::Held(Hold::UpstreamCuts(r))),
        };
        if registry.can_apply||registry.closed_shell_proven ||
            registry.arcs.len()!=4 ||registry.physical_event_count!=4{
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnexpectedSourceCuts));
        }
        let mut all_arcs=BTreeSet::new();
        let mut cap_counts=[0usize;3];
        let mut degree=BTreeMap::<GlobalEventId,usize>::new();
        for arc in &registry.arcs{
            if !all_arcs.insert(arc.source) || arc.second_face!=side ||
                arc.virtual_periodic_anchor {
                return Ok(TiltedOneWallSourceResult::Held(Hold::SourceCarrierMismatch));
            }
            let kind=if arc.first_face==ylo {0}
                else if arc.first_face==yhi {1}
                else if arc.first_face==xlo {2}
                else{return Ok(TiltedOneWallSourceResult::Held(Hold::UnexpectedSourceCuts));};
            let harmonic=matches!(arc.carrier,Curve3::Harmonic(_)) &&
                matches!(arc.face_pcurves[0],Curve2::Harmonic(_)) &&
                matches!(arc.face_pcurves[1],Curve2::CylinderSection(_));
            let generator=matches!(arc.carrier,Curve3::Line(_)) &&
                matches!(arc.face_pcurves[0],Curve2::Line(_)) &&
                matches!(arc.face_pcurves[1],Curve2::Line(_));
            if !(harmonic || (kind==2 && generator)) {
                return Ok(TiltedOneWallSourceResult::Held(Hold::SourceCarrierMismatch));
            }
            cap_counts[kind]+=1;
            let Some([start,end])=arc.physical_endpoints else {
                return Ok(TiltedOneWallSourceResult::Held(Hold::UnreconciledPhysicalCycle));
            };
            if start==end {return Ok(TiltedOneWallSourceResult::Held(
                Hold::UnreconciledPhysicalCycle));}
            *degree.entry(start).or_default()+=1;
            *degree.entry(end).or_default()+=1;

            let crate::ParameterRange::Bounded{start:lo,end:hi}=arc.parameter_domain
                else{return Err(invalid("G8a unbounded source cut arc"));};
            if kind<2 && hi-lo<=std::f64::consts::PI+tolerance.angular*32.0 {
                return Ok(TiltedOneWallSourceResult::Held(Hold::SourceCarrierMismatch));
            }
            let pa=a.face_surface(arc.first_face)
                .ok_or_else(||invalid("G8a missing original box surface"))?;
            let pb=b.face_surface(arc.second_face)
                .ok_or_else(||invalid("G8a missing cylinder support"))?;
            for f in [0.0,0.125,0.25,0.5,0.75,0.875,1.0]{
                let t=lo+(hi-lo)*f;
                let xyz=arc.carrier.evaluate(t)?;
                let uv0=arc.face_pcurves[0].evaluate(t)?;
                let uv1=arc.face_pcurves[1].evaluate(t)?;
                let lhs=pa.evaluate(uv0.u,uv0.v)?.point;
                let rhs=pb.evaluate(uv1.u,uv1.v)?.point;
                if distance(lhs,xyz)>margin/8.0 ||
                    distance(rhs,xyz)>margin/8.0 ||
                    xyz.y<bounds.min.y-margin/8.0 ||
                    xyz.y>bounds.max.y+margin/8.0 ||
                    xyz.z<=bounds.min.z+margin/8.0 ||
                    xyz.z>=bounds.max.z-margin/8.0 ||
                    xyz.x<bounds.min.x-margin/8.0 ||
                    xyz.x>=bounds.max.x-margin/8.0 ||
                    (kind==0&&(xyz.y-bounds.min.y).abs()>margin/8.0) ||
                    (kind==1&&(xyz.y-bounds.max.y).abs()>margin/8.0) ||
                    (kind==2&&(xyz.x-bounds.min.x).abs()>margin/8.0){
                    return Ok(TiltedOneWallSourceResult::Held(Hold::SourceCarrierMismatch));
                }
            }
        }
        if cap_counts!=[1,1,2] ||degree.len()!=4 ||
            degree.values().any(|count|*count!=2){
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnreconciledPhysicalCycle));
        }
        let events=registry.global_nodes.events.iter()
            .filter(|e|e.kind==GlobalEventKind::BoundaryEndpoint)
            .collect::<Vec<_>>();
        if events.len()!=4 {
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnreconciledPhysicalCycle));
        }
        let expected=[
            bottom.wall_points[0],bottom.wall_points[1],
            top.wall_points[0],top.wall_points[1],
        ];
        let mut seen=[false;4];
        for e in &events{
            let indices=expected.iter().enumerate()
                .filter(|(_,p)|distance(e.point,**p)<=margin/8.0)
                .map(|(index,_)|index).collect::<Vec<_>>();
            if indices.len()!=1 ||seen[indices[0]] || !degree.contains_key(&e.id){
                return Ok(TiltedOneWallSourceResult::Held(Hold::UnreconciledPhysicalCycle));
            }
            seen[indices[0]]=true;
        }
        // A valence-two graph on four events can be two disjoint double
        // edges. Require a single connected 4-cycle by source event IDs.
        let mut reached=BTreeSet::new();
        let mut frontier=vec![events[0].id];
        while let Some(id)=frontier.pop() {
            if !reached.insert(id){continue;}
            for arc in &registry.arcs {
                let [a,b]=arc.physical_endpoints.expect("checked above");
                if a==id {frontier.push(b);}
                if b==id {frontier.push(a);}
            }
        }
        if reached.len()!=4{
            return Ok(TiltedOneWallSourceResult::Held(Hold::UnreconciledPhysicalCycle));
        }
        let refs=registry.arcs.iter().map(|a|a.source)
            .collect::<Vec<_>>().try_into()
            .map_err(|_|invalid("G8a missing four source arcs"))?;
        let physical=events.iter().map(|e|e.id)
            .collect::<Vec<_>>().try_into()
            .map_err(|_|invalid("G8a missing four source events"))?;
        Ok(TiltedOneWallSourceResult::Ready(Self{
            bottom,top,box_bottom:ylo,box_top:yhi,box_min_x:xlo,
            cutter_side:side,cut_arcs:refs,events:physical,
            source_cuts:registry,can_apply:false,closed_shell_proven:false,
        }))
    }
}
