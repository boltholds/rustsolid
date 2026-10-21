//! G8b independent SOURCE-ONLY bounded volume for a tilted one-wall notch.
//! No Boolean candidate, BrepModel, mesh, output mass, or G8a/G8c plan is consulted.
//! The derived lower/upper bounds include engineering f64 guards, NOT directed-rounded proof.
use rustsolid_foundation::{GeometryTolerance, Point3, GeometryError};
use rustsolid_geometry::{Frame3, ParameterRange};
use rustsolid_source::{BooleanKind, BooleanIntent, SourceGeometry, SourceRef, SourceSnapshot, SourceStamp};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeInterval {
    pub lower:f64, pub upper:f64,
    pub integration_panels:usize, pub numerical_guard:f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct SourceVolumeEvidence {
    retained:VolumeInterval,
    removed_inside_box:[f64;2],
    original_box_volume:f64,
    analytic_removed_cross_check:f64,
    numerical_guard:f64,
    angular_panels:usize,
    height_panels:usize,
    source_stamp:SourceStamp,
    source_refs:[SourceRef;2],
    intent_stamp:SourceStamp,
}
impl SourceVolumeEvidence {
    pub fn interval(&self)->VolumeInterval {self.retained}
    pub fn removed_interval(&self)->[f64;2] {self.removed_inside_box}
    pub fn source_stamp(&self)->&SourceStamp{&self.source_stamp}
    pub fn intent_stamp(&self)->&SourceStamp{&self.intent_stamp}
    pub fn source_refs(&self)->[SourceRef;2]{self.source_refs}
    pub fn original_box_volume(&self)->f64{self.original_box_volume}
    pub fn angular_panels(&self)->usize{self.angular_panels}
    pub fn height_panels(&self)->usize{self.height_panels}
    pub fn numerical_guard(&self)->f64{self.numerical_guard}
    pub fn analytic_removed_cross_check(&self)->f64{self.analytic_removed_cross_check}
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TiltedVolumeHold {
    UnsupportedOperand,
    InconsistentUnits,
    InvalidOriginalBoxMass,
    InvalidSourceSection,
    NumericalBudgetExceeded,
    UnresolvedPolygonEnclosure,
    UnresolvedHeightQuadrature,
    SourceConservationFailed,
    IndependentAreaCrossCheckMismatch,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TiltedSectionAreaBracket {
    pub height_y:f64,
    pub lower_area:f64,
    pub upper_area:f64,
    pub absolute_jacobian:f64,
    pub normalized_wall_offset:f64,
    pub major_arc_span:f64,
    pub angular_panels:usize,
}
#[derive(Clone,Copy)]
struct SourceEllipse {
    origin_offset:f64, // original cylinder origin.x - original MIN-X
    origin_y:f64,
    x_slope:f64,
    x_amplitude:f64,
    jac:f64, // absolute determinant of the horizontal radial projection
    ymin:f64,ymax:f64,
}
impl SourceEllipse {
    fn from_sources(min:Point3,max:Point3,frame:Frame3,radius:f64)->Option<Self>{
        let a=frame.normal;
        if a.y<=0.90||!a.y.is_finite(){return None;}
        let inv=1./a.y;
        let e1=frame.x.sub(a.scale(frame.x.y*inv)).scale(radius);
        let e2=frame.y.sub(a.scale(frame.y.y*inv)).scale(radius);
        let source=Self{
            origin_offset:frame.origin.x-min.x,
            origin_y:frame.origin.y,
            x_slope:a.x*inv,
            x_amplitude:e1.x.hypot(e2.x),
            jac:(e1.x*e2.z-e2.x*e1.z).abs(),
            ymin:min.y,ymax:max.y,
        };
        let d0=source.d(source.ymin);
        let d1=source.d(source.ymax);
        let expected=radius*radius*inv;
        if ![d0,d1,expected,source.x_amplitude,source.jac,source.x_slope]
            .iter().all(|x|x.is_finite())
            ||source.x_amplitude<=0.||source.jac<=0.
            ||!(0.0..1.0).contains(&d0)||!(0.0..1.0).contains(&d1)
            ||(source.jac-expected).abs()>512.*f64::EPSILON*expected.max(1.){
            return None;
        }
        Some(source)
    }
    fn d(self,y:f64)->f64{
        (self.origin_offset+self.x_slope*(y-self.origin_y))/self.x_amplitude
    }
    fn height(self)->f64{self.ymax-self.ymin}
    fn span(self,y:f64)->Option<f64>{
        let d=self.d(y);
        if !(0.0..1.0).contains(&d){return None;}
        let span=2.*(-d).acos();
        if span>PI&&span<2.*PI&&span.is_finite(){Some(span)}else{None}
    }
    fn fourth_bound(self)->Option<f64>{
        // Area(d)=J*(PI/2+asin(d)+d sqrt(1-d^2))
        // A''''_y=-6J*d*(d'_y)^4/(1-d^2)^(5/2).
        // d>0 and affine in y, so use its maximum endpoint value.
        let dmax=self.d(self.ymin).max(self.d(self.ymax));
        let denominator=((1.-dmax)*(1.+dmax)).powf(2.5);
        let slope=self.x_slope/self.x_amplitude;
        let result=6.*self.jac*dmax*slope.powi(4)/denominator;
        if result.is_finite()&&result>=0.{Some(result)}else{None}
    }
}
fn polygon(
    source:SourceEllipse,y:f64,n:usize,
)->Result<TiltedSectionAreaBracket,TiltedVolumeHold>{
    if n<8||!n.is_power_of_two(){
        return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
    }
    let span=source.span(y).ok_or(TiltedVolumeHold::InvalidSourceSection)?;
    let step=span/(n as f64);
    let half=step/2.;
    if step<=0.||step>=PI/2.||half.cos()<=0.{
        return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
    }
    // Shoelace of N source ellipse chords and their true wall chord:
    // inner = |J|/2 * [N sin(step)-sin(span)].
    // Tangent-intersection polygon adds the strictly positive
    //   |J| * N * sin(step/2)^3 / cos(step/2).
    // No circular-segment antiderivative or sampled mesh is consulted.
    let inner=source.jac*0.5*((n as f64)*step.sin()-span.sin());
    let gap=source.jac*(n as f64)*half.sin().powi(3)/half.cos();
    let outer=inner+gap;
    if !inner.is_finite()||!outer.is_finite()||inner<=0.||
        outer<=inner||outer>=source.jac*PI {
        return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
    }
    Ok(TiltedSectionAreaBracket{
        height_y:y,lower_area:inner,upper_area:outer,
        absolute_jacobian:source.jac,normalized_wall_offset:source.d(y),
        major_arc_span:span,angular_panels:n,
    })
}
#[derive(Default)]
struct Kahan{sum:f64,carry:f64}
impl Kahan{
    fn add(&mut self,v:f64){
        let x=v-self.carry;let t=self.sum+x;
        self.carry=(t-self.sum)-x;self.sum=t;
    }
}
struct SweptBracket {
    retained:VolumeInterval,
    removed:[f64;2],
    sections:[TiltedSectionAreaBracket;3],
    n:usize,m:usize,error:f64,guard:f64,
}
fn enclose(
    source:SourceEllipse,box_volume:f64,target:f64,
)->Result<SweptBracket,TiltedVolumeHold>{
    let height=source.height();
    if height<=0.||target<=0.||!target.is_finite(){
        return Err(TiltedVolumeHold::InvalidSourceSection);
    }
    let fourth=source.fourth_bound()
        .ok_or(TiltedVolumeHold::InvalidSourceSection)?;
    let mut m=8usize;
    let remainder=loop{
        let h=height/(m as f64);
        let e=height*h.powi(4)*fourth/180.; // composite Simpson remainder
        if !e.is_finite(){return Err(TiltedVolumeHold::NumericalBudgetExceeded);}
        if e<=target/16.{break e;}
        if m>=1024 {return Err(TiltedVolumeHold::UnresolvedHeightQuadrature);}
        m*=2;
    };
    let max_span=source.span(source.ymin)
        .ok_or(TiltedVolumeHold::InvalidSourceSection)?
        .max(source.span(source.ymax)
            .ok_or(TiltedVolumeHold::InvalidSourceSection)?);
    let mut n=1024usize;
    let gap=loop{
        let u=max_span/(2.*n as f64);
        let g=source.jac*(n as f64)*u.sin().powi(3)/u.cos();
        if !g.is_finite()||g<=0. {
            return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
        }
        if height*g<=target/2.{break g;}
        if n>=262144 {
            return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
        }
        n*=2;
    };
    // Explicit practical f64 guard. These are conservative engineering
    // margins and NOT directed-rounding certified floating point intervals.
    let guard=4096.*f64::EPSILON*
        (box_volume+2.*source.jac*PI*height)
        +64.*f64::EPSILON*(n as f64)*source.jac*height;
    if !guard.is_finite()||guard<0.||guard>=target/8. {
        return Err(TiltedVolumeHold::NumericalBudgetExceeded);
    }
    if height*gap+2.*(remainder+guard)>target{
        return Err(TiltedVolumeHold::UnresolvedPolygonEnclosure);
    }
    let h=height/m as f64;
    let mut lower=Kahan::default();
    let mut upper=Kahan::default();
    for i in 0..=m {
        let y=if i==m{source.ymax}else{source.ymin+i as f64*h};
        let area=polygon(source,y,n)?;
        let weight=if i==0||i==m{1.}else if i%2==0{2.}else{4.};
        lower.add(weight*area.lower_area);
        upper.add(weight*area.upper_area);
    }
    let removed_lo=h*lower.sum/3.-remainder-guard;
    let removed_hi=h*upper.sum/3.+remainder+guard;
    let retained_lo=box_volume-removed_hi;
    let retained_hi=box_volume-removed_lo;
    if ![removed_lo,removed_hi,retained_lo,retained_hi].iter().all(|x|x.is_finite())
        ||removed_lo<=0.||removed_hi<=removed_lo||retained_lo<=0.
        ||retained_hi<=retained_lo||retained_hi>=box_volume
        ||removed_hi-removed_lo>target {
        return Err(TiltedVolumeHold::SourceConservationFailed);
    }
    Ok(SweptBracket{
        retained:VolumeInterval{
            lower:retained_lo,upper:retained_hi,
            integration_panels:m,numerical_guard:guard,
        },
        removed:[removed_lo,removed_hi],
        sections:[
            polygon(source,source.ymin,n)?,
            polygon(source,0.5*source.ymin+0.5*source.ymax,n)?,
            polygon(source,source.ymax,n)?,
        ],
        n,m,error:remainder,guard,
    })
}
fn analytic_cross_check(source:SourceEllipse)->Result<(f64,f64),TiltedVolumeHold>{
    // Independent analytic area antiderivative for the *source* only.
    // P'(d) = pi/2 + asin(d) + d sqrt(1-d^2).
    let d0=source.d(source.ymin);
    let d1=source.d(source.ymax);
    let h=source.height();
    let jac=source.jac;
    let f=|d:f64|PI/2.+d.asin()+d*((1.-d)*(1.+d)).sqrt();
    let p=|d:f64|{
        let root=((1.-d)*(1.+d)).sqrt();
        PI*d/2.+d*d.asin()+root-root.powi(3)/3.
    };
    let (integral,error)=if (d1-d0).abs()>1e-4{
        (jac*h*(p(d1)-p(d0))/(d1-d0),0.)
    }else{
        // Avoid subtractive cancellation for nearly constant sections.
        const M:usize=32;
        let dy=h/(M as f64);
        let mut sum=Kahan::default();
        for i in 0..=M{
            let y=if i==M{source.ymax}else{source.ymin+i as f64*dy};
            let weight=if i==0||i==M{1.}else if i%2==0{2.}else{4.};
            sum.add(weight*f(source.d(y)));
        }
        (jac*sum.sum*dy/3.,
            h*dy.powi(4)*source.fourth_bound()
                .ok_or(TiltedVolumeHold::InvalidSourceSection)?/180.)
    };
    if !integral.is_finite()||integral<=0.||!error.is_finite(){
        return Err(TiltedVolumeHold::IndependentAreaCrossCheckMismatch);
    }
    Ok((integral,error))
}

/// Source-only G8b witness. G8a graph is not required; this function independently
/// rejects any cylinder clipped by axial end caps or crossing Xmax/Z walls.
pub fn evaluate_tilted_one_wall(
    snapshot:&SourceSnapshot,intent:&BooleanIntent,
)->Result<SourceVolumeEvidence,TiltedVolumeHold>{
    if snapshot.check_intent(intent).is_err() || intent.operation()!=BooleanKind::Difference {
        return Err(TiltedVolumeHold::UnsupportedOperand);
    }
    let (left,right)=intent.operands();
    let block=snapshot.resolve(left).map_err(|_|TiltedVolumeHold::UnsupportedOperand)?;
    let cyl=snapshot.resolve(right).map_err(|_|TiltedVolumeHold::UnsupportedOperand)?;
    if block.units()!=cyl.units() || !matches!(intent.unit_policy(),rustsolid_source::UnitPolicy::RequireIdentical) {
        return Err(TiltedVolumeHold::InconsistentUnits);
    }
    let SourceGeometry::AxisAlignedBlock{min,max}=block.geometry() else {
        return Err(TiltedVolumeHold::UnsupportedOperand);
    };
    let SourceGeometry::FiniteCylinder{frame,radius,axial}=cyl.geometry() else {
        return Err(TiltedVolumeHold::UnsupportedOperand);
    };
    let ParameterRange::Bounded{start,end}=axial else {
        return Err(TiltedVolumeHold::InvalidSourceSection);
    };
    let delta=max.sub(*min);
    let box_volume=delta.x*delta.y*delta.z;
    let height=end-start;
    let scale=delta.x.hypot(delta.y).hypot(delta.z).max(*radius).max(height).max(1.);
    let tol=intent.tolerance();
    let budget=8.*tol.length_at(scale)*scale*scale;
    if !box_volume.is_finite() || box_volume<=0. || !budget.is_finite()||budget<=0. {
        return Err(TiltedVolumeHold::InvalidOriginalBoxMass);
    }
    if frame.normal.y<=0.90 {return Err(TiltedVolumeHold::InvalidSourceSection);}
    // Projection along finite cylinder axis, for all angles and all box heights.
    // This enforces absence of cutter end caps in the body being cut.
    let axis=frame.normal;
    let v_amp=radius*frame.x.y.hypot(frame.y.y)/axis.y;
    let v_min=(min.y-frame.origin.y)/axis.y-v_amp;
    let v_max=(max.y-frame.origin.y)/axis.y+v_amp;
    let margin=tol.length_at(scale)*64.;
    if !v_min.is_finite() || !v_max.is_finite() ||v_min<=start+margin||v_max>=end-margin{
        return Err(TiltedVolumeHold::InvalidSourceSection);
    }
    // G8 source is limited to the one-wall case, and must remain clear of
    // Xmax, Zmin and Zmax for every horizontal cutter section.
    let e1=frame.x.sub(axis.scale(frame.x.y/axis.y)).scale(*radius);
    let e2=frame.y.sub(axis.scale(frame.y.y/axis.y)).scale(*radius);
    let x_amp=e1.x.hypot(e2.x);
    let z_amp=e1.z.hypot(e2.z);
    for y in [min.y,max.y]{
        let cx=frame.origin.x+(y-frame.origin.y)*axis.x/axis.y;
        let cz=frame.origin.z+(y-frame.origin.y)*axis.z/axis.y;
        if cx+x_amp>=max.x-margin || cz-z_amp<=min.z+margin||cz+z_amp>=max.z-margin{
            return Err(TiltedVolumeHold::InvalidSourceSection);
        }
    }
    let shape=SourceEllipse::from_sources(*min,*max,*frame,*radius)
        .ok_or(TiltedVolumeHold::InvalidSourceSection)?;
    let bound=enclose(shape,box_volume,budget/8.)?;
    let (analytic,analytic_error)=analytic_cross_check(shape)?;
    if analytic<bound.removed[0]-analytic_error-bound.guard||
       analytic>bound.removed[1]+analytic_error+bound.guard {
        return Err(TiltedVolumeHold::IndependentAreaCrossCheckMismatch);
    }
    // Only the finite-height box intersection contributes material removal.
    // The old G8 implementation incorrectly compared against the entire
    // cutter cylinder's full volume as a source proxy.
    if (bound.retained.lower+bound.removed[1]-box_volume).abs()>budget/4.||
       (bound.retained.upper+bound.removed[0]-box_volume).abs()>budget/4. {
        return Err(TiltedVolumeHold::SourceConservationFailed);
    }
    Ok(SourceVolumeEvidence{
        retained:bound.retained,removed_inside_box:bound.removed,
        original_box_volume:box_volume,analytic_removed_cross_check:analytic,
        numerical_guard:bound.guard,angular_panels:bound.n,height_panels:bound.m,
        source_stamp:snapshot.stamp().clone(),source_refs:[left,right],
        intent_stamp:intent.intent_stamp().clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsolid_source::{SourceId,SourceRevision,SourceSpec,UnitPolicy};
    use rustsolid_geometry::{ModelUnits,Frame3,ParameterRange};
    fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
    fn new_source(i:u128,geom:SourceGeometry)->SourceSpec {
        SourceSpec::new(SourceId::new(i).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),geom,GeometryTolerance::default()).unwrap()
    }
    fn sample(angle:f64)->(SourceSnapshot,BooleanIntent) {
        let frame=Frame3::from_axes(p(0.,-2.,0.),
            p(angle.cos(),-angle.sin(),0.),p(0.,0.,-1.),GeometryTolerance::default()).unwrap();
        let block=new_source(1,SourceGeometry::AxisAlignedBlock{
            min:p(-1.,1.,-3.),max:p(3.,4.,3.)
        });
        let cyl=new_source(2,SourceGeometry::FiniteCylinder{
            frame,radius:1.5,axial:ParameterRange::bounded(0.,10.).unwrap()
        });
        let snap=SourceSnapshot::new([block.clone(),cyl.clone()]).unwrap();
        let intent=snap.intent(BooleanKind::Difference,block.reference(),cyl.reference(),
            GeometryTolerance::default(),UnitPolicy::RequireIdentical).unwrap();
        (snap,intent)
    }
    #[test]
    fn tilted_source_interval_is_finite_and_nontrivial(){
        for tilt in [0.04,-0.03,0.055]{
            let (snap,intent)=sample(tilt);
            let result=evaluate_tilted_one_wall(&snap,&intent).unwrap();
            let v=result.interval();
            assert!(v.lower>0.&&v.lower<v.upper&&v.upper<72.);
            assert!(v.upper-v.lower<1e-4,
                "expected tightly bounded G8b source enclosure: {v:?}");
            assert!(result.analytic_removed_cross_check()>0.);
        }
    }
    #[test]
    fn rejects_stale_source_or_improper_boolean(){
        let (snap,intent)=sample(0.04);
        let stale=snap.with_replacement(
            snap.get(SourceId::new(1).unwrap()).unwrap().successor(
                SourceGeometry::AxisAlignedBlock{min:p(-1.,1.,-3.),max:p(4.,4.,3.)},
                GeometryTolerance::default()).unwrap()).unwrap();
        assert!(evaluate_tilted_one_wall(&stale,&intent).is_err());
        let union=snap.intent(BooleanKind::Union,intent.operands().0,intent.operands().1,
            GeometryTolerance::default(),UnitPolicy::RequireIdentical).unwrap();
        assert!(evaluate_tilted_one_wall(&snap,&union).is_err());
    }
    #[test]
    fn rejects_cylinder_caps_inside_box(){
        let (snap,intent)=sample(0.04);
        let right=SourceId::new(2).unwrap();
        let spec=snap.get(right).unwrap();
        let SourceGeometry::FiniteCylinder{frame,radius,..}=spec.geometry() else{unreachable!()};
        let invalid=spec.successor(SourceGeometry::FiniteCylinder{
            frame:*frame,radius:*radius,axial:ParameterRange::bounded(2.,4.).unwrap()
        },GeometryTolerance::default()).unwrap();
        let next=snap.with_replacement(invalid).unwrap();
        let changed=next.intent(BooleanKind::Difference,
            next.get(SourceId::new(1).unwrap()).unwrap().reference(),
            next.get(right).unwrap().reference(),GeometryTolerance::default(),
            UnitPolicy::RequireIdentical).unwrap();
        assert!(evaluate_tilted_one_wall(&next,&changed).is_err());
    }
}
