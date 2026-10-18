//! G8b: source-only independent volume enclosure for the certified G8a
//! tilted one-wall notch. No candidate B-rep, display mesh or flux is read.
//!
//! In a horizontal plane the cutter forms an affine ellipse, clipped by the
//! original MIN-X wall. Chord and tangent polygon shoelace areas bound the
//! removed section. Composite Simpson with a bounded fourth derivative then
//! encloses the height integral. All bounds are exact in real arithmetic;
//! f64 includes an explicit practical guard, not formal directed rounding.

use std::f64::consts::PI;
use crate::{
    BrepBody, CylindricalBrep, GeometryError, GeometryTolerance,
    SourceVolumeBracket, TiltedOneWallSourcePlan,
};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}

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
#[derive(Debug,Clone,PartialEq)]
pub struct IndependentTiltedVolumeOracle {
    /// Source-only retained material, not a verified result shell volume.
    pub retained:SourceVolumeBracket,
    pub removed_inside_box:[f64;2],
    pub original_box_volume:f64,
    /// A second check from the elliptic segment antiderivative, not used
    /// to construct either polygon/Simpson interval endpoint.
    pub analytic_removed_cross_check:f64,
    pub analytic_cross_check_error:f64,
    /// Ymin, midpoint, Ymax: all source-derived.
    pub section_witnesses:[TiltedSectionAreaBracket;3],
    pub angular_panels:usize,
    pub height_panels:usize,
    pub simpson_remainder_bound:f64,
    pub numerical_guard:f64,
    pub source_revisions:[Option<u64>;2],
    pub tolerance:GeometryTolerance,
    pub source_bracket_established:bool,
    pub result_shell_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum IndependentTiltedVolumeResult {
    Ready(IndependentTiltedVolumeOracle),
    Held(TiltedVolumeHold),
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
    fn from_sources(bounds:crate::BoundingBox,cyl:&CylindricalBrep)->Option<Self>{
        let a=cyl.frame.normal;
        if a.y<=0.90||!a.y.is_finite(){return None;}
        let inv=1./a.y;
        let e1=cyl.frame.x.sub(a.scale(cyl.frame.x.y*inv)).scale(cyl.radius);
        let e2=cyl.frame.y.sub(a.scale(cyl.frame.y.y*inv)).scale(cyl.radius);
        let source=Self{
            origin_offset:cyl.frame.origin.x-bounds.min.x,
            origin_y:cyl.frame.origin.y,
            x_slope:a.x*inv,
            x_amplitude:e1.x.hypot(e2.x),
            jac:(e1.x*e2.z-e2.x*e1.z).abs(),
            ymin:bounds.min.y,ymax:bounds.max.y,
        };
        let d0=source.d(source.ymin);
        let d1=source.d(source.ymax);
        let expected=cyl.radius*cyl.radius*inv;
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
    retained:SourceVolumeBracket,
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
        retained:SourceVolumeBracket{
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
impl IndependentTiltedVolumeOracle {
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        source:&TiltedOneWallSourcePlan,tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        if self.tolerance!=tol ||Self::build(first,second,source,tol)?
            !=IndependentTiltedVolumeResult::Ready(self.clone()){
            return Err(invalid("G8b source-derived volume witness changed"));
        }
        Ok(())
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        source:&TiltedOneWallSourcePlan,tol:GeometryTolerance,
    )->Result<IndependentTiltedVolumeResult,GeometryError>{
        tol.validate()?;
        source.validate_against(first,second,tol)?;
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cyl))=(first,second)
            else{return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::UnsupportedOperand));};
        if first.shared()?.units!=second.shared()?.units{
            return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::InconsistentUnits));
        }
        let extent=block.bbox.max.sub(block.bbox.min);
        let box_volume=extent.x*extent.y*extent.z;
        let scale=extent.length_squared().sqrt()
            .max(cyl.radius).max(cyl.height).max(1.);
        let budget=8.*tol.length_at(scale)*scale*scale;
        if !box_volume.is_finite()||box_volume<=0.||
            !budget.is_finite()||budget<=0.||
            (block.mass.volume-box_volume).abs()>budget {
            return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::InvalidOriginalBoxMass));
        }
        let shape=match SourceEllipse::from_sources(block.bbox,cyl){
            Some(s)=>s,None=>return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::InvalidSourceSection)),
        };
        let bound=match enclose(shape,box_volume,budget/8.){
            Ok(v)=>v,Err(e)=>return Ok(IndependentTiltedVolumeResult::Held(e)),
        };
        let (analytic,analytic_error)=match analytic_cross_check(shape){
            Ok(v)=>v,Err(e)=>return Ok(IndependentTiltedVolumeResult::Held(e)),
        };
        if analytic<bound.removed[0]-analytic_error-bound.guard ||
            analytic>bound.removed[1]+analytic_error+bound.guard {
            return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::IndependentAreaCrossCheckMismatch));
        }
        if bound.removed[1]>cyl.mass.volume+budget ||
            (bound.retained.lower+bound.removed[1]-box_volume).abs()>budget/4. ||
            (bound.retained.upper+bound.removed[0]-box_volume).abs()>budget/4. {
            return Ok(IndependentTiltedVolumeResult::Held(
                TiltedVolumeHold::SourceConservationFailed));
        }
        Ok(IndependentTiltedVolumeResult::Ready(Self{
            retained:bound.retained,removed_inside_box:bound.removed,
            original_box_volume:box_volume,
            analytic_removed_cross_check:analytic,
            analytic_cross_check_error:analytic_error,
            section_witnesses:bound.sections,
            angular_panels:bound.n,height_panels:bound.m,
            simpson_remainder_bound:bound.error,numerical_guard:bound.guard,
            source_revisions:[Some(block.revision),None],tolerance:tol,
            source_bracket_established:true,
            result_shell_certified:false,can_apply:false,
        }))
    }
}
