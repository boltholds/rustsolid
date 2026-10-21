//! G8b independent source-only volume enclosure for an axis-compatible tilted
//! circular cylinder cutting exactly MIN-X of a finite block.
//! No candidate B-rep, source graph, cached mass, mesh or G8c result is read.
//! f64 with explicit engineering guard: NOT formal directed-rounding intervals.
use std::f64::consts::PI;
use rustsolid_foundation::{GeometryTolerance, Point3};
use rustsolid_geometry::{Frame3, ParameterRange};
use rustsolid_source::{BooleanIntent, BooleanKind, SourceError, SourceGeometry,
    SourceRef, SourceSnapshot, SourceStamp, UnitPolicy};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G8SourceHold {
    UnsupportedOperation,
    UnsupportedSources,
    MixedUnits,
    CylinderCapIntersectsBox,
    UnsupportedWallFamily,
    InvalidSourceSection,
    NumericalBudgetExceeded,
    UnresolvedPolygonEnclosure,
    UnresolvedHeightQuadrature,
    SourceConservationFailed,
    IndependentAreaCrossCheckMismatch,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuardedSourceVolume {
    pub lower:f64,
    pub upper:f64,
    pub integration_panels:usize,
    pub numerical_guard:f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TiltedSectionAreaBracket {
    pub height_y:f64,
    pub lower_area:f64,
    pub upper_area:f64,
    pub absolute_jacobian:f64,
    pub normalized_wall_offset:f64,
    pub major_arc_span:f64,
    pub angular_panels:usize,
}
#[derive(Debug, Clone, PartialEq)]
pub struct G8SourceVolumeEvidence {
    pub retained:GuardedSourceVolume,
    pub removed_inside_box:[f64;2],
    pub original_box_volume:f64,
    pub independent_cross_check:f64,
    pub cross_check_error:f64,
    pub sections:[TiltedSectionAreaBracket;3],
    pub angular_panels:usize,
    pub height_panels:usize,
    pub simpson_remainder:f64,
    pub numerical_guard:f64,
    source_stamp:SourceStamp,
    operands:[SourceRef;2],
}
impl G8SourceVolumeEvidence {
    pub fn source_stamp(&self)->&SourceStamp{&self.source_stamp}
    pub fn operands(&self)->[SourceRef;2]{self.operands}
}
#[derive(Debug, Clone, PartialEq)]
pub enum G8SourceVolumeOutcome {
    Ready(G8SourceVolumeEvidence),
    Held(G8SourceHold),
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
)->Result<TiltedSectionAreaBracket,G8SourceHold>{
    if n<8||!n.is_power_of_two(){
        return Err(G8SourceHold::UnresolvedPolygonEnclosure);
    }
    let span=source.span(y).ok_or(G8SourceHold::InvalidSourceSection)?;
    let step=span/(n as f64);
    let half=step/2.;
    if step<=0.||step>=PI/2.||half.cos()<=0.{
        return Err(G8SourceHold::UnresolvedPolygonEnclosure);
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
        return Err(G8SourceHold::UnresolvedPolygonEnclosure);
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
    retained:GuardedSourceVolume,
    removed:[f64;2],
    sections:[TiltedSectionAreaBracket;3],
    n:usize,m:usize,error:f64,guard:f64,
}
fn enclose(
    source:SourceEllipse,box_volume:f64,target:f64,
)->Result<SweptBracket,G8SourceHold>{
    let height=source.height();
    if height<=0.||target<=0.||!target.is_finite(){
        return Err(G8SourceHold::InvalidSourceSection);
    }
    let fourth=source.fourth_bound()
        .ok_or(G8SourceHold::InvalidSourceSection)?;
    let mut m=8usize;
    let remainder=loop{
        let h=height/(m as f64);
        let e=height*h.powi(4)*fourth/180.; // composite Simpson remainder
        if !e.is_finite(){return Err(G8SourceHold::NumericalBudgetExceeded);}
        if e<=target/16.{break e;}
        if m>=1024 {return Err(G8SourceHold::UnresolvedHeightQuadrature);}
        m*=2;
    };
    let max_span=source.span(source.ymin)
        .ok_or(G8SourceHold::InvalidSourceSection)?
        .max(source.span(source.ymax)
            .ok_or(G8SourceHold::InvalidSourceSection)?);
    let mut n=1024usize;
    let gap=loop{
        let u=max_span/(2.*n as f64);
        let g=source.jac*(n as f64)*u.sin().powi(3)/u.cos();
        if !g.is_finite()||g<=0. {
            return Err(G8SourceHold::UnresolvedPolygonEnclosure);
        }
        if height*g<=target/2.{break g;}
        if n>=262144 {
            return Err(G8SourceHold::UnresolvedPolygonEnclosure);
        }
        n*=2;
    };
    // Explicit practical f64 guard. These are conservative engineering
    // margins and NOT directed-rounding certified floating point intervals.
    let guard=4096.*f64::EPSILON*
        (box_volume+2.*source.jac*PI*height)
        +64.*f64::EPSILON*(n as f64)*source.jac*height;
    if !guard.is_finite()||guard<0.||guard>=target/8. {
        return Err(G8SourceHold::NumericalBudgetExceeded);
    }
    if height*gap+2.*(remainder+guard)>target{
        return Err(G8SourceHold::UnresolvedPolygonEnclosure);
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
        return Err(G8SourceHold::SourceConservationFailed);
    }
    Ok(SweptBracket{
        retained:GuardedSourceVolume{
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
fn analytic_cross_check(source:SourceEllipse)->Result<(f64,f64),G8SourceHold>{
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
                .ok_or(G8SourceHold::InvalidSourceSection)?/180.)
    };
    if !integral.is_finite()||integral<=0.||!error.is_finite(){
        return Err(G8SourceHold::IndependentAreaCrossCheckMismatch);
    }
    Ok((integral,error))
}

/// Computes only from original immutable source parameters. No body constructors
/// are invoked. Unsupported source, unit and cap families explicitly fail closed.
pub fn g8_source_volume(
    snapshot:&SourceSnapshot, intent:&BooleanIntent,
)->Result<G8SourceVolumeOutcome,SourceError>{
    snapshot.check_intent(intent)?;
    if intent.operation()!=BooleanKind::Difference {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedOperation));
    }
    let (lhs,rhs)=intent.operands();
    let left=snapshot.resolve(lhs)?;
    let right=snapshot.resolve(rhs)?;
    if left.units()!=right.units() || !matches!(intent.unit_policy(),UnitPolicy::RequireIdentical) {
        // Explicit conversions require a normalized physical source snapshot;
        // never compute with mixed-units raw coordinates.
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::MixedUnits));
    }
    let (
        SourceGeometry::AxisAlignedBlock{min,max},
        SourceGeometry::FiniteCylinder{frame,radius,axial},
    )=(left.geometry(),right.geometry()) else {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedSources));
    };
    let (ParameterRange::Bounded{start,end})=*axial else {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedSources));
    };
    let min=*min;let max=*max;let frame=*frame;let radius=*radius;
    let tol=intent.tolerance();
    let extent=max.sub(min);
    let box_volume=extent.x*extent.y*extent.z;
    let scale=extent.length_squared().sqrt()
        .max(radius).max(end-start).max(1.0);
    let budget=8.0*tol.length_at(scale)*scale*scale;
    if !box_volume.is_finite() || box_volume<=0.0 ||
        !budget.is_finite() || budget<=0.0 {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::NumericalBudgetExceeded));
    }
    let eps=tol.length_at(scale)*128.0;
    // This bounded family excludes all intrusion of cylinder caps into the box;
    // the eight corners are sufficient for an affine axial coordinate.
    for x in [min.x,max.x] {
        for y in [min.y,max.y] {
            for z in [min.z,max.z] {
                let axial_coordinate=Point3{x,y,z}.sub(frame.origin).dot(frame.normal);
                if axial_coordinate<=start+eps || axial_coordinate>=end-eps {
                    return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::CylinderCapIntersectsBox));
                }
            }
        }
    }
    if frame.normal.y<=0.90 {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedWallFamily));
    }
    let inv=1.0/frame.normal.y;
    let e1=frame.x.sub(frame.normal.scale(frame.x.y*inv)).scale(radius);
    let e2=frame.y.sub(frame.normal.scale(frame.y.y*inv)).scale(radius);
    let rx=e1.x.hypot(e2.x);
    let rz=e1.z.hypot(e2.z);
    if !rx.is_finite() || !rz.is_finite() || rx<=eps || rz<=eps {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::InvalidSourceSection));
    }
    // At every horizontal section the ellipse crosses MIN-X only, with
    // positive clearance from XMAX, ZMIN and ZMAX (affine centres, fixed support).
    for y in [min.y,max.y] {
        let delta=(y-frame.origin.y)*inv;
        let cx=frame.origin.x+frame.normal.x*delta;
        let cz=frame.origin.z+frame.normal.z*delta;
        if cx-rx>=min.x-eps || cx+rx<=min.x+eps ||
           cx+rx>=max.x-eps || cz-rz<=min.z+eps ||
           cz+rz>=max.z-eps {
            return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedWallFamily));
        }
    }
    let Some(shape)=SourceEllipse::from_sources(min,max,frame,radius) else {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::InvalidSourceSection));
    };
    let bound=match enclose(shape,box_volume,budget/8.0) {
        Ok(v)=>v, Err(e)=>return Ok(G8SourceVolumeOutcome::Held(e)),
    };
    let (analytic,analytic_error)=match analytic_cross_check(shape){
        Ok(v)=>v, Err(e)=>return Ok(G8SourceVolumeOutcome::Held(e)),
    };
    if analytic<bound.removed[0]-analytic_error-bound.guard ||
       analytic>bound.removed[1]+analytic_error+bound.guard {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::IndependentAreaCrossCheckMismatch));
    }
    let cylinder_volume=PI*radius*radius*(end-start);
    if !cylinder_volume.is_finite() || bound.removed[1]>cylinder_volume+budget ||
       (bound.retained.lower+bound.removed[1]-box_volume).abs()>budget/4.0 ||
       (bound.retained.upper+bound.removed[0]-box_volume).abs()>budget/4.0 {
        return Ok(G8SourceVolumeOutcome::Held(G8SourceHold::SourceConservationFailed));
    }
    Ok(G8SourceVolumeOutcome::Ready(G8SourceVolumeEvidence {
        retained:bound.retained,
        removed_inside_box:bound.removed,
        original_box_volume:box_volume,
        independent_cross_check:analytic,
        cross_check_error:analytic_error,
        sections:bound.sections,
        angular_panels:bound.n,
        height_panels:bound.m,
        simpson_remainder:bound.error,
        numerical_guard:bound.guard,
        source_stamp:snapshot.stamp().clone(),
        operands:[lhs,rhs],
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsolid_geometry::ModelUnits;
    use rustsolid_source::{SourceId,SourceRevision,SourceSpec};
    fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
    fn source(i:u128,geom:SourceGeometry)->SourceSpec {
        SourceSpec::new(SourceId::new(i).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),geom,GeometryTolerance::default()).unwrap()
    }
    fn snapshot(theta:f64)->(SourceSnapshot,BooleanIntent){
        let a=source(1,SourceGeometry::AxisAlignedBlock {
            min:p(-1.,1.,-3.),max:p(3.,4.,3.),
        });
        let frame=Frame3::from_axes(p(0.,-2.,0.),
            p(theta.cos(),-theta.sin(),0.),p(0.,0.,-1.),
            GeometryTolerance::default()).unwrap();
        let b=source(2,SourceGeometry::FiniteCylinder {
            frame,radius:1.5,axial:ParameterRange::bounded(0.,10.).unwrap(),
        });
        let snap=SourceSnapshot::new([a.clone(),b.clone()]).unwrap();
        let intent=snap.intent(BooleanKind::Difference,a.reference(),b.reference(),
            GeometryTolerance::default(),UnitPolicy::RequireIdentical).unwrap();
        (snap,intent)
    }
    #[test]
    fn g8_source_only_ellipse_bracket_tracks_tilts(){
        for theta in [0.04,-0.03,0.055] {
            let (snap,intent)=snapshot(theta);
            let G8SourceVolumeOutcome::Ready(w)=g8_source_volume(&snap,&intent).unwrap()
                else { panic!("G8b independent source interval held"); };
            assert!(w.retained.lower>0. && w.retained.upper<w.original_box_volume);
            assert!(w.retained.lower<w.retained.upper);
            assert!(w.independent_cross_check>=w.removed_inside_box[0]-1e-7);
            assert!(w.independent_cross_check<=w.removed_inside_box[1]+1e-7);
            assert_eq!(w.source_stamp(),snap.stamp());
            assert_eq!(w.operands()[0].id,intent.operands().0.id);
        }
    }
    #[test]
    fn source_revisions_are_checked_before_volume(){
        let (snap,intent)=snapshot(0.04);
        let original=snap.get(intent.operands().0.id).unwrap();
        let newer=original.successor(original.geometry().clone(),original.tolerance()).unwrap();
        let changed=snap.with_replacement(newer).unwrap();
        assert!(g8_source_volume(&changed,&intent).is_err());
    }
    #[test]
    fn reject_cap_intrusion_and_unsupported_operation(){
        let (snap,intent)=snapshot(0.04);
        let refs=intent.operands();
        let union=snap.intent(BooleanKind::Union,refs.0,refs.1,
            intent.tolerance(),UnitPolicy::RequireIdentical).unwrap();
        assert!(matches!(g8_source_volume(&snap,&union).unwrap(),
            G8SourceVolumeOutcome::Held(G8SourceHold::UnsupportedOperation)));
        let cyl=snap.get(refs.1.id).unwrap();
        let SourceGeometry::FiniteCylinder{frame,radius,..}=cyl.geometry() else{unreachable!()};
        let changed=cyl.successor(SourceGeometry::FiniteCylinder {
            frame:*frame,radius:*radius,axial:ParameterRange::bounded(0.,4.).unwrap(),
        },cyl.tolerance()).unwrap();
        let next=snap.with_replacement(changed).unwrap();
        let intent=next.intent(BooleanKind::Difference,
            next.get(refs.0.id).unwrap().reference(),
            next.get(refs.1.id).unwrap().reference(),
            GeometryTolerance::default(),UnitPolicy::RequireIdentical).unwrap();
        assert!(matches!(g8_source_volume(&next,&intent).unwrap(),
            G8SourceVolumeOutcome::Held(G8SourceHold::CylinderCapIntersectsBox)));
    }
}
