//! G8c independent oriented flux of the ACTUAL tilted result BrepModel.
//! This module never integrates the original cylinder/ellipse footprints,
//! invokes G8b's analytic area formula or trusts staged candidate mass.
//!
//! F=(x-xmin,0,0), div F=1. Exterior Xmax contributes its real B-rep
//! oriented planar coedge-loop area. Other box faces have zero flux.
//! The inward result-cylinder flux comes from a 1D Green integral along
//! the four original physical cylinder-face UV coedges. Their exact
//! CylinderSection2 carriers are trigonometric polynomials; multiply
//! Fourier coefficients and integrate harmonics without mesh facets.

use std::f64::consts::TAU;
use crate::{
    BrepBody,BrepModel,Curve2,FaceId,FaceSense,GeometryError,
    GeometryTolerance,IndependentTiltedVolumeOracle,ParameterRange,
    Point3,Surface3,TiltedOneWallShellPlan,
};

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn bound(r:ParameterRange)->Result<(f64,f64),GeometryError>{
    match r{
        ParameterRange::Bounded{start,end}
            if start.is_finite()&&end.is_finite()&&start<end=>Ok((start,end)),
        _=>Err(invalid("G8c flux requires bounded physical UV carriers")),
    }
}
fn magnitude(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn distance(a:Point3,b:Point3)->f64{magnitude(a.sub(b))}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TiltedFluxHold {
    InvalidShell,
    IncorrectFaceSense,
    MissingPhysicalCoedge,
    OpenOrPeriodicUvLoop,
    UnsupportedCurve,
    IncorrectFluxOrientation,
    SourceFluxDisagreement,
    ResultMassDisagreement,
    NumericalBudgetExceeded,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TiltedBoundaryFluxWitness {
    pub exterior_x_face:FaceId,
    pub retained_x_faces:[FaceId;2],
    pub inward_cylinder_face:FaceId,
    pub original_x_face_flux:f64,
    pub inward_cylinder_flux:f64,
    pub result_volume:f64,
    pub physical_boundary_uv_integral:f64,
    pub harmonic_integrated_uses:usize,
    pub vertical_generator_uses:usize,
    pub numerical_guard:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct IndependentTiltedBoundaryFluxOracle {
    pub witness:TiltedBoundaryFluxWitness,
    pub source_retained:[f64;2],
    pub candidate_staged_volume:f64,
    pub tolerance:GeometryTolerance,
    pub volume_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum IndependentTiltedBoundaryFluxResult{
    Ready(IndependentTiltedBoundaryFluxOracle),
    Held(TiltedFluxHold),
}

// Small finite Fourier algebra for harmonics k=-3..3.
// Input n(u) has degree 1, v(u) has degree 1, and
// G(u,v(u))=R n[(Ox-xmin+R n)*v + Ax*v^2/2] has degree <=3.
// No quadrature samples or source-section antiderivatives are used.
#[derive(Debug,Clone,Copy,Default)]
struct Complex{re:f64,im:f64}
impl Complex{
    fn mul(self,other:Self)->Self{
        Self{
            re:self.re*other.re-self.im*other.im,
            im:self.re*other.im+self.im*other.re,
        }
    }
}
#[derive(Debug,Clone,Copy,Default)]
struct Fourier{a:[Complex;7]}
impl Fourier{
    fn scalar(x:f64)->Self{
        let mut z=Self::default();z.a[3].re=x;z
    }
    fn harmonic(a:f64,cos:f64,sin:f64)->Self{
        let mut z=Self::scalar(a);
        z.a[4]=Complex{re:cos/2.,im:-sin/2.};
        z.a[2]=Complex{re:cos/2.,im:sin/2.};
        z
    }
    fn scale(mut self,s:f64)->Self{
        for c in &mut self.a{c.re*=s;c.im*=s;}
        self
    }
    fn add(mut self,b:Self)->Self{
        for i in 0..7{
            self.a[i].re+=b.a[i].re;
            self.a[i].im+=b.a[i].im;
        }
        self
    }
    fn mul(self,b:Self)->Self{
        let mut z=Self::default();
        for i in -3i32..=3{
            for j in -3i32..=3{
                if (-3..=3).contains(&(i+j)){
                    let product=self.a[(i+3) as usize]
                        .mul(b.a[(j+3) as usize]);
                    let out=&mut z.a[(i+j+3) as usize];
                    out.re+=product.re;
                    out.im+=product.im;
                }
            }
        }
        z
    }
    fn integral(self,lo:f64,hi:f64)->Result<f64,GeometryError>{
        if !lo.is_finite()||!hi.is_finite(){
            return Err(invalid("G8c nonfinite harmonic integration span"));
        }
        let mut out=self.a[3].re*(hi-lo);
        let mid=0.5*lo+0.5*hi;
        let half_delta=0.5*(hi-lo);
        for k in 1..=3{
            let kf=k as f64;
            let twice=2.*(kf*half_delta).sin();
            let ds=twice*(kf*mid).cos();
            let dc=-twice*(kf*mid).sin();
            let c=self.a[3+k];
            out+=2.*(c.re*ds+c.im*dc)/kf;
        }
        if !out.is_finite(){
            return Err(invalid("G8c harmonic flux evaluation nonfinite"));
        }
        Ok(out)
    }
}
fn cylinder_coedge_green(
    model:&BrepModel,face:FaceId,xmin:f64,
)->Result<(f64,usize,usize),GeometryError>{
    let Surface3::Cylinder(cylinder)=model.face_surface(face)
        .ok_or_else(||invalid("G8c missing actual cylinder face support"))?
        else{return Err(invalid("G8c inner face is not analytic cylinder"));};
    if model.face_sense(face)!=Some(FaceSense::Reversed){
        return Err(invalid("G8c result cylinder must be inward-oriented"));
    }
    let radial=Fourier::harmonic(
        0.,cylinder.frame.x.x,cylinder.frame.y.x);
    let r=cylinder.radius;
    let xoffset=cylinder.frame.origin.x-xmin;
    let ax=cylinder.frame.normal.x;
    let mut current=model.loops[face.0 as usize].first_coedge;
    let first=current;
    let mut sum=0.;
    let mut harmonic_count=0usize;
    let mut generator_count=0usize;
    let mut points=Vec::new();
    for _ in 0..4{
        let geom=model.geometry.coedges.get(current.0 as usize)
            .ok_or_else(||invalid("G8c physical cylinder coedge absent"))?;
        let (lo,hi)=bound(geom.domain)?;
        let curve=model.coedge_pcurve(current)
            .ok_or_else(||invalid("G8c missing actual result pcurve"))?;
        let uva=curve.evaluate(lo)?;
        let uvb=curve.evaluate(hi)?;
        points.push((uva,uvb));
        match curve{
            Curve2::CylinderSection(section)=>{
                if section.direction.abs()!=1.0 {
                    return Err(invalid("G8c nonunit angular pcurve direction"));
                }
                let v=Fourier::harmonic(
                    section.axial_center,section.axial_cosine,section.axial_sine);
                let n_v=radial.mul(v);
                let n2_v=radial.mul(radial).mul(v);
                let n_v2=radial.mul(v.mul(v));
                let integrand=n_v.scale(r*xoffset)
                    .add(n2_v.scale(r*r))
                    .add(n_v2.scale(0.5*r*ax));
                sum+=integrand.integral(uva.u,uvb.u)?;
                harmonic_count+=1;
            },
            Curve2::Line(line)=>{
                // Physical generator at constant cylinder angle: du=0,
                // so it contributes zero to G(u,v)*du, regardless of tilt.
                if line.delta.u.abs()>1e-12||(uvb.u-uva.u).abs()>1e-10{
                    return Err(invalid("G8c wall line is not an axial generator"));
                }
                generator_count+=1;
            },
            _=>return Err(invalid("G8c cylinder UV uses unverified curve")),
        }
        current=model.coedges[current.0 as usize].next;
    }
    if current!=first||points.len()!=4||
        harmonic_count+generator_count!=4{
        return Err(invalid("G8c result cylinder has incomplete boundary"));
    }
    for i in 0..4{
        let a=points[i].1;
        let b=points[(i+1)%4].0;
        if (a.u-b.u).abs()>2e-7||(a.v-b.v).abs()>2e-7{
            return Err(invalid("G8c physical cylinder UV boundary does not close"));
        }
    }
    Ok((sum,harmonic_count,generator_count))
}
fn physical_planar_flux(m:&BrepModel,face:FaceId,cx:f64)
    ->Result<f64,GeometryError>
{
    let model=m.faces.get(face.0 as usize)
        .ok_or_else(||invalid("G8c missing exterior result face"))?;
    if model.loops.len()!=1{
        return Err(invalid("G8c exterior planar face has nontrivial loops"));
    }
    let normal=m.oriented_face_normal(face,0.,0.)?;
    if normal.x.abs()<1.-1e-9 ||normal.y.abs()>1e-9||normal.z.abs()>1e-9{
        return Err(invalid("G8c unexpected orientation of X wall"));
    }
    let initial=m.loops[model.loops[0].0 as usize].first_coedge;
    let mut at=initial;
    let mut yz=Vec::new();
    loop {
        let co=&m.coedges[at.0 as usize];
        let edge=&m.edges[co.edge.0 as usize];
        let id=if co.reversed{edge.end}else{edge.start};
        let q=m.vertices[id.0 as usize].position;
        yz.push(q);
        at=co.next;
        if at==initial{break;}
        if yz.len()>12{return Err(invalid("G8c unbounded planar loop"));}
    }
    if yz.len()!=4 ||yz.iter().any(|q|(q.x-yz[0].x).abs()>1e-8){
        return Err(invalid("G8c exterior source wall lacks physical rectangle"));
    }
    let area=yz.iter().enumerate().map(|(i,q)|{
        let next=yz[(i+1)%yz.len()];
        q.y*next.z-next.y*q.z
    }).sum::<f64>().abs()*0.5;
    if area<=0.||!area.is_finite(){
        return Err(invalid("G8c physical wall area invalid"));
    }
    Ok((yz[0].x-cx)*normal.x*area)
}
fn integrate(
    model:&BrepModel,tol:GeometryTolerance,
)->Result<TiltedBoundaryFluxWitness,TiltedFluxHold>{
    model.validate().map_err(|_|TiltedFluxHold::InvalidShell)?;
    if model.origin!=crate::BrepOrigin::AnalyticTiltedOneWallCylinderDifference
        ||model.shells.len()!=1||model.faces.len()!=8
        ||model.vertices.len()!=12||model.edges.len()!=18{
        return Err(TiltedFluxHold::InvalidShell);
    }
    let xmin=model.vertices[0].position.x;
    let eps=tol.length_at(model.bbox.max.sub(model.bbox.min)
        .length_squared().sqrt().max(1.))*128.;
    for index in [0u32,1,2,4,5,6]{
        let n=model.oriented_face_normal(FaceId(index),0.,0.)
            .map_err(|_|TiltedFluxHold::IncorrectFaceSense)?;
        if n.x.abs()>1e-8 &&index!=5&&index!=6 ||
            (index==5||index==6)&&n.x> -1.+1e-8 {
            return Err(TiltedFluxHold::IncorrectFaceSense);
        }
    }
    let cylinder=FaceId(7);
    let cyl_surface=model.face_surface(cylinder)
        .ok_or(TiltedFluxHold::InvalidShell)?;
    let Surface3::Cylinder(cyl)=cyl_surface else{
        return Err(TiltedFluxHold::InvalidShell);
    };
    let n=model.oriented_face_normal(cylinder,0.,0.)
        .map_err(|_|TiltedFluxHold::IncorrectFaceSense)?;
    if n.x> -0.90 ||model.face_sense(cylinder)!=Some(FaceSense::Reversed){
        return Err(TiltedFluxHold::IncorrectFaceSense);
    }
    let outer=physical_planar_flux(model,FaceId(3),xmin)
        .map_err(|_|TiltedFluxHold::MissingPhysicalCoedge)?;
    let (boundary,harmonic,lines)=
        cylinder_coedge_green(model,cylinder,xmin)
            .map_err(|_|TiltedFluxHold::MissingPhysicalCoedge)?;
    if boundary<=eps ||!boundary.is_finite(){
        return Err(TiltedFluxHold::IncorrectFluxOrientation);
    }
    let flux=outer-boundary;
    if !flux.is_finite()||flux<=0.{
        return Err(TiltedFluxHold::IncorrectFluxOrientation);
    }
    let extent=model.bbox.max.sub(model.bbox.min);
    let scale=extent.length_squared().sqrt()
        .max(cyl.radius).max(1.);
    let guard=512.*f64::EPSILON*
        (outer.abs()+boundary.abs()+cyl.radius.powi(2)*scale);
    if !guard.is_finite()||guard>tol.length_at(scale)*scale*scale{
        return Err(TiltedFluxHold::NumericalBudgetExceeded);
    }
    Ok(TiltedBoundaryFluxWitness{
        exterior_x_face:FaceId(3),retained_x_faces:[FaceId(5),FaceId(6)],
        inward_cylinder_face:cylinder,original_x_face_flux:outer,
        inward_cylinder_flux:-boundary,result_volume:flux,
        physical_boundary_uv_integral:boundary,
        harmonic_integrated_uses:harmonic,vertical_generator_uses:lines,
        numerical_guard:guard,
    })
}
impl IndependentTiltedBoundaryFluxOracle {
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        shell:&TiltedOneWallShellPlan,
        tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        shell.validate_against(first,second,tol)?;
        if self.tolerance!=tol ||
            Self::build(first,second,shell,tol)?
                !=IndependentTiltedBoundaryFluxResult::Ready(self.clone()){
            return Err(invalid("G8c independent B-rep flux evidence changed"));
        }
        Ok(())
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        shell:&TiltedOneWallShellPlan,tol:GeometryTolerance,
    )->Result<IndependentTiltedBoundaryFluxResult,GeometryError>{
        tol.validate()?;
        shell.validate_against(first,second,tol)?;
        let model=shell.brep();
        let witness=match integrate(model,tol){
            Ok(x)=>x,
            Err(hold)=>return Ok(IndependentTiltedBoundaryFluxResult::Held(hold)),
        };
        let IndependentTiltedVolumeOracle{retained,..}=&shell.source_volume;
        let span=model.bbox.max.sub(model.bbox.min)
            .length_squared().sqrt().max(1.);
        let budget=4.*tol.length_at(span)*span*span;
        let slack=budget+witness.numerical_guard;
        if !budget.is_finite()||budget<=0.{
            return Ok(IndependentTiltedBoundaryFluxResult::Held(
                TiltedFluxHold::NumericalBudgetExceeded));
        }
        if witness.result_volume<retained.lower-slack ||
            witness.result_volume>retained.upper+slack{
            return Ok(IndependentTiltedBoundaryFluxResult::Held(
                TiltedFluxHold::SourceFluxDisagreement));
        }
        if model.mass.volume<retained.lower ||
            model.mass.volume>retained.upper ||
            (model.mass.volume-witness.result_volume).abs()>slack{
            return Ok(IndependentTiltedBoundaryFluxResult::Held(
                TiltedFluxHold::ResultMassDisagreement));
        }
        Ok(IndependentTiltedBoundaryFluxResult::Ready(Self{
            witness,source_retained:[retained.lower,retained.upper],
            candidate_staged_volume:model.mass.volume,
            tolerance:tol,volume_certified:true,can_apply:false,
        }))
    }
}
