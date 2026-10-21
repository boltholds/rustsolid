//! Independent G8c oriented flux from the checked, actual analytic BrepModel.
//! The result is not a Boolean certificate until separately compared with
//! source-only bounds and complete material/provenance evidence.
use rustsolid_foundation::{GeometryError,GeometryTolerance,Point3};
use rustsolid_geometry::{Curve2,ParameterRange,Surface3};
use rustsolid_topology::{BrepModel,FaceId,FaceSense,EdgeKind};

#[derive(Debug,Clone,PartialEq)]
pub enum BoundaryFluxHold {
    InvalidModel, UnsupportedFamily, InvalidOrientation, MissingCoedge,
    InvalidCarrier, OpenUvLoop, InvalidFlux, NumericalBudgetExceeded,
}
impl From<GeometryError> for BoundaryFluxHold {
    fn from(_:GeometryError)->Self{Self::InvalidCarrier}
}
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn bound(r:ParameterRange)->Result<(f64,f64),BoundaryFluxHold>{
    match r {
        ParameterRange::Bounded{start,end} if start.is_finite()&&end.is_finite()&&end>start=>Ok((start,end)),
        _=>Err(BoundaryFluxHold::InvalidCarrier),
    }
}
fn length(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
#[derive(Debug,Clone,PartialEq)]
pub struct BoundaryFluxEvidence {
    volume:f64,
    guard:f64,
    exterior_x_face:FaceId,
    inner_cylinder_face:FaceId,
    outward_x_flux:f64,
    inward_cylinder_flux:f64,
    harmonics:usize,
    generators:usize,
}
impl BoundaryFluxEvidence {
    pub fn result_volume(&self)->f64{self.volume}
    pub fn numerical_guard(&self)->f64{self.guard}
    pub fn exterior_x_face(&self)->FaceId{self.exterior_x_face}
    pub fn inner_cylinder_face(&self)->FaceId{self.inner_cylinder_face}
    pub fn outward_x_flux(&self)->f64{self.outward_x_flux}
    pub fn inward_cylinder_flux(&self)->f64{self.inward_cylinder_flux}
    pub fn harmonic_uses(&self)->usize{self.harmonics}
    pub fn generator_uses(&self)->usize{self.generators}
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

fn cylinder_green(
    m:&BrepModel,face:FaceId,xmin:f64,
)->Result<(f64,usize,usize),BoundaryFluxHold>{
    let f=m.faces().get(face.0 as usize).ok_or(BoundaryFluxHold::InvalidModel)?;
    let Surface3::Cylinder(cyl)=m.carriers().surfaces[f.surface.0 as usize]
        else{return Err(BoundaryFluxHold::UnsupportedFamily);};
    if f.sense!=FaceSense::Reversed || f.loops.len()!=1{
        return Err(BoundaryFluxHold::InvalidOrientation);
    }
    let radial=Fourier::harmonic(0.,cyl.frame.x.x,cyl.frame.y.x);
    let r=cyl.radius;
    let xoffset=cyl.frame.origin.x-xmin;
    let ax=cyl.frame.normal.x;
    let first=m.loops()[f.loops[0].0 as usize].first;
    let mut cid=first;
    let mut sum=0.0;
    let mut harmonic=0usize;
    let mut generator=0usize;
    let mut uv_ends=Vec::new();
    for _ in 0..4 {
        let co=m.coedges().get(cid.0 as usize).ok_or(BoundaryFluxHold::MissingCoedge)?;
        if co.face!=face||co.loop_id!=f.loops[0]{
            return Err(BoundaryFluxHold::MissingCoedge);
        }
        let edge=m.edges().get(co.edge.0 as usize).ok_or(BoundaryFluxHold::MissingCoedge)?;
        if edge.kind!=EdgeKind::Physical{return Err(BoundaryFluxHold::MissingCoedge);}
        let (lo,hi)=bound(edge.domain)?;
        let start=if co.reversed{hi}else{lo};
        let end=if co.reversed{lo}else{hi};
        let t0=co.edge_to_pcurve.map(start).map_err(|_|BoundaryFluxHold::InvalidCarrier)?;
        let t1=co.edge_to_pcurve.map(end).map_err(|_|BoundaryFluxHold::InvalidCarrier)?;
        let trim=*m.carriers().curves2.get(co.pcurve.0 as usize).ok_or(BoundaryFluxHold::InvalidCarrier)?;
        let a=trim.evaluate(t0)?;
        let b=trim.evaluate(t1)?;
        uv_ends.push((a,b));
        match trim {
            Curve2::CylinderSection(section)=>{
                if section.direction.abs()!=1. {return Err(BoundaryFluxHold::InvalidCarrier);}
                let v=Fourier::harmonic(section.axial_center,section.axial_cosine,section.axial_sine);
                let n_v=radial.mul(v);
                let n2_v=radial.mul(radial).mul(v);
                let n_v2=radial.mul(v.mul(v));
                let integrand=n_v.scale(r*xoffset)
                    .add(n2_v.scale(r*r))
                    .add(n_v2.scale(0.5*r*ax));
                sum+=integrand.integral(a.u,b.u)?;
                harmonic+=1;
            },
            Curve2::Line(line)=>{
                if line.delta.u.abs()>1e-12||(b.u-a.u).abs()>1e-10{
                    return Err(BoundaryFluxHold::InvalidCarrier);
                }
                generator+=1;
            },
            _=>return Err(BoundaryFluxHold::InvalidCarrier),
        }
        cid=co.next;
    }
    if cid!=first || harmonic+generator!=4 {return Err(BoundaryFluxHold::OpenUvLoop);}
    for i in 0..4{
        let b=uv_ends[i].1;
        let a=uv_ends[(i+1)%4].0;
        if (a.u-b.u).abs()>2e-7 ||(a.v-b.v).abs()>2e-7{
            return Err(BoundaryFluxHold::OpenUvLoop);
        }
    }
    Ok((sum,harmonic,generator))
}

fn physical_planar_area(m:&BrepModel, face:FaceId)->Result<(f64,f64),BoundaryFluxHold>{
    let f=m.faces().get(face.0 as usize).ok_or(BoundaryFluxHold::InvalidModel)?;
    if f.loops.len()!=1{return Err(BoundaryFluxHold::UnsupportedFamily);}
    let first=m.loops()[f.loops[0].0 as usize].first;
    let mut at=first;
    let mut points=Vec::new();
    loop{
        let co=m.coedges().get(at.0 as usize).ok_or(BoundaryFluxHold::MissingCoedge)?;
        let edge=m.edges().get(co.edge.0 as usize).ok_or(BoundaryFluxHold::MissingCoedge)?;
        if co.face!=face ||co.loop_id!=f.loops[0]{return Err(BoundaryFluxHold::MissingCoedge);}
        let vertex=if co.reversed{edge.vertices[1]}else{edge.vertices[0]};
        points.push(m.vertices()[vertex.0 as usize].point);
        at=co.next;
        if at==first{break;}
        if points.len()>12{return Err(BoundaryFluxHold::UnsupportedFamily);}
    }
    if points.len()!=4 ||points.iter().any(|p|(p.x-points[0].x).abs()>1e-8){
        return Err(BoundaryFluxHold::UnsupportedFamily);
    }
    let area=points.iter().enumerate().map(|(i,p)|{
        let q=points[(i+1)%points.len()];
        p.y*q.z-q.y*p.z
    }).sum::<f64>().abs()*0.5;
    if !area.is_finite()||area<=0.{return Err(BoundaryFluxHold::InvalidFlux);}
    Ok((points[0].x,area))
}

/// Integrate the ACTUAL coedge curves of a source-qualified G8 shape.
/// x_min is the source reference plane for F=(x-x_min, 0, 0), div F=1.
/// The function never reads an expected result mass or G8b polygon integral.
pub fn evaluate_tilted_boundary_flux(
    m:&BrepModel,x_min:f64,tol:GeometryTolerance,
)->Result<BoundaryFluxEvidence,BoundaryFluxHold>{
    tol.validate().map_err(|_|BoundaryFluxHold::NumericalBudgetExceeded)?;
    m.validate().map_err(|_|BoundaryFluxHold::InvalidModel)?;
    if !x_min.is_finite()||m.shells().len()!=1||m.faces().len()!=8
        ||m.vertices().len()!=12||m.edges().len()!=18||m.coedges().len()!=36 {
        return Err(BoundaryFluxHold::UnsupportedFamily);
    }
    let mut vmin=Point3{x:f64::INFINITY,y:f64::INFINITY,z:f64::INFINITY};
    let mut vmax=Point3{x:f64::NEG_INFINITY,y:f64::NEG_INFINITY,z:f64::NEG_INFINITY};
    for v in m.vertices(){
        let q=v.point;
        vmin.x=vmin.x.min(q.x);vmin.y=vmin.y.min(q.y);vmin.z=vmin.z.min(q.z);
        vmax.x=vmax.x.max(q.x);vmax.y=vmax.y.max(q.y);vmax.z=vmax.z.max(q.z);
    }
    let scale=length(vmax.sub(vmin)).max(1.);
    let eps=tol.length_at(scale)*128.;
    if (vmin.x-x_min).abs()>eps{return Err(BoundaryFluxHold::UnsupportedFamily);}
    let mut cylinder=None;
    let mut exterior=None;
    let mut xmin_faces=0usize;
    let mut outward_flux=0.;
    for f in m.faces(){
        let surface=m.carriers().surfaces[f.surface.0 as usize];
        match surface{
            Surface3::Plane(plane)=>{
                let mut n=plane.frame.normal;
                if f.sense==FaceSense::Reversed{n=n.scale(-1.);}
                if n.x.abs()<1e-8{continue;}
                if (n.x.abs()-1.).abs()>1e-8||n.y.abs()>1e-8||n.z.abs()>1e-8{
                    return Err(BoundaryFluxHold::UnsupportedFamily);
                }
                let (x,area)=physical_planar_area(m,f.id)?;
                if n.x>0. {
                    if exterior.replace(f.id).is_some()||x<=x_min+eps{
                        return Err(BoundaryFluxHold::UnsupportedFamily);
                    }
                    outward_flux+=(x-x_min)*area;
                }else if (x-x_min).abs()<=eps{
                    xmin_faces+=1;
                }else{
                    return Err(BoundaryFluxHold::InvalidOrientation);
                }
            }
            Surface3::Cylinder(cyl)=>{
                if cylinder.replace(f.id).is_some() || f.sense!=FaceSense::Reversed{
                    return Err(BoundaryFluxHold::InvalidOrientation);
                }
                let d=Surface3::Cylinder(cyl).evaluate(0.,0.)?;
                let n=d.du.cross(d.dv).scale(-1.);
                if length(n)<=0.||n.x/length(n)>-0.90{
                    return Err(BoundaryFluxHold::InvalidOrientation);
                }
            }
            Surface3::Sphere(_)=>return Err(BoundaryFluxHold::UnsupportedFamily),
        }
    }
    if exterior.is_none()||cylinder.is_none()||xmin_faces!=2{
        return Err(BoundaryFluxHold::UnsupportedFamily);
    }
    let c=cylinder.expect("checked");
    let (boundary,harmonics,generators)=cylinder_green(m,c,x_min)?;
    if !boundary.is_finite()||boundary<=eps{return Err(BoundaryFluxHold::InvalidFlux);}
    let volume=outward_flux-boundary;
    if !volume.is_finite()||volume<=0.{return Err(BoundaryFluxHold::InvalidFlux);}
    let Surface3::Cylinder(support)=m.carriers().surfaces[m.faces()[c.0 as usize].surface.0 as usize]
        else{unreachable!()};
    let guard=512.*f64::EPSILON*(outward_flux.abs()+boundary.abs()+support.radius.powi(2)*scale);
    if !guard.is_finite()||guard<0.||guard>tol.length_at(scale)*scale*scale{
        return Err(BoundaryFluxHold::NumericalBudgetExceeded);
    }
    Ok(BoundaryFluxEvidence{
        volume,guard,exterior_x_face:exterior.expect("checked"),
        inner_cylinder_face:c,outward_x_flux:outward_flux,
        inward_cylinder_flux:-boundary,harmonics,generators,
    })
}
