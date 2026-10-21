//! G8c: actual-result boundary flux over a checked analytic B-rep v2.
//! This module consumes NO SourceSpec, G8b brackets, old shell plan or staged mass.
//! The flux is a narrow, source-independent G8 MIN-X family witness, not a
//! complete arbitrary Boolean material/lineage certification.
use rustsolid_foundation::{GeometryTolerance, Point3};
use rustsolid_geometry::{Curve2, Curve3, ParameterRange, Surface3};
use rustsolid_topology::{BrepModel, FaceId, FaceSense, BrepIssue};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G8FluxHold {
    InvalidShell,
    UnsupportedResultFamily,
    MissingPhysicalCoedge,
    IncorrectFaceSense,
    OpenOrPeriodicUvLoop,
    UnsupportedCurve,
    IncorrectFluxOrientation,
    NumericalBudgetExceeded,
}
#[derive(Debug, Clone, PartialEq)]
pub struct G8BoundaryFlux {
    pub retained_volume:f64,
    pub numerical_guard:f64,
    pub original_x_face_flux:f64,
    pub inward_cylinder_flux:f64,
    pub physical_boundary_uv_integral:f64,
    pub harmonic_integrated_uses:usize,
    pub axial_generator_uses:usize,
    pub exterior_x_face:FaceId,
    pub cylinder_face:FaceId,
    pub min_x_faces:[FaceId;2],
    pub model_incarnation:u64,
    /// Length extent bounds numerical checks; not a stored B-rep AABB.
    pub length_scale:f64,
}
#[derive(Debug, Clone, PartialEq)]
pub enum G8FluxOutcome {
    Ready(G8BoundaryFlux),
    Held(G8FluxHold),
}
fn magnitude(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn dist(a:Point3,b:Point3)->f64{magnitude(a.sub(b))}
fn normal(model:&BrepModel, face:FaceId)->Result<Point3,G8FluxHold>{
    let f=model.faces().get(face.0 as usize).ok_or(G8FluxHold::InvalidShell)?;
    let support=model.carriers().surfaces.get(f.surface.0 as usize)
        .ok_or(G8FluxHold::InvalidShell)?;
    let eval=support.evaluate(0.,0.).map_err(|_|G8FluxHold::InvalidShell)?;
    let mut n=eval.du.cross(eval.dv);
    if f.sense==FaceSense::Reversed {n=n.scale(-1.);}
    let norm=magnitude(n);
    if !norm.is_finite()||norm<=0.0 {return Err(G8FluxHold::InvalidShell);}
    Ok(n.scale(1./norm))
}
fn physical_planar_area(model:&BrepModel,face:FaceId)
    ->Result<(f64,f64),G8FluxHold>
{
    let f=model.faces().get(face.0 as usize).ok_or(G8FluxHold::InvalidShell)?;
    if f.loops.len()!=1 {return Err(G8FluxHold::MissingPhysicalCoedge);}
    let lp=model.loops().get(f.loops[0].0 as usize).ok_or(G8FluxHold::MissingPhysicalCoedge)?;
    let first=lp.first;let mut at=first;
    let mut points=Vec::new();
    loop {
        let c=model.coedges().get(at.0 as usize).ok_or(G8FluxHold::MissingPhysicalCoedge)?;
        let edge=model.edges().get(c.edge.0 as usize).ok_or(G8FluxHold::MissingPhysicalCoedge)?;
        let vertex=if c.reversed {edge.vertices[1]} else {edge.vertices[0]};
        points.push(model.vertices().get(vertex.0 as usize)
            .ok_or(G8FluxHold::InvalidShell)?.point);
        at=c.next;
        if at==first {break;}
        if points.len()>12 {return Err(G8FluxHold::UnsupportedResultFamily);}
    }
    if points.len()!=4 {return Err(G8FluxHold::UnsupportedResultFamily);}
    let x=points[0].x;
    if points.iter().any(|q|!q.x.is_finite()||(q.x-x).abs()>1e-8) {
        return Err(G8FluxHold::UnsupportedResultFamily);
    }
    let area=points.iter().enumerate().map(|(i,q)|{
        let b=points[(i+1)%points.len()];
        q.y*b.z-b.y*q.z
    }).sum::<f64>().abs()*0.5;
    if !area.is_finite()||area<=0.0 {return Err(G8FluxHold::UnsupportedResultFamily);}
    Ok((x,area))
}
fn edge_span(v:ParameterRange)->Result<(f64,f64),G8FluxHold>{
    match v {
        ParameterRange::Bounded{start,end}
            if start.is_finite()&&end.is_finite()&&start<end=>Ok((start,end)),
        _=>Err(G8FluxHold::MissingPhysicalCoedge),
    }
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
    fn integral(self,lo:f64,hi:f64)->Result<f64,G8FluxHold>{
        if !lo.is_finite()||!hi.is_finite(){
            return Err(G8FluxHold::NumericalBudgetExceeded);
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
            return Err(G8FluxHold::NumericalBudgetExceeded);
        }
        Ok(out)
    }
}

fn cylinder_coedge_green(
    model:&BrepModel, face:FaceId, xmin:f64,
)->Result<(f64,usize,usize),G8FluxHold>{
    let f=model.faces().get(face.0 as usize).ok_or(G8FluxHold::InvalidShell)?;
    if f.loops.len()!=1 || f.sense!=FaceSense::Reversed {
        return Err(G8FluxHold::IncorrectFaceSense);
    }
    let Surface3::Cylinder(cyl)=model.carriers().surfaces[f.surface.0 as usize]
        else{return Err(G8FluxHold::UnsupportedResultFamily)};
    let radial=Fourier::harmonic(0.,cyl.frame.x.x,cyl.frame.y.x);
    let xoffset=cyl.frame.origin.x-xmin;
    let ax=cyl.frame.normal.x;
    let r=cyl.radius;
    let mut current=model.loops()[f.loops[0].0 as usize].first;
    let first=current;
    let mut sum=0.0;
    let mut harmonic=0usize;let mut generators=0usize;
    let mut points=Vec::new();
    for _ in 0..4{
        let co=model.coedges().get(current.0 as usize)
            .ok_or(G8FluxHold::MissingPhysicalCoedge)?;
        let edge=model.edges().get(co.edge.0 as usize)
            .ok_or(G8FluxHold::MissingPhysicalCoedge)?;
        let (lo,hi)=edge_span(edge.domain)?;
        let curve=*model.carriers().curves2.get(co.pcurve.0 as usize)
            .ok_or(G8FluxHold::MissingPhysicalCoedge)?;
        let start_t=if co.reversed{hi}else{lo};
        let end_t=if co.reversed{lo}else{hi};
        let uv_a=curve.evaluate(co.edge_to_pcurve.map(start_t)
            .map_err(|_|G8FluxHold::MissingPhysicalCoedge)?)
            .map_err(|_|G8FluxHold::MissingPhysicalCoedge)?;
        let uv_b=curve.evaluate(co.edge_to_pcurve.map(end_t)
            .map_err(|_|G8FluxHold::MissingPhysicalCoedge)?)
            .map_err(|_|G8FluxHold::MissingPhysicalCoedge)?;
        points.push((uv_a,uv_b));
        match curve {
            Curve2::CylinderSection(section)=>{
                if section.direction.abs()!=1.0 {
                    return Err(G8FluxHold::UnsupportedCurve);
                }
                let v=Fourier::harmonic(section.axial_center,
                    section.axial_cosine,section.axial_sine);
                let n_v=radial.mul(v);
                let n2_v=radial.mul(radial).mul(v);
                let n_v2=radial.mul(v.mul(v));
                let integrand=n_v.scale(r*xoffset)
                    .add(n2_v.scale(r*r))
                    .add(n_v2.scale(0.5*r*ax));
                sum+=integrand.integral(uv_a.u,uv_b.u)?;
                harmonic+=1;
            }
            Curve2::Line(line)=>{
                if line.delta.u.abs()>1e-12 ||(uv_b.u-uv_a.u).abs()>1e-10 {
                    return Err(G8FluxHold::UnsupportedCurve);
                }
                generators+=1;
            }
            _=>return Err(G8FluxHold::UnsupportedCurve),
        }
        current=co.next;
    }
    if current!=first||harmonic+generators!=4 {
        return Err(G8FluxHold::OpenOrPeriodicUvLoop);
    }
    for i in 0..4{
        let a=points[i].1;let b=points[(i+1)%4].0;
        if (a.u-b.u).abs()>2e-7||(a.v-b.v).abs()>2e-7 {
            return Err(G8FluxHold::OpenOrPeriodicUvLoop);
        }
    }
    Ok((sum,harmonic,generators))
}
pub fn g8_boundary_flux(model:&BrepModel)->G8FluxOutcome {
    use G8FluxHold as Hold;
    if model.validate().is_err() {
        return G8FluxOutcome::Held(Hold::InvalidShell);
    }
    // Narrow physical family: six ordinary box planes + two MIN-X patches,
    // represented as 7 planar faces and one inward cutter face.
    if model.shells().len()!=1||model.faces().len()!=8 ||
        model.vertices().len()!=12||model.edges().len()!=18 ||
        model.coedges().len()!=36 {
        return G8FluxOutcome::Held(Hold::UnsupportedResultFamily);
    }
    let mut minx_faces=Vec::new();
    let mut outer_faces=Vec::new();
    let mut cylinder_faces=Vec::new();
    for face in model.faces(){
        let Some(surface)=model.carriers().surfaces.get(face.surface.0 as usize) else {
            return G8FluxOutcome::Held(Hold::InvalidShell);
        };
        match surface {
            Surface3::Plane(_)=>{
                let Ok(n)=normal(model,face.id) else {
                    return G8FluxOutcome::Held(Hold::IncorrectFaceSense);
                };
                if n.x > 1.0-1e-9 {
                    let Ok((x,area))=physical_planar_area(model,face.id) else {
                        return G8FluxOutcome::Held(Hold::MissingPhysicalCoedge);
                    };
                    outer_faces.push((face.id,x,area));
                }else if n.x < -1.0+1e-9 {
                    let Ok((x,_))=physical_planar_area(model,face.id) else {
                        return G8FluxOutcome::Held(Hold::MissingPhysicalCoedge);
                    };
                    minx_faces.push((face.id,x));
                }else if n.x.abs()>1e-9 {
                    return G8FluxOutcome::Held(Hold::UnsupportedResultFamily);
                }
            },
            Surface3::Cylinder(_)=>cylinder_faces.push(face.id),
            Surface3::Sphere(_)=>return G8FluxOutcome::Held(Hold::UnsupportedResultFamily),
        }
    }
    if outer_faces.len()!=1||minx_faces.len()!=2||cylinder_faces.len()!=1 {
        return G8FluxOutcome::Held(Hold::UnsupportedResultFamily);
    }
    let xmin=minx_faces[0].1;
    if (minx_faces[1].1-xmin).abs()>1e-8 {
        return G8FluxOutcome::Held(Hold::UnsupportedResultFamily);
    }
    let (exterior,xmax,area)=outer_faces[0];
    if xmax<=xmin || !area.is_finite() {return G8FluxOutcome::Held(Hold::IncorrectFluxOrientation);}
    let cylinder=cylinder_faces[0];
    let Ok(n)=normal(model,cylinder) else {
        return G8FluxOutcome::Held(Hold::IncorrectFaceSense);
    };
    if n.x> -0.90 || model.faces()[cylinder.0 as usize].sense!=FaceSense::Reversed{
        return G8FluxOutcome::Held(Hold::IncorrectFaceSense);
    }
    let outer=(xmax-xmin)*area;
    let Ok((boundary,harmonic,generators))=cylinder_coedge_green(model,cylinder,xmin) else {
        return G8FluxOutcome::Held(Hold::MissingPhysicalCoedge);
    };
    let span={
        let mut lo=Point3{x:f64::INFINITY,y:f64::INFINITY,z:f64::INFINITY};
        let mut hi=Point3{x:f64::NEG_INFINITY,y:f64::NEG_INFINITY,z:f64::NEG_INFINITY};
        for v in model.vertices(){
            let q=v.point;
            lo.x=lo.x.min(q.x);lo.y=lo.y.min(q.y);lo.z=lo.z.min(q.z);
            hi.x=hi.x.max(q.x);hi.y=hi.y.max(q.y);hi.z=hi.z.max(q.z);
        }
        magnitude(hi.sub(lo)).max(1.)
    };
    let tol=model.tolerance();
    let eps=tol.length_at(span)*128.0;
    if boundary<=eps || !boundary.is_finite(){
        return G8FluxOutcome::Held(Hold::IncorrectFluxOrientation);
    }
    let volume=outer-boundary;
    if volume<=0.0 || !volume.is_finite() {
        return G8FluxOutcome::Held(Hold::IncorrectFluxOrientation);
    }
    let Surface3::Cylinder(cyl)=model.carriers().surfaces[model.faces()[cylinder.0 as usize].surface.0 as usize]
        else {unreachable!()};
    let scale=span.max(cyl.radius);
    let guard=512.0*f64::EPSILON*(outer.abs()+boundary.abs()+cyl.radius.powi(2)*scale);
    if !guard.is_finite() || guard>tol.length_at(scale)*scale*scale {
        return G8FluxOutcome::Held(Hold::NumericalBudgetExceeded);
    }
    G8FluxOutcome::Ready(G8BoundaryFlux{
        retained_volume:volume,numerical_guard:guard,
        original_x_face_flux:outer,inward_cylinder_flux:-boundary,
        physical_boundary_uv_integral:boundary,
        harmonic_integrated_uses:harmonic,axial_generator_uses:generators,
        exterior_x_face:exterior,cylinder_face:cylinder,
        min_x_faces:[minx_faces[0].0,minx_faces[1].0],
        model_incarnation:model.incarnation().get(),
        length_scale:scale,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finite_fourier_integrates_one_cycle_without_false_dc() {
        let u=Fourier::harmonic(2.0,0.7,-0.9);
        let integral=u.integral(0.,std::f64::consts::TAU).unwrap();
        assert!((integral-2.0*std::f64::consts::TAU).abs()<1e-12);
    }
}
