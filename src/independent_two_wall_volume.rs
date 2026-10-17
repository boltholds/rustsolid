//! G6b: independently certify each of two curved material shell volumes.
//!
//! Path A integrates strictly convex source material widths with explicit
//! midpoint-lower and trapezoid-upper numerical brackets. Path B traverses
//! the actual result B-rep coedges and evaluates oriented boundary flux
//! with F=(0,0,z-cylinder_center_z) (div F=1). No triangulation, candidate
//! analytic antiderivative or stored mass is used for either path.
//! The resulting limited-family certificate is NOT a formal rounded-interval
//! arithmetic proof for arbitrary geometry.

use std::collections::BTreeSet;
use crate::{
    AnalyticTwoWallPartitionPlan,BrepBody,BrepModel,BrepOrigin,
    GeometryError,GeometryTolerance,SourceVolumeBracket,
    Curve3,FaceId,FaceSense,ParameterRange,Point3,ShellId,Surface3,
};
fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn bounds(range:ParameterRange)->Option<[f64;2]>{
    match range{
        ParameterRange::Bounded{start,end} if start.is_finite()&&end.is_finite()
            &&end>start=>Some([start,end]),
        _=>None,
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TwoWallVolumeHold {
    UnsupportedOperand,
    UnsupportedCandidate,
    InvalidSourceMass,
    UnsupportedBoundary,
    MissingPhysicalTrim,
    IncorrectFaceOrientation,
    InvalidNumericalBudget,
    UnresolvedSourceIntegration,
    ShellFluxDisagreement,
    CandidateVolumeDisagreement,
    GlobalMassDisagreement,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TwoWallFluxWitness {
    pub shell:ShellId,
    pub result_volume:f64,
    pub z_wall_face:FaceId,
    pub cylinder_face:FaceId,
    pub source_z_wall_flux:f64,
    pub inward_cylinder_flux:f64,
    pub angular_span:f64,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TwoWallVolumeComponent {
    pub shell:ShellId,
    pub z_side:i8,
    pub source:SourceVolumeBracket,
    pub boundary:TwoWallFluxWitness,
    pub candidate_volume:f64,
    pub allowed_error:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct AnalyticTwoWallVolumeOracle {
    pub components:Vec<TwoWallVolumeComponent>,
    pub source_total:SourceVolumeBracket,
    pub boundary_total:f64,
    pub stored_mass:f64,
    pub original_box_volume:f64,
    pub removed_inside_box:[f64;2],
    pub source_revisions:[Option<u64>;2],
    pub tolerance:GeometryTolerance,
    pub volume_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum AnalyticTwoWallVolumeResult {
    Ready(AnalyticTwoWallVolumeOracle),
    Held(TwoWallVolumeHold),
}
#[derive(Default)]
struct Compensated{sum:f64,error:f64}
impl Compensated{
    fn add(&mut self,value:f64){
        let term=value-self.error;let next=self.sum+term;
        self.error=(next-self.sum)-term;self.sum=next;
    }
}
/// w(x)=b-sqrt(r^2-x^2), |x|<r, b>r.
/// w''(x)=r^2/(r^2-x^2)^(3/2)>0: midpoint lower/trapezoid upper
/// in exact arithmetic. A conservative f64 margin must also fit.
fn source_bracket(
    r:f64,xa:f64,xb:f64,b:f64,height:f64,target:f64,
)->Result<SourceVolumeBracket,TwoWallVolumeHold>{
    if ![r,xa,xb,b,height,target].iter().all(|x|x.is_finite())
        ||r<=0.||height<=0.||target<=0.||b<=r
        ||xa>=xb||xa.abs()>=r||xb.abs()>=r {
        return Err(TwoWallVolumeHold::UnsupportedOperand);
    }
    let span=xb-xa;
    let max_x=xa.abs().max(xb.abs());
    let min_radius=((r-max_x)*(r+max_x)).sqrt();
    if !min_radius.is_finite()||min_radius<=0. {
        return Err(TwoWallVolumeHold::InvalidNumericalBudget);
    }
    let margin=256.0*f64::EPSILON
        *(r+b+r*r/min_radius)*span*height;
    if !margin.is_finite()||margin>=target/4. {
        return Err(TwoWallVolumeHold::InvalidNumericalBudget);
    }
    let eval=|x:f64|->Result<f64,TwoWallVolumeHold>{
        let radial=(r-x)*(r+x);
        if !radial.is_finite()||radial<0. {
            return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
        }
        let width=b-radial.sqrt();
        if !width.is_finite()||width<=0. {
            return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
        }
        Ok(width)
    };
    let mut lower=Compensated::default();
    let mut upper=Compensated::default();
    let mut slices=vec![(xa,xb,0usize)];
    let mut panels=0usize;
    const MAX_PANELS:usize=131072;
    while let Some((lo,hi,depth))=slices.pop(){
        let mid=0.5*lo+0.5*hi;
        let length=hi-lo;
        if mid<=lo||mid>=hi||!length.is_finite(){
            return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
        }
        let left=eval(lo)?;
        let center=eval(mid)?;
        let right=eval(hi)?;
        let minimum=length*center;
        let maximum=0.5*length*(left+right);
        if maximum+margin/height<minimum{
            return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
        }
        let local_budget=(target-margin)*0.5*length/span/height;
        if maximum-minimum<=local_budget {
            lower.add(minimum);
            upper.add(maximum.max(minimum));
            panels+=1;
            if panels>MAX_PANELS {
                return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
            }
        }else{
            if depth>=28||panels+slices.len()>=MAX_PANELS{
                return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
            }
            slices.push((mid,hi,depth+1));
            slices.push((lo,mid,depth+1));
        }
    }
    let lo=(lower.sum*height-margin).max(0.);
    let hi=upper.sum*height+margin;
    if !lo.is_finite()||!hi.is_finite()||hi<=lo||hi-lo>target{
        return Err(TwoWallVolumeHold::UnresolvedSourceIntegration);
    }
    Ok(SourceVolumeBracket{
        lower:lo,upper:hi,integration_panels:panels,
        numerical_guard:margin,
    })
}
fn face_uses(m:&BrepModel,id:FaceId)->Result<Vec<usize>,TwoWallVolumeHold>{
    let face=m.faces.get(id.0 as usize).filter(|f|f.id==id&&f.loops.len()==1)
        .ok_or(TwoWallVolumeHold::UnsupportedBoundary)?;
    let lp=m.loops.get(face.loops[0].0 as usize)
        .ok_or(TwoWallVolumeHold::MissingPhysicalTrim)?;
    let mut ids=Vec::new();
    let mut current=lp.first_coedge;
    loop {
        let co=m.coedges.get(current.0 as usize)
            .filter(|c|c.face==id&&c.loop_id==lp.id)
            .ok_or(TwoWallVolumeHold::MissingPhysicalTrim)?;
        if ids.contains(&(current.0 as usize)){
            return Err(TwoWallVolumeHold::MissingPhysicalTrim);
        }
        ids.push(current.0 as usize);
        current=co.next;
        if current==lp.first_coedge {break;}
        if ids.len()>16{return Err(TwoWallVolumeHold::UnsupportedBoundary);}
    }
    if ids.len()<3{return Err(TwoWallVolumeHold::MissingPhysicalTrim);}
    Ok(ids)
}
fn shell_flux(
    m:&BrepModel,shell:ShellId,center_z:f64,eps:f64,
)->Result<TwoWallFluxWitness,TwoWallVolumeHold>{
    let owning=m.shells.get(shell.0 as usize)
        .filter(|s|s.id==shell&&s.faces.len()==6&&s.closed)
        .ok_or(TwoWallVolumeHold::UnsupportedBoundary)?;
    let mut z_wall=None;let mut cylinder_wall=None;
    let mut other=0usize;
    let mut plane_flux=0.;
    let mut cylinder_flux=0.;
    let mut angle_span=0.;
    for &fid in &owning.faces{
        let support=m.face_surface(fid)
            .ok_or(TwoWallVolumeHold::UnsupportedBoundary)?;
        let uses=face_uses(m,fid)?;
        let points=uses.iter().map(|&u|{
            let co=&m.coedges[u];
            let e=&m.edges[co.edge.0 as usize];
            let id=if co.reversed{e.end}else{e.start};
            m.vertices[id.0 as usize].position
        }).collect::<Vec<_>>();
        match support{
            Surface3::Plane(_)=>{
                let n=m.oriented_face_normal(fid,0.,0.)
                    .map_err(|_|TwoWallVolumeHold::IncorrectFaceOrientation)?;
                if n.z.abs()>1.-1e-9{
                    if uses.len()!=4||z_wall.replace(fid).is_some()
                        ||n.x.abs()>1e-9||n.y.abs()>1e-9
                        ||points.iter().any(|p|(p.z-points[0].z).abs()>eps){
                        return Err(TwoWallVolumeHold::IncorrectFaceOrientation);
                    }
                    let twice_area=points.iter().enumerate().map(|(i,p)|{
                        let next=points[(i+1)%points.len()];
                        p.x*next.y-next.x*p.y
                    }).sum::<f64>();
                    let area=0.5*twice_area.abs();
                    if !area.is_finite()||area<=eps*eps{
                        return Err(TwoWallVolumeHold::MissingPhysicalTrim);
                    }
                    plane_flux+=(points[0].z-center_z)*n.z*area;
                }else{
                    if n.z.abs()>1e-9 {
                        return Err(TwoWallVolumeHold::IncorrectFaceOrientation);
                    }
                    other+=1;
                }
            },
            Surface3::Cylinder(cyl)=>{
                if uses.len()!=4||cylinder_wall.replace(fid).is_some()
                    ||m.face_sense(fid)!=Some(FaceSense::Reversed)
                    ||(cyl.frame.origin.z-center_z).abs()>eps
                    ||distance(cyl.frame.x,Point3{x:1.,y:0.,z:0.})>1e-9
                    ||distance(cyl.frame.y,Point3{x:0.,y:0.,z:-1.})>1e-9 {
                    return Err(TwoWallVolumeHold::IncorrectFaceOrientation);
                }
                let mut arcs=Vec::<([f64;2],f64)>::new();
                let mut generators=0usize;
                for &idx in &uses {
                    let eid=m.coedges[idx].edge;
                    match m.edge_curve(eid){
                        Some(Curve3::Circle(c))=>{
                            let range=bounds(m.geometry.edges[eid.0 as usize].domain)
                                .ok_or(TwoWallVolumeHold::MissingPhysicalTrim)?;
                            if (c.radius-cyl.radius).abs()>eps
                                ||(c.frame.origin.x-cyl.frame.origin.x).abs()>eps
                                ||(c.frame.origin.z-cyl.frame.origin.z).abs()>eps
                                ||distance(c.frame.x,cyl.frame.x)>1e-9
                                ||distance(c.frame.y,cyl.frame.y)>1e-9 {
                                return Err(TwoWallVolumeHold::MissingPhysicalTrim);
                            }
                            arcs.push((range,c.frame.origin.y));
                        },
                        Some(Curve3::Line(_))=>generators+=1,
                        _=>return Err(TwoWallVolumeHold::MissingPhysicalTrim),
                    }
                }
                if arcs.len()!=2||generators!=2 {
                    return Err(TwoWallVolumeHold::MissingPhysicalTrim);
                }
                let [lo,hi]=arcs[0].0;
                if (arcs[1].0[0]-lo).abs()>eps
                    ||(arcs[1].0[1]-hi).abs()>eps{
                    return Err(TwoWallVolumeHold::MissingPhysicalTrim);
                }
                let height=(arcs[0].1-arcs[1].1).abs();
                let middle=0.5*(hi+lo);
                let n=m.oriented_face_normal(fid,middle,
                    0.5*(arcs[0].1+arcs[1].1)-cyl.frame.origin.y)
                    .map_err(|_|TwoWallVolumeHold::IncorrectFaceOrientation)?;
                if (n.z-middle.sin()).abs()>1e-9
                    ||height<=eps||hi-lo<=1e-9 {
                    return Err(TwoWallVolumeHold::IncorrectFaceOrientation);
                }
                // F=(0,0,z-cz), divergence 1. The inward cylinder
                // n_z=sin(theta), z-cz=-r*sin(theta), dS=r dtheta dy.
                let sine_sq=0.5*(hi-lo)
                    -0.25*((2.*hi).sin()-(2.*lo).sin());
                if sine_sq<=0.||!sine_sq.is_finite(){
                    return Err(TwoWallVolumeHold::InvalidNumericalBudget);
                }
                cylinder_flux-=cyl.radius*cyl.radius*height*sine_sq;
                angle_span=hi-lo;
            },
            Surface3::Sphere(_)=>return Err(TwoWallVolumeHold::UnsupportedBoundary),
        }
    }
    let total=plane_flux+cylinder_flux;
    if other!=4||z_wall.is_none()||cylinder_wall.is_none()
        ||!total.is_finite()||total<=0.{
        return Err(TwoWallVolumeHold::ShellFluxDisagreement);
    }
    Ok(TwoWallFluxWitness{
        shell,result_volume:total,z_wall_face:z_wall.unwrap(),
        cylinder_face:cylinder_wall.unwrap(),
        source_z_wall_flux:plane_flux,inward_cylinder_flux:cylinder_flux,
        angular_span:angle_span,
    })
}
impl AnalyticTwoWallVolumeOracle{
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticTwoWallPartitionPlan,tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        if tol!=self.tolerance{
            return Err(invalid("G6b source tolerance changed"));
        }
        let rebuilt=Self::build(first,second,candidate,tol)?;
        if rebuilt!=AnalyticTwoWallVolumeResult::Ready(self.clone()){
            return Err(invalid("G6b independent volume certificate differs from sources"));
        }
        Ok(())
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticTwoWallPartitionPlan,tol:GeometryTolerance,
    )->Result<AnalyticTwoWallVolumeResult,GeometryError>{
        tol.validate()?;
        let (BrepBody::Polyhedral(a),BrepBody::Cylinder(b))=(first,second)
            else{
                return Ok(AnalyticTwoWallVolumeResult::Held(
                    TwoWallVolumeHold::UnsupportedOperand));
            };
        candidate.validate_against(first,second,tol)?;
        let m=candidate.brep();
        if m.origin!=BrepOrigin::AnalyticTwoWallCylinderDifference
            ||m.shells.len()!=2||candidate.components.len()!=2
            ||candidate.can_apply||!candidate.closed_shell_proven{
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::UnsupportedCandidate));
        }
        if m.units!=first.shared()?.units||m.units!=second.shared()?.units{
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::UnsupportedOperand));
        }
        let extent=a.bbox.max.sub(a.bbox.min);
        let box_volume=extent.x*extent.y*extent.z;
        let scale=extent.length_squared().sqrt().max(b.radius)
            .max(b.height).max(1.);
        let eps=tol.length_at(scale);
        let budget=8.*eps*scale*scale;
        if !budget.is_finite()||budget<=0.||!box_volume.is_finite()
            ||box_volume<=0.||(a.mass.volume-box_volume).abs()>budget {
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::InvalidSourceMass));
        }
        let mut per_shell=Vec::new();
        let mut lower=Compensated::default();
        let mut upper=Compensated::default();
        let mut flux=Compensated::default();
        let mut seen=BTreeSet::<i8>::new();
        for part in &candidate.components{
            if !seen.insert(part.z_side)||![-1_i8,1_i8].contains(&part.z_side){
                return Ok(AnalyticTwoWallVolumeResult::Held(
                    TwoWallVolumeHold::UnsupportedCandidate));
            }
            let z=if part.z_side<0{a.bbox.min.z}else{a.bbox.max.z};
            let cylinder_center=b.frame.origin;
            let b_z=(z-cylinder_center.z).abs();
            let bracket=match source_bracket(
                b.radius,a.bbox.min.x-cylinder_center.x,
                a.bbox.max.x-cylinder_center.x,b_z,extent.y,budget/16.,
            ){
                Ok(x)=>x,Err(e)=>return Ok(AnalyticTwoWallVolumeResult::Held(e)),
            };
            let result=match shell_flux(m,part.shell,cylinder_center.z,64.*eps){
                Ok(x)=>x,Err(e)=>return Ok(AnalyticTwoWallVolumeResult::Held(e)),
            };
            if result.result_volume<bracket.lower-budget/8.
                ||result.result_volume>bracket.upper+budget/8. {
                return Ok(AnalyticTwoWallVolumeResult::Held(
                    TwoWallVolumeHold::ShellFluxDisagreement));
            }
            if part.analytic_volume<bracket.lower-budget/8.
                ||part.analytic_volume>bracket.upper+budget/8.
                ||(part.analytic_volume-result.result_volume).abs()>budget/8. {
                return Ok(AnalyticTwoWallVolumeResult::Held(
                    TwoWallVolumeHold::CandidateVolumeDisagreement));
            }
            lower.add(bracket.lower);
            upper.add(bracket.upper);
            flux.add(result.result_volume);
            per_shell.push(TwoWallVolumeComponent{
                shell:part.shell,z_side:part.z_side,
                source:bracket,boundary:result,
                candidate_volume:part.analytic_volume,
                allowed_error:budget/8.,
            });
        }
        if seen.len()!=2 {
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::UnsupportedCandidate));
        }
        let total=SourceVolumeBracket{
            lower:lower.sum,upper:upper.sum,
            integration_panels:per_shell.iter().map(|x|x.source.integration_panels).sum(),
            numerical_guard:per_shell.iter().map(|x|x.source.numerical_guard).sum(),
        };
        if !m.mass.volume.is_finite()
            ||(m.mass.volume-flux.sum).abs()>budget
            ||flux.sum<total.lower-budget||flux.sum>total.upper+budget{
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::GlobalMassDisagreement));
        }
        let removed=[box_volume-total.upper,box_volume-total.lower];
        if removed[0]<=0.||removed[1]<=removed[0]
            ||removed[1]>=box_volume
            ||removed[1]>b.mass.volume+budget{
            return Ok(AnalyticTwoWallVolumeResult::Held(
                TwoWallVolumeHold::GlobalMassDisagreement));
        }
        Ok(AnalyticTwoWallVolumeResult::Ready(Self{
            components:per_shell,source_total:total,boundary_total:flux.sum,
            stored_mass:m.mass.volume,original_box_volume:box_volume,
            removed_inside_box:removed,source_revisions:[Some(a.revision),None],
            tolerance:tol,volume_certified:true,can_apply:false,
        }))
    }
}
