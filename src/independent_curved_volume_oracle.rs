//! G4c: independent source-volume integration and oriented B-rep flux.
//! A bounded, exact-geometry four-corner Box minus Cylinder certificate.
//! No mesh, candidate area antiderivative, or candidate mass enters either
//! oracle. Bounds follow convexity in real arithmetic with explicit f64 guard;
//! they are not arbitrary-surface directed-rounding interval arithmetic.
//! No G5 commit is enabled.

use crate::{AnalyticFacePatchPartitionPlan,BrepBody,BrepModel,BrepOrigin,
    Curve3,EdgeId,FaceId,FaceSense,GeometryError,GeometryTolerance,
    ParameterRange,Point3,ShellId,Surface3};

fn invalid(x:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(x.into())}
fn distance(a:Point3,b:Point3)->f64{let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)}
fn interval(x:ParameterRange)->Option<[f64;2]>{
    match x{ParameterRange::Bounded{start,end} if start.is_finite()
        &&end.is_finite()&&end>start=>Some([start,end]),_=>None}
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum IndependentVolumeHold {
    UnsupportedOperand,InvalidSourceGeometry,InconsistentUnits,
    IncorrectCandidateOrigin,UnexpectedShellStructure,
    MissingPhysicalFaceCarrier,IncompatibleCylinderChart,
    InconsistentBoundaryNormals,NumericalBudgetExceeded,
    SourceQuadratureUnresolved,ShellFluxMismatch,
    ComponentVolumeMismatch,ResultMassMismatch,
    SourceMassConservationFailed,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct SourceVolumeBracket {
    pub lower:f64,pub upper:f64,
    pub integration_panels:usize,pub numerical_guard:f64,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ShellFluxWitness {
    pub result_volume:f64,
    pub original_x_wall_flux:f64,
    pub inward_cylinder_flux:f64,
    pub cylinder_angle_span:f64,
    pub x_wall:FaceId,pub cylinder_wall:FaceId,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct IndependentVolumeComponent {
    pub shell:ShellId,pub corner_sides:[i8;2],
    pub source_volume:SourceVolumeBracket,
    pub result_flux:ShellFluxWitness,
    pub candidate_component_volume:f64,
    pub permitted_absolute_error:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct IndependentVolumeOraclePlan {
    pub components:Vec<IndependentVolumeComponent>,
    pub source_total:SourceVolumeBracket,
    pub result_flux_total:f64,
    pub candidate_mass_volume:f64,
    pub original_box_volume:f64,
    /// Only intersection WITHIN the box; not the entire cutter's volume.
    pub removed_volume_bounds:[f64;2],
    /// Cylinder has no revision counter; candidate revalidation checks its full geometry.
    pub source_revisions:[Option<u64>;2],
    pub construction_tolerance:GeometryTolerance,
    pub volume_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum IndependentVolumeOracleResult{
    Ready(IndependentVolumeOraclePlan),Held(IndependentVolumeHold),
}
#[derive(Default)]
struct Kahan{sum:f64,c:f64}
impl Kahan{fn add(&mut self,x:f64){
    let y=x-self.c;let t=self.sum+y;
    self.c=(t-self.sum)-y;self.sum=t;
}}
#[derive(Clone,Copy)]
struct Slice{lo:f64,hi:f64,depth:usize}

/// Independent source-only numerical area integral of
/// w(u)=b-sqrt(r^2-u^2), u in [sqrt(r^2-b^2),a].
/// w''(u)=r^2/(r^2-u^2)^(3/2)>0: midpoint is an area LOWER bound,
/// trapezoid an UPPER bound. Adaptive partition gives a finite bracket.
fn source_quadrature(
    r:f64,a:f64,b:f64,height:f64,target:f64,
)->Result<SourceVolumeBracket,IndependentVolumeHold>{
    if ![r,a,b,height,target].iter().all(|x|x.is_finite()&&*x>0.0)
       ||a>=r||b>=r||a.hypot(b)<=r{
        return Err(IndependentVolumeHold::InvalidSourceGeometry);
    }
    let u0=((r-b)*(r+b)).sqrt();
    let min_y=((r-a)*(r+a)).sqrt();
    let span=a-u0;
    if !u0.is_finite()||span<=0.0||min_y<=0.0{
        return Err(IndependentVolumeHold::InvalidSourceGeometry);
    }
    let guard=256.0*f64::EPSILON
        *(r+b+r*r/min_y)*span*height;
    if !guard.is_finite()||guard>=target*0.25{
        return Err(IndependentVolumeHold::NumericalBudgetExceeded);
    }
    let width=|u:f64|->Result<f64,IndependentVolumeHold>{
        if !u.is_finite()||u<u0||u>a {
            return Err(IndependentVolumeHold::InvalidSourceGeometry);
        }
        let arg=(r-u)*(r+u);
        if arg<0.0||!arg.is_finite(){
            return Err(IndependentVolumeHold::InvalidSourceGeometry);
        }
        let w=b-arg.sqrt();
        if !w.is_finite()||w< -64.0*f64::EPSILON*(r+b){
            return Err(IndependentVolumeHold::InvalidSourceGeometry);
        }
        Ok(w.max(0.0))
    };
    let mut lower=Kahan::default();let mut upper=Kahan::default();
    let mut stack=vec![Slice{lo:u0,hi:a,depth:0}];
    let mut panels=0usize;
    const LIMIT:usize=131072;
    while let Some(s)=stack.pop(){
        let mid=0.5*s.lo+0.5*s.hi;
        let dx=s.hi-s.lo;
        if dx<=0.0||mid<=s.lo||mid>=s.hi{
            return Err(IndependentVolumeHold::SourceQuadratureUnresolved);
        }
        let w0=width(s.lo)?;let wm=width(mid)?;let w1=width(s.hi)?;
        let lo=dx*wm;let hi=0.5*dx*(w0+w1);
        if !lo.is_finite()||!hi.is_finite()||hi+guard/height<lo{
            return Err(IndependentVolumeHold::SourceQuadratureUnresolved);
        }
        let local=(target-guard)*0.5*dx/span/height;
        if hi-lo<=local{
            lower.add(lo);upper.add(hi.max(lo));panels+=1;
            if panels>LIMIT {
                return Err(IndependentVolumeHold::SourceQuadratureUnresolved);
            }
        }else{
            if s.depth>=28||panels+stack.len()>=LIMIT{
                return Err(IndependentVolumeHold::SourceQuadratureUnresolved);
            }
            stack.push(Slice{lo:mid,hi:s.hi,depth:s.depth+1});
            stack.push(Slice{lo:s.lo,hi:mid,depth:s.depth+1});
        }
    }
    let low=(lower.sum*height-guard).max(0.0);
    let high=upper.sum*height+guard;
    if !low.is_finite()||!high.is_finite()||high<=low
        ||high-low>target{
        return Err(IndependentVolumeHold::SourceQuadratureUnresolved);
    }
    Ok(SourceVolumeBracket{lower:low,upper:high,
        integration_panels:panels,numerical_guard:guard})
}
fn face_uses(m:&BrepModel,id:FaceId)->Result<Vec<usize>,IndependentVolumeHold>{
    let face=m.faces.get(id.0 as usize).filter(|f|f.id==id&&f.loops.len()==1)
        .ok_or(IndependentVolumeHold::UnexpectedShellStructure)?;
    let loop_id=face.loops[0];
    let lp=m.loops.get(loop_id.0 as usize)
        .ok_or(IndependentVolumeHold::MissingPhysicalFaceCarrier)?;
    let mut ids=Vec::new();let mut at=lp.first_coedge;
    loop{
        let co=m.coedges.get(at.0 as usize)
            .filter(|co|co.face==id&&co.loop_id==loop_id)
            .ok_or(IndependentVolumeHold::MissingPhysicalFaceCarrier)?;
        if ids.contains(&(at.0 as usize)){
            return Err(IndependentVolumeHold::MissingPhysicalFaceCarrier);
        }
        ids.push(at.0 as usize);at=co.next;
        if at==lp.first_coedge{break;}
        if ids.len()>16{
            return Err(IndependentVolumeHold::UnexpectedShellStructure);
        }
    }
    if ids.len()<3{return Err(IndependentVolumeHold::UnexpectedShellStructure);}
    Ok(ids)
}
fn shell_flux(m:&BrepModel,shell:ShellId,cx:f64,eps:f64)
    ->Result<ShellFluxWitness,IndependentVolumeHold>{
    let owned=m.shells.get(shell.0 as usize)
        .filter(|s|s.id==shell&&s.closed&&s.faces.len()==5)
        .ok_or(IndependentVolumeHold::UnexpectedShellStructure)?;
    let mut x_wall=None;let mut cylinder_wall=None;
    let mut other_planes=0usize;
    let mut plane_flux=0.0;let mut cyl_flux=0.0;let mut angle_span=0.0;
    for &fid in &owned.faces{
        let surface=m.face_surface(fid)
            .ok_or(IndependentVolumeHold::MissingPhysicalFaceCarrier)?;
        let uses=face_uses(m,fid)?;
        let points=uses.iter().map(|&idx|{
            let co=&m.coedges[idx];let edge=&m.edges[co.edge.0 as usize];
            let v=if co.reversed{edge.end}else{edge.start};
            m.vertices[v.0 as usize].position
        }).collect::<Vec<_>>();
        match surface{
            Surface3::Plane(_)=>{
                let n=m.oriented_face_normal(fid,0.,0.)
                    .map_err(|_|IndependentVolumeHold::InconsistentBoundaryNormals)?;
                if n.x.abs()>1.-1e-9 {
                    if uses.len()!=4||x_wall.replace(fid).is_some()
                        ||n.y.abs()>1e-9||n.z.abs()>1e-9
                        ||points.iter().any(|p|(p.x-points[0].x).abs()>eps){
                        return Err(IndependentVolumeHold::InconsistentBoundaryNormals);
                    }
                    let double_area=points.iter().enumerate().map(|(i,p)|{
                        let q=points[(i+1)%points.len()];
                        p.y*q.z-q.y*p.z
                    }).sum::<f64>();
                    let area=0.5*double_area.abs();
                    if !area.is_finite()||area<=eps*eps {
                        return Err(IndependentVolumeHold::MissingPhysicalFaceCarrier);
                    }
                    plane_flux+=(points[0].x-cx)*n.x*area;
                }else{
                    if n.x.abs()>1e-9 {
                        return Err(IndependentVolumeHold::InconsistentBoundaryNormals);
                    }
                    other_planes+=1;
                }
            },
            Surface3::Cylinder(cyl)=>{
                if cylinder_wall.replace(fid).is_some()||uses.len()!=4
                    ||m.face_sense(fid)!=Some(FaceSense::Reversed)
                    ||(cyl.frame.origin.x-cx).abs()>eps
                    ||distance(cyl.frame.x,Point3{x:1.,y:0.,z:0.})>1e-9
                    ||distance(cyl.frame.normal,Point3{x:0.,y:1.,z:0.})>1e-9{
                    return Err(IndependentVolumeHold::IncompatibleCylinderChart);
                }
                let mut circles=Vec::<([f64;2],f64)>::new();
                let mut straight=0usize;
                for &id in &uses {
                    let eid=m.coedges[id].edge;
                    match m.edge_curve(eid) {
                        Some(Curve3::Circle(c))=>{
                            let dom=interval(m.geometry.edges[eid.0 as usize].domain)
                                .ok_or(IndependentVolumeHold::MissingPhysicalFaceCarrier)?;
                            if (c.radius-cyl.radius).abs()>eps
                                ||(c.frame.origin.x-cyl.frame.origin.x).abs()>eps
                                ||(c.frame.origin.z-cyl.frame.origin.z).abs()>eps
                                ||distance(c.frame.x,cyl.frame.x)>1e-9
                                ||distance(c.frame.y,cyl.frame.y)>1e-9 {
                                return Err(IndependentVolumeHold::IncompatibleCylinderChart);
                            }
                            circles.push((dom,c.frame.origin.y));
                        },
                        Some(Curve3::Line(_))=>straight+=1,
                        _=>return Err(IndependentVolumeHold::MissingPhysicalFaceCarrier),
                    }
                }
                if circles.len()!=2||straight!=2 {
                    return Err(IndependentVolumeHold::MissingPhysicalFaceCarrier);
                }
                let [lo,hi]=circles[0].0;
                if (lo-circles[1].0[0]).abs()>eps
                    ||(hi-circles[1].0[1]).abs()>eps {
                    return Err(IndependentVolumeHold::IncompatibleCylinderChart);
                }
                let height=(circles[0].1-circles[1].1).abs();
                if height<=eps||hi-lo<=1e-9 {
                    return Err(IndependentVolumeHold::MissingPhysicalFaceCarrier);
                }
                let middle=0.5*(lo+hi);
                let normal=m.oriented_face_normal(fid,middle,
                    0.5*(circles[0].1+circles[1].1)-cyl.frame.origin.y)
                    .map_err(|_|IndependentVolumeHold::InconsistentBoundaryNormals)?;
                if (normal.x+middle.cos()).abs()>1e-9 {
                    return Err(IndependentVolumeHold::InconsistentBoundaryNormals);
                }
                // Divergence theorem with F=(x-cx,0,0), div F=1.
                // Actual inward oriented cylindrical B-rep contributes
                // -radius^2 * height * integral(cos^2(u)du).
                let cosine_integral=0.5*(hi-lo)
                    +0.25*((2.0*hi).sin()-(2.0*lo).sin());
                if !cosine_integral.is_finite()||cosine_integral<=0.0 {
                    return Err(IndependentVolumeHold::NumericalBudgetExceeded);
                }
                cyl_flux-=cyl.radius*cyl.radius*height*cosine_integral;
                angle_span=hi-lo;
            },
            Surface3::Sphere(_)=>return Err(IndependentVolumeHold::UnsupportedOperand),
        }
    }
    let result_volume=plane_flux+cyl_flux;
    if other_planes!=3||x_wall.is_none()||cylinder_wall.is_none()
        ||!result_volume.is_finite()||result_volume<=0.0{
        return Err(IndependentVolumeHold::ShellFluxMismatch);
    }
    Ok(ShellFluxWitness{result_volume,original_x_wall_flux:plane_flux,
        inward_cylinder_flux:cyl_flux,cylinder_angle_span:angle_span,
        x_wall:x_wall.unwrap(),cylinder_wall:cylinder_wall.unwrap()})
}
impl IndependentVolumeOraclePlan{
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticFacePatchPartitionPlan,tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        if self.construction_tolerance!=tol{
            return Err(invalid("G4c tolerance changed"));
        }
        if Self::build(first,second,candidate,tol)?
            !=IndependentVolumeOracleResult::Ready(self.clone()){
            return Err(invalid("G4c independently reconstructed evidence differs"));
        }
        Ok(())
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticFacePatchPartitionPlan,tol:GeometryTolerance,
    )->Result<IndependentVolumeOracleResult,GeometryError>{
        tol.validate()?;
        let (BrepBody::Polyhedral(a),BrepBody::Cylinder(b))=(first,second)
            else{return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::UnsupportedOperand));};
        candidate.validate_against(first,second,tol)?;
        let m=candidate.brep();
        if m.origin!=BrepOrigin::AnalyticPartialBoxCylinderDifference
            ||candidate.can_apply||!candidate.closed_shell_proven {
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::IncorrectCandidateOrigin));
        }
        if m.units!=first.shared()?.units ||m.units!=second.shared()?.units{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::InconsistentUnits));
        }
        if m.shells.len()!=4||candidate.corners.len()!=4{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::UnexpectedShellStructure));
        }
        let bounds=a.bbox;
        let sizes=bounds.max.sub(bounds.min);
        let h=sizes.y;
        let box_volume=sizes.x*sizes.y*sizes.z;
        let center=b.frame.origin;
        let scale=sizes.length_squared().sqrt().max(b.radius).max(b.height).max(1.0);
        let precision=tol.length_at(scale);
        let budget=8.0*precision*scale*scale;
        if !box_volume.is_finite()||box_volume<=0.0
            ||!budget.is_finite()||budget<=0.0{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::InvalidSourceGeometry));
        }
        if (a.mass.volume-box_volume).abs()>budget{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::SourceMassConservationFailed));
        }
        let mut components=Vec::new();
        let mut low=Kahan::default();let mut high=Kahan::default();
        let mut flux_total=Kahan::default();let mut declared_total=Kahan::default();
        let mut seen=std::collections::BTreeSet::new();
        for corner in &candidate.corners {
            if !seen.insert((corner.x_side,corner.z_side))
                ||![-1_i8,1_i8].contains(&corner.x_side)
                ||![-1_i8,1_i8].contains(&corner.z_side)
                ||corner.shell.0 as usize>=m.shells.len(){
                return Ok(IndependentVolumeOracleResult::Held(
                    IndependentVolumeHold::UnexpectedShellStructure));
            }
            let x=if corner.x_side<0{bounds.min.x}else{bounds.max.x};
            let z=if corner.z_side<0{bounds.min.z}else{bounds.max.z};
            let a_local=(x-center.x).abs();let b_local=(z-center.z).abs();
            let source=match source_quadrature(b.radius,a_local,b_local,h,budget/16.0){
                Ok(v)=>v,Err(reason)=>
                    return Ok(IndependentVolumeOracleResult::Held(reason)),
            };
            let result=match shell_flux(m,corner.shell,center.x,precision*64.0){
                Ok(v)=>v,Err(reason)=>
                    return Ok(IndependentVolumeOracleResult::Held(reason)),
            };
            if result.result_volume<source.lower-budget/8.0
                ||result.result_volume>source.upper+budget/8.0{
                return Ok(IndependentVolumeOracleResult::Held(
                    IndependentVolumeHold::ShellFluxMismatch));
            }
            if corner.exact_volume<source.lower-budget/8.0
                ||corner.exact_volume>source.upper+budget/8.0
                ||(corner.exact_volume-result.result_volume).abs()>budget/8.0{
                return Ok(IndependentVolumeOracleResult::Held(
                    IndependentVolumeHold::ComponentVolumeMismatch));
            }
            low.add(source.lower);high.add(source.upper);
            flux_total.add(result.result_volume);declared_total.add(corner.exact_volume);
            components.push(IndependentVolumeComponent{
                shell:corner.shell,corner_sides:[corner.x_side,corner.z_side],
                source_volume:source,result_flux:result,
                candidate_component_volume:corner.exact_volume,
                permitted_absolute_error:budget/8.0,
            });
        }
        if seen.len()!=4{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::UnexpectedShellStructure));
        }
        let total=SourceVolumeBracket{
            lower:low.sum,upper:high.sum,
            integration_panels:components.iter().map(|c|c.source_volume.integration_panels).sum(),
            numerical_guard:components.iter().map(|c|c.source_volume.numerical_guard).sum(),
        };
        if (declared_total.sum-flux_total.sum).abs()>budget
            ||flux_total.sum<total.lower-budget
            ||flux_total.sum>total.upper+budget{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::ComponentVolumeMismatch));
        }
        if m.mass.volume<total.lower-budget
            ||m.mass.volume>total.upper+budget
            ||(m.mass.volume-flux_total.sum).abs()>budget{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::ResultMassMismatch));
        }
        let removed=[box_volume-total.upper,box_volume-total.lower];
        if removed[0]<=0.0||removed[1]<=removed[0]
            ||removed[1]>=box_volume
            ||removed[1]>b.mass.volume+budget{
            return Ok(IndependentVolumeOracleResult::Held(
                IndependentVolumeHold::SourceMassConservationFailed));
        }
        Ok(IndependentVolumeOracleResult::Ready(Self{
            components,source_total:total,result_flux_total:flux_total.sum,
            candidate_mass_volume:m.mass.volume,original_box_volume:box_volume,
            removed_volume_bounds:removed,
            source_revisions:[Some(a.revision),None],
            construction_tolerance:tol,volume_certified:true,can_apply:false,
        }))
    }
}
