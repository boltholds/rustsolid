//! G7b: independent source polygon enclosure vs actual curved B-rep flux.
//!
//! Path A bounds the source cylinder disk inside the box by inscribed and
//! tangent-circumscribed polygons. Path B integrates the divergent field
//! F=(x-cx,0,0) through the signed analytic RESULT shell coedges.
//! Neither path reads G7a's analytic circular-segment antiderivative,
//! centroid, or declared mass as an input. Numerical f64 guards are explicit;
//! this is not formal directed-rounding interval arithmetic for all CAD curves.

use crate::{
    AnalyticOneWallPlan,BrepBody,BrepModel,BrepOrigin,Curve3,
    FaceId,FaceSense,GeometryError,GeometryTolerance,ParameterRange,
    Point3,ShellId,SourceVolumeBracket,Surface3,
};
fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(message.into())
}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn bounded(range:ParameterRange)->Option<[f64;2]>{
    match range{
        ParameterRange::Bounded{start,end} if start.is_finite()
            &&end.is_finite()&&end>start=>Some([start,end]),
        _=>None,
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum OneWallVolumeHold{
    UnsupportedOperand,
    UnsupportedCandidate,
    InvalidSourceGeometry,
    UnresolvedPolygonEnclosure,
    InvalidNumericalBudget,
    MissingPhysicalCarrier,
    IncorrectBoundaryOrientation,
    SourceFluxDisagreement,
    CandidateMassDisagreement,
    SourceVolumeConservationFailed,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct OneWallFluxWitness {
    pub shell:ShellId,
    pub volume:f64,
    pub exterior_x_wall:FaceId,
    pub split_x_wall_faces:[FaceId;2],
    pub inner_cylinder_wall:FaceId,
    pub outer_flux:f64,
    pub inner_flux:f64,
    pub true_angle_span:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct AnalyticOneWallVolumeOracle{
    pub source:SourceVolumeBracket,
    pub result_flux:OneWallFluxWitness,
    pub candidate_volume:f64,
    pub original_box_volume:f64,
    pub removed_inside_box:[f64;2],
    pub source_revisions:[Option<u64>;2],
    pub tolerance:GeometryTolerance,
    pub volume_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub enum AnalyticOneWallVolumeResult{
    Ready(AnalyticOneWallVolumeOracle),
    Held(OneWallVolumeHold),
}
#[derive(Default)]
struct Kahan{sum:f64,carry:f64}
impl Kahan {
    fn add(&mut self,x:f64){
        let correction=x-self.carry;
        let next=self.sum+correction;
        self.carry=(next-self.sum)-correction;
        self.sum=next;
    }
}
/// Exact-in-real-arithmetic geometric enclosure of a circle clipped by
/// one vertical box wall. Chords lie INSIDE the disk; consecutive tangent
/// intersections lie OUTSIDE. Both close on the physical clipping chord.
/// The method uses only points, cross products and monotone refinement,
/// not the G7a circular segment mass antiderivative.
fn source_polygon_enclosure(
    radius:f64,d:f64,box_area:f64,height:f64,target:f64,
)->Result<SourceVolumeBracket,OneWallVolumeHold>{
    if ![radius,d,box_area,height,target].iter().all(|x|x.is_finite()&&*x>0.)
        ||d>=radius ||box_area<=std::f64::consts::PI*radius*radius/2. {
        return Err(OneWallVolumeHold::InvalidSourceGeometry);
    }
    let alpha=(-d/radius).acos();
    let ta=-alpha;let tb=alpha;
    let span=tb-ta;
    if !span.is_finite()||span<=std::f64::consts::PI
        ||span>=std::f64::consts::TAU{
        return Err(OneWallVolumeHold::InvalidSourceGeometry);
    }
    let point=|theta:f64,rad:f64|->(f64,f64){
        (rad*theta.cos(),-rad*theta.sin())
    };
    const MAX_N:usize=131072;
    let mut n=128usize;
    loop {
        let step=span/n as f64;
        let outer_radius=radius/(step/2.).cos();
        if !outer_radius.is_finite()||outer_radius<=radius {
            return Err(OneWallVolumeHold::UnresolvedPolygonEnclosure);
        }
        let start=point(ta,radius);
        let end=point(tb,radius);
        let mut inside=Kahan::default();
        let mut outside=Kahan::default();
        let cross=|a:(f64,f64),b:(f64,f64)|->f64{
            a.0*b.1-b.0*a.1
        };
        // Signed orientation is clockwise in (x,z).
        let mut before_inner=start;
        let mut before_outer=start;
        for i in 0..n {
            let next_angle=ta+(i+1) as f64*step;
            let to_inner=point(next_angle,radius);
            inside.add(cross(before_inner,to_inner));
            before_inner=to_inner;
            let mid=ta+(i as f64+0.5)*step;
            let to_outer=point(mid,outer_radius);
            outside.add(cross(before_outer,to_outer));
            before_outer=to_outer;
        }
        outside.add(cross(before_outer,end));
        inside.add(cross(before_inner,start));
        outside.add(cross(end,start));
        let lower_removed=0.5*inside.sum.abs();
        let upper_removed=0.5*outside.sum.abs();
        // Bound model roundoff from compensated accumulation; do not claim
        // mathematically certified IEEE-754 interval endpoints.
        let guard=32.*f64::EPSILON
            *(radius*radius*(n as f64+1.)+box_area)*height;
        if !guard.is_finite()||guard>=target/4. {
            return Err(OneWallVolumeHold::InvalidNumericalBudget);
        }
        if !lower_removed.is_finite()||!upper_removed.is_finite()
            ||lower_removed<=0.||upper_removed<lower_removed {
            return Err(OneWallVolumeHold::UnresolvedPolygonEnclosure);
        }
        let source_lower=(box_area-upper_removed)*height-guard;
        let source_upper=(box_area-lower_removed)*height+guard;
        if source_lower<=0.||!source_lower.is_finite()
            ||!source_upper.is_finite(){
            return Err(OneWallVolumeHold::InvalidSourceGeometry);
        }
        if source_upper-source_lower<=target{
            return Ok(SourceVolumeBracket {
                lower:source_lower,upper:source_upper,
                integration_panels:n,numerical_guard:guard,
            });
        }
        if n>=MAX_N{
            return Err(OneWallVolumeHold::UnresolvedPolygonEnclosure);
        }
        n*=2;
    }
}
fn coedge_loop(m:&BrepModel,face:FaceId)->Result<Vec<usize>,OneWallVolumeHold>{
    let selected=m.faces.get(face.0 as usize)
        .filter(|f|f.id==face&&f.loops.len()==1)
        .ok_or(OneWallVolumeHold::MissingPhysicalCarrier)?;
    let loop_id=selected.loops[0];
    let lp=m.loops.get(loop_id.0 as usize)
        .ok_or(OneWallVolumeHold::MissingPhysicalCarrier)?;
    let mut uses=Vec::new();
    let mut cursor=lp.first_coedge;
    loop {
        let co=m.coedges.get(cursor.0 as usize)
            .filter(|c|c.face==face&&c.loop_id==loop_id)
            .ok_or(OneWallVolumeHold::MissingPhysicalCarrier)?;
        if uses.contains(&(cursor.0 as usize)){
            return Err(OneWallVolumeHold::MissingPhysicalCarrier);
        }
        uses.push(cursor.0 as usize);
        cursor=co.next;
        if cursor==lp.first_coedge{break;}
        if uses.len()>16{return Err(OneWallVolumeHold::MissingPhysicalCarrier);}
    }
    if uses.len()<3{return Err(OneWallVolumeHold::MissingPhysicalCarrier);}
    Ok(uses)
}
fn oriented_flux(
    m:&BrepModel,shell:ShellId,cx:f64,eps:f64,
)->Result<OneWallFluxWitness,OneWallVolumeHold>{
    let s=m.shells.get(shell.0 as usize)
        .filter(|sh|sh.id==shell&&sh.faces.len()==8&&sh.closed)
        .ok_or(OneWallVolumeHold::UnsupportedCandidate)?;
    let mut right_wall=None;
    let mut left_walls=Vec::<FaceId>::new();
    let mut cylinder=None;
    let mut other_planes=0usize;
    let mut outer_flux=Kahan::default();
    let mut inner_flux=0.;
    let mut angular_span=0.;
    for &fid in &s.faces{
        let support=m.face_surface(fid)
            .ok_or(OneWallVolumeHold::MissingPhysicalCarrier)?;
        let uses=coedge_loop(m,fid)?;
        let points=uses.iter().map(|&idx|{
            let co=&m.coedges[idx];
            let edge=&m.edges[co.edge.0 as usize];
            let vertex=if co.reversed{edge.end}else{edge.start};
            m.vertices[vertex.0 as usize].position
        }).collect::<Vec<_>>();
        match support {
            Surface3::Plane(_)=>{
                let normal=m.oriented_face_normal(fid,0.,0.)
                    .map_err(|_|OneWallVolumeHold::IncorrectBoundaryOrientation)?;
                if normal.x.abs()>1.-1e-9 {
                    if uses.len()!=4||normal.y.abs()>1e-9||normal.z.abs()>1e-9
                        ||points.iter().any(|p|(p.x-points[0].x).abs()>eps) {
                        return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
                    }
                    if normal.x>0. {
                        if right_wall.replace(fid).is_some() {
                            return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
                        }
                    }else{left_walls.push(fid);}
                    let twice_area=points.iter().enumerate().map(|(i,p)|{
                        let next=points[(i+1)%points.len()];
                        p.y*next.z-next.y*p.z
                    }).sum::<f64>();
                    let area=0.5*twice_area.abs();
                    if !area.is_finite()||area<=eps*eps {
                        return Err(OneWallVolumeHold::MissingPhysicalCarrier);
                    }
                    outer_flux.add((points[0].x-cx)*normal.x*area);
                }else{
                    if normal.x.abs()>1e-9{
                        return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
                    }
                    other_planes+=1;
                }
            },
            Surface3::Cylinder(cyl)=>{
                if cylinder.replace(fid).is_some()
                    ||m.face_sense(fid)!=Some(FaceSense::Reversed)
                    ||uses.len()!=4
                    ||(cyl.frame.origin.x-cx).abs()>eps
                    ||distance(cyl.frame.x,Point3{x:1.,y:0.,z:0.})>1e-9
                    ||distance(cyl.frame.y,Point3{x:0.,y:0.,z:-1.})>1e-9{
                    return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
                }
                let mut arcs=Vec::<([f64;2],f64)>::new();
                let mut generators=0usize;
                for &idx in &uses{
                    let edgeid=m.coedges[idx].edge;
                    match m.edge_curve(edgeid){
                        Some(Curve3::Circle(c))=>{
                            let domain=bounded(m.geometry.edges[edgeid.0 as usize].domain)
                                .ok_or(OneWallVolumeHold::MissingPhysicalCarrier)?;
                            if (c.radius-cyl.radius).abs()>eps
                                ||(c.frame.origin.x-cyl.frame.origin.x).abs()>eps
                                ||(c.frame.origin.z-cyl.frame.origin.z).abs()>eps
                                ||distance(c.frame.x,cyl.frame.x)>1e-9
                                ||distance(c.frame.y,cyl.frame.y)>1e-9 {
                                return Err(OneWallVolumeHold::MissingPhysicalCarrier);
                            }
                            arcs.push((domain,c.frame.origin.y));
                        },
                        Some(Curve3::Line(_))=>generators+=1,
                        _=>return Err(OneWallVolumeHold::MissingPhysicalCarrier),
                    }
                }
                if arcs.len()!=2||generators!=2{
                    return Err(OneWallVolumeHold::MissingPhysicalCarrier);
                }
                let [ta,tb]=arcs[0].0;
                if (arcs[1].0[0]-ta).abs()>eps
                    ||(arcs[1].0[1]-tb).abs()>eps
                    ||tb-ta<=std::f64::consts::PI
                    ||tb-ta>=std::f64::consts::TAU {
                    return Err(OneWallVolumeHold::MissingPhysicalCarrier);
                }
                let height=(arcs[0].1-arcs[1].1).abs();
                let mid=(ta+tb)/2.;
                let n=m.oriented_face_normal(fid,mid,
                    0.5*(arcs[0].1+arcs[1].1)-cyl.frame.origin.y)
                    .map_err(|_|OneWallVolumeHold::IncorrectBoundaryOrientation)?;
                if height<=eps||(n.x+mid.cos()).abs()>1e-9 {
                    return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
                }
                let cos2=0.5*(tb-ta)+0.25*((2.*tb).sin()-(2.*ta).sin());
                if cos2<=0.||!cos2.is_finite(){
                    return Err(OneWallVolumeHold::InvalidNumericalBudget);
                }
                inner_flux=-cyl.radius*cyl.radius*height*cos2;
                angular_span=tb-ta;
            },
            Surface3::Sphere(_)=>return Err(OneWallVolumeHold::UnsupportedCandidate),
        }
    }
    if right_wall.is_none()||left_walls.len()!=2||other_planes!=4
        ||cylinder.is_none(){
        return Err(OneWallVolumeHold::IncorrectBoundaryOrientation);
    }
    let volume=outer_flux.sum+inner_flux;
    if !volume.is_finite()||volume<=0. {
        return Err(OneWallVolumeHold::SourceFluxDisagreement);
    }
    Ok(OneWallFluxWitness{
        shell,volume,exterior_x_wall:right_wall.unwrap(),
        split_x_wall_faces:[left_walls[0],left_walls[1]],
        inner_cylinder_wall:cylinder.unwrap(),
        outer_flux:outer_flux.sum,inner_flux,
        true_angle_span:angular_span,
    })
}
impl AnalyticOneWallVolumeOracle {
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticOneWallPlan,tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        if self.tolerance!=tol{
            return Err(invalid("G7b source tolerance changed"));
        }
        if Self::build(first,second,candidate,tol)?
            !=AnalyticOneWallVolumeResult::Ready(self.clone()){
            return Err(invalid("G7b independent volume proof differs from sources"));
        }
        Ok(())
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        candidate:&AnalyticOneWallPlan,tol:GeometryTolerance,
    )->Result<AnalyticOneWallVolumeResult,GeometryError>{
        tol.validate()?;
        let (BrepBody::Polyhedral(a),BrepBody::Cylinder(b))=(first,second)
            else{return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::UnsupportedOperand));};
        candidate.validate_against(first,second,tol)?;
        let m=candidate.brep();
        if m.origin!=BrepOrigin::AnalyticOneWallCylinderDifference
            ||m.shells.len()!=1||candidate.can_apply||!candidate.closed_shell_proven{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::UnsupportedCandidate));
        }
        if m.units!=first.shared()?.units || m.units!=second.shared()?.units{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::UnsupportedOperand));
        }
        let sizes=a.bbox.max.sub(a.bbox.min);
        let extent=sizes.length_squared().sqrt().max(b.radius)
            .max(b.height).max(1.);
        let eps=tol.length_at(extent);
        let budget=8.*eps*extent*extent;
        let box_volume=sizes.x*sizes.y*sizes.z;
        if !budget.is_finite()||budget<=0.||!box_volume.is_finite()
            ||box_volume<=0.||(a.mass.volume-box_volume).abs()>budget {
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::InvalidSourceGeometry));
        }
        let d=b.frame.origin.x-a.bbox.min.x;
        let area=sizes.x*sizes.z;
        let bracket=match source_polygon_enclosure(
            b.radius,d,area,sizes.y,budget/16.,
        ){
            Ok(x)=>x,Err(err)=>return Ok(AnalyticOneWallVolumeResult::Held(err)),
        };
        let flux=match oriented_flux(m,ShellId(0),b.frame.origin.x,64.*eps){
            Ok(x)=>x,Err(err)=>return Ok(AnalyticOneWallVolumeResult::Held(err)),
        };
        if flux.volume<bracket.lower-budget/8.
            ||flux.volume>bracket.upper+budget/8.{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::SourceFluxDisagreement));
        }
        if candidate.component.exact_volume<bracket.lower-budget/8.
            ||candidate.component.exact_volume>bracket.upper+budget/8.
            ||(candidate.component.exact_volume-flux.volume).abs()>budget/8.{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::CandidateMassDisagreement));
        }
        if !m.mass.volume.is_finite()
            ||(m.mass.volume-flux.volume).abs()>budget
            ||m.mass.volume<bracket.lower-budget
            ||m.mass.volume>bracket.upper+budget{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::CandidateMassDisagreement));
        }
        let removed=[
            box_volume-bracket.upper,
            box_volume-bracket.lower,
        ];
        if removed[0]<=0.||removed[1]<=removed[0]
            ||removed[1]>=box_volume
            ||removed[1]>b.mass.volume+budget{
            return Ok(AnalyticOneWallVolumeResult::Held(
                OneWallVolumeHold::SourceVolumeConservationFailed));
        }
        Ok(AnalyticOneWallVolumeResult::Ready(Self{
            source:bracket,result_flux:flux,
            candidate_volume:m.mass.volume,
            original_box_volume:box_volume,
            removed_inside_box:removed,
            source_revisions:[Some(a.revision),None],
            tolerance:tol,volume_certified:true,can_apply:false,
        }))
    }
}
