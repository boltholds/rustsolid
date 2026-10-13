//! Source-certified 3D halfedge stitching PLAN for one connected three-disk
//! chain across two distinct blind-depth events.
//!
//! Unlike v0.33 this is not component grafting: the middle cylinder joins
//! the through cutter at the deeper floor and the third cylinder joins the
//! already connected contour at a shallower floor. Global pair-junction IDs
//! and circular-arc atoms are shared across ALL three height slabs.
//!
//! This builds physical horizontal/vertical edge-use incidence, oriented
//! source-supported wall patches, planar floors and cap contour loops.
//! It deliberately does NOT publish a BrepModel or a Boolean commit.
use std::collections::{BTreeMap,BTreeSet};
use std::f64::consts::{PI,TAU};
use crate::{BlindHoleParameters,GeometryError,GeometryTolerance,HoleParameters,
    MultiCutPlanState,MultiCutRegionPlanner,Point3,Solid};

fn invalid(msg:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(msg.into())}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TripleStitchHold {
    UnresolvedSourceArrangement,
    UnsupportedCircleContact,
    NearTangency,
    AmbiguousArcExposure,
    OpenPlanarContour,
    NonmanifoldPhysicalEdge,
    InconsistentAnalyticArea,
    InconsistentEulerCharacteristic,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TripleStitchState { Ready, Held(TripleStitchHold) }

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TripleChainParameters {
    pub through:HoleParameters,
    pub deeper_blind:BlindHoleParameters,
    pub shallower_blind:BlindHoleParameters,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleJunction {
    pub id:usize,
    pub source_pair:[usize;2],
    pub position_xz:[f64;2],
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleArc {
    pub id:usize,
    pub source:usize,
    pub start_junction:usize,
    pub end_junction:usize,
    pub start_angle:f64,
    pub end_angle:f64,
    pub arc_length:f64,
    pub planar_area_contribution:f64,
    pub visible_in_slab:[bool;3],
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct DirectedTripleArc {
    pub arc:usize,
    pub reversed:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleContour {
    pub parts:Vec<DirectedTripleArc>,
    pub signed_area:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleSlab {
    pub y_lower:f64,
    pub y_upper:f64,
    pub active_sources:Vec<usize>,
    pub exposed_contour:TripleContour,
    pub union_area:f64,
    pub removed_volume:f64,
}
/// Interface between lower and upper void slabs, oriented CCW in world XZ
/// for planar area integration. The physical +Y floor face has REVERSED UV
/// edge directions relative to this contour.
#[derive(Debug,Clone,PartialEq)]
pub struct TripleFloorInterface {
    pub y:f64,
    pub contour:TripleContour,
    pub floor_area:f64,
    pub incoming_lower_arcs:Vec<usize>,
    pub born_upper_arcs:Vec<usize>,
}
/// One genuine cylindrical surface patch to be constructed from one source
/// angular interval, with V pointing DOWNWARD from the source box top plane.
/// Its material-facing normal must be -radial, not +radial.
#[derive(Debug,Clone,PartialEq)]
pub struct TripleWallPatch {
    pub source:usize,
    pub arc:usize,
    pub slab:usize,
    pub y_lower:f64,
    pub y_upper:f64,
    pub uv_bounds:[f64;4], // [u0,u1,v0,v1]; v = topY-y
    pub corners:[Point3;4], // lower start/end, upper start/end
    pub analytic_area:f64,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord)]
pub enum TriplePhysicalEdgeKey {
    CircleArc {arc:usize,level:usize},
    VerticalTrace {junction:usize,slab:usize},
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleEdgePair {
    pub key:TriplePhysicalEdgeKey,
    pub opposite_use_signs:[i8;2],
}
#[derive(Debug,Clone,PartialEq)]
pub struct TripleStitchPlan {
    pub state:TripleStitchState,
    pub input:TripleChainParameters,
    pub source_revision:u64,
    pub source_bounds:crate::BoundingBox,
    pub source_volume:f64,
    pub tolerance:GeometryTolerance,
    pub length_budget:f64,
    pub levels:Vec<f64>,
    pub junctions:Vec<TripleJunction>,
    pub atomic_arcs:Vec<TripleArc>,
    pub slabs:Vec<TripleSlab>,
    pub floors:Vec<TripleFloorInterface>,
    pub bottom_cap:Option<TripleContour>,
    pub top_cap:Option<TripleContour>,
    pub walls:Vec<TripleWallPatch>,
    pub shared_physical_edges:Vec<TripleEdgePair>,
    pub predicted_shell_vertices:Option<usize>,
    pub predicted_shell_edges:Option<usize>,
    pub predicted_shell_faces:Option<usize>,
    pub predicted_euler:Option<isize>,
    pub predicted_genus:Option<usize>,
    pub analytic_surface_area:Option<f64>,
    pub predicted_remaining_volume:Option<f64>,
}
impl TripleStitchPlan{
    /// A valid 3D halfedge plan is not yet a topology-validated BrepModel.
    pub fn can_apply(&self)->bool{false}
    pub fn validate_against(&self,source:&Solid,input:TripleChainParameters)
        ->Result<(),GeometryError>{
        let expected=ConnectedTripleStitchPlanner::plan(source,input)?;
        if *self!=expected{return Err(invalid("stitch plan differs from source-rebuilt witnesses"));}
        Ok(())
    }
}
#[derive(Debug,Clone,Copy,Default)]
pub struct ConnectedTripleStitchPlanner;

#[derive(Debug,Clone,Copy)]
struct Disk{cx:f64,cz:f64,r:f64}
impl Disk{
    fn point(self,angle:f64,y:f64)->Point3{
        Point3{x:self.cx+self.r*angle.cos(),y,
            z:self.cz+self.r*angle.sin()}
    }
    fn gap(self,p:Point3)->f64{
        (p.x-self.cx).hypot(p.z-self.cz)-self.r
    }
    fn area(self)->f64{PI*self.r*self.r}
}
#[derive(Debug,Clone,Copy)]
struct Event {angle:f64,junction:usize}
fn angles(d:Disk,p:[f64;2])->f64{
    (p[1]-d.cz).atan2(p[0]-d.cx).rem_euclid(TAU)
}
fn safe_near(a:f64,b:f64,eps:f64)->bool{(a-b).abs()<=eps}
fn kahan(values:impl IntoIterator<Item=f64>)->f64{
    let (mut total,mut correction)=(0.0,0.0);
    for x in values{
        let y=x-correction;
        let t=total+y;
        correction=(t-total)-y;
        total=t;
    }
    total
}
fn green(d:Disk,u0:f64,u1:f64,origin:(f64,f64))->f64{
    let x=d.cx-origin.0;
    let z=d.cz-origin.1;
    0.5*(d.r*d.r*(u1-u0)+
        d.r*x*(u1.sin()-u0.sin())+
        d.r*z*(u0.cos()-u1.cos()))
}
fn pair_intersections(a:Disk,b:Disk,y:f64,eps:f64)
    ->Result<[[f64;2];2],TripleStitchHold>{
    let dx=b.cx-a.cx;let dz=b.cz-a.cz;
    let separation=dx.hypot(dz);
    let diff=(a.r-b.r).abs();
    let sum=a.r+b.r;
    if !separation.is_finite() || !sum.is_finite()||
        !separation.is_normal() && separation!=0.0 {
        return Err(TripleStitchHold::UnsupportedCircleContact);
    }
    if (separation-diff).abs()<=eps*16.0 ||
        (separation-sum).abs()<=eps*16.0 {
        return Err(TripleStitchHold::NearTangency);
    }
    if separation<=diff+eps*16.0 || separation>=sum-eps*16.0 {
        return Err(TripleStitchHold::UnsupportedCircleContact);
    }
    let along=(separation*separation+a.r*a.r-b.r*b.r)/(2.0*separation);
    let sq=a.r*a.r-along*along;
    if !sq.is_finite() || sq<=eps*eps*1024.0{
        return Err(TripleStitchHold::NearTangency);
    }
    let h=sq.sqrt();
    let ux=dx/separation;let uz=dz/separation;
    let pts=[
        [a.cx+along*ux-h*uz,a.cz+along*uz+h*ux],
        [a.cx+along*ux+h*uz,a.cz+along*uz-h*ux]
    ];
    for p in pts{
        let candidate=Point3{x:p[0],z:p[1],y};
        if !candidate.is_finite()||
            a.gap(candidate).abs()>eps*8.0 ||
            b.gap(candidate).abs()>eps*8.0{
            return Err(TripleStitchHold::AmbiguousArcExposure);
        }
    }
    Ok(pts)
}
fn trace(
    uses:Vec<DirectedTripleArc>,arcs:&[TripleArc],
    empty_reason:TripleStitchHold
)->Result<Vec<TripleContour>,TripleStitchHold>{
    if uses.is_empty(){return Err(empty_reason);}
    let mut successor=BTreeMap::new();
    let mut incoming=BTreeMap::<usize,usize>::new();
    for (i,use_) in uses.iter().enumerate(){
        let arc=arcs.get(use_.arc).ok_or(TripleStitchHold::OpenPlanarContour)?;
        let start=if use_.reversed{arc.end_junction}else{arc.start_junction};
        let end=if use_.reversed{arc.start_junction}else{arc.end_junction};
        if successor.insert(start,i).is_some(){
            return Err(TripleStitchHold::OpenPlanarContour);
        }
        *incoming.entry(end).or_default()+=1;
    }
    if incoming.len()!=successor.len() ||
        incoming.values().any(|&count|count!=1) ||
        successor.keys().any(|key|incoming.get(key)!=Some(&1)){
        return Err(TripleStitchHold::OpenPlanarContour);
    }
    let mut used=vec![false;uses.len()];
    let mut loops=Vec::new();
    for initial in 0..uses.len(){
        if used[initial]{continue;}
        let mut members=Vec::new();
        let mut curr=initial;
        loop{
            if used[curr]{
                if curr!=initial{return Err(TripleStitchHold::OpenPlanarContour);}
                break;
            }
            used[curr]=true;
            let part=uses[curr];
            let arc=&arcs[part.arc];
            let ending=if part.reversed{arc.start_junction}else{arc.end_junction};
            members.push(part);
            curr=*successor.get(&ending)
                .ok_or(TripleStitchHold::OpenPlanarContour)?;
            if members.len()>uses.len(){return Err(TripleStitchHold::OpenPlanarContour);}
        }
        let signed=kahan(members.iter().map(|m|{
            if m.reversed{-arcs[m.arc].planar_area_contribution}
            else{arcs[m.arc].planar_area_contribution}
        }));
        if !signed.is_finite(){return Err(TripleStitchHold::InconsistentAnalyticArea);}
        loops.push(TripleContour{parts:members,signed_area:signed});
    }
    Ok(loops)
}
fn put_use(map:&mut BTreeMap<TriplePhysicalEdgeKey,Vec<i8>>,
           key:TriplePhysicalEdgeKey,sign:i8) {
    map.entry(key).or_default().push(sign);
}
fn held(source:&Solid,input:TripleChainParameters,eps:f64,reason:TripleStitchHold)
    ->TripleStitchPlan {
    TripleStitchPlan{
        state:TripleStitchState::Held(reason),input,
        source_revision:source.revision,source_bounds:source.bbox,
        source_volume:source.mass.volume,tolerance:source.tolerance,
        length_budget:eps,levels:Vec::new(),junctions:Vec::new(),
        atomic_arcs:Vec::new(),slabs:Vec::new(),floors:Vec::new(),
        bottom_cap:None,top_cap:None,walls:Vec::new(),
        shared_physical_edges:Vec::new(),
        predicted_shell_vertices:None,predicted_shell_edges:None,
        predicted_shell_faces:None,predicted_euler:None,
        predicted_genus:None,analytic_surface_area:None,
        predicted_remaining_volume:None
    }
}
impl ConnectedTripleStitchPlanner{
    pub fn plan(source:&Solid,input:TripleChainParameters)
        ->Result<TripleStitchPlan,GeometryError>{
        crate::boolean::certify_box(source)
            .map_err(|e|invalid(format!("triple stitch requires true certified box: {e}")))?;
        source.tolerance.validate()?;
        let b=source.bbox;
        let height=b.max.y-b.min.y;
        let span=b.max.sub(b.min);
        let extent=span.x.max(span.y).max(span.z);
        let t=input.through;
        let a=input.deeper_blind;
        let c=input.shallower_blind;
        let disks=[
            Disk{cx:t.center_x,cz:t.center_z,r:t.radius},
            Disk{cx:a.center_x,cz:a.center_z,r:a.radius},
            Disk{cx:c.center_x,cz:c.center_z,r:c.radius},
        ];
        let coords=[b.min.x,b.max.x,b.min.y,b.max.y,b.min.z,b.max.z,
            t.center_x,t.center_z,a.center_x,a.center_z,c.center_x,c.center_z];
        if !coords.into_iter().all(f64::is_finite) ||
            !a.depth.is_finite()||!c.depth.is_finite(){
            return Err(invalid("nonfinite triple chain source geometry"));
        }
        let largest=coords.into_iter().map(f64::abs).fold(1.0_f64,f64::max);
        let eps=source.tolerance.length_at(extent)
            .max(64.0*f64::EPSILON*largest);
        for d in disks {
            if !d.r.is_finite()||d.r<=eps||
                d.cx-d.r<=b.min.x+eps||d.cx+d.r>=b.max.x-eps||
                d.cz-d.r<=b.min.z+eps||d.cz+d.r>=b.max.z-eps {
                return Err(invalid("all three disks need certified box wall clearance"));
            }
        }
        if !height.is_finite()||a.depth<=eps*32.0||
            c.depth<=eps*32.0||height-a.depth<=eps*32.0||
            height-c.depth<=eps*32.0{
            return Err(invalid("blind depths must be strictly interior"));
        }
        if a.depth-c.depth<=eps*32.0 {
            return Ok(held(source,input,eps,TripleStitchHold::UnsupportedCircleContact));
        }
        let tb=(disks[0].cx-disks[2].cx)
            .hypot(disks[0].cz-disks[2].cz);
        if (tb-disks[0].r-disks[2].r).abs()<=eps*32.0 {
            return Ok(held(source,input,eps,TripleStitchHold::NearTangency));
        }
        if tb<=disks[0].r+disks[2].r+eps*32.0 {
            return Ok(held(source,input,eps,TripleStitchHold::UnsupportedCircleContact));
        }
        let ta=match pair_intersections(disks[0],disks[1],b.max.y,eps){
            Ok(v)=>v,Err(why)=>return Ok(held(source,input,eps,why))
        };
        let ab=match pair_intersections(disks[1],disks[2],b.max.y,eps){
            Ok(v)=>v,Err(why)=>return Ok(held(source,input,eps,why))
        };
        let planning=MultiCutRegionPlanner::plan(source,&[t],&[a,c])?;
        if planning.state!=MultiCutPlanState::Ready||planning.slabs.len()!=3 {
            return Ok(held(source,input,eps,TripleStitchHold::UnresolvedSourceArrangement));
        }
        let levels=vec![b.min.y,b.max.y-a.depth,b.max.y-c.depth,b.max.y];
        if !levels.windows(2).all(|w|w[1]-w[0]>eps*32.0){
            return Ok(held(source,input,eps,TripleStitchHold::UnresolvedSourceArrangement));
        }
        let mut junctions=Vec::new();
        for (pair,points) in [([0,1],ta),([1,2],ab)] {
            for point in points {
                let id=junctions.len();
                junctions.push(TripleJunction{
                    id,source_pair:pair,position_xz:point
                });
            }
        }
        let mut events=[Vec::<Event>::new(),Vec::<Event>::new(),Vec::<Event>::new()];
        for j in &junctions {
            for &i in &j.source_pair {
                events[i].push(Event{
                    angle:angles(disks[i],j.position_xz),
                    junction:j.id
                });
            }
        }
        let mut atomic_arcs=Vec::new();
        let anchor=(t.center_x,t.center_z);
        for disk_index in 0..3 {
            let disk=disks[disk_index];
            events[disk_index].sort_by(|x,y|x.angle.total_cmp(&y.angle));
            if events[disk_index].len()!=if disk_index==1{4}else{2}{
                return Ok(held(source,input,eps,TripleStitchHold::OpenPlanarContour));
            }
            for k in 0..events[disk_index].len(){
                let prev=events[disk_index][k];
                let next=events[disk_index][(k+1)%events[disk_index].len()];
                let end=if k+1==events[disk_index].len(){
                    next.angle+TAU
                }else{next.angle};
                let delta=end-prev.angle;
                if !delta.is_finite()||delta*disk.r<=eps*32.0{
                    return Ok(held(source,input,eps,TripleStitchHold::AmbiguousArcExposure));
                }
                let midpoint=disk.point(prev.angle+delta*0.5,b.max.y);
                let mut visible=[false;3];
                for slab in 0..3 {
                    let active=match slab{0=>&[0][..],1=>&[0,1][..],
                        _=>&[0,1,2][..]};
                    if !active.contains(&disk_index){continue;}
                    let mut covered=false;
                    for &other in active{
                        if other==disk_index{continue;}
                        let gap=disks[other].gap(midpoint);
                        if !gap.is_finite()||gap.abs()<=eps*8.0{
                            return Ok(held(source,input,eps,
                                TripleStitchHold::AmbiguousArcExposure));
                        }
                        if gap<0.0{covered=true;}
                    }
                    visible[slab]=!covered;
                }
                atomic_arcs.push(TripleArc{
                    id:atomic_arcs.len(),source:disk_index,
                    start_junction:prev.junction,end_junction:next.junction,
                    start_angle:prev.angle,end_angle:end,
                    arc_length:delta*disk.r,
                    planar_area_contribution:green(disk,prev.angle,end,anchor),
                    visible_in_slab:visible
                });
            }
        }
        let mut slabs=Vec::new();
        let area_budget=eps*(disks.iter().map(|d|d.r).sum::<f64>())*256.0+
            span.x*span.z*1024.0*f64::EPSILON;
        for slab in 0..3 {
            let uses=atomic_arcs.iter().filter(|a|a.visible_in_slab[slab])
                .map(|a|DirectedTripleArc{arc:a.id,reversed:false})
                .collect::<Vec<_>>();
            let loops=match trace(uses,&atomic_arcs,TripleStitchHold::OpenPlanarContour){
                Ok(x)=>x,Err(e)=>return Ok(held(source,input,eps,e))
            };
            if loops.len()!=1||loops[0].signed_area<=0.0{
                return Ok(held(source,input,eps,TripleStitchHold::OpenPlanarContour));
            }
            let area=loops[0].signed_area;
            if (area-planning.slabs[slab].union_area).abs()>area_budget{
                return Ok(held(source,input,eps,TripleStitchHold::InconsistentAnalyticArea));
            }
            let removed=area*(levels[slab+1]-levels[slab]);
            slabs.push(TripleSlab{
                y_lower:levels[slab],y_upper:levels[slab+1],
                active_sources:(0..=slab).collect(),
                exposed_contour:loops[0].clone(),union_area:area,
                removed_volume:removed
            });
        }
        let mut floors=Vec::new();
        for interface in 1..3{
            let mut floor=Vec::new();
            let mut lost=Vec::new();
            let mut born=Vec::new();
            for a in &atomic_arcs {
                let lower=a.visible_in_slab[interface-1];
                let upper=a.visible_in_slab[interface];
                match (lower,upper) {
                    (true,false)=>{
                        floor.push(DirectedTripleArc{arc:a.id,reversed:true});
                        lost.push(a.id);
                    },
                    (false,true)=>{
                        floor.push(DirectedTripleArc{arc:a.id,reversed:false});
                        born.push(a.id);
                    },
                    _=>{}
                }
            }
            let loops=match trace(floor,&atomic_arcs,TripleStitchHold::OpenPlanarContour){
                Ok(x)=>x,Err(reason)=>return Ok(held(source,input,eps,reason))
            };
            if loops.len()!=1||loops[0].signed_area<=area_budget {
                return Ok(held(source,input,eps,TripleStitchHold::OpenPlanarContour));
            }
            let floor_area=slabs[interface].union_area-slabs[interface-1].union_area;
            if floor_area<=area_budget ||
                (loops[0].signed_area-floor_area).abs()>area_budget {
                return Ok(held(source,input,eps,TripleStitchHold::InconsistentAnalyticArea));
            }
            floors.push(TripleFloorInterface{
                y:levels[interface],contour:loops[0].clone(),
                floor_area,incoming_lower_arcs:lost,born_upper_arcs:born
            });
        }
        let mut walls=Vec::new();
        let mut edge_uses=BTreeMap::<TriplePhysicalEdgeKey,Vec<i8>>::new();
        for slab in 0..3{
            for a in &atomic_arcs{
                if !a.visible_in_slab[slab]{continue;}
                let disk=disks[a.source];
                let corners=[
                    disk.point(a.start_angle,levels[slab]),
                    disk.point(a.end_angle,levels[slab]),
                    disk.point(a.start_angle,levels[slab+1]),
                    disk.point(a.end_angle,levels[slab+1])
                ];
                let v0=b.max.y-levels[slab+1];
                let v1=b.max.y-levels[slab];
                let face=TripleWallPatch{
                    source:a.source,arc:a.id,slab,
                    y_lower:levels[slab],y_upper:levels[slab+1],
                    uv_bounds:[a.start_angle,a.end_angle,v0,v1],
                    corners,analytic_area:a.arc_length*(levels[slab+1]-levels[slab])
                };
                for p in &face.corners{
                    if !p.is_finite() || disk.gap(*p).abs()>eps*8.0{
                        return Ok(held(source,input,eps,TripleStitchHold::AmbiguousArcExposure));
                    }
                }
                walls.push(face);
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                    arc:a.id,level:slab},1);
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                    arc:a.id,level:slab+1},-1);
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::VerticalTrace{
                    junction:a.start_junction,slab},-1);
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::VerticalTrace{
                    junction:a.end_junction,slab},1);
            }
        }
        for a in &atomic_arcs{
            if a.visible_in_slab[0]{
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                    arc:a.id,level:0},-1); // canonical bottom cap inner loop
            }
            if a.visible_in_slab[2]{
                put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                    arc:a.id,level:3},1); // canonical top cap inner loop
            }
            for interface in 1..3{
                let (lower,upper)=(a.visible_in_slab[interface-1],
                    a.visible_in_slab[interface]);
                if lower&&!upper {
                    put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                        arc:a.id,level:interface},1); // physical floor, +Y
                }
                if upper&&!lower {
                    put_use(&mut edge_uses,TriplePhysicalEdgeKey::CircleArc{
                        arc:a.id,level:interface},-1);
                }
            }
        }
        if edge_uses.values().any(|uses|uses.len()!=2||
            uses[0]!=-uses[1]){
            return Ok(held(source,input,eps,TripleStitchHold::NonmanifoldPhysicalEdge));
        }
        let mut vertices=BTreeSet::<(usize,usize)>::new();
        for key in edge_uses.keys(){
            match *key{
                TriplePhysicalEdgeKey::CircleArc{arc,level}=>{
                    let a=&atomic_arcs[arc];
                    vertices.insert((a.start_junction,level));
                    vertices.insert((a.end_junction,level));
                }
                TriplePhysicalEdgeKey::VerticalTrace{junction,slab}=>{
                    vertices.insert((junction,slab));
                    vertices.insert((junction,slab+1));
                }
            }
        }
        // Base canonical box has V=8, E=12, F=6, and two physical
        // cap-inner loops must be subtracted from Euler's relation.
        let count_v=8+vertices.len();
        let count_e=12+edge_uses.len();
        let count_f=6+walls.len()+floors.len();
        let chi=count_v as isize-count_e as isize+count_f as isize-2;
        if chi!=0 {
            return Ok(held(source,input,eps,
                TripleStitchHold::InconsistentEulerCharacteristic));
        }
        let remaining=source.mass.volume-
            kahan(slabs.iter().map(|s|s.removed_volume));
        let from_existing=planning.predicted_remaining_volume
            .ok_or_else(||invalid("missing v0.32 independent mass witness"))?;
        let volume_budget=area_budget*height*16.0;
        if remaining<=0.0||!remaining.is_finite()||
            (remaining-from_existing).abs()>volume_budget{
            return Ok(held(source,input,eps,TripleStitchHold::InconsistentAnalyticArea));
        }
        let perimeter_area=kahan(walls.iter().map(|p|p.analytic_area));
        let floor_area=kahan(floors.iter().map(|f|f.floor_area));
        let surface_area=source.mass.surface_area-slabs[0].union_area-
            slabs[2].union_area+perimeter_area+floor_area;
        if !surface_area.is_finite()||surface_area<=0.0{
            return Ok(held(source,input,eps,TripleStitchHold::InconsistentAnalyticArea));
        }
        let shared_physical_edges=edge_uses.into_iter().map(|(key,uses)|
            TripleEdgePair{key,opposite_use_signs:[uses[0],uses[1]]}
        ).collect();
        Ok(TripleStitchPlan{
            state:TripleStitchState::Ready,input,
            source_revision:source.revision,source_bounds:source.bbox,
            source_volume:source.mass.volume,tolerance:source.tolerance,
            length_budget:eps,levels,junctions,atomic_arcs,
            bottom_cap:Some(slabs[0].exposed_contour.clone()),
            top_cap:Some(slabs[2].exposed_contour.clone()),
            slabs,floors,walls,shared_physical_edges,
            predicted_shell_vertices:Some(count_v),
            predicted_shell_edges:Some(count_e),
            predicted_shell_faces:Some(count_f),
            predicted_euler:Some(chi),
            predicted_genus:Some(1),
            analytic_surface_area:Some(surface_area),
            predicted_remaining_volume:Some(remaining)
        })
    }
}
