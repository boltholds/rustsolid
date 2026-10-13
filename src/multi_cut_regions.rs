//! Read-only planar circle-union arrangement, evaluated independently in
//! horizontal constant-active-cutter slabs for multiple through/blind bores.
//! The signed circular-arc Green integral is analytical, not polygonized.
//!
//! All commits to B-rep topology remain disabled. Near tangencies, triple
//! source junctions and numerically unresolved events are explicitly held.
use std::collections::BTreeMap;
use std::f64::consts::{PI,TAU};

use crate::{
    BlindHoleParameters,BoundingBox,GeometryError,GeometryTolerance,
    HoleParameters,Point3,Solid,
};

const MAX_CUTTERS:usize=12;
fn invalid(msg:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(msg.into())}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiCutSource { Through(usize), Blind(usize) }

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiCutHold {
    NearTangencyOrCoincidence,
    NearTripleJunction,
    UnresolvedArc,
    IncompleteBoundaryCycles,
    UnresolvedDepthEvent,
    NumericalAreaAmbiguity,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MultiCutPlanState { Ready, Held(MultiCutHold) }

/// Source-bound intersection or virtual (un-cut) full-circle seam.
#[derive(Debug,Clone,PartialEq)]
pub struct CircleUnionVertex {
    pub id:usize,
    pub disk_a:usize,
    pub disk_b:usize,
    pub point:Point3,
    pub is_virtual_seam:bool,
}

/// Directed CCW portion of a real input circle on the union boundary.
/// The interior of the union is always on the left of this directed arc.
/// start_angle may be negative/positive; end_angle > start_angle, and
/// the interval may wrap past 2*pi without altering its source geometry.
#[derive(Debug,Clone,PartialEq)]
pub struct CircleUnionArc {
    pub disk_index:usize,
    pub start_vertex:usize,
    pub end_vertex:usize,
    pub start_angle:f64,
    pub end_angle:f64,
    pub midpoint:Point3,
    pub signed_area_contribution:f64,
}

#[derive(Debug,Clone,PartialEq)]
pub struct CircleUnionLoop {
    pub arcs:Vec<CircleUnionArc>,
    pub signed_area:f64,
    pub is_hole:bool,
}

/// One vertical slab over which all active circular cutters are constant.
/// Loops with positive signed area bound connected occupied disk unions;
/// negative loops bound voids in those planar unions.
#[derive(Debug,Clone,PartialEq)]
pub struct MultiCutRegionSlab {
    pub y_lower:f64,
    pub y_upper:f64,
    /// Indices into the source sequence: through first, blind second.
    pub active_disks:Vec<usize>,
    pub vertices:Vec<CircleUnionVertex>,
    pub boundary_loops:Vec<CircleUnionLoop>,
    pub connected_union_components:usize,
    pub planar_holes:usize,
    pub union_area:f64,
    pub removed_volume:f64,
    /// A numerical resolution scale, NOT a mathematically rigorous
    /// interval-arithmetic bound on all floating-point error.
    pub area_resolution_budget:f64,
}

#[derive(Debug,Clone,PartialEq)]
pub struct MultiCutRegionPlan {
    pub state:MultiCutPlanState,
    pub slabs:Vec<MultiCutRegionSlab>,
    pub removed_volume:Option<f64>,
    pub predicted_remaining_volume:Option<f64>,
    pub source_bounds:BoundingBox,
    pub source_volume:f64,
    pub source_revision:u64,
    pub tolerance:GeometryTolerance,
    pub length_budget:f64,
    pub through:Vec<HoleParameters>,
    pub blind:Vec<BlindHoleParameters>,
}
impl MultiCutRegionPlan {
    /// Source-backed mathematical geometry is NOT executable curved B-rep.
    pub fn can_apply(&self)->bool{false}
    pub fn validate_against(&self,source:&Solid,through:&[HoleParameters],
                            blind:&[BlindHoleParameters])->Result<(),GeometryError>{
        let expected=MultiCutRegionPlanner::plan(source,through,blind)?;
        if *self!=expected{return Err(invalid("multi-cut plan differs from its source"));}
        Ok(())
    }
}

#[derive(Debug,Clone,Copy,Default)]
pub struct MultiCutRegionPlanner;

#[derive(Debug,Clone,Copy)]
struct Disk { x:f64,z:f64,r:f64,source:MultiCutSource }
impl Disk {
    fn on(self,theta:f64,y:f64)->Point3{
        Point3{x:self.x+self.r*theta.cos(),y,z:self.z+self.r*theta.sin()}
    }
    fn gap(self,point:Point3)->f64{
        (point.x-self.x).hypot(point.z-self.z)-self.r
    }
    fn area(self)->f64{PI*self.r*self.r}
}
#[derive(Debug,Clone,Copy)]
struct Event {theta:f64,vertex:usize}
fn angle(point:Point3,disk:Disk)->f64{
    (point.z-disk.z).atan2(point.x-disk.x).rem_euclid(TAU)
}
fn sum_compensated(values:impl IntoIterator<Item=f64>)->f64{
    let (mut total,mut correction)=(0.0,0.0);
    for value in values{
        let adjusted=value-correction;
        let next=total+adjusted;
        correction=(next-total)-adjusted;
        total=next;
    }
    total
}
fn arc_integral(d:Disk,start:f64,end:f64,anchor:(f64,f64))->f64{
    // Green theorem: A = 1/2 integral (x dz - z dx). The shared
    // anchor prevents catastrophic cancellation under large translation.
    let cx=d.x-anchor.0;
    let cz=d.z-anchor.1;
    let delta=end-start;
    0.5*(d.r*d.r*delta+
        d.r*cx*(end.sin()-start.sin())+
        d.r*cz*(start.cos()-end.cos()))
}
fn hold()->MultiCutHold{MultiCutHold::NumericalAreaAmbiguity}

/// A closed planar contour graph reconstructed from source circle-circle
/// intersection events. Zero-length/seam fake edges are never committed.
fn union_slab(disks:&[Disk],active:&[usize],y:f64,eps:f64,area_scale:f64)
    ->Result<(Vec<CircleUnionVertex>,Vec<CircleUnionLoop>,f64,usize,usize,f64),MultiCutHold>
{
    let mut vertices=Vec::<CircleUnionVertex>::new();
    let mut events=vec![Vec::<Event>::new();disks.len()];
    for (offset,&i) in active.iter().enumerate(){
        let a=disks[i];
        for &j in active.iter().skip(offset+1){
            let b=disks[j];
            let dx=b.x-a.x;let dz=b.z-a.z;
            let distance=dx.hypot(dz);
            let sum=a.r+b.r;let difference=(a.r-b.r).abs();
            if !distance.is_finite()||!sum.is_finite(){return Err(hold());}
            if (distance<=eps && difference<=eps) ||
                (distance-sum).abs()<=eps*4.0 ||
                (distance-difference).abs()<=eps*4.0 {
                return Err(MultiCutHold::NearTangencyOrCoincidence);
            }
            if distance>sum+eps*4.0 || distance<difference-eps*4.0{
                continue;
            }
            if distance<=eps{return Err(MultiCutHold::NearTangencyOrCoincidence);}
            let along=(distance*distance+a.r*a.r-b.r*b.r)/(2.0*distance);
            let opposing=distance-along;
            let square=a.r*a.r-along*along;
            let sq_b=b.r*b.r-opposing*opposing;
            if !square.is_finite()||!sq_b.is_finite()||square<=0.0||sq_b<=0.0{
                return Err(MultiCutHold::UnresolvedArc);
            }
            let h=square.sqrt();let h_b=sq_b.sqrt();
            if h<=eps*32.0||h_b<=eps*32.0||(h-h_b).abs()>eps*8.0{
                return Err(MultiCutHold::NearTangencyOrCoincidence);
            }
            let ux=dx/distance;let uz=dz/distance;
            let cx=a.x+along*ux;let cz=a.z+along*uz;
            for sign in [-1.0,1.0] {
                let point=Point3{x:cx-sign*h*uz,y,z:cz+sign*h*ux};
                if !point.is_finite() ||
                    a.gap(point).abs()>eps*8.0 ||
                    b.gap(point).abs()>eps*8.0 {
                    return Err(MultiCutHold::UnresolvedArc);
                }
                // Near triple circle crossings are an unresolved topology
                // transition, even when two crossing sources are certified.
                for &k in active{
                    if k!=i&&k!=j && disks[k].gap(point).abs()<=eps*32.0{
                        return Err(MultiCutHold::NearTripleJunction);
                    }
                }
                let id=vertices.len();
                vertices.push(CircleUnionVertex{
                    id,disk_a:i,disk_b:j,point,is_virtual_seam:false
                });
                events[i].push(Event{theta:angle(point,a),vertex:id});
                events[j].push(Event{theta:angle(point,b),vertex:id});
            }
        }
    }
    let mut arcs=Vec::<CircleUnionArc>::new();
    let anchor=(disks[active[0]].x,disks[active[0]].z);
    for &i in active{
        let disk=disks[i];
        let owns=|theta:f64|->Result<bool,MultiCutHold>{
            let point=disk.on(theta,y);
            for &j in active{
                if j==i{continue;}
                let gap=disks[j].gap(point);
                if !gap.is_finite()||gap.abs()<=eps*8.0{
                    return Err(MultiCutHold::UnresolvedArc);
                }
                if gap<0.0{return Ok(false);}
            }
            Ok(true)
        };
        if events[i].is_empty(){
            if owns(0.0)?{
                let id=vertices.len();
                let point=disk.on(0.0,y);
                vertices.push(CircleUnionVertex{
                    id,disk_a:i,disk_b:i,point,is_virtual_seam:true
                });
                arcs.push(CircleUnionArc{
                    disk_index:i,start_vertex:id,end_vertex:id,
                    start_angle:0.0,end_angle:TAU,
                    midpoint:disk.on(PI,y),
                    signed_area_contribution:arc_integral(disk,0.0,TAU,anchor),
                });
            }
            continue;
        }
        events[i].sort_by(|a,b|a.theta.total_cmp(&b.theta));
        if events[i].len()%2!=0{return Err(MultiCutHold::IncompleteBoundaryCycles);}
        for k in 0..events[i].len(){
            let from=events[i][k];
            let to=events[i][(k+1)%events[i].len()];
            let end=if k+1==events[i].len(){to.theta+TAU}else{to.theta};
            let span=end-from.theta;
            if !span.is_finite()||span<=0.0||
                span*disk.r<=eps*32.0 {
                return Err(MultiCutHold::UnresolvedArc);
            }
            if owns(from.theta+span*0.5)?{
                arcs.push(CircleUnionArc{
                    disk_index:i,start_vertex:from.vertex,end_vertex:to.vertex,
                    start_angle:from.theta,end_angle:end,
                    midpoint:disk.on(from.theta+span*0.5,y),
                    signed_area_contribution:arc_integral(disk,from.theta,end,anchor),
                });
            }
        }
    }
    if arcs.is_empty(){return Err(MultiCutHold::IncompleteBoundaryCycles);}
    let mut directed=BTreeMap::<usize,usize>::new();
    let mut incoming=BTreeMap::<usize,usize>::new();
    for (idx,arc) in arcs.iter().enumerate(){
        if directed.insert(arc.start_vertex,idx).is_some(){
            return Err(MultiCutHold::IncompleteBoundaryCycles);
        }
        *incoming.entry(arc.end_vertex).or_default()+=1;
    }
    for (&vertex,_) in &directed {
        if incoming.get(&vertex)!=Some(&1){
            return Err(MultiCutHold::IncompleteBoundaryCycles);
        }
    }
    if incoming.len()!=directed.len(){
        return Err(MultiCutHold::IncompleteBoundaryCycles);
    }
    let mut used=vec![false;arcs.len()];
    let mut loops=Vec::<CircleUnionLoop>::new();
    let area_budget=eps*active.iter().map(|&i|disks[i].r).sum::<f64>()*
        128.0 + area_scale*256.0*f64::EPSILON;
    for start in 0..arcs.len(){
        if used[start]{continue;}
        let mut entries=Vec::new();
        let mut cursor=start;
        loop{
            if used[cursor] {
                if cursor!=start {
                    return Err(MultiCutHold::IncompleteBoundaryCycles);
                }
                break;
            }
            used[cursor]=true;
            let arc=arcs[cursor].clone();
            let next=*directed.get(&arc.end_vertex)
                .ok_or(MultiCutHold::IncompleteBoundaryCycles)?;
            entries.push(arc);
            cursor=next;
            if entries.len()>arcs.len(){
                return Err(MultiCutHold::IncompleteBoundaryCycles);
            }
        }
        let signed_area=sum_compensated(
            entries.iter().map(|a|a.signed_area_contribution)
        );
        if !signed_area.is_finite()||signed_area.abs()<=area_budget{
            return Err(MultiCutHold::NumericalAreaAmbiguity);
        }
        loops.push(CircleUnionLoop{
            arcs:entries,signed_area,is_hole:signed_area<0.0
        });
    }
    if used.iter().any(|&u|!u){return Err(MultiCutHold::IncompleteBoundaryCycles);}
    let union_area=sum_compensated(loops.iter().map(|l|l.signed_area));
    let component_count=loops.iter().filter(|l|!l.is_hole).count();
    let holes=loops.len()-component_count;
    let sum_areas: f64=active.iter().map(|&i|disks[i].area()).sum();
    let largest=active.iter().map(|&i|disks[i].area()).fold(0.0_f64,f64::max);
    if !union_area.is_finite()||union_area<largest-area_budget ||
        union_area>sum_areas+area_budget ||
        union_area>area_scale+area_budget ||
        component_count==0 {
        return Err(MultiCutHold::NumericalAreaAmbiguity);
    }
    Ok((vertices,loops,union_area,component_count,holes,area_budget))
}
impl MultiCutRegionPlanner{
    pub fn plan(source:&Solid,through:&[HoleParameters],blind:&[BlindHoleParameters])
        ->Result<MultiCutRegionPlan,GeometryError>
    {
        crate::boolean::certify_box(source)
            .map_err(|e|invalid(format!("multi-cut source must be exact box: {e}")))?;
        let tol=source.tolerance;
        tol.validate()?;
        if through.is_empty()||blind.is_empty()||
            through.len()+blind.len()>MAX_CUTTERS{
            return Err(invalid("multi-cut planner requires 1..=12 total cutters and both families"));
        }
        let b=source.bbox;
        let width=b.max.x-b.min.x;
        let height=b.max.y-b.min.y;
        let depth=b.max.z-b.min.z;
        let extent=width.max(height).max(depth);
        let mut largest=[b.min.x,b.max.x,b.min.z,b.max.z,b.min.y,b.max.y]
            .into_iter().map(f64::abs).fold(1.0_f64,f64::max);
        let mut disks=Vec::<Disk>::new();
        for (i,t) in through.iter().enumerate(){
            if !t.center_x.is_finite()||!t.center_z.is_finite()||
                !t.radius.is_finite(){return Err(invalid("nonfinite through cutter"));}
            largest=largest.max(t.center_x.abs()).max(t.center_z.abs());
            disks.push(Disk{x:t.center_x,z:t.center_z,r:t.radius,
                source:MultiCutSource::Through(i)});
        }
        let mut floors=Vec::new();
        for (i,cut) in blind.iter().enumerate(){
            if !cut.center_x.is_finite()||!cut.center_z.is_finite()||
                !cut.radius.is_finite()||!cut.depth.is_finite() {
                return Err(invalid("nonfinite blind cutter"));
            }
            largest=largest.max(cut.center_x.abs()).max(cut.center_z.abs());
            disks.push(Disk{x:cut.center_x,z:cut.center_z,r:cut.radius,
                source:MultiCutSource::Blind(i)});
            floors.push(b.max.y-cut.depth);
        }
        if !largest.is_finite()||!extent.is_finite()||height<=0.0{
            return Err(invalid("nonfinite source bounds"));
        }
        let eps=tol.length_at(extent)
            .max(64.0*f64::EPSILON*largest);
        for d in &disks{
            if d.r<=eps || !d.area().is_finite() ||
                d.x-d.r<=b.min.x+eps || d.x+d.r>=b.max.x-eps ||
                d.z-d.r<=b.min.z+eps || d.z+d.r>=b.max.z-eps {
                return Err(invalid("each disk must be strictly within the box sides"));
            }
        }
        for (i,cut) in blind.iter().enumerate(){
            if cut.depth<=eps||height-cut.depth<=eps ||
                floors[i]<=b.min.y+eps || floors[i]>=b.max.y-eps {
                return Err(invalid("blind depth must be strictly interior"));
            }
        }
        let mut plan=MultiCutRegionPlan {
            state:MultiCutPlanState::Ready,slabs:Vec::new(),
            removed_volume:None,predicted_remaining_volume:None,
            source_bounds:b,source_volume:source.mass.volume,
            source_revision:source.revision,tolerance:tol,length_budget:eps,
            through:through.to_vec(),blind:blind.to_vec()
        };
        let mut levels=Vec::with_capacity(floors.len()+2);
        levels.push(b.min.y);levels.extend(floors);levels.push(b.max.y);
        levels.sort_by(f64::total_cmp);
        for w in levels.windows(2) {
            if w[1]-w[0]>0.0 && w[1]-w[0]<=eps*4.0 {
                plan.state=MultiCutPlanState::Held(MultiCutHold::UnresolvedDepthEvent);
                return Ok(plan);
            }
        }
        levels.dedup_by(|a,c|*a==*c);
        let area_scale=width*depth;
        let mut removed_by_slab=Vec::new();
        for pair in levels.windows(2){
            let (y0,y1)=(pair[0],pair[1]);
            let midpoint=y0+(y1-y0)*0.5;
            let mut active=(0..through.len()).collect::<Vec<_>>();
            for (i,cut) in blind.iter().enumerate(){
                if midpoint>b.max.y-cut.depth{
                    active.push(through.len()+i);
                }
            }
            let (vertices,loops,area,components,holes,budget)=
                match union_slab(&disks,&active,midpoint,eps,area_scale){
                    Ok(value)=>value,
                    Err(reason)=>{
                        plan.state=MultiCutPlanState::Held(reason);
                        plan.slabs.clear();
                        return Ok(plan);
                    }
                };
            let volume=area*(y1-y0);
            if !volume.is_finite(){
                plan.state=MultiCutPlanState::Held(MultiCutHold::NumericalAreaAmbiguity);
                return Ok(plan);
            }
            removed_by_slab.push(volume);
            plan.slabs.push(MultiCutRegionSlab{
                y_lower:y0,y_upper:y1,active_disks:active,vertices,
                boundary_loops:loops,connected_union_components:components,
                planar_holes:holes,union_area:area,removed_volume:volume,
                area_resolution_budget:budget
            });
        }
        let removed=sum_compensated(removed_by_slab);
        let remaining=source.mass.volume-removed;
        let volume_budget=eps*area_scale*512.0+
            source.mass.volume*1024.0*f64::EPSILON;
        if !removed.is_finite()||!remaining.is_finite()||
            removed<=0.0||remaining<=0.0 ||
            remaining>source.mass.volume+volume_budget {
            plan.state=MultiCutPlanState::Held(MultiCutHold::NumericalAreaAmbiguity);
            plan.slabs.clear();
            return Ok(plan);
        }
        plan.removed_volume=Some(removed);
        plan.predicted_remaining_volume=Some(remaining);
        Ok(plan)
    }
}
