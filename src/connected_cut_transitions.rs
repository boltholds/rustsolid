//! Source-certified transitions between adjacent analytical circular-union
//! slabs. These are *geometric face-splitting instructions*, not a committed
//! B-rep. A global curved topology builder must validate coedges and UV charts
//! before constructing a shell.
use std::f64::consts::{PI,TAU};
use crate::{
    BlindHoleParameters,CircleUnionArc,CircleUnionLoop,GeometryError,
    HoleParameters,MultiCutPlanState,MultiCutRegionPlan,MultiCutRegionPlanner,
    Point3,Solid,
};
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TransitionArcKind {
    /// Both sides expose precisely the same source-owned circular arc.
    Continuous,
    /// Exposed cylindrical wall terminates at this blind floor.
    Terminating,
    /// A new cylindrical wall begins above this blind floor.
    Emerging,
}
#[derive(Debug,Clone,PartialEq)]
pub struct TransitionArc {
    pub disk_index:usize,
    pub theta_start:f64,
    pub theta_end:f64,
    pub kind:TransitionArcKind,
    pub start:Point3,
    pub end:Point3,
    /// This is the *physical supporting cylinder*, not a new topology edge.
    pub source_radius:f64,
}
#[derive(Debug,Clone,PartialEq)]
pub struct DepthTransition {
    pub height_y:f64,
    pub activated_blind:Vec<usize>,
    pub lower_slab:usize,
    pub upper_slab:usize,
    pub arcs:Vec<TransitionArc>,
    /// Area of new void opening above this floor. The actual planar material
    /// floor is a subset of this region with upward material normal.
    pub new_void_area:f64,
    pub floor_material_area:f64,
    pub component_count_before:usize,
    pub component_count_after:usize,
    pub planar_holes_before:usize,
    pub planar_holes_after:usize,
}
#[derive(Debug,Clone,PartialEq)]
pub struct ConnectedCutTransitionPlan {
    pub source_plan:MultiCutRegionPlan,
    pub transitions:Vec<DepthTransition>,
    pub predicted_volume:f64,
    pub can_apply:bool,
}
impl ConnectedCutTransitionPlan {
    pub fn validate_against(&self,box_:&Solid,
        through:&[HoleParameters],blind:&[BlindHoleParameters])
        ->Result<(),GeometryError>{
        if self!=&ConnectedCutTransitionPlanner::plan(box_,through,blind)?{
            return Err(invalid("transition witnesses differ from independent source replay"));
        }
        Ok(())
    }
}
#[derive(Debug,Clone,Copy,Default)]
pub struct ConnectedCutTransitionPlanner;

fn source_arc_intervals(loops:&[CircleUnionLoop],disk:usize)->Vec<(f64,f64)>{
    let mut spans=loops.iter().flat_map(|lp|lp.arcs.iter())
        .filter(|a|a.disk_index==disk)
        .map(|a|(a.start_angle,a.end_angle)).collect::<Vec<_>>();
    spans.sort_by(|a,b|a.0.total_cmp(&b.0));
    spans
}
fn intervals_on_circle(loops:&[CircleUnionLoop],disk:usize)
    ->Vec<(f64,f64)>{
    // Split wrapped intervals into the canonical [0,2*pi] angular chart.
    // Endpoints from different slabs need not have identical subdivision.
    let mut result=Vec::new();
    for (start,end) in source_arc_intervals(loops,disk){
        if !start.is_finite()||!end.is_finite()||end<=start {continue;}
        let a=start.rem_euclid(TAU);
        let span=end-start;
        if span>=TAU-32.0*f64::EPSILON{
            return vec![(0.0,TAU)];
        }
        if a+span<=TAU {
            result.push((a,a+span));
        }else{
            result.push((a,TAU));
            result.push((0.0,a+span-TAU));
        }
    }
    result.sort_by(|a,b|a.0.total_cmp(&b.0));
    result
}
fn source_disk(index:usize,through:&[HoleParameters],
    blind:&[BlindHoleParameters])->(f64,f64,f64){
    if index<through.len(){
        let x=through[index];(x.center_x,x.center_z,x.radius)
    }else{
        let x=blind[index-through.len()];
        (x.center_x,x.center_z,x.radius)
    }
}
fn same_interval_value(mid:f64,spans:&[(f64,f64)],angular_eps:f64)->bool{
    spans.iter().any(|&(a,b)|mid>a+angular_eps && mid<b-angular_eps)
}
impl ConnectedCutTransitionPlanner {
    pub fn plan(box_:&Solid,through:&[HoleParameters],
        blind:&[BlindHoleParameters])->Result<ConnectedCutTransitionPlan,GeometryError>
    {
        let plan=MultiCutRegionPlanner::plan(box_,through,blind)?;
        if plan.state!=MultiCutPlanState::Ready{
            return Err(invalid("global source arrangement must be Ready before cross-layer transitions"));
        }
        plan.validate_against(box_,through,blind)?;
        if plan.slabs.len()<2||blind.len()<2{
            return Err(invalid("connected triple planning requires at least two blind cuts and height slabs"));
        }
        let eps=plan.length_budget;
        let mut transitions=Vec::new();
        let total=through.len()+blind.len();
        let mut found_connected_three=false;
        for i in 0..plan.slabs.len()-1{
            let lower=&plan.slabs[i];
            let upper=&plan.slabs[i+1];
            let height=lower.y_upper;
            if (height-upper.y_lower).abs()>eps*4.0{
                return Err(invalid("consecutive depth slabs have a gap"));
            }
            let activated=upper.active_disks.iter().copied()
                .filter(|j|!lower.active_disks.contains(j))
                .collect::<Vec<_>>();
            if activated.is_empty(){
                return Err(invalid("depth transition activates no source cutter"));
            }
            if activated.iter().any(|&j|j<through.len()) {
                return Err(invalid("full through cutters cannot activate mid-height"));
            }
            let mut arcs=Vec::new();
            let angle_eps=32.0*f64::EPSILON;
            for disk in 0..total {
                let before=if lower.active_disks.contains(&disk){
                    intervals_on_circle(&lower.boundary_loops,disk)
                }else{Vec::new()};
                let after=if upper.active_disks.contains(&disk){
                    intervals_on_circle(&upper.boundary_loops,disk)
                }else{Vec::new()};
                let mut cuts=vec![0.0,TAU];
                for &(a,b) in before.iter().chain(&after){
                    cuts.push(a);cuts.push(b);
                }
                cuts.sort_by(f64::total_cmp);
                cuts.dedup_by(|a,b|(*a-*b).abs()<=angle_eps);
                let (x,z,r)=source_disk(disk,through,blind);
                let angular_tol=(eps/r).max(angle_eps)*32.0;
                for pair in cuts.windows(2){
                    let a=pair[0];let b=pair[1];
                    if b-a<=angular_tol {
                        // A boundary chart too fine to certify would result
                        // in a false topological edge if published.
                        return Err(invalid("depth transition has sub-resolution angular interval"));
                    }
                    let mid=a+(b-a)*0.5;
                    let below=same_interval_value(mid,&before,angle_eps);
                    let above=same_interval_value(mid,&after,angle_eps);
                    if !(below||above){continue;}
                    let kind=match (below,above){
                        (true,true)=>TransitionArcKind::Continuous,
                        (true,false)=>TransitionArcKind::Terminating,
                        (false,true)=>TransitionArcKind::Emerging,
                        (false,false)=>continue,
                    };
                    let on=|theta:f64|Point3{
                        x:x+r*theta.cos(),y:height,z:z+r*theta.sin()
                    };
                    arcs.push(TransitionArc{
                        disk_index:disk,theta_start:a,theta_end:b,
                        kind,start:on(a),end:on(b),source_radius:r
                    });
                }
            }
            let area=upper.union_area-lower.union_area;
            let scale=box_.mass.volume/(box_.bbox.max.y-box_.bbox.min.y);
            let budget=(eps*scale.sqrt()*256.0).max(
                1024.0*f64::EPSILON*scale);
            if !area.is_finite()||area<=budget ||
                area>scale+budget{
                return Err(invalid("new floor material area is unprovable"));
            }
            // Connected three-cutter event: at least 3 active circles,
            // nonvirtual crossing arcs from 3 source circles and a source
            // contour from one connected void component.
            if upper.active_disks.len()>=3 &&
                upper.connected_union_components==1 {
                let mut observed=std::collections::BTreeSet::new();
                for lp in &upper.boundary_loops{
                    for arc in &lp.arcs{observed.insert(arc.disk_index);}
                }
                if observed.len()>=3 {found_connected_three=true;}
            }
            transitions.push(DepthTransition{
                height_y:height,
                activated_blind:activated,
                lower_slab:i,upper_slab:i+1,arcs,
                new_void_area:area,floor_material_area:area,
                component_count_before:lower.connected_union_components,
                component_count_after:upper.connected_union_components,
                planar_holes_before:lower.planar_holes,
                planar_holes_after:upper.planar_holes,
            });
        }
        if !found_connected_three{
            return Err(invalid("a connected three-circle exterior is not evidenced"));
        }
        let volume=plan.predicted_remaining_volume
            .ok_or_else(||invalid("source plan lacks predicted remaining volume"))?;
        Ok(ConnectedCutTransitionPlan{
            source_plan:plan,transitions,
            predicted_volume:volume,
            can_apply:false,
        })
    }
}
