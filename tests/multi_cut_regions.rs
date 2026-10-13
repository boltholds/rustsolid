//! First-party multi-circle union oracle and depth-sweep regression corpus.
//! The oracle samples UNION OF VERTICAL INTERVALS, not exposed circle arcs.
//! The active-source arrangement may have multiple components and true holes.
//! None of these plans are executable B-rep or faceted approximation.
use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters,HoleParameters,MultiCutRegionPlanner,MultiCutPlanState,
    MultiCutHold,MixedCrossingBrep,Point2,Point3,Solid,
};
fn point(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(point(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn b(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn close(a:f64,e:f64,tol:f64){
    assert!((a-e).abs()<=tol,
        "got {a:.14}, expected {e:.14}; difference {:.4e}",(a-e).abs());
}

/// Evaluate union of independently derived vertical circle slice intervals.
fn interval_width(x:f64,disks:&[(f64,f64,f64)])->f64{
    let mut spans=disks.iter().filter_map(|&(cx,cz,r)|{
        let squared=r*r-(x-cx).powi(2);
        if squared>0.0{
            let d=squared.sqrt();
            Some((cz-d,cz+d))
        }else{None}
    }).collect::<Vec<_>>();
    if spans.is_empty(){return 0.0;}
    spans.sort_by(|a,b|a.0.total_cmp(&b.0));
    let mut lower=spans[0].0;
    let mut upper=spans[0].1;
    let mut sum=0.0;
    for &(start,end) in spans.iter().skip(1){
        if start>upper{
            sum+=upper-lower;
            lower=start;upper=end;
        }else{upper=upper.max(end);}
    }
    sum+upper-lower
}

/// Cosine-transformed Simpson integral; its algorithm and numerical
/// error mechanism are independent of the production Green arc integral.
fn numeric_union_area(disks:&[(f64,f64,f64)],n:usize)->f64{
    assert!(n%2==0 && n>=1024);
    let left=disks.iter().map(|&(x,_,r)|x-r).fold(f64::INFINITY,f64::min);
    let right=disks.iter().map(|&(x,_,r)|x+r).fold(f64::NEG_INFINITY,f64::max);
    let step=PI/n as f64;
    let mut sum=0.0;
    for i in 0..=n{
        let theta=i as f64*step;
        let x=left+(right-left)*0.5*(1.0-theta.cos());
        let dx=0.5*(right-left)*theta.sin();
        let integrand=interval_width(x,disks)*dx;
        let weight=if i==0||i==n{1.0}else if i%2==0{2.0}else{4.0};
        sum+=weight*integrand;
    }
    sum*step/3.0
}
fn sources(through:&[HoleParameters],blind:&[BlindHoleParameters])->Vec<(f64,f64,f64)>{
    through.iter().map(|c|(c.center_x,c.center_z,c.radius))
        .chain(blind.iter().map(|c|(c.center_x,c.center_z,c.radius)))
        .collect()
}

#[test]
fn one_crossing_pair_matches_existing_v031_executable_analytic_mass(){
    let source=block();
    let through=[t(0.,0.,0.5)];
    let blind=[b(0.5,0.,0.5,1.5)];
    let plan=MultiCutRegionPlanner::plan(&source,&through,&blind).unwrap();
    assert_eq!(plan.state,MultiCutPlanState::Ready);
    assert!(!plan.can_apply());
    assert_eq!(plan.slabs.len(),2);
    assert_eq!(plan.slabs[0].active_disks,[0]);
    assert_eq!(plan.slabs[1].active_disks,[0,1]);
    let actual=MixedCrossingBrep::from_box(&source,through[0],blind[0]).unwrap();
    close(plan.predicted_remaining_volume.unwrap(),actual.model().mass.volume,1e-10);
    close(plan.predicted_remaining_volume.unwrap(),140.14094873829234,1e-10);
    for slab in &plan.slabs{
        assert_eq!(slab.connected_union_components,1);
        assert_eq!(slab.planar_holes,0);
        assert_eq!(slab.boundary_loops.len(),1);
        close(slab.union_area,
              numeric_union_area(&sources(&through,&blind)
                .into_iter().take(slab.active_disks.len()).collect::<Vec<_>>(),16384),
              3e-5);
    }
    plan.validate_against(&source,&through,&blind).unwrap();
}

#[test]
fn three_interacting_cutters_with_two_blind_depths_are_certified_by_slices(){
    let source=block();
    let through=[t(0.,0.,0.75)];
    let blind=[b(0.7,0.,0.75,2.0),b(1.4,0.25,0.6,1.0)];
    let disks=sources(&through,&blind);
    let plan=MultiCutRegionPlanner::plan(&source,&through,&blind).unwrap();
    assert_eq!(plan.state,MultiCutPlanState::Ready);
    assert_eq!(plan.slabs.len(),3);
    assert_eq!(plan.slabs[0].active_disks,[0]);
    assert_eq!(plan.slabs[1].active_disks,[0,1]);
    assert_eq!(plan.slabs[2].active_disks,[0,1,2]);
    let mut numeric_removed=0.0;
    for slab in &plan.slabs{
        assert!(slab.y_upper>slab.y_lower);
        let active=slab.active_disks.iter().map(|&i|disks[i])
            .collect::<Vec<_>>();
        let independent_area=numeric_union_area(&active,32768);
        close(slab.union_area,independent_area,4e-5);
        close(slab.removed_volume,
            slab.union_area*(slab.y_upper-slab.y_lower),1e-11);
        assert_eq!(slab.connected_union_components,1);
        assert_eq!(slab.planar_holes,0);
        for lp in &slab.boundary_loops{
            assert!(lp.signed_area>0.0);
            for a in &lp.arcs{
                assert!(a.end_angle>a.start_angle);
                let c=disks[a.disk_index];
                close((a.midpoint.x-c.0).hypot(a.midpoint.z-c.1),c.2,1e-10);
            }
        }
        numeric_removed+=independent_area*(slab.y_upper-slab.y_lower);
    }
    close(plan.predicted_remaining_volume.unwrap(),
        144.0-numeric_removed,1e-4);
    assert!(!plan.can_apply());
    plan.validate_against(&source,&through,&blind).unwrap();
}

#[test]
fn three_pairwise_overlapping_disks_create_exposed_multi_source_arcs(){
    let through=[t(0.,0.,0.7)];
    let blind=[b(0.75,0.,0.7,1.8),b(0.35,0.65,0.7,1.8)];
    let plan=MultiCutRegionPlanner::plan(&block(),&through,&blind).unwrap();
    assert_eq!(plan.state,MultiCutPlanState::Ready);
    let top=plan.slabs.last().unwrap();
    assert_eq!(top.active_disks,[0,1,2]);
    let represented=top.boundary_loops.iter().flat_map(|lp|lp.arcs.iter())
        .map(|arc|arc.disk_index).collect::<std::collections::BTreeSet<_>>();
    assert_eq!(represented.len(),3);
    assert_eq!(top.connected_union_components,1);
    let disks=sources(&through,&blind);
    close(top.union_area,numeric_union_area(&disks,32768),4e-5);
}

#[test]
fn ring_of_three_disks_detects_true_undrilled_planar_hole(){
    // Side length about 1.8, disk radius 0.95: each pair intersects,
    // but the central circumcircle radius exceeds the disk radius.
    let through=[t(-0.9,-0.52,0.95)];
    let blind=[b(0.9,-0.52,0.95,2.),b(0.,1.04,0.95,2.)];
    let plan=MultiCutRegionPlanner::plan(&block(),&through,&blind).unwrap();
    assert_eq!(plan.state,MultiCutPlanState::Ready);
    let top=plan.slabs.last().unwrap();
    assert_eq!(top.connected_union_components,1);
    assert_eq!(top.planar_holes,1);
    assert_eq!(top.boundary_loops.len(),2);
    assert_eq!(top.boundary_loops.iter().filter(|lp|lp.is_hole).count(),1);
    assert!(top.boundary_loops.iter().any(|lp|lp.signed_area<0.0));
    close(top.union_area,
        numeric_union_area(&sources(&through,&blind),32768),5e-5);
}

#[test]
fn nested_and_disjoint_cutters_never_invent_intersection_vertices(){
    let through=[t(-1.3,0.,0.8)];
    let blind=[b(-1.3,0.,0.2,2.),b(1.35,0.1,0.35,1.)];
    let plan=MultiCutRegionPlanner::plan(&block(),&through,&blind).unwrap();
    assert_eq!(plan.state,MultiCutPlanState::Ready);
    let top=plan.slabs.last().unwrap();
    assert!(top.vertices.iter().all(|v|v.is_virtual_seam));
    assert_eq!(top.connected_union_components,2);
    assert_eq!(top.planar_holes,0);
    let area=PI*(0.8f64.powi(2)+0.35f64.powi(2));
    close(top.union_area,area,1e-9);
    assert!(top.boundary_loops.iter().all(|lp|lp.arcs.len()==1));
}

#[test]
fn unsafe_external_tangency_and_triple_junction_are_held_without_mass(){
    let source=block();
    for offset in [1.4,1.4-1e-10,1.4+1e-10]{
        let p=MultiCutRegionPlanner::plan(
            &source,&[t(0.,0.,0.7)],&[b(offset,0.,0.7,1.5)]
        ).unwrap();
        assert_eq!(p.state,MultiCutPlanState::Held(
            MultiCutHold::NearTangencyOrCoincidence));
        assert!(p.slabs.is_empty()&&p.predicted_remaining_volume.is_none());
        assert!(!p.can_apply());
    }
    let z=3.0_f64.sqrt()*0.5-1.0;
    let triple=MultiCutRegionPlanner::plan(
        &source,&[t(0.,0.,1.)],&[b(1.,0.,1.,1.5),b(0.5,z,1.,1.5)]
    ).unwrap();
    assert_eq!(triple.state,MultiCutPlanState::Held(
        MultiCutHold::NearTripleJunction));
    assert!(triple.slabs.is_empty()&&triple.removed_volume.is_none());
}

#[test]
fn near_equal_depth_events_hold_but_exact_equal_depths_are_valid(){
    let source=block();
    let through=[t(0.,0.,0.5)];
    let cuts=[b(0.5,0.,0.5,1.2),b(-1.5,1.0,0.4,1.2+1e-10)];
    let held=MultiCutRegionPlanner::plan(&source,&through,&cuts).unwrap();
    assert_eq!(held.state,MultiCutPlanState::Held(
        MultiCutHold::UnresolvedDepthEvent));
    let equal=[cuts[0],b(-1.5,1.0,0.4,1.2)];
    let p=MultiCutRegionPlanner::plan(&source,&through,&equal).unwrap();
    assert_eq!(p.state,MultiCutPlanState::Ready);
    assert_eq!(p.slabs.len(),2);
}

#[test]
fn depth_monotonicity_permutation_and_large_translation_preserve_source_mass(){
    let through=[t(0.,0.,0.75)];
    let a=[b(0.7,0.,0.75,2.0),b(1.4,0.25,0.6,1.0)];
    let p=MultiCutRegionPlanner::plan(&block(),&through,&a).unwrap();
    let swap=[a[1],a[0]];
    let q=MultiCutRegionPlanner::plan(&block(),&through,&swap).unwrap();
    assert_eq!(p.state,MultiCutPlanState::Ready);
    assert_eq!(q.state,MultiCutPlanState::Ready);
    close(p.predicted_remaining_volume.unwrap(),
        q.predicted_remaining_volume.unwrap(),1e-10);
    let deeper=[a[0],b(a[1].center_x,a[1].center_z,a[1].radius,1.4)];
    let r=MultiCutRegionPlanner::plan(&block(),&through,&deeper).unwrap();
    assert_eq!(r.state,MultiCutPlanState::Ready);
    assert!(r.predicted_remaining_volume.unwrap()<
        p.predicted_remaining_volume.unwrap());
    let offset=10000.0;
    let shifted=Solid::block(point(-3.+offset,-1.,-3.-offset),6.,4.,6.).unwrap();
    let t_shift=[HoleParameters{center_x:offset,center_z:-offset,radius:through[0].radius}];
    let b_shift=a.map(|cut|b(cut.center_x+offset,cut.center_z-offset,cut.radius,cut.depth));
    let translated=MultiCutRegionPlanner::plan(&shifted,&t_shift,&b_shift).unwrap();
    assert_eq!(translated.state,MultiCutPlanState::Ready);
    close(translated.predicted_remaining_volume.unwrap(),
        p.predicted_remaining_volume.unwrap(),1e-8);
}

#[test]
fn source_replay_rejects_tampering_and_never_modifies_topology(){
    let src=block();
    let t=[t(0.,0.,0.75)];
    let blinds=[b(0.7,0.,0.75,2.),b(1.4,0.25,0.6,1.)];
    let before=(src.revision,src.mass,src.edges.len(),src.faces.len());
    let plan=MultiCutRegionPlanner::plan(&src,&t,&blinds).unwrap();
    plan.validate_against(&src,&t,&blinds).unwrap();
    let mut tampered=plan.clone();
    tampered.slabs[2].union_area+=0.01;
    assert!(tampered.validate_against(&src,&t,&blinds).is_err());
    let mut tampered=plan.clone();
    tampered.slabs[2].boundary_loops[0].arcs[0].midpoint.x+=0.01;
    assert!(tampered.validate_against(&src,&t,&blinds).is_err());
    assert!(plan.validate_against(&src,&t,&[blinds[0],b(1.4,0.25,0.6,1.2)]).is_err());
    assert_eq!((src.revision,src.mass,src.edges.len(),src.faces.len()),before);
    assert!(!plan.can_apply());
}

#[test]
fn rejects_nonbox_bad_radii_missing_families_and_unbounded_work(){
    let src=block();
    let not_box=Solid::extrude_xz(&[
        Point2{x:-3.,z:-3.},Point2{x:3.,z:-3.},
        Point2{x:2.,z:3.},Point2{x:-3.,z:3.}],4.).unwrap();
    assert!(MultiCutRegionPlanner::plan(&not_box,
        &[t(0.,0.,0.5)],&[b(0.5,0.,0.5,2.)]).is_err());
    assert!(MultiCutRegionPlanner::plan(&src,&[],&[b(0.5,0.,0.5,2.)]).is_err());
    assert!(MultiCutRegionPlanner::plan(&src,&[t(0.,0.,0.5)],&[]).is_err());
    for bad in [
        t(2.9,0.,0.5),t(0.,0.,0.),t(0.,0.,f64::NAN),
        t(0.,0.,-0.1),
    ]{
        assert!(MultiCutRegionPlanner::plan(&src,&[bad],&[b(0.5,0.,0.5,2.)]).is_err());
    }
    for bad in [
        b(0.5,0.,0.5,0.),b(0.5,0.,0.5,4.),b(0.5,0.,0.5,-1.),
        b(0.5,0.,f64::INFINITY,2.),b(2.8,0.,0.5,2.),
    ]{
        assert!(MultiCutRegionPlanner::plan(&src,&[t(0.,0.,0.5)],&[bad]).is_err());
    }
    let many=(0..12).map(|i|t(-2.5+i as f64*0.4,0.,0.1))
        .collect::<Vec<_>>();
    assert!(MultiCutRegionPlanner::plan(&src,&many,&[b(0.5,0.,0.5,2.)]).is_err());
}
