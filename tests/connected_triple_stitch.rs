//! Independent chained three-cylinder 3D halfedge stitching contracts.
//! No vendor code, tessellated area or proprietary CAD output is an oracle.
use rustsolid::{
    BlindHoleParameters,ConnectedTripleStitchPlanner,HoleParameters,
    Point3,Solid,TripleChainParameters,TripleStitchHold,TripleStitchState,
    TriplePhysicalEdgeKey,
};
use std::f64::consts::PI;
use std::collections::BTreeSet;
fn pos(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(pos(-4.,-2.,-3.),8.,5.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn b(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn fixture()->TripleChainParameters{
    TripleChainParameters{
        through:t(-1.15,0.,0.8),
        deeper_blind:b(0.,0.,0.8,2.8),
        shallower_blind:b(1.15,0.,0.8,1.2)
    }
}
fn approx(x:f64,y:f64,eps:f64){
    assert!((x-y).abs()<=eps,
        "observed {x:.14}, expected {y:.14}, error {:.3e}",(x-y).abs())
}
fn pair_lens(r1:f64,r2:f64,d:f64)->f64{
    let a=(d*d+r1*r1-r2*r2)/(2.*d);
    let other=d-a;
    let h=(r1*r1-a*a).sqrt();
    r1*r1*h.atan2(a)+r2*r2*h.atan2(other)-d*h
}
fn width(x:f64,disks:&[(f64,f64,f64)])->f64{
    let mut intervals=disks.iter().filter_map(|&(cx,cz,r)|{
        let v=r*r-(x-cx)*(x-cx);
        if v>0.0{
            let h=v.sqrt();
            Some((cz-h,cz+h))
        }else{None}
    }).collect::<Vec<_>>();
    intervals.sort_by(|a,b|a.0.total_cmp(&b.0));
    if intervals.is_empty(){return 0.0;}
    let (mut lo,mut hi)=intervals[0];
    let mut length=0.0;
    for (a,b) in intervals.into_iter().skip(1){
        if a>hi{
            length+=hi-lo;lo=a;hi=b;
        }else{hi=hi.max(b);}
    }
    length+hi-lo
}
/// Numerical interval-union quadrature is mathematically independent of
/// circle Green integrals, seam orientation and event-graph decomposition.
fn oracle_area(disks:&[(f64,f64,f64)])->f64{
    let a=disks.iter().map(|d|d.0-d.2).fold(f64::INFINITY,f64::min);
    let b=disks.iter().map(|d|d.0+d.2).fold(f64::NEG_INFINITY,f64::max);
    let n=32768;
    let delta=PI/n as f64;
    let mut sum=0.0;
    for i in 0..=n{
        let th=i as f64*delta;
        let x=a+0.5*(b-a)*(1.0-th.cos());
        let jac=0.5*(b-a)*th.sin();
        let f=width(x,disks)*jac;
        let w=if i==0||i==n{1.0}else if i%2==0{2.0}else{4.0};
        sum+=w*f;
    }
    sum*delta/3.0
}
#[test]
fn genuinely_connected_three_cutters_stitch_closed_source_atom_contours(){
    let source=block();let input=fixture();
    let plan=ConnectedTripleStitchPlanner::plan(&source,input).unwrap();
    assert_eq!(plan.state,TripleStitchState::Ready);
    assert!(!plan.can_apply());
    assert_eq!(plan.levels.len(),4);
    for (&level,&expected) in plan.levels.iter().zip(&[-2.0,0.2,1.8,3.0]) {
        approx(level,expected,1e-12);
    }
    assert_eq!(plan.junctions.len(),4);
    assert_eq!(plan.atomic_arcs.len(),8);
    assert_eq!(plan.slabs.len(),3);
    assert_eq!(plan.floors.len(),2);
    assert_eq!(plan.slabs[0].active_sources,[0]);
    assert_eq!(plan.slabs[1].active_sources,[0,1]);
    assert_eq!(plan.slabs[2].active_sources,[0,1,2]);
    assert_eq!(plan.bottom_cap.as_ref().unwrap().parts.len(),2);
    assert_eq!(plan.top_cap.as_ref().unwrap().parts.len(),4);
    assert_eq!(plan.floors[0].contour.parts.len(),4);
    assert_eq!(plan.floors[1].contour.parts.len(),2);
    assert_eq!(plan.predicted_genus,Some(1));
    assert_eq!(plan.predicted_euler,Some(0));
    assert_eq!(plan.predicted_shell_vertices,Some(22));
    assert_eq!(plan.predicted_shell_edges,Some(38));
    assert_eq!(plan.predicted_shell_faces,Some(18));
    assert_eq!(plan.walls.len(),10);
    assert_eq!(plan.shared_physical_edges.len(),26);
    assert!(plan.shared_physical_edges.iter().all(|e|{
        e.opposite_use_signs[0]==-e.opposite_use_signs[1]
    }));
    plan.validate_against(&source,input).unwrap();
}
#[test]
fn two_real_depth_floor_interfaces_are_connected_crescents_with_source_arcs(){
    let plan=ConnectedTripleStitchPlanner::plan(&block(),fixture()).unwrap();
    assert_eq!(plan.state,TripleStitchState::Ready);
    let r=0.8;let d=1.15;let circle=PI*r*r;
    let lens=pair_lens(r,r,d);
    let expected_floor=circle-lens;
    for floor in &plan.floors{
        approx(floor.floor_area,expected_floor,1e-10);
        approx(floor.contour.signed_area,expected_floor,1e-10);
        assert!(!floor.incoming_lower_arcs.is_empty());
        assert!(!floor.born_upper_arcs.is_empty());
    }
    // At first floor, B is not yet active and its intersections only
    // subdivide existing source A circular arcs; the physical crescent
    // still has a valid closed contour.
    let deep_sources=plan.floors[0].contour.parts.iter()
        .map(|a|plan.atomic_arcs[a.arc].source).collect::<BTreeSet<_>>();
    let shallow_sources=plan.floors[1].contour.parts.iter()
        .map(|a|plan.atomic_arcs[a.arc].source).collect::<BTreeSet<_>>();
    assert_eq!(deep_sources,BTreeSet::from([0,1]));
    assert_eq!(shallow_sources,BTreeSet::from([1,2]));
}
#[test]
fn analytical_slab_volume_matches_independent_pair_lenses_and_strip_union(){
    let plan=ConnectedTripleStitchPlanner::plan(&block(),fixture()).unwrap();
    let r=0.8;let lens=pair_lens(r,r,1.15);
    let disk=PI*r*r;
    let expected=240.0-disk*5.0-(disk-lens)*2.8-(disk-lens)*1.2;
    approx(plan.predicted_remaining_volume.unwrap(),expected,1e-10);
    let disks=[
        (-1.15,0.,r),(0.,0.,r),(1.15,0.,r)
    ];
    let mut removed=0.0;
    for (idx,slab) in plan.slabs.iter().enumerate(){
        let numeric=oracle_area(&disks[..idx+1]);
        approx(slab.union_area,numeric,5e-5);
        removed+=numeric*(slab.y_upper-slab.y_lower);
    }
    approx(plan.predicted_remaining_volume.unwrap(),240.0-removed,2e-4);
    assert!(plan.analytic_surface_area.unwrap()>block().mass.surface_area-20.0);
}
#[test]
fn source_supported_uv_wall_patches_share_pairwise_physical_junctions(){
    let plan=ConnectedTripleStitchPlanner::plan(&block(),fixture()).unwrap();
    assert_eq!(plan.state,TripleStitchState::Ready);
    let cut=fixture();
    let sources=[
        (cut.through.center_x,cut.through.center_z,cut.through.radius),
        (cut.deeper_blind.center_x,cut.deeper_blind.center_z,cut.deeper_blind.radius),
        (cut.shallower_blind.center_x,cut.shallower_blind.center_z,cut.shallower_blind.radius),
    ];
    for w in &plan.walls{
        let (cx,cz,r)=sources[w.source];
        let [u0,u1,v0,v1]=w.uv_bounds;
        approx(u1-u0,plan.atomic_arcs[w.arc].arc_length/r,1e-10);
        approx(v1-v0,w.y_upper-w.y_lower,1e-12);
        assert!(w.analytic_area>0.0);
        for (idx,p) in w.corners.iter().enumerate(){
            let u=if idx%2==0{u0}else{u1};
            let y=if idx<2{w.y_lower}else{w.y_upper};
            approx(p.x,cx+r*u.cos(),1e-12);
            approx(p.z,cz+r*u.sin(),1e-12);
            approx(p.y,y,1e-12);
            approx((p.x-cx).hypot(p.z-cz),r,1e-12);
        }
    }
    let n_circle=plan.shared_physical_edges.iter()
        .filter(|e|matches!(e.key,TriplePhysicalEdgeKey::CircleArc{..})).count();
    let n_vertical=plan.shared_physical_edges.iter()
        .filter(|e|matches!(e.key,TriplePhysicalEdgeKey::VerticalTrace{..})).count();
    assert_eq!(n_circle,16);
    assert_eq!(n_vertical,10);
}
#[test]
fn rotated_and_far_translated_chains_preserve_genus_area_and_volume(){
    let first=ConnectedTripleStitchPlanner::plan(&block(),fixture()).unwrap();
    let r=0.8;
    let rotated=TripleChainParameters{
        through:t(0.,-1.15,r),
        deeper_blind:b(0.,0.,r,2.8),
        shallower_blind:b(0.,1.15,r,1.2),
    };
    let rot=ConnectedTripleStitchPlanner::plan(&block(),rotated).unwrap();
    assert_eq!(rot.state,TripleStitchState::Ready);
    approx(rot.predicted_remaining_volume.unwrap(),
        first.predicted_remaining_volume.unwrap(),1e-9);
    approx(rot.analytic_surface_area.unwrap(),
        first.analytic_surface_area.unwrap(),1e-9);
    let dx=10000.;let dz=-10000.;
    let moved=Solid::block(pos(-4.+dx,-2.,-3.+dz),8.,5.,6.).unwrap();
    let i=fixture();
    let shifted=TripleChainParameters{
        through:t(i.through.center_x+dx,i.through.center_z+dz,r),
        deeper_blind:b(i.deeper_blind.center_x+dx,i.deeper_blind.center_z+dz,r,2.8),
        shallower_blind:b(i.shallower_blind.center_x+dx,
            i.shallower_blind.center_z+dz,r,1.2),
    };
    let plan=ConnectedTripleStitchPlanner::plan(&moved,shifted).unwrap();
    assert_eq!(plan.state,TripleStitchState::Ready);
    approx(plan.predicted_remaining_volume.unwrap(),
        first.predicted_remaining_volume.unwrap(),1e-8);
}
#[test]
fn unsupported_near_tangency_interacting_through_shallow_and_depth_order_hold(){
    let src=block();
    let mut p=fixture();
    p.shallower_blind.center_x=1.6;
    let held=ConnectedTripleStitchPlanner::plan(&src,p).unwrap();
    assert_eq!(held.state,TripleStitchState::Held(TripleStitchHold::NearTangency));
    assert!(held.predicted_remaining_volume.is_none());
    assert!(held.atomic_arcs.is_empty());
    assert!(!held.can_apply());
    let mut p=fixture();
    p.shallower_blind.center_x=-0.4; // third cutter contacts through too
    let held=ConnectedTripleStitchPlanner::plan(&src,p).unwrap();
    assert_eq!(held.state,TripleStitchState::Held(
        TripleStitchHold::UnsupportedCircleContact));
    let mut p=fixture();
    p.shallower_blind.depth=p.deeper_blind.depth;
    let held=ConnectedTripleStitchPlanner::plan(&src,p).unwrap();
    assert_eq!(held.state,TripleStitchState::Held(
        TripleStitchHold::UnsupportedCircleContact));
}
#[test]
fn source_tampering_is_rejected_and_no_input_topology_is_modified(){
    let src=block();let i=fixture();
    let baseline=(src.revision,src.mass,src.faces.len(),src.edges.len());
    let plan=ConnectedTripleStitchPlanner::plan(&src,i).unwrap();
    plan.validate_against(&src,i).unwrap();
    let mut broken=plan.clone();
    broken.walls[0].corners[0].x+=0.05;
    assert!(broken.validate_against(&src,i).is_err());
    let mut broken=plan.clone();
    broken.shared_physical_edges[0].opposite_use_signs=[1,1];
    assert!(broken.validate_against(&src,i).is_err());
    let mut broken=plan.clone();
    broken.floors[1].contour.signed_area+=0.5;
    assert!(broken.validate_against(&src,i).is_err());
    let mut altered=i;
    altered.shallower_blind.depth=1.4;
    assert!(plan.validate_against(&src,altered).is_err());
    assert_eq!(baseline,(src.revision,src.mass,src.faces.len(),src.edges.len()));
}
#[test]
fn invalid_nonbox_bad_radii_and_excessive_depth_refuse_without_side_effects(){
    let s=block();
    let mut i=fixture();
    i.shallower_blind.radius=-1.0;
    assert!(ConnectedTripleStitchPlanner::plan(&s,i).is_err());
    let mut i=fixture();
    i.shallower_blind.depth=5.0;
    assert!(ConnectedTripleStitchPlanner::plan(&s,i).is_err());
    let mut i=fixture();
    i.shallower_blind.center_x=f64::NAN;
    assert!(ConnectedTripleStitchPlanner::plan(&s,i).is_err());
    let no_box=Solid::extrude_xz(&[
        rustsolid::Point2{x:-4.,z:-3.},rustsolid::Point2{x:4.,z:-3.},
        rustsolid::Point2{x:3.,z:3.},rustsolid::Point2{x:-4.,z:3.}
    ],5.).unwrap();
    assert!(ConnectedTripleStitchPlanner::plan(&no_box,fixture()).is_err());
}
