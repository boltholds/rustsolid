use rustsolid::{
    BlindHoleParameters,HoleParameters,MixedCrossingSurfacePlanner,
    MixedFootprintOverlapPlanner,MixedFootprintContact,MixedPatchKind,
    Point3,Solid,MixedThroughBlindBrep,MixedThroughBlindParameters,
};
use std::f64::consts::{PI,TAU};

fn box_()->Solid{Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn close(a:f64,c:f64,e:f64){assert!((a-c).abs()<e,"got {a}; expected {c}");}
fn sample()->(HoleParameters,BlindHoleParameters){(t(0.,0.,0.5),b(0.5,0.,0.5,1.5))}

#[test]
fn real_analytic_support_partition_preserves_the_through_tunnel_genus(){
    let source=box_();let (through,blind)=sample();
    let plan=MixedCrossingSurfacePlanner::plan(&source,through,blind).unwrap();
    assert_eq!(plan.expected_genus,1);
    assert_eq!(plan.expected_euler_characteristic,0);
    assert_eq!(plan.top_boundary_arc_count,2);
    assert_eq!(plan.bottom_boundary_arc_count,1);
    close(plan.splitting_height_y,1.5,1e-12);
    assert!(!plan.can_apply());
    assert_eq!(plan.footprint.contact,MixedFootprintContact::TwoTransverseCrossings);
    plan.validate_against(&source,through,blind).unwrap();
    close(plan.predicted_remaining_volume,140.14094873829234,1e-10);
}

#[test]
fn support_areas_match_independent_sector_geometry_not_mesh(){
    let source=box_();let (through,blind)=sample();
    let plan=MixedCrossingSurfacePlanner::plan(&source,through,blind).unwrap();
    let overlap=plan.footprint.overlap_area.unwrap();
    let arc=plan.footprint.exposed_arcs.unwrap();
    let box_top_area=36.;
    let upper_removed=PI*(through.radius*through.radius+blind.radius*blind.radius)-overlap;
    let expected=[
        box_top_area-upper_removed,
        box_top_area-PI*through.radius.powi(2),
        96.,
        TAU*through.radius*(4.-blind.depth),
        arc[0].sweep_radians*through.radius*blind.depth,
        arc[1].sweep_radians*blind.radius*blind.depth,
        PI*blind.radius.powi(2)-overlap,
    ];
    for (part,value) in plan.patches.iter().zip(expected) {
        close(part.surface_area,value,1e-11);
        assert!(part.surface_area>0.);
    }
    assert_eq!(plan.patches[0].kind,MixedPatchKind::TopCapComplement);
    assert_eq!(plan.patches[1].kind,MixedPatchKind::BottomCapComplement);
    assert_eq!(plan.patches[3].kind,MixedPatchKind::ThroughLowerWall);
    assert_eq!(plan.patches[4].kind,MixedPatchKind::ThroughUpperWall);
    assert_eq!(plan.patches[5].kind,MixedPatchKind::BlindUpperWall);
    assert_eq!(plan.patches[6].kind,MixedPatchKind::BlindFloorCrescent);
    close(plan.predicted_surface_area,expected.iter().sum(),1e-9);
}

#[test]
fn through_upper_wall_disappears_where_blind_void_already_exists(){
    let source=box_();let through=t(0.,0.,1.);let blind=b(1.1,0.,0.8,2.);
    let plan=MixedCrossingSurfacePlanner::plan(&source,through,blind).unwrap();
    let exterior=plan.footprint.exposed_arcs.unwrap()[0].sweep_radians;
    let upper=plan.patches[4].surface_area;
    close(upper,through.radius*exterior*blind.depth,1e-10);
    assert!(upper<TAU*through.radius*blind.depth);
    // Physical exterior arcs belong to their own source support cylinders.
    let mid=plan.footprint.exposed_arcs.unwrap()[0].midpoint;
    assert!((mid.x-blind.center_x).hypot(mid.z-blind.center_z)>blind.radius);
}

#[test]
fn read_only_surface_plan_refuses_disjoint_tangent_and_contained_cases(){
    let source=box_();let through=t(0.,0.,0.5);
    for blind in [
        b(2.,0.,0.5,1.),
        b(1.,0.,0.5,1.),
        b(0.,0.,0.3,1.),
        b(0.,0.,0.5,1.),
    ]{
        assert!(MixedCrossingSurfacePlanner::plan(&source,through,blind).is_err());
    }
}

#[test]
fn tampering_cannot_create_a_valid_analytic_support_plan(){
    let source=box_();let (through,blind)=sample();
    let p=MixedCrossingSurfacePlanner::plan(&source,through,blind).unwrap();
    let mut q=p.clone();
    q.patches[4].surface_area*=0.75;
    assert!(q.validate_against(&source,through,blind).is_err());
    let mut q=p.clone();
    q.expected_genus=0;
    assert!(q.validate_against(&source,through,blind).is_err());
    let mut q=p.clone();
    q.splitting_height_y+=0.1;
    assert!(q.validate_against(&source,through,blind).is_err());
    assert!(p.validate_against(&source,through,b(0.6,0.,0.5,1.5)).is_err());
}

#[test]
fn old_mixed_boolean_remains_fail_closed_until_coedges_are_real(){
    let source=box_();let (through,blind)=sample();
    let original=(source.revision,source.mass.volume);
    let p=MixedCrossingSurfacePlanner::plan(&source,through,blind).unwrap();
    assert!(!p.can_apply());
    assert!(MixedThroughBlindBrep::from_box(&source,MixedThroughBlindParameters{
        through:vec![through],blind:vec![blind]
    }).is_err());
    assert_eq!(original,(source.revision,source.mass.volume));
}
