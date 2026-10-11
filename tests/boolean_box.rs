//! Independent Boolean acceptance: mass, symmetry, genus boundaries, lineage,
//! fresh identities, contact topology, near-coincidence and read-only carriers.
use rustsolid::{
    BooleanError, BooleanKernel, BooleanOperand, BooleanOperation, BooleanOutcome,
    BrepArena, BrepBody, CylindricalBrep, GeometryTolerance, Point3, Point2,
    Solid, TopologyEntity,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_(x:f64,y:f64,z:f64,w:f64,h:f64,d:f64)->Solid{
    Solid::block(p(x,y,z),w,h,d).unwrap()
}
fn a()->Solid{box_(0.,0.,0.,2.,3.,4.)}
fn b()->Solid{box_(1.,1.,1.,2.,3.,4.)}
fn run(a:&Solid,b:&Solid,op:BooleanOperation)->rustsolid::BooleanResult{
    BooleanKernel::apply_solids(a,b,op).unwrap()
}
fn expect_solid(result:&rustsolid::BooleanResult,volume:f64,area:f64,centroid:[f64;3]){
    let solid=result.solid().expect("expected one closed result shell");
    solid.validate().unwrap();
    assert!((solid.mass.volume-volume).abs()<1e-8,
        "unexpected volume: {}",solid.mass.volume);
    assert!((solid.mass.surface_area-area).abs()<1e-8,
        "unexpected area: {}",solid.mass.surface_area);
    for (u,v) in solid.mass.centroid.array().into_iter().zip(centroid) {
        assert!((u-v).abs()<1e-8,"unexpected centroid {u} != {v}");
    }
    assert_eq!(solid.euler_characteristic(),2);
    assert_eq!(solid.shells.len(),1);
    assert_eq!(result.face_sources.len(),solid.faces.len());
}

#[test]
fn intersecting_boxes_are_exact_and_preserve_volume_inclusion_exclusion(){
    let (a,b)=(a(),b());
    let ab=run(&a,&b,BooleanOperation::Intersection);
    expect_solid(&ab,6.,22.,[1.5,2.,2.5]);
    let union=run(&a,&b,BooleanOperation::Union);
    expect_solid(&union,42.,82.,[1.5,2.,2.5]);
    let difference=run(&a,&b,BooleanOperation::Difference);
    expect_solid(&difference,18.,52.,[15./18.,24./18.,33./18.]);
    assert!((union.solid().unwrap().mass.volume -
        (a.mass.volume+b.mass.volume-ab.solid().unwrap().mass.volume)).abs()<1e-9);
    assert!((difference.solid().unwrap().mass.volume +
        ab.solid().unwrap().mass.volume-a.mass.volume).abs()<1e-9);
    assert!(difference.face_sources.iter().any(|f|
        f.operand==BooleanOperand::Right && f.reversed));
    for (result,left,right) in [(&ab,&a,&b),(&union,&a,&b),(&difference,&a,&b)] {
        for face in &result.face_sources {
            let source=match face.operand{BooleanOperand::Left=>left,BooleanOperand::Right=>right};
            assert!(source.faces.get(face.source_face.0 as usize).is_some());
        }
    }
    // Read-only operands are never modified by Boolean evaluation.
    a.validate().unwrap();
    b.validate().unwrap();
    assert_eq!((a.revision,b.revision),(0,0));
}

#[test]
fn union_intersection_commute_but_difference_has_direction(){
    let (a,b)=(a(),b());
    for op in [BooleanOperation::Intersection,BooleanOperation::Union]{
        let one=run(&a,&b,op);let other=run(&b,&a,op);
        let x=one.solid().unwrap();let y=other.solid().unwrap();
        assert!((x.mass.volume-y.mass.volume).abs()<1e-9);
        assert!((x.mass.surface_area-y.mass.surface_area).abs()<1e-9);
        assert_eq!(x.bbox,y.bbox);
    }
    let ab=run(&a,&b,BooleanOperation::Difference);
    let ba=run(&b,&a,BooleanOperation::Difference);
    assert_eq!(ab.solid().unwrap().mass.volume,ba.solid().unwrap().mass.volume);
    assert_ne!(ab.solid().unwrap().bbox,ba.solid().unwrap().bbox);
}

#[test]
fn exact_coincident_boxes_and_whole_subtraction(){
    let a=a();let b=a.clone();
    expect_solid(&run(&a,&b,BooleanOperation::Union),24.,52.,[1.,1.5,2.]);
    expect_solid(&run(&a,&b,BooleanOperation::Intersection),24.,52.,[1.,1.5,2.]);
    let d=run(&a,&b,BooleanOperation::Difference);
    assert!(matches!(d.outcome,BooleanOutcome::Empty));
    assert!(d.face_sources.is_empty());
}

#[test]
fn disjoint_intersection_is_empty_and_single_shell_union_fails(){
    let a=a();let b=box_(10.,10.,10.,1.,1.,1.);
    assert!(matches!(run(&a,&b,BooleanOperation::Intersection).outcome,
        BooleanOutcome::Empty));
    expect_solid(&run(&a,&b,BooleanOperation::Difference),24.,52.,[1.,1.5,2.]);
    assert!(matches!(
        BooleanKernel::apply_solids(&a,&b,BooleanOperation::Union),
        Err(BooleanError::UnsupportedResultTopology(_))
    ));
}

#[test]
fn face_touch_union_fuses_and_touch_intersection_is_empty(){
    let a=box_(0.,0.,0.,1.,1.,1.);
    let b=box_(1.,0.,0.,1.,1.,1.);
    expect_solid(&run(&a,&b,BooleanOperation::Union),2.,10.,[1.,0.5,0.5]);
    expect_solid(&run(&a,&b,BooleanOperation::Difference),1.,6.,[0.5,0.5,0.5]);
    assert!(matches!(run(&a,&b,BooleanOperation::Intersection).outcome,
        BooleanOutcome::Empty));
}

#[test]
fn edge_or_point_contact_cannot_fake_a_manifold() {
    let a=box_(0.,0.,0.,1.,1.,1.);
    for b in [box_(1.,1.,0.,1.,1.,1.),box_(1.,1.,1.,1.,1.,1.)] {
        assert!(matches!(
            BooleanKernel::apply_solids(&a,&b,BooleanOperation::Union),
            Err(BooleanError::UnsupportedResultTopology(_))
        ));
    }
}

#[test]
fn cavity_and_split_components_are_explicitly_unsupported(){
    let large=box_(0.,0.,0.,3.,3.,3.);
    let inner=box_(1.,1.,1.,1.,1.,1.);
    let through=box_(1.,0.,0.,1.,3.,3.);
    for b in [&inner,&through]{
        assert!(matches!(
            BooleanKernel::apply_solids(&large,b,BooleanOperation::Difference),
            Err(BooleanError::UnsupportedResultTopology(_))
        ));
    }
    expect_solid(&run(&large,&inner,BooleanOperation::Intersection),
        1.,6.,[1.5,1.5,1.5]);
}

#[test]
fn near_coincident_boundaries_are_not_silently_snapped(){
    let a=box_(0.,0.,0.,1.,1.,1.);
    let b=box_(1.0+5e-10,0.,0.,1.,1.,1.);
    assert!(matches!(
        BooleanKernel::apply_solids(&a,&b,BooleanOperation::Union),
        Err(BooleanError::AmbiguousBoundary)
    ));
}

#[test]
fn nonbox_or_analytic_source_cannot_be_replaced_by_its_bbox(){
    let solid=a();
    let concave=Solid::extrude_xz(&[
        Point2{x:0.,z:0.},Point2{x:3.,z:0.},Point2{x:3.,z:1.},
        Point2{x:1.,z:1.},Point2{x:1.,z:3.},Point2{x:0.,z:3.}
    ],2.).unwrap();
    assert!(matches!(
        BooleanKernel::apply_solids(&solid,&concave,BooleanOperation::Union),
        Err(BooleanError::UnsupportedOperand(_))
    ));
    let cylinder=CylindricalBrep::upright(
        p(0.,0.,0.),1.,2.,GeometryTolerance::default()).unwrap();
    assert!(matches!(
        BooleanKernel::apply(&BrepBody::Polyhedral(solid),
            &BrepBody::Cylinder(cylinder),BooleanOperation::Difference),
        Err(BooleanError::UnsupportedOperand(_))
    ));
}

#[test]
fn arena_boolean_is_nondestructive_and_returns_new_identity_incarnation(){
    let left=BrepArena::from_solid(a()).unwrap();
    let right=BrepArena::from_solid(b()).unwrap();
    let first=left.snapshot().unwrap();let second=right.snapshot().unwrap();
    let output=left.boolean(&right,BooleanOperation::Intersection).unwrap();
    expect_solid(&output,6.,22.,[1.5,2.,2.5]);
    let output_arena=BrepArena::from_solid(output.into_solid().unwrap()).unwrap();
    assert_ne!(left.arena_token(),output_arena.arena_token());
    assert_ne!(right.arena_token(),output_arena.arena_token());
    assert_ne!(left.handle(TopologyEntity::Body),output_arena.handle(TopologyEntity::Body));
    left.check_snapshot(&first).unwrap();
    right.check_snapshot(&second).unwrap();
    assert_eq!(left.undo_depth(),0);
    assert_eq!(right.undo_depth(),0);
    output_arena.validate().unwrap();
}

#[test]
fn lattice_translation_metamorphic_volume_contracts(){
    // The overlap volume is an independent closed-form oracle, not another
    // CAD implementation. All 27 intersections are strictly three-dimensional.
    let origin=p(0.,0.,0.);
    let a=Solid::block(origin,3.,3.,3.).unwrap();
    for dx in [-2.0,0.0,2.0] {
        for dy in [-2.0,0.0,2.0] {
            for dz in [-2.0,0.0,2.0] {
                let b=box_(dx,dy,dz,3.,3.,3.);
                let intersection=run(&a,&b,BooleanOperation::Intersection);
                let union=run(&a,&b,BooleanOperation::Union);
                let difference=run(&a,&b,BooleanOperation::Difference);
                let exact_overlap=(3.0-f64::abs(dx))*
                    (3.0-f64::abs(dy))*(3.0-f64::abs(dz));
                assert!((intersection.solid().unwrap().mass.volume-exact_overlap).abs()<1e-8);
                assert!((union.solid().unwrap().mass.volume-(54.0-exact_overlap)).abs()<1e-8);
                if exact_overlap==27.0 {
                    assert!(matches!(difference.outcome,BooleanOutcome::Empty));
                } else {
                    assert!((difference.solid().unwrap().mass.volume -
                        (27.0-exact_overlap)).abs()<1e-8);
                    difference.solid().unwrap().validate().unwrap();
                }
                union.solid().unwrap().validate().unwrap();
                intersection.solid().unwrap().validate().unwrap();
            }
        }
    }
}

#[test]
fn large_world_coordinates_do_not_destroy_local_mass(){
    let a=box_(1.0e9,1.0e9,1.0e9,2.,3.,4.);
    let b=box_(1.0e9+1.,1.0e9+1.,1.0e9+1.,2.,3.,4.);
    let result=run(&a,&b,BooleanOperation::Intersection);
    assert!((result.solid().unwrap().mass.volume-6.0).abs()<1e-8);
    let p=result.solid().unwrap().mass.centroid;
    assert!((p.x-(1.0e9+1.5)).abs()<1e-6);
    assert!((p.y-(1.0e9+2.0)).abs()<1e-6);
    assert!((p.z-(1.0e9+2.5)).abs()<1e-6);
}
