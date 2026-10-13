//! Independent source/geometry contracts for a read-only mixed through/blind
//! circle overlap plan. No third-party CAD engine or B-rep geometry is used
//! as a truth oracle. The current mixed-body executor must stay fail-closed.

use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters, FootprintSource, HoleParameters,
    MixedFootprintContact, MixedFootprintOverlapPlanner,
    MixedThroughBlindBrep, MixedThroughBlindParameters, Point2, Point3, Solid,
};

fn origin(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(origin(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn through(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn blind(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn approx(value:f64,expected:f64,eps:f64){
    assert!((value-expected).abs()<=eps,
        "expected {expected:.14}, received {value:.14}, error {:.6e}",(value-expected).abs());
}

/// Deliberately independent numerical area oracle. For disks centered at 0
/// and D on one axis, each vertical strip's intersected width equals twice
/// the smaller chord height. Integrate with a cosine-variable transformation
/// to regularize the square-root singularities at support endpoints.
/// The exact analytic lens formula is NOT used here.
fn strip_integral_overlap(r1:f64,r2:f64,d:f64)->f64 {
    assert!(d>0.0);
    let start=(-r1).max(d-r2);
    let end=r1.min(d+r2);
    if start>=end{return 0.0;}
    let crossing_chord=(d*d+r1*r1-r2*r2)/(2.0*d);
    let n=4096usize;
    let integrate=|left:f64,right:f64|->f64{
        if right<=left{return 0.0;}
        let step=PI/(n as f64);
        let f=|theta:f64|->f64{
            let x=left+0.5*(right-left)*(1.0-theta.cos());
            let ch1=(r1*r1-x*x).max(0.0).sqrt();
            let ch2=(r2*r2-(x-d)*(x-d)).max(0.0).sqrt();
            2.0*ch1.min(ch2)*0.5*(right-left)*theta.sin()
        };
        let mut integral=f(0.0)+f(PI);
        for k in 1..n {
            integral+=if k%2==0{2.0}else{4.0}*f(step*k as f64);
        }
        integral*step/3.0
    };
    if crossing_chord>start && crossing_chord<end {
        integrate(start,crossing_chord)+integrate(crossing_chord,end)
    } else {
        integrate(start,end)
    }
}

#[test]
fn published_first_party_crossing_vector_matches_independent_strip_quadrature(){
    let a=through(0.,0.,0.5);
    let b=blind(0.5,0.,0.5,1.5);
    let box_=box_();
    let p=MixedFootprintOverlapPlanner::plan(&box_,a,b).unwrap();
    assert_eq!(p.contact,MixedFootprintContact::TwoTransverseCrossings);
    assert!(!p.can_apply());
    let exact_overlap=0.3070924246521893_f64;
    approx(p.overlap_area.unwrap(),exact_overlap,1e-12);
    approx(p.overlap_area.unwrap(),strip_integral_overlap(0.5,0.5,0.5),2e-9);
    approx(p.predicted_remaining_volume.unwrap(),140.14094873829234,1e-10);
    p.validate_against(&box_,a,b).unwrap();
}

#[test]
fn crossing_points_are_on_both_analytic_circles_and_arcs_are_exterior(){
    let a=through(0.,0.,0.5);
    let b=blind(0.5,0.,0.5,1.5);
    let p=MixedFootprintOverlapPlanner::plan(&box_(),a,b).unwrap();
    let points=p.crossing_points.unwrap();
    assert_ne!(points[0],points[1]);
    for v in points {
        approx((v.x-a.center_x).hypot(v.z-a.center_z),a.radius,1e-10);
        approx((v.x-b.center_x).hypot(v.z-b.center_z),b.radius,1e-10);
        approx(v.y,3.0,1e-12);
    }
    let arcs=p.exposed_arcs.unwrap();
    assert_eq!(arcs[0].source,FootprintSource::Through);
    assert_eq!(arcs[1].source,FootprintSource::Blind);
    assert_eq!(arcs[0].end_crossing,arcs[1].start_crossing);
    assert_eq!(arcs[1].end_crossing,arcs[0].start_crossing);
    for arc in arcs {
        assert!(arc.sweep_radians>0.0 && arc.sweep_radians<2.0*PI);
        let (ox,oz,r)=match arc.source {
            FootprintSource::Through=>(b.center_x,b.center_z,b.radius),
            FootprintSource::Blind=>(a.center_x,a.center_z,a.radius),
        };
        assert!((arc.midpoint.x-ox).hypot(arc.midpoint.z-oz)>r+p.length_budget);
    }
}

#[test]
fn disjoint_formula_is_identical_to_v029_source_geometry(){
    let a=through(-1.5,0.,0.45);
    let b=blind(1.5,0.,0.5,2.0);
    let box_=box_();
    let p=MixedFootprintOverlapPlanner::plan(&box_,a,b).unwrap();
    assert_eq!(p.contact,MixedFootprintContact::Disjoint);
    assert_eq!(p.overlap_area,Some(0.0));
    let expected=144.0-PI*a.radius*a.radius*4.0-PI*b.radius*b.radius*b.depth;
    approx(p.predicted_remaining_volume.unwrap(),expected,1e-11);
    let full=MixedThroughBlindBrep::from_box(&box_,MixedThroughBlindParameters{
        through:vec![a],blind:vec![b],
    }).unwrap();
    approx(p.predicted_remaining_volume.unwrap(),full.model().mass.volume,1e-10);
}

#[test]
fn nested_disks_obey_area_inclusion_and_distinct_cut_depths(){
    let box_=box_();
    let a=through(0.,0.,0.9);
    let b=blind(0.1,0.,0.3,1.8);
    let inside=MixedFootprintOverlapPlanner::plan(&box_,a,b).unwrap();
    assert_eq!(inside.contact,MixedFootprintContact::BlindInsideThrough);
    approx(inside.overlap_area.unwrap(),PI*b.radius*b.radius,1e-12);
    approx(inside.predicted_remaining_volume.unwrap(),144.0-PI*0.81*4.0,1e-10);

    let a=through(0.0,0.0,0.3);
    let b=blind(0.1,0.,0.9,2.0);
    let outer=MixedFootprintOverlapPlanner::plan(&box_,a,b).unwrap();
    assert_eq!(outer.contact,MixedFootprintContact::ThroughInsideBlind);
    approx(outer.overlap_area.unwrap(),PI*0.09,1e-12);
    approx(outer.predicted_remaining_volume.unwrap(),
        144.0-PI*0.09*4.0-PI*(0.81-0.09)*2.0,1e-10);

    let same=MixedFootprintOverlapPlanner::plan(&box_,
        through(0.,0.,0.5),blind(0.,0.,0.5,1.5)).unwrap();
    assert_eq!(same.contact,MixedFootprintContact::Coincident);
    approx(same.overlap_area.unwrap(),PI*0.25,1e-12);
    approx(same.predicted_remaining_volume.unwrap(),144.0-PI*0.25*4.0,1e-10);
    assert!(same.exposed_arcs.is_none());
    assert!(!same.can_apply());
}

#[test]
fn external_and_internal_near_tangency_hold_instead_of_inventing_crossings(){
    let box_=box_();
    for d in [1.0,1.0-2e-10,1.0+2e-10] {
        let p=MixedFootprintOverlapPlanner::plan(&box_,
            through(0.,0.,0.5),blind(d,0.,0.5,1.5)).unwrap();
        assert_eq!(p.contact,MixedFootprintContact::HeldNearTangency);
        assert!(p.overlap_area.is_none() && p.predicted_remaining_volume.is_none());
        assert!(p.crossing_points.is_none() && p.exposed_arcs.is_none());
        assert!(!p.can_apply());
    }
    for d in [0.5,0.5+1e-10,0.5-1e-10] {
        let p=MixedFootprintOverlapPlanner::plan(&box_,
            through(0.,0.,0.9),blind(d,0.,0.4,1.5)).unwrap();
        assert_eq!(p.contact,MixedFootprintContact::HeldNearTangency);
    }
    let nearly_coincident=MixedFootprintOverlapPlanner::plan(&box_,
        through(0.,0.,0.5),blind(1e-11,0.,0.5000000001,1.5)).unwrap();
    assert_eq!(nearly_coincident.contact,MixedFootprintContact::HeldNearTangency);
}

#[test]
fn independent_numeric_oracle_matches_unequal_radii_in_both_circle_orders(){
    let box_=box_();
    for (r1,r2,d) in [
        (0.5,0.5,0.5), (0.7,0.5,0.6), (0.95,0.4,0.8),
        (0.65,0.9,0.6), (0.9,0.75,0.3), (0.45,0.95,1.1),
    ]{
        let through=through(0.,0.,r1);
        let blind=blind(d,0.,r2,1.3);
        let p=MixedFootprintOverlapPlanner::plan(&box_,through,blind).unwrap();
        assert_eq!(p.contact,MixedFootprintContact::TwoTransverseCrossings,
            "failed for {r1}, {r2}, {d}");
        let integral=strip_integral_overlap(r1,r2,d);
        approx(p.overlap_area.unwrap(),integral,1e-7);
        let reverse=MixedFootprintOverlapPlanner::plan(&box_,
            through(d,0.,r2),blind(0.,0.,r1,1.3)).unwrap();
        approx(reverse.overlap_area.unwrap(),p.overlap_area.unwrap(),1e-10);
        approx(reverse.predicted_remaining_volume.unwrap(),
            144.0-PI*r2*r2*4.0-(PI*r1*r1-integral)*1.3,1e-7);
    }
}

#[test]
fn changing_blind_depth_is_monotonic_for_crossing_disks(){
    let box_=box_();
    let through=through(0.,0.,0.5);
    let shallow=MixedFootprintOverlapPlanner::plan(&box_,through,blind(0.5,0.,0.5,0.5)).unwrap();
    let deep=MixedFootprintOverlapPlanner::plan(&box_,through,blind(0.5,0.,0.5,2.8)).unwrap();
    assert!(deep.predicted_remaining_volume.unwrap() < shallow.predicted_remaining_volume.unwrap());
    let outside=PI*0.25-shallow.overlap_area.unwrap();
    approx(shallow.predicted_remaining_volume.unwrap()-deep.predicted_remaining_volume.unwrap(),
        outside*(2.8-0.5),1e-10);
}

#[test]
fn rotating_and_translating_source_preserves_analytic_volume(){
    let sample=MixedFootprintOverlapPlanner::plan(&box_(),
        through(0.,0.,0.5),blind(0.5,0.,0.5,1.5)).unwrap();
    let rotated=MixedFootprintOverlapPlanner::plan(&box_(),
        through(0.,0.,0.5),blind(0.,0.5,0.5,1.5)).unwrap();
    approx(sample.overlap_area.unwrap(),rotated.overlap_area.unwrap(),1e-12);
    approx(sample.predicted_remaining_volume.unwrap(),
        rotated.predicted_remaining_volume.unwrap(),1e-10);

    let dx=1e6;
    let dz=-1e6;
    let translated=Solid::block(origin(dx-3.,-1.,dz-3.),6.,4.,6.).unwrap();
    let shifted=MixedFootprintOverlapPlanner::plan(&translated,
        through(dx,dz,0.5),blind(dx+0.5,dz,0.5,1.5)).unwrap();
    assert_eq!(shifted.contact,MixedFootprintContact::TwoTransverseCrossings);
    approx(sample.overlap_area.unwrap(),shifted.overlap_area.unwrap(),1e-10);
    approx(sample.predicted_remaining_volume.unwrap(),
        shifted.predicted_remaining_volume.unwrap(),1e-9);
}

#[test]
fn source_or_result_tampering_is_detected_without_editing_source(){
    let source=box_();
    let through=through(0.,0.,0.5);
    let blind=blind(0.5,0.,0.5,1.5);
    let p=MixedFootprintOverlapPlanner::plan(&source,through,blind).unwrap();
    let before=(source.revision,source.mass,source.faces.len(),source.vertices.len());
    p.validate_against(&source,through,blind).unwrap();

    let mut tampered=p.clone();
    tampered.overlap_area=Some(tampered.overlap_area.unwrap()+0.01);
    assert!(tampered.validate_against(&source,through,blind).is_err());
    let mut tampered=p.clone();
    tampered.crossing_points.as_mut().unwrap()[0].x+=0.05;
    assert!(tampered.validate_against(&source,through,blind).is_err());
    let mut tampered=p.clone();
    tampered.length_budget*=2.0;
    assert!(tampered.validate_against(&source,through,blind).is_err());
    assert!(p.validate_against(&source,through,blind(0.6,0.,0.5,1.5)).is_err());

    assert_eq!(before,(source.revision,source.mass,source.faces.len(),source.vertices.len()));
    assert!(!p.can_apply());
}

#[test]
fn existing_mixed_body_still_rejects_through_blind_crossing(){
    let source=box_();
    let through=through(0.,0.,0.5);
    let blind=blind(0.5,0.,0.5,1.5);
    let p=MixedFootprintOverlapPlanner::plan(&source,through,blind).unwrap();
    assert_eq!(p.contact,MixedFootprintContact::TwoTransverseCrossings);
    assert!(!p.can_apply());
    let before=source.mass;
    let attempt=MixedThroughBlindBrep::from_box(&source,MixedThroughBlindParameters{
        through:vec![through],blind:vec![blind]
    });
    assert!(attempt.is_err(),"geometry-only overlap plan must not unlock B-rep commit");
    assert_eq!(source.mass,before);
    assert_eq!(source.revision,0);
}

#[test]
fn nonbox_and_invalid_cutter_inputs_fail_closed(){
    let invalid_solid=Solid::extrude_xz(&[
        Point2{x:-3.,z:-3.},Point2{x:3.,z:-3.},
        Point2{x:2.,z:3.},Point2{x:-3.,z:3.}],4.0).unwrap();
    assert!(MixedFootprintOverlapPlanner::plan(&invalid_solid,
        through(0.,0.,0.5),blind(0.5,0.,0.5,1.0)).is_err());
    let source=box_();
    for a in [
        through(2.9,0.,0.5),through(0.,0.,0.),
        through(0.,0.,-1.),through(f64::NAN,0.,0.5),
    ] {
        assert!(MixedFootprintOverlapPlanner::plan(&source,a,blind(0.5,0.,0.5,1.5)).is_err());
    }
    for b in [
        blind(0.5,0.,0.5,0.),blind(0.5,0.,0.5,4.),
        blind(0.5,0.,0.5,4.1),blind(0.5,0.,0.5,-0.01),
        blind(2.8,0.,0.5,1.5),blind(0.5,0.,f64::INFINITY,1.5),
    ] {
        assert!(MixedFootprintOverlapPlanner::plan(&source,through(0.,0.,0.5),b).is_err());
    }
}

#[test]
fn under_resolved_contact_never_reports_false_crossing(){
    let source=Solid::block_with_tolerance(origin(-3.,-1.,-3.),6.,4.,6.,
        rustsolid::GeometryTolerance{
            absolute_length:1e-3,relative_length:0.,angular:1e-8
        }).unwrap();
    let p=MixedFootprintOverlapPlanner::plan(&source,
        through(0.,0.,0.5),blind(0.9995,0.,0.5,1.5)).unwrap();
    assert_eq!(p.contact,MixedFootprintContact::HeldNearTangency);
    assert!(p.overlap_area.is_none());
}
