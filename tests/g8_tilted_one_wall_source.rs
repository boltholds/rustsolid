use std::collections::{BTreeMap,BTreeSet};
use rustsolid::{
    AnalyticOneWallPlan,AnalyticOneWallResult,BooleanOperation,
    BrepBody,Curve2,Curve3,CylindricalBrep,Frame3,GeometryTolerance,
    GlobalEventId,Point3,Solid,TiltedOneWallHold,
    TiltedOneWallSourcePlan,TiltedOneWallSourceResult,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn body()->BrepBody{
    BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap())
}
fn cutter_x_tilt(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),
        p(theta.cos(),-theta.sin(),0.),p(0.,0.,-1.),
        tol(),
    ).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn cutter_z_tilt(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(1.,0.,0.),
        p(0.,-theta.sin(),-theta.cos()),tol(),
    ).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn ready(a:&BrepBody,b:&BrepBody)->TiltedOneWallSourcePlan{
    let result=TiltedOneWallSourcePlan::build(
        a,b,BooleanOperation::Difference,tol()).unwrap();
    let TiltedOneWallSourceResult::Ready(plan)=result else{
        panic!("G8a source geometry Held: {result:?}");
    };
    plan
}
fn assert_source_cycle(plan:&TiltedOneWallSourcePlan,body:&BrepBody,cutter:&BrepBody){
    assert!(!plan.can_apply&&!plan.closed_shell_proven);
    assert_eq!(plan.source_cuts.arcs.len(),4);
    assert_eq!(plan.source_cuts.physical_event_count,4);
    assert_eq!(plan.cut_arcs.len(),4);
    assert_eq!(plan.events.len(),4);
    let cuts=plan.cut_arcs.iter().copied().collect::<BTreeSet<_>>();
    assert_eq!(cuts.len(),4);
    let unique=plan.events.iter().copied().collect::<BTreeSet<_>>();
    assert_eq!(unique.len(),4);
    let mut counts=BTreeMap::<GlobalEventId,usize>::new();
    for a in &plan.source_cuts.arcs{
        let [u,v]=a.physical_endpoints.expect("physical G1 endpoints");
        *counts.entry(u).or_default()+=1;
        *counts.entry(v).or_default()+=1;
        let first=body.shared().unwrap();
        let second=cutter.shared().unwrap();
        let pa=first.face_surface(a.first_face).unwrap();
        let pb=second.face_surface(a.second_face).unwrap();
        let rustsolid::ParameterRange::Bounded{start,end}=a.parameter_domain
            else{panic!("unbounded G8a curve")};
        for f in [0.,0.1,0.25,0.5,0.75,0.9,1.]{
            let t=start+(end-start)*f;
            let point=a.carrier.evaluate(t).unwrap();
            let uv0=a.face_pcurves[0].evaluate(t).unwrap();
            let uv1=a.face_pcurves[1].evaluate(t).unwrap();
            assert!(distance(point,pa.evaluate(uv0.u,uv0.v).unwrap().point)<2e-7);
            assert!(distance(point,pb.evaluate(uv1.u,uv1.v).unwrap().point)<2e-7);
        }
    }
    assert_eq!(counts.len(),4);
    assert!(counts.values().all(|degree|*degree==2));
    plan.validate_against(body,cutter,tol()).unwrap();
}
#[test]
fn g8a_tilted_x_axis_certifies_four_source_arcs_two_major_ellipses_and_one_physical_cycle(){
    let (a,b)=(body(),cutter_x_tilt(0.04));
    let old_a=format!("{a:?}");let old_b=format!("{b:?}");
    let plan=ready(&a,&b);
    assert_source_cycle(&plan,&a,&b);
    for section in [plan.bottom,plan.top]{
        assert!(section.retained_major_span>std::f64::consts::PI);
        assert!(section.retained_major_span<std::f64::consts::TAU);
        assert!(section.wall_points[0].z<section.wall_points[1].z ||
                section.wall_points[1].z<section.wall_points[0].z);
    }
    let b_sections=plan.source_cuts.arcs.iter().filter(|arc|
        arc.first_face==plan.box_bottom||arc.first_face==plan.box_top
    ).collect::<Vec<_>>();
    assert_eq!(b_sections.len(),2);
    assert!(b_sections.iter().all(|arc|
        matches!(arc.carrier,Curve3::Harmonic(_)) &&
        matches!(arc.face_pcurves[1],Curve2::CylinderSection(_))));
    assert_eq!(format!("{a:?}"),old_a);
    assert_eq!(format!("{b:?}"),old_b);
}
#[test]
fn g8a_ellipse_implicit_distance_and_face_events_independent_of_g1(){
    let (a,b)=(body(),cutter_x_tilt(-0.03));
    let plan=ready(&a,&b);
    let BrepBody::Cylinder(cyl)=b else{unreachable!()};
    let origin=cyl.frame.origin;
    let axis=cyl.frame.normal;
    for section in [plan.bottom,plan.top]{
        assert!(section.cosine.y.abs()<1e-12);
        assert!(section.sine.y.abs()<1e-12);
        assert!((section.center.y-section.height_y).abs()<1e-12);
        for i in 0..=64 {
            let u=std::f64::consts::TAU*(i as f64)/64.;
            let v=section.evaluate(u).unwrap();
            let delta=v.sub(origin);
            let radial=delta.sub(axis.scale(delta.dot(axis)));
            assert!((radial.length_squared()-cyl.radius*cyl.radius).abs()<1e-10);
            assert!((v.y-section.height_y).abs()<1e-10);
        }
        for &v in &section.wall_points{
            assert!((v.x+1.).abs()<1e-10);
            assert!((v.y-section.height_y).abs()<1e-10);
            assert!(v.z>-3.&&v.z<3.);
        }
    }
}
#[test]
fn g8a_source_revision_and_provenance_tampering_are_detected(){
    let (a,b)=(body(),cutter_x_tilt(0.04));
    let plan=ready(&a,&b);
    let mut changed=a.clone();
    if let BrepBody::Polyhedral(s)=&mut changed{s.revision+=1;}
    assert!(plan.validate_against(&changed,&b,tol()).is_err());
    let mut tampered=plan.clone();
    tampered.top.retained_major_span-=0.05;
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut tampered=plan.clone();
    tampered.source_cuts.arcs[0].physical_endpoints=None;
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
}
#[test]
fn g8a_conservative_rejections_preserve_g7_and_only_preview_eligible_tilts(){
    let a=body();
    let b=cutter_x_tilt(0.04);
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&b,
        BooleanOperation::Union,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(TiltedOneWallHold::UnsupportedOperation)));
    // Upright G7 remains fully executable through its *own* bounded
    // document, never through this read-only G8a source preview.
    let upright=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),1.5,5.,tol()).unwrap());
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&upright,
        BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(TiltedOneWallHold::NotPositiveYTransverseTilt)));
    assert!(matches!(AnalyticOneWallPlan::build(&a,&upright,
        BooleanOperation::Difference,tol()).unwrap(),
        AnalyticOneWallResult::Ready(_)));
    let short_frame=Frame3::from_axes(
        p(0.,0.9,0.),
        p(0.04_f64.cos(),-0.04_f64.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let short=BrepBody::from(CylindricalBrep::new(
        short_frame,1.5,3.3,tol()).unwrap());
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&short,
        BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(TiltedOneWallHold::CutterCapsReachBox)));
    // Cylinder reaches Z walls or the MAX-X face, so not a one-wall cut.
    let broad_frame=Frame3::from_axes(
        p(0.,-2.,0.),p(0.04_f64.cos(),-0.04_f64.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let broad=BrepBody::from(CylindricalBrep::new(
        broad_frame,4.0,10.,tol()).unwrap());
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&broad,
        BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(TiltedOneWallHold::NotIsolatedMinXWall)));
}
#[test]
fn g8a_translated_offcenter_and_z_tilt_are_source_owned(){
    let a=BrepBody::from(Solid::block(p(4.,2.,-6.),5.,4.,8.).unwrap());
    let theta=0.03_f64;
    let frame=Frame3::from_axes(
        p(4.8,-3.,-2.5),p(theta.cos(),-theta.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let b=BrepBody::from(CylindricalBrep::new(frame,1.8,14.,tol()).unwrap());
    assert_source_cycle(&ready(&a,&b),&a,&b);
    // Perpendicular tilted direction: MIN-X wall cut is represented by
    // two original analytic line generators rather than conic curves.
    let a=body();
    let b=cutter_z_tilt(0.04);
    let plan=ready(&a,&b);
    assert_source_cycle(&plan,&a,&b);
    let wall=plan.source_cuts.arcs.iter().filter(|arc|
        arc.first_face==plan.box_min_x).collect::<Vec<_>>();
    assert_eq!(wall.len(),2);
    assert!(wall.iter().all(|arc|
        matches!(arc.carrier,Curve3::Line(_)) &&
        matches!(arc.face_pcurves[1],Curve2::Line(_))));
}
