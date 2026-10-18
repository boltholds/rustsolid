use rustsolid::{
    AnalyticOneWallPlan,AnalyticOneWallResult,BooleanOperation,
    BrepBody,CylindricalBrep,Frame3,GeometryTolerance,
    IndependentTiltedVolumeOracle,IndependentTiltedVolumeResult,
    Point3,Solid,TiltedOneWallSourcePlan,TiltedOneWallSourceResult,
    TiltedVolumeHold,
};
use std::f64::consts::{PI,TAU};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn box_body()->BrepBody{
    BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap())
}
fn tilted_x(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(theta.cos(),-theta.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn tilted_z(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(1.,0.,0.),
        p(0.,-theta.sin(),-theta.cos()),tol()).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn ready(a:&BrepBody,b:&BrepBody)
    ->(TiltedOneWallSourcePlan,IndependentTiltedVolumeOracle)
{
    let TiltedOneWallSourceResult::Ready(source)=
        TiltedOneWallSourcePlan::build(a,b,BooleanOperation::Difference,tol())
        .unwrap() else{panic!("not G8a Ready")};
    let result=IndependentTiltedVolumeOracle::build(a,b,&source,tol()).unwrap();
    let IndependentTiltedVolumeResult::Ready(oracle)=result else{
        panic!("G8b Held: {result:?}");
    };
    (source,oracle)
}
fn assert_valid(
    a:&BrepBody,b:&BrepBody,
    source:&TiltedOneWallSourcePlan,oracle:&IndependentTiltedVolumeOracle,
){
    assert!(!oracle.can_apply&&!oracle.result_shell_certified);
    assert!(oracle.source_bracket_established);
    assert!(oracle.angular_panels.is_power_of_two());
    assert!(oracle.angular_panels>=1024);
    assert_eq!(oracle.height_panels%2,0);
    assert_eq!(oracle.retained.integration_panels,oracle.height_panels);
    assert_eq!(oracle.retained.numerical_guard,oracle.numerical_guard);
    let [lo,hi]=oracle.removed_inside_box;
    assert!(lo>0.&&hi>lo);
    assert!((lo-oracle.analytic_removed_cross_check).abs()
        <= hi-lo+oracle.analytic_cross_check_error+oracle.numerical_guard);
    assert!(oracle.analytic_removed_cross_check>=lo-oracle.numerical_guard
        -oracle.analytic_cross_check_error);
    assert!(oracle.analytic_removed_cross_check<=hi+oracle.numerical_guard
        +oracle.analytic_cross_check_error);
    assert!(oracle.retained.lower>0.);
    assert!(oracle.retained.upper<oracle.original_box_volume);
    assert!(oracle.retained.upper-oracle.retained.lower<1e-5);
    assert!((oracle.retained.lower+hi-oracle.original_box_volume).abs()<1e-10);
    assert!((oracle.retained.upper+lo-oracle.original_box_volume).abs()<1e-10);
    for section in &oracle.section_witnesses{
        assert!(section.lower_area<section.upper_area);
        assert!(section.normalized_wall_offset>0.);
        assert!(section.normalized_wall_offset<1.);
        assert!(section.major_arc_span>PI);
        assert!(section.major_arc_span<TAU);
        assert_eq!(section.angular_panels,oracle.angular_panels);
    }
    oracle.validate_against(a,b,source,tol()).unwrap();
}

#[test]
fn g8b_swept_ellipse_has_narrow_source_only_enclosure_at_both_x_tilt_directions(){
    for theta in [0.04,-0.03,0.11]{
        let (a,b)=(box_body(),tilted_x(theta));
        let before_a=format!("{a:?}");let before_b=format!("{b:?}");
        let (source,oracle)=ready(&a,&b);
        assert_valid(&a,&b,&source,&oracle);
        assert!(oracle.retained.upper-oracle.retained.lower<2e-7);
        assert_eq!(format!("{a:?}"),before_a);
        assert_eq!(format!("{b:?}"),before_b);
    }
}
#[test]
fn g8b_tilt_in_yz_plane_integrates_constant_projected_ellipse_area(){
    let (a,b)=(box_body(),tilted_z(0.04));
    let (source,oracle)=ready(&a,&b);
    assert_valid(&a,&b,&source,&oracle);
    let BrepBody::Cylinder(cyl)=b else{panic!("cylinder expected")};
    let jac=cyl.radius*cyl.radius/cyl.frame.normal.y;
    let d=(cyl.frame.origin.x-(-1.))/cyl.radius;
    let removed=3.*jac*(PI/2.+d.asin()+d*(1.-d*d).sqrt());
    assert!(oracle.removed_inside_box[0]<removed);
    assert!(oracle.removed_inside_box[1]>removed);
    assert!((oracle.analytic_removed_cross_check-removed).abs()<2e-9);
    assert_eq!(oracle.height_panels,8);
    assert!(oracle.simpson_remainder_bound<1e-20);
    let [a,b,c]=oracle.section_witnesses;
    assert!((a.lower_area-b.lower_area).abs()<1e-11);
    assert!((b.lower_area-c.lower_area).abs()<1e-11);
}
#[test]
fn g8b_polygon_formula_matches_actual_ellipse_vertex_shoelace(){
    let (a,b)=(box_body(),tilted_x(0.04));
    let (_,oracle)=ready(&a,&b);
    let BrepBody::Cylinder(cyl)=b else{panic!("expected cylinder")};
    let s=oracle.section_witnesses[1]; // interior Y
    let axis=cyl.frame.normal;let frame=cyl.frame;let r=cyl.radius;
    let c=p(
        frame.origin.x+axis.x*(s.height_y-frame.origin.y)/axis.y,
        s.height_y,
        frame.origin.z+axis.z*(s.height_y-frame.origin.y)/axis.y);
    let ux=p(
        r*(frame.x.x-axis.x*frame.x.y/axis.y),0.,
        r*(frame.x.z-axis.z*frame.x.y/axis.y));
    let uy=p(
        r*(frame.y.x-axis.x*frame.y.y/axis.y),0.,
        r*(frame.y.z-axis.z*frame.y.y/axis.y));
    let phase=uy.x.atan2(ux.x);
    let lo=phase-s.major_arc_span/2.;
    let extent=s.major_arc_span;
    let n=128usize;
    let step=extent/(n as f64);
    let q=|angle:f64,radial:f64|->(f64,f64){
        let (si,co)=angle.sin_cos();
        (c.x+radial*(ux.x*co+uy.x*si),
         c.z+radial*(ux.z*co+uy.z*si))
    };
    let shoelace=|points:&[(f64,f64)]|->f64{
        0.5*points.iter().enumerate().map(|(i,a)|{
            let b=points[(i+1)%points.len()];
            a.0*b.1-b.0*a.1
        }).sum::<f64>().abs()
    };
    let inner=(0..=n).map(|i|q(lo+i as f64*step,1.))
        .collect::<Vec<_>>();
    let tangent=(0..n).map(|i|
        q(lo+(i as f64+0.5)*step,1./(step/2.).cos()))
        .collect::<Vec<_>>();
    let mut outside=vec![inner[0]];
    outside.extend(tangent);
    outside.push(inner[n]);
    let area_inner=shoelace(&inner);
    let area_outer=shoelace(&outside);
    let jac=s.absolute_jacobian;
    let compact_inner=0.5*jac*((n as f64)*step.sin()-extent.sin());
    let compact_gap=jac*(n as f64)*(step/2.).sin().powi(3)/(step/2.).cos();
    assert!((area_inner-compact_inner).abs()<2e-11);
    assert!((area_outer-(compact_inner+compact_gap)).abs()<2e-11);
    assert!(area_outer>area_inner);
    assert!((inner[0].0-(-1.)).abs()<1e-10);
    assert!((inner[n].0-(-1.)).abs()<1e-10);
}
#[test]
fn g8b_shifted_offcenter_box_does_not_change_interval_semantics(){
    let theta=0.03_f64;
    let a=BrepBody::from(Solid::block(p(4.,2.,-6.),5.,4.,8.).unwrap());
    let frame=Frame3::from_axes(
        p(4.8,-3.,-2.5),p(theta.cos(),-theta.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let b=BrepBody::from(CylindricalBrep::new(frame,1.8,14.,tol()).unwrap());
    let (source,oracle)=ready(&a,&b);
    assert_valid(&a,&b,&source,&oracle);
    assert_eq!(oracle.original_box_volume,160.);
    assert!(oracle.analytic_removed_cross_check>0.);
    assert!(oracle.analytic_removed_cross_check<160.);
}
#[test]
fn g8b_rejects_corrupted_brackets_revisions_and_source_provenance(){
    let (a,b)=(box_body(),tilted_x(0.04));
    let (source,oracle)=ready(&a,&b);
    let mut tampered=oracle.clone();
    tampered.retained.lower+=0.1;
    assert!(tampered.validate_against(&a,&b,&source,tol()).is_err());
    let mut tampered=oracle.clone();
    tampered.section_witnesses[0].major_arc_span+=0.02;
    assert!(tampered.validate_against(&a,&b,&source,tol()).is_err());
    let mut tampered=oracle.clone();
    tampered.result_shell_certified=true;
    assert!(tampered.validate_against(&a,&b,&source,tol()).is_err());
    let mut changed=a.clone();
    if let BrepBody::Polyhedral(s)=&mut changed{s.revision+=1;}
    assert!(oracle.validate_against(&changed,&b,&source,tol()).is_err());
    let mut bad_source=source.clone();
    bad_source.top.center.x+=0.01;
    assert!(IndependentTiltedVolumeOracle::build(&a,&b,&bad_source,tol()).is_err());
    let mut changed_tol=tol();
    changed_tol.relative_length+=1e-11;
    assert!(oracle.validate_against(&a,&b,&source,changed_tol).is_err());
}
#[test]
fn g8b_keeps_upright_g7_contract_outside_the_tilted_oracle(){
    let a=box_body();
    let b=BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),1.5,5.,tol()).unwrap());
    let r=AnalyticOneWallPlan::build(&a,&b,BooleanOperation::Difference,tol()).unwrap();
    assert!(matches!(r,AnalyticOneWallResult::Ready(_)));
    assert!(matches!(TiltedOneWallSourcePlan::build(
        &a,&b,BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(_)));
}
#[test]
fn g8b_very_strict_f64_tolerance_never_promotes_unresolved_volume_to_approved(){
    // Only a genuine G8a-ready input can be evaluated. If strict
    // upstream geometric tolerances are unresolved, the source stage
    // must Held. Otherwise, G8b must Held on its numeric budget or
    // return a true interval while NEVER approving a Boolean commit.
    let a=box_body();
    let b=tilted_x(0.04);
    let strict=GeometryTolerance{
        absolute_length:1e-13,relative_length:1e-16,angular:1e-8,
    };
    match TiltedOneWallSourcePlan::build(&a,&b,BooleanOperation::Difference,strict){
        Err(_)|Ok(TiltedOneWallSourceResult::Held(_))=>{},
        Ok(TiltedOneWallSourceResult::Ready(source))=>{
            match IndependentTiltedVolumeOracle::build(&a,&b,&source,strict).unwrap(){
                IndependentTiltedVolumeResult::Held(
                    TiltedVolumeHold::NumericalBudgetExceeded|
                    TiltedVolumeHold::UnresolvedPolygonEnclosure|
                    TiltedVolumeHold::UnresolvedHeightQuadrature)=>{},
                IndependentTiltedVolumeResult::Ready(oracle)=>{
                    assert!(oracle.source_bracket_established);
                    assert!(!oracle.can_apply&&!oracle.result_shell_certified);
                },
                other=>panic!("unexpected precision outcome: {other:?}"),
            }
        }
    }
}
