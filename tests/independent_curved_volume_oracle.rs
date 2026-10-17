//! G4c: source quadrature vs independent oriented curved-shell flux.
use rustsolid::{
    AnalyticFacePatchPartitionPlan,AnalyticFacePatchPartitionResult,
    BooleanOperation,BrepBody,CylindricalBrep,GeometryTolerance,
    IndependentVolumeOraclePlan,IndependentVolumeOracleResult,
    Point3,Solid,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn cube()->BrepBody{
    BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
}
fn cylinder()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}
fn prepare(a:&BrepBody,b:&BrepBody)->AnalyticFacePatchPartitionPlan{
    let s=AnalyticFacePatchPartitionPlan::build(
        a,b,BooleanOperation::Difference,tol()).unwrap();
    let AnalyticFacePatchPartitionResult::Ready(plan)=s else {
        panic!("expected bounded candidate: {s:?}");
    };plan
}
fn certify(a:&BrepBody,b:&BrepBody,
    candidate:&AnalyticFacePatchPartitionPlan)->IndependentVolumeOraclePlan{
    let result=IndependentVolumeOraclePlan::build(a,b,candidate,tol()).unwrap();
    let IndependentVolumeOracleResult::Ready(plan)=result else{
        panic!("independent source vs flux proof not ready: {result:?}");
    };plan
}
#[test]
fn four_corner_shells_have_independent_positive_source_and_flux_volume(){
    let (a,b)=(cube(),cylinder());
    let candidate=prepare(&a,&b);
    let proof=certify(&a,&b,&candidate);
    assert!(proof.volume_certified && !proof.can_apply);
    assert_eq!(proof.components.len(),4);
    assert_eq!(candidate.brep().shells.len(),4);
    assert!((proof.original_box_volume-27.).abs()<1e-11);
    if let BrepBody::Polyhedral(ref original)=a {
        assert_eq!(proof.source_revisions[0],Some(original.revision));
    }
    assert!(proof.source_total.integration_panels>4);
    assert!(proof.source_total.numerical_guard>0.);
    assert!(proof.source_total.lower>0.);
    assert!(proof.source_total.upper<proof.original_box_volume);
    assert!(proof.source_total.upper-proof.source_total.lower<1e-7);
    let mut total=0.;
    for c in &proof.components{
        assert!(c.source_volume.lower>0.);
        assert!(c.source_volume.upper>c.source_volume.lower);
        assert!(c.result_flux.result_volume>0.);
        assert!(c.result_flux.original_x_wall_flux>0.);
        assert!(c.result_flux.inward_cylinder_flux<0.);
        assert!(c.result_flux.cylinder_angle_span>0.);
        assert!(c.result_flux.result_volume>=c.source_volume.lower-c.permitted_absolute_error);
        assert!(c.result_flux.result_volume<=c.source_volume.upper+c.permitted_absolute_error);
        assert!((c.result_flux.result_volume-c.candidate_component_volume).abs()
            <c.permitted_absolute_error);
        total+=c.result_flux.result_volume;
    }
    assert!((total-proof.result_flux_total).abs()<1e-10);
    assert!((candidate.brep().mass.volume-proof.result_flux_total).abs()<1e-7);
    assert!((proof.source_total.lower+
        proof.removed_volume_bounds[1]-proof.original_box_volume).abs()<1e-9);
    assert!((proof.source_total.upper+
        proof.removed_volume_bounds[0]-proof.original_box_volume).abs()<1e-9);
    proof.validate_against(&a,&b,&candidate,tol()).unwrap();
}
#[test]
fn translated_off_center_non_square_box_certifies_all_components_separately(){
    let a=BrepBody::from(Solid::block(p(2.,6.,-8.),4.5,2.5,3.5).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        p(4.,1.,-6.2),2.55,10.,tol()).unwrap());
    let candidate=prepare(&a,&b);
    let proof=certify(&a,&b,&candidate);
    assert_eq!(proof.components.len(),4);
    let volumes=proof.components.iter()
        .map(|p|p.result_flux.result_volume).collect::<Vec<_>>();
    assert!(volumes.iter().all(|x|*x>0.));
    assert!(volumes.iter().any(|x|(volumes[0]-x).abs()>1e-5));
    assert!((proof.original_box_volume-4.5*2.5*3.5).abs()<1e-10);
    assert!(proof.removed_volume_bounds[0]>0.);
    proof.validate_against(&a,&b,&candidate,tol()).unwrap();
}
#[test]
fn oracle_and_candidate_tampering_rejected_by_independent_source_rebuild(){
    let (a,b)=(cube(),cylinder());
    let candidate=prepare(&a,&b);
    let proof=certify(&a,&b,&candidate);
    let mut damaged=proof.clone();
    damaged.components[0].result_flux.result_volume+=0.01;
    assert!(damaged.validate_against(&a,&b,&candidate,tol()).is_err());
    let mut damaged=proof.clone();
    damaged.components[1].source_volume.upper+=0.02;
    assert!(damaged.validate_against(&a,&b,&candidate,tol()).is_err());
    let mut damaged=proof.clone();
    damaged.candidate_mass_volume+=1.;
    assert!(damaged.validate_against(&a,&b,&candidate,tol()).is_err());
    let mut damaged=proof.clone();
    damaged.removed_volume_bounds.swap(0,1);
    assert!(damaged.validate_against(&a,&b,&candidate,tol()).is_err());
    let mut corrupted=candidate.clone();
    corrupted.corners[0].exact_volume+=1.;
    assert!(IndependentVolumeOraclePlan::build(&a,&b,&corrupted,tol()).is_err());
    let mut corrupted=candidate.clone();
    corrupted.material_witnesses[0].retained_point.x+=0.2;
    assert!(IndependentVolumeOraclePlan::build(&a,&b,&corrupted,tol()).is_err());
}
#[test]
fn old_source_revision_and_cylinder_geometry_invalidate_all_volume_evidence(){
    let (a,b)=(cube(),cylinder());
    let candidate=prepare(&a,&b);
    let proof=certify(&a,&b,&candidate);
    let mut edited=a.clone();
    if let BrepBody::Polyhedral(ref mut s)=edited{s.revision+=1;}
    assert!(proof.validate_against(&edited,&b,&candidate,tol()).is_err());
    let other=BrepBody::from(CylindricalBrep::upright(
        p(0.1,0.,0.),2.,5.,tol()).unwrap());
    assert!(proof.validate_against(&a,&other,&candidate,tol()).is_err());
    let changed=GeometryTolerance{
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8,
    };
    assert!(proof.validate_against(&a,&b,&candidate,changed).is_err());
}
#[test]
fn can_apply_remains_false_even_when_both_independent_volumes_pass(){
    let (a,b)=(cube(),cylinder());
    let candidate=prepare(&a,&b);
    let proof=certify(&a,&b,&candidate);
    assert!(proof.volume_certified);
    assert!(!proof.can_apply && !candidate.can_apply);
    assert!(proof.components.iter().all(|part|
        part.result_flux.original_x_wall_flux
            > part.result_flux.result_volume));
}
