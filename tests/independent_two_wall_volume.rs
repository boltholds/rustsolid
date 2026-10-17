use rustsolid::{
    AnalyticTwoWallPartitionPlan,AnalyticTwoWallPartitionResult,
    AnalyticTwoWallVolumeOracle,AnalyticTwoWallVolumeResult,
    BooleanOperation,BrepBody,CylindricalBrep,GeometryTolerance,
    Point3,Solid,
};
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tolerance()->GeometryTolerance{GeometryTolerance::default()}
fn a()->BrepBody{
    BrepBody::from(Solid::block(q(-1.2,1.,-3.),2.4,3.,6.).unwrap())
}
fn b()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(q(0.,0.,0.),1.5,5.,tolerance()).unwrap())
}
fn ready(a:&BrepBody,b:&BrepBody)
    ->(AnalyticTwoWallPartitionPlan,AnalyticTwoWallVolumeOracle){
    let AnalyticTwoWallPartitionResult::Ready(candidate)=
        AnalyticTwoWallPartitionPlan::build(
            a,b,BooleanOperation::Difference,tolerance()).unwrap()
    else{panic!("two-wall geometry source unsupported")};
    let r=AnalyticTwoWallVolumeOracle::build(a,b,&candidate,tolerance()).unwrap();
    let AnalyticTwoWallVolumeResult::Ready(proof)=r else{
        panic!("two-wall G6b independent volume Held: {r:?}");
    };
    (candidate,proof)
}
#[test]
fn two_material_shells_have_independent_source_brackets_and_signed_boundary_flux(){
    let (a,b)=(a(),b());
    let (candidate,proof)=ready(&a,&b);
    assert!(proof.volume_certified&&!proof.can_apply);
    assert_eq!(proof.components.len(),2);
    assert_eq!(candidate.brep().shells.len(),2);
    assert!((proof.original_box_volume-43.2).abs()<1e-10);
    assert!(proof.source_total.lower>0.);
    assert!(proof.source_total.upper-proof.source_total.lower<1e-6);
    assert!(proof.source_total.integration_panels>2);
    let mut sum=0.;
    for shell in &proof.components{
        assert!(shell.source.lower>0.);
        assert!(shell.source.upper>shell.source.lower);
        assert!(shell.boundary.result_volume>0.);
        assert!(shell.boundary.source_z_wall_flux>0.);
        assert!(shell.boundary.inward_cylinder_flux<0.);
        assert!(shell.boundary.angular_span>0.);
        assert!(shell.boundary.result_volume>=shell.source.lower-shell.allowed_error);
        assert!(shell.boundary.result_volume<=shell.source.upper+shell.allowed_error);
        assert!((shell.boundary.result_volume-shell.candidate_volume).abs()<shell.allowed_error);
        sum+=shell.boundary.result_volume;
    }
    assert!((sum-proof.boundary_total).abs()<1e-9);
    assert!((proof.stored_mass-proof.boundary_total).abs()<1e-7);
    assert!((proof.removed_inside_box[0]+proof.source_total.upper
        -proof.original_box_volume).abs()<1e-8);
    assert!((proof.removed_inside_box[1]+proof.source_total.lower
        -proof.original_box_volume).abs()<1e-8);
    proof.validate_against(&a,&b,&candidate,tolerance()).unwrap();
}
#[test]
fn asymmetric_and_translated_two_wall_volume_witnesses_are_distinct(){
    let a=BrepBody::from(Solid::block(q(1.,2.,-5.),3.,4.,9.).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        q(2.4,0.,-1.),1.8,8.,tolerance()).unwrap());
    let (candidate,proof)=ready(&a,&b);
    assert!(proof.components[0].boundary.result_volume
        !=proof.components[1].boundary.result_volume);
    assert_eq!(proof.components.len(),2);
    assert!(proof.boundary_total>0.);
    assert_eq!(proof.source_revisions[1],None);
    proof.validate_against(&a,&b,&candidate,tolerance()).unwrap();
}
#[test]
fn tampering_with_volume_stamps_and_G1_carriers_is_caught_before_any_commit(){
    let (a,b)=(a(),b());
    let (candidate,proof)=ready(&a,&b);
    let mut broken=proof.clone();
    broken.components[0].boundary.result_volume+=0.3;
    assert!(broken.validate_against(&a,&b,&candidate,tolerance()).is_err());
    let mut broken=proof.clone();
    broken.source_total.lower*=0.5;
    assert!(broken.validate_against(&a,&b,&candidate,tolerance()).is_err());
    let mut broken=proof.clone();
    broken.stored_mass+=1.;
    assert!(broken.validate_against(&a,&b,&candidate,tolerance()).is_err());
    let mut tampered=candidate.clone();
    tampered.components[0].global_arcs[0].pair=rustsolid::FacePairId(u32::MAX);
    assert!(AnalyticTwoWallVolumeOracle::build(&a,&b,&tampered,tolerance()).is_err());
}
#[test]
fn revisions_geometry_changes_and_different_tolerances_invalidate_certificate(){
    let (a,b)=(a(),b());
    let (candidate,proof)=ready(&a,&b);
    let mut original=a.clone();
    if let BrepBody::Polyhedral(ref mut s)=original{s.revision+=1;}
    assert!(proof.validate_against(&original,&b,&candidate,tolerance()).is_err());
    let changed=BrepBody::from(CylindricalBrep::upright(
        q(0.1,0.,0.),1.5,5.,tolerance()).unwrap());
    assert!(proof.validate_against(&a,&changed,&candidate,tolerance()).is_err());
    let custom=GeometryTolerance {
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8,
    };
    assert!(proof.validate_against(&a,&b,&candidate,custom).is_err());
}
#[test]
fn two_wall_result_approval_never_mutates_original_source_or_enables_general_boolean(){
    let (a,b)=(a(),b());
    let sa=format!("{a:?}");
    let sb=format!("{b:?}");
    let (candidate,proof)=ready(&a,&b);
    assert!(!proof.can_apply&&!candidate.can_apply);
    assert_eq!(format!("{a:?}"),sa);
    assert_eq!(format!("{b:?}"),sb);
}
