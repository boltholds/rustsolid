use rustsolid::{
    AnalyticOneWallPlan,AnalyticOneWallResult,
    AnalyticOneWallVolumeOracle,AnalyticOneWallVolumeResult,
    BooleanOperation,BrepBody,CylindricalBrep,
    GeometryTolerance,Point3,Solid,
};
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn a()->BrepBody{
    BrepBody::from(Solid::block(q(-1.,1.,-3.),4.,3.,6.).unwrap())
}
fn b()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(
        q(0.,0.,0.),1.5,5.,tol()).unwrap())
}
fn ready(left:&BrepBody,right:&BrepBody)
    ->(AnalyticOneWallPlan,AnalyticOneWallVolumeOracle){
    let AnalyticOneWallResult::Ready(candidate)=
        AnalyticOneWallPlan::build(
            left,right,BooleanOperation::Difference,tol()).unwrap()
    else{panic!("one-wall source cutoff not supported")};
    let proof=AnalyticOneWallVolumeOracle::build(
        left,right,&candidate,tol()).unwrap();
    let AnalyticOneWallVolumeResult::Ready(proof)=proof else {
        panic!("G7b independent proof Held: {proof:?}");
    };
    (candidate,proof)
}
#[test]
fn independent_source_polygons_and_real_boundary_flux_certify_one_shell(){
    let (left,right)=(a(),b());
    let (model,proof)=ready(&left,&right);
    assert!(proof.volume_certified&&!proof.can_apply);
    assert!(model.closed_shell_proven&&!model.can_apply);
    assert_eq!(proof.source.integration_panels>=16384,true);
    assert!(proof.source.lower>0.);
    assert!(proof.source.upper>proof.source.lower);
    assert!(proof.source.upper-proof.source.lower<1e-6);
    assert!(proof.source.numerical_guard>0.);
    assert_eq!(proof.result_flux.shell.0,0);
    assert!(proof.result_flux.volume>0.);
    assert!(proof.result_flux.outer_flux>proof.result_flux.volume);
    assert!(proof.result_flux.inner_flux<0.);
    assert!(proof.result_flux.true_angle_span>std::f64::consts::PI);
    assert_eq!(proof.result_flux.split_x_wall_faces.len(),2);
    assert!(proof.result_flux.volume>=proof.source.lower-1e-6);
    assert!(proof.result_flux.volume<=proof.source.upper+1e-6);
    assert!((proof.result_flux.volume-proof.candidate_volume).abs()<1e-6);
    assert!((proof.original_box_volume-72.).abs()<1e-10);
    assert!((proof.removed_inside_box[0]+proof.source.upper
        -proof.original_box_volume).abs()<1e-8);
    assert!((proof.removed_inside_box[1]+proof.source.lower
        -proof.original_box_volume).abs()<1e-8);
    proof.validate_against(&left,&right,&model,tol()).unwrap();
}
#[test]
fn translated_nonsymmetric_notch_has_valid_independent_mass_balance(){
    let left=BrepBody::from(Solid::block(q(4.,2.,-6.),5.,4.,8.).unwrap());
    let right=BrepBody::from(CylindricalBrep::upright(
        q(5.2,0.,-2.5),1.8,8.,tol()).unwrap());
    let (candidate,oracle)=ready(&left,&right);
    assert!(oracle.result_flux.volume>0.);
    assert!(oracle.original_box_volume>oracle.result_flux.volume);
    assert_eq!(candidate.brep().shells.len(),1);
    assert!(oracle.removed_inside_box[0]>0.);
    oracle.validate_against(&left,&right,&candidate,tol()).unwrap();
}
#[test]
fn falsified_polygon_bounds_flux_or_saved_volume_are_rejected(){
    let (left,right)=(a(),b());
    let (candidate,oracle)=ready(&left,&right);
    let mut damaged=oracle.clone();
    damaged.source.upper+=1.;
    assert!(damaged.validate_against(&left,&right,&candidate,tol()).is_err());
    let mut damaged=oracle.clone();
    damaged.result_flux.inner_flux-=1.;
    assert!(damaged.validate_against(&left,&right,&candidate,tol()).is_err());
    let mut damaged=oracle.clone();
    damaged.candidate_volume+=1.;
    assert!(damaged.validate_against(&left,&right,&candidate,tol()).is_err());
    let mut changed=candidate.clone();
    changed.component.new_cut_arcs[0].pair=rustsolid::FacePairId(u32::MAX);
    assert!(AnalyticOneWallVolumeOracle::build(
        &left,&right,&changed,tol()).is_err());
}
#[test]
fn source_revision_cylinder_changes_and_tolerance_changes_invalidate_evidence(){
    let (left,right)=(a(),b());
    let (candidate,oracle)=ready(&left,&right);
    let mut stale=left.clone();
    if let BrepBody::Polyhedral(ref mut s)=stale{s.revision+=1;}
    assert!(oracle.validate_against(&stale,&right,&candidate,tol()).is_err());
    let other=BrepBody::from(CylindricalBrep::upright(
        q(0.1,0.,0.),1.5,5.,tol()).unwrap());
    assert!(oracle.validate_against(&left,&other,&candidate,tol()).is_err());
    let changed=GeometryTolerance {
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8,
    };
    assert!(oracle.validate_against(&left,&right,&candidate,changed).is_err());
}
#[test]
fn independent_g7b_never_commits_mutates_sources_or_promotes_general_csg(){
    let (left,right)=(a(),b());
    let x=format!("{left:?}");let y=format!("{right:?}");
    let (candidate,proof)=ready(&left,&right);
    assert!(!candidate.can_apply);
    assert!(!proof.can_apply);
    assert_eq!(format!("{left:?}"),x);
    assert_eq!(format!("{right:?}"),y);
    assert_eq!(proof.source_revisions[1],None);
}
