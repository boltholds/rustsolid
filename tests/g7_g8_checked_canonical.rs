//! The production G7/G8 plan exposes checked canonical topology during migration.
//! Checked topology is NOT permission for Boolean commit.
use rustsolid::{AnalyticOneWallPlan,AnalyticOneWallResult,BooleanOperation,BrepBody,
    CylindricalBrep,Frame3,GeometryTolerance,IndependentTiltedVolumeOracle,
    IndependentTiltedVolumeResult,Point3,Solid,TiltedOneWallShellPlan,
    TiltedOneWallShellResult,TiltedOneWallSourcePlan,TiltedOneWallSourceResult};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn box_body()->BrepBody {BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap())}
#[test]
fn g7_builder_exports_checked_canonical_model(){
    let a=box_body();
    let b=BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),1.5,5.,tol()).unwrap());
    let AnalyticOneWallResult::Ready(g7)=AnalyticOneWallPlan::build(&a,&b,BooleanOperation::Difference,tol()).unwrap()
        else{panic!("G7 failed")};
    let canonical=g7.checked_canonical().unwrap();
    assert_eq!((canonical.vertices().len(),canonical.edges().len(),canonical.faces().len()),(12,18,8));
    assert_eq!(canonical.carriers().curves3,g7.brep().geometry.curves3);
}
#[test]
fn g8_builder_exports_checked_canonical_model(){
    let a=box_body();let angle=0.04_f64;
    let frame=Frame3::from_axes(p(0.,-2.,0.),p(angle.cos(),-angle.sin(),0.),p(0.,0.,-1.),tol()).unwrap();
    let b=BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap());
    let TiltedOneWallSourceResult::Ready(src)=TiltedOneWallSourcePlan::build(&a,&b,BooleanOperation::Difference,tol()).unwrap()
        else{panic!("G8a failed")};
    let IndependentTiltedVolumeResult::Ready(volume)=IndependentTiltedVolumeOracle::build(&a,&b,&src,tol()).unwrap()
        else{panic!("G8b failed")};
    let TiltedOneWallShellResult::Ready(g8)=TiltedOneWallShellPlan::build(&a,&b,&src,&volume,tol()).unwrap()
        else{panic!("G8c failed")};
    let canonical=g8.checked_canonical().unwrap();
    assert_eq!((canonical.vertices().len(),canonical.edges().len(),canonical.faces().len()),(12,18,8));
    assert_eq!(canonical.carriers().curves3,g8.brep().geometry.curves3);
}
