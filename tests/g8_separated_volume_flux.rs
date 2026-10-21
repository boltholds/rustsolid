//! Qualify source-only G8b and actual v2-coedge G8c against the existing G8.
use rustsolid::{BrepBody,BooleanOperation,CylindricalBrep,Frame3,GeometryTolerance,
    IndependentTiltedVolumeOracle,IndependentTiltedVolumeResult,
    IndependentTiltedBoundaryFluxOracle,IndependentTiltedBoundaryFluxResult,
    Point3,Solid,TiltedOneWallSourcePlan,TiltedOneWallSourceResult,
    TiltedOneWallShellPlan,TiltedOneWallShellResult};
use rustsolid_source::{BooleanKind,SourceGeometry,SourceId,SourceRevision,
    SourceSnapshot,SourceSpec,UnitPolicy};
use rustsolid_geometry::{ModelUnits,ParameterRange};
use rustsolid_validation::{evaluate_tilted_one_wall,evaluate_tilted_boundary_flux};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn run(angle:f64,tilt_x:bool){
    let frame=if tilt_x {
        Frame3::from_axes(p(0.,-2.,0.),p(angle.cos(),-angle.sin(),0.),
            p(0.,0.,-1.),tol()).unwrap()
    }else{
        Frame3::from_axes(p(0.,-2.,0.),p(1.,0.,0.),
            p(0.,-angle.sin(),-angle.cos()),tol()).unwrap()
    };
    let a=BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
    let b=BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap());
    let TiltedOneWallSourceResult::Ready(s)=TiltedOneWallSourcePlan::build(
        &a,&b,BooleanOperation::Difference,tol()).unwrap() else{panic!("G8a held")};
    let IndependentTiltedVolumeResult::Ready(old_a)=IndependentTiltedVolumeOracle::build(
        &a,&b,&s,tol()).unwrap() else{panic!("G8b held")};
    let TiltedOneWallShellResult::Ready(shell)=TiltedOneWallShellPlan::build(
        &a,&b,&s,&old_a,tol()).unwrap() else{panic!("G8c held")};
    let IndependentTiltedBoundaryFluxResult::Ready(old_b)=
        IndependentTiltedBoundaryFluxOracle::build(&a,&b,&shell,tol()).unwrap()
        else{panic!("old boundary flux held")};
    let box_source=SourceSpec::new(SourceId::new(1).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),SourceGeometry::AxisAlignedBlock{
            min:p(-1.,1.,-3.),max:p(3.,4.,3.)},tol()).unwrap();
    let cyl_source=SourceSpec::new(SourceId::new(2).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),SourceGeometry::FiniteCylinder{
            frame,radius:1.5,axial:ParameterRange::bounded(0.,10.).unwrap()},tol()).unwrap();
    let snapshot=SourceSnapshot::new([box_source.clone(),cyl_source.clone()]).unwrap();
    let intent=snapshot.intent(BooleanKind::Difference,box_source.reference(),
        cyl_source.reference(),tol(),UnitPolicy::RequireIdentical).unwrap();
    let source=evaluate_tilted_one_wall(&snapshot,&intent).unwrap();
    let canonical=shell.checked_canonical().unwrap();
    let flux=evaluate_tilted_boundary_flux(&canonical,-1.,tol()).unwrap();
    let limits=source.interval();
    assert!((limits.lower-old_a.retained.lower).abs()<1e-7);
    assert!((limits.upper-old_a.retained.upper).abs()<1e-7);
    assert!((flux.result_volume()-old_b.witness.result_volume).abs()<1e-7,
        "new flux {} old {}",flux.result_volume(),old_b.witness.result_volume);
    let slack=source.numerical_guard()+flux.numerical_guard()+4.*tol().length_at(10.)*100.;
    assert!(flux.result_volume()>=limits.lower-slack);
    assert!(flux.result_volume()<=limits.upper+slack);
    assert_eq!(flux.harmonic_uses()+flux.generator_uses(),4);
    assert!(evaluate_tilted_boundary_flux(&canonical,-2.,tol()).is_err());
}
#[test]
fn x_axis_tilt_regression(){for angle in [0.04,-0.03,0.055]{run(angle,true);}}
#[test]
fn z_axis_tilt_regression(){run(0.04,false);}
