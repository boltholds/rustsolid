//! Production G7 canonical constructor regression against prior source-certified
//! G7 analytical shell and INDEPENDENT source vs actual boundary flux.
//! The new producer does not import the old BrepBody or old PatchBuilder.

use rustsolid as old;
use rustsolid_construct::{g7_one_wall_direct,G7BuildOutcome,G7BuildHold};
use rustsolid_geometry::{Frame3,ModelUnits,ParameterRange};
use rustsolid_source::{
    BooleanKind,SourceGeometry,SourceId,SourceRevision,SourceSnapshot,
    SourceSpec,UnitPolicy,
};
use rustsolid_validation::{compare_g8_volume,g8_boundary_flux,G8ComparisonOutcome,G8FluxOutcome};

fn p(x:f64,y:f64,z:f64)->old::Point3{old::Point3{x,y,z}}
fn tol()->old::GeometryTolerance{old::GeometryTolerance::default()}
fn pair()->(SourceSnapshot,rustsolid_source::BooleanIntent){
    let a=SourceSpec::new(SourceId::new(1).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),
        SourceGeometry::AxisAlignedBlock{min:p(-1.,1.,-3.),max:p(3.,4.,3.)},
        tol()).unwrap();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(1.,0.,0.),p(0.,0.,-1.),tol()).unwrap();
    let b=SourceSpec::new(SourceId::new(2).unwrap(),SourceRevision::new(1).unwrap(),
        ModelUnits::default(),
        SourceGeometry::FiniteCylinder{frame,radius:1.5,
            axial:ParameterRange::bounded(0.,5.).unwrap()},
        tol()).unwrap();
    let snapshot=SourceSnapshot::new([a.clone(),b.clone()]).unwrap();
    let intent=snapshot.intent(BooleanKind::Difference,
        a.reference(),b.reference(),tol(),UnitPolicy::RequireIdentical).unwrap();
    (snapshot,intent)
}
#[test]
fn source_only_producer_directly_builds_same_g7_brep_and_independent_volume(){
    let (snapshot,intent)=pair();
    let G7BuildOutcome::Staged(staged)=g7_one_wall_direct(&snapshot,&intent).unwrap()
        else{panic!("direct source-bound constructor held")};
    let new=staged.model();
    assert_eq!((new.vertices().len(),new.edges().len(),new.faces().len()),(12,18,8));
    new.validate().unwrap();
    let a=old::BrepBody::from(old::Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
    let b=old::BrepBody::from(old::CylindricalBrep::upright(
        p(0.,0.,0.),1.5,5.,tol()).unwrap());
    let old::AnalyticOneWallResult::Ready(baseline)=
        old::AnalyticOneWallPlan::build(&a,&b,old::BooleanOperation::Difference,tol()).unwrap()
        else{panic!("reference G7 producer unexpectedly held")};
    let old_model=baseline.brep();
    assert_eq!(new.edges().len(),old_model.edges.len());
    assert_eq!(new.faces().len(),old_model.faces.len());
    assert_eq!(new.coedges().len(),old_model.coedges.len());
    let G8FluxOutcome::Ready(actual)=g8_boundary_flux(new)
        else{panic!("G7 direct coedges rejected by independent flux oracle")};
    let G8ComparisonOutcome::VolumeCompatible(report)=
        compare_g8_volume(&snapshot,&intent,new).unwrap()
        else{panic!("G7 source/flux independent comparisons incompatible")};
    assert!((actual.retained_volume-old_model.mass.volume).abs()<2e-7,
        "direct G7 flux {} vs old material volume {}",
        actual.retained_volume,old_model.mass.volume);
    assert!((report.actual_boundary_volume-old_model.mass.volume).abs()<2e-7);
    assert_eq!(report.assurance_level,"GuardedF64");
    assert_eq!(staged.source_stamp(),snapshot.stamp());
    assert_eq!(staged.intent_stamp(),intent.intent_stamp());
}
#[test]
fn direct_g7_rejects_tangency_and_stale_source_before_any_topology_commit(){
    let (snapshot,intent)=pair();
    let cylinder=snapshot.get(intent.operands().1.id).unwrap();
    let SourceGeometry::FiniteCylinder{frame,axial,..}=cylinder.geometry()
        else{unreachable!()};
    let replaced=cylinder.successor(SourceGeometry::FiniteCylinder{
        frame:*frame,radius:1.0,axial:*axial,
    },cylinder.tolerance()).unwrap();
    let changed=snapshot.with_replacement(replaced).unwrap();
    assert!(g7_one_wall_direct(&changed,&intent).is_err());
    let new_intent=changed.intent(BooleanKind::Difference,
        changed.get(intent.operands().0.id).unwrap().reference(),
        changed.get(intent.operands().1.id).unwrap().reference(),
        tol(),UnitPolicy::RequireIdentical).unwrap();
    assert!(matches!(g7_one_wall_direct(&changed,&new_intent).unwrap(),
        G7BuildOutcome::Held(G7BuildHold::NotStrictOneWallNotch)));
}
