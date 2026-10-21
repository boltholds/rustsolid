//! Independent v2 source-only G8b witness cross-check against the historical
//! G8b regression family. This test is migration evidence, not a production
//! dependency from rustsolid-validation to the old root modeling kernel.

use rustsolid as old;
use rustsolid_geometry::{Frame3, ModelUnits, ParameterRange};
use rustsolid_source::{BooleanIntent, BooleanKind, SourceGeometry, SourceId,
    SourceRevision, SourceSnapshot, SourceSpec, UnitPolicy};
use rustsolid_validation::{g8_source_volume, G8SourceVolumeOutcome};

fn p(x:f64,y:f64,z:f64)->old::Point3{old::Point3{x,y,z}}
fn tol()->old::GeometryTolerance{old::GeometryTolerance::default()}

fn pair(theta:f64,yz:bool)->(
    old::BrepBody,old::BrepBody,
    SourceSnapshot,BooleanIntent
){
    let block=old::Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap();
    let frame=if yz {
        Frame3::from_axes(p(0.,-2.,0.),p(1.,0.,0.),
            p(0.,-theta.sin(),-theta.cos()),tol()).unwrap()
    }else{
        Frame3::from_axes(p(0.,-2.,0.),
            p(theta.cos(),-theta.sin(),0.),p(0.,0.,-1.),tol()).unwrap()
    };
    let cutter=old::CylindricalBrep::new(frame,1.5,10.,tol()).unwrap();
    let sources=[
        SourceSpec::new(SourceId::new(11).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),
            SourceGeometry::AxisAlignedBlock{min:block.bbox.min,max:block.bbox.max},
            tol()).unwrap(),
        SourceSpec::new(SourceId::new(17).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),
            SourceGeometry::FiniteCylinder{
                frame,radius:1.5,axial:ParameterRange::bounded(0.,10.).unwrap(),
            },
            tol()).unwrap(),
    ];
    let snap=SourceSnapshot::new(sources.clone()).unwrap();
    let intent=snap.intent(BooleanKind::Difference,
        sources[0].reference(),sources[1].reference(),
        tol(),UnitPolicy::RequireIdentical).unwrap();
    (block.into(),cutter.into(),snap,intent)
}

#[test]
fn g8_source_only_brackets_agree_with_historical_g8b_math(){
    for (theta,yz) in [(0.04,false),(-0.03,false),(0.055,false),(0.04,true)] {
        let (a,b,snap,intent)=pair(theta,yz);
        let old::TiltedOneWallSourceResult::Ready(src)=
            old::TiltedOneWallSourcePlan::build(&a,&b,old::BooleanOperation::Difference,tol()).unwrap()
            else {panic!("source geometry G8a held on pinned case theta={theta}, yz={yz}")};
        let old::IndependentTiltedVolumeResult::Ready(reference)=
            old::IndependentTiltedVolumeOracle::build(&a,&b,&src,tol()).unwrap()
            else {panic!("baseline G8b held")};
        let G8SourceVolumeOutcome::Ready(next)=g8_source_volume(&snap,&intent).unwrap()
            else {panic!("new standalone independent G8b held")};
        let max_diff=(next.retained.lower-reference.retained.lower).abs()
            .max((next.retained.upper-reference.retained.upper).abs());
        assert!(max_diff<5e-7,
            "theta={theta}, yz={yz}, baseline={:?}, next={:?}",
            reference.retained,next.retained);
        assert_eq!(next.angular_panels,reference.angular_panels);
        assert_eq!(next.height_panels,reference.height_panels);
        assert!(next.independent_cross_check>=next.removed_inside_box[0]
            -next.cross_check_error-next.numerical_guard);
        assert!(next.independent_cross_check<=next.removed_inside_box[1]
            +next.cross_check_error+next.numerical_guard);
        assert_eq!(next.source_stamp(),snap.stamp());
    }
}
