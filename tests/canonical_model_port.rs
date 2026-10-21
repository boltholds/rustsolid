//! Transitional regression-only qualification of the new canonical B-rep.
//! This adapter must NOT become a public compatibility layer: migration should
//! construct topology::BrepDraft directly from source/imprint/assembly contracts.

use rustsolid as old;
use rustsolid_geometry::{Curve2, Curve3, ParameterRange, Surface3};
use rustsolid_topology as topo;

fn interval(range:ParameterRange)->(f64,f64){
    match range {
        ParameterRange::Bounded{start,end} if start.is_finite() && end.is_finite() && end>start=>
            (start,end),
        _=>panic!("source requires bounded edge/trim parameter"),
    }
}

fn conversion_for(
    model:&old::BrepModel, coedge:&old::Coedge,
    edge:&old::Edge, new:&topo::BrepEdge,
)->topo::AffineParameterMap{
    let id=edge.id.0 as usize;
    let cid=coedge.id.0 as usize;
    let (t0,t1)=interval(model.geometry.edges[id].domain);
    let (u0,u1)=interval(model.geometry.coedges[cid].domain);
    let curve=model.edge_curve(edge.id).expect("existing 3D carrier");
    let chart=model.coedge_pcurve(coedge.id).expect("existing UV trim");
    let surface=model.face_surface(coedge.face).expect("existing surface");
    let slope=(u1-u0)/(t1-t0);
    let tries=[
        topo::AffineParameterMap{scale:slope,offset:u0-slope*t0},
        topo::AffineParameterMap{scale:-slope,offset:u1+slope*t0},
    ];
    let mut candidates=tries.map(|map|{
        let mut peak=0.0f64;
        for fraction in [0.0,0.25,0.5,0.75,1.0]{
            let t=t0+(t1-t0)*fraction;
            let from_curve=curve.evaluate(t).unwrap();
            let uv=chart.evaluate(map.map(t).unwrap()).unwrap();
            let from_surface=surface.evaluate(uv.u,uv.v).unwrap().point;
            let diff=from_curve.sub(from_surface);
            peak=peak.max(diff.x.hypot(diff.y).hypot(diff.z));
        }
        (map,peak)
    });
    candidates.sort_by(|a,b|a.1.total_cmp(&b.1));
    assert!(candidates[0].1<1e-7,
        "source coedge {} edge {} cannot be represented by affine pcurve mapping: {:?}",
        cid,id,candidates);
    assert_eq!(new.id,edge.id);
    candidates[0].0
}

fn test_only_from_transitional(model:&old::BrepModel)->topo::BrepModel {
    model.validate().expect("original analytic source B-rep");
    let vertices=model.vertices.iter().map(|v|topo::BrepVertex{
        id:v.id,point:v.position,
    }).collect();
    let edges=model.edges.iter().map(|e|{
        topo::BrepEdge{
            id:e.id,vertices:[e.start,e.end],
            curve:topo::Curve3Id(model.geometry.edges[e.id.0 as usize].curve.0),
            domain:model.geometry.edges[e.id.0 as usize].domain,
            uses:e.coedges,
            kind:if model.coedges[e.coedges[0].0 as usize].face==
                    model.coedges[e.coedges[1].0 as usize].face {
                topo::EdgeKind::PeriodicSeam
            }else{topo::EdgeKind::Physical},
        }
    }).collect::<Vec<_>>();
    let coedges=model.coedges.iter().map(|c|{
        let edge=&model.edges[c.edge.0 as usize];
        let map=conversion_for(model,c,edge,&edges[c.edge.0 as usize]);
        topo::BrepCoedge{
            id:c.id,edge:c.edge,face:c.face,loop_id:c.loop_id,
            pcurve:topo::Curve2Id(model.geometry.coedges[c.id.0 as usize].pcurve.0),
            edge_to_pcurve:map,
            reversed:c.reversed,next:c.next,prev:c.prev,twin:c.twin,
        }
    }).collect();
    let loops=model.loops.iter().map(|l|topo::BrepLoop{
        id:l.id,face:l.face,role:l.role,first:l.first_coedge,
    }).collect();
    let faces=model.faces.iter().map(|f|topo::BrepFace{
        id:f.id,surface:topo::SurfaceId(f.geometry.surface.0),loops:f.loops.clone(),
        sense:match f.sense{old::FaceSense::Forward=>topo::FaceSense::Forward,
                           old::FaceSense::Reversed=>topo::FaceSense::Reversed},
    }).collect();
    let shells=model.shells.iter().map(|s|topo::BrepShell{
        id:s.id,faces:s.faces.clone(),sense:topo::ShellSense::Outward,
    }).collect();
    let carriers=topo::GeometryCatalog{
        curves3:model.geometry.curves3.clone(),
        curves2:model.geometry.curves2.clone(),
        surfaces:model.geometry.surfaces.clone(),
    };
    topo::BrepDraft{
        incarnation:topo::LocalBrepIncarnation::new(17).unwrap(),
        units:model.units,tolerance:model.tolerance,
        vertices,edges,coedges,loops,faces,shells,carriers
    }.check().expect("v2 canonical topology and 3D/UV agreement")
}

fn p(x:f64,y:f64,z:f64)->old::Point3{old::Point3{x,y,z}}
fn tol()->old::GeometryTolerance{old::GeometryTolerance::default()}

#[test]
fn block_and_closed_cylinder_retain_source_carrier_topology() {
    let box_brep=old::BrepBody::from(old::Solid::block(p(0.,0.,0.),3.,4.,5.).unwrap())
        .shared().unwrap();
    let box_model=test_only_from_transitional(&box_brep);
    assert_eq!((box_model.faces().len(),box_model.edges().len()),(6,12));
    assert!(box_model.carriers().curves3.iter().all(|c|matches!(c,Curve3::Line(_))));

    let cylinder=old::CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap();
    let old_model=old::BrepBody::from(cylinder).shared().unwrap();
    let new=test_only_from_transitional(&old_model);
    assert_eq!((new.faces().len(),new.edges().len(),new.coedges().len()),(3,3,6));
    assert_eq!(new.edges().iter().filter(|e|e.kind==topo::EdgeKind::PeriodicSeam).count(),1);
    assert_eq!(new.edges().iter().filter(|e|e.vertices[0]==e.vertices[1]).count(),2);
    assert_eq!(new.carriers().surfaces.iter().filter(|s|matches!(s,Surface3::Cylinder(_))).count(),1);
}

#[test]
fn g7_one_wall_produces_valid_meshless_v2_analytic_brep() {
    let a=old::BrepBody::from(old::Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
    let b=old::BrepBody::from(old::CylindricalBrep::upright(
        p(0.,0.,0.),1.5,5.,tol()).unwrap());
    let old::AnalyticOneWallResult::Ready(plan)=
        old::AnalyticOneWallPlan::build(&a,&b,old::BooleanOperation::Difference,tol()).unwrap()
    else {panic!("G7 source is not supported")};
    let original=plan.brep();
    let checked=test_only_from_transitional(original);
    assert_eq!((checked.vertices().len(),checked.edges().len(),checked.coedges().len(),
        checked.faces().len()),(12,18,36,8));
    assert_eq!(checked.carriers().curves3,original.geometry.curves3);
}

#[test]
fn g8_tilted_one_wall_produces_valid_meshless_v2_analytic_brep() {
    let a=old::BrepBody::from(old::Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
    let theta=0.04_f64;
    let frame=old::Frame3::from_axes(p(0.,-2.,0.),
        p(theta.cos(),-theta.sin(),0.),p(0.,0.,-1.),tol()).unwrap();
    let b=old::BrepBody::from(old::CylindricalBrep::new(frame,1.5,10.,tol()).unwrap());
    let old::TiltedOneWallSourceResult::Ready(source)=
        old::TiltedOneWallSourcePlan::build(&a,&b,old::BooleanOperation::Difference,tol()).unwrap()
    else{panic!("G8a held")};
    let old::IndependentTiltedVolumeResult::Ready(volume)=
        old::IndependentTiltedVolumeOracle::build(&a,&b,&source,tol()).unwrap()
    else{panic!("G8b held")};
    let old::TiltedOneWallShellResult::Ready(shell)=
        old::TiltedOneWallShellPlan::build(&a,&b,&source,&volume,tol()).unwrap()
    else{panic!("G8c held")};
    let actual=shell.brep();
    let checked=test_only_from_transitional(actual);
    assert_eq!((checked.vertices().len(),checked.edges().len(),checked.coedges().len(),
        checked.faces().len()),(12,18,36,8));
    assert_eq!(checked.carriers().curves3,actual.geometry.curves3);
}


/// Compares G8c raw flux from the ACTUAL canonical coedges against the old
/// independently computed G8c witness and the new source-only G8b enclosure.
/// This is a migration qualification, not a production B-rep adapter.
#[test]
fn g8_independent_source_and_actual_result_flux_agree_on_x_and_yz_tilts() {
    use rustsolid_geometry::{Frame3,ModelUnits,ParameterRange};
    use rustsolid_source::{
        BooleanKind,SourceGeometry,SourceId,SourceRevision,SourceSnapshot,
        SourceSpec,UnitPolicy,
    };
    use rustsolid_validation::{
        compare_g8_volume,g8_boundary_flux,G8ComparisonOutcome,G8FluxOutcome,
    };
    for (theta,yz) in [(0.04_f64,false),(-0.03,false),(0.04,true)] {
        let a=old::BrepBody::from(old::Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap());
        let frame=if yz {
            Frame3::from_axes(p(0.,-2.,0.),p(1.,0.,0.),
                p(0.,-theta.sin(),-theta.cos()),tol()).unwrap()
        }else{
            Frame3::from_axes(p(0.,-2.,0.),
                p(theta.cos(),-theta.sin(),0.),p(0.,0.,-1.),tol()).unwrap()
        };
        let b=old::BrepBody::from(old::CylindricalBrep::new(frame,1.5,10.,tol()).unwrap());
        let old::TiltedOneWallSourceResult::Ready(source)=
            old::TiltedOneWallSourcePlan::build(&a,&b,old::BooleanOperation::Difference,tol()).unwrap()
            else{panic!("old G8a source Held")};
        let old::IndependentTiltedVolumeResult::Ready(vol)=
            old::IndependentTiltedVolumeOracle::build(&a,&b,&source,tol()).unwrap()
            else{panic!("old G8b source Held")};
        let old::TiltedOneWallShellResult::Ready(shell)=
            old::TiltedOneWallShellPlan::build(&a,&b,&source,&vol,tol()).unwrap()
            else{panic!("old G8c shell Held")};
        let old::IndependentTiltedBoundaryFluxResult::Ready(baseline)=
            old::IndependentTiltedBoundaryFluxOracle::build(&a,&b,&shell,tol()).unwrap()
            else{panic!("old G8c flux Held")};
        let canonical=test_only_from_transitional(shell.brep());
        let G8FluxOutcome::Ready(actual)=g8_boundary_flux(&canonical)
            else{panic!("new actual result flux Held for theta={theta}, yz={yz}")};
        assert!(
            (actual.retained_volume-baseline.witness.result_volume).abs()<2e-7,
            "independent actual flux disagreement for theta={theta}, yz={yz}: {:?}",
            actual,
        );
        assert_eq!(actual.harmonic_integrated_uses,
            baseline.witness.harmonic_integrated_uses);
        assert_eq!(actual.axial_generator_uses,
            baseline.witness.vertical_generator_uses);
        let spec_a=SourceSpec::new(
            SourceId::new(101).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),
            SourceGeometry::AxisAlignedBlock{min:p(-1.,1.,-3.),max:p(3.,4.,3.)},
            tol(),
        ).unwrap();
        let spec_b=SourceSpec::new(
            SourceId::new(102).unwrap(),SourceRevision::new(1).unwrap(),
            ModelUnits::default(),
            SourceGeometry::FiniteCylinder{
                frame,radius:1.5,axial:ParameterRange::bounded(0.,10.).unwrap(),
            },tol(),
        ).unwrap();
        let snapshot=SourceSnapshot::new([spec_a.clone(),spec_b.clone()]).unwrap();
        let intent=snapshot.intent(BooleanKind::Difference,
            spec_a.reference(),spec_b.reference(),tol(),UnitPolicy::RequireIdentical).unwrap();
        let G8ComparisonOutcome::VolumeCompatible(agreement)=
            compare_g8_volume(&snapshot,&intent,&canonical).unwrap()
            else{panic!("new source/actual boundary independent comparison Held")};
        assert!(agreement.source_interval[0]<=agreement.actual_boundary_volume
            + agreement.numerical_slack);
        assert!(agreement.source_interval[1]>=agreement.actual_boundary_volume
            - agreement.numerical_slack);
        assert_eq!(agreement.assurance_level,"GuardedF64");
        assert_eq!(agreement.model_incarnation,canonical.incarnation().get());
        assert_eq!(&agreement.source_stamp,snapshot.stamp());
    }
}
