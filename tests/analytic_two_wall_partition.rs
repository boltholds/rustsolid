use rustsolid::{
    AnalyticTwoWallHold,AnalyticTwoWallPartitionPlan,
    AnalyticTwoWallPartitionResult,BrepBody,BrepOrigin,BooleanOperation,
    CylindricalBrep,FaceSense,GeometryTolerance,GraphOperand,
    Point3,Solid,
};
use std::collections::{BTreeMap,BTreeSet};
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tolerance()->GeometryTolerance{GeometryTolerance::default()}
fn block()->BrepBody{
    BrepBody::from(Solid::block(q(-1.2,1.,-3.),2.4,3.,6.).unwrap())
}
fn cylinder()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(q(0.,0.,0.),1.5,5.,tolerance()).unwrap())
}
fn ready(a:&BrepBody,b:&BrepBody)->AnalyticTwoWallPartitionPlan{
    let r=AnalyticTwoWallPartitionPlan::build(
        a,b,BooleanOperation::Difference,tolerance()).unwrap();
    let AnalyticTwoWallPartitionResult::Ready(plan)=r else{
        panic!("strict two-wall difference unexpectedly Held: {r:?}");
    };
    plan
}
fn norm(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

#[test]
fn transverse_two_wall_cutter_produces_two_real_closed_genus_zero_shells(){
    let (a,b)=(block(),cylinder());
    let unchanged_a=format!("{a:?}");
    let unchanged_b=format!("{b:?}");
    let p=ready(&a,&b);
    let m=p.brep();
    assert!(p.closed_shell_proven&&!p.can_apply);
    assert_eq!(m.origin,BrepOrigin::AnalyticTwoWallCylinderDifference);
    assert_eq!(m.shells.len(),2);
    assert_eq!(m.faces.len(),12);
    assert_eq!(m.loops.len(),12);
    assert_eq!(m.vertices.len(),16);
    assert_eq!(m.edges.len(),24);
    assert_eq!(m.coedges.len(),48);
    assert_eq!(m.euler_characteristic(),4);
    assert_eq!(m.summary().genus,0);
    assert_eq!(p.radial.shell_components,2);
    assert_eq!(p.radial.total_genus,0);
    assert_eq!(p.radial.euler_characteristic,4);
    assert_eq!(p.components.len(),2);
    assert_eq!(p.source_faces.len(),12);
    assert_eq!(p.material_sides.len(),12);
    m.validate().unwrap();
    p.validate_against(&a,&b,tolerance()).unwrap();
    assert_eq!(format!("{a:?}"),unchanged_a);
    assert_eq!(format!("{b:?}"),unchanged_b);
}

#[test]
fn every_curved_boundary_has_opposite_real_coedges_and_two_original_uv_carriers(){
    let (a,b)=(block(),cylinder());
    let p=ready(&a,&b);
    let m=p.brep();
    let mut arcs=BTreeSet::new();
    let mut events=BTreeMap::new();
    for c in &p.components {
        assert_eq!(m.shells[c.shell.0 as usize].faces,c.faces);
        assert_eq!(c.faces.len(),6);
        let senses=c.faces.iter().map(|fid|m.faces[fid.0 as usize].sense).collect::<Vec<_>>();
        assert_eq!(senses,vec![FaceSense::Forward;5]
            .into_iter().chain([FaceSense::Reversed]).collect::<Vec<_>>());
        for (index,&source) in c.global_arcs.iter().enumerate(){
            assert!(arcs.insert(source));
            let edge=c.edges[[1,5,9,10][index]];
            let e=&m.edges[edge.0 as usize];
            let [u,v]=e.coedges;
            let first=&m.coedges[u.0 as usize];
            let second=&m.coedges[v.0 as usize];
            assert_ne!(first.face,second.face);
            assert_ne!(first.reversed,second.reversed);
            assert_eq!(first.twin,v);
            assert_eq!(second.twin,u);
            for (event,vertex) in [
                (c.global_endpoints[index][0],e.start),
                (c.global_endpoints[index][1],e.end),
            ]{
                if let Some(previous)=events.insert(event,vertex) {
                    assert_eq!(previous,vertex);
                }
            }
        }
        for id in &c.edges {
            let e=&m.edges[id.0 as usize];
            for &cid in &e.coedges {
                let co=&m.coedges[cid.0 as usize];
                let edge=m.edge_curve(*id).unwrap();
                let trim=m.coedge_pcurve(cid).unwrap();
                let support=m.face_surface(co.face).unwrap();
                let rustsolid::ParameterRange::Bounded {start:lo,end:hi}=
                    m.geometry.edges[id.0 as usize].domain
                else{panic!("unbounded Edge3")};
                let rustsolid::ParameterRange::Bounded {start:tlo,end:thi}=
                    m.geometry.coedges[cid.0 as usize].domain
                else{panic!("unbounded trim")};
                for f in [0.0,0.125,0.5,0.875,1.0]{
                    let uv=trim.evaluate(tlo+(thi-tlo)*f).unwrap();
                    let xyz=support.evaluate(uv.u,uv.v).unwrap().point;
                    let ef=if co.reversed{1.-f}else{f};
                    assert!(norm(xyz,edge.evaluate(lo+(hi-lo)*ef).unwrap())<1e-7);
                }
            }
        }
    }
    assert_eq!(arcs.len(),8);
    assert_eq!(arcs.len(),p.split_sources.registry.arcs.len());
    assert_eq!(events.len(),p.split_sources.registry.physical_event_count);
    assert_eq!(p.source_faces.iter().filter(|x|
        x.operand==GraphOperand::Second&&x.reversed_for_difference).count(),2);
}

#[test]
fn every_retained_wall_and_curved_patch_has_strict_material_side_witnesses(){
    let plan=ready(&block(),&cylinder());
    let b=block().shared().unwrap().bbox;
    let mut seen=BTreeSet::new();
    let inside=|v:Point3|->bool{
        v.x>b.min.x&&v.x<b.max.x&&v.y>b.min.y&&v.y<b.max.y
            &&v.z>b.min.z&&v.z<b.max.z&&v.x.hypot(v.z)>1.5
    };
    for witness in &plan.material_sides{
        assert!(seen.insert(witness.face));
        assert!(inside(witness.retained));
        assert!(!inside(witness.removed));
    }
    assert_eq!(seen.len(),12);
    let mut tampered=plan.clone();
    tampered.material_sides[0].removed.x+=0.5;
    assert!(tampered.validate_against(&block(),&cylinder(),tolerance()).is_err());
}

#[test]
fn source_specific_result_volume_matches_independent_numerical_integral(){
    let plan=ready(&block(),&cylinder());
    // A Simpson computation independently of the builder antiderivative.
    let n=32768usize;
    let dx=2.4/n as f64;
    let r=1.5_f64;
    let mut sum=0.;
    for i in 0..=n {
        let x=-1.2+i as f64*dx;
        let weight=if i==0||i==n{1.}else if i%2==0{2.}else{4.};
        sum+=weight*(3.-(r*r-x*x).sqrt());
    }
    let expected=sum*dx/3.*3.;
    assert!(expected>0.);
    for part in &plan.components{
        assert!((part.analytic_volume-expected).abs()<1e-9);
    }
    assert!((plan.brep().mass.volume-2.*expected).abs()<2e-9);
    assert!(plan.brep().mass.centroid.x.abs()<1e-10);
    assert!((plan.brep().mass.centroid.y-2.5).abs()<1e-10);
    assert!(plan.brep().mass.centroid.z.abs()<1e-10);
}
#[test]
fn shifted_asymmetric_rectangular_box_still_has_real_analytic_components(){
    let a=BrepBody::from(Solid::block(q(1.,2.,-5.),3.,4.,9.).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        q(2.4,0.,-1.),1.8,8.,tolerance()).unwrap());
    let plan=ready(&a,&b);
    assert_eq!(plan.brep().shells.len(),2);
    assert_ne!(plan.components[0].analytic_volume,plan.components[1].analytic_volume);
    plan.brep().validate().unwrap();
    plan.validate_against(&a,&b,tolerance()).unwrap();
}
#[test]
fn manifold_failure_and_source_revision_are_rejected_without_partial_result(){
    let (a,b)=(block(),cylinder());
    let plan=ready(&a,&b);
    let mut m=plan.brep().clone();
    m.shells[1].faces[0]=m.shells[0].faces[0];
    assert!(m.validate().is_err());
    let mut m=plan.brep().clone();
    let id=plan.components[0].faces[5];
    m.faces[id.0 as usize].sense=FaceSense::Forward;
    assert!(m.validate().is_err());
    let mut m=plan.brep().clone();
    m.shells.pop();
    assert!(m.validate().is_err());
    let mut a2=a.clone();
    if let BrepBody::Polyhedral(ref mut solid)=a2 {
        solid.revision+=1;
    }
    assert!(plan.validate_against(&a2,&b,tolerance()).is_err());
    let b2=BrepBody::from(CylindricalBrep::upright(
        q(0.1,0.,0.),1.5,5.,tolerance()).unwrap());
    assert!(plan.validate_against(&a,&b2,tolerance()).is_err());
}
#[test]
fn tangent_opposite_wall_wrong_radius_or_operation_stays_held(){
    let a=block();
    for r in [1.2,3.0,3.1,0.5]{
        let b=BrepBody::from(CylindricalBrep::upright(
            q(0.,0.,0.),r,5.,tolerance()).unwrap());
        assert!(matches!(
            AnalyticTwoWallPartitionPlan::build(
                &a,&b,BooleanOperation::Difference,tolerance()).unwrap(),
            AnalyticTwoWallPartitionResult::Held(_)
        ));
    }
    let b=cylinder();
    assert!(matches!(
        AnalyticTwoWallPartitionPlan::build(
            &a,&b,BooleanOperation::Union,tolerance()).unwrap(),
        AnalyticTwoWallPartitionResult::Held(AnalyticTwoWallHold::UnsupportedOperation)
    ));
}
