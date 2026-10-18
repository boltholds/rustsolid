use rustsolid::{
    AnalyticOneWallPlan,AnalyticOneWallResult,AnalyticOneWallHold,
    BooleanOperation,BrepBody,BrepOrigin,CylindricalBrep,
    FaceSense,GeometryTolerance,GraphOperand,Point3,Solid,Surface3,
};
use std::collections::{BTreeMap,BTreeSet};
fn q(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn box_body()->BrepBody{
    BrepBody::from(Solid::block(q(-1.,1.,-3.),4.,3.,6.).unwrap())
}
fn cutter()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(q(0.,0.,0.),1.5,5.,tol()).unwrap())
}
fn ready(a:&BrepBody,b:&BrepBody)->AnalyticOneWallPlan{
    let result=AnalyticOneWallPlan::build(
        a,b,BooleanOperation::Difference,tol()).unwrap();
    let AnalyticOneWallResult::Ready(plan)=result else{
        panic!("one-wall source shape Held: {result:?}");
    };
    plan
}
fn dist(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}
#[test]
fn one_wall_notch_is_exact_single_closed_curved_shell_with_two_planar_wall_patches(){
    let (a,b)=(box_body(),cutter());
    let before_a=format!("{a:?}");
    let before_b=format!("{b:?}");
    let plan=ready(&a,&b);
    let m=plan.brep();
    assert_eq!(m.origin,BrepOrigin::AnalyticOneWallCylinderDifference);
    assert_eq!(m.shells.len(),1);
    assert_eq!(m.vertices.len(),12);
    assert_eq!(m.edges.len(),18);
    assert_eq!(m.coedges.len(),36);
    assert_eq!(m.faces.len(),8);
    assert_eq!(m.loops.len(),8);
    assert_eq!(m.euler_characteristic(),2);
    assert_eq!(m.summary().genus,0);
    assert!(plan.closed_shell_proven&&!plan.can_apply);
    assert_eq!(plan.source_faces.len(),8);
    let matching=plan.source_faces.iter().filter(|x|
        x.source_operand==GraphOperand::First
            &&x.source_face==plan.source_faces[5].source_face
    ).collect::<Vec<_>>();
    assert_eq!(matching.len(),2,
        "two physically distinct material faces must share original x-wall source");
    assert_eq!(plan.radial.shell_components,1);
    assert_eq!(plan.radial.euler_characteristic,2);
    assert_eq!(plan.radial.total_genus,0);
    m.validate().unwrap();
    plan.validate_against(&a,&b,tol()).unwrap();
    assert_eq!(format!("{a:?}"),before_a);
    assert_eq!(format!("{b:?}"),before_b);
}
#[test]
fn major_circle_arc_and_generators_have_exact_twin_physical_coedges_and_uv(){
    let (a,b)=(box_body(),cutter());
    let plan=ready(&a,&b);
    let m=plan.brep();
    assert!(plan.component.major_arc_parameter_span>std::f64::consts::PI);
    assert!(plan.component.major_arc_parameter_span<std::f64::consts::TAU);
    assert_eq!(m.face_sense(plan.component.faces[7]),Some(FaceSense::Reversed));
    let mut cuts=BTreeSet::new();
    let mut events=BTreeMap::new();
    for (i,&arc) in plan.component.new_cut_arcs.iter().enumerate(){
        assert!(cuts.insert(arc));
        let eid=plan.component.edges[[4,10,16,17][i]];
        let e=&m.edges[eid.0 as usize];
        let [a,b]=e.coedges;
        let lhs=&m.coedges[a.0 as usize];
        let rhs=&m.coedges[b.0 as usize];
        assert_ne!(lhs.face,rhs.face);
        assert_ne!(lhs.reversed,rhs.reversed);
        assert_eq!(lhs.twin,b);
        assert_eq!(rhs.twin,a);
        for (event,vertex) in [
            (plan.component.new_cut_events[i][0],e.start),
            (plan.component.new_cut_events[i][1],e.end),
        ]{
            if let Some(before)=events.insert(event,vertex){
                assert_eq!(before,vertex);
            }
        }
        for &cid in &e.coedges{
            let edge=m.edge_curve(eid).unwrap();
            let co=m.coedges[cid.0 as usize].clone();
            let trim=m.coedge_pcurve(cid).unwrap();
            let surface=m.face_surface(co.face).unwrap();
            let rustsolid::ParameterRange::Bounded{start:e0,end:e1}=
                m.geometry.edges[eid.0 as usize].domain else{panic!("3D unbounded")};
            let rustsolid::ParameterRange::Bounded{start:t0,end:t1}=
                m.geometry.coedges[cid.0 as usize].domain else{panic!("UV unbounded")};
            for t in [0.,0.05,0.25,0.5,0.75,0.95,1.]{
                let uv=trim.evaluate(t0+(t1-t0)*t).unwrap();
                let from_surface=surface.evaluate(uv.u,uv.v).unwrap().point;
                let frac=if co.reversed{1.-t}else{t};
                let from_edge=edge.evaluate(e0+(e1-e0)*frac).unwrap();
                assert!(dist(from_surface,from_edge)<1e-7);
            }
        }
    }
    assert_eq!(cuts.len(),4);
    assert_eq!(cuts.len(),plan.edge_splits.registry.arcs.len());
    assert_eq!(events.len(),plan.edge_splits.registry.physical_event_count);
    assert_eq!(plan.source_faces.iter().filter(|f|
        f.reversed_for_difference
            &&matches!(m.face_surface(f.result_face),Some(Surface3::Cylinder(_)))
    ).count(),1);
}
#[test]
fn every_bounded_cap_and_split_plane_has_strict_inside_outside_material_samples(){
    let plan=ready(&box_body(),&cutter());
    let m=plan.brep();
    assert_eq!(plan.material_witnesses.len(),8);
    let bbox=box_body().shared().unwrap().bbox;
    let material=|p:Point3|->bool{
        p.x>bbox.min.x&&p.x<bbox.max.x
            &&p.y>bbox.min.y&&p.y<bbox.max.y
            &&p.z>bbox.min.z&&p.z<bbox.max.z
            &&p.x.hypot(p.z)>1.5
    };
    let mut seen=BTreeSet::new();
    for w in &plan.material_witnesses{
        assert!(seen.insert(w.face));
        assert!(material(w.retained_point),"{:?}",w.face);
        assert!(!material(w.removed_point));
    }
    for c in &m.faces{
        assert!(seen.contains(&c.id));
    }
}
#[test]
fn mass_is_source_exact_and_independently_matched_to_numerical_circle_clipping(){
    let plan=ready(&box_body(),&cutter());
    let r=1.5_f64;
    let n=65536usize;
    let a=-1.;let b=r;
    let dx=(b-a)/n as f64;
    let mut sum=0.;
    let mut moment=0.;
    for i in 0..=n{
        let x=a+dx*(i as f64);
        let w=2.*((r-x)*(r+x)).max(0.).sqrt();
        let factor=if i==0||i==n{1.}else if i%2==0{2.}else{4.};
        sum+=factor*w;
        moment+=factor*x*w;
    }
    let disk_area=sum*dx/3.;
    let disk_x=moment*dx/3.;
    let expected_area=24.-disk_area;
    let expected_volume=expected_area*3.;
    let expected_x=(24.*1.-disk_x)/expected_area;
    assert!((plan.component.exact_area_xz-expected_area).abs()<1e-6);
    assert!((plan.component.exact_volume-expected_volume).abs()<3e-6);
    assert!((plan.brep().mass.volume-expected_volume).abs()<3e-6);
    assert!((plan.brep().mass.centroid.x-expected_x).abs()<1e-6);
    assert!((plan.brep().mass.centroid.z).abs()<1e-9);
    assert!((plan.brep().mass.centroid.y-2.5).abs()<1e-9);
}
#[test]
fn translated_and_offcenter_one_wall_notch_remains_source_owned(){
    let a=BrepBody::from(Solid::block(q(4.,2.,-6.),5.,4.,8.).unwrap());
    let b=BrepBody::from(CylindricalBrep::upright(
        q(5.2,0.,-2.5),1.8,8.,tol()).unwrap());
    let result=ready(&a,&b);
    assert_eq!(result.brep().shells.len(),1);
    assert_eq!(result.component.edges.len(),18);
    assert!(result.component.exact_volume>0.);
    assert_eq!(result.brep().summary().genus,0);
    result.validate_against(&a,&b,tol()).unwrap();
}
#[test]
fn topology_and_material_orientation_corruption_rejected_prior_to_any_commit(){
    let (a,b)=(box_body(),cutter());
    let plan=ready(&a,&b);
    let mut corrupted=plan.brep().clone();
    corrupted.faces[7].sense=FaceSense::Forward;
    assert!(corrupted.validate().is_err());
    let mut corrupted=plan.brep().clone();
    let eid=plan.component.edges[4].0 as usize;
    let twin=corrupted.edges[eid].coedges[1];
    corrupted.edges[eid].coedges[0]=twin;
    assert!(corrupted.validate().is_err());
    let mut source_changed=a.clone();
    if let BrepBody::Polyhedral(ref mut body)=source_changed {
        body.revision+=1;
    }
    assert!(plan.validate_against(&source_changed,&b,tol()).is_err());
    let mut altered=plan.clone();
    altered.component.major_arc_parameter_span-=0.1;
    assert!(altered.validate_against(&a,&b,tol()).is_err());
}
#[test]
fn unsupported_tangent_two_wall_four_corner_and_union_are_held(){
    let a=box_body();
    for r in [1.0,3.0,3.2,0.5]{
        let b=BrepBody::from(CylindricalBrep::upright(
            q(0.,0.,0.),r,5.,tol()).unwrap());
        assert!(matches!(AnalyticOneWallPlan::build(
            &a,&b,BooleanOperation::Difference,tol()).unwrap(),
            AnalyticOneWallResult::Held(_)
        ));
    }
    assert!(matches!(AnalyticOneWallPlan::build(
        &a,&cutter(),BooleanOperation::Union,tol()).unwrap(),
        AnalyticOneWallResult::Held(AnalyticOneWallHold::UnsupportedOperation)
    ));
}
