//! G4b2b3b: four REAL disconnected closed analytic B-rep corner shells.
use rustsolid::{
    AnalyticFacePatchPartitionPlan,AnalyticFacePatchPartitionResult,
    AnalyticFacePatchHold,BrepBody,BrepOrigin,BooleanOperation,
    CylindricalBrep,FaceSense,GeometryTolerance,
    GlobalEventId,GraphOperand,Point3,Solid,Surface3,
};
use std::collections::{BTreeMap,BTreeSet};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn cutter(origin:Point3,r:f64,height:f64)->BrepBody {
    BrepBody::from(CylindricalBrep::upright(origin,r,height,tol()).unwrap())
}
fn square()->BrepBody {
    BrepBody::from(Solid::block(p(-1.5,1.,-1.5),3.,3.,3.).unwrap())
}
fn cylinder()->BrepBody{cutter(p(0.,0.,0.),2.,5.)}
fn result(a:&BrepBody,b:&BrepBody)->AnalyticFacePatchPartitionResult {
    AnalyticFacePatchPartitionPlan::build(
        a,b,BooleanOperation::Difference,tol()).unwrap()
}
fn ready(a:&BrepBody,b:&BrepBody)->AnalyticFacePatchPartitionPlan{
    let outcome=result(a,b);
    let AnalyticFacePatchPartitionResult::Ready(plan)=outcome else {
        panic!("expected a real 4-component curved shell: {outcome:?}");
    };
    plan
}
fn distance(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

#[test]
fn strict_partial_box_cylinder_builds_four_closed_independent_genus_zero_shells(){
    let (a,b)=(square(),cylinder());
    let source_a=format!("{a:?}");
    let source_b=format!("{b:?}");
    let plan=ready(&a,&b);
    let brep=plan.brep();
    assert_eq!(brep.origin,BrepOrigin::AnalyticPartialBoxCylinderDifference);
    assert!(plan.closed_shell_proven&&!plan.can_apply);
    assert_eq!(brep.shells.len(),4);
    assert_eq!(plan.corners.len(),4);
    assert_eq!(plan.material_witnesses.len(),20);
    assert_eq!(brep.vertices.len(),24);
    assert_eq!(brep.edges.len(),36);
    assert_eq!(brep.coedges.len(),72);
    assert_eq!(brep.loops.len(),20);
    assert_eq!(brep.faces.len(),20);
    assert_eq!(brep.euler_characteristic(),8);
    assert_eq!(brep.summary().genus,0);
    assert_eq!(plan.radial.shell_components,4);
    assert_eq!(plan.radial.total_genus,0);
    assert_eq!(plan.radial.euler_characteristic,8);
    assert_eq!(plan.radial.vertices,24);
    assert_eq!(plan.radial.edges,36);
    assert_eq!(plan.radial.material_face_cells,20);
    brep.validate().unwrap();
    assert_eq!(format!("{a:?}"),source_a);
    assert_eq!(format!("{b:?}"),source_b);
    plan.validate_against(&a,&b,tol()).unwrap();
}

#[test]
fn true_cut_physical_edges_have_opposite_coedges_and_independent_uv_carriers(){
    let (a,b)=(square(),cylinder());
    let plan=ready(&a,&b);
    let model=plan.brep();
    let mut arcs=BTreeSet::new();
    let mut events=BTreeMap::<GlobalEventId,BTreeSet<u32>>::new();
    for corner in &plan.corners {
        assert_eq!(model.shells[corner.shell.0 as usize].faces.to_vec(),corner.faces);
        for (i,&arc) in corner.cut_arcs.iter().enumerate(){
            assert!(arcs.insert(arc));
            let e=corner.edges[[2,5,7,8][i]];
            let edge=&model.edges[e.0 as usize];
            let [first,second]=edge.coedges;
            let x=&model.coedges[first.0 as usize];
            let y=&model.coedges[second.0 as usize];
            assert_ne!(x.face,y.face);
            assert_ne!(x.reversed,y.reversed);
            assert_eq!(x.twin,y.id);
            assert_eq!(y.twin,x.id);
            assert!(model.face_surface(x.face).is_some());
            assert!(model.coedge_pcurve(x.id).is_some());
            assert!(model.coedge_pcurve(y.id).is_some());
            for (event,vertex) in [
                (corner.global_endpoints[i][0],edge.start),
                (corner.global_endpoints[i][1],edge.end),
            ]{
                events.entry(event).or_default().insert(vertex.0);
            }
            let [s,t]=edge.coedges;
            for co_id in [s,t] {
                let co=&model.coedges[co_id.0 as usize];
                let carrier=model.edge_curve(e).unwrap();
                let eg=model.geometry.edges[e.0 as usize];
                let pc=model.coedge_pcurve(co_id).unwrap();
                let uv_range=model.geometry.coedges[co_id.0 as usize].domain;
                let rustsolid::ParameterRange::Bounded{
                    start:tr0,end:tr1}=uv_range else{panic!("unbounded cut UV")};
                let rustsolid::ParameterRange::Bounded{
                    start:ed0,end:ed1}=eg.domain else{panic!("unbounded cut Curve3")};
                let surface=model.face_surface(co.face).unwrap();
                for f in [0.0,0.125,0.5,0.875,1.0] {
                    let uv=pc.evaluate(tr0+(tr1-tr0)*f).unwrap();
                    let on_surface=surface.evaluate(uv.u,uv.v).unwrap().point;
                    let t=if co.reversed{ed1-(ed1-ed0)*f}else{ed0+(ed1-ed0)*f};
                    let on_curve=carrier.evaluate(t).unwrap();
                    assert!(distance(on_surface,on_curve)<1e-7);
                }
            }
        }
        let faces:Vec<_>=corner.faces.iter().map(|f|model.faces[f.0 as usize].sense).collect();
        assert_eq!(faces,vec![FaceSense::Forward,FaceSense::Forward,
                              FaceSense::Forward,FaceSense::Forward,
                              FaceSense::Reversed]);
    }
    assert_eq!(arcs.len(),16);
    assert_eq!(arcs.len(),plan.edge_splits.registry.arcs.len());
    assert_eq!(events.len(),plan.edge_splits.registry.physical_event_count);
    assert!(events.values().all(|set|set.len()==1));
    assert_eq!(plan.source_faces.len(),20);
    assert_eq!(plan.source_faces.iter().filter(|p|
        p.source_operand==GraphOperand::Second && p.flipped_for_difference).count(),4);
}

#[test]
fn source_specific_analytic_mass_is_consistent_with_independent_quadrature(){
    let (a,b)=(square(),cylinder());
    let plan=ready(&a,&b);
    let brep=plan.brep();
    let n=16384;
    let r=2.0_f64;
    let a_side=1.5_f64;
    let u0=(r*r-a_side*a_side).sqrt();
    let h=(a_side-u0)/n as f64;
    let mut area=0.;
    let mut moment_x=0.;
    let mut moment_z=0.;
    for i in 0..=n {
        let u=u0+i as f64*h;
        let lo=(r*r-u*u).sqrt();
        let strip=a_side-lo;
        let weight=if i==0 || i==n {1.}else if i%2==0 {2.}else{4.};
        area+=weight*strip;
        moment_x+=weight*u*strip;
        moment_z+=weight*0.5*(a_side*a_side-lo*lo);
    }
    area*=h/3.;
    moment_x*=h/3.;
    moment_z*=h/3.;
    assert_eq!(plan.corners.len(),4);
    for patch in &plan.corners{
        assert!((patch.exact_area_xy-area).abs()<5e-10);
        assert!((patch.exact_volume-area*3.).abs()<2e-9);
    }
    assert!((brep.mass.volume-4.*area*3.).abs()<5e-9);
    assert!(brep.mass.surface_area>2.*brep.mass.volume/3.);
    assert!(brep.mass.centroid.x.abs()<1e-10);
    assert!((brep.mass.centroid.y-2.5).abs()<1e-10);
    assert!(brep.mass.centroid.z.abs()<1e-10);
    assert!(moment_x>0. && moment_z>0.);
}

#[test]
fn translated_asymmetric_box_yields_four_source_owned_components(){
    let a=BrepBody::from(Solid::block(
        p(2.,6.,-8.),4.5,2.5,3.5).unwrap());
    // X in [2,6.5], Z in [-8,-4.5]; center at (4,-6.2).
    // Choose a radius > all four wall distances and < every corner radius:
    // dx=[2,2.5], dz=[1.8,1.7] -> impossible: radius >2.5
    // while closest corner has distance sqrt(2^2+1.7^2)=2.62.
    let b=cutter(p(4.,1.,-6.2),2.55,10.);
    let plan=ready(&a,&b);
    assert_eq!(plan.brep().shells.len(),4);
    plan.brep().validate().unwrap();
    plan.validate_against(&a,&b,tol()).unwrap();
    assert!(plan.brep().mass.centroid.x>2.);
    assert!(plan.brep().mass.centroid.z< -4.5);
}

#[test]
fn disconnected_shell_membership_and_material_sense_fail_closed_when_tampered(){
    let (a,b)=(square(),cylinder());
    let original=ready(&a,&b);
    let mut modified=original.brep().clone();
    modified.shells[1].faces[0]=modified.shells[0].faces[0];
    assert!(modified.validate().is_err());
    let mut modified=original.brep().clone();
    let face=original.corners[0].faces[4];
    modified.faces[face.0 as usize].sense=FaceSense::Forward;
    assert!(modified.validate().is_err());
    let mut modified=original.brep().clone();
    let e=original.corners[0].edges[2];
    let c=modified.edges[e.0 as usize].coedges[0];
    modified.coedges[c.0 as usize].twin=rustsolid::CoedgeId(u32::MAX);
    assert!(modified.validate().is_err());
    let mut modified=original.clone();
    modified.corners[0].cut_arcs[0].pair=rustsolid::FacePairId(u32::MAX);
    assert!(modified.validate_against(&a,&b,tol()).is_err());
}

#[test]
fn near_tangent_corner_moved_cylinder_and_unsupported_operation_hold(){
    let (a,b)=(square(),cylinder());
    assert!(matches!(AnalyticFacePatchPartitionPlan::build(
        &a,&b,BooleanOperation::Union,tol()).unwrap(),
        AnalyticFacePatchPartitionResult::Held(
            AnalyticFacePatchHold::UnsupportedOperation)));
    let plan=ready(&a,&b);
    let moved=cutter(p(0.1,0.,0.),2.,5.);
    assert!(plan.validate_against(&a,&moved,tol()).is_err());
    let almost_touch=cutter(p(0.,0.,0.),(1.5_f64*2.0_f64.sqrt())-1e-10,5.);
    let attempt=result(&a,&almost_touch);
    assert!(matches!(attempt,AnalyticFacePatchPartitionResult::Held(_)));
    let no_intersection=cutter(p(10.,0.,0.),1.0,5.0);
    assert!(matches!(result(&a,&no_intersection),
        AnalyticFacePatchPartitionResult::Held(_)));
}


#[test]
fn every_material_patch_has_two_strict_source_membership_side_witnesses(){
    let (a,b)=(square(),cylinder());
    let plan=ready(&a,&b);
    let source=a.shared().unwrap();
    let bounds=source.bbox;
    let r=2.0_f64;
    let center=p(0.,0.,0.);
    let contains=|point:Point3|->bool {
        point.x>bounds.min.x&&point.x<bounds.max.x
            &&point.y>bounds.min.y&&point.y<bounds.max.y
            &&point.z>bounds.min.z&&point.z<bounds.max.z
            &&(point.x-center.x).hypot(point.z-center.z)>r
    };
    let mut unique=BTreeSet::new();
    for witness in &plan.material_witnesses {
        assert!(unique.insert(witness.face));
        assert!(contains(witness.retained_point));
        assert!(!contains(witness.removed_point));
        let signed=|point:Point3| {
            (point.x-witness.boundary_sample.x)*witness.outward.x
            +(point.y-witness.boundary_sample.y)*witness.outward.y
            +(point.z-witness.boundary_sample.z)*witness.outward.z
        };
        assert!((signed(witness.retained_point)+signed(witness.removed_point)).abs()<1e-10);
    }
    assert_eq!(unique.len(),plan.brep().faces.len());
    let mut tampered=plan.clone();
    tampered.material_witnesses[0].retained_point.x+=0.4;
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
}
