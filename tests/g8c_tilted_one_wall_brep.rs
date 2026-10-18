use std::collections::{BTreeMap,BTreeSet};
use rustsolid::{
    AnalyticOneWallPlan,AnalyticOneWallResult,BooleanOperation,BrepBody,
    BrepOrigin,Curve2,Curve3,CylindricalBrep,FaceSense,Frame3,
    GeometryTolerance,GlobalEventId,GraphOperand,
    IndependentTiltedBoundaryFluxOracle,IndependentTiltedBoundaryFluxResult,
    IndependentTiltedVolumeOracle,IndependentTiltedVolumeResult,
    Point3,Solid,Surface3,TiltedOneWallShellPlan,TiltedOneWallShellResult,
    TiltedOneWallSourcePlan,TiltedOneWallSourceResult,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn dist(a:Point3,b:Point3)->f64{
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}
fn box_body()->BrepBody{
    BrepBody::from(Solid::block(p(-1.,1.,-3.),4.,3.,6.).unwrap())
}
fn tilted_x(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(theta.cos(),-theta.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn tilted_z(theta:f64)->BrepBody{
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(1.,0.,0.),
        p(0.,-theta.sin(),-theta.cos()),tol()).unwrap();
    BrepBody::from(CylindricalBrep::new(frame,1.5,10.,tol()).unwrap())
}
fn make(a:&BrepBody,b:&BrepBody)
    ->(TiltedOneWallShellPlan,rustsolid::IndependentTiltedBoundaryFluxOracle)
{
    let TiltedOneWallSourceResult::Ready(s)=
        TiltedOneWallSourcePlan::build(a,b,BooleanOperation::Difference,tol()).unwrap()
        else{panic!("not G8a source")};
    let IndependentTiltedVolumeResult::Ready(volume)=
        IndependentTiltedVolumeOracle::build(a,b,&s,tol()).unwrap()
        else{panic!("not G8b volume")};
    let shape=TiltedOneWallShellPlan::build(a,b,&s,&volume,tol()).unwrap();
    let TiltedOneWallShellResult::Ready(shell)=shape else{
        panic!("G8c shell Held: {shape:?}");
    };
    let flux=IndependentTiltedBoundaryFluxOracle::build(a,b,&shell,tol()).unwrap();
    let IndependentTiltedBoundaryFluxResult::Ready(cert)=flux else{
        panic!("G8c actual coedge flux Held: {flux:?}");
    };
    (shell,cert)
}
fn assert_candidate(
    a:&BrepBody,b:&BrepBody,
    shell:&TiltedOneWallShellPlan,
    flux:&rustsolid::IndependentTiltedBoundaryFluxOracle,
){
    let m=shell.brep();
    m.validate().unwrap();
    assert_eq!(m.origin,BrepOrigin::AnalyticTiltedOneWallCylinderDifference);
    assert_eq!(m.shells.len(),1);
    assert_eq!(m.vertices.len(),12);
    assert_eq!(m.edges.len(),18);
    assert_eq!(m.coedges.len(),36);
    assert_eq!(m.faces.len(),8);
    assert_eq!(m.loops.len(),8);
    assert_eq!(m.euler_characteristic(),2);
    assert_eq!(m.summary().genus,0);
    assert_eq!(m.summary().cylindrical_faces,1);
    assert!(shell.closed_shell_proven&&!shell.boundary_flux_certified);
    assert!(!shell.can_apply);
    assert_eq!(shell.radial.shell_components,1);
    assert_eq!(shell.radial.euler_characteristic,2);
    assert_eq!(shell.radial.total_genus,0);
    assert_eq!(flux.witness.harmonic_integrated_uses+
        flux.witness.vertical_generator_uses,4);
    assert!(flux.volume_certified&&!flux.can_apply);
    assert_eq!(m.face_sense(shell.faces[7]),Some(FaceSense::Reversed));
    assert_eq!(shell.source_faces.len(),8);
    assert_eq!(shell.source_faces.iter().filter(|f|
        f.operand==GraphOperand::First &&
        f.original_face==shell.source_faces[5].original_face).count(),2);
    assert!(flux.witness.inward_cylinder_flux<0.);
    assert!(flux.witness.original_x_face_flux>0.);
    assert!(flux.witness.result_volume>0.);
    assert!(flux.witness.result_volume<
        flux.witness.original_x_face_flux);
    assert!((flux.witness.result_volume
        -(flux.witness.original_x_face_flux+
          flux.witness.inward_cylinder_flux)).abs()<1e-10);
    assert!(flux.witness.result_volume>=shell.source_volume.retained.lower-1e-7);
    assert!(flux.witness.result_volume<=shell.source_volume.retained.upper+1e-7);
    assert!((flux.witness.result_volume-m.mass.volume).abs()<2e-7);
    shell.validate_against(a,b,tol()).unwrap();
    flux.validate_against(a,b,shell,tol()).unwrap();
}
#[test]
fn g8c_x_tilt_constructs_closed_analytic_shell_and_certifies_actual_coedge_flux(){
    for theta in [0.04,-0.03,0.055]{
        let (a,b)=(box_body(),tilted_x(theta));
        let old_a=format!("{a:?}");let old_b=format!("{b:?}");
        let (shell,flux)=make(&a,&b);
        assert_candidate(&a,&b,&shell,&flux);
        assert_eq!(flux.witness.harmonic_integrated_uses,4);
        assert_eq!(flux.witness.vertical_generator_uses,0);
        assert_eq!(format!("{a:?}"),old_a);
        assert_eq!(format!("{b:?}"),old_b);
    }
}
#[test]
fn g8c_yz_tilt_keeps_true_source_generators_and_elliptic_cap_carriers(){
    let a=box_body();let b=tilted_z(0.04);
    let (shell,flux)=make(&a,&b);
    assert_candidate(&a,&b,&shell,&flux);
    assert_eq!(flux.witness.vertical_generator_uses,2);
    assert_eq!(flux.witness.harmonic_integrated_uses,2);
    let m=shell.brep();
    for i in [4usize,10]{
        assert!(matches!(m.edge_curve(shell.edges[i]),Some(Curve3::Harmonic(_))));
    }
    for i in [16usize,17]{
        assert!(matches!(m.edge_curve(shell.edges[i]),Some(Curve3::Line(_))));
    }
}
#[test]
fn g8c_every_cut_edge_and_dual_uv_trim_matches_unique_original_g1_arc(){
    let (a,b)=(box_body(),tilted_x(0.04));
    let (shell,_)=make(&a,&b);
    let m=shell.brep();
    let registry=&shell.source_volume; // source proof stays separate
    assert!(registry.source_bracket_established);
    let source=match TiltedOneWallSourcePlan::build(&a,&b,
        BooleanOperation::Difference,tol()).unwrap(){
        TiltedOneWallSourceResult::Ready(x)=>x,
        _=>panic!("invalid sources"),
    };
    let mut ids=BTreeSet::new();
    let mut events=BTreeMap::<GlobalEventId,rustsolid::VertexId>::new();
    for cut in &shell.physical_cuts{
        assert!(ids.insert(cut.source_arc));
        let original=source.source_cuts.arcs.iter()
            .find(|arc|arc.source==cut.source_arc).unwrap();
        let e=&m.edges[cut.result_edge.0 as usize];
        assert_eq!(cut.first_face,original.first_face);
        assert_eq!(cut.second_face,original.second_face);
        assert_eq!(m.edge_curve(cut.result_edge),Some(original.carrier));
        assert_eq!(m.geometry.edges[cut.result_edge.0 as usize].domain,
            original.parameter_domain);
        for (event,vertex) in cut.source_endpoints.iter()
            .zip(cut.result_vertices.iter()){
            if let Some(previous)=events.insert(*event,*vertex) {
                assert_eq!(previous,*vertex);
            }
        }
        let [aa,bb]=e.coedges;
        let lhs=&m.coedges[aa.0 as usize];
        let rhs=&m.coedges[bb.0 as usize];
        assert_ne!(lhs.reversed,rhs.reversed);
        assert_eq!(lhs.twin,bb);
        assert_eq!(rhs.twin,aa);
        for &co_id in &[aa,bb]{
            let co=&m.coedges[co_id.0 as usize];
            let pcurve=m.coedge_pcurve(co_id).unwrap();
            let support=m.face_surface(co.face).unwrap();
            let edge=m.edge_curve(co.edge).unwrap();
            let rustsolid::ParameterRange::Bounded{start,end}=
                m.geometry.edges[co.edge.0 as usize].domain
                else{panic!("unbounded source edge")};
            let rustsolid::ParameterRange::Bounded{start:lo,end:hi}=
                m.geometry.coedges[co.id.0 as usize].domain
                else{panic!("unbounded pcurve")};
            for f in [0.,0.08,0.3,0.5,0.7,0.92,1.]{
                let uv=pcurve.evaluate(lo+(hi-lo)*f).unwrap();
                let q=support.evaluate(uv.u,uv.v).unwrap().point;
                let t=if co.reversed{1.-f}else{f};
                let direct=edge.evaluate(start+(end-start)*t).unwrap();
                assert!(dist(q,direct)<5e-8);
            }
        }
    }
    assert_eq!(events.len(),4);
    assert_eq!(ids.len(),4);
    for i in 0..m.faces.len(){
        let count=m.coedges.iter().filter(|c|c.face.0 as usize==i).count();
        assert_eq!(count,if i<=1{6}else{4});
    }
}
#[test]
fn g8c_periodic_cylinder_uv_loop_is_physically_closed_without_virtual_vertices(){
    let (a,b)=(box_body(),tilted_x(-0.03));
    let (shell,_)=make(&a,&b);
    let m=shell.brep();
    let face=shell.faces[7];
    let mut at=m.loops[m.faces[face.0 as usize].loops[0].0 as usize].first_coedge;
    let first=at;
    let mut points=Vec::new();
    loop {
        let g=m.geometry.coedges[at.0 as usize];
        let rustsolid::ParameterRange::Bounded{start,end}=g.domain
            else{panic!("G8c must be bounded")};
        let c=m.coedge_pcurve(at).unwrap();
        points.push((c.evaluate(start).unwrap(),c.evaluate(end).unwrap()));
        at=m.coedges[at.0 as usize].next;
        if at==first{break;}
    }
    assert_eq!(points.len(),4);
    for i in 0..4{
        let a=points[i].1;let b=points[(i+1)%4].0;
        assert!((a.u-b.u).abs()<1e-8);
        assert!((a.v-b.v).abs()<1e-8);
    }
    assert!(points.iter().any(|(a,b)|(a.u-b.u).abs()>3.));
}
#[test]
fn g8c_translated_offcenter_sources_keep_volume_and_material_orientation(){
    let a=BrepBody::from(Solid::block(p(4.,2.,-6.),5.,4.,8.).unwrap());
    let angle=0.03_f64;
    let frame=Frame3::from_axes(
        p(4.8,-3.,-2.5),
        p(angle.cos(),-angle.sin(),0.),p(0.,0.,-1.),tol()).unwrap();
    let b=BrepBody::from(CylindricalBrep::new(frame,1.8,14.,tol()).unwrap());
    let (shell,flux)=make(&a,&b);
    assert_candidate(&a,&b,&shell,&flux);
    assert!((flux.witness.original_x_face_flux-160.).abs()<1e-9);
}
#[test]
fn g8c_invalid_face_sense_twins_pcurves_source_volume_and_revisions_fail_closed(){
    let a=box_body();let b=tilted_x(0.04);
    let (shell,flux)=make(&a,&b);
    let m=shell.brep();
    let mut broken=m.clone();
    broken.faces[7].sense=FaceSense::Forward;
    assert!(broken.validate().is_err());
    let mut broken=m.clone();
    let edge=broken.edges[shell.edges[4].0 as usize].coedges;
    broken.coedges[edge[0].0 as usize].twin=edge[0];
    assert!(broken.validate().is_err());
    let mut broken=m.clone();
    let cut=broken.edges[shell.edges[4].0 as usize].coedges;
    let co=broken.geometry.coedges[cut[0].0 as usize];
    if let Curve2::Harmonic(ref mut c)=
        broken.geometry.curves2[co.pcurve.0 as usize] {
        c.center.u+=0.04;
    }else{panic!("expected harmonic plane trim")};
    assert!(broken.validate().is_err());
    let mut tampered=shell.clone();
    tampered.source_volume.retained.lower+=0.1;
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut tampered=shell.clone();
    tampered.physical_cuts[0].source_endpoints[0]=GlobalEventId(u32::MAX);
    assert!(tampered.validate_against(&a,&b,tol()).is_err());
    let mut changed=a.clone();
    if let BrepBody::Polyhedral(s)=&mut changed{s.revision+=1;}
    assert!(flux.validate_against(&changed,&b,&shell,tol()).is_err());
    let mut false_flux=flux.clone();
    false_flux.witness.inward_cylinder_flux+=0.5;
    assert!(false_flux.validate_against(&a,&b,&shell,tol()).is_err());
}
#[test]
fn g8c_cannot_promote_unsupported_upright_or_multiwall_boolean(){
    let a=box_body();
    let upright=BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),1.5,5.,tol()).unwrap());
    assert!(matches!(AnalyticOneWallPlan::build(&a,&upright,
        BooleanOperation::Difference,tol()).unwrap(),
        AnalyticOneWallResult::Ready(_)));
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&upright,
        BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(_)));
    let frame=Frame3::from_axes(
        p(0.,-2.,0.),p(0.04_f64.cos(),-0.04_f64.sin(),0.),
        p(0.,0.,-1.),tol()).unwrap();
    let wide=BrepBody::from(CylindricalBrep::new(frame,4.,10.,tol()).unwrap());
    assert!(matches!(TiltedOneWallSourcePlan::build(&a,&wide,
        BooleanOperation::Difference,tol()).unwrap(),
        TiltedOneWallSourceResult::Held(_)));
}
