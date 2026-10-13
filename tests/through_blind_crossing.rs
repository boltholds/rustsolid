//! Certified analytical B-rep for one crossing through/blind cylinder pair.
//! Mathematical, topology and transaction oracles are original first-party
//! specifications. No decompiler/source translation or vendor reference used.

use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters, BrepOrigin, Curve2, Curve3, FaceId,
    FaceSense, HoleParameters, MixedFootprintContact, Point3, Solid,
    ThroughBlindCrossingBrep,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn source()->Solid {Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn through(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn blind(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn sample()->ThroughBlindCrossingBrep{
    ThroughBlindCrossingBrep::from_box(&source(),
        through(0.,0.,0.5),blind(0.5,0.,0.5,1.5)).unwrap()
}
fn close(actual:f64,wanted:f64,eps:f64){
    assert!((actual-wanted).abs()<=eps,
        "actual={actual}, wanted={wanted}, residual={}",(actual-wanted).abs());
}
fn check(body:&ThroughBlindCrossingBrep){
    body.validate().unwrap();
    assert_eq!(body.overlap().contact,MixedFootprintContact::TwoTransverseCrossings);
    assert!(!body.overlap().can_apply());
    let model=body.model();
    model.validate().unwrap();
    assert_eq!(model.origin,BrepOrigin::AnalyticThroughBlindCrossing);
    let s=model.summary();
    assert_eq!(s.vertices,13);
    assert_eq!(s.edges,21);
    assert_eq!(s.coedges,42);
    assert_eq!(s.faces,10);
    assert_eq!(s.loops,12);
    assert_eq!(s.shells,1);
    assert_eq!(s.euler_characteristic,0);
    assert_eq!(s.genus,1);
    assert_eq!(model.faces[0].loops.len(),2);
    assert_eq!(model.faces[1].loops.len(),2);
    assert_eq!(model.faces[8].loops.len(),1);
    assert_eq!(model.faces[9].loops.len(),1);
    assert_eq!(model.loops[7].first_coedge.0,26);
    assert_eq!(model.loops[10].face,FaceId(0));
    assert_eq!(model.loops[11].face,FaceId(9));
    for id in 0..6 {
        assert_eq!(model.face_sense(FaceId(id)),Some(FaceSense::Forward));
    }
    for id in [6,7,9] {
        assert_eq!(model.face_sense(FaceId(id)),Some(FaceSense::Reversed));
    }
    assert_eq!(model.face_sense(FaceId(8)),Some(FaceSense::Forward));
    assert_eq!(model.edges[18].start.0,9);
    assert_eq!(model.edges[18].end.0,8);
    assert_eq!(model.edges[19].start,model.edges[19].end);
    assert!(matches!(model.edge_curve(rustsolid::EdgeId(18)),Some(Curve3::Circle(_))));
    assert!(matches!(model.coedge_pcurve(rustsolid::CoedgeId(26)),
        Some(Curve2::ClockwiseCircle(_))));
    assert!(matches!(model.coedge_pcurve(rustsolid::CoedgeId(36)),
        Some(Curve2::Circle(_))));
    assert_eq!(model.coedges[24].face,FaceId(9));
    assert_eq!(model.coedges[26].face,FaceId(8));
    assert_eq!(model.edges[20].coedges,
        [rustsolid::CoedgeId(40),rustsolid::CoedgeId(41)]);
}

#[test]
fn true_crossed_through_blind_builds_one_genus_one_shell_with_split_height_wall(){
    let body=sample();
    check(&body);
    let m=body.model();
    let through=body.parameters().0;
    let blind=body.parameters().1;
    let ai=body.overlap().overlap_area.unwrap();
    let at=PI*through.radius.powi(2);
    let ab=PI*blind.radius.powi(2);
    let predicted=144.0-at*4.0-(ab-ai)*blind.depth;
    close(m.mass.volume,predicted,1e-10);
    close(m.mass.volume,140.14094873829234,1e-10);
    assert!(m.mass.surface_area.is_finite()&&m.mass.surface_area>0.0);
    let upper=m.oriented_face_normal(FaceId(6),2.0,0.75).unwrap();
    let lower=m.oriented_face_normal(FaceId(9),2.0,1.0).unwrap();
    let floor=m.oriented_face_normal(FaceId(8),0.0,0.0).unwrap();
    assert!(upper.y.abs()<1e-12);
    assert!(lower.y.abs()<1e-12);
    assert!(floor.y>0.0);
}

#[test]
fn translated_and_rotated_pair_reconstruct_valid_analytic_closed_shells(){
    for (dx,dz) in [(0.0,0.5),(-0.5,0.0),(0.5,0.0)] {
        let model=ThroughBlindCrossingBrep::from_box(&source(),
            through(0.,0.,0.5),blind(dx,dz,0.5,1.5)).unwrap();
        check(&model);
        close(model.model().mass.volume,140.14094873829234,1e-9);
    }
}

#[test]
fn uneven_radii_have_correct_volume_and_crescent_area(){
    let body=ThroughBlindCrossingBrep::from_box(&source(),
        through(0.,0.,0.7),blind(0.6,0.,0.5,2.2)).unwrap();
    check(&body);
    let t=PI*0.49;
    let b=PI*0.25;
    let ai=body.overlap().overlap_area.unwrap();
    close(body.model().mass.volume,144.0-t*4.0-(b-ai)*2.2,1e-10);
}

#[test]
fn crossing_source_edits_are_single_commit_and_replayable(){
    let mut body=sample();let first=body.model().mass;
    {
        let mut tx=body.begin_edit();
        tx.replace_through(through(-0.1,0.05,0.6)).unwrap();
        tx.replace_blind(blind(0.55,0.05,0.55,2.0)).unwrap();
        assert_eq!(tx.commit().unwrap(),1);
    }
    check(&body);
    assert_ne!(body.model().mass.volume,first.volume);
    assert_eq!(body.undo_depth(),1);
    assert_eq!(body.redo_depth(),0);
    assert!(body.undo().unwrap());
    assert_eq!(body.revision(),2);
    check(&body);
    assert_eq!(body.model().mass,first);
    assert!(body.redo().unwrap());
    assert_eq!(body.revision(),3);
    check(&body);
}

#[test]
fn failed_staging_and_failed_late_commit_do_not_modify_history_or_brep(){
    let mut body=sample();
    let original=body.model().mass;
    {
        let mut tx=body.begin_edit();
        tx.replace_through(through(-0.1,0.0,0.5)).unwrap();
        assert!(tx.replace_blind(blind(0.5,0.,-1.0,1.5)).is_err());
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=body.begin_edit();
        tx.replace_through(through(2.9,0.0,0.5)).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(body.model().mass,original);
    assert_eq!(body.revision(),0);
    assert_eq!(body.undo_depth(),0);
    check(&body);
}

#[test]
fn tangent_nested_disjoint_and_invalid_inputs_do_not_create_intersection_topology(){
    let block=source();
    for h in [
        blind(1.,0.,0.5,1.5),
        blind(2.5,0.,0.2,1.5),
        blind(0.,0.,0.3,1.5),
        blind(0.,0.,0.5,1.5),
        blind(0.5,0.,0.5,0.0),
        blind(0.5,0.,0.5,4.0),
    ]{
        assert!(ThroughBlindCrossingBrep::from_box(
            &block,through(0.,0.,0.5),h).is_err());
    }
    assert_eq!(block.mass.volume,144.0);
}
