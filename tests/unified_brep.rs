use rustsolid::{
    BrepBody, BrepModel, BrepOrigin, CoedgeId, Curve2, Curve3, CylindricalBrep,
    EdgeId, FaceId, Frame3, GeometryTolerance, Point3, Solid, Surface3,
};
use std::f64::consts::{PI, TAU};

fn cube() -> Solid {
    Solid::block(Point3{x:0.0,y:0.0,z:0.0},3.0,4.0,5.0).unwrap()
}
fn cylinder() -> CylindricalBrep {
    CylindricalBrep::upright(Point3{x:0.0,y:0.0,z:0.0},2.0,5.0,
        GeometryTolerance::default()).unwrap()
}

#[test]
fn both_sources_share_identical_brep_queries_without_faceting_curves() {
    let body = BrepBody::from(cube());
    let model = body.shared().unwrap();
    assert_eq!(model.origin,BrepOrigin::Polyhedral);
    assert_eq!(model.summary().plane_faces,6);
    assert_eq!(model.summary().cylindrical_faces,0);
    assert_eq!(model.summary().straight_edges,12);
    assert_eq!(model.summary().closed_edges,0);
    assert_eq!(model.summary().genus,0);
    assert!(matches!(model.face_surface(FaceId(0)),Some(Surface3::Plane(_))));
    assert!(matches!(model.edge_curve(EdgeId(0)),Some(Curve3::Line(_))));
    assert!(matches!(model.coedge_pcurve(CoedgeId(0)),Some(Curve2::Line(_))));
    assert_eq!(body.tessellate(8).unwrap().triangles.len(),12);

    let body = BrepBody::from(cylinder());
    let model = body.shared().unwrap();
    assert_eq!(model.origin,BrepOrigin::AnalyticCylinder);
    assert_eq!(model.summary().plane_faces,2);
    assert_eq!(model.summary().cylindrical_faces,1);
    assert_eq!(model.summary().straight_edges,1);
    assert_eq!(model.summary().circular_edges,2);
    assert_eq!(model.summary().closed_edges,2);
    assert_eq!(model.summary().seam_edges,1);
    assert_eq!(model.summary().genus,0);
    assert_eq!(model.euler_characteristic(),2);
    assert!(matches!(model.face_surface(FaceId(2)),Some(Surface3::Cylinder(_))));
    assert!(matches!(model.edge_curve(EdgeId(0)),Some(Curve3::Circle(_))));
    assert!(matches!(model.coedge_pcurve(CoedgeId(1)),Some(Curve2::Circle(_))));
    assert_eq!(model.coedges_of_edge(EdgeId(2)), Some([CoedgeId(4), CoedgeId(5)]));
    assert_eq!(model.next_of_edge(CoedgeId(4)),Some(CoedgeId(5)));
    assert_eq!(model.next_in_loop(CoedgeId(0)),Some(CoedgeId(4)));
    assert_eq!(model.loops_of_face(FaceId(0)).unwrap().len(),1);
    assert_eq!(body.tessellate(8).unwrap().triangles.len(),32);
}

#[test]
fn tessellation_resolution_is_independent_from_shared_analytic_brep() {
    let body=BrepBody::from(cylinder());
    let first=body.shared().unwrap();
    let coarse=body.tessellate(8).unwrap();
    let fine=body.tessellate(128).unwrap();
    let second=body.shared().unwrap();
    assert_eq!((coarse.triangles.len(),fine.triangles.len()),(32,512));
    assert_eq!(first.summary(),second.summary());
    assert_eq!(first.geometry.surfaces,second.geometry.surfaces);
    assert_eq!(first.geometry.curves3,second.geometry.curves3);
    assert_eq!(first.mass.volume,PI*20.0);
}

#[test]
fn shared_validator_detects_bad_twin_and_seam_uv_even_with_valid_source() {
    let mut m=BrepModel::from_cylinder(&cylinder()).unwrap();
    m.coedges[4].twin=CoedgeId(0);
    assert!(m.validate().is_err());
    let mut m=BrepModel::from_cylinder(&cylinder()).unwrap();
    m.geometry.curves2[4]=m.geometry.curves2[5];
    assert!(m.validate().is_err());
    let mut m=BrepModel::from_cylinder(&cylinder()).unwrap();
    m.faces[2].geometry.bounds.u_periodic=false;
    assert!(m.validate().is_err());
    let mut m=BrepModel::from_cylinder(&cylinder()).unwrap();
    m.geometry.curves3[0]=m.geometry.curves3[2];
    assert!(m.validate().is_err());
}

#[test]
fn the_generic_validator_handles_one_coedge_circle_loop_and_shared_face_seam() {
    let m=BrepModel::from_cylinder(&cylinder()).unwrap();
    m.validate().unwrap();
    let cap = m.loops_of_face(FaceId(0)).unwrap()[0];
    assert_eq!(m.loops[cap.0 as usize].first_coedge,CoedgeId(1));
    assert_eq!(m.next_in_loop(CoedgeId(1)),Some(CoedgeId(1)));
    let pair = m.coedges_of_edge(EdgeId(2)).unwrap();
    assert_eq!(m.coedges[pair[0].0 as usize].face,m.coedges[pair[1].0 as usize].face);
    let a=m.coedge_pcurve(pair[0]).unwrap().evaluate(0.5).unwrap();
    let b=m.coedge_pcurve(pair[1]).unwrap().evaluate(0.5).unwrap();
    assert!((a.u-b.u).abs()>TAU-1e-12);
}

#[test]
fn offset_analytic_geometry_and_polyhedral_local_edits_convert_to_common_graph() {
    let tol=GeometryTolerance::default();
    let frame=Frame3::from_axes(Point3{x:1.0e9,y:-5.0,z:8.0},
        Point3{x:0.0,y:1.0,z:0.0},Point3{x:0.0,y:0.0,z:1.0},tol).unwrap();
    let c=CylindricalBrep::new(frame,2.0,6.0,tol).unwrap();
    BrepModel::from_cylinder(&c).unwrap().validate().unwrap();

    let mut poly=cube();
    poly.edit_atomic(|tx|{tx.split_edge(EdgeId(0),0.25)?;Ok(())}).unwrap();
    let generic=BrepModel::from_polyhedral(&poly,rustsolid::ModelUnits::default()).unwrap();
    assert_eq!(generic.summary().vertices,9);
    assert_eq!(generic.summary().edges,13);
    assert_eq!(generic.summary().cylindrical_faces,0);
    generic.validate().unwrap();
}
