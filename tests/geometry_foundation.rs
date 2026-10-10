use rustsolid::{
    AngleUnit, Circle2, Circle3, CoedgeId, Curve2, Curve3, CylinderSurface, EdgeId,
    FaceId, Frame3, GeometryStore, GeometryTolerance, LengthUnit, Line2, Line3,
    ModelUnits, ParameterRange, PlaneSurface, Point2Param, Point3, Solid,
    SphereSurface, Surface3,
};
use std::f64::consts::{FRAC_PI_2, TAU};

fn pt(x: f64, y: f64, z: f64) -> Point3 { Point3 {x,y,z} }
fn assert_close(a:f64,b:f64,eps:f64) { assert!((a-b).abs()<=eps,"{a} != {b}"); }
fn assert_pt(actual:Point3,expected:Point3,eps:f64){
    assert_close(actual.x,expected.x,eps);
    assert_close(actual.y,expected.y,eps);
    assert_close(actual.z,expected.z,eps);
}
fn frame() -> Frame3 {
    Frame3::from_axes(pt(1.0,2.0,3.0), pt(2.0,0.0,0.0), pt(1.0,3.0,0.0), GeometryTolerance::default()).unwrap()
}
fn cube() -> Solid {
    Solid::block(pt(5.0,-7.0,10.0), 2.0, 3.0, 4.0).unwrap()
}

#[test]
fn frame_orthonormalizes_axes_and_projects_point() {
    let frame=frame();
    assert_pt(frame.x,pt(1.0,0.0,0.0),1e-12);
    assert_pt(frame.y,pt(0.0,1.0,0.0),1e-12);
    assert_pt(frame.normal,pt(0.0,0.0,1.0),1e-12);
    let p=frame.position(4.0,-2.0).unwrap();
    assert_pt(p,pt(5.0,0.0,3.0),1e-12);
    assert_eq!(frame.project(p).unwrap(),Point2Param{u:4.0,v:-2.0});
    assert!(Frame3::from_axes(pt(0.0,0.0,0.0),pt(1.0,0.0,0.0),pt(2.0,0.0,0.0),GeometryTolerance::default()).is_err());
}
#[test]
fn analytic_curves_keep_parameter_contracts() {
    let tol=GeometryTolerance::default();
    let (line,interval)=Line3::between(pt(1.0,2.0,3.0),pt(1.0,5.0,3.0),tol).unwrap();
    assert_eq!(interval,ParameterRange::Bounded {start:0.0,end:3.0});
    assert_pt(Curve3::Line(line).evaluate(3.0).unwrap(),pt(1.0,5.0,3.0),1e-12);
    assert_pt(Curve3::Line(line).derivative(1.0).unwrap(),pt(0.0,1.0,0.0),1e-12);
    let circle=Circle3::new(frame(),2.0,tol).unwrap();
    assert_pt(Curve3::Circle(circle).evaluate(FRAC_PI_2).unwrap(),pt(1.0,4.0,3.0),1e-12);
    assert_pt(Curve3::Circle(circle).derivative(0.0).unwrap(),pt(0.0,2.0,0.0),1e-12);
    assert!(Circle3::new(frame(),0.0,tol).is_err());
    assert!(ParameterRange::bounded(1.0,1.0).is_err());
    assert!(!(ParameterRange::Bounded{start:0.0,end:2.0}).includes(4.0,0.0));
    let line2=Curve2::Line(Line2{origin:Point2Param{u:1.0,v:1.0},delta:Point2Param{u:2.0,v:-3.0}});
    assert_eq!(line2.evaluate(0.5).unwrap(),Point2Param{u:2.0,v:-0.5});
    let circle2=Curve2::Circle(Circle2::new(Point2Param{u:2.0,v:3.0},1.0,tol).unwrap());
    let uv=circle2.evaluate(FRAC_PI_2).unwrap();
    assert_close(uv.u,2.0,1e-12);
    assert_close(uv.v,4.0,1e-12);
}
#[test]
fn analytic_surface_evaluation_and_derivatives() {
    let t=GeometryTolerance::default();
    let p=Surface3::Plane(PlaneSurface{frame:frame()});
    let e=p.evaluate(2.0,5.0).unwrap();
    assert_pt(e.point,pt(3.0,7.0,3.0),1e-12);
    assert_pt(e.du,pt(1.0,0.0,0.0),1e-12);
    assert_pt(e.dv,pt(0.0,1.0,0.0),1e-12);
    assert_eq!(p.domain().u,ParameterRange::Unbounded);
    let cyl=Surface3::Cylinder(CylinderSurface::new(frame(),2.0,t).unwrap());
    assert_pt(cyl.evaluate(0.0,4.0).unwrap().point,pt(3.0,2.0,7.0),1e-12);
    assert_pt(cyl.evaluate(0.0,4.0).unwrap().du,pt(0.0,2.0,0.0),1e-12);
    assert_pt(cyl.evaluate(0.0,4.0).unwrap().dv,pt(0.0,0.0,1.0),1e-12);
    assert!(cyl.domain().u_periodic);
    assert_eq!(cyl.domain().u,ParameterRange::Bounded{start:0.0,end:TAU});
    let sphere=Surface3::Sphere(SphereSurface::new(frame(),2.0,t).unwrap());
    assert_pt(sphere.evaluate(0.0,FRAC_PI_2).unwrap().point,pt(1.0,2.0,5.0),1e-12);
    assert_pt(sphere.evaluate(FRAC_PI_2,0.0).unwrap().point,pt(1.0,4.0,3.0),1e-12);
    assert_eq!(sphere.domain().v,ParameterRange::Bounded{start:-FRAC_PI_2,end:FRAC_PI_2});
    assert!(sphere.evaluate(f64::INFINITY,0.0).is_err());
}
#[test]
fn derived_store_maps_every_brep_entity_to_analytic_geometry() {
    let body=cube();
    let store=GeometryStore::from_solid(&body,ModelUnits::default()).unwrap();
    assert_eq!((store.curves3.len(),store.surfaces.len(),store.curves2.len()),
        (body.edges.len(),body.faces.len(),body.coedges.len()));
    assert_eq!(store.units.length,LengthUnit::Millimeter);
    assert_eq!(store.units.angle,AngleUnit::Radian);
    store.validate_bindings(&body).unwrap();
    for e in &body.edges {
        assert!(matches!(store.edge_curve(e.id),Some(Curve3::Line(_))));
        let edge=store.edges.get(&e.id).unwrap();
        let ParameterRange::Bounded {start,end}=edge.domain else {panic!("edge needs bounded parameter")};
        let curve=store.edge_curve(e.id).unwrap();
        assert_pt(curve.evaluate(start).unwrap(),body.vertices[e.start.0 as usize].position,1e-10);
        assert_pt(curve.evaluate(end).unwrap(),body.vertices[e.end.0 as usize].position,1e-10);
    }
    for face in &body.faces {
        assert!(matches!(store.face_surface(face.id),Some(Surface3::Plane(_))));
        let binding=store.faces.get(&face.id).unwrap();
        assert!(binding.bounds.u.includes(0.0,5.0));
    }
    for coedge in &body.coedges {
        let edge=&body.edges[coedge.edge.0 as usize];
        let carrier=store.face_surface(coedge.face).unwrap();
        let pcurve=store.trim_curve(coedge.id).unwrap();
        let start=body.vertices[coedge.start_vertex(edge).0 as usize].position;
        let end=body.vertices[coedge.end_vertex(edge).0 as usize].position;
        for (t,vertex) in [(0.0,start),(1.0,end)] {
            let uv=pcurve.evaluate(t).unwrap();
            assert_pt(carrier.evaluate(uv.u,uv.v).unwrap().point,vertex,1e-9);
        }
    }
    assert!(store.trim_curve(CoedgeId(999)).is_none());
    assert!(store.face_surface(FaceId(999)).is_none());
}
#[test]
fn scoped_tolerances_are_checked_at_brep_bindings() {
    let body=cube();
    let mut store=GeometryStore::from_solid(&body,ModelUnits::default()).unwrap();
    let invalid=GeometryTolerance{absolute_length:10.0,..Default::default()};
    store.set_edge_tolerance(EdgeId(0),invalid).unwrap();
    assert!(store.validate_bindings(&body).is_err());
    store.set_edge_tolerance(EdgeId(0),body.tolerance).unwrap();
    store.validate_bindings(&body).unwrap();
    let vertex=body.edges[0].start;
    store.set_vertex_tolerance(vertex,invalid).unwrap();
    assert!(store.validate_bindings(&body).is_err());
    assert!(store.set_edge_tolerance(EdgeId(999),body.tolerance).is_err());
    assert!(store.set_vertex_tolerance(rustsolid::VertexId(999),body.tolerance).is_err());
    assert_close(LengthUnit::Millimeter.convert(25.4,LengthUnit::Inch).unwrap(),1.0,1e-12);
    assert_close(LengthUnit::Inch.convert(1.0,LengthUnit::Meter).unwrap(),0.0254,1e-12);
}
#[test]
fn stores_are_bound_to_body_revision_and_incarnation() {
    let mut body=cube();
    let store=GeometryStore::from_solid(&body,ModelUnits::default()).unwrap();
    let other=cube();
    assert!(store.check_current(&other).is_err());
    body.edit_atomic(|tx| {tx.split_edge(EdgeId(0),0.5)?;Ok(())}).unwrap();
    assert!(store.check_current(&body).is_err());
    let current=GeometryStore::from_solid(&body,ModelUnits::default()).unwrap();
    current.validate_bindings(&body).unwrap();
}
#[test]
fn binding_survives_small_part_far_from_origin() {
    let body=Solid::block(pt(1e12,1e12,-1e12),1.0,1.0,1.0).unwrap();
    let store=GeometryStore::from_solid(&body,ModelUnits::default()).unwrap();
    store.validate_bindings(&body).unwrap();
}
