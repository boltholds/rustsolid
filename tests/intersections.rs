use rustsolid::{
    CylinderSurface, Frame3, GeometryError, GeometryTolerance, IntersectionAmbiguity,
    IntersectionContact, IntersectionCurve, Line3, PlaneSurface, Point3,
    Surface3, SurfaceIntersection, SurfaceIntersectionEngine,
};
use std::f64::consts::{PI,TAU};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn engine()->SurfaceIntersectionEngine{SurfaceIntersectionEngine::new(tol()).unwrap()}
fn cylinder(origin:Point3,r:f64)->Surface3{
    Surface3::Cylinder(CylinderSurface::new(Frame3::from_axes(origin,p(1.,0.,0.),p(0.,0.,-1.),tol()).unwrap(),r,tol()).unwrap())
}
fn plane(origin:Point3,x:Point3,y:Point3)->Surface3{
    Surface3::Plane(PlaneSurface{frame:Frame3::from_axes(origin,x,y,tol()).unwrap()})
}
fn plane_x(x:f64)->Surface3{plane(p(x,0.,0.),p(0.,1.,0.),p(0.,0.,1.))}
fn plane_y(y:f64)->Surface3{plane(p(0.,y,0.),p(1.,0.,0.),p(0.,0.,-1.))}
fn plane_z(z:f64)->Surface3{plane(p(0.,0.,z),p(1.,0.,0.),p(0.,1.,0.))}
fn norm(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn sub(a:Point3,b:Point3)->Point3{p(a.x-b.x,a.y-b.y,a.z-b.z)}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn close(a:f64,b:f64,eps:f64){assert!((a-b).abs()<=eps,"expected {b}, got {a}, diff={}",(a-b).abs());}
fn xyz(a:Point3,b:Point3,eps:f64){close(norm(sub(a,b)),0.,eps);}
fn certify(result:&SurfaceIntersection,a:Surface3,b:Surface3){
    let proof=result.certify(a,b,tol()).unwrap();
    assert!(proof.samples>=5);
    assert!(proof.max_residual_a<=proof.budget);
    assert!(proof.max_residual_b<=proof.budget);
}
#[test]
fn orthogonal_planes_intersect_in_exact_y_axis_and_swap_as_sets(){
    let a=plane_x(0.);let b=plane_z(0.);
    let result=engine().intersect(a,b).unwrap();
    assert_eq!(result.kind(),"crossing");
    assert_eq!(result.curves().len(),1);
    certify(&result,a,b);
    let reversed=engine().intersect(b,a).unwrap();
    certify(&reversed,b,a);
    for t in [-12.,0.,7.] {
        let q=result.curves()[0].evaluate(t).unwrap();
        close(q.x,0.,1e-10);close(q.z,0.,1e-10);
        let line=reversed.curves()[0];
        let IntersectionCurve::Line(reverse)=line else {panic!("expected line")};
        // Point-to-line distance independent of orientation/choice of origin.
        let v=sub(q,reverse.origin);
        let projected=dot(v,reverse.direction);
        xyz(v,p(reverse.direction.x*projected,reverse.direction.y*projected,reverse.direction.z*projected),1e-9);
    }
}
#[test]
fn parallel_planes_and_coincident_planes_are_distinct(){
    assert!(matches!(engine().intersect(plane_x(0.),plane_x(2.)).unwrap(),SurfaceIntersection::Empty));
    assert!(matches!(engine().intersect(plane_x(0.),plane_x(0.)).unwrap(),SurfaceIntersection::Coincident));
}
#[test]
fn cylinder_perpendicular_plane_is_exact_periodic_circle(){
    let a=cylinder(p(0.,0.,0.),2.);let b=plane_y(3.);
    let out=engine().intersect(a,b).unwrap();
    assert_eq!(out.kind(),"crossing");
    assert_eq!(out.curves().len(),1);
    assert!(matches!(out.curves()[0],IntersectionCurve::PlaneCylinderSection(_)));
    certify(&out,a,b);
    for t in [0.,0.37,PI,TAU]{
        let q=out.curves()[0].evaluate(t).unwrap();
        close(q.y,3.,1e-11);close(q.x*q.x+q.z*q.z,4.,1e-10);
    }
    xyz(out.curves()[0].evaluate(0.).unwrap(),out.curves()[0].evaluate(TAU).unwrap(),1e-10);
}
#[test]
fn oblique_plane_section_satisfies_two_independent_implicit_equations(){
    let a=cylinder(p(0.,0.,0.),2.);
    let b=plane(p(0.,1.,0.),p(0.,0.,1.),p(2.,-1.,0.));
    let out=engine().intersect(a,b).unwrap();
    assert_eq!(out.curves().len(),1);
    certify(&out,a,b);
    for t in [0.04,0.7,1.9,3.7,5.3]{
        let q=out.curves()[0].evaluate(t).unwrap();
        // plane x + 2(y-1) = 0
        close(q.x+2.*(q.y-1.),0.,1e-9);
        close(q.x*q.x+q.z*q.z,4.,1e-9);
    }
}
#[test]
fn longitudinal_plane_creates_exactly_two_generatrix_lines(){
    let a=plane_x(0.);let b=cylinder(p(0.,0.,0.),2.);
    let out=engine().intersect(a,b).unwrap();
    assert_eq!(out.curves().len(),2);
    certify(&out,a,b);
    let mut z=Vec::new();
    for branch in out.curves(){
        assert!(matches!(branch,IntersectionCurve::Line(_)));
        let q=branch.evaluate(0.5).unwrap();
        close(q.x,0.,1e-10);close(q.y,0.5,1e-10);
        close(q.z.abs(),2.,1e-10);
        z.push(q.z);
    }
    assert!(z[0]*z[1]<0.0);
}
#[test]
fn tangent_is_not_secant_and_near_tangent_is_unresolved(){
    let c=cylinder(p(0.,0.,0.),2.);
    let at=engine().intersect(plane_x(2.),c).unwrap();
    assert_eq!(at.kind(),"tangent");
    assert_eq!(at.curves().len(),1);
    certify(&at,plane_x(2.),c);
    assert!(matches!(engine().intersect(plane_x(2.000001),c).unwrap(),SurfaceIntersection::Empty));
    assert_eq!(engine().intersect(plane_x(1.999999),c).unwrap().curves().len(),2);
    for x in [2.-4e-10,2.+4e-10]{
        let out=engine().intersect(plane_x(x),c).unwrap();
        assert!(matches!(out,SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearTangency,..}),"{out:?}");
    }
}
#[test]
fn near_parallel_orientation_is_reported_not_silently_classified(){
    // Normal is (-1,+1e-10,0), nearly perpendicular to +Y axis.
    let almost_axis_parallel=plane(p(0.,0.,0.),p(0.,0.,1.),p(1e-10,1.,0.));
    let c=cylinder(p(0.,0.,0.),2.);
    let result=engine().intersect(almost_axis_parallel,c).unwrap();
    assert!(matches!(result,SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearParallel,..}),"{result:?}");
}
#[test]
fn operand_symmetry_and_cylinder_reparameterization_preserve_points(){
    let origin=p(2.,-5.,3.);
    let cylinder=cylinder(origin,2.);
    let plane=plane(p(0.,1.,0.),p(0.,0.,1.),p(2.,-1.,0.));
    let ab=engine().intersect(plane,cylinder).unwrap();
    let ba=engine().intersect(cylinder,plane).unwrap();
    assert_eq!(ab,ba);
    certify(&ab,plane,cylinder);
    let Surface3::Cylinder(body)=cylinder else {unreachable!()};
    let phi=0.39_f64;
    let new_x=p(body.frame.x.x*phi.cos()+body.frame.y.x*phi.sin(),
        body.frame.x.y*phi.cos()+body.frame.y.y*phi.sin(),
        body.frame.x.z*phi.cos()+body.frame.y.z*phi.sin());
    let new_y=p(body.frame.y.x*phi.cos()-body.frame.x.x*phi.sin(),
        body.frame.y.y*phi.cos()-body.frame.x.y*phi.sin(),
        body.frame.y.z*phi.cos()-body.frame.x.z*phi.sin());
    let rotated=Surface3::Cylinder(CylinderSurface::new(
        Frame3::from_axes(origin,new_x,new_y,tol()).unwrap(),2.,tol()).unwrap());
    let reparam=engine().intersect(plane,rotated).unwrap();
    certify(&reparam,plane,rotated);
    for t in [0.1,0.7,2.1,4.3,5.9]{
        xyz(ab.curves()[0].evaluate(t).unwrap(),reparam.curves()[0].evaluate(t-phi).unwrap(),1e-10);
    }
}
#[test]
fn translated_analytic_sections_are_stable_when_representable(){
    let shift=p(1e9,-1e9,1e9);
    let original=engine().intersect(plane_y(3.),cylinder(p(0.,0.,0.),2.)).unwrap();
    let relocated_plane=plane(p(shift.x,shift.y+3.,shift.z),p(1.,0.,0.),p(0.,0.,-1.));
    let relocated_cylinder=cylinder(shift,2.);
    let relocated=engine().intersect(relocated_plane,relocated_cylinder).unwrap();
    certify(&relocated,relocated_plane,relocated_cylinder);
    for t in [0.3,1.4,3.5,5.7]{
        let before=original.curves()[0].evaluate(t).unwrap();
        let after=relocated.curves()[0].evaluate(t).unwrap();
        xyz(sub(after,shift),before,1e-6);
    }
}
#[test]
fn invalid_inputs_unsupported_pairs_and_tampered_curve_are_rejected(){
    let bad=Surface3::Cylinder(CylinderSurface{
        frame:Frame3::from_axes(p(0.,0.,0.),p(1.,0.,0.),p(0.,0.,-1.),tol()).unwrap(),
        radius:-1.,
    });
    assert!(engine().intersect(plane_x(0.),bad).is_err());
    let sphere=Surface3::Sphere(rustsolid::SphereSurface{
        frame:Frame3::from_axes(p(0.,0.,0.),p(1.,0.,0.),p(0.,1.,0.),tol()).unwrap(),radius:1.,
    });
    assert!(matches!(engine().intersect(sphere,plane_x(0.)),Err(GeometryError::InvalidEdit(_))));
    let swapped=engine().intersect(cylinder(p(0.,0.,0.),2.),plane_x(0.)).unwrap();
    certify(&swapped,cylinder(p(0.,0.,0.),2.),plane_x(0.));
    let forged=SurfaceIntersection::Curves{contact:IntersectionContact::Crossing,
        curves:vec![IntersectionCurve::Line(Line3{origin:p(20.,20.,20.),direction:p(0.,1.,0.)})]};
    assert!(forged.certify(plane_x(0.),cylinder(p(0.,0.,0.),2.),tol()).is_err());
}
