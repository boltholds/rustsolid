use rustsolid::{
    GeometryTolerance, Orientation, Point2, Point3, PredicateKernel, PredicateClassification,
    Solid, orient2d,
};

#[test]
fn near_cancelled_determinant_is_classified_by_robust_predicate() {
    // Naive (b-a)x(c-a) in f64 loses the unit determinant to cancellation.
    let a = [0.0, 0.0];
    let b = [100_000_000.0, 100_000_001.0];
    let c = [100_000_001.0, 100_000_002.0];
    assert_eq!(orient2d(a,b,c).unwrap(), Orientation::Negative);
    let kernel=PredicateKernel::new(GeometryTolerance::default()).unwrap();
    // The exact sign is negative, but this very thin shape is within
    // default CAD tolerance at 1e8 model-unit extent.
    assert_eq!(kernel.orient2d(a,b,c).unwrap(),PredicateClassification::WithinTolerance);
    let strict=PredicateKernel::new(GeometryTolerance {
        absolute_length:1e-16, relative_length:0.0, angular:1e-8,
    }).unwrap();
    assert_eq!(strict.orient2d(a,b,c).unwrap(),PredicateClassification::Negative);
    assert_eq!(strict.orient2d(a,c,b).unwrap(),PredicateClassification::Positive);
}

#[test]
fn geometric_tolerance_has_a_separate_near_degenerate_class() {
    let kernel=PredicateKernel::new(GeometryTolerance{
        absolute_length:1e-4, relative_length:0.0, angular:1e-8,
    }).unwrap();
    let a=[0.0,0.0]; let b=[1.0,0.0];
    assert_eq!(kernel.orient2d(a,b,[0.5,1e-5]).unwrap(),PredicateClassification::WithinTolerance);
    assert_eq!(kernel.orient2d(a,b,[0.5,1e-2]).unwrap(),PredicateClassification::Positive);
    assert!(kernel.orient2d(a,b,[f64::NAN,0.0]).is_err());
}

#[test]
fn exact_sign_predicates_are_used_by_profile_and_face_editing() {
    // A tiny but resolvable concavity must not silently reverse triangulation.
    let poly=[Point2{x:0.0,z:0.0},Point2{x:10.0,z:0.0},
        Point2{x:10.0,z:10.0},Point2{x:5.0,z:9.9},Point2{x:0.0,z:10.0}];
    let mut body=Solid::extrude_xz(&poly,2.0).unwrap();
    body.validate().unwrap();
    let area=99.5;
    assert!((body.mass.volume - area*2.0).abs() < 1e-8);
    let old=body.mesh.triangles.len();
    body.edit_atomic(|tx| {tx.split_edge(body_dummy_edge(),0.5)?;Ok(())}).unwrap();
    assert_eq!(body.mesh.triangles.len(),old+2);
    body.validate().unwrap();
}

fn body_dummy_edge()->rustsolid::EdgeId { rustsolid::EdgeId(0) }

#[test]
fn translate_a_small_polyhedron_far_from_origin() {
    let body=Solid::block(Point3{x:1e12,y:1e12,z:1e12},1.0,1.0,1.0).unwrap();
    body.validate().unwrap();
}
