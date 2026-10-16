//! G4b2a: all original axis-aligned box faces accounted for, with a
//! real validated canonical B-rep, physical twin coedges and vertex fans.
//! This is a restricted preview, not a general curved Boolean or G5 commit.
use rustsolid::{
    BrepBody,BooleanOperation,BooleanOperand,CanonicalBoxCandidate,
    CanonicalBoxCandidateHold,CanonicalBoxCandidateResult,
    GeometryTolerance,Point2,Point3,Solid,SourceFaceDisposition,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_(x:f64,y:f64,z:f64,w:f64,h:f64,d:f64)->BrepBody{
    BrepBody::from(Solid::block(p(x,y,z),w,h,d).unwrap())
}
fn ab()->(BrepBody,BrepBody){
    (box_(0.,0.,0.,2.,3.,4.),box_(1.,1.,1.,2.,3.,4.))
}
fn ready(result:CanonicalBoxCandidateResult)->CanonicalBoxCandidate{
    let CanonicalBoxCandidateResult::Ready(candidate)=result else{
        panic!("expected a certified box-only B-rep candidate: {result:?}");
    };
    candidate
}
fn candidate(
    a:&BrepBody,b:&BrepBody,op:BooleanOperation,
)->CanonicalBoxCandidate{
    ready(CanonicalBoxCandidate::build(a,b,op,GeometryTolerance::default()).unwrap())
}
fn eq(a:f64,b:f64){assert!((a-b).abs()<=1e-8,"{a} vs {b}");}

#[test]
fn three_csg_operations_generate_real_closed_canonical_brep_candidates(){
    let (a,b)=ab();
    let before_a=format!("{a:?}");
    let before_b=format!("{b:?}");
    for (op,volume) in [
        (BooleanOperation::Union,42.),
        (BooleanOperation::Intersection,6.),
        (BooleanOperation::Difference,18.),
    ]{
        let preview=candidate(&a,&b,op);
        assert!(preview.closed_shell_proven());
        assert!(!preview.can_apply());
        preview.brep().validate().unwrap();
        eq(preview.brep().mass.volume,volume);
        assert_eq!(preview.brep().shells.len(),1);
        assert!(preview.brep().shells[0].closed);
        assert_eq!(preview.source_faces().len(),12);
        assert_eq!(preview.face_provenance().len(),preview.brep().faces.len());
        assert_eq!(preview.vertices().len(),preview.brep().vertices.len());
        assert_eq!(preview.radial_topology().vertices,preview.brep().vertices.len());
        assert_eq!(preview.radial_topology().edges,preview.brep().edges.len());
        assert_eq!(preview.radial_topology().material_face_cells,preview.brep().faces.len());
        assert_eq!(preview.radial_topology().shell_components,1);
        assert_eq!(preview.radial_topology().euler_characteristic,2);
        assert_eq!(preview.radial_topology().total_genus,0);
        preview.validate_against(&a,&b,GeometryTolerance::default()).unwrap();
    }
    assert_eq!(format!("{a:?}"),before_a);
    assert_eq!(format!("{b:?}"),before_b);
}

#[test]
fn untouched_deleted_and_partitioned_faces_have_exact_area_coverage(){
    let (a,b)=ab();
    let plan=candidate(&a,&b,BooleanOperation::Difference);
    let coverage=plan.source_faces();
    assert!(coverage.iter().any(|face|
        face.disposition==SourceFaceDisposition::RetainedUncut));
    assert!(coverage.iter().any(|face|
        face.disposition==SourceFaceDisposition::Discarded));
    assert!(coverage.iter().any(|face|
        face.disposition==SourceFaceDisposition::Partitioned));
    let mut area=0.0;
    for face in coverage {
        assert!(face.original_area>0.0);
        assert!(face.total_tiles>=face.exposed_tiles);
        assert!(face.exposed_area>=0.0);
        assert!(face.exposed_area<=face.original_area+1e-9);
        area+=face.exposed_area;
    }
    eq(area,plan.brep().mass.surface_area);
    assert!(plan.face_provenance().iter().any(|face|
        face.operand==BooleanOperand::Right&&face.reversed));
    assert!(plan.face_provenance().iter().any(|face|
        face.operand==BooleanOperand::Left&&!face.reversed));
}

#[test]
fn canonical_result_vertex_ids_are_exact_source_plane_lattice_intersections(){
    let (a,b)=ab();
    let plan=candidate(&a,&b,BooleanOperation::Intersection);
    let axes=[
        [0.0,1.0,2.0,3.0],
        [0.0,1.0,3.0,4.0],
        [0.0,1.0,4.0,5.0],
    ];
    let mut grid=std::collections::BTreeSet::new();
    for certificate in plan.vertices(){
        let v=&plan.brep().vertices[certificate.result_vertex.0 as usize];
        assert_eq!(v.position.x,axes[0][certificate.coordinate_index[0] as usize]);
        assert_eq!(v.position.y,axes[1][certificate.coordinate_index[1] as usize]);
        assert_eq!(v.position.z,axes[2][certificate.coordinate_index[2] as usize]);
        assert!(grid.insert(certificate.coordinate_index),
            "no coordinate snapping or duplicate topological vertex");
    }
}

#[test]
fn unsupported_nonbox_polyhedral_cannot_be_replaced_by_its_bounding_box(){
    let (a,_) = ab();
    let concave=Solid::extrude_xz(&[
        Point2{x:0.,z:0.},Point2{x:3.,z:0.},Point2{x:3.,z:1.},
        Point2{x:1.,z:1.},Point2{x:1.,z:3.},Point2{x:0.,z:3.},
    ],2.).unwrap();
    let response=CanonicalBoxCandidate::build(
        &a,&BrepBody::from(concave),BooleanOperation::Difference,
        GeometryTolerance::default()).unwrap();
    assert!(matches!(response,CanonicalBoxCandidateResult::Held(
        CanonicalBoxCandidateHold::UnsupportedOperand)));
}

#[test]
fn near_coincident_and_exactly_coincident_support_planes_are_held(){
    let a=box_(0.,0.,0.,2.,3.,4.);
    for b in [
        box_(2.,1.,1.,2.,3.,4.),
        box_(2.+1e-10,1.,1.,2.,3.,4.),
    ]{
        let response=CanonicalBoxCandidate::build(
            &a,&b,BooleanOperation::Difference,
            GeometryTolerance::default()).unwrap();
        assert!(matches!(response,CanonicalBoxCandidateResult::Held(
            CanonicalBoxCandidateHold::NearCoincidentSourcePlanes)));
    }
}

#[test]
fn full_disjoint_intersection_is_a_typed_empty_held_result(){
    let a=box_(0.,0.,0.,1.,1.,1.);
    let b=box_(5.,5.,5.,1.,1.,1.);
    let response=CanonicalBoxCandidate::build(
        &a,&b,BooleanOperation::Intersection,
        GeometryTolerance::default()).unwrap();
    assert!(matches!(response,CanonicalBoxCandidateResult::Held(
        CanonicalBoxCandidateHold::EmptyResult)));
}

#[test]
fn source_or_tolerance_change_invalidates_built_geometry_candidate(){
    let (a,b)=ab();
    let plan=candidate(&a,&b,BooleanOperation::Union);
    let shifted=box_(1.2,1.,1.,2.,3.,4.);
    assert!(plan.validate_against(&a,&shifted,GeometryTolerance::default()).is_err());
    // Recomputed source lattice/graph changes and invalidates the draft.
    let stricter=GeometryTolerance{
        absolute_length:1e-7,
        relative_length:1e-10,
        angular:1e-8,
    };
    // The same exact geometry remains buildable under a new tolerance,
    // but the older source-bound candidate cannot be silently reused.
    assert!(plan.validate_against(&a,&b,stricter).is_err());
    assert_eq!(plan.construction_tolerance(),GeometryTolerance::default());
    assert_eq!(plan.source_revisions(),[0,0]);
    let mut stale=b.clone();
    if let BrepBody::Polyhedral(ref mut solid)=stale {
        solid.revision+=1;
    }
    assert!(plan.validate_against(&a,&stale,GeometryTolerance::default()).is_err(),
        "same geometry with a changed edit revision cannot reuse this preview");
    let fresh=CanonicalBoxCandidate::build(
        &a,&b,BooleanOperation::Union,stricter).unwrap();
    assert!(matches!(fresh,CanonicalBoxCandidateResult::Ready(_)));
}

#[test]
fn one_intersecting_axis_aligned_fixture_has_a_complete_result_shell(){
    let a=box_(0.,0.,0.,4.,4.,4.);
    let b=box_(2.,1.,1.,4.,2.,2.);
    let plan=candidate(&a,&b,BooleanOperation::Difference);
    eq(plan.brep().mass.volume,64.-8.);
    assert!(plan.closed_shell_proven());
    assert!(!plan.can_apply());
    assert_eq!(plan.radial_topology().total_genus,0);
    plan.validate_against(&a,&b,GeometryTolerance::default()).unwrap();
}
