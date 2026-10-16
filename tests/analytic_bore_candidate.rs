//! G4b2b1 acceptance: real curved B-reps, multiple source face-pair
//! closed loops, physical circle coedges and a nonmaterial periodic seam.
use rustsolid::{
    AnalyticBoreCandidate,AnalyticBoreCandidateResult,AnalyticBoreKind,
    AnalyticBoreHold,AnalyticBoreFaceUse,BrepBody,BrepOrigin,
    BooleanOperation,CoedgeId,Curve2,Curve3,CylindricalBrep,EdgeId,
    FaceId,Frame3,GeometryTolerance,GraphNodeKind,GraphState,
    Point3,Solid,Surface3,
};
use std::f64::consts::TAU;

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->BrepBody {
    BrepBody::from(Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap())
}
fn cylinder(bottom:f64,height:f64,radius:f64,x:f64,z:f64)->BrepBody {
    BrepBody::from(CylindricalBrep::upright(
        p(x,bottom,z),radius,height,GeometryTolerance::default()).unwrap())
}
fn through()->(BrepBody,BrepBody){
    (block(),cylinder(-1.,4.,1.,0.,0.))
}
fn blind()->(BrepBody,BrepBody){
    (block(),cylinder(1.,2.,1.,0.,0.))
}
fn build(a:&BrepBody,b:&BrepBody)->AnalyticBoreCandidateResult {
    AnalyticBoreCandidate::build(a,b,BooleanOperation::Difference,
        GeometryTolerance::default()).unwrap()
}
fn ready(result:AnalyticBoreCandidateResult)->AnalyticBoreCandidate{
    let AnalyticBoreCandidateResult::Ready(c)=result else{
        panic!("expected source-certified curved result: {result:?}");
    };
    c
}
fn close(a:f64,b:f64){
    assert!((a-b).abs()<1e-8,"actual {a} expected {b}");
}

#[test]
fn real_through_bore_has_two_independent_source_circle_graphs_and_genus_one(){
    let (a,b)=through();
    let original_a=format!("{a:?}");
    let original_b=format!("{b:?}");
    let candidate=ready(build(&a,&b));
    assert_eq!(candidate.kind(),AnalyticBoreKind::Through);
    assert_eq!(candidate.brep().origin,BrepOrigin::AnalyticThroughHole);
    let m=candidate.brep();
    m.validate().unwrap();
    assert_eq!(m.summary().genus,1);
    assert_eq!(m.euler_characteristic(),0);
    assert_eq!(m.faces.len(),7);
    assert_eq!(m.summary().cylindrical_faces,1);
    assert_eq!(m.summary().seam_edges,1);
    assert_eq!(m.summary().closed_edges,2);
    assert_eq!(candidate.face_coverage().len(),9);
    assert_eq!(candidate.circle_graphs().len(),2);
    for (i,graph) in candidate.circle_graphs().iter().enumerate(){
        assert_eq!(graph.state,GraphState::Ready);
        assert_eq!(graph.source_faces[0].face,FaceId(i as u32));
        assert_eq!(graph.source_faces[1].face,FaceId(2));
        assert_eq!(graph.nodes.len(),1);
        assert_eq!(graph.arcs.len(),1);
        assert_eq!(graph.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
        assert!(graph.arcs[0].closed);
        assert_eq!(graph.arcs[0].from,graph.arcs[0].to);
        assert_eq!(graph.arcs[0].parameter_end-graph.arcs[0].parameter_start,TAU);
    }
    let rings=candidate.rings();
    assert_eq!(rings[0].result_edge,EdgeId(12));
    assert_eq!(rings[1].result_edge,EdgeId(13));
    assert_eq!(rings[0].source_box_face,Some(FaceId(0)));
    assert_eq!(rings[1].source_box_face,Some(FaceId(1)));
    assert!(rings.iter().all(|ring|ring.closed_graph_witness));
    for ring in rings{
        let edge=&m.edges[ring.result_edge.0 as usize];
        assert_eq!(edge.start,edge.end);
        assert!(matches!(m.edge_curve(edge.id),Some(Curve3::Circle(_))));
        for &c in &ring.result_coedges{
            assert!(matches!(m.coedge_pcurve(c),Some(Curve2::Circle(_)|Curve2::Line(_))));
        }
    }
    assert!(!candidate.can_apply());
    assert!(candidate.closed_shell_proven());
    close(m.mass.volume,144.-4.*std::f64::consts::PI);
    close(m.mass.surface_area,168.+6.*std::f64::consts::PI);
    candidate.validate_against(&a,&b,GeometryTolerance::default()).unwrap();
    assert_eq!(format!("{a:?}"),original_a);
    assert_eq!(format!("{b:?}"),original_b);
}

#[test]
fn real_blind_bore_has_partial_depth_wall_and_source_cylinder_floor(){
    let (a,b)=blind();
    let candidate=ready(build(&a,&b));
    let m=candidate.brep();m.validate().unwrap();
    assert_eq!(candidate.kind(),AnalyticBoreKind::Blind);
    assert_eq!(m.origin,BrepOrigin::AnalyticBlindHole);
    assert_eq!(m.summary().genus,0);
    assert_eq!(m.euler_characteristic(),2);
    assert_eq!((m.faces.len(),m.edges.len(),m.coedges.len()),(8,15,30));
    assert_eq!(candidate.circle_graphs().len(),1);
    assert_eq!(candidate.circle_graphs()[0].source_faces[0].face,FaceId(1));
    assert_eq!(candidate.rings()[0].planar_result_face,FaceId(7));
    assert_eq!(candidate.rings()[0].source_box_face,None);
    assert_eq!(candidate.rings()[0].source_cylinder_face,FaceId(0));
    assert!(!candidate.rings()[0].closed_graph_witness);
    assert_eq!(candidate.rings()[1].planar_result_face,FaceId(1));
    assert!(candidate.rings()[1].closed_graph_witness);
    assert!(matches!(m.face_surface(FaceId(7)),Some(Surface3::Plane(_))));
    assert!(candidate.face_coverage().iter().any(|f|
        f.source_face==FaceId(0)&&f.result_face==Some(FaceId(7))&&
        f.use_kind==AnalyticBoreFaceUse::BoreFloorFromCutterCap));
    assert_eq!(candidate.face_coverage().iter().filter(|f|
        f.use_kind==AnalyticBoreFaceUse::OriginalWithCircularOpening).count(),1);
    close(m.mass.volume,144.-2.*std::f64::consts::PI);
    close(m.mass.surface_area,168.+4.*std::f64::consts::PI);
    candidate.validate_against(&a,&b,GeometryTolerance::default()).unwrap();
}

#[test]
fn periodic_cylinder_seam_has_two_distinct_uv_charts_but_one_physical_edge(){
    for (a,b) in [through(),blind()]{
        let c=ready(build(&a,&b));
        let seam=c.seam();
        let m=c.brep();
        assert_eq!(seam.source_cylinder_face,FaceId(2));
        assert_eq!(seam.result_edge,EdgeId(14));
        assert_eq!(seam.uses,[CoedgeId(28),CoedgeId(29)]);
        assert!((seam.uv_u[0]-seam.uv_u[1]).abs()-TAU<1e-8);
        let edge=&m.edges[seam.result_edge.0 as usize];
        assert_ne!(edge.start,edge.end);
        assert_eq!(m.coedges[28].face,FaceId(6));
        assert_eq!(m.coedges[29].face,FaceId(6));
        assert_eq!(m.coedges[28].twin,CoedgeId(29));
        assert_eq!(m.coedges[29].twin,CoedgeId(28));
        assert!(!c.can_apply());
    }
}

#[test]
fn off_center_through_and_partial_blind_bores_preserve_analytic_mass(){
    let a=block();
    for (bottom,height) in [(-1.0,4.0),(0.5,2.5)]{
        let b=cylinder(bottom,height,0.75,0.4,-0.6);
        let c=ready(build(&a,&b));
        let expected=144.-std::f64::consts::PI*0.75_f64.powi(2)*height;
        close(c.brep().mass.volume,expected);
        assert!(!c.can_apply());
    }
}

#[test]
fn rotated_cylinder_uv_frame_does_not_change_true_physical_circle_geometry(){
    let a=block();
    let phi=0.37_f64;
    let frame=Frame3::from_axes(
        p(0.2,-1.,-0.1),
        p(phi.cos(),0.,-phi.sin()),
        p(-phi.sin(),0.,-phi.cos()),
        GeometryTolerance::default(),
    ).unwrap();
    let b=BrepBody::from(CylindricalBrep::new(
        frame,0.75,4.,GeometryTolerance::default()).unwrap());
    let c=ready(build(&a,&b));
    assert_eq!(c.circle_graphs().len(),2);
    assert_eq!(c.brep().summary().genus,1);
    c.validate_against(&a,&b,GeometryTolerance::default()).unwrap();
}

#[test]
fn nonthrough_nonblind_or_nearly_touching_side_is_never_a_valid_candidate(){
    let a=block();
    for b in [
        cylinder(-2.,3.,1.,0.,0.),
        cylinder(-1.,4.,1.,2.1,0.),
        cylinder(1.,2.,1.,2.1,0.),
    ]{
        assert!(matches!(build(&a,&b),AnalyticBoreCandidateResult::Held(
            AnalyticBoreHold::UnsupportedPartialIntersection)));
    }
}
#[test]
fn unsupported_operation_or_source_kind_fails_closed(){
    let (a,b)=through();
    let union=AnalyticBoreCandidate::build(&a,&b,BooleanOperation::Union,
        GeometryTolerance::default()).unwrap();
    assert!(matches!(union,AnalyticBoreCandidateResult::Held(
        AnalyticBoreHold::UnsupportedOperation)));
    let b_box=BrepBody::from(Solid::block(p(-1.,0.,-1.),2.,1.,2.).unwrap());
    assert!(matches!(build(&a,&b_box),AnalyticBoreCandidateResult::Held(
        AnalyticBoreHold::UnsupportedOperand)));
}

#[test]
fn stale_source_revision_geometry_and_tolerance_do_not_reuse_curved_preview(){
    let (a,b)=through();
    let c=ready(build(&a,&b));
    let different=cylinder(-1.,4.,1.1,0.,0.);
    assert!(c.validate_against(&a,&different,GeometryTolerance::default()).is_err());
    let mut edited=a.clone();
    if let BrepBody::Polyhedral(ref mut src)=edited{src.revision+=1;}
    assert!(c.validate_against(&edited,&b,GeometryTolerance::default()).is_err());
    let policy=GeometryTolerance{
        absolute_length:1e-8,relative_length:1e-12,angular:1e-8,
    };
    assert!(c.validate_against(&a,&b,policy).is_err());
    assert!(matches!(AnalyticBoreCandidate::build(
        &a,&b,BooleanOperation::Difference,policy).unwrap(),
        AnalyticBoreCandidateResult::Held(AnalyticBoreHold::TolerancePolicyMismatch)));
}
