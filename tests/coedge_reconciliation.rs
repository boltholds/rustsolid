//! A physical intersection arc reconciles exactly two source-owned
//! analytic pcurves. The candidate is not yet assigned committed B-rep IDs.
use std::f64::consts::TAU;
use rustsolid::{
    BrepBody, CoedgeReconciliationPlan, Curve2, Curve3, CylindricalBrep,
    EdgeReconciliationResult, FaceId, Frame3, GeometryTolerance,
    GraphNodeKind, GraphOperand, GraphState, HoleParameters,
    IntersectionGraphBuilder, Point3, Solid, Surface3, ThroughHoleBrep,
    UvFaceArrangement, UvFaceRegionPlan,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn dist(a:Point3,b:Point3)->f64{(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn plane_face(a:&rustsolid::BrepModel,y:f64)->FaceId{
    a.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=a.face_surface(f.id)? else{return None};
        ((s.frame.origin.y-y).abs()<1e-9&&s.frame.normal.y.abs()>0.9)
            .then_some(f.id)
    }).unwrap()
}
#[test]
fn two_intersection_arcs_reconcile_owning_face_pcurves_and_both_uv_cells(){
    let source=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    let annulus=BrepBody::from(ThroughHoleBrep::from_box(&source,
        HoleParameters{center_x:0.,center_z:0.,radius:1.}).unwrap());
    let plane=BrepBody::from(Solid::block(p(0.,4.,-4.),2.,2.,8.).unwrap());
    let left=annulus.shared().unwrap();let right=plane.shared().unwrap();
    let x_face=right.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=right.face_surface(f.id)? else{return None};
        ((s.frame.origin.x).abs()<1e-9&&s.frame.normal.x.abs()>0.9)
            .then_some(f.id)
    }).unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&left,FaceId(1),&right,x_face).unwrap();
    let EdgeReconciliationResult::Ready(plan)=
        CoedgeReconciliationPlan::build(&left,&right,&graph,tol()).unwrap()
        else {panic!("line segments should have two physical pcurves")};
    assert_eq!(plan.edges.len(),2);
    assert!(!plan.can_apply);
    for edge in &plan.edges{
        assert!(matches!(edge.curve,Curve3::Line(_)));
        assert_eq!(edge.uses[0].source.operand,GraphOperand::First);
        assert_eq!(edge.uses[1].source.operand,GraphOperand::Second);
        assert!(!edge.uses[0].reversed&&edge.uses[1].reversed);
        assert!(!edge.virtual_periodic_anchor);
        assert!(matches!(edge.uses[0].curve2,Curve2::Line(_)));
        assert!(matches!(edge.uses[1].curve2,Curve2::Line(_)));
        for f in [0.,0.125,0.5,0.875,1.] {
            let a=edge.uses[0].point(f).unwrap();
            let b=edge.uses[1].point(1.-f).unwrap();
            let original=edge.curve.evaluate(match edge.parameter_domain{
                rustsolid::ParameterRange::Bounded{start,end}=>start+(end-start)*f,
                _=>panic!("unexpected unbounded curve"),
            }).unwrap();
            assert!(dist(a,b)<1e-8);
            assert!(dist(a,original)<1e-8);
        }
    }
    plan.validate_against(&left,&right,&graph,tol()).unwrap();
    let UvFaceArrangement::Ready(cells)=UvFaceRegionPlan::build(
        &left,&right,&graph,GraphOperand::First,tol()).unwrap()
        else{panic!("missing full planar UV arrangement");};
    assert_eq!(cells.cut_twins.len(),plan.edges.len());
    for twin in &cells.cut_twins{
        let edge=plan.edges.iter().find(|e|e.source_arc==twin.arc).unwrap();
        let surface=left.face_surface(FaceId(1)).unwrap();
        let p0=surface.evaluate(twin.positive_start.u,twin.positive_start.v).unwrap().point;
        let p1=surface.evaluate(twin.positive_end.u,twin.positive_end.v).unwrap().point;
        assert!(
            (dist(p0,edge.start)<1e-8 && dist(p1,edge.end)<1e-8) ||
            (dist(p1,edge.start)<1e-8 && dist(p0,edge.end)<1e-8)
        );
    }
    let mut forged=plan.clone();
    forged.edges[0].start.x+=0.1;
    assert!(forged.validate_against(&left,&right,&graph,tol()).is_err());
}
#[test]
fn closed_periodic_ellipse_has_a_virtual_parameter_anchor_not_material_vertex(){
    let tilt=0.3_f64;
    let (s,c)=tilt.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(c,s,0.),p(0.,0.,-1.),tol()).unwrap();
    let cylinder=BrepBody::from(CylindricalBrep::new(frame,2.,5.,tol()).unwrap());
    let horizontal=BrepBody::from(Solid::block(p(-7.,2.5,-7.),14.,1.,14.).unwrap());
    let left=horizontal.shared().unwrap();let right=cylinder.shared().unwrap();
    let face=plane_face(&left,2.5);
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&left,face,&right,FaceId(2)).unwrap();
    assert_eq!(graph.state,GraphState::Ready);
    assert_eq!(graph.arcs.len(),1);
    assert_eq!(graph.nodes.len(),1);
    assert_eq!(graph.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
    let EdgeReconciliationResult::Ready(owners)=
        CoedgeReconciliationPlan::build(&left,&right,&graph,tol()).unwrap()
        else{panic!("closed section not reconciled")};
    assert_eq!(owners.edges.len(),1);
    let edge=&owners.edges[0];
    assert!(edge.virtual_periodic_anchor);
    assert!(matches!(edge.curve,Curve3::Harmonic(_)));
    assert!(matches!(edge.uses[0].curve2,Curve2::Harmonic(_)));
    assert!(matches!(edge.uses[1].curve2,Curve2::CylinderSection(_)));
    assert_eq!(edge.origin_node,edge.destination_node);
    assert!(dist(edge.start,edge.end)<1e-8);
    assert!(!owners.can_apply);
}
#[test]
fn open_seam_crossing_keeps_one_physical_arc_and_unwrapped_uv(){
    let block=BrepBody::from(Solid::block(p(-2.5,2.,-2.2),5.,1.,3.2).unwrap());
    let left=block.shared().unwrap();let face=plane_face(&left,2.);
    let tilt=0.4_f64.atan();let (si,co)=tilt.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(co,si,0.),p(0.,0.,-1.),tol()).unwrap();
    let cylinder=BrepBody::from(CylindricalBrep::new(frame,2.,2.4,tol()).unwrap());
    let right=cylinder.shared().unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&left,face,&right,FaceId(2)).unwrap();
    let EdgeReconciliationResult::Ready(plan)=
        CoedgeReconciliationPlan::build(&left,&right,&graph,tol()).unwrap()
        else{panic!("periodic crossing not reconciled")};
    assert_eq!(plan.edges.len(),1);
    assert!(!plan.edges[0].virtual_periodic_anchor);
    let edge=&plan.edges[0];
    let front=edge.uses[1].evaluate(1.).unwrap();
    let back=edge.uses[1].evaluate(0.).unwrap();
    assert!(front.u<TAU && back.u>TAU);
    assert_eq!(graph.nodes.len(),2);
}
