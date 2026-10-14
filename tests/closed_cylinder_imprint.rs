//! Independently derived closed oblique imprint on a cylindrical B-rep.
//! The imprint changes only its face partition and carries no Boolean material
//! deletion/creation. The input TopologyDeltaPlan remains preview-only.
use std::f64::consts::TAU;
use rustsolid::{
    BrepModel, ClosedCylinderImprint, CoedgeId, Curve2, Curve3,
    CylindricalBrep, FaceId, Frame3, GeometryTolerance, GraphNodeKind,
    GraphState, ImprintHold, ImprintResult, IntersectionGraph,
    IntersectionGraphBuilder, ModelUnits, PlanRequirement, PlanState,
    Point3, Solid, Surface3, TopologyDeltaPlan, TopologyDeltaPlanner,
};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tol()->GeometryTolerance {GeometryTolerance::default()}
fn distance(a:Point3,b:Point3)->f64 {
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)
}

fn setup(plane_y:f64)->(
    BrepModel,FaceId,BrepModel,CylindricalBrep,IntersectionGraph,TopologyDeltaPlan
){
    let angle=0.3_f64;
    let (s,c)=angle.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(c,s,0.),p(0.,0.,-1.),tol()).unwrap();
    let source=CylindricalBrep::new(frame,2.,5.,tol()).unwrap();
    let side=BrepModel::from_cylinder(&source).unwrap();
    let block=Solid::block(p(-7.,plane_y,-7.),14.,1.,14.).unwrap();
    let plane=BrepModel::from_polyhedral(&block,ModelUnits::default()).unwrap();
    let face=plane.faces.iter().find_map(|f|{
        let Surface3::Plane(surface)=plane.face_surface(f.id)? else{return None};
        ((surface.frame.origin.y-plane_y).abs()<1e-10 &&
            surface.frame.normal.y.abs()>0.9).then_some(f.id)
    }).unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&plane,face,&side,FaceId(2)).unwrap();
    let plan=TopologyDeltaPlanner::new().plan(&graph,&plane,face,&side,FaceId(2)).unwrap();
    (plane,face,side,source,graph,plan)
}

#[test]
fn closed_oblique_imprint_creates_one_physical_ellipse_and_two_opposite_uv_coedges(){
    let (plane,_,side,source,graph,plan)=setup(2.5);
    assert_eq!(graph.state,GraphState::Ready);
    assert_eq!(graph.arcs.len(),1);
    assert!(graph.arcs[0].closed);
    assert_eq!(graph.nodes[0].kind,GraphNodeKind::PeriodicAnchor);
    assert_eq!(plan.state,PlanState::PreviewOnly);
    assert!(!plan.can_apply);
    assert!(plan.requirements.contains(&PlanRequirement::GeneralizedBrepEditor));

    // A single intersection curve has independent, *persistable* surface UV
    // carriers. A pcurve on a plane face is not a coedge on the cylinder.
    let [planar,cylindrical]=graph.arcs[0].pcurves(&plane,&side).unwrap();
    let cp=planar.to_curve2().unwrap();
    let cc=cylindrical.to_curve2().unwrap();
    assert!(matches!(cp,Curve2::Harmonic(_)));
    assert!(matches!(cc,Curve2::CylinderSection(_)));
    for t in [0.,0.125,0.25,0.375,0.5,0.75,1.] {
        let u=t*TAU;
        let first=cp.evaluate(u).unwrap();
        let second=cc.evaluate(u).unwrap();
        let q=graph.arcs[0].evaluate(t).unwrap();
        let a=planar.surface.evaluate(first.u,first.v).unwrap().point;
        let b=cylindrical.surface.evaluate(second.u,second.v).unwrap().point;
        assert!(distance(q,a)<1e-8);
        assert!(distance(q,b)<1e-8);
    }

    let mut editor=ClosedCylinderImprint::new(source.clone()).unwrap();
    let before=editor.model().mass;
    let commit=editor.imprint(&graph,&plan,&plane,&side).unwrap();
    let ImprintResult::Committed(changes)=commit else{panic!("{commit:?}");};
    assert_eq!(changes.revision_before,0);
    assert_eq!(changes.revision_after,1);
    assert_eq!(editor.revision(),1);
    let result=editor.model();
    result.validate().unwrap();
    assert_eq!(result.mass,before);
    assert_eq!(result.bbox,side.bbox);
    assert_eq!(result.euler_characteristic(),2);
    assert_eq!((result.vertices.len(),result.edges.len(),
        result.coedges.len(),result.loops.len(),result.faces.len()),(3,5,10,4,4));
    assert_eq!(changes.new_ellipse_edge.0,4);
    assert_eq!(changes.new_upper_face.0,3);
    assert!(matches!(result.edge_curve(changes.new_ellipse_edge),
        Some(Curve3::Harmonic(_))));
    let [c0,c1]=changes.ellipse_coedges;
    assert_eq!(result.coedges[c0.0 as usize].twin,c1);
    assert_eq!(result.coedges[c1.0 as usize].twin,c0);
    assert_ne!(result.coedges[c0.0 as usize].face,
        result.coedges[c1.0 as usize].face);
    assert_ne!(result.coedges[c0.0 as usize].reversed,
        result.coedges[c1.0 as usize].reversed);
    assert_eq!(result.coedges[CoedgeId(4).0 as usize].face,FaceId(2));
    assert_eq!(result.coedges[CoedgeId(6).0 as usize].face,FaceId(3));

    // Check actual canonical coedge pcurves, not preview-only witnesses.
    let edge=result.edge_curve(changes.new_ellipse_edge).unwrap();
    for fraction in [0.,0.13,0.31,0.6,0.84,1.] {
        let q=edge.evaluate(fraction*TAU).unwrap();
        for &co in &[c0,c1] {
            let g=&result.geometry.coedges[co.0 as usize];
            let Curve2::CylinderSection(pc)=result.coedge_pcurve(co).unwrap()
                else{panic!("missing persisted cylindrical section trim");};
            let t=if co==c1 {TAU*(1.-fraction)}else{fraction*TAU};
            let uv=pc.evaluate(t).unwrap();
            let surf=result.face_surface(result.coedges[co.0 as usize].face).unwrap();
            let other=surf.evaluate(uv.u,uv.v).unwrap().point;
            assert!(distance(q,other)<1e-8);
            assert_eq!(g.domain,rustsolid::ParameterRange::Bounded{start:0.,end:TAU});
        }
    }
    // As a source-preserving imprint there is no Boolean material change.
    assert_eq!(result.mass.volume,source.mass.volume);
    assert!(!plan.can_apply);

    // A harmonic lower/upper side boundary is not a rectangular trim
    // domain. The old preview clipper must fail closed, never classify it
    // using only its bounding UV box or number of coedges.
    assert!(rustsolid::TrimmedIntersectionEngine::new(tol()).unwrap()
        .intersect_faces(&plane,graph.source_faces[0].face,
            result,FaceId(2)).is_err());
    // The carrier editor validates the persisted nonlinear trims and
    // their true 3D curve owners independently of the original preview.
    rustsolid::BrepCarrierEditor::new(result.clone()).unwrap();

    let mut damaged=result.clone();
    let id=damaged.geometry.coedges[c0.0 as usize].pcurve.0 as usize;
    let Curve2::CylinderSection(mut pc)=damaged.geometry.curves2[id]
        else{panic!("ellipse pcurve has wrong type");};
    pc.axial_cosine+=0.1;
    damaged.geometry.curves2[id]=Curve2::CylinderSection(pc);
    assert!(damaged.validate().is_err(),
        "tampering with an owned analytic trim must fail B-rep validation");
}

#[test]
fn closed_imprint_undo_redo_are_atomic_and_revisioned(){
    let (plane,_,side,source,graph,plan)=setup(2.5);
    let mut editor=ClosedCylinderImprint::new(source).unwrap();
    let base=format!("{:?}",editor.model());
    assert!(matches!(editor.imprint(&graph,&plan,&plane,&side).unwrap(),ImprintResult::Committed(_)));
    let imprinted=format!("{:?}",editor.model());
    assert_ne!(base,imprinted);
    assert_eq!((editor.undo_depth(),editor.redo_depth()),(1,0));
    assert_eq!(editor.imprint(&graph,&plan,&plane,&side).unwrap(),
        ImprintResult::Held(ImprintHold::AlreadyImprinted));
    assert_eq!(editor.revision(),1);
    assert!(editor.undo().unwrap());
    assert_eq!(editor.revision(),2);
    assert_eq!(format!("{:?}",editor.model()),base);
    assert_eq!((editor.undo_depth(),editor.redo_depth()),(0,1));
    assert!(editor.redo().unwrap());
    assert_eq!(editor.revision(),3);
    assert_eq!(format!("{:?}",editor.model()),imprinted);
    assert_eq!((editor.undo_depth(),editor.redo_depth()),(1,0));
    assert!(!editor.redo().unwrap());
    editor.model().validate().unwrap();
}

#[test]
fn unsupported_and_mutated_intent_hold_or_fail_without_partial_topology(){
    let (plane,fa,side,source,graph,plan)=setup(2.5);
    let mut editor=ClosedCylinderImprint::new(source).unwrap();
    let original=format!("{:?}",editor.model());
    // Reverse source ordering is not within the bounded executor contract.
    let reverse=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&side,FaceId(2),&plane,fa).unwrap();
    let inverse=TopologyDeltaPlanner::new().plan(&reverse,&side,FaceId(2),&plane,fa).unwrap();
    assert_eq!(editor.imprint(&reverse,&inverse,&side,&plane).unwrap(),
        ImprintResult::Held(ImprintHold::UnsupportedSource));

    // A forged topology delta must fail provenance validation, never commit.
    let mut tampered=plan.clone();
    tampered.arcs[0].parameter_end-=0.02;
    assert!(editor.imprint(&graph,&tampered,&plane,&side).is_err());
    assert_eq!(editor.revision(),0);
    assert_eq!(editor.undo_depth(),0);
    assert_eq!(format!("{:?}",editor.model()),original);
}

#[test]
fn near_cap_and_unsupported_geometry_do_not_publish_partial_imprint(){
    let (plane,_,side,source,graph,plan)=setup(2.5);
    let mut editor=ClosedCylinderImprint::new(source.clone()).unwrap();
    let original=format!("{:?}",editor.model());
    // Source is pinned: a different cylinder cannot silently reuse a graph.
    let other=CylindricalBrep::new(source.frame,2.,6.,tol()).unwrap();
    let outsider=ClosedCylinderImprint::new(other).unwrap();
    assert_ne!(format!("{:?}",outsider.model()),original);
    let mut tampered=side.clone();
    tampered.geometry.curves2[0]=Curve2::Line(rustsolid::Line2{
        origin:rustsolid::Point2Param{u:0.0,v:0.0},
        delta:rustsolid::Point2Param{u:0.0,v:0.0}});
    assert!(editor.imprint(&graph,&plan,&plane,&tampered).is_err());
    assert_eq!(format!("{:?}",editor.model()),original);

    // Near-cap plane has one section outside the finite side; it is not
    // eligible for closed-imprint material preservation.
    let (plane2,_,side2,_,graph2,plan2)=setup(4.95);
    let held=editor.imprint(&graph2,&plan2,&plane2,&side2).unwrap();
    assert!(matches!(held,ImprintResult::Held(_)));
    assert_eq!(editor.revision(),0);
    assert_eq!(editor.undo_depth(),0);
}
