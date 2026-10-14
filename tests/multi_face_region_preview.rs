//! Multiple source-face pair preview with strict all-or-nothing Held states.
//! This is an independent graph/region witness, NEVER a Boolean commit.
use rustsolid::{
    BooleanOperation, BrepBody, CylindricalBrep, FaceId,
    GeometryTolerance, HoleParameters, MultiFaceRegionHold, MultiFaceRegionPlan,
    MultiFaceRegionState, Point3, Solid, Surface3, ThroughHoleBrep,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn block(origin:[f64;3],dims:[f64;3])->BrepBody{
    BrepBody::from(Solid::block(
        p(origin[0],origin[1],origin[2]),dims[0],dims[1],dims[2]).unwrap())
}
fn cylinder()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap())
}

#[test]
fn disjoint_source_bodies_have_no_false_crossing_pairs(){
    let left=block([8.,8.,8.],[2.,2.,2.]);
    let right=cylinder();
    let plan=MultiFaceRegionPlan::plan(&left,&right,
        BooleanOperation::Intersection,tol()).unwrap();
    assert_eq!(plan.state,MultiFaceRegionState::Empty,"{plan:?}");
    assert!(plan.pairs.is_empty());
    assert_eq!(plan.total_analytic_arcs,0);
    assert!(!plan.can_apply);
    plan.validate_against(&left,&right,tol()).unwrap();
}

#[test]
fn multiple_finite_faces_produce_source_owned_local_region_witnesses(){
    let left=block([-1.5,1.0,-1.5],[3.,3.,3.]);
    let right=cylinder();
    let original=format!("{left:?}");
    let plan=MultiFaceRegionPlan::plan(
        &left,&right,BooleanOperation::Difference,tol()).unwrap();
    assert_eq!(plan.state,MultiFaceRegionState::ClassifiedPreview,"{plan:?}");
    assert!(plan.pairs.len()>=4,"{plan:?}");
    assert!(plan.total_analytic_arcs>=4);
    assert_eq!(plan.total_side_samples,plan.total_analytic_arcs*6);
    assert!(!plan.can_apply);
    for pair in &plan.pairs{
        assert_eq!(pair.graph.source_faces[0].face,pair.left_face);
        assert_eq!(pair.graph.source_faces[1].face,pair.right_face);
        assert!(!pair.graph.arcs.is_empty());
        assert!(pair.regions.samples.iter().all(|s|s.keep_positive!=s.keep_negative));
    }
    plan.validate_against(&left,&right,tol()).unwrap();
    assert_eq!(format!("{left:?}"),original,"read-only preview modified the source");

    let mut broken=plan.clone();
    broken.total_side_samples+=1;
    assert!(broken.validate_against(&left,&right,tol()).is_err());
}

#[test]
fn native_circular_inner_loops_are_no_longer_globally_rejected(){
    let source=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    let pierced=BrepBody::from(ThroughHoleBrep::from_box(
        &source,HoleParameters{center_x:0.,center_z:0.,radius:0.75}).unwrap());
    let other=block([-1.,1.,-1.],[2.,2.,2.]);
    let plan=MultiFaceRegionPlan::plan(
        &pierced,&other,BooleanOperation::Union,tol()).unwrap();
    assert_eq!(plan.state,MultiFaceRegionState::ClassifiedPreview,"{plan:?}");
    assert!(!plan.pairs.is_empty());
    assert!(plan.total_analytic_arcs>=1);
    assert_eq!(plan.total_reconciled_edges,plan.total_analytic_arcs);
    assert!(!plan.can_apply);
    plan.validate_against(&pierced,&other,tol()).unwrap();
}

#[test]
fn near_tangent_or_boundary_overlap_is_never_represented_as_classified_cut(){
    let source=block([2.0+4e-10,1.0,-3.0],[2.0,2.0,6.0]);
    let plan=MultiFaceRegionPlan::plan(
        &source,&cylinder(),BooleanOperation::Difference,tol()).unwrap();
    assert!(matches!(plan.state,MultiFaceRegionState::Held(_)),"{plan:?}");
    assert!(plan.pairs.is_empty());
    assert_eq!(plan.total_side_samples,0);
    assert!(!plan.can_apply);
}

#[test]
fn selected_annular_cap_and_crossing_plane_produce_verified_face_cells(){
    let base=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    let ring=BrepBody::from(ThroughHoleBrep::from_box(
        &base,HoleParameters{center_x:0.,center_z:0.,radius:1.}).unwrap());
    let crossing=block([0.,4.,-4.],[2.,2.,8.]);
    let mesh=crossing.shared().unwrap();
    let xface=mesh.faces.iter().find_map(|face|{
        let Surface3::Plane(plane)=mesh.face_surface(face.id)? else{return None};
        ((plane.frame.origin.x).abs()<1e-9 && plane.frame.normal.x.abs()>0.9)
            .then_some(face.id)
    }).expect("physical cutter plane");
    let selection=[(FaceId(1),xface)];
    let plan=MultiFaceRegionPlan::plan_selected(&ring,&crossing,
        BooleanOperation::Difference,tol(),&selection).unwrap();
    assert_eq!(plan.state,MultiFaceRegionState::ClassifiedPreview,"{plan:?}");
    assert_eq!(plan.pairs.len(),1);
    assert_eq!(plan.total_analytic_arcs,2);
    assert_eq!(plan.total_reconciled_edges,2);
    assert_eq!(plan.total_uv_cells,2);
    let face=plan.pairs[0].first_uv_cells.as_ref().expect("certified annular UV cells");
    assert!((face.cells[0].area+face.cells[1].area-
        (36.0-std::f64::consts::PI)).abs()<1e-7);
    assert_eq!(face.cut_twins.len(),plan.pairs[0].reconciled.edges.len());
    assert!(!plan.can_apply);
    plan.validate_against(&ring,&crossing,tol()).unwrap();
    let mut tampered=plan.clone();
    tampered.pairs[0].reconciled.edges[0].start.x+=0.1;
    assert!(tampered.validate_against(&ring,&crossing,tol()).is_err());
}
