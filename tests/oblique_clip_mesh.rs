//! Certified display proof for the *volume-changing* oblique halfspace clip.
//! Independent facet volume uses a regular-polygon cross-section oracle.
use std::f64::consts::{PI,TAU};
use rustsolid::{
    BrepModel,ClosedCylinderImprint,Curve3,CylindricalBrep,FaceId,Frame3,
    GeometryTolerance,GraphState,ImprintResult,
    IntersectionGraphBuilder,ModelUnits,ObliqueHalfspaceEditor,
    ObliqueHalfspaceSide,ObliqueClipOutcome,Point3,Solid,Surface3,
    TopologyDeltaPlanner,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn close(a:f64,b:f64,eps:f64){assert!((a-b).abs()<eps,"{a} versus {b}");}
fn fixture(angle:f64,height:f64)->ClosedCylinderImprint{
    let (s,c)=angle.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(c,s,0.),p(0.,0.,-1.),tol()).unwrap();
    let original=CylindricalBrep::new(frame,2.,5.,tol()).unwrap();
    let side=BrepModel::from_cylinder(&original).unwrap();
    let plane=BrepModel::from_polyhedral(&
        Solid::block(p(-9.,height,-9.),18.,1.,18.).unwrap(),ModelUnits::default()).unwrap();
    let fa=plane.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=plane.face_surface(f.id)? else{return None};
        ((s.frame.origin.y-height).abs()<1e-9 &&
         s.frame.normal.y.abs()>0.9).then_some(f.id)
    }).unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&plane,fa,&side,FaceId(2)).unwrap();
    assert_eq!(graph.state,GraphState::Ready);
    let plan=TopologyDeltaPlanner::new().plan(&graph,&plane,fa,&side,FaceId(2)).unwrap();
    let mut imp=ClosedCylinderImprint::new(original).unwrap();
    assert!(matches!(imp.imprint(&graph,&plan,&plane,&side).unwrap(),
        ImprintResult::Committed(_)));
    imp
}

#[test]
fn lower_and_upper_meshes_are_closed_own_each_face_and_partition_prism_volume(){
    for (angle,height) in [(0.,2.5),(0.12,2.2),(0.32,2.5),(0.48,2.)]{
        let imp=fixture(angle,height);
        let mut below=ObliqueHalfspaceEditor::new(&imp).unwrap();
        let mut above=ObliqueHalfspaceEditor::new(&imp).unwrap();
        assert!(matches!(below.clip(ObliqueHalfspaceSide::Negative).unwrap(),ObliqueClipOutcome::Committed(_)));
        assert!(matches!(above.clip(ObliqueHalfspaceSide::Positive).unwrap(),ObliqueClipOutcome::Committed(_)));
        for n in [12usize,24,64,128] {
            let a=below.tessellate(n).unwrap();
            let b=above.tessellate(n).unwrap();
            below.validate_display(&a).unwrap();
            above.validate_display(&b).unwrap();
            for c in [&a,&b]{
                assert_eq!(c.mesh.vertices.len(),2*n+2);
                assert_eq!(c.mesh.triangles.len(),4*n);
                assert_eq!(c.certificate.mesh_euler,2);
                assert_eq!(c.certificate.triangle_faces,[n,n,2*n]);
                assert!(c.certificate.area_by_face.iter().all(|v|*v>0.0&&v.is_finite()));
                assert!(c.certificate.max_edge_residual<1e-8);
                assert!(c.certificate.facet_volume>0.0);
            }
            let polygon=0.5*n as f64*4.0*(TAU/n as f64).sin();
            close(a.certificate.facet_volume+b.certificate.facet_volume,
                polygon*5.0,1e-7);
            close(a.certificate.facet_volume,a.certificate.source_polygon_volume,1e-8);
            close(b.certificate.facet_volume,b.certificate.source_polygon_volume,1e-8);
            let Curve3::Harmonic(e)=below.model().edge_curve(rustsolid::EdgeId(1)).unwrap()
                else{panic!("analytic ellipse missing in negative cut");};
            let amplitude=(e.cosine.y*e.sine.z-e.cosine.z*e.sine.y).powi(2)+
                (e.cosine.z*e.sine.x-e.cosine.x*e.sine.z).powi(2)+
                (e.cosine.x*e.sine.y-e.cosine.y*e.sine.x).powi(2);
            let polygon_ellipse_area=0.5*n as f64*(TAU/n as f64).sin()*amplitude.sqrt();
            close(a.certificate.area_by_face[1],polygon_ellipse_area,1e-8);
            close(b.certificate.area_by_face[0],polygon_ellipse_area,1e-8);
            let original_side_faceted=
                2.0*n as f64*2.0*(std::f64::consts::PI/n as f64).sin()*5.0;
            close(a.certificate.area_by_face[2]+b.certificate.area_by_face[2],
                original_side_faceted,1e-8);
            // Strict faceted body volume must not exceed closed-form solid.
            assert!(a.certificate.facet_volume<=below.model().mass.volume*(1.+1e-8));
            assert!(b.certificate.facet_volume<=above.model().mass.volume*(1.+1e-8));
        }
        close(below.model().mass.volume+above.model().mass.volume,
            PI*4.0*5.0,1e-8);
    }
}

#[test]
fn triangle_owner_geometry_winding_and_stale_revision_are_fail_closed(){
    let imp=fixture(0.3,2.5);
    let mut e=ObliqueHalfspaceEditor::new(&imp).unwrap();
    assert!(e.tessellate(24).is_err(),
        "a four-face imprint cannot impersonate an already clipped body");
    assert!(matches!(e.clip(ObliqueHalfspaceSide::Positive).unwrap(),
        ObliqueClipOutcome::Committed(_)));
    let mesh=e.tessellate(48).unwrap();
    let mut t=mesh.clone();
    t.mesh.triangle_faces[3]=FaceId(0);
    assert!(e.validate_display(&t).is_err());
    let mut t=mesh.clone();
    t.mesh.triangles[1].swap(0,1);
    assert!(e.validate_display(&t).is_err());
    let mut t=mesh.clone();
    t.mesh.vertices[0].z+=0.01;
    assert!(e.validate_display(&t).is_err());
    let mut t=mesh.clone();
    t.certificate.facet_volume+=0.4;
    assert!(e.validate_display(&t).is_err());
    let rev=e.revision();
    assert!(e.tessellate(11).is_err());
    assert!(e.tessellate(4097).is_err());
    assert_eq!(e.revision(),rev);
    assert!(e.undo().unwrap());
    assert!(e.validate_display(&mesh).is_err());
    assert!(e.redo().unwrap());
    assert!(e.validate_display(&mesh).is_err(),
        "redo requires a fresh display certificate");
    e.validate_display(&e.tessellate(48).unwrap()).unwrap();
}
