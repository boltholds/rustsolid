//! First-party analytic ownership and watertightness checks of the v0.38
//! display. A mesh is NEVER the source of exact geometry or the original mass.
use std::f64::consts::TAU;
use rustsolid::{
    BrepModel, ClosedCylinderImprint, CylindricalBrep, FaceId, Frame3,
    GeometryTolerance, IntersectionGraphBuilder, ImprintResult,
    ModelUnits, Point3, Solid, Surface3, TopologyDeltaPlanner,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn dist(a:Point3,b:Point3)->f64{(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn close(a:f64,b:f64,eps:f64){assert!((a-b).abs()<eps,"actual {a} expected {b} tolerance {eps}");}
fn fixture(tilt:f64,plane_y:f64)->ClosedCylinderImprint {
    let (sin,cos)=tilt.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(cos,sin,0.),p(0.,0.,-1.),tol()).unwrap();
    let original=CylindricalBrep::new(frame,2.,5.,tol()).unwrap();
    let side=BrepModel::from_cylinder(&original).unwrap();
    let block=Solid::block(p(-8.,plane_y,-8.),16.,1.,16.).unwrap();
    let plane=BrepModel::from_polyhedral(&block,ModelUnits::default()).unwrap();
    let face=plane.faces.iter().find_map(|f|{
        let Surface3::Plane(s)=plane.face_surface(f.id)? else{return None};
        ((s.frame.origin.y-plane_y).abs()<1e-10 &&
            s.frame.normal.y.abs()>0.9).then_some(f.id)
    }).unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&plane,face,&side,FaceId(2)).unwrap();
    let plan=TopologyDeltaPlanner::new().plan(&graph,&plane,face,&side,FaceId(2)).unwrap();
    let mut editor=ClosedCylinderImprint::new(original).unwrap();
    assert!(matches!(editor.imprint(&graph,&plan,&plane,&side).unwrap(),
        ImprintResult::Committed(_)));
    editor
}

#[test]
fn all_four_faces_are_meshed_with_a_single_shared_elliptic_ring() {
    for n in [12usize,16,24,64,128,256] {
        let editor=fixture(0.3,2.5);
        let display=editor.tessellate(n).unwrap();
        editor.validate_display(&display).unwrap();
        let mesh=&display.mesh;
        let cert=display.certificate;
        assert_eq!(cert.angular_segments,n);
        assert_eq!(cert.source_revision,1);
        assert_eq!(mesh.vertices.len(),3*n+2);
        assert_eq!(mesh.triangles.len(),6*n);
        assert_eq!(mesh.triangle_faces.len(),6*n);
        assert_eq!(cert.mesh_euler,2);
        assert_eq!(cert.triangle_count,6*n);
        assert!(cert.max_analytic_edge_residual<1e-8);
        let mut owner_counts=[0usize;4];
        for face in &mesh.triangle_faces {owner_counts[face.0 as usize]+=1;}
        assert_eq!(owner_counts,[n,n,2*n,2*n]);
        let source=editor.model();
        let ellipse=source.edge_curve(rustsolid::EdgeId(4)).unwrap();
        for i in 0..n {
            let u=TAU*i as f64/n as f64;
            assert!(dist(mesh.vertices[n+i],ellipse.evaluate(u).unwrap())<1e-9);
        }
        assert!(dist(mesh.vertices[n],source.vertices[2].position)<1e-9);
        let expected=0.5*n as f64*4.0*(TAU/n as f64).sin()*5.0;
        close(cert.inscribed_prism_volume,expected,1e-8);
        close(cert.facet_volume,expected,1e-8);
        assert!(cert.facet_volume<source.mass.volume);
        // The analytic area of the original cylinder is not replaced by
        // its faceted area, and both split cylindrical faces are represented.
        let side_area=cert.face_area[2]+cert.face_area[3];
        close(side_area,2.0*(n as f64)*2.0*(std::f64::consts::PI/n as f64).sin()*5.0,1e-8);
        assert!(cert.face_area.iter().all(|a|a.is_finite() && *a>0.0));
    }
}

#[test]
fn rotated_cylinders_keep_face_ownership_mesh_volume_and_source_mass(){
    for tilt in [0.0f64,0.1,0.21,0.4,0.52] {
        let editor=fixture(tilt,2.5);
        let m=editor.tessellate(96).unwrap();
        editor.validate_display(&m).unwrap();
        let expected=(96.0*2.0)* (TAU/96.0).sin()*5.0;
        close(m.certificate.facet_volume,expected,1e-8);
        assert!(m.certificate.face_area[2]>0.0);
        assert!(m.certificate.face_area[3]>0.0);
        editor.model().validate().unwrap();
    }
}

#[test]
fn tampering_with_triangle_owner_winding_ring_or_certificate_is_rejected(){
    let editor=fixture(0.3,2.5);
    let actual=editor.tessellate(32).unwrap();
    let mut broken=actual.clone();
    broken.mesh.triangle_faces[4]=FaceId(0);
    assert!(editor.validate_display(&broken).is_err(),
        "a side triangle must never claim ownership of a cap");
    let mut broken=actual.clone();
    broken.mesh.triangles[8].swap(0,1);
    assert!(editor.validate_display(&broken).is_err(),
        "incorrect mesh winding must be rejected");
    let mut broken=actual.clone();
    broken.mesh.vertices[32].x+=0.01;
    assert!(editor.validate_display(&broken).is_err(),
        "new ellipse vertices must stay on the physical edge");
    let mut broken=actual.clone();
    broken.certificate.facet_volume+=5.0;
    assert!(editor.validate_display(&broken).is_err(),
        "certificate must be derived from actual triangles");
    let mut broken=actual.clone();
    broken.mesh.triangles.pop();
    broken.mesh.triangle_faces.pop();
    assert!(editor.validate_display(&broken).is_err(),
        "a display with a boundary hole must fail");
}

#[test]
fn failed_sampling_and_stale_display_never_change_source_revision(){
    let original=CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap();
    let editor=ClosedCylinderImprint::new(original).unwrap();
    assert!(editor.tessellate(24).is_err(),
        "unimprinted canonical cylinder is outside split-side mesher scope");
    let mut editor=fixture(0.3,2.5);
    let original=format!("{:?}",editor.model());
    let rev=editor.revision();
    assert!(editor.tessellate(11).is_err());
    assert!(editor.tessellate(4097).is_err());
    assert_eq!(editor.revision(),rev);
    assert_eq!(format!("{:?}",editor.model()),original);
    let mesh=editor.tessellate(64).unwrap();
    assert!(editor.undo().unwrap());
    assert!(editor.validate_display(&mesh).is_err());
    assert!(editor.redo().unwrap());
    assert!(editor.validate_display(&mesh).is_err(),
        "a cached mesh must carry the current revision, not past topology");
    let new=editor.tessellate(64).unwrap();
    editor.validate_display(&new).unwrap();
}
