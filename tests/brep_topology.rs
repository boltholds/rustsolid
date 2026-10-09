use rustsolid::{
    CoedgeId, GeometryError, GeometryTolerance, LoopRole, Point2, Point3, Solid,
    TopologyEntity,
};

fn cube() -> Solid {
    Solid::block(Point3 { x: 0.0, y: 0.0, z: 0.0 }, 2.0, 3.0, 4.0).unwrap()
}

#[test]
fn cube_has_consistent_coedges_loops_shell_and_twins() {
    let solid = cube();
    assert_eq!((solid.shells.len(), solid.faces.len(), solid.loops.len()), (1, 6, 6));
    assert_eq!((solid.vertices.len(), solid.edges.len(), solid.coedges.len()), (8, 12, 24));
    for edge in &solid.edges {
        let a = &solid.coedges[edge.coedges[0].0 as usize];
        let b = &solid.coedges[edge.coedges[1].0 as usize];
        assert_eq!(a.twin, b.id);
        assert_eq!(b.twin, a.id);
        assert_ne!(a.reversed, b.reversed);
        assert_ne!(a.face, b.face);
        assert_eq!(a.start_vertex(edge), b.end_vertex(edge));
        assert_eq!(a.end_vertex(edge), b.start_vertex(edge));
    }
    for face in &solid.faces {
        let loop_data = &solid.loops[face.loops[0].0 as usize];
        assert_eq!(loop_data.role, LoopRole::Outer);
        let first = loop_data.first_coedge;
        let mut current = first;
        let mut ring = Vec::new();
        loop {
            let c = &solid.coedges[current.0 as usize];
            assert_eq!(c.face, face.id);
            assert_eq!(c.loop_id, loop_data.id);
            assert_eq!(solid.coedges[c.next.0 as usize].prev, c.id);
            ring.push(c.start_vertex(&solid.edges[c.edge.0 as usize]));
            current = c.next;
            if current == first { break; }
            assert!(ring.len() < 100);
        }
        assert_eq!(ring, face.boundary);
    }
    solid.check_topology().unwrap();
    solid.validate().unwrap();
}

#[test]
fn broken_twin_has_a_typed_witness() {
    let mut solid = cube();
    solid.coedges[0].twin = CoedgeId(0);
    let issue = solid.check_topology().unwrap_err();
    assert_eq!(issue.code, "coedge.inconsistent_twin");
    assert_eq!(issue.entity, TopologyEntity::Coedge(CoedgeId(0)));
    assert!(matches!(solid.validate(), Err(GeometryError::InvalidTopology(_))));
}

#[test]
fn broken_loop_and_face_cache_are_detected() {
    let mut solid = cube();
    solid.coedges[0].next = CoedgeId(0);
    assert_eq!(solid.check_topology().unwrap_err().code, "coedge.broken_links");

    let mut solid = cube();
    solid.faces[0].boundary.reverse();
    assert_eq!(solid.check_topology().unwrap_err().code, "face.boundary_mismatch");
}

#[test]
fn shell_membership_is_checked() {
    let mut solid = cube();
    solid.shells[0].faces[1] = solid.shells[0].faces[0];
    assert_eq!(solid.check_topology().unwrap_err().code, "shell.face_ownership");
}

#[test]
fn non_planar_face_is_rejected_with_vertex_witness() {
    let mut solid = cube();
    // Move a top-cap corner away from its face plane.
    let id = solid.faces[1].boundary[0];
    solid.vertices[id.0 as usize].position.y += 0.02;
    let problem = solid.check_topology().unwrap_err();
    assert_eq!(problem.code, "face.nonplanar");
    assert!(matches!(problem.entity, TopologyEntity::Vertex(_)));
}

#[test]
fn mesh_ownership_and_geometry_cannot_silently_drift() {
    let mut solid = cube();
    solid.mesh.vertices[0].x += 0.1;
    assert!(matches!(solid.validate(), Err(GeometryError::InvalidTopology(_))));

    let mut solid = cube();
    solid.mesh.triangle_faces[0] = solid.faces[2].id;
    assert!(matches!(solid.validate(), Err(GeometryError::InvalidTopology(_))));
}

#[test]
fn invalid_tolerances_are_rejected_and_tiny_features_are_configurable() {
    let p = [Point2{x:0.0,z:0.0}, Point2{x:1.0e-8,z:0.0},
             Point2{x:1.0e-8,z:1.0}, Point2{x:0.0,z:1.0}];
    let default_result = Solid::extrude_xz(&p, 1.0);
    assert!(default_result.is_ok());
    let strict = GeometryTolerance { absolute_length: 1.0e-6, ..Default::default() };
    assert!(Solid::extrude_xz_with_tolerance(&p, 1.0, strict).is_err());
    let custom = GeometryTolerance { absolute_length: 1.0e-12, relative_length: 1.0e-13, ..Default::default() };
    let body = Solid::extrude_xz_with_tolerance(&p, 1.0, custom).unwrap();
    assert!((body.mass.volume - 1.0e-8).abs() < 1.0e-16);
    assert_eq!(body.coedges.len(), 24);

    let invalid = GeometryTolerance { absolute_length: f64::NAN, ..Default::default() };
    assert!(matches!(Solid::extrude_xz_with_tolerance(&p,1.0,invalid), Err(GeometryError::InvalidTolerance(_))));
}

#[test]
fn translated_small_profile_at_large_coordinates_keeps_topology() {
    let origin = Point3 { x: 1.0e12, y: 1.0e12, z: -1.0e12 };
    let body = Solid::block(origin, 1.0, 1.0, 1.0).unwrap();
    assert_eq!(body.mass.volume, 1.0);
    assert_eq!(body.mass.centroid.x, origin.x + 0.5);
    assert_eq!(body.mass.centroid.y, origin.y + 0.5);
    body.check_topology().unwrap();
}

#[test]
fn input_extrusion_rejects_unresolvable_thickness_and_retains_tolerance() {
    let p = [Point2 {x:0.0,z:0.0}, Point2 {x:1.0,z:0.0},
             Point2 {x:1.0,z:1.0}, Point2 {x:0.0,z:1.0}];
    assert!(Solid::extrude_xz(&p, 1.0e-12).is_err());
    let t = GeometryTolerance { absolute_length: 1.0e-14, relative_length: 1.0e-14, ..Default::default() };
    let solid = Solid::extrude_xz_with_tolerance(&p, 1.0e-12, t).unwrap();
    assert_eq!(solid.tolerance, t);
    solid.validate().unwrap();
}

#[test]
fn self_touching_face_is_rejected_after_editing_geometry() {
    let mut solid = cube();
    // Make a crossed side quad by interchanging two geometric corners; the
    // graph still references the old vertices and must fail geometric checks.
    let a = solid.faces[2].boundary[1];
    let b = solid.faces[2].boundary[2];
    let pa = solid.vertices[a.0 as usize].position;
    let pb = solid.vertices[b.0 as usize].position;
    solid.vertices[a.0 as usize].position = pb;
    solid.vertices[b.0 as usize].position = pa;
    assert!(solid.check_topology().is_err());
}
