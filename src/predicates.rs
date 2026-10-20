//! Topology-aware predicate adapters. Numerical code lives in rustsolid-predicates.
use crate::{FaceId, GeometryError, Solid};
pub use rustsolid_predicates::{Orientation, PredicateClassification, PredicateKernel,
    orient2d, orient3d, Plane3, PlaneIntersection, SegmentPlaneIntersection,
    intersect_segment_plane};
pub(crate) use rustsolid_predicates::determinant2d;

/// Plane supporting a planar B-rep face (not its trimmed polygon domain).
pub fn face_support_plane(solid:&Solid, face: FaceId) -> Result<Plane3,GeometryError> {
    solid.validate()?;
    let ring=&solid.faces.get(face.0 as usize)
        .ok_or_else(||GeometryError::InvalidEdit(format!("face {:?} does not exist",face)))?.boundary;
    let origin=solid.vertices[ring[0].0 as usize].position;
    for i in 1..ring.len()-1 {
        let b=solid.vertices[ring[i].0 as usize].position;
        let c=solid.vertices[ring[i+1].0 as usize].position;
        if let Ok(p)=Plane3::through(origin,b,c,solid.tolerance){return Ok(p);}
    }
    Err(GeometryError::InvalidTopology("face has no noncollinear support triple".into()))
}
