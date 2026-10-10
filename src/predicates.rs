//! Robust orientation signs and tolerance-aware intersections of *analytic planes*.
//!
//! `robust` supplies adaptive/exact-sign orientation predicates for finite
//! IEEE-754 coordinates. Constructed intersection coordinates are still f64
//! approximations; these routines do not clip to trimmed faces or NURBS.
use crate::{FaceId, GeometryError, GeometryTolerance, Point3, Solid};
use robust::{Coord, Coord3D};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation { Positive, Negative, Zero }

fn sign(value: f64) -> Result<Orientation, GeometryError> {
    if !value.is_finite() {
        return Err(GeometryError::InvalidDimension("orientation calculation overflowed"));
    }
    Ok(if value > 0.0 { Orientation::Positive }
        else if value < 0.0 { Orientation::Negative } else { Orientation::Zero })
}

pub fn orient2d(a: [f64;2], b: [f64;2], c: [f64;2]) -> Result<Orientation, GeometryError> {
    if ![a,b,c].iter().flatten().all(|v| v.is_finite()) {
        return Err(GeometryError::InvalidDimension("orientation coordinates must be finite"));
    }
    sign(determinant2d(a,b,c))
}


/// A single adaptive two-dimensional determinant shared by the polygon
/// builder, B-rep loop validator, triangulator and sketch face editor.
/// Preconditions for internal callers: all coordinates are finite.
pub(crate) fn determinant2d(a:[f64;2], b:[f64;2], c:[f64;2]) -> f64 {
    robust::orient2d(Coord{x:a[0],y:a[1]}, Coord{x:b[0],y:b[1]}, Coord{x:c[0],y:c[1]})
}

/// Engineering-tolerance classification is explicitly separate from the exact
/// orientation sign of representable IEEE-754 input coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PredicateClassification { Positive, Negative, WithinTolerance }

#[derive(Debug, Clone, Copy)]
pub struct PredicateKernel {
    pub tolerance: GeometryTolerance,
}

impl PredicateKernel {
    pub fn new(tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        Ok(Self {tolerance})
    }

    pub fn orient2d(self, a:[f64;2], b:[f64;2], c:[f64;2])
        -> Result<PredicateClassification, GeometryError>
    {
        if ![a,b,c].iter().flatten().all(|v|v.is_finite()) {
            return Err(GeometryError::InvalidDimension("orientation coordinates must be finite"));
        }
        let determinant=determinant2d(a,b,c);
        if !determinant.is_finite() {
            return Err(GeometryError::InvalidDimension("orientation calculation overflowed"));
        }
        // Derive a LOCAL geometric extent to avoid dependence on a part's
        // distance from global origin; determinant sign is adaptive-exact,
        // but "near-zero" remains a user-controlled engineering decision.
        let extent = (b[0]-a[0]).hypot(b[1]-a[1])
            .max((c[0]-a[0]).hypot(c[1]-a[1]))
            .max((c[0]-b[0]).hypot(c[1]-b[1]));
        if !extent.is_finite() {
            return Err(GeometryError::InvalidDimension("orientation extent overflowed"));
        }
        if determinant.abs() <= self.tolerance.area_at(extent) {
            Ok(PredicateClassification::WithinTolerance)
        } else if determinant > 0.0 {
            Ok(PredicateClassification::Positive)
        } else {
            Ok(PredicateClassification::Negative)
        }
    }
}

pub fn orient3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Result<Orientation, GeometryError> {
    if ![a,b,c,d].iter().all(|v| v.is_finite()) {
        return Err(GeometryError::InvalidDimension("orientation coordinates must be finite"));
    }
    fn coord(p: Point3) -> Coord3D<f64> { Coord3D{x:p.x,y:p.y,z:p.z} }
    sign(robust::orient3d(coord(a),coord(b),coord(c),coord(d)))
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane3 {
    pub origin: Point3,
    pub normal: Point3, // unit-length, oriented by the construction triple
}

fn difference(a: Point3, b: Point3) -> Point3 { Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z} }
fn dot(a:Point3,b:Point3)->f64 { a.x*b.x+a.y*b.y+a.z*b.z }
fn cross(a:Point3,b:Point3)->Point3 {Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}}
fn mul(a:Point3,t:f64)->Point3 { Point3{x:a.x*t,y:a.y*t,z:a.z*t} }
fn add(a:Point3,b:Point3)->Point3 { Point3{x:a.x+b.x,y:a.y+b.y,z:a.z+b.z} }
fn norm(a:Point3)->f64 { a.x.hypot(a.y).hypot(a.z) }

impl Plane3 {
    pub fn through(a: Point3, b: Point3, c: Point3, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if ![a,b,c].iter().all(|p|p.is_finite()) {
            return Err(GeometryError::InvalidDimension("plane points must be finite"));
        }
        let u=difference(b,a); let v=difference(c,a);
        let extent=norm(u).max(norm(v)).max(norm(difference(c,b)));
        if !extent.is_finite() || extent <= tolerance.length_at(extent) {
            return Err(GeometryError::InvalidDimension("plane points are indistinguishable"));
        }
        let cross_uv=cross(u,v); let area=norm(cross_uv);
        if !area.is_finite() || area<=tolerance.area_at(extent) {
            return Err(GeometryError::InvalidDimension("plane points are collinear or below tolerance"));
        }
        Ok(Self{origin:a,normal:mul(cross_uv,1.0/area)})
    }

    pub fn signed_distance(self, p:Point3)->Result<f64,GeometryError> {
        if !p.is_finite(){return Err(GeometryError::InvalidDimension("point must be finite"));}
        let distance=dot(self.normal,difference(p,self.origin));
        if !distance.is_finite(){return Err(GeometryError::InvalidDimension("plane distance overflow"));}
        Ok(distance)
    }

    /// Computes the supporting-plane intersection. A Line point is an f64
    /// approximation and must not be used as an exact Boolean vertex witness.
    pub fn intersect(self, other:Self, tolerance: GeometryTolerance) -> Result<PlaneIntersection,GeometryError> {
        tolerance.validate()?;
        let direction=cross(self.normal,other.normal);
        let d=norm(direction);
        if !d.is_finite() {return Err(GeometryError::InvalidDimension("invalid plane normals"));}
        if d<=tolerance.angular.sin() {
            let offset=other.signed_distance(self.origin)?.abs();
            return Ok(if offset<=tolerance.absolute_length { PlaneIntersection::Coincident }
                else { PlaneIntersection::Parallel });
        }
        // Shift the origin to the first plane so even very large global
        // coordinates do not enter the cross-products of the intersection.
        let h=dot(other.normal,difference(other.origin,self.origin));
        let relative=mul(cross(direction,self.normal),h/(d*d));
        let point=add(self.origin,relative);
        if !point.is_finite(){return Err(GeometryError::InvalidDimension("plane intersection overflow"));}
        Ok(PlaneIntersection::Line { point, direction:mul(direction,1.0/d) })
    }
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub enum PlaneIntersection {
    Parallel,
    Coincident,
    Line {point:Point3,direction:Point3},
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub enum SegmentPlaneIntersection {
    Disjoint,
    Coplanar,
    Point { point:Point3, parameter:f64 },
}

pub fn intersect_segment_plane(
    start: Point3, end: Point3, plane:Plane3, tolerance:GeometryTolerance,
) -> Result<SegmentPlaneIntersection,GeometryError> {
    tolerance.validate()?;
    if !start.is_finite() || !end.is_finite() {
        return Err(GeometryError::InvalidDimension("segment coordinates must be finite"));
    }
    let delta=difference(end,start);let length=norm(delta);
    if !length.is_finite() || length<=tolerance.absolute_length {
        return Err(GeometryError::InvalidDimension("segment length below tolerance"));
    }
    let a=plane.signed_distance(start)?;let b=plane.signed_distance(end)?;
    let epsilon=tolerance.length_at(length);
    if a.abs()<=epsilon && b.abs()<=epsilon {return Ok(SegmentPlaneIntersection::Coplanar);}
    if a.abs()<=epsilon {return Ok(SegmentPlaneIntersection::Point{point:start,parameter:0.0});}
    if b.abs()<=epsilon {return Ok(SegmentPlaneIntersection::Point{point:end,parameter:1.0});}
    if (a>0.0)==(b>0.0) {return Ok(SegmentPlaneIntersection::Disjoint);}
    let t=a/(a-b);
    if !t.is_finite() || !(0.0..=1.0).contains(&t) {
        return Err(GeometryError::InvalidDimension("segment intersection not representable"));
    }
    let point=add(start,mul(delta,t));
    if !point.is_finite(){return Err(GeometryError::InvalidDimension("segment intersection overflow"));}
    Ok(SegmentPlaneIntersection::Point { point, parameter:t })
}

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
