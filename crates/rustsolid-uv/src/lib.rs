//! Polyline-only UV operations backed by geo.
//!
//! This crate deliberately rejects non-linear analytic input at its boundary:
//! curves, source IDs, and coedge provenance belong to the higher-level
//! arrangement kernel. A polygonal Boolean output is NOT a certified B-rep.
//! It is appropriate for polygonal patches, preview, independent checks and
//! for algorithms where a separate accuracy certificate is available.

use geo::{BooleanOps, Coord, LineString, Polygon, Validation};
pub use rustsolid_geometry::Point2Param;


/// Rings are open (first point is not duplicated at the end).
#[derive(Debug, Clone, PartialEq)]
pub struct UvPolygon {
    pub outer: Vec<Point2Param>,
    pub holes: Vec<Vec<Point2Param>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvBoolean {
    Union,
    Difference,
    Intersection,
    Xor,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UvError {
    InvalidRing(&'static str),
    InvalidPolygon(String),
}

impl std::fmt::Display for UvError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRing(msg) => write!(f, "invalid UV ring: {msg}"),
            Self::InvalidPolygon(msg) => write!(f, "invalid UV polygon: {msg}"),
        }
    }
}
impl std::error::Error for UvError {}

fn ring(points: &[Point2Param]) -> Result<LineString<f64>, UvError> {
    if points.len() < 3 { return Err(UvError::InvalidRing("fewer than three vertices")); }
    if points.iter().any(|p| !p.u.is_finite() || !p.v.is_finite()) {
        return Err(UvError::InvalidRing("non-finite coordinates"));
    }
    if points.windows(2).any(|w| w[0] == w[1]) {
        return Err(UvError::InvalidRing("consecutive duplicate points"));
    }
    if points[0] == *points.last().unwrap() {
        return Err(UvError::InvalidRing("supply an open ring; closure is automatic"));
    }
    let mut coords: Vec<Coord<f64>> = points.iter().map(|p| Coord { x:p.u, y:p.v }).collect();
    coords.push(coords[0]);
    Ok(LineString::new(coords))
}

fn polygon(p: &UvPolygon) -> Result<Polygon<f64>, UvError> {
    let shape = Polygon::new(ring(&p.outer)?,
        p.holes.iter().map(|points|ring(points)).collect::<Result<Vec<_>,_>>()?);
    shape.check_validation().map_err(|e| UvError::InvalidPolygon(e.to_string()))?;
    Ok(shape)
}

fn open_ring(ls: &LineString<f64>) -> Vec<Point2Param> {
    let mut result: Vec<_> = ls.0.iter().map(|p| Point2Param{u:p.x,v:p.y}).collect();
    if result.len() > 1 && result[0] == *result.last().unwrap() { result.pop(); }
    result
}

/// Executes a *polygonal* UV Boolean; the output carries no analytic
/// curves, coedge ownership, or rigorous error certificate.
pub fn polygon_boolean(
    left: &UvPolygon, right: &UvPolygon, op: UvBoolean,
) -> Result<Vec<UvPolygon>, UvError> {
    let a = polygon(left)?;
    let b = polygon(right)?;
    let result = match op {
        UvBoolean::Union => a.union(&b),
        UvBoolean::Difference => a.difference(&b),
        UvBoolean::Intersection => a.intersection(&b),
        UvBoolean::Xor => a.xor(&b),
    };
    for p in &result.0 {
        p.check_validation().map_err(|e| UvError::InvalidPolygon(e.to_string()))?;
    }
    Ok(result.0.iter().map(|p| UvPolygon {
        outer: open_ring(p.exterior()),
        holes: p.interiors().iter().map(open_ring).collect(),
    }).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::Area;
    fn rectangle(a:f64,b:f64,c:f64,d:f64) -> UvPolygon {
        UvPolygon{outer:vec![
            Point2Param{u:a,v:b},Point2Param{u:c,v:b},
            Point2Param{u:c,v:d},Point2Param{u:a,v:d}
        ],holes:vec![]}
    }
    #[test]
    fn difference_produces_hole() {
        let a=rectangle(0.,0.,10.,10.);
        let b=rectangle(2.,2.,8.,8.);
        let out=polygon_boolean(&a,&b,UvBoolean::Difference).unwrap();
        assert_eq!(out.len(),1);
        assert_eq!(out[0].holes.len(),1);
        let area=polygon(&out[0]).unwrap().unsigned_area();
        assert!((area-64.0).abs()<1e-8);
    }
    #[test]
    fn intersection_and_disjoint() {
        let a=rectangle(0.,0.,2.,2.);
        let b=rectangle(1.,1.,3.,3.);
        let out=polygon_boolean(&a,&b,UvBoolean::Intersection).unwrap();
        assert_eq!(out.len(),1);
        assert!((polygon(&out[0]).unwrap().unsigned_area()-1.0).abs()<1e-8);
        assert!(polygon_boolean(&a,&rectangle(4.,4.,5.,5.),UvBoolean::Intersection).unwrap().is_empty());
    }
    #[test]
    fn reject_self_intersecting_polygon() {
        let p=UvPolygon{outer:vec![
            Point2Param{u:0.,v:0.},Point2Param{u:2.,v:2.},
            Point2Param{u:0.,v:2.},Point2Param{u:2.,v:0.},
        ],holes:vec![]};
        assert!(polygon_boolean(&p,&rectangle(0.,0.,1.,1.),UvBoolean::Union).is_err());
    }
    #[test]
    fn reject_invalid_polygon() {
        let p=UvPolygon{outer:vec![Point2Param{u:0.,v:0.},Point2Param{u:1.,v:f64::NAN},Point2Param{u:2.,v:0.}],holes:vec![]};
        assert!(polygon_boolean(&p,&rectangle(0.,0.,1.,1.),UvBoolean::Union).is_err());
    }
}
