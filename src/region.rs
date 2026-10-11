//! Independent solid-region predicates for planar polyhedra and exact
//! analytic right circular cylinders. Surface points are an explicit state:
//! callers cannot silently turn a tangent/contact into a cut.
use crate::{BrepBody, GeometryError, GeometryTolerance, Point3, Solid, CylindricalBrep};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointRegion {
    Interior,
    Boundary,
    Exterior,
}

/// Both model-space and floating-point coordinates are honored, while the
/// geometric tolerance remains based on the local feature extent.
#[derive(Debug, Clone, Copy)]
pub struct RegionClassifier {
    pub tolerance: GeometryTolerance,
}
impl RegionClassifier {
    pub fn new(tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        Ok(Self {tolerance})
    }
    pub fn classify(&self, body:&BrepBody, point:Point3)
        -> Result<PointRegion,GeometryError>
    {
        match body {
            BrepBody::Polyhedral(solid) => self.classify_solid(solid,point),
            BrepBody::Cylinder(cylinder) => self.classify_cylinder(cylinder,point),
        }
    }
    pub fn classify_cylinder(&self, cylinder:&CylindricalBrep, point:Point3)
        -> Result<PointRegion,GeometryError>
    {
        cylinder.validate()?;
        if !point.is_finite() {
            return Err(GeometryError::InvalidEdit("region query requires finite coordinates".into()));
        }
        let delta=point.sub(cylinder.frame.origin);
        let axial=delta.dot(cylinder.frame.normal);
        let radial=delta.dot(cylinder.frame.x).hypot(delta.dot(cylinder.frame.y));
        if !radial.is_finite() || !axial.is_finite() {
            return Err(GeometryError::InvalidEdit("region query coordinate overflow".into()));
        }
        let eps=self.tolerance.length_at(cylinder.radius.max(cylinder.height))
            .max(cylinder.tolerance.length_at(cylinder.radius.max(cylinder.height)));
        if radial>cylinder.radius+eps || axial< -eps || axial>cylinder.height+eps {
            return Ok(PointRegion::Exterior);
        }
        if (radial-cylinder.radius).abs()<=eps || axial.abs()<=eps ||
            (axial-cylinder.height).abs()<=eps {
            return Ok(PointRegion::Boundary);
        }
        Ok(PointRegion::Interior)
    }
    pub fn classify_solid(&self, solid:&Solid, point:Point3)
        -> Result<PointRegion,GeometryError>
    {
        solid.validate()?;
        if !point.is_finite() {
            return Err(GeometryError::InvalidEdit("region query requires finite coordinates".into()));
        }
        let extent=solid.bbox.max.sub(solid.bbox.min).length_squared().sqrt();
        let eps=self.tolerance.length_at(extent).max(solid.tolerance.length_at(extent));
        if point.x<solid.bbox.min.x-eps || point.y<solid.bbox.min.y-eps ||
            point.z<solid.bbox.min.z-eps || point.x>solid.bbox.max.x+eps ||
            point.y>solid.bbox.max.y+eps || point.z>solid.bbox.max.z+eps {
            return Ok(PointRegion::Exterior);
        }
        // Compute the oriented solid-angle sum of the *validated closed mesh*.
        // This works for concave prisms as well as convex bodies, independent
        // of ray directions, triangulation diagonals and polygonal face splits.
        // All vector arithmetic is point-relative for translation stability.
        let mut sum=0.0_f64;
        let mut correction=0.0_f64;
        for tri in &solid.mesh.triangles {
            let a=solid.mesh.vertices[tri[0] as usize].sub(point);
            let b=solid.mesh.vertices[tri[1] as usize].sub(point);
            let c=solid.mesh.vertices[tri[2] as usize].sub(point);
            let la=a.length_squared().sqrt();
            let lb=b.length_squared().sqrt();
            let lc=c.length_squared().sqrt();
            if !la.is_finite() || !lb.is_finite() || !lc.is_finite() {
                return Err(GeometryError::InvalidEdit("region solid-angle overflow".into()));
            }
            if la<=eps || lb<=eps || lc<=eps {return Ok(PointRegion::Boundary);}
            let n=b.sub(a).cross(c.sub(a));
            let nlen=n.length_squared().sqrt();
            if nlen<=0.0 || !nlen.is_finite() {
                return Err(GeometryError::InvalidTopology("degenerate region triangle".into()));
            }
            let plane_dist=a.dot(n).abs()/nlen;
            if plane_dist<=eps {
                // Project the query onto the plane and test barycentric
                // inclusion. Exact boundary classification includes edges,
                // interior face points and sites near triangulation seams.
                let ab=b.sub(a);let ac=c.sub(a);let ap=a.scale(-1.0);
                let d00=ab.dot(ab);let d01=ab.dot(ac);let d11=ac.dot(ac);
                let d20=ap.dot(ab);let d21=ap.dot(ac);
                let den=d00*d11-d01*d01;
                if den<=0.0 || !den.is_finite(){
                    return Err(GeometryError::InvalidTopology("invalid region triangle area".into()));
                }
                let u=(d11*d20-d01*d21)/den;
                let v=(d00*d21-d01*d20)/den;
                let slack=eps/(d00.max(d11).sqrt());
                if u>=-slack && v>=-slack && u+v<=1.0+slack {
                    return Ok(PointRegion::Boundary);
                }
            }
            let numerator=a.dot(b.cross(c));
            let denominator=la*lb*lc+a.dot(b)*lc+b.dot(c)*la+c.dot(a)*lb;
            if !numerator.is_finite() || !denominator.is_finite() {
                return Err(GeometryError::InvalidEdit("region solid-angle accumulator overflow".into()));
            }
            let angle=2.0*numerator.atan2(denominator);
            let kahan=angle-correction;
            let next=sum+kahan;
            correction=(next-sum)-kahan;
            sum=next;
        }
        if sum.abs()>3.0*PI {Ok(PointRegion::Interior)}
        else if sum.abs()<PI {Ok(PointRegion::Exterior)}
        else {
            Err(GeometryError::InvalidEdit(
                "region solid-angle winding is numerically ambiguous".into()))
        }
    }
}

/// A Boolean truth table is defined only for proper 3D interior/exterior
/// samples. Contact points are kept as Boundary and must be classified from
/// neighborhoods before editing source topology.
pub fn boolean_region(
    a:PointRegion,b:PointRegion,op:crate::BooleanOperation,
) -> Option<bool> {
    use PointRegion::{Interior,Exterior,Boundary};
    if a==Boundary || b==Boundary {return None;}
    let a=a==Interior;
    let b=b==Interior;
    Some(match op {
        crate::BooleanOperation::Union=>a||b,
        crate::BooleanOperation::Intersection=>a&&b,
        crate::BooleanOperation::Difference=>a&&!b,
    })
}
