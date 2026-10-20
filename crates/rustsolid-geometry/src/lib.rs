//! Analytic geometry independent of B-rep storage and topological ownership.
//! All coordinates use explicit model units; source analytic curves are not tessellated.
use rustsolid_foundation::{GeometryError, GeometryTolerance, Point3};
use std::f64::consts::{FRAC_PI_2, TAU};

fn invalid(reason: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(reason.into())
}
fn invalid_dimension(reason: &'static str) -> GeometryError {
    GeometryError::InvalidDimension(reason)
}
fn length(p: Point3) -> f64 { p.x.hypot(p.y).hypot(p.z) }
fn unit(v: Point3, epsilon: f64) -> Result<Point3, GeometryError> {
    let magnitude = length(v);
    if !magnitude.is_finite() || magnitude <= epsilon {
        return Err(invalid_dimension("geometry direction below tolerance"));
    }
    Ok(v.scale(1.0 / magnitude))
}
fn distance(a: Point3, b: Point3) -> f64 { length(a.sub(b)) }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit { Millimeter, Meter, Inch }
impl LengthUnit {
    pub fn meters_per_unit(self) -> f64 {
        match self { Self::Millimeter => 0.001, Self::Meter => 1.0, Self::Inch => 0.0254 }
    }
    pub fn convert(self, amount: f64, target: LengthUnit) -> Result<f64, GeometryError> {
        let converted = amount * (self.meters_per_unit() / target.meters_per_unit());
        if !converted.is_finite() { return Err(invalid_dimension("length unit conversion overflow")); }
        Ok(converted)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AngleUnit { Radian }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModelUnits { pub length: LengthUnit, pub angle: AngleUnit }
impl Default for ModelUnits {
    fn default() -> Self { Self { length: LengthUnit::Millimeter, angle: AngleUnit::Radian } }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParameterRange {
    Unbounded,
    Bounded { start: f64, end: f64 },
}
impl ParameterRange {
    pub fn bounded(start: f64, end: f64) -> Result<Self, GeometryError> {
        if !start.is_finite() || !end.is_finite() || start >= end {
            return Err(invalid_dimension("parameter interval must have finite ordered endpoints"));
        }
        Ok(Self::Bounded { start, end })
    }
    pub fn includes(self, parameter: f64, slack: f64) -> bool {
        if !parameter.is_finite() || !slack.is_finite() || slack < 0.0 { return false; }
        match self {
            Self::Unbounded => true,
            Self::Bounded { start, end } => parameter >= start - slack && parameter <= end + slack,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParameterDomain2 {
    pub u: ParameterRange,
    pub v: ParameterRange,
    pub u_periodic: bool,
    pub v_periodic: bool,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2Param { pub u: f64, pub v: f64 }
impl Point2Param {
    pub fn new(u: f64, v: f64) -> Result<Self, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("UV parameters must be finite")); }
        Ok(Self { u, v })
    }
}

/// Right-handed orthonormal frame. The stored X/Y directions are explicit;
/// `normal` is derived, not independently guessed from a surface's mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame3 {
    pub origin: Point3,
    pub x: Point3,
    pub y: Point3,
    pub normal: Point3,
}
impl Frame3 {
    pub fn from_axes(origin: Point3, x_axis: Point3, y_axis: Point3,
                     tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !origin.is_finite() || !x_axis.is_finite() || !y_axis.is_finite() {
            return Err(invalid_dimension("frame coordinates must be finite"));
        }
        let x = unit(x_axis, tolerance.absolute_length)?;
        let n = unit(x.cross(y_axis), tolerance.absolute_length * length(y_axis).max(1.0))?;
        let y = unit(n.cross(x), tolerance.absolute_length)?;
        if !x.is_finite() || !y.is_finite() || !n.is_finite() {
            return Err(invalid_dimension("invalid frame basis"));
        }
        Ok(Self { origin, x, y, normal: n })
    }
    pub fn position(self, u: f64, v: f64) -> Result<Point3, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("surface coordinates must be finite")); }
        let p = self.origin.add(self.x.scale(u)).add(self.y.scale(v));
        if !p.is_finite() { return Err(invalid_dimension("surface evaluation overflow")); }
        Ok(p)
    }
    pub fn project(self, point: Point3) -> Result<Point2Param, GeometryError> {
        if !point.is_finite() { return Err(invalid_dimension("projected point must be finite")); }
        let delta = point.sub(self.origin);
        Point2Param::new(delta.dot(self.x), delta.dot(self.y))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line3 { pub origin: Point3, pub direction: Point3 }
impl Line3 {
    pub fn between(a: Point3, b: Point3, tolerance: GeometryTolerance)
        -> Result<(Self, ParameterRange), GeometryError>
    {
        tolerance.validate()?;
        if !a.is_finite() || !b.is_finite() { return Err(invalid_dimension("line coordinates must be finite")); }
        let delta = b.sub(a);
        let span = length(delta);
        if !span.is_finite() || span <= tolerance.length_at(span) {
            return Err(invalid_dimension("line length below tolerance"));
        }
        Ok((Self { origin: a, direction: delta.scale(1.0 / span) }, ParameterRange::bounded(0.0, span)?))
    }
    pub fn evaluate(self, t: f64) -> Result<Point3, GeometryError> {
        if !t.is_finite() { return Err(invalid_dimension("line parameter must be finite")); }
        let p = self.origin.add(self.direction.scale(t));
        if !p.is_finite() { return Err(invalid_dimension("line evaluation overflow")); }
        Ok(p)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle3 { pub frame: Frame3, pub radius: f64 }
impl Circle3 {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("circle radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
    pub fn evaluate(self, angle: f64) -> Result<Point3, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("circle angle must be finite")); }
        let p = self.frame.origin.add(self.frame.x.scale(self.radius * angle.cos()))
            .add(self.frame.y.scale(self.radius * angle.sin()));
        if !p.is_finite() { return Err(invalid_dimension("circle evaluation overflow")); }
        Ok(p)
    }
    pub fn derivative(self, angle: f64) -> Result<Point3, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("circle angle must be finite")); }
        Ok(self.frame.x.scale(-self.radius * angle.sin())
            .add(self.frame.y.scale(self.radius * angle.cos())))
    }
}
/// Parametric spatial ellipse with the source-cylinder angular parameter.
/// The two basis vectors need not be orthogonal: oblique plane sections are
/// stored analytically, without fitting sampled 3D points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicCurve3 {
    pub center:Point3,
    pub cosine:Point3,
    pub sine:Point3,
}
impl HarmonicCurve3 {
    pub fn evaluate(self,t:f64)->Result<Point3,GeometryError>{
        if !t.is_finite(){return Err(invalid_dimension("nonfinite spatial ellipse parameter"));}
        let (s,c)=t.sin_cos();
        let p=self.center.add(self.cosine.scale(c)).add(self.sine.scale(s));
        if !p.is_finite(){return Err(invalid_dimension("spatial ellipse overflow"));}
        Ok(p)
    }
    pub fn derivative(self,t:f64)->Result<Point3,GeometryError>{
        if !t.is_finite(){return Err(invalid_dimension("nonfinite spatial ellipse parameter"));}
        let (s,c)=t.sin_cos();
        let p=self.sine.scale(c).sub(self.cosine.scale(s));
        if !p.is_finite(){return Err(invalid_dimension("spatial ellipse tangent overflow"));}
        Ok(p)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve3 { Line(Line3), Circle(Circle3), Harmonic(HarmonicCurve3) }
impl Curve3 {
    pub fn evaluate(self, t: f64) -> Result<Point3, GeometryError> {
        match self { Self::Line(line) => line.evaluate(t), Self::Circle(circle) => circle.evaluate(t),
            Self::Harmonic(curve)=>curve.evaluate(t) }
    }
    pub fn derivative(self, t: f64) -> Result<Point3, GeometryError> {
        match self {
            Self::Line(line) => Ok(line.direction),
            Self::Circle(circle) => circle.derivative(t),
            Self::Harmonic(curve)=>curve.derivative(t),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line2 { pub origin: Point2Param, pub delta: Point2Param }
impl Line2 {
    pub fn evaluate(self, t: f64) -> Result<Point2Param, GeometryError> {
        Point2Param::new(self.origin.u + t * self.delta.u, self.origin.v + t * self.delta.v)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2 { pub center: Point2Param, pub radius: f64 }
impl Circle2 {
    pub fn new(center: Point2Param, radius: f64, tol: GeometryTolerance) -> Result<Self, GeometryError> {
        tol.validate()?;
        if !radius.is_finite() || radius <= tol.absolute_length {
            return Err(invalid_dimension("parametric circle radius below tolerance"));
        }
        Ok(Self {center, radius})
    }
    pub fn evaluate(self, angle: f64) -> Result<Point2Param, GeometryError> {
        if !angle.is_finite() { return Err(invalid_dimension("parametric circle angle must be finite")); }
        Point2Param::new(self.center.u + self.radius * angle.cos(), self.center.v + self.radius * angle.sin())
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedCircle2 {
    pub center: Point2Param,
    pub radius: f64,
    /// UV angle is phase + direction * edge_parameter.
    pub phase: f64,
    /// Exactly +1 or -1; mirrors the source 3D angular parameter if needed.
    pub direction: f64,
}
impl OrientedCircle2 {
    pub fn new(center: Point2Param, radius: f64, phase: f64,
               direction: f64, tol: GeometryTolerance) -> Result<Self, GeometryError> {
        tol.validate()?;
        if !center.u.is_finite() || !center.v.is_finite() ||
            !radius.is_finite() || radius <= tol.absolute_length ||
            !phase.is_finite() || (direction != 1.0 && direction != -1.0) {
            return Err(invalid_dimension("invalid oriented planar circular trim"));
        }
        Ok(Self { center, radius, phase, direction })
    }
    pub fn evaluate(self, t: f64) -> Result<Point2Param, GeometryError> {
        if !t.is_finite() { return Err(invalid_dimension("nonfinite circular trim parameter")); }
        let angle = self.phase + self.direction * t;
        Point2Param::new(
            self.center.u + self.radius * angle.cos(),
            self.center.v + self.radius * angle.sin()
        )
    }
}

/// Exact trigonometric image of an oblique circular-cylinder intersection
/// on a planar UV chart, parameterized by the original source-cylinder angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicCurve2 {
    pub center:Point2Param,
    pub cosine:Point2Param,
    pub sine:Point2Param,
}
impl HarmonicCurve2 {
    pub fn evaluate(self,t:f64)->Result<Point2Param,GeometryError>{
        if !t.is_finite(){return Err(invalid_dimension("nonfinite harmonic UV parameter"));}
        let (s,c)=t.sin_cos();
        Point2Param::new(self.center.u+self.cosine.u*c+self.sine.u*s,
            self.center.v+self.cosine.v*c+self.sine.v*s)
    }
}

/// Cylindrical UV trace of an oblique section. Keeps unwrapped U and allows
/// reversed coedge use without changing the underlying physical ellipse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSection2 {
    pub angle_offset:f64,
    pub direction:f64,
    pub axial_center:f64,
    pub axial_cosine:f64,
    pub axial_sine:f64,
}
impl CylinderSection2 {
    pub fn evaluate(self,t:f64)->Result<Point2Param,GeometryError>{
        if !t.is_finite() || !self.angle_offset.is_finite() ||
            (self.direction!=1.0 && self.direction!= -1.0) ||
            ![self.axial_center,self.axial_cosine,self.axial_sine].iter().all(|v|v.is_finite()) {
            return Err(invalid_dimension("invalid oblique-cylinder UV carrier"));
        }
        let angle=self.angle_offset+self.direction*t;
        let (s,c)=angle.sin_cos();
        Point2Param::new(angle,self.axial_center+self.axial_cosine*c+self.axial_sine*s)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve2 {
    Line(Line2), Circle(Circle2), OrientedCircle(OrientedCircle2),
    Harmonic(HarmonicCurve2), CylinderSection(CylinderSection2),
}
impl Curve2 {
    pub fn evaluate(self, t: f64) -> Result<Point2Param, GeometryError> {
        match self {
            Self::Line(l) => l.evaluate(t), Self::Circle(c) => c.evaluate(t),
            Self::OrientedCircle(c) => c.evaluate(t),
            Self::Harmonic(c) => c.evaluate(t),
            Self::CylinderSection(c) => c.evaluate(t),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneSurface { pub frame: Frame3 }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSurface { pub frame: Frame3, pub radius: f64 }
impl CylinderSurface {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("cylinder radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereSurface { pub frame: Frame3, pub radius: f64 }
impl SphereSurface {
    pub fn new(frame: Frame3, radius: f64, tolerance: GeometryTolerance) -> Result<Self, GeometryError> {
        tolerance.validate()?;
        if !radius.is_finite() || radius <= tolerance.absolute_length {
            return Err(invalid_dimension("sphere radius below tolerance"));
        }
        Ok(Self { frame, radius })
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface3 {
    Plane(PlaneSurface),
    Cylinder(CylinderSurface),
    Sphere(SphereSurface),
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceDerivatives {
    pub point: Point3, pub du: Point3, pub dv: Point3,
}
impl Surface3 {
    pub fn domain(self) -> ParameterDomain2 {
        let periodic = ParameterRange::Bounded { start: 0.0, end: TAU };
        match self {
            Self::Plane(_) => ParameterDomain2 { u: ParameterRange::Unbounded, v: ParameterRange::Unbounded, u_periodic:false, v_periodic:false },
            Self::Cylinder(_) => ParameterDomain2 { u: periodic, v: ParameterRange::Unbounded, u_periodic:true, v_periodic:false },
            Self::Sphere(_) => ParameterDomain2 { u: periodic, v: ParameterRange::Bounded { start:-FRAC_PI_2, end:FRAC_PI_2 }, u_periodic:true, v_periodic:false },
        }
    }
    pub fn evaluate(self, u: f64, v: f64) -> Result<SurfaceDerivatives, GeometryError> {
        if !u.is_finite() || !v.is_finite() { return Err(invalid_dimension("surface parameters must be finite")); }
        let result = match self {
            Self::Plane(surface) => SurfaceDerivatives {
                point: surface.frame.position(u, v)?, du: surface.frame.x, dv: surface.frame.y,
            },
            Self::Cylinder(surface) => {
                let (su, cu) = u.sin_cos();
                let radial = surface.frame.x.scale(cu).add(surface.frame.y.scale(su));
                SurfaceDerivatives {
                    point: surface.frame.origin.add(radial.scale(surface.radius)).add(surface.frame.normal.scale(v)),
                    du: surface.frame.x.scale(-surface.radius * su).add(surface.frame.y.scale(surface.radius * cu)),
                    dv: surface.frame.normal,
                }
            }
            Self::Sphere(surface) => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let radial = surface.frame.x.scale(cu).add(surface.frame.y.scale(su));
                let tangent = surface.frame.x.scale(-su).add(surface.frame.y.scale(cu));
                SurfaceDerivatives {
                    point: surface.frame.origin.add(radial.scale(surface.radius * cv))
                        .add(surface.frame.normal.scale(surface.radius * sv)),
                    du: tangent.scale(surface.radius * cv),
                    dv: radial.scale(-surface.radius * sv)
                        .add(surface.frame.normal.scale(surface.radius * cv)),
                }
            }
        };
        if !result.point.is_finite() || !result.du.is_finite() || !result.dv.is_finite() {
            return Err(invalid_dimension("surface evaluation overflow"));
        }
        Ok(result)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve3Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Curve2Id(pub u32);
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SurfaceId(pub u32);

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn projected_frame_roundtrips() {
        let f=Frame3::from_axes(
            Point3{x:1.0,y:2.0,z:-4.0},
            Point3{x:1.0,y:0.0,z:0.0},
            Point3{x:0.0,y:1.0,z:0.0},
            GeometryTolerance::default(),
        ).unwrap();
        let uv=f.project(Point3{x:4.0,y:6.0,z:-4.0}).unwrap();
        let back=f.position(uv.u,uv.v).unwrap();
        assert!((back.x-4.0).abs()<1e-12);
        assert!((back.y-6.0).abs()<1e-12);
        assert!((back.z+4.0).abs()<1e-12);
    }
    #[test]
    fn parameter_ranges_are_well_formed() {
        assert!(ParameterRange::bounded(1.0,1.0).is_err());
        assert!(ParameterRange::bounded(0.0,2.0).unwrap().includes(1.0,0.0));
        assert!(!ParameterRange::bounded(0.0,2.0).unwrap().includes(3.0,0.0));
    }
}
