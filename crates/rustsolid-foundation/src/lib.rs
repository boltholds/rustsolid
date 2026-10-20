//! Kernel fundamentals, shared without dependency on topology or modeling operations.
//! Length coordinates retain the existing coordinate axes; the API may evolve.
use std::fmt;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeometryError {
    InvalidProfile(&'static str),
    InvalidDimension(&'static str),
    InvalidTopology(String),
    InvalidTolerance(&'static str),
    InvalidEdit(String),
    TriangulationFailed,
}

impl fmt::Display for GeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidProfile(reason) => write!(f, "invalid extrusion profile: {reason}"),
            Self::InvalidDimension(reason) => write!(f, "invalid dimension: {reason}"),
            Self::InvalidTopology(reason) => write!(f, "invalid topology: {reason}"),
            Self::InvalidTolerance(reason) => write!(f, "invalid tolerance: {reason}"),
            Self::InvalidEdit(reason) => write!(f, "invalid edit: {reason}"),
            Self::TriangulationFailed => write!(f, "polygon triangulation failed"),
        }
    }
}

impl std::error::Error for GeometryError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    pub x: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub fn array(self) -> [f64; 3] { [self.x, self.y, self.z] }
    pub fn is_finite(self) -> bool { self.x.is_finite() && self.y.is_finite() && self.z.is_finite() }
    pub fn add(self, v: Point3) -> Self {
        Self { x: self.x + v.x, y: self.y + v.y, z: self.z + v.z }
    }
    pub fn sub(self, v: Point3) -> Self {
        Self { x: self.x - v.x, y: self.y - v.y, z: self.z - v.z }
    }
    pub fn scale(self, s: f64) -> Self {
        Self { x: self.x * s, y: self.y * s, z: self.z * s }
    }
    pub fn dot(self, v: Point3) -> f64 { self.x * v.x + self.y * v.y + self.z * v.z }
    pub fn cross(self, v: Point3) -> Self {
        Self {
            x: self.y * v.z - self.z * v.y,
            y: self.z * v.x - self.x * v.z,
            z: self.x * v.y - self.y * v.x,
        }
    }
    pub fn length_squared(self) -> f64 { self.dot(self) }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GeometryTolerance {
    /// Minimum positional resolution, in model length units.
    pub absolute_length: f64,
    /// Positional resolution relative to the measured feature extent.
    pub relative_length: f64,
    /// Maximum allowed normal deviation in radians.
    pub angular: f64,
}

impl Default for GeometryTolerance {
    fn default() -> Self {
        Self { absolute_length: 1.0e-9, relative_length: 1.0e-12, angular: 1.0e-8 }
    }
}

impl GeometryTolerance {
    pub fn validate(self) -> Result<(), GeometryError> {
        if !self.absolute_length.is_finite() || self.absolute_length <= 0.0 {
            return Err(GeometryError::InvalidTolerance("absolute_length must be finite and positive"));
        }
        if !self.relative_length.is_finite() || !(0.0..1.0).contains(&self.relative_length) {
            return Err(GeometryError::InvalidTolerance("relative_length must be finite, nonnegative and < 1"));
        }
        if !self.angular.is_finite() || self.angular <= 0.0 || self.angular >= std::f64::consts::FRAC_PI_2 {
            return Err(GeometryError::InvalidTolerance("angular must be finite and between 0 and pi/2"));
        }
        Ok(())
    }

    pub fn length_at(self, feature_extent: f64) -> f64 {
        self.absolute_length.max(self.relative_length * feature_extent)
    }

    pub fn area_at(self, feature_extent: f64) -> f64 {
        self.length_at(feature_extent) * feature_extent
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn point_and_tol() {
        let a = Point3 { x: 2.0, y: 0.0, z: -1.0 };
        let b = Point3 { x: -2.0, y: 1.0, z: 1.0 };
        assert_eq!(a.add(b), Point3{x:0.0,y:1.0,z:0.0});
        assert!(GeometryTolerance::default().validate().is_ok());
        assert!(GeometryTolerance{absolute_length:-1.0,..Default::default()}.validate().is_err());
    }
}
