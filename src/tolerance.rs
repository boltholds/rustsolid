//! Model-space geometric tolerances. Distances are in the caller's model units.
//! Relative tolerance is evaluated against local *feature extent*, never global
//! distance from the origin; small parts translated far away remain valid.
use crate::GeometryError;
use serde::{Deserialize, Serialize};

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
