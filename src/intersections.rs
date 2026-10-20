//! Public analytic supporting-surface intersection kernel lives in rustsolid-intersect.
//! Topology-aware trimming stays in the modeling/Boolean integration layer.
pub use rustsolid_intersect::{
    SurfaceIntersectionEngine, SurfaceIntersection, IntersectionCurve,
    PlaneCylinderSection, IntersectionContact, IntersectionAmbiguity,
    IntersectionCertificate, IntersectionPair, ALL_INTERSECTION_PAIRS,
};
