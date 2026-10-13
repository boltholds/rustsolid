//! Source-certified read-only geometry for an upright through cylinder and
//! a finite-depth top-entry blind cylinder in an axis-aligned box.
//! No B-rep mutation. The mixed-body executor still rejects touching bores.

use std::f64::consts::{PI, TAU};
use crate::{BlindHoleParameters, BoundingBox, GeometryError, GeometryTolerance,
    HoleParameters, Point3, Solid};

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidEdit(message.into())
}

/// Mathematical relationship between two disk footprints in world XZ.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MixedFootprintContact {
    Disjoint,
    TwoTransverseCrossings,
    BlindInsideThrough,
    ThroughInsideBlind,
    Coincident,
    HeldNearTangency,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FootprintSource { Through, Blind }

/// Source-bound CCW arc of the exposed two-disk union contour.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FootprintExposedArc {
    pub source: FootprintSource,
    pub start_crossing: usize,
    pub end_crossing: usize,
    pub start_angle: f64,
    pub sweep_radians: f64,
    pub midpoint: Point3,
}

/// Source-bound mathematical intent. Public fields may be inspected or
/// modified, but validate_against rejects any tampering and can_apply is
/// unconditionally false. No source or topology mutation is possible.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedFootprintOverlapPlan {
    pub contact: MixedFootprintContact,
    pub overlap_area: Option<f64>,
    pub predicted_remaining_volume: Option<f64>,
    pub crossing_points: Option<[Point3; 2]>,
    pub exposed_arcs: Option<[FootprintExposedArc; 2]>,
    pub length_budget: f64,
    pub source_bounds: BoundingBox,
    pub source_volume: f64,
    pub source_revision: u64,
    pub through: HoleParameters,
    pub blind: BlindHoleParameters,
    pub tolerance: GeometryTolerance,
}
impl MixedFootprintOverlapPlan {
    pub fn can_apply(&self) -> bool { false }
    pub fn validate_against(&self, box_: &Solid, through: HoleParameters,
                            blind: BlindHoleParameters) -> Result<(), GeometryError> {
        let expected = MixedFootprintOverlapPlanner::plan(box_, through, blind)?;
        if *self != expected { return Err(invalid("overlap plan/source witness mismatch")); }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct MixedFootprintOverlapPlanner;

fn bounded(bounds: BoundingBox, x: f64, z: f64, radius: f64, eps: f64) -> bool {
    x-radius > bounds.min.x+eps && x+radius < bounds.max.x-eps &&
    z-radius > bounds.min.z+eps && z+radius < bounds.max.z-eps
}
fn point_on_circle(x: f64, z: f64, radius: f64, y: f64, angle: f64) -> Point3 {
    Point3 { x: x+radius*angle.cos(), y, z: z+radius*angle.sin() }
}

/// Select one CCW circle arc with midpoint outside the other disk, and
/// prove that the opposite arc midpoint lies inside it.
fn exposed_arc(role: FootprintSource, x: f64, z: f64, radius: f64,
               ox: f64, oz: f64, other_radius: f64,
               points: [Point3;2], eps: f64) -> Option<FootprintExposedArc> {
    let angle = points.map(|p| (p.z-z).atan2(p.x-x).rem_euclid(TAU));
    let sample = |i: usize, j: usize| {
        let sweep = (angle[j]-angle[i]).rem_euclid(TAU);
        let midpoint = point_on_circle(x,z,radius,points[0].y,angle[i]+sweep*0.5);
        let separation = (midpoint.x-ox).hypot(midpoint.z-oz);
        (sweep,midpoint,separation)
    };
    let first = sample(0,1);
    let second = sample(1,0);
    let span_valid = |s:f64| s.is_finite() && s>32.0*f64::EPSILON &&
        s<TAU-32.0*f64::EPSILON;
    if span_valid(first.0) && first.2>other_radius+eps*4.0 &&
        second.2<other_radius-eps*4.0 {
        Some(FootprintExposedArc { source:role, start_crossing:0,end_crossing:1,
            start_angle:angle[0],sweep_radians:first.0,midpoint:first.1 })
    } else if span_valid(second.0) && second.2>other_radius+eps*4.0 &&
        first.2<other_radius-eps*4.0 {
        Some(FootprintExposedArc { source:role, start_crossing:1,end_crossing:0,
            start_angle:angle[1],sweep_radians:second.0,midpoint:second.1 })
    } else { None }
}
impl MixedFootprintOverlapPlanner {
    pub fn plan(box_:&Solid, through:HoleParameters, blind:BlindHoleParameters)
        -> Result<MixedFootprintOverlapPlan,GeometryError> {
        // A generic polyhedron cannot be replaced by its bounding box.
        crate::boolean::certify_box(box_)
            .map_err(|e|invalid(format!("source must be an exact axis-aligned box: {e}")))?;
        let tol = box_.tolerance;
        tol.validate()?;
        let bounds=box_.bbox;
        let width=bounds.max.x-bounds.min.x;
        let height=bounds.max.y-bounds.min.y;
        let depth=bounds.max.z-bounds.min.z;
        let extent=width.max(height).max(depth);
        if !extent.is_finite() || !height.is_finite() || height<=0.0 {
            return Err(invalid("degenerate box extent"));
        }
        let coords=[bounds.min.x,bounds.max.x,bounds.min.z,bounds.max.z,
            through.center_x,through.center_z,blind.center_x,blind.center_z];
        if !coords.into_iter().all(f64::is_finite) ||
            !through.radius.is_finite() || !blind.radius.is_finite() ||
            !blind.depth.is_finite() {return Err(invalid("nonfinite cutter geometry"));}
        // Respect the tolerance and floating-point representability at
        // large translated coordinates, without snapping physical shapes.
        let largest=coords.into_iter().map(f64::abs).fold(1.0_f64,f64::max);
        let eps=tol.length_at(extent).max(64.0*f64::EPSILON*largest);
        if !eps.is_finite() ||
            through.radius<=eps || blind.radius<=eps ||
            blind.depth<=eps || height-blind.depth<=eps ||
            !bounded(bounds,through.center_x,through.center_z,through.radius,eps) ||
            !bounded(bounds,blind.center_x,blind.center_z,blind.radius,eps) {
            return Err(invalid("cutters must lie inside box with positive radius/depth and clearance"));
        }
        let dx=blind.center_x-through.center_x;
        let dz=blind.center_z-through.center_z;
        let distance=dx.hypot(dz);
        let sum=through.radius+blind.radius;
        let diff=(through.radius-blind.radius).abs();
        if !distance.is_finite() || !sum.is_finite() {
            return Err(invalid("nonfinite cutter separation"));
        }
        let area_through=PI*through.radius*through.radius;
        let area_blind=PI*blind.radius*blind.radius;
        if !area_through.is_finite() || !area_blind.is_finite() {
            return Err(invalid("cutter area overflow"));
        }
        let mut plan=MixedFootprintOverlapPlan {
            contact:MixedFootprintContact::HeldNearTangency,
            overlap_area:None, predicted_remaining_volume:None,
            crossing_points:None, exposed_arcs:None, length_budget:eps,
            source_bounds:bounds, source_volume:box_.mass.volume,
            source_revision:box_.revision, through,blind,tolerance:tol,
        };
        let (contact,area)=if distance==0.0 && through.radius==blind.radius {
            (MixedFootprintContact::Coincident,area_through)
        } else if (distance<=eps && diff<=eps) ||
            (distance-sum).abs()<=eps || (distance-diff).abs()<=eps {
            return Ok(plan);
        } else if distance>sum+eps {
            (MixedFootprintContact::Disjoint,0.0)
        } else if distance<diff-eps {
            if through.radius>blind.radius {
                (MixedFootprintContact::BlindInsideThrough,area_blind)
            } else {
                (MixedFootprintContact::ThroughInsideBlind,area_through)
            }
        } else {
            if distance<=eps {return Ok(plan);}
            let ux=dx/distance;
            let uz=dz/distance;
            let x1=(distance*distance+through.radius*through.radius-
                blind.radius*blind.radius)/(2.0*distance);
            let x2=distance-x1;
            let square=through.radius*through.radius-x1*x1;
            if !square.is_finite()||square<=0.0 {return Ok(plan);}
            let h=square.sqrt();
            let second_square=blind.radius*blind.radius-x2*x2;
            if h<=eps*32.0 || !second_square.is_finite()||second_square<=0.0 ||
                (h-second_square.sqrt()).abs()>eps*8.0 {
                return Ok(plan);
            }
            let cx=through.center_x+x1*ux;
            let cz=through.center_z+x1*uz;
            let points=[
                Point3{x:cx-h*uz,y:bounds.max.y,z:cz+h*ux},
                Point3{x:cx+h*uz,y:bounds.max.y,z:cz-h*ux}
            ];
            for p in points {
                let d1=(p.x-through.center_x).hypot(p.z-through.center_z);
                let d2=(p.x-blind.center_x).hypot(p.z-blind.center_z);
                if !d1.is_finite()||!d2.is_finite()||
                    (d1-through.radius).abs()>eps*8.0||
                    (d2-blind.radius).abs()>eps*8.0 {return Ok(plan);}
            }
            let arc_t=exposed_arc(FootprintSource::Through,
                through.center_x,through.center_z,through.radius,
                blind.center_x,blind.center_z,blind.radius,points,eps);
            let arc_b=exposed_arc(FootprintSource::Blind,
                blind.center_x,blind.center_z,blind.radius,
                through.center_x,through.center_z,through.radius,points,eps);
            let (Some(arc_t),Some(arc_b))=(arc_t,arc_b) else {return Ok(plan)};
            if arc_t.end_crossing!=arc_b.start_crossing ||
                arc_b.end_crossing!=arc_t.start_crossing {return Ok(plan);}
            // Sector angles minus chord triangles, valid for unequal radii.
            let lens=through.radius.powi(2)*h.atan2(x1) +
                blind.radius.powi(2)*h.atan2(x2)-distance*h;
            let area_budget=eps*sum*16.0;
            if !lens.is_finite() || lens<0.0 ||
                lens>area_through.min(area_blind)+area_budget {return Ok(plan);}
            plan.crossing_points=Some(points);
            plan.exposed_arcs=Some([arc_t,arc_b]);
            (MixedFootprintContact::TwoTransverseCrossings,lens)
        };
        let outside=area_blind-area;
        let area_budget=eps*sum*16.0;
        if !outside.is_finite()||outside< -area_budget {
            plan.crossing_points=None;plan.exposed_arcs=None;
            return Ok(plan);
        }
        let removed=area_through*height+outside.max(0.0)*blind.depth;
        let remaining=box_.mass.volume-removed;
        if !removed.is_finite()||!remaining.is_finite()||
            remaining< -area_budget*height ||
            remaining>box_.mass.volume+area_budget*height {
            return Err(invalid("overlap volume violates source bounds"));
        }
        plan.contact=contact;
        plan.overlap_area=Some(area);
        plan.predicted_remaining_volume=Some(remaining.max(0.0));
        Ok(plan)
    }
}
