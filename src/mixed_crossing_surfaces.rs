//! Read-only analytic decomposition of the boundary created by one transverse
//! through/blind circular footprint overlap. Not a B-rep topology commit.
//!
//! The result describes actual mathematical planar/cylindrical support regions,
//! their source ownership and expected areas. Constructing shared edges,
//! coedges, UV trim loops and a watertight shell remains a separate milestone.

use std::f64::consts::{PI,TAU};
use crate::{BlindHoleParameters,GeometryError,HoleParameters,Point3,Solid,
    MixedFootprintContact,MixedFootprintOverlapPlan,MixedFootprintOverlapPlanner};

fn invalid(why:&str)->GeometryError{
    GeometryError::InvalidEdit(why.to_string())
}

/// One source support patch, representing part of the *actual* difference
/// boundary as a mathematical region, not yet a B-rep face.
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum MixedPatchKind {
    /// Exposed upper box cap, excluding the union of both cutter disks.
    TopCapComplement,
    /// Exposed lower box cap, excluding the through disk only.
    BottomCapComplement,
    /// Four unchanged vertical box sides.
    OuterBoxSides,
    /// Through cylinder: full circle, bottom to blind floor.
    ThroughLowerWall,
    /// Through cylinder: angular subset outside blind cutter, upper segment.
    ThroughUpperWall,
    /// Blind cylinder: angular subset outside through cutter, upper segment.
    BlindUpperWall,
    /// Horizontal blind floor minus overlapping through-disk footprint.
    BlindFloorCrescent,
}

/// Analytic support-region area; this is not a tessellation or generated face.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct MixedSurfacePatch {
    pub kind:MixedPatchKind,
    pub surface_area:f64,
    pub is_cylindrical:bool,
    pub reversed_material_sense:bool,
}

#[derive(Debug,Clone,PartialEq)]
pub struct MixedCrossingSurfacePlan {
    pub footprint:MixedFootprintOverlapPlan,
    pub patches:[MixedSurfacePatch;7],
    pub predicted_surface_area:f64,
    pub predicted_remaining_volume:f64,
    pub expected_genus:usize,
    pub expected_euler_characteristic:isize,
    pub top_boundary_arc_count:usize,
    pub bottom_boundary_arc_count:usize,
    pub splitting_height_y:f64,
}

impl MixedCrossingSurfacePlan {
    /// This contract is read-only; the actual shared B-rep is not built.
    pub fn can_apply(&self)->bool{false}
    pub fn validate_against(&self,source:&Solid,t:HoleParameters,b:BlindHoleParameters)
        ->Result<(),GeometryError>
    {
        let expected=MixedCrossingSurfacePlanner::plan(source,t,b)?;
        if *self!=expected{return Err(invalid("surface plan/source witness mismatch"));}
        Ok(())
    }
}

#[derive(Debug,Clone,Copy,Default)]
pub struct MixedCrossingSurfacePlanner;

impl MixedCrossingSurfacePlanner {
    pub fn plan(source:&Solid,t:HoleParameters,b:BlindHoleParameters)
        ->Result<MixedCrossingSurfacePlan,GeometryError>
    {
        let fp=MixedFootprintOverlapPlanner::plan(source,t,b)?;
        if fp.contact!=MixedFootprintContact::TwoTransverseCrossings {
            return Err(invalid("surface planner requires exactly two certified transverse crossings"));
        }
        fp.validate_against(source,t,b)?;
        let arcs=fp.exposed_arcs.ok_or_else(||invalid("missing exposed source arcs"))?;
        let overlap=fp.overlap_area.ok_or_else(||invalid("missing certified area"))?;
        let remaining=fp.predicted_remaining_volume
            .ok_or_else(||invalid("missing expected source volume"))?;
        let h=source.bbox.max.y-source.bbox.min.y;
        let d=b.depth;
        let lower=h-d;
        if lower<=source.tolerance.length_at(h) ||
            !lower.is_finite(){return Err(invalid("unresolved blind-floor level"));}

        // Along the upper through-cylinder wall, only the arc outside the
        // blind footprint is part of the physical material boundary. Conversely,
        // the blind-cylinder wall exists only on the arc outside the through.
        let exposed_t=arcs[0].sweep_radians;
        let exposed_b=arcs[1].sweep_radians;
        if exposed_t<=0.0 || exposed_t>=TAU ||
            exposed_b<=0.0||exposed_b>=TAU {
            return Err(invalid("invalid crossing angular partition"));
        }
        let at=PI*t.radius*t.radius;
        let ab=PI*b.radius*b.radius;
        let top= (source.bbox.max.x-source.bbox.min.x)*
            (source.bbox.max.z-source.bbox.min.z) - (at+ab-overlap);
        let bottom=(source.bbox.max.x-source.bbox.min.x)*
            (source.bbox.max.z-source.bbox.min.z)-at;
        let sides=source.mass.surface_area-
            2.0*(source.bbox.max.x-source.bbox.min.x)*
            (source.bbox.max.z-source.bbox.min.z);
        let through_lower=TAU*t.radius*lower;
        let through_upper=exposed_t*t.radius*d;
        let blind_upper=exposed_b*b.radius*d;
        let floor=ab-overlap;

        let patches=[
            MixedSurfacePatch{kind:MixedPatchKind::TopCapComplement,surface_area:top,
                is_cylindrical:false,reversed_material_sense:false},
            MixedSurfacePatch{kind:MixedPatchKind::BottomCapComplement,surface_area:bottom,
                is_cylindrical:false,reversed_material_sense:false},
            MixedSurfacePatch{kind:MixedPatchKind::OuterBoxSides,surface_area:sides,
                is_cylindrical:false,reversed_material_sense:false},
            MixedSurfacePatch{kind:MixedPatchKind::ThroughLowerWall,
                surface_area:through_lower,is_cylindrical:true,reversed_material_sense:true},
            MixedSurfacePatch{kind:MixedPatchKind::ThroughUpperWall,
                surface_area:through_upper,is_cylindrical:true,reversed_material_sense:true},
            MixedSurfacePatch{kind:MixedPatchKind::BlindUpperWall,
                surface_area:blind_upper,is_cylindrical:true,reversed_material_sense:true},
            MixedSurfacePatch{kind:MixedPatchKind::BlindFloorCrescent,
                surface_area:floor,is_cylindrical:false,reversed_material_sense:true},
        ];
        if patches.iter().any(|p|!p.surface_area.is_finite()||p.surface_area<=0.0){
            return Err(invalid("surface partition has degenerate or nonfinite measure"));
        }
        // One full-height tunnel survives this transverse mixing;
        // a blind bite into its upper wall cannot close the tunnel.
        let surface_area=patches.iter().map(|p|p.surface_area).sum::<f64>();
        Ok(MixedCrossingSurfacePlan{
            footprint:fp,patches,predicted_surface_area:surface_area,
            predicted_remaining_volume:remaining,
            expected_genus:1,expected_euler_characteristic:0,
            top_boundary_arc_count:2,bottom_boundary_arc_count:1,
            splitting_height_y:source.bbox.max.y-b.depth,
        })
    }
}
