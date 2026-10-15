//! G3: whole-solid membership of complete G2 planar UV face cells.
//!
//! A classification is Ready only when inequalities over ALL analytic
//! boundary spans (including true circular arcs) bound the entire connected
//! cell on one side of the opposite solid's material. An interior point probe
//! independently checks the claimed label. Sampling alone is never accepted
//! as proof that a cell is uniformly inside or outside.
//!
//! This bounded implementation supports canonical axis-aligned polyhedral
//! boxes and analytically oriented finite cylinders as opposite solids.
//! Other bodies, near-boundary uncertainty, or incomplete opposing face
//! coverage yield Held. No topology is mutated; can_apply is always false.

use crate::{
    BooleanOperation, BooleanRegionClassifier, BrepBody, BrepModel, FaceId,
    GeometryError, GeometryTolerance, GraphOperand, MultiTraceUvPlan,
    MultiUvBoundarySpan, MultiUvFaceCell, MultiUvBoundaryProvenance, Point2Param, Point3,
    Surface3, UvCycleKind, VolumePosition, FacePairGraph, TrimPosition,
};
use crate::trimmed_intersections::classify_trim_point;
use crate::multi_trace_uv_arrangement::contains;

fn invalid(message:impl Into<String>)->GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn uv_distance(a:Point2Param,b:Point2Param)->f64{
    (a.u-b.u).hypot(a.v-b.v)
}
fn magnitude(p:Point3)->f64 {p.x.hypot(p.y).hypot(p.z)}
fn extent(body:&BrepModel)->f64{
    magnitude(body.bbox.max.sub(body.bbox.min)).max(1.0)
}
fn max_abs(a:ScalarInterval)->f64 {a.min.abs().max(a.max.abs())}
fn min_abs(a:ScalarInterval)->f64{
    if a.min>0.0 {a.min}
    else if a.max<0.0 {-a.max}
    else {0.0}
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ScalarInterval {
    pub min:f64,
    pub max:f64,
}
impl ScalarInterval {
    fn new()->Self {
        Self {min:f64::INFINITY,max:f64::NEG_INFINITY}
    }
    fn include(&mut self,value:f64)->Result<(),GeometryError>{
        if !value.is_finite(){return Err(invalid(
            "G3 projected cell interval is not finite"));}
        self.min=self.min.min(value);
        self.max=self.max.max(value);
        Ok(())
    }
    fn is_valid(self)->bool{
        self.min.is_finite()&&self.max.is_finite()&&self.min<=self.max
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum CertifiedCellLocation { Inside, Outside }

#[derive(Debug,Clone,PartialEq)]
pub enum WholeSolidCellCertificate {
    AxisAlignedBox {
        intervals:[ScalarInterval;3],
        /// For lower and upper slab limits of each XYZ axis. A true
        /// G2 cut on the corresponding opposing B-rep face certifies
        /// representability-level endpoint rounding.
        cut_plane_support:[[bool;2];3],
    },
    AnalyticCylinder {
        axial:ScalarInterval,
        transverse:[ScalarInterval;2],
        radius_lower_bound:f64,
        radius_upper_bound:f64,
    },
}
#[derive(Debug,Clone,PartialEq)]
pub struct WholeSolidCellWitness {
    pub cell_index:usize,
    pub uv:Point2Param,
    pub point:Point3,
    pub location:CertifiedCellLocation,
    pub keep:bool,
    /// Only B-side patches of Difference have their normal reversed.
    pub reverse_orientation:bool,
    pub certificate:WholeSolidCellCertificate,
}
#[derive(Debug,Clone,PartialEq)]
pub struct WholeSolidRegionPlan {
    pub operation:BooleanOperation,
    pub source_face:FaceId,
    pub source_operand:GraphOperand,
    pub witnesses:Vec<WholeSolidCellWitness>,
    pub kept_cells:usize,
    /// G4 topology reconstruction and G5 atomic commit are still required.
    pub can_apply:bool,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum WholeSolidRegionHold {
    IncompleteOppositeFaceCoverage,
    UnsupportedOppositeSolid,
    UnsupportedSourceSurface,
    SourceFaceBoundary,
    NoInteriorWitness,
    NonUniformCellBounds,
    UnresolvedNumericBoundary,
    ContradictoryMembership,
    BudgetExceeded,
}
#[derive(Debug,Clone,PartialEq)]
pub enum WholeSolidRegionResult {
    Ready(WholeSolidRegionPlan),
    Held(WholeSolidRegionHold),
}

fn keep(operation:BooleanOperation,operand:GraphOperand,location:CertifiedCellLocation)->bool{
    match(operation,operand,location) {
        (BooleanOperation::Union,_,CertifiedCellLocation::Outside)=>true,
        (BooleanOperation::Intersection,_,CertifiedCellLocation::Inside)=>true,
        (BooleanOperation::Difference,GraphOperand::First,CertifiedCellLocation::Outside)=>true,
        (BooleanOperation::Difference,GraphOperand::Second,CertifiedCellLocation::Inside)=>true,
        _=>false,
    }
}

/// Exact-in-form extrema of one linear 3D functional over a G2 boundary
/// span. A circle is never replaced with chords.
fn projected_range(
    span:MultiUvBoundarySpan,origin:Point3,
    x:Point3,y:Point3,normal:Point3,
)->Result<ScalarInterval,GeometryError>{
    let mut out=ScalarInterval::new();
    let base=normal.dot(origin);
    let a=normal.dot(x);
    let b=normal.dot(y);
    let eval=|uv:Point2Param|base+a*uv.u+b*uv.v;
    out.include(eval(span.start))?;
    out.include(eval(span.end))?;
    if let Some(circle)=span.circular {
        let amplitude=a.hypot(b);
        if !amplitude.is_finite()||!circle.radius.is_finite()||circle.radius<=0.0 {
            return Err(invalid("G3 circle has invalid analytic radius"));
        }
        if amplitude>0.0 {
            let theta=b.atan2(a);
            // A circle arc may be unwrapped negatively across a chart seam.
            // Evaluate extrema only when the critical angle lies on this arc.
            let lo=circle.angle_start.min(circle.angle_end);
            let hi=circle.angle_start.max(circle.angle_end);
            let phase=[theta,theta+std::f64::consts::PI];
            for critical in phase {
                for wrap in -3..=3 {
                    let t=critical+(wrap as f64)*std::f64::consts::TAU;
                    if t>=lo&&t<=hi {
                        out.include(base+a*(circle.center.u+circle.radius*t.cos())
                            +b*(circle.center.v+circle.radius*t.sin()))?;
                    }
                }
            }
        }
    }
    if !out.is_valid(){return Err(invalid("G3 cell range is inverted"));}
    Ok(out)
}

fn cell_range(cell:&MultiUvFaceCell,source:Surface3,axis:Point3)
    ->Result<ScalarInterval,GeometryError>{
    let Surface3::Plane(plane)=source else{
        return Err(invalid("G3 only supports source planar UV cells"));
    };
    let mut interval=ScalarInterval::new();
    for span in cell.cycles.iter().flat_map(|cycle|&cycle.spans) {
        let bound=projected_range(
            *span,plane.frame.origin,plane.frame.x,plane.frame.y,axis)?;
        interval.include(bound.min)?;
        interval.include(bound.max)?;
    }
    if !interval.is_valid(){return Err(invalid("empty G3 UV cell boundary"));}
    Ok(interval)
}

/// A point on an independently validated 3D face-pair cut may acquire a
/// tiny coordinate overshoot when projected through the source UV frame.
/// Only such physical cut curves, lying on the REAL opposing box face,
/// authorize a representability-sized adjustment of the interval decision.
fn cut_plane_support(
    cell:&MultiUvFaceCell, source_surface:Surface3,
    source_operand:GraphOperand, pairs:&[FacePairGraph],
    opposite:&BrepModel, lower:[f64;3],upper:[f64;3],
    precision:f64,
)->Result<[[bool;2];3],GeometryError>{
    let mut supported=[[false;2];3];
    let Surface3::Plane(source_plane)=source_surface else{
        return Err(invalid("G3 cut support requires a planar source face"));
    };
    for span in cell.cycles.iter().flat_map(|c|&c.spans) {
        let MultiUvBoundaryProvenance::Cut{source,..}=span.provenance else{
            continue;
        };
        let pair=pairs.iter().find(|pair|pair.id==source.pair)
            .ok_or_else(||invalid("G3 cut has no source face pair"))?;
        if pair.graph.arcs.get(source.local.0 as usize)
            .map(|arc|arc.id)!=Some(source.local) {
            return Err(invalid("G3 cut has no source analytic arc"));
        }
        let opposing_face=if source_operand==GraphOperand::First {
            pair.second_face
        }else{pair.first_face};
        let Some(Surface3::Plane(plane))=opposite.face_surface(opposing_face) else {
            continue;
        };
        let normal=plane.frame.normal.array();
        let origin=plane.frame.origin.array();
        let first=source_plane.frame.position(span.start.u,span.start.v)?.array();
        let last=source_plane.frame.position(span.end.u,span.end.v)?.array();
        for axis in 0..3 {
            if normal[axis].abs()<1.0-1e-10 ||
                (0..3).any(|i|i!=axis && normal[i].abs()>1e-10) {
                continue;
            }
            for (side,bound) in [lower[axis],upper[axis]].into_iter().enumerate(){
                if (origin[axis]-bound).abs()<=precision &&
                    (first[axis]-bound).abs()<=precision &&
                    (last[axis]-bound).abs()<=precision {
                    supported[axis][side]=true;
                }
            }
        }
    }
    Ok(supported)
}

fn mapped_location(
    cell:&MultiUvFaceCell,source:Surface3,other:&BrepBody,
    source_operand:GraphOperand,pairs:&[FacePairGraph],
    other_model:&BrepModel,precision:f64,
)->Result<Result<(CertifiedCellLocation,WholeSolidCellCertificate),WholeSolidRegionHold>,
    GeometryError>
{
    match other {
        BrepBody::Polyhedral(box_solid)=>{
            // The legacy Boolean kernel verifies the *actual* six planar
            // faces, orientation and axis alignment, rather than trusting
            // the bounding box of an arbitrary polyhedron.
            if crate::boolean::certify_box(box_solid).is_err() {
                return Ok(Err(WholeSolidRegionHold::UnsupportedOppositeSolid));
            }
            let axes=[
                Point3{x:1.0,y:0.0,z:0.0},
                Point3{x:0.0,y:1.0,z:0.0},
                Point3{x:0.0,y:0.0,z:1.0},
            ];
            let mut intervals=[ScalarInterval{min:0.0,max:0.0};3];
            for (i,&axis) in axes.iter().enumerate(){
                intervals[i]=cell_range(cell,source,axis)?;
            }
            let lower=box_solid.bbox.min.array();
            let upper=box_solid.bbox.max.array();
            let cut_support=cut_plane_support(
                cell,source,source_operand,pairs,other_model,
                lower,upper,precision)?;
            // Keep the geometric inequality STRICT. Floating point budget
            // distinguishes a numeric Held from a wider unresolved region;
            // it must never swallow a real thin inside/outside strip.
            let above=|actual:f64,limit:f64,allow:bool|
                actual>=limit ||
                (allow && actual>=limit-precision);
            let below=|actual:f64,limit:f64,allow:bool|
                actual<=limit ||
                (allow && actual<=limit+precision);
            let outside=(0..3).any(|i|
                below(intervals[i].max,lower[i],cut_support[i][0]) ||
                above(intervals[i].min,upper[i],cut_support[i][1]));
            let inside=(0..3).all(|i|
                above(intervals[i].min,lower[i],cut_support[i][0]) &&
                below(intervals[i].max,upper[i],cut_support[i][1]));
            let near=(0..3).any(|i|[
                (intervals[i].max-lower[i]).abs(),
                (intervals[i].min-upper[i]).abs(),
                (intervals[i].min-lower[i]).abs(),
                (intervals[i].max-upper[i]).abs(),
            ].into_iter().any(|d|d<=precision));
            let result=if outside {
                Some(CertifiedCellLocation::Outside)
            }else if inside {
                Some(CertifiedCellLocation::Inside)
            }else{None};
            Ok(result.map(|location|(
                location,WholeSolidCellCertificate::AxisAlignedBox{
                    intervals,cut_plane_support:cut_support,
                }
            )).ok_or(if near {
                WholeSolidRegionHold::UnresolvedNumericBoundary
            }else {
                WholeSolidRegionHold::NonUniformCellBounds
            }))
        }
        BrepBody::Cylinder(cylinder)=>{
            let frame=cylinder.frame;
            let axis=cell_range(cell,source,frame.normal)?;
            let x=cell_range(cell,source,frame.x)?;
            let y=cell_range(cell,source,frame.y)?;
            let axial=ScalarInterval {
                min:axis.min-frame.normal.dot(frame.origin),
                max:axis.max-frame.normal.dot(frame.origin),
            };
            let px=ScalarInterval {
                min:x.min-frame.x.dot(frame.origin),
                max:x.max-frame.x.dot(frame.origin),
            };
            let py=ScalarInterval {
                min:y.min-frame.y.dot(frame.origin),
                max:y.max-frame.y.dot(frame.origin),
            };
            let lower=min_abs(px).hypot(min_abs(py));
            let upper=max_abs(px).hypot(max_abs(py));
            if !lower.is_finite() || !upper.is_finite() {
                return Err(invalid("nonfinite G3 cylinder interval radius"));
            }
            // Conservative enclosure: a Cartesian rectangle contains the
            // entire curved cell; it may be wider and cause a safe Held.
            let outside=axial.max<=0.0 ||
                axial.min>=cylinder.height ||
                lower>=cylinder.radius;
            let inside=axial.min>=0.0 &&
                axial.max<=cylinder.height &&
                upper<=cylinder.radius;
            let near=[
                axial.min.abs(),axial.max.abs(),
                (axial.min-cylinder.height).abs(),
                (axial.max-cylinder.height).abs(),
                (lower-cylinder.radius).abs(),
                (upper-cylinder.radius).abs(),
            ].into_iter().any(|distance|distance<=precision);
            let location=if outside {
                Some(CertifiedCellLocation::Outside)
            }else if inside {
                Some(CertifiedCellLocation::Inside)
            }else{None};
            Ok(location.map(|location|(
                location,
                WholeSolidCellCertificate::AnalyticCylinder {
                    axial,transverse:[px,py],
                    radius_lower_bound:lower,radius_upper_bound:upper,
                },
            )).ok_or(if near {
                WholeSolidRegionHold::UnresolvedNumericBoundary
            }else {
                WholeSolidRegionHold::NonUniformCellBounds
            }))
        }
        _=>Ok(Err(WholeSolidRegionHold::UnsupportedOppositeSolid)),
    }
}

/// The source G2 cycle keeps material on the left of its directed boundary.
/// Finding a point slightly left of several different analytic spans avoids
/// cutting through holes or sampling exactly on the opposite solid's shell.
fn interior_probe(
    cell:&MultiUvFaceCell,source_face:FaceId,
    source_model:&BrepModel,source_surface:Surface3,
    other:&BrepBody,tolerance:GeometryTolerance,
)->Result<Result<(Point2Param,Point3,VolumePosition),WholeSolidRegionHold>,GeometryError>{
    let Surface3::Plane(plane)=source_surface else {
        return Ok(Err(WholeSolidRegionHold::UnsupportedSourceSurface));
    };
    let eps=tolerance.length_at(extent(source_model));
    for cycle in &cell.cycles {
        if cycle.kind!=UvCycleKind::Outer {continue;}
        for &span in &cycle.spans {
            let midpoint=span.midpoint();
            let tangent=if let Some(circle)=span.circular {
                let angle=(circle.angle_start+circle.angle_end)*0.5;
                let sign=if circle.angle_end>=circle.angle_start {1.0}else{-1.0};
                Point2Param {u:-angle.sin()*sign,v:angle.cos()*sign}
            }else{
                Point2Param {
                    u:span.end.u-span.start.u,
                    v:span.end.v-span.start.v,
                }
            };
            let size=tangent.u.hypot(tangent.v);
            let distance=if let Some(c)=span.circular {
                c.radius*(c.angle_end-c.angle_start).abs()
            }else {uv_distance(span.start,span.end)};
            if !size.is_finite()||size<=0.0 ||
                !distance.is_finite()||distance<=eps*32.0 {continue;}
            for factor in [0.1,0.025,0.00625,0.0015625] {
                let shift=distance*factor;
                if shift<=eps*64.0 {continue;}
                let sample=Point2Param {
                    u:midpoint.u-tangent.v/size*shift,
                    v:midpoint.v+tangent.u/size*shift,
                };
                let within_outer=contains(cycle,sample,eps);
                if within_outer!=Ok(true){continue;}
                let mut in_hole=false;
                for hole in &cell.cycles {
                    if hole.kind!=UvCycleKind::Inner{continue;}
                    if contains(hole,sample,eps)!=Ok(false) {
                        in_hole=true;break;
                    }
                }
                if in_hole {continue;}
                let point=plane.frame.position(sample.u,sample.v)?;
                if classify_trim_point(source_model,source_face,point,tolerance)?
                    !=TrimPosition::Interior {continue;}
                let location=BooleanRegionClassifier::classify_point(other,point)?;
                if location==VolumePosition::Boundary{continue;}
                return Ok(Ok((sample,point,location)));
            }
        }
    }
    Ok(Err(WholeSolidRegionHold::NoInteriorWitness))
}

fn precision_budget(source:&BrepModel,other:&BrepModel)->Result<f64,GeometryError>{
    let max_coord=[
        source.bbox.min,source.bbox.max,other.bbox.min,other.bbox.max,
    ].into_iter()
        .flat_map(|point|[point.x.abs(),point.y.abs(),point.z.abs()])
        .fold(1.0_f64,f64::max);
    let precision=64.0*f64::EPSILON*max_coord;
    if !precision.is_finite(){return Err(invalid("G3 floating point budget overflow"));}
    Ok(precision)
}

impl WholeSolidRegionPlan {
    /// Reject even a source-valid UV preview when a changed body or an
    /// edited operation would alter the regenerated material classification.
    pub fn validate_against(
        &self,left:&BrepBody,right:&BrepBody,
        pairs:&[FacePairGraph],tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        let uv=MultiTraceUvPlan::build_selected(
            &left.shared()?,&right.shared()?,
            pairs,self.source_operand,self.source_face,tolerance)?;
        let crate::MultiTraceUvResult::Ready(uv)=uv else{
            return Err(invalid("G3 source no longer constructs a ready UV face"));
        };
        let rebuilt=Self::build(left,right,pairs,&uv,self.operation,tolerance)?;
        if rebuilt!=WholeSolidRegionResult::Ready(self.clone()) {
            return Err(invalid("G3 classified plan differs from current source solids"));
        }
        Ok(())
    }

    /// Selected source pairs MUST cover every face of the opposite body,
    /// including pairs whose intersection graph is Empty. Coverage cannot
    /// be inferred from the pairs that happen to have analytic cut arcs.
    pub fn build(
        left:&BrepBody,right:&BrepBody,
        pairs:&[FacePairGraph],uv:&MultiTraceUvPlan,
        operation:BooleanOperation,tolerance:GeometryTolerance,
    )->Result<WholeSolidRegionResult,GeometryError>{
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let lm=left.shared()?;
        let rm=right.shared()?;
        uv.validate_against(&lm,&rm,pairs,tolerance)?;

        let (model,other,other_model)=if uv.operand==GraphOperand::First {
            (&lm,right,&rm)
        }else{(&rm,left,&lm)};
        // Selected-pair completeness prevents a collection of partial cuts
        // from masquerading as a partition of the whole opposite solid.
        let mut covered=vec![false;other_model.faces.len()];
        for pair in pairs {
            let selected=if uv.operand==GraphOperand::First {
                pair.first_face==uv.face
            }else{pair.second_face==uv.face};
            if !selected{continue;}
            let target=if uv.operand==GraphOperand::First {
                pair.second_face
            }else{pair.first_face};
            let Some(entry)=covered.get_mut(target.0 as usize) else {
                return Err(invalid("G3 selected opposite face is missing"));
            };
            if *entry {
                return Err(invalid("G3 duplicate opposite face source pair"));
            }
            *entry=true;
        }
        if !covered.iter().all(|x|*x) {
            return Ok(WholeSolidRegionResult::Held(
                WholeSolidRegionHold::IncompleteOppositeFaceCoverage));
        }
        if uv.cells.len()>256 {
            return Ok(WholeSolidRegionResult::Held(
                WholeSolidRegionHold::BudgetExceeded));
        }
        let Some(source_surface)=model.face_surface(uv.face) else{
            return Err(invalid("G3 missing source surface"));
        };
        if !matches!(source_surface,Surface3::Plane(_)) {
            return Ok(WholeSolidRegionResult::Held(
                WholeSolidRegionHold::UnsupportedSourceSurface));
        }
        let precision=precision_budget(model,other_model)?;
        let mut witnesses=Vec::with_capacity(uv.cells.len());
        for (cell_index,cell) in uv.cells.iter().enumerate() {
            let (label,certificate)=match mapped_location(
                cell,source_surface,other,uv.operand,pairs,other_model,precision)? {
                Ok(proof)=>proof,
                Err(why)=>return Ok(WholeSolidRegionResult::Held(why)),
            };
            let (point_uv,point,observed)=match interior_probe(
                cell,uv.face,model,source_surface,other,tolerance)? {
                Ok(data)=>data,
                Err(why)=>return Ok(WholeSolidRegionResult::Held(why)),
            };
            let expected=match label {
                CertifiedCellLocation::Inside=>VolumePosition::Inside,
                CertifiedCellLocation::Outside=>VolumePosition::Outside,
            };
            if observed!=expected {
                return Ok(WholeSolidRegionResult::Held(
                    WholeSolidRegionHold::ContradictoryMembership));
            }
            witnesses.push(WholeSolidCellWitness {
                cell_index,uv:point_uv,point,location:label,
                keep:keep(operation,uv.operand,label),
                reverse_orientation:operation==BooleanOperation::Difference &&
                    uv.operand==GraphOperand::Second,
                certificate,
            });
        }
        let kept_cells=witnesses.iter().filter(|cell|cell.keep).count();
        Ok(WholeSolidRegionResult::Ready(Self {
            operation,source_face:uv.face,source_operand:uv.operand,
            witnesses,kept_cells,can_apply:false,
        }))
    }
}

#[cfg(test)]
mod g3_unit_tests {
    use super::*;
    use crate::{Frame3, UvCircularSpan, MultiUvBoundaryProvenance,
        GlobalArcRef, FacePairId, GraphArcId};

    fn q(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
    fn xyz(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}

    #[test]
    fn circular_arc_extrema_are_exact_in_form_not_endpoint_only(){
        let arc=MultiUvBoundarySpan {
            start:q(1.,0.),end:q(-1.,0.),
            provenance:MultiUvBoundaryProvenance::Cut {
                segment:0,source:GlobalArcRef {
                    pair:FacePairId(0),local:GraphArcId(0),
                },
            },
            circular:Some(UvCircularSpan{
                center:q(0.,0.),radius:1.,
                angle_start:0.0,angle_end:std::f64::consts::PI,
            }),
        };
        let range=projected_range(
            arc,xyz(0.,0.,0.),xyz(1.,0.,0.),xyz(0.,1.,0.),
            xyz(0.,1.,0.)).unwrap();
        assert!((range.min-0.0).abs()<1e-14);
        assert!((range.max-1.0).abs()<1e-14,
            "upper semicircle must reach y=1 internally");
    }

    #[test]
    fn entire_negative_oriented_circle_has_both_radial_extrema(){
        let span=MultiUvBoundarySpan {
            start:q(2.,0.),end:q(2.,0.),
            provenance:MultiUvBoundaryProvenance::Cut {
                segment:0,source:GlobalArcRef{
                    pair:FacePairId(0),local:GraphArcId(0),
                },
            },
            circular:Some(UvCircularSpan {
                center:q(0.,0.),radius:2.,
                angle_start:0.0,angle_end:-std::f64::consts::TAU,
            }),
        };
        let range=projected_range(
            span,xyz(0.,0.,0.),xyz(1.,0.,0.),xyz(0.,1.,0.),
            xyz(1.,0.,0.)).unwrap();
        assert!((range.min+2.).abs()<1e-14);
        assert!((range.max-2.).abs()<1e-14);
    }

    #[test]
    fn difference_reverses_only_second_operand_inside_patches(){
        assert!(keep(BooleanOperation::Difference,GraphOperand::First,
            CertifiedCellLocation::Outside));
        assert!(!keep(BooleanOperation::Difference,GraphOperand::First,
            CertifiedCellLocation::Inside));
        assert!(keep(BooleanOperation::Difference,GraphOperand::Second,
            CertifiedCellLocation::Inside));
        assert!(!keep(BooleanOperation::Difference,GraphOperand::Second,
            CertifiedCellLocation::Outside));
        assert!(keep(BooleanOperation::Intersection,GraphOperand::Second,
            CertifiedCellLocation::Inside));
        assert!(keep(BooleanOperation::Union,GraphOperand::Second,
            CertifiedCellLocation::Outside));
    }

    #[test]
    fn source_space_interval_is_translation_invariant_in_exact_coordinates(){
        let frame=Frame3::from_axes(
            xyz(100_000.,200_000.,-300_000.),
            xyz(1.,0.,0.),xyz(0.,0.,1.),
            crate::GeometryTolerance::default(),
        ).unwrap();
        let interval=projected_range(
            MultiUvBoundarySpan {
                start:q(-3.,-2.),end:q(4.,5.),
                provenance:MultiUvBoundaryProvenance::Cut {
                    segment:0,source:GlobalArcRef{
                        pair:FacePairId(0),local:GraphArcId(0),
                    },
                },
                circular:None,
            },
            frame.origin,frame.x,frame.y,xyz(1.,0.,0.),
        ).unwrap();
        assert_eq!(interval.min,99_997.);
        assert_eq!(interval.max,100_004.);
    }
}
