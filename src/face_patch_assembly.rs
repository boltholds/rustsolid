//! G4a: provenance-carrying, read-only local face patch/coedge draft.
//!
//! G1 reconciles source intersection events, G2 subdivides one planar source
//! face into analytic UV cells, and G3 proves the complete cells' membership
//! against the opposite *whole* solid. G4a prepares the directed uses needed
//! for a later global topology builder. A G2 cut is not a physical B-rep EdgeId.
//!
//! A complete closed shell additionally requires *all* source faces, true
//! cross-face/cooperand cut twins, shared physical vertices, an independent
//! mass oracle and an atomic commit (G4b/G5). None is claimed here.

use crate::{
    BrepBody, BrepModel, BooleanOperation, Curve2, FaceId, FacePairGraph,
    FaceSense, GeometryError, GeometryTolerance, GlobalArcRef, GraphOperand,
    IntersectionCurve, MultiTraceUvPlan, MultiTraceUvResult,
    MultiUvBoundaryProvenance, MultiUvBoundarySpan, Point2Param, Point3,
    Surface3, UvCycleKind, WholeSolidRegionPlan, WholeSolidRegionResult,
};

const MAX_PATCH_USES: usize = 16_384;
const MAX_PATCH_CUTS: usize = 2_048;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn length(a: Point3, b: Point3) -> f64 {
    let d = a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}
fn uv_distance(a: Point2Param, b: Point2Param) -> f64 {
    (a.u - b.u).hypot(a.v - b.v)
}
fn spatial_span(
    surface: Surface3,
    span: MultiUvBoundarySpan,
) -> Result<(Point3, Point3), GeometryError> {
    let a = surface.evaluate(span.start.u, span.start.v)?.point;
    let b = surface.evaluate(span.end.u, span.end.v)?.point;
    if !a.is_finite() || !b.is_finite() {
        return Err(invalid("G4 surface evaluated to nonfinite coordinates"));
    }
    Ok((a,b))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FacePatchAssemblyHold {
    BudgetExceeded,
    UnsupportedSourceSurface,
    UnsupportedCutCarrier,
    SourceTrimMismatch,
    CutCarrierMismatch,
    CutIncidenceMismatch,
    NonSeparatingCut,
    OpenOrMisorientedCycle,
}

/// A preview-local directed source edge or cut use. Neither `id` nor any
/// index in `FacePatchAssemblyPlan` is a persistent topological name.
#[derive(Debug, Clone, PartialEq)]
pub struct FacePatchDraftUse {
    pub cell_index: usize,
    pub cycle_index: usize,
    pub span_index: usize,
    pub provenance: MultiUvBoundaryProvenance,
    pub uv_span: MultiUvBoundarySpan,
    pub start: Point3,
    pub end: Point3,
    pub kept: bool,
    /// Includes the original source face sense and the G3 Difference flip.
    pub reversed_for_shell: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FacePatchDraftLoop {
    pub cell_index: usize,
    pub cycle_index: usize,
    pub kind: UvCycleKind,
    pub signed_area_uv: f64,
    pub uses: Vec<usize>,
    pub reversed_for_shell: bool,
}

/// Two opposite local uses of one G2 cut. Exactly one is selected by G3.
/// Its real physical twin on the *other source face* is still outstanding.
#[derive(Debug, Clone, PartialEq)]
pub struct FacePatchCutPair {
    pub segment: u32,
    pub source_arc: GlobalArcRef,
    pub counterpart_operand: GraphOperand,
    pub counterpart_face: FaceId,
    pub uses: [usize;2],
    pub selected_use: usize,
    pub discarded_use: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FacePatchAssemblyPlan {
    pub source_operand: GraphOperand,
    pub source_face: FaceId,
    pub operation: BooleanOperation,
    /// All input spans are preserved for audit, including discarded ones.
    pub uses: Vec<FacePatchDraftUse>,
    pub retained_loops: Vec<FacePatchDraftLoop>,
    pub cut_pairs: Vec<FacePatchCutPair>,
    pub kept_area_uv: f64,
    /// Source coedge spans retained on the selected patches; full source
    /// adjacency is required before these may become shared physical edges.
    pub open_source_uses: usize,
    /// Every locally selected physical cut needs a cross-face coedge twin.
    pub open_cut_uses: usize,
    /// No temporary patch may be swapped into a persisted B-rep.
    pub can_apply: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FacePatchAssemblyResult {
    Ready(FacePatchAssemblyPlan),
    Held(FacePatchAssemblyHold),
}

fn check_source_span(
    model: &BrepModel,
    face: FaceId,
    span: MultiUvBoundarySpan,
    epsilon: f64,
) -> Result<bool, GeometryError> {
    let MultiUvBoundaryProvenance::Source {
        edge, coedge, loop_id, parameter_start, parameter_end,
    } = span.provenance else { return Ok(false); };

    let Some(source_use) = model.coedges.get(coedge.0 as usize) else {
        return Ok(false);
    };
    if source_use.id != coedge || source_use.face != face ||
        source_use.edge != edge || source_use.loop_id != loop_id ||
        model.edge_curve(edge).is_none()
    {
        return Ok(false);
    }
    let Some(curve) = model.coedge_pcurve(coedge) else {return Ok(false);};
    if !matches!(curve, Curve2::Line(_) | Curve2::Circle(_)) ||
        !parameter_start.is_finite() || !parameter_end.is_finite()
    {
        return Ok(false);
    }
    let a = curve.evaluate(parameter_start)?;
    let b = curve.evaluate(parameter_end)?;
    Ok(uv_distance(a, span.start) <= epsilon &&
       uv_distance(b, span.end) <= epsilon)
}

/// Verify the *bounded* physical 3D intersection carrier. UV proximity to
/// another source face is never enough to invent an EdgeId or merge vertices.
fn check_cut_span(
    arc: crate::GraphArc,
    a: Point3,
    b: Point3,
    midpoint: Point3,
    budget: f64,
) -> Result<bool, GeometryError> {
    let IntersectionCurve::Line(line) = arc.curve else { return Ok(false); };
    let len2 = line.direction.length_squared();
    if !len2.is_finite() || len2 <= 0.0 ||
        !arc.parameter_start.is_finite() || !arc.parameter_end.is_finite() ||
        arc.parameter_start >= arc.parameter_end || !budget.is_finite() || budget <= 0.0
    {
        return Ok(false);
    }
    let slack = budget / len2.sqrt();
    for point in [a,b,midpoint] {
        let parameter = point.sub(line.origin).dot(line.direction)/len2;
        if !parameter.is_finite() ||
            parameter < arc.parameter_start - slack ||
            parameter > arc.parameter_end + slack ||
            length(point, line.evaluate(parameter)?) > budget
        {
            return Ok(false);
        }
    }
    Ok(true)
}

impl FacePatchAssemblyPlan {
    pub fn validate_against(
        &self,
        left: &BrepBody, right: &BrepBody,
        pairs: &[FacePairGraph], tolerance: GeometryTolerance,
    ) -> Result<(), GeometryError> {
        let a = left.shared()?;
        let b = right.shared()?;
        let rebuilt_uv = MultiTraceUvPlan::build_selected(
            &a, &b, pairs, self.source_operand, self.source_face, tolerance
        )?;
        let MultiTraceUvResult::Ready(uv) = rebuilt_uv else {
            return Err(invalid("G4 original source no longer constructs G2"));
        };
        let rebuilt_g3 = WholeSolidRegionPlan::build(
            left, right, pairs, &uv, self.operation, tolerance
        )?;
        let WholeSolidRegionResult::Ready(g3) = rebuilt_g3 else {
            return Err(invalid("G4 original source no longer constructs G3"));
        };
        if Self::build(left, right, pairs, &uv, &g3, tolerance)? !=
            FacePatchAssemblyResult::Ready(self.clone())
        {
            return Err(invalid("G4 draft is not identical to source reconstruction"));
        }
        Ok(())
    }

    pub fn build(
        left: &BrepBody, right: &BrepBody,
        pairs: &[FacePairGraph], uv: &MultiTraceUvPlan,
        g3: &WholeSolidRegionPlan, tolerance: GeometryTolerance,
    ) -> Result<FacePatchAssemblyResult, GeometryError> {
        tolerance.validate()?;
        left.validate()?;
        right.validate()?;
        let a = left.shared()?;
        let b = right.shared()?;
        uv.validate_against(&a, &b, pairs, tolerance)?;
        g3.validate_against(left, right, pairs, tolerance)?;
        if uv.face != g3.source_face || uv.operand != g3.source_operand ||
            g3.witnesses.len() != uv.cells.len() ||
            g3.kept_cells != g3.witnesses.iter().filter(|w| w.keep).count()
        {
            return Err(invalid("G4 face/cell classification provenance mismatch"));
        }
        let model = if uv.operand == GraphOperand::First {&a} else {&b};
        let Some(surface) = model.face_surface(uv.face) else {
            return Err(invalid("G4 source face has no analytic support"));
        };
        if !matches!(surface, Surface3::Plane(_)) {
            return Ok(FacePatchAssemblyResult::Held(
                FacePatchAssemblyHold::UnsupportedSourceSurface
            ));
        }
        let Some(sense) = model.face_sense(uv.face) else {
            return Err(invalid("G4 source face sense is absent"));
        };
        let scale = model.bbox.max.sub(model.bbox.min).length_squared().sqrt();
        if !scale.is_finite() {return Err(invalid("G4 invalid source scale"));}
        let eps = tolerance.length_at(scale.max(1.0)) * 16.0;
        let mut uses = Vec::<FacePatchDraftUse>::new();
        let mut loops = Vec::<FacePatchDraftLoop>::new();
        let mut cut_uses = vec![Vec::<usize>::new(); uv.cuts.len()];
        if uv.cuts.len() > MAX_PATCH_CUTS {
            return Ok(FacePatchAssemblyResult::Held(FacePatchAssemblyHold::BudgetExceeded));
        }
        let mut kept_area_uv = 0.0;
        let mut open_source_uses = 0usize;
        for (cell_index,cell) in uv.cells.iter().enumerate() {
            let witness = &g3.witnesses[cell_index];
            if witness.cell_index != cell_index {
                return Err(invalid("G4 witness index does not match G2 source cell"));
            }
            let reversed = (sense == FaceSense::Reversed) ^ witness.reverse_orientation;
            if witness.keep {
                kept_area_uv += cell.area;
            }
            for (cycle_index,cycle) in cell.cycles.iter().enumerate() {
                if cycle.spans.is_empty() || uses.len().saturating_add(cycle.spans.len()) > MAX_PATCH_USES {
                    return Ok(FacePatchAssemblyResult::Held(FacePatchAssemblyHold::BudgetExceeded));
                }
                let mut cycle_uses = Vec::with_capacity(cycle.spans.len());
                for (span_index,&span) in cycle.spans.iter().enumerate() {
                    let next = cycle.spans[(span_index+1)%cycle.spans.len()];
                    if uv_distance(span.end,next.start) > eps {
                        return Ok(FacePatchAssemblyResult::Held(
                            FacePatchAssemblyHold::OpenOrMisorientedCycle
                        ));
                    }
                    let (mut start,mut end) = spatial_span(surface,span)?;
                    match span.provenance {
                        MultiUvBoundaryProvenance::Source{..} => {
                            if !check_source_span(model,uv.face,span,eps)? {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::SourceTrimMismatch
                                ));
                            }
                            if witness.keep {open_source_uses+=1;}
                        },
                        MultiUvBoundaryProvenance::Cut{segment,source} => {
                            let Some(cut) = uv.cuts.get(segment as usize) else {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutIncidenceMismatch
                                ));
                            };
                            if cut.segment != segment || cut.source != source || span.circular.is_some() {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutIncidenceMismatch
                                ));
                            }
                            let Some(pair) = pairs.iter().find(|pair|pair.id==source.pair) else {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutIncidenceMismatch
                                ));
                            };
                            if (uv.operand==GraphOperand::First && pair.first_face!=uv.face) ||
                               (uv.operand==GraphOperand::Second && pair.second_face!=uv.face)
                            {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutIncidenceMismatch
                                ));
                            }
                            let Some(arc) = pair.graph.arcs.iter().find(|arc|arc.id==source.local) else{
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutIncidenceMismatch
                                ));
                            };
                            if !matches!(arc.curve,IntersectionCurve::Line(_)) {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::UnsupportedCutCarrier
                                ));
                            }
                            let midpoint = surface.evaluate(
                                (span.start.u+span.end.u)*0.5,
                                (span.start.v+span.end.v)*0.5,
                            )?.point;
                            if !check_cut_span(*arc,start,end,midpoint,
                                arc.tolerance.evaluation_budget.max(eps))?
                            {
                                return Ok(FacePatchAssemblyResult::Held(
                                    FacePatchAssemblyHold::CutCarrierMismatch
                                ));
                            }
                            cut_uses[segment as usize].push(uses.len());
                        },
                    }
                    if reversed {
                        std::mem::swap(&mut start,&mut end);
                    }
                    let id = uses.len();
                    uses.push(FacePatchDraftUse {
                        cell_index, cycle_index,span_index,
                        provenance:span.provenance,uv_span:span,
                        start,end,kept:witness.keep,reversed_for_shell:reversed,
                    });
                    cycle_uses.push(id);
                }
                if witness.keep {
                    if reversed {cycle_uses.reverse();}
                    loops.push(FacePatchDraftLoop {
                        cell_index,cycle_index,kind:cycle.kind,
                        signed_area_uv:cycle.signed_area,
                        uses:cycle_uses,reversed_for_shell:reversed,
                    });
                }
            }
        }
        if !kept_area_uv.is_finite() || kept_area_uv < 0.0 {
            return Err(invalid("G4 invalid kept face area"));
        }

        let mut cut_pairs = Vec::with_capacity(uv.cuts.len());
        for (index,cut) in uv.cuts.iter().enumerate() {
            if cut.segment as usize != index || cut_uses[index].len()!=2 {
                return Ok(FacePatchAssemblyResult::Held(
                    FacePatchAssemblyHold::CutIncidenceMismatch
                ));
            }
            let one=cut_uses[index][0];
            let two=cut_uses[index][1];
            let x=&uses[one];
            let y=&uses[two];
            if x.cell_index == y.cell_index ||
                x.uv_span.provenance!=y.uv_span.provenance ||
                uv_distance(x.uv_span.start,y.uv_span.end)>eps ||
                uv_distance(x.uv_span.end,y.uv_span.start)>eps
            {
                return Ok(FacePatchAssemblyResult::Held(
                    FacePatchAssemblyHold::CutIncidenceMismatch
                ));
            }
            // A true volume-separating G2 cut must change the G3 material
            // label. Equal labels cannot be treated as an implicit welded seam.
            if x.kept == y.kept {
                return Ok(FacePatchAssemblyResult::Held(
                    FacePatchAssemblyHold::NonSeparatingCut
                ));
            }
            let (selected,discarded) = if x.kept {(one,two)}else{(two,one)};
            let Some(pair)=pairs.iter().find(|pair|pair.id==cut.source.pair) else {
                return Ok(FacePatchAssemblyResult::Held(
                    FacePatchAssemblyHold::CutIncidenceMismatch
                ));
            };
            let (other_operand,counterpart_face) =
                if uv.operand==GraphOperand::First {
                    (GraphOperand::Second,pair.second_face)
                }else{(GraphOperand::First,pair.first_face)};
            cut_pairs.push(FacePatchCutPair {
                segment:cut.segment,source_arc:cut.source,
                counterpart_operand:other_operand,counterpart_face,
                uses:[one,two],selected_use:selected,discarded_use:discarded,
            });
        }
        Ok(FacePatchAssemblyResult::Ready(Self {
            source_operand:uv.operand,source_face:uv.face,
            operation:g3.operation,uses,retained_loops:loops,
            cut_pairs,kept_area_uv,open_source_uses,
            open_cut_uses:uv.cuts.len(),can_apply:false,
        }))
    }
}
