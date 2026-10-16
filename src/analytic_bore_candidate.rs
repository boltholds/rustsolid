//! G4b2b1: source-checked curved B-rep candidates for strict vertical bores.
//!
//! The existing ThroughHoleBrep/BlindHoleBrep producers already construct
//! actual closed analytic B-reps. This adapter certifies original A/B face
//! lineage, multiple source face-pair circle graphs, two independent pcurves
//! on each physical circular rim, and a true same-face periodic seam.
//! No tessellation is promoted to topology and no G5 document edit is made.
//!
//! This IS NOT a general multi-cylinder / partially intersecting surface
//! Boolean: unsupported geometry returns Held.

use crate::{
    BooleanError, BooleanFaceSource, BooleanKernel, BooleanOperand, BooleanOperation,
    BooleanOutcome, BrepBody, BrepModel, BrepOrigin, CoedgeId, Curve2, Curve3,
    EdgeId, FaceId, FaceSense, Frame3, GeometryError, GeometryTolerance,
    GraphNodeKind, GraphOperand, GraphState, IntersectionGraph,
    IntersectionGraphBuilder, Point3, Surface3, TrimContact,
};
use std::f64::consts::{PI,TAU};
use std::collections::BTreeSet;

fn invalid(message:impl Into<String>)->GeometryError {
    GeometryError::InvalidTopology(message.into())
}
fn distance(a:Point3,b:Point3)->f64{
    let q=a.sub(b);
    q.x.hypot(q.y).hypot(q.z)
}
fn near(a:f64,b:f64,eps:f64)->bool{
    a.is_finite() && b.is_finite() && (a-b).abs()<=eps
}
fn precision(a:Point3,b:Point3,c:Point3,tolerance:GeometryTolerance)->f64{
    let scale=[a,b,c].into_iter()
        .flat_map(|q|[q.x.abs(),q.y.abs(),q.z.abs()])
        .fold(1.0_f64,f64::max);
    tolerance.absolute_length.max(64.0*f64::EPSILON*scale)
}
fn merge_tolerance(a:GeometryTolerance,b:GeometryTolerance)->GeometryTolerance{
    GeometryTolerance{
        absolute_length:a.absolute_length.max(b.absolute_length),
        relative_length:a.relative_length.max(b.relative_length),
        angular:a.angular.max(b.angular),
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticBoreKind { Through, Blind }
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticBoreFaceUse {
    OriginalUncut,
    OriginalWithCircularOpening,
    BoreWallInverted,
    BoreFloorFromCutterCap,
    DiscardedSourceCap,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct AnalyticBoreFaceCoverage{
    pub operand:GraphOperand,
    pub source_face:FaceId,
    pub result_face:Option<FaceId>,
    pub use_kind:AnalyticBoreFaceUse,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum AnalyticBoreHold{
    UnsupportedOperand,
    UnsupportedOperation,
    UnsupportedPartialIntersection,
    IncompleteClosedIntersectionGraph,
    MissingIndependentPcurve,
    IncorrectCarrierOrTopology,
    IncorrectSourceLineage,
    IncorrectAnalyticMass,
    TolerancePolicyMismatch,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct AnalyticRingWitness{
    pub result_edge:EdgeId,
    pub result_coedges:[CoedgeId;2],
    pub planar_result_face:FaceId,
    pub cylindrical_result_face:FaceId,
    /// Through: both box bottom/top planes. Blind: only the opening uses
    /// the box top plane; the bottom circle is on the source cylinder cap.
    pub source_box_face:Option<FaceId>,
    pub source_cylinder_face:FaceId,
    /// An original closed graph arc has one synthetic *chart* anchor,
    /// which is NOT a physical source B-rep vertex identity.
    pub closed_graph_witness:bool,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct AnalyticSeamWitness {
    pub source_cylinder_face:FaceId,
    pub result_edge:EdgeId,
    pub uses:[CoedgeId;2],
    /// Separate, unwrapped U coordinates on the same periodic face.
    pub uv_u:[f64;2],
}
/// Two different source-face UV charts for the SAME closed 3D graph arc.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct AnalyticSourcePcurvePair {
    pub source_box_face:FaceId,
    pub source_cylinder_face:FaceId,
    pub planar:Curve2,
    pub cylinder:Curve2,
    pub graph_interval:crate::ParameterRange,
    /// The source graph's single periodic anchor is not a material vertex.
    pub virtual_periodic_anchor:bool,
}
#[derive(Debug,Clone)]
pub struct AnalyticBoreCandidate {
    model:BrepModel,
    operation:BooleanOperation,
    kind:AnalyticBoreKind,
    lineage:Vec<AnalyticBoreFaceCoverage>,
    source_faces:Vec<BooleanFaceSource>,
    rings:[AnalyticRingWitness;2],
    source_circle_graphs:Vec<IntersectionGraph>,
    source_pcurves:Vec<AnalyticSourcePcurvePair>,
    seam:AnalyticSeamWitness,
    source_revision:u64,
    input_cylinder_frame:Frame3,
    input_cylinder_radius:f64,
    input_cylinder_height:f64,
    construction_tolerance:GeometryTolerance,
    can_apply:bool,
}
#[derive(Debug,Clone)]
pub enum AnalyticBoreCandidateResult{
    Ready(AnalyticBoreCandidate),
    Held(AnalyticBoreHold),
}

impl AnalyticBoreCandidate {
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn kind(&self)->AnalyticBoreKind{self.kind}
    pub fn face_coverage(&self)->&[AnalyticBoreFaceCoverage]{&self.lineage}
    pub fn face_provenance(&self)->&[BooleanFaceSource]{&self.source_faces}
    pub fn circle_graphs(&self)->&[IntersectionGraph]{&self.source_circle_graphs}
    pub fn source_pcurves(&self)->&[AnalyticSourcePcurvePair]{&self.source_pcurves}
    pub fn rings(&self)->&[AnalyticRingWitness;2]{&self.rings}
    pub fn seam(&self)->AnalyticSeamWitness{self.seam}
    pub fn source_revision(&self)->u64{self.source_revision}
    pub fn can_apply(&self)->bool{self.can_apply}
    pub fn closed_shell_proven(&self)->bool{true}

    pub fn validate_against(&self,left:&BrepBody,right:&BrepBody,
        tolerance:GeometryTolerance)->Result<(),GeometryError>
    {
        let (BrepBody::Polyhedral(a),BrepBody::Cylinder(b))=(left,right) else{
            return Err(invalid("analytic bore source type changed"));
        };
        if a.revision!=self.source_revision ||
            b.frame!=self.input_cylinder_frame ||
            b.radius!=self.input_cylinder_radius ||
            b.height!=self.input_cylinder_height ||
            self.construction_tolerance!=tolerance
        {return Err(invalid("analytic bore sources, revision or tolerance changed"));}
        let rebuilt=Self::build(left,right,BooleanOperation::Difference,tolerance)?;
        let AnalyticBoreCandidateResult::Ready(expected)=rebuilt else{
            return Err(invalid("source no longer produces a supported analytic bore"));
        };
        if self.kind!=expected.kind ||
            self.operation!=expected.operation ||
            self.lineage!=expected.lineage ||
            self.source_faces!=expected.source_faces ||
            self.rings!=expected.rings ||
            self.seam!=expected.seam ||
            self.source_circle_graphs!=expected.source_circle_graphs ||
            self.source_pcurves!=expected.source_pcurves ||
            format!("{:?}",self.model)!=format!("{:?}",expected.model)
        {return Err(invalid("analytic bore carrier/topology differs from recomputed source"));}
        self.model.validate()
    }

    pub fn build(left:&BrepBody,right:&BrepBody,operation:BooleanOperation,
        tolerance:GeometryTolerance)->Result<AnalyticBoreCandidateResult,GeometryError>
    {
        tolerance.validate()?;left.validate()?;right.validate()?;
        if operation!=BooleanOperation::Difference {
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::UnsupportedOperation));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cylinder))=(left,right) else{
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::UnsupportedOperand));
        };
        if merge_tolerance(block.tolerance,cylinder.tolerance)!=tolerance{
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::TolerancePolicyMismatch));
        }
        // This executor is deliberately restricted to actual AABB +Y bores
        // with a whole through cut or a strictly interior blind floor.
        match crate::boolean::certify_box(block){
            Ok(())=>{},
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::UnsupportedOperand)),
        }
        let result=match BooleanKernel::apply(left,right,operation){
            Ok(s)=>s,
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::UnsupportedPartialIntersection)),
        };
        let (kind,model)=match result.outcome {
            BooleanOutcome::Analytic(a)=>{
                a.validate()?;
                (AnalyticBoreKind::Through,a.model().clone())
            },
            BooleanOutcome::BlindHole(a)=>{
                a.validate()?;
                (AnalyticBoreKind::Blind,a.model().clone())
            },
            _=>return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::UnsupportedPartialIntersection)),
        };
        model.validate()?;
        if model.origin!=match kind{
            AnalyticBoreKind::Through=>BrepOrigin::AnalyticThroughHole,
            AnalyticBoreKind::Blind=>BrepOrigin::AnalyticBlindHole,
        }{
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::IncorrectCarrierOrTopology));
        }
        let source_a=left.shared()?;
        let source_b=right.shared()?;
        let mut closed_graphs=Vec::<IntersectionGraph>::new();
        let mut closed_pcurves=Vec::<AnalyticSourcePcurvePair>::new();
        // The *two* face pairs for a through-hole are independently built:
        // bottom box/cylinder-wall and top box/cylinder-wall.
        // A blind hole has only the top opening; the floor comes directly
        // from the cylinder's lower source cap.
        let cap_faces=match kind {
            AnalyticBoreKind::Through=>vec![FaceId(0),FaceId(1)],
            AnalyticBoreKind::Blind=>vec![FaceId(1)],
        };
        for face in cap_faces {
            let graph=IntersectionGraphBuilder::new(tolerance)?
                .intersect_faces(&source_a,face,&source_b,FaceId(2))?;
            graph.validate_against(&source_a,face,&source_b,FaceId(2))?;
            if graph.state!=GraphState::Ready || graph.nodes.len()!=1 ||
               graph.arcs.len()!=1 || graph.nodes[0].kind!=GraphNodeKind::PeriodicAnchor ||
               graph.arcs[0].from!=graph.arcs[0].to ||
               !graph.arcs[0].closed ||
               graph.arcs[0].contact==TrimContact::Tangent ||
               !near(graph.arcs[0].parameter_end-graph.arcs[0].parameter_start,
                     TAU,tolerance.angular*16.0)
            {
                return Ok(AnalyticBoreCandidateResult::Held(
                    AnalyticBoreHold::IncompleteClosedIntersectionGraph));
            }
            // The source plane and source cylinder have different UV
            // parameter charts. Reconstruct BOTH source-owned analytic
            // curves on their original support surfaces, at the original
            // closed-arc parameter values.
            let arc=graph.arcs[0];
            let [on_plane,on_cylinder]=arc.pcurves(&source_a,&source_b)?;
            let [a_curve,b_curve]=[on_plane.to_curve2()?,on_cylinder.to_curve2()?];
            let Some(plane)=source_a.face_surface(face) else {
                return Err(invalid("analytic source cap face missing"));
            };
            let Some(side)=source_b.face_surface(FaceId(2)) else {
                return Err(invalid("analytic source cutter wall missing"));
            };
            let eps=precision(source_a.bbox.min,source_a.bbox.max,
                cylinder.frame.origin,tolerance)*16.0;
            for fraction in [0.0_f64,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0] {
                let parameter=arc.parameter_start+
                    (arc.parameter_end-arc.parameter_start)*fraction;
                let xyz=arc.evaluate(fraction)?;
                let uv0=a_curve.evaluate(parameter)?;
                let uv1=b_curve.evaluate(parameter)?;
                if distance(plane.evaluate(uv0.u,uv0.v)?.point,xyz)>eps ||
                   distance(side.evaluate(uv1.u,uv1.v)?.point,xyz)>eps {
                    return Ok(AnalyticBoreCandidateResult::Held(
                        AnalyticBoreHold::MissingIndependentPcurve));
                }
            }
            closed_pcurves.push(AnalyticSourcePcurvePair{
                source_box_face:face,source_cylinder_face:FaceId(2),
                planar:a_curve,cylinder:b_curve,
                graph_interval:crate::ParameterRange::bounded(
                    arc.parameter_start,arc.parameter_end)?,
                virtual_periodic_anchor:true,
            });
            closed_graphs.push(graph);
        }
        let expected_faces=if kind==AnalyticBoreKind::Through{7}else{8};
        let expected_genus=if kind==AnalyticBoreKind::Through{1}else{0};
        let sum=model.summary();
        if model.faces.len()!=expected_faces || model.edges.len()!=15 ||
            model.coedges.len()!=30 || model.loops.len()!=9 ||
            model.vertices.len()!=10 || model.shells.len()!=1 ||
            !model.shells[0].closed ||
            sum.circular_edges!=2 || sum.closed_edges!=2 ||
            sum.seam_edges!=1 || sum.cylindrical_faces!=1 ||
            sum.genus!=expected_genus ||
            sum.euler_characteristic!=2-2*(expected_genus as isize)
        {
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::IncorrectCarrierOrTopology));
        }
        // Compare original 12 physical box edges, 24 coedges and six outer
        // plane carriers against their actual source, not its loose bbox.
        for i in 0..12 {
            let edge=EdgeId(i as u32);
            let original=&source_a.edges[i];let generated=&model.edges[i];
            if original.id!=generated.id || original.start!=generated.start ||
                original.end!=generated.end ||
                original.coedges!=generated.coedges ||
                source_a.edge_curve(edge)!=model.edge_curve(edge)
            {
                return Ok(AnalyticBoreCandidateResult::Held(
                    AnalyticBoreHold::IncorrectSourceLineage));
            }
        }
        for i in 0..24 {
            let old=&source_a.coedges[i];
            let new=&model.coedges[i];
            let key=CoedgeId(i as u32);
            if old.id!=new.id || old.edge!=new.edge ||
                old.face!=new.face || old.loop_id!=new.loop_id ||
                old.twin!=new.twin || old.reversed!=new.reversed ||
                source_a.coedge_pcurve(key)!=model.coedge_pcurve(key)
            {
                return Ok(AnalyticBoreCandidateResult::Held(
                    AnalyticBoreHold::IncorrectSourceLineage));
            }
        }
        for i in 0..6{
            if source_a.face_surface(FaceId(i))!=model.face_surface(FaceId(i)){
                return Ok(AnalyticBoreCandidateResult::Held(
                    AnalyticBoreHold::IncorrectSourceLineage));
            }
        }
        let lineage=match certify_lineage(&model,&result.face_sources,kind){
            Some(rows)=>rows,
            None=>return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::IncorrectSourceLineage)),
        };
        let lower_y=match kind{
            AnalyticBoreKind::Through=>block.bbox.min.y,
            AnalyticBoreKind::Blind=>cylinder.frame.origin.y,
        };
        let upper_y=block.bbox.max.y;
        let lower_face=if kind==AnalyticBoreKind::Through{
            FaceId(0)
        }else{FaceId(7)};
        let ring_spec=[
            (EdgeId(12),lower_face,lower_y,
             if kind==AnalyticBoreKind::Through{Some(FaceId(0))}else{None}),
            (EdgeId(13),FaceId(1),upper_y,Some(FaceId(1))),
        ];
        let mut ring_vec=Vec::new();
        for (edge,face,y,original_box_face) in ring_spec {
            let graph=original_box_face.and_then(|f|
                closed_graphs.iter().find(|g|g.source_faces[0].face==f));
            let witness=match certify_ring(&model,edge,face,y,cylinder,graph,tolerance)?{
                Some(w)=>w,
                None=>return Ok(AnalyticBoreCandidateResult::Held(
                    AnalyticBoreHold::MissingIndependentPcurve)),
            };
            ring_vec.push(AnalyticRingWitness{
                source_box_face:original_box_face,
                source_cylinder_face:FaceId(if original_box_face.is_some(){2}else{0}),
                ..witness
            });
        }
        let rings=[ring_vec[0],ring_vec[1]];
        let seam=match certify_periodic_seam(&model,cylinder,tolerance)?{
            Some(proof)=>proof,
            None=>return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::MissingIndependentPcurve)),
        };
        if !certify_analytic_mass(&model,block,cylinder,kind,tolerance){
            return Ok(AnalyticBoreCandidateResult::Held(
                AnalyticBoreHold::IncorrectAnalyticMass));
        }
        Ok(AnalyticBoreCandidateResult::Ready(Self{
            model,operation,kind,lineage,source_faces:result.face_sources,
            rings,source_circle_graphs:closed_graphs,
            source_pcurves:closed_pcurves,seam,
            source_revision:block.revision,
            input_cylinder_frame:cylinder.frame,
            input_cylinder_radius:cylinder.radius,
            input_cylinder_height:cylinder.height,
            construction_tolerance:tolerance,can_apply:false,
        }))
    }
}

fn certify_lineage(
    m:&BrepModel,face_sources:&[BooleanFaceSource],kind:AnalyticBoreKind,
)->Option<Vec<AnalyticBoreFaceCoverage>>{
    if face_sources.len()!=m.faces.len(){return None;}
    let mut seen=BTreeSet::new();
    for source in face_sources{
        if !seen.insert(source.result_face){return None;}
        if source.result_face.0<6 {
            if source.operand!=BooleanOperand::Left ||
               source.result_face!=source.source_face ||source.reversed {return None;}
        }else if source.result_face==FaceId(6){
            if source.operand!=BooleanOperand::Right ||
               source.source_face!=FaceId(2) || !source.reversed {return None;}
        }else if kind==AnalyticBoreKind::Blind && source.result_face==FaceId(7){
            if source.operand!=BooleanOperand::Right ||
               source.source_face!=FaceId(0) ||!source.reversed {return None;}
        }else{return None;}
    }
    if seen.len()!=m.faces.len() || (0..m.faces.len()).any(|i|!seen.contains(&FaceId(i as u32))) {
        return None;
    }
    let mut out=Vec::with_capacity(9);
    for i in 0..6{
        let face=FaceId(i as u32);
        let cut=i==1 || (i==0 && kind==AnalyticBoreKind::Through);
        let mface=m.faces.get(i)?;
        if mface.id!=face || mface.loops.len()!=if cut{2}else{1}{
            return None;
        }
        out.push(AnalyticBoreFaceCoverage{
            operand:GraphOperand::First,source_face:face,
            result_face:Some(face),
            use_kind:if cut{AnalyticBoreFaceUse::OriginalWithCircularOpening}
                else{AnalyticBoreFaceUse::OriginalUncut},
        });
    }
    for i in 0..3{
        let face=FaceId(i as u32);
        let use_kind=match i{
            0 if kind==AnalyticBoreKind::Blind=>AnalyticBoreFaceUse::BoreFloorFromCutterCap,
            2=>AnalyticBoreFaceUse::BoreWallInverted,
            _=>AnalyticBoreFaceUse::DiscardedSourceCap,
        };
        out.push(AnalyticBoreFaceCoverage{
            operand:GraphOperand::Second,source_face:face,
            result_face:match use_kind{
                AnalyticBoreFaceUse::BoreFloorFromCutterCap=>Some(FaceId(7)),
                AnalyticBoreFaceUse::BoreWallInverted=>Some(FaceId(6)),
                _=>None,
            },
            use_kind,
        });
    }
    if m.face_sense(FaceId(6))!=Some(FaceSense::Reversed) ||
        (kind==AnalyticBoreKind::Blind && m.face_sense(FaceId(7))!=Some(FaceSense::Reversed))
    {
        return None;
    }
    Some(out)
}

fn certify_ring(
    m:&BrepModel,edge_id:EdgeId,planar_face:FaceId,y:f64,
    cylinder:&crate::CylindricalBrep,
    graph:Option<&IntersectionGraph>,tolerance:GeometryTolerance,
)->Result<Option<AnalyticRingWitness>,GeometryError>{
    let Some(edge)=m.edges.get(edge_id.0 as usize).filter(|edge|
        edge.id==edge_id && edge.start==edge.end) else{return Ok(None);};
    let Some(Curve3::Circle(circle))=m.edge_curve(edge_id) else{return Ok(None);};
    let coedges=edge.coedges;
    if coedges[0]==coedges[1] {return Ok(None);}
    let Some(a)=m.coedges.get(coedges[0].0 as usize) else{return Ok(None);};
    let Some(b)=m.coedges.get(coedges[1].0 as usize) else{return Ok(None);};
    let planar=if a.face==planar_face{a}else{b};
    let lateral=if a.face==FaceId(6){a}else{b};
    if planar.face!=planar_face || lateral.face!=FaceId(6) ||
       planar.twin!=lateral.id || lateral.twin!=planar.id ||
       planar.edge!=edge_id || lateral.edge!=edge_id ||
       planar.reversed==lateral.reversed ||
       !matches!(m.face_surface(planar_face),Some(Surface3::Plane(_))) ||
       !matches!(m.face_surface(FaceId(6)),Some(Surface3::Cylinder(_))) ||
       !matches!(m.coedge_pcurve(planar.id),Some(Curve2::Circle(_))) ||
       !matches!(m.coedge_pcurve(lateral.id),Some(Curve2::Line(_)))
    {return Ok(None);}
    let q=cylinder.frame.origin;
    let eps=precision(q,m.bbox.min,m.bbox.max,tolerance)*16.0;
    for f in [0.0_f64,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0]{
        let sample=circle.evaluate(f*TAU)?;
        let rho=(sample.x-q.x).hypot(sample.z-q.z);
        if !near(sample.y,y,eps) || !near(rho,cylinder.radius,eps) {
            return Ok(None);
        }
        for co in [planar,lateral]{
            let Some(trim)=m.coedge_pcurve(co.id)else{return Ok(None);};
            let data=&m.geometry.coedges[co.id.0 as usize];
            let domain=match data.domain{
                crate::ParameterRange::Bounded{start,end}=>(start,end),
                _=>return Ok(None),
            };
            let fraction=if co.reversed{1.0-f}else{f};
            let uv=trim.evaluate(domain.0+(domain.1-domain.0)*fraction)?;
            let Some(surface)=m.face_surface(co.face) else{return Ok(None);};
            let surface_sample=surface.evaluate(uv.u,uv.v)?.point;
            if distance(sample,surface_sample)>eps {return Ok(None);}
        }
        if let Some(g)=graph{
            let Some(arc)=g.arcs.first()else{return Ok(None);};
            let radial=sample.sub(q);
            let angle=radial.dot(cylinder.frame.y)
                .atan2(radial.dot(cylinder.frame.x)).rem_euclid(TAU);
            let param=(angle-arc.parameter_start).rem_euclid(TAU)+arc.parameter_start;
            if param>arc.parameter_end+eps ||
                distance(arc.curve.evaluate(param)?,sample)>eps
            {return Ok(None);}
        }
    }
    Ok(Some(AnalyticRingWitness{
        result_edge:edge_id,result_coedges:coedges,
        planar_result_face:planar_face,
        cylindrical_result_face:FaceId(6),
        source_box_face:None,source_cylinder_face:FaceId(2),
        closed_graph_witness:graph.is_some(),
    }))
}
fn certify_periodic_seam(
    m:&BrepModel,cylinder:&crate::CylindricalBrep,tolerance:GeometryTolerance,
)->Result<Option<AnalyticSeamWitness>,GeometryError>{
    let edge_id=EdgeId(14);
    let Some(edge)=m.edges.get(14)else{return Ok(None);};
    let [a,b]=edge.coedges;
    if a==b||edge.start==edge.end{return Ok(None);}
    let Some(x)=m.coedges.get(a.0 as usize)else{return Ok(None);};
    let Some(y)=m.coedges.get(b.0 as usize)else{return Ok(None);};
    if x.face!=FaceId(6)||y.face!=FaceId(6)||
        x.twin!=b||y.twin!=a ||x.reversed==y.reversed ||
        !matches!(m.edge_curve(edge_id),Some(Curve3::Line(_)))
    {return Ok(None);}
    let Some(Curve2::Line(c0))=m.coedge_pcurve(a)else{return Ok(None);};
    let Some(Curve2::Line(c1))=m.coedge_pcurve(b)else{return Ok(None);};
    let Some(surface)=m.face_surface(FaceId(6))else{return Ok(None);};
    let Some(Curve3::Line(line))=m.edge_curve(edge_id) else{return Ok(None);};
    let u0=c0.evaluate(0.5)?.u;
    let u1=c1.evaluate(0.5)?.u;
    if !near((u0-u1).abs(),TAU,tolerance.angular*16.0) {
        return Ok(None);
    }
    let eps=precision(cylinder.frame.origin,m.bbox.min,m.bbox.max,tolerance)*16.0;
    for t in [0.0_f64,0.25,0.5,0.75,1.0]{
        let world=line.evaluate(t*cylinder.height)?;
        let uv0=c0.evaluate(if x.reversed{1.0-t}else{t})?;
        let uv1=c1.evaluate(if y.reversed{1.0-t}else{t})?;
        if distance(surface.evaluate(uv0.u,uv0.v)?.point,world)>eps ||
           distance(surface.evaluate(uv1.u,uv1.v)?.point,world)>eps
        {return Ok(None);}
    }
    Ok(Some(AnalyticSeamWitness {
        source_cylinder_face:FaceId(2),result_edge:edge_id,
        uses:[a,b],uv_u:[u0,u1],
    }))
}
fn certify_analytic_mass(m:&BrepModel,block:&crate::Solid,
    cylinder:&crate::CylindricalBrep,kind:AnalyticBoreKind,
    tolerance:GeometryTolerance)->bool
{
    let bb=block.bbox;
    let dx=bb.max.x-bb.min.x;
    let dy=bb.max.y-bb.min.y;
    let dz=bb.max.z-bb.min.z;
    let box_volume=dx*dy*dz;
    let bore_height=if kind==AnalyticBoreKind::Through{dy}else{cylinder.height};
    let removed=PI*cylinder.radius*cylinder.radius*bore_height;
    let remaining=box_volume-removed;
    let area=2.0*(dx*dy+dy*dz+dx*dz)+
        TAU*cylinder.radius*bore_height-
        if kind==AnalyticBoreKind::Through{2.0*PI*cylinder.radius.powi(2)}else{0.0};
    let box_center=Point3{x:(bb.min.x+bb.max.x)*0.5,
        y:(bb.min.y+bb.max.y)*0.5,z:(bb.min.z+bb.max.z)*0.5};
    let bore_center=Point3{x:cylinder.frame.origin.x,
        y:cylinder.frame.origin.y+0.5*bore_height,
        z:cylinder.frame.origin.z};
    let centroid=Point3{
        x:(box_center.x*box_volume-bore_center.x*removed)/remaining,
        y:(box_center.y*box_volume-bore_center.y*removed)/remaining,
        z:(box_center.z*box_volume-bore_center.z*removed)/remaining,
    };
    let eps=precision(bb.min,bb.max,cylinder.frame.origin,tolerance)*16.0;
    remaining.is_finite()&&remaining>0.0&&area.is_finite()&&area>0.0&&
        near(m.mass.volume,remaining,remaining*1e-10+eps)&&
        near(m.mass.surface_area,area,area*1e-10+eps)&&
        distance(m.mass.centroid,centroid)<=eps
}
