//! Read-only, source-bound Boolean *region* witnesses for curved intersections.
//!
//! Point-in-volume is independent of face clipping. For validated polyhedral
//! solids, use triangle solid angles (not a bbox as a geometric stand-in).
//! For the upright or rotated analytic cylinder use its exact radial and
//! axial inequalities. Graph decisions sample both sides of every crossing
//! trace *within each source face*, then classify against the opposite volume.
//! Nothing here cuts, mutates, or authorizes applying TopologyDeltaPlan.
use crate::{
    BrepBody, BrepModel, FaceId, GeometryError, GeometryTolerance, GraphArcId,
    GraphOperand, GraphState, IntersectionGraph, Point3, Solid, Surface3,
    TrimContact, TrimPosition, BooleanOperation,
};
use crate::trimmed_intersections::classify_trim_point;
use std::f64::consts::PI;

fn invalid(msg:impl Into<String>) -> GeometryError { GeometryError::InvalidTopology(msg.into()) }
fn length(p:Point3)->f64 { p.x.hypot(p.y).hypot(p.z) }
fn normalize(p:Point3)->Result<Point3,GeometryError>{
    let l=length(p);
    if !l.is_finite() || l<=0.0 {return Err(invalid("degenerate Boolean region probe direction"));}
    Ok(p.scale(1.0/l))
}
fn close(a:f64,b:f64,eps:f64)->bool{(a-b).abs()<=eps}
fn position_eps(body:&BrepBody)->f64{
    let (tol,extent)=match body {
        BrepBody::Polyhedral(s)=>(s.tolerance,
            length(s.bbox.max.sub(s.bbox.min))),
        BrepBody::Cylinder(c)=>(c.tolerance,c.radius.max(c.height)),
    };
    tol.length_at(extent)
}
fn position_scale(point:Point3)->f64{
    point.x.abs().max(point.y.abs()).max(point.z.abs()).max(1.0)
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum VolumePosition { Inside, Outside, Boundary }
impl VolumePosition {
    pub fn is_strict(self)->bool { self!=Self::Boundary }
}

/// A certified input volume is required; boundary points never silently
/// become 'inside' or 'outside' through an epsilon convention.
#[derive(Debug,Clone,Copy,Default)]
pub struct BooleanRegionClassifier;
impl BooleanRegionClassifier {
    pub fn classify_point(body:&BrepBody,point:Point3)->Result<VolumePosition,GeometryError>{
        body.validate()?;
        Self::classify_validated(body,point)
    }
    fn classify_validated(body:&BrepBody,point:Point3)
        ->Result<VolumePosition,GeometryError>
    {
        if !point.is_finite(){return Err(invalid("nonfinite Boolean membership query"));}
        let eps=position_eps(body).max(8.0*f64::EPSILON*position_scale(point));
        match body {
            BrepBody::Polyhedral(s)=>classify_polyhedron(s,point,eps),
            BrepBody::Cylinder(c)=>{
                let offset=point.sub(c.frame.origin);
                let along=offset.dot(c.frame.normal);
                let x=offset.dot(c.frame.x);
                let y=offset.dot(c.frame.y);
                let radius=x.hypot(y);
                if !along.is_finite() || !radius.is_finite(){
                    return Err(invalid("analytic cylinder membership overflow"));
                }
                if along < -eps || along > c.height+eps || radius > c.radius+eps {
                    return Ok(VolumePosition::Outside);
                }
                if close(along,0.0,eps) || close(along,c.height,eps) ||
                   close(radius,c.radius,eps) {
                    return Ok(VolumePosition::Boundary);
                }
                Ok(VolumePosition::Inside)
            }
        }
    }

    /// Produces immutable, source-revalidated evidence for a selected
    /// two-face IntersectionGraph. A classified crossing is a **preview**:
    /// full face partitioning, multiple traces, curved trims, and region
    /// reconstruction are separate topology-editing gates.
    pub fn classify_graph(graph:&IntersectionGraph,
        left:&BrepBody,left_face:FaceId,
        right:&BrepBody,right_face:FaceId,
        operation:BooleanOperation)->Result<GraphRegionPlan,GeometryError>
    {
        let lm=left.shared()?;
        let rm=right.shared()?;
        graph.validate_against(&lm,left_face,&rm,right_face)?;
        let mut result=GraphRegionPlan{
            operation,state:GraphRegionState::NoCut,hold:None,
            source_faces:[left_face,right_face],samples:Vec::new(),
        };
        match graph.state {
            GraphState::Empty=>return Ok(result),
            GraphState::Ambiguous(_)=>{
                result.state=GraphRegionState::Held;
                result.hold=Some(RegionHold::AmbiguousIntersection);
                return Ok(result);
            }
            GraphState::Ready=>{},
        }
        if graph.arcs.is_empty(){return Ok(result);}
        let bodies=[left,right];
        let models=[&lm,&rm];
        for arc in &graph.arcs {
            if arc.contact!=TrimContact::Crossing {
                result.state=GraphRegionState::Held;
                result.hold=Some(match arc.contact {
                    TrimContact::Tangent=>RegionHold::Tangency,
                    TrimContact::BoundaryOverlap=>RegionHold::BoundaryOverlap,
                    TrimContact::Crossing=>unreachable!(),
                });
                result.samples.clear();
                return Ok(result);
            }
            for sample_index in 0..3 {
                let fraction=[0.2,0.5,0.8][sample_index];
                let point=arc.evaluate(fraction)?;
                let parameter=arc.parameter_start+
                    (arc.parameter_end-arc.parameter_start)*fraction;
                let tangent=normalize(arc.curve.derivative(parameter)?)?;
                for source_index in 0..2 {
                    let opposite_index=1-source_index;
                    let source_face=graph.source_faces[source_index].face;
                    let opposite_face=graph.source_faces[opposite_index].face;
                    let model=models[source_index];
                    let own_surface=model.face_surface(source_face)
                        .ok_or_else(||invalid("missing source surface for region classification"))?;
                    let other_surface=models[opposite_index].face_surface(opposite_face)
                        .ok_or_else(||invalid("missing opposite surface for region classification"))?;
                    let normal=outward_surface_normal(own_surface,point)?;
                    let other_normal=outward_surface_normal(other_surface,point)?;
                    let across=normalize(normal.cross(tangent))?;
                    // If a tiny probe in the source's tangent plane does not
                    // traverse the other surface transversally, mark held.
                    if across.dot(other_normal).abs() <= 1e-3 {
                        return Ok(GraphRegionPlan::held(operation,[left_face,right_face],
                            RegionHold::NearlyTangent));
                    }
                    let body=bodies[source_index];
                    let other=bodies[opposite_index];
                    let diag=length(model.bbox.max.sub(model.bbox.min));
                    let other_diag=length(models[opposite_index].bbox.max.sub(models[opposite_index].bbox.min));
                    let feature=diag.min(other_diag);
                    let precision=position_eps(body).max(position_eps(other))
                        .max(graph.tolerance.length)
                        .max(graph.tolerance.evaluation_budget);
                    let base=(feature*1e-4).max(precision*128.0);
                    if !base.is_finite() || base<=0.0 {
                        return Err(invalid("invalid Boolean sampling scale"));
                    }
                    let mut chosen=None;
                    // Re-run at increasingly large scales if the probe fell
                    // within an ambiguous 'on' band; never fuse coordinates.
                    for factor in [1.0,2.0,4.0,8.0] {
                        let distance=base*factor;
                        if distance>=feature*0.03 {break;}
                        let positive=project_to_surface(
                            own_surface,point.add(across.scale(distance)))?;
                        let negative=project_to_surface(
                            own_surface,point.sub(across.scale(distance)))?;
                        if classify_trim_point(model,source_face,positive,body_tolerance(body))?
                            !=TrimPosition::Interior ||
                           classify_trim_point(model,source_face,negative,body_tolerance(body))?
                            !=TrimPosition::Interior {
                            break;
                        }
                        let pos=Self::classify_validated(other,positive)?;
                        let neg=Self::classify_validated(other,negative)?;
                        if pos==VolumePosition::Boundary || neg==VolumePosition::Boundary {
                            continue;
                        }
                        if pos==neg {
                            return Ok(GraphRegionPlan::held(operation,[left_face,right_face],
                                RegionHold::NonSeparatingTrace));
                        }
                        chosen=Some((positive,negative,pos,neg));
                        break;
                    }
                    let Some((positive,negative,pos,neg))=chosen else {
                        return Ok(GraphRegionPlan::held(operation,[left_face,right_face],
                            RegionHold::UnresolvedProbe));
                    };
                    let operand=if source_index==0 {GraphOperand::First}else{GraphOperand::Second};
                    let keep=|status|keep_face_piece(operation,operand,status);
                    result.samples.push(ArcRegionSample{
                        arc:arc.id,operand,face:source_face,
                        fraction,point,positive,negative,
                        positive_other:pos,negative_other:neg,
                        keep_positive:keep(pos),keep_negative:keep(neg),
                        reverse_orientation:operation==BooleanOperation::Difference &&
                            operand==GraphOperand::Second,
                    });
                }
            }
        }
        result.state=GraphRegionState::ClassifiedPreview;
        Ok(result)
    }
}
fn body_tolerance(body:&BrepBody)->GeometryTolerance{
    match body {BrepBody::Polyhedral(s)=>s.tolerance,BrepBody::Cylinder(c)=>c.tolerance}
}
fn keep_face_piece(op:BooleanOperation,source:GraphOperand,other:VolumePosition)->bool{
    match (op,source,other) {
        (BooleanOperation::Union,_,VolumePosition::Outside)=>true,
        (BooleanOperation::Intersection,_,VolumePosition::Inside)=>true,
        (BooleanOperation::Difference,GraphOperand::First,VolumePosition::Outside)=>true,
        (BooleanOperation::Difference,GraphOperand::Second,VolumePosition::Inside)=>true,
        _=>false,
    }
}
fn outward_surface_normal(surface:Surface3,p:Point3)->Result<Point3,GeometryError>{
    match surface {
        Surface3::Plane(plane)=>normalize(plane.frame.normal),
        Surface3::Cylinder(c)=>{
            let radial=p.sub(c.frame.origin);
            let axial=radial.dot(c.frame.normal);
            normalize(radial.sub(c.frame.normal.scale(axial)))
        },
        Surface3::Sphere(_)=>Err(invalid("sphere region probes are not supported")),
    }
}
fn project_to_surface(surface:Surface3,p:Point3)->Result<Point3,GeometryError>{
    match surface {
        Surface3::Plane(plane)=>{
            let n=plane.frame.normal;
            Ok(p.sub(n.scale(p.sub(plane.frame.origin).dot(n))))
        },
        Surface3::Cylinder(c)=>{
            let q=p.sub(c.frame.origin);
            let axial=q.dot(c.frame.normal);
            let radial=q.sub(c.frame.normal.scale(axial));
            Ok(c.frame.origin.add(c.frame.normal.scale(axial))
                .add(normalize(radial)?.scale(c.radius)))
        },
        Surface3::Sphere(_)=>Err(invalid("sphere region probes are not supported")),
    }
}
fn classify_polyhedron(s:&Solid,p:Point3,eps:f64)
    ->Result<VolumePosition,GeometryError>
{
    let min=s.bbox.min.array();let max=s.bbox.max.array();let pos=p.array();
    if (0..3).any(|i|pos[i]<min[i]-eps || pos[i]>max[i]+eps){
        return Ok(VolumePosition::Outside);
    }
    let mut omega=0.0;
    for triangle in &s.mesh.triangles {
        let a=s.mesh.vertices[triangle[0] as usize].sub(p);
        let b=s.mesh.vertices[triangle[1] as usize].sub(p);
        let c=s.mesh.vertices[triangle[2] as usize].sub(p);
        let edge0=b.sub(a);let edge1=c.sub(a);
        let cross=edge0.cross(edge1);
        let twice_area=length(cross);
        if !twice_area.is_finite() || twice_area<=0.0 {
            return Err(invalid("degenerate polyhedral triangle in region classifier"));
        }
        let signed_distance=-a.dot(cross)/twice_area;
        if signed_distance.abs()<=eps {
            let projected=a.scale(-1.0).add(cross.scale(signed_distance/twice_area));
            let aa=edge0.dot(edge0);
            let ab=edge0.dot(edge1);
            let bb=edge1.dot(edge1);
            let ac=projected.dot(edge0);
            let bc=projected.dot(edge1);
            let denominator=aa*bb-ab*ab;
            if !denominator.is_finite() || denominator<=0.0 {
                return Err(invalid("degenerate triangle barycentric frame"));
            }
            let u=(bb*ac-ab*bc)/denominator;
            let v=(aa*bc-ab*ac)/denominator;
            let slack=4.0*eps/
                length(edge0).min(length(edge1)).min(length(c.sub(b)));
            if u>=-slack && v>=-slack && u+v<=1.0+slack {
                return Ok(VolumePosition::Boundary);
            }
        }
        let la=length(a);let lb=length(b);let lc=length(c);
        if la<=eps || lb<=eps || lc<=eps {
            return Ok(VolumePosition::Boundary);
        }
        let numerator=a.dot(b.cross(c));
        let denominator=la*lb*lc+a.dot(b)*lc+b.dot(c)*la+c.dot(a)*lb;
        omega+=2.0*numerator.atan2(denominator);
    }
    // Solid-angle winding counts the actual closed triangle boundary and
    // remains valid for nonconvex polyhedra. These are numerical certificates,
    // not a claim of exact arithmetic.
    if !omega.is_finite(){return Err(invalid("nonfinite solid-angle sum"));}
    if (omega.abs()-4.0*PI).abs()<1e-5 {
        Ok(VolumePosition::Inside)
    } else if omega.abs()<1e-5 {
        Ok(VolumePosition::Outside)
    } else {
        Err(invalid("ambiguous polyhedral winding: cannot assign Boolean region"))
    }
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum GraphRegionState { NoCut, ClassifiedPreview, Held }
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum RegionHold {
    AmbiguousIntersection,Tangency,BoundaryOverlap,
    NearlyTangent,NonSeparatingTrace,UnresolvedProbe,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ArcRegionSample {
    pub arc:GraphArcId,
    pub operand:GraphOperand,
    pub face:FaceId,
    /// Relative position along the validated analytic intersection arc.
    pub fraction:f64,
    pub point:Point3,
    /// Two points on the source surface on opposite sides of the cut.
    pub positive:Point3,
    pub negative:Point3,
    pub positive_other:VolumePosition,
    pub negative_other:VolumePosition,
    pub keep_positive:bool,
    pub keep_negative:bool,
    pub reverse_orientation:bool,
}
#[derive(Debug,Clone,PartialEq)]
pub struct GraphRegionPlan {
    pub operation:BooleanOperation,
    pub state:GraphRegionState,
    pub hold:Option<RegionHold>,
    pub source_faces:[FaceId;2],
    pub samples:Vec<ArcRegionSample>,
}
impl GraphRegionPlan {
    fn held(operation:BooleanOperation,source_faces:[FaceId;2],why:RegionHold)->Self{
        Self {operation,state:GraphRegionState::Held,hold:Some(why),
            source_faces,samples:Vec::new()}
    }
    /// A tampered plan cannot masquerade as a source-bound region proof.
    /// The graph is rechecked against both *current* input bodies.
    pub fn validate_against(&self,graph:&IntersectionGraph,
        left:&BrepBody,left_face:FaceId,right:&BrepBody,right_face:FaceId)
        ->Result<(),GeometryError>
    {
        let regenerated=BooleanRegionClassifier::classify_graph(
            graph,left,left_face,right,right_face,self.operation)?;
        if &regenerated!=self {
            return Err(invalid("Boolean region plan differs from source-bound proof"));
        }
        Ok(())
    }
    pub fn can_apply(&self)->bool {false}
}
