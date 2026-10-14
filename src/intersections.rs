//! Analytic supporting-surface intersections with explicit uncertainty near
//! topology-changing degeneracies. Geometry is independent of render facets.
//!
//! Current support: plane/plane; infinite plane/right circular cylinder,
//! including parallel-axis secants/tangencies and oblique closed sections.
//! These are UNTRIMMED intersections and are not yet Boolean edges.
use crate::{CylinderSurface, Frame3, GeometryError, GeometryTolerance, Line3,
    ParameterRange, Plane3, PlaneIntersection, PlaneSurface, Point3, Surface3};
use std::f64::consts::TAU;

fn invalid(message: impl Into<String>) -> GeometryError {
    GeometryError::InvalidEdit(message.into())
}
fn norm(v: Point3) -> f64 { v.x.hypot(v.y).hypot(v.z) }
fn distance(a: Point3, b: Point3) -> f64 { norm(a.sub(b)) }
fn checked(p: Point3) -> Result<Point3, GeometryError> {
    if p.is_finite() {Ok(p)} else {Err(invalid("non-finite intersection construction"))}
}

/// Stable registry for the independently contracted analytical pair families.
/// New surface pairs require an executable first-party mathematical witness.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IntersectionPair { PlanePlane, PlaneCylinder }
pub const ALL_INTERSECTION_PAIRS: [IntersectionPair; 2] = [
    IntersectionPair::PlanePlane, IntersectionPair::PlaneCylinder,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectionContact { Crossing, Tangent }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectionAmbiguity { NearParallel, NearTangency }

/// An exact analytic parameterization in the mathematical model, with f64
/// evaluations. This is a support curve, not a clipped B-rep edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntersectionCurve {
    Line(Line3),
    PlaneCylinderSection(PlaneCylinderSection),
}
impl IntersectionCurve {
    pub fn evaluate(self, parameter: f64) -> Result<Point3, GeometryError> {
        match self {
            Self::Line(c) => c.evaluate(parameter),
            Self::PlaneCylinderSection(c) => c.evaluate(parameter),
        }
    }
    pub fn derivative(self, parameter: f64) -> Result<Point3, GeometryError> {
        match self {
            Self::Line(c) => {
                if !parameter.is_finite() {return Err(invalid("non-finite line parameter"));}
                Ok(c.direction)
            }
            Self::PlaneCylinderSection(c) => c.derivative(parameter),
        }
    }
    pub fn parameter_domain(self) -> ParameterRange {
        match self {
            Self::Line(_) => ParameterRange::Unbounded,
            Self::PlaneCylinderSection(_) => ParameterRange::Bounded{start:0.0,end:TAU},
        }
    }
}

/// v(u) is solved from n·(S(u,v)-plane.origin)=0, not fitted from facets.
/// Circular and tilted sections share the same analytic representation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneCylinderSection {
    cylinder: CylinderSurface,
    /// n·(plane.origin-cylinder.origin)
    projected_offset: f64,
    /// n·cylinder.axis. Must be resolvably nonzero.
    axial: f64,
    /// Dot products of the plane normal with the cylinder basis.
    coefficients: [f64;2],
}
impl PlaneCylinderSection {
    /// The cylinder chart in which this analytic section is parameterized.
    /// The section parameter is the *unwrapped* angular U of this chart.
    pub(crate) fn cylinder(self) -> CylinderSurface { self.cylinder }

    /// Coefficients v(u)=a+b cos(u)+c sin(u), valid only for the oblique
    /// (non-generator) family constructed by SurfaceIntersectionEngine.
    pub fn axial_coefficients(self)->Result<[f64;3],GeometryError>{
        if !self.axial.is_finite() || self.axial==0.0 {
            return Err(invalid("oblique section has degenerate axial coefficient"));
        }
        let k=self.cylinder.radius/self.axial;
        let values=[self.projected_offset/self.axial,
            -k*self.coefficients[0],-k*self.coefficients[1]];
        if !values.iter().all(|x|x.is_finite()) {
            return Err(invalid("oblique section axial coefficients overflow"));
        }
        Ok(values)
    }

    /// Exact-in-form spatial harmonic carrier of the plane/cylinder section.
    /// No fitting, tessellation or sampling is involved.
    pub fn harmonic_curve(self)->Result<crate::HarmonicCurve3,GeometryError>{
        let [mean,c,s]=self.axial_coefficients()?;
        let frame=self.cylinder.frame;
        let center=frame.origin.add(frame.normal.scale(mean));
        let cosine=frame.x.scale(self.cylinder.radius).add(frame.normal.scale(c));
        let sine=frame.y.scale(self.cylinder.radius).add(frame.normal.scale(s));
        if !center.is_finite() || !cosine.is_finite() || !sine.is_finite(){
            return Err(invalid("nonfinite oblique spatial curve coefficients"));
        }
        Ok(crate::HarmonicCurve3{center,cosine,sine})
    }

    pub fn axial_parameter(self, u: f64) -> Result<f64, GeometryError> {
        if !u.is_finite() {return Err(invalid("non-finite circular parameter"));}
        let (sin,cos)=u.sin_cos();
        let radial_dot=self.cylinder.radius*(self.coefficients[0]*cos+self.coefficients[1]*sin);
        let v=(self.projected_offset-radial_dot)/self.axial;
        if !v.is_finite() {return Err(invalid("plane/cylinder section parameter overflow"));}
        Ok(v)
    }
    pub fn evaluate(self,u:f64)->Result<Point3,GeometryError> {
        let v=self.axial_parameter(u)?;
        checked(Surface3::Cylinder(self.cylinder).evaluate(u,v)?.point)
    }
    pub fn derivative(self,u:f64)->Result<Point3,GeometryError> {
        if !u.is_finite(){return Err(invalid("non-finite circular parameter"));}
        let (sin,cos)=u.sin_cos();
        let dv=self.cylinder.radius*(self.coefficients[0]*sin-self.coefficients[1]*cos)/self.axial;
        let d=Surface3::Cylinder(self.cylinder).evaluate(u,self.axial_parameter(u)?)?;
        checked(d.du.add(d.dv.scale(dv)))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceIntersection {
    Empty,
    /// Infinite supporting planes coincide; overlapping trimmed regions are
    /// not calculated by this API.
    Coincident,
    /// A tiny perturbation could change the number/type of intersection
    /// branches; downstream modeling MUST NOT use this as a Boolean cut.
    Ambiguous { reason:IntersectionAmbiguity, gap:f64 },
    Curves {contact:IntersectionContact, curves:Vec<IntersectionCurve>},
}
impl SurfaceIntersection {
    pub fn curves(&self) -> &[IntersectionCurve] {
        match self {Self::Curves{curves,..}=>curves.as_slice(), _=>&[]}
    }
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Empty=>"empty", Self::Coincident=>"coincident",
            Self::Ambiguous{..}=>"ambiguous",
            Self::Curves{contact:IntersectionContact::Tangent,..}=>"tangent",
            Self::Curves{..}=>"crossing",
        }
    }
    pub fn certify(self: &Self, a:Surface3, b:Surface3, tolerance:GeometryTolerance)
        -> Result<IntersectionCertificate,GeometryError>
    {
        tolerance.validate()?;
        let Self::Curves{curves,..}=self else {
            return Err(invalid("only constructed intersection branches can be residue-certified"));
        };
        if curves.is_empty(){return Err(invalid("intersection has no branches"));}
        let mut sample_count=0;
        let mut max_residual_a: f64=0.0;
        let mut max_residual_b: f64=0.0;
        let mut max_coordinate: f64=1.0;
        for curve in curves {
            let domain=curve.parameter_domain();
            let params: Vec<f64>=match domain {
                ParameterRange::Unbounded=>vec![-3.0,-1.0,0.0,1.0,3.0],
                ParameterRange::Bounded{start,end}=>(0..=16).map(|i|start+(end-start)*i as f64/16.0).collect(),
            };
            for t in params {
                let p=curve.evaluate(t)?;
                let derivative=curve.derivative(t)?;
                if !p.is_finite() || !derivative.is_finite() || norm(derivative)==0.0 {
                    return Err(invalid("intersection evaluation or tangent is degenerate"));
                }
                max_coordinate=max_coordinate.max(p.x.abs().max(p.y.abs()).max(p.z.abs()));
                max_residual_a=max_residual_a.max(implicit_residual(a,p)?);
                max_residual_b=max_residual_b.max(implicit_residual(b,p)?);
                sample_count+=1;
            }
        }
        // Error budget includes rounding at large translated coordinates.
        // These samples are evidence, NOT a completeness/exactness proof.
        let scale=characteristic_scale(a).max(characteristic_scale(b));
        let budget=tolerance.length_at(scale).max(64.0*f64::EPSILON*max_coordinate);
        if max_residual_a>budget || max_residual_b>budget {
            return Err(invalid(format!("intersection surface residual exceeds budget: {max_residual_a:e}, {max_residual_b:e} > {budget:e}")));
        }
        Ok(IntersectionCertificate{samples:sample_count,max_residual_a,max_residual_b,budget})
    }
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct IntersectionCertificate {
    pub samples:usize,
    pub max_residual_a:f64,
    pub max_residual_b:f64,
    pub budget:f64,
}

fn characteristic_scale(s:Surface3)->f64{
    match s {
        Surface3::Plane(_)=>1.0,
        Surface3::Cylinder(c)=>c.radius,
        Surface3::Sphere(s)=>s.radius,
    }
}
fn implicit_residual(s:Surface3,p:Point3)->Result<f64,GeometryError>{
    let value=match s {
        Surface3::Plane(plane)=>plane.frame.normal.dot(p.sub(plane.frame.origin)).abs(),
        Surface3::Cylinder(cylinder)=>{
            let v=p.sub(cylinder.frame.origin);
            let axial=cylinder.frame.normal.dot(v);
            (norm(v.sub(cylinder.frame.normal.scale(axial)))-cylinder.radius).abs()
        }
        Surface3::Sphere(_)=>return Err(invalid("sphere intersection certification is not implemented")),
    };
    if value.is_finite(){Ok(value)}else{Err(invalid("surface residual is not finite"))}
}
fn valid_frame(frame:Frame3, tolerance:GeometryTolerance)->Result<(),GeometryError>{
    if !frame.origin.is_finite() || !frame.x.is_finite() || !frame.y.is_finite() || !frame.normal.is_finite(){
        return Err(invalid("non-finite surface frame"));
    }
    let eps=tolerance.angular.min(1e-5);
    for v in [frame.x,frame.y,frame.normal] {
        if (norm(v)-1.0).abs()>eps {return Err(invalid("surface frame axes must be normalized"));}
    }
    if frame.x.dot(frame.y).abs()>eps || frame.x.dot(frame.normal).abs()>eps ||
       frame.y.dot(frame.normal).abs()>eps || norm(frame.x.cross(frame.y).sub(frame.normal))>eps {
        return Err(invalid("surface frame must be right handed and orthogonal"));
    }
    Ok(())
}

/// Solves untrimmed implicit surface intersections; no arbitrary faceting or
/// fallback Boolean and no assumption that every near-degenerate input has a
/// determinable topology at the chosen tolerance.
#[derive(Debug,Clone,Copy)]
pub struct SurfaceIntersectionEngine {pub tolerance:GeometryTolerance}
impl SurfaceIntersectionEngine {
    pub fn new(tolerance:GeometryTolerance)->Result<Self,GeometryError>{
        tolerance.validate()?;
        Ok(Self{tolerance})
    }
    pub fn intersect(&self,a:Surface3,b:Surface3)->Result<SurfaceIntersection,GeometryError>{
        match (a,b) {
            (Surface3::Plane(first),Surface3::Plane(second))=>self.plane_plane(first,second),
            (Surface3::Plane(plane),Surface3::Cylinder(cylinder)) |
            (Surface3::Cylinder(cylinder),Surface3::Plane(plane))=>self.plane_cylinder(plane,cylinder),
            _=>Err(invalid("surface pair is not supported by this analytic intersection slice")),
        }
    }
    fn plane_plane(&self,a:PlaneSurface,b:PlaneSurface)->Result<SurfaceIntersection,GeometryError>{
        valid_frame(a.frame,self.tolerance)?;
        valid_frame(b.frame,self.tolerance)?;
        let d=norm(a.frame.normal.cross(b.frame.normal));
        let epsilon=32.0*f64::EPSILON;
        if d>epsilon && d<=self.tolerance.angular.sin() {
            return Ok(SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearParallel,gap:d});
        }
        let first=Plane3{origin:a.frame.origin,normal:a.frame.normal};
        let second=Plane3{origin:b.frame.origin,normal:b.frame.normal};
        Ok(match first.intersect(second,self.tolerance)? {
            PlaneIntersection::Parallel=>SurfaceIntersection::Empty,
            PlaneIntersection::Coincident=>SurfaceIntersection::Coincident,
            PlaneIntersection::Line{point,direction}=>SurfaceIntersection::Curves{
                contact:IntersectionContact::Crossing,
                curves:vec![IntersectionCurve::Line(Line3{origin:point,direction})],
            },
        })
    }
    fn plane_cylinder(&self,plane:PlaneSurface,cylinder:CylinderSurface)
        ->Result<SurfaceIntersection,GeometryError>
    {
        valid_frame(plane.frame,self.tolerance)?;
        valid_frame(cylinder.frame,self.tolerance)?;
        if !cylinder.radius.is_finite() || cylinder.radius<=self.tolerance.absolute_length {
            return Err(invalid("cylinder radius below tolerance"));
        }
        let n=plane.frame.normal;
        let axis=cylinder.frame.normal;
        let axial=n.dot(axis);
        let relative=plane.frame.origin.sub(cylinder.frame.origin);
        let offset=n.dot(relative);
        let coefficients=[n.dot(cylinder.frame.x),n.dot(cylinder.frame.y)];
        if !axial.is_finite() || !offset.is_finite() || !coefficients.iter().all(|x|x.is_finite()) {
            return Err(invalid("non-finite cylinder/plane orientation or offset"));
        }
        if axial.abs()>self.tolerance.angular.sin(){
            return Ok(SurfaceIntersection::Curves{
                contact:IntersectionContact::Crossing,
                curves:vec![IntersectionCurve::PlaneCylinderSection(PlaneCylinderSection{
                    cylinder,projected_offset:offset,axial,coefficients,
                })],
            });
        }
        if axial.abs()>32.0*f64::EPSILON {
            return Ok(SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearParallel,gap:axial.abs()});
        }
        // Plane is (at the prescribed angular tolerance) parallel to the
        // cylinder axis. Its circular cross-section has 0, 1 or 2 roots.
        let projection=coefficients[0].hypot(coefficients[1]);
        let radius=cylinder.radius*projection;
        if !radius.is_finite() || radius<=self.tolerance.absolute_length {
            return Err(invalid("cylinder cross-section is ill-conditioned"));
        }
        let gap=offset.abs()-radius;
        let rounding_floor=16.0*f64::EPSILON*offset.abs().max(radius).max(1.0);
        let distance_tolerance=self.tolerance.length_at(cylinder.radius);
        if gap.abs()>rounding_floor && gap.abs()<=distance_tolerance {
            return Ok(SurfaceIntersection::Ambiguous{reason:IntersectionAmbiguity::NearTangency,gap});
        }
        if gap>distance_tolerance {return Ok(SurfaceIntersection::Empty);}
        let phase=coefficients[1].atan2(coefficients[0]);
        let tangent=gap.abs()<=rounding_floor;
        let roots=if tangent {vec![if offset>=0.0 {phase}else{phase+std::f64::consts::PI}]}
            else {let theta=(offset/radius).clamp(-1.0,1.0).acos();vec![phase-theta,phase+theta]};
        let mut branches=Vec::with_capacity(roots.len());
        for u in roots {
            let origin=Surface3::Cylinder(cylinder).evaluate(u,0.0)?.point;
            branches.push(IntersectionCurve::Line(Line3{origin,direction:axis}));
        }
        Ok(SurfaceIntersection::Curves{
            contact:if tangent {IntersectionContact::Tangent}else{IntersectionContact::Crossing},
            curves:branches,
        })
    }
}
