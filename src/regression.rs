//! Independent, versioned CAD regression fixtures: closed-form golden values,
//! metamorphic invariants, and deterministic shrinking of failing profiles.
//!
//! A passing kernel self-validation is NOT an oracle. Fixtures record mathematical
//! expectations independent of OCCT or any other CAD engine. Every failure
//! includes the original model and (where possible) a reduced reproducer.
use crate::{
    BrepBody, GeometryError, Point3, ProbeCase, ProbeExpectation, ProbeModel,
};
use crate::accuracy::{build, evaluate_probe_case, mesh_signed_volume};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

pub const REGRESSION_SCHEMA: &str = "rustsolid.regressions.v1";
pub const BUILTIN_REGRESSIONS: &str = include_str!("../tests/fixtures/regressions.v1.json");

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RegressionCorpus {
    pub schema: String,
    pub title: String,
    pub fixtures: Vec<RegressionFixture>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RegressionFixture {
    pub id: String,
    pub description: String,
    pub expectation: ProbeExpectation,
    pub model: ProbeModel,
    #[serde(default, skip_serializing_if="Option::is_none")]
    pub golden: Option<AnalyticGolden>,
    #[serde(default)]
    pub relations: Vec<MetamorphicRelation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalyticGolden {
    pub volume: f64,
    pub surface_area: f64,
    pub centroid: [f64; 3],
    pub bbox_min: [f64; 3],
    pub bbox_max: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub enum MetamorphicRelation {
    Translate { offset: [f64; 3] },
    Scale { factor: f64 },
    ReverseWinding,
    FacetRefinement { coarse: usize, fine: usize },
}
impl MetamorphicRelation {
    fn label(&self) -> &'static str {
        match self {
            Self::Translate{..} => "translation",
            Self::Scale{..} => "uniform_scale",
            Self::ReverseWinding => "reverse_winding",
            Self::FacetRefinement{..} => "facet_refinement",
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RegressionCheck {
    pub name: String,
    pub passed: bool,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct RegressionOutcome {
    pub fixture: RegressionFixture,
    pub passed: bool,
    pub checks: Vec<RegressionCheck>,
    /// Present only if shrinking preserved a base-probe failure signature.
    #[serde(skip_serializing_if="Option::is_none")]
    pub minimized_probe: Option<ProbeCase>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RegressionReport {
    pub schema: String,
    pub title: String,
    pub total: usize,
    pub passed: usize,
    pub failed: usize,
    pub outcomes: Vec<RegressionOutcome>,
}
impl RegressionReport {
    pub fn all_passed(&self) -> bool { self.failed == 0 }
}

fn invalid(reason: impl Into<String>) -> GeometryError {
    GeometryError::InvalidEdit(reason.into())
}
fn ensure_finite(values: impl IntoIterator<Item=f64>, name: &str) -> Result<(),GeometryError> {
    if values.into_iter().all(f64::is_finite) { Ok(()) }
    else { Err(invalid(format!("{name} has non-finite numeric values"))) }
}
fn valid_fixture(f: &RegressionFixture) -> Result<(), GeometryError> {
    if f.id.is_empty() || f.id.len()>128 ||
        !f.id.bytes().all(|b|b.is_ascii_alphanumeric() || matches!(b,b'_'|b'-'|b'.'|b'/')) {
        return Err(invalid("fixture id must be 1..128 ASCII letters/digits/._-/"));
    }
    if f.description.trim().is_empty() || f.description.len()>500 {
        return Err(invalid("fixture description must be 1..500 chars"));
    }
    if f.relations.len()>16 {return Err(invalid("fixture has more than 16 metamorphic relations"));}
    if f.golden.is_some() && f.expectation==ProbeExpectation::Reject {
        return Err(invalid("expected-invalid model must not include analytic golden values"));
    }
    if let Some(g) = &f.golden {
        ensure_finite([g.volume,g.surface_area].into_iter().chain(g.centroid)
            .chain(g.bbox_min).chain(g.bbox_max),"golden values")?;
        if g.volume <= 0.0 || g.surface_area <= 0.0 ||
            (0..3).any(|i|g.bbox_max[i]<=g.bbox_min[i]) {
            return Err(invalid("golden mass or bbox is non-positive"));
        }
    }
    for relation in &f.relations {
        if f.expectation != ProbeExpectation::Accept {
            return Err(invalid("metamorphic transformations need an accepted base model"));
        }
        match relation {
            MetamorphicRelation::Translate{offset} => {
                ensure_finite(*offset,"translation")?;
                if matches!(f.model, ProbeModel::EditBlock {..}) ||
                    (matches!(f.model, ProbeModel::Extrusion{..}) && offset[1]!=0.0) {
                    return Err(invalid("translation not supported for this model representation"));
                }
            }
            MetamorphicRelation::Scale{factor} => {
                if !factor.is_finite() || !(1.0e-4..=1.0e4).contains(factor) ||
                    matches!(f.model, ProbeModel::EditBlock{..}) {
                    return Err(invalid("scale relation requires a positive representable non-edit model"));
                }
            }
            MetamorphicRelation::ReverseWinding => {
                if !matches!(f.model, ProbeModel::Extrusion{..}) {
                    return Err(invalid("winding relation requires an extrusion"));
                }
            }
            MetamorphicRelation::FacetRefinement{coarse,fine} => {
                if !matches!(f.model, ProbeModel::Cylinder{..}) ||
                    !(8..=*fine).contains(coarse) || *fine>16384 || coarse==fine {
                    return Err(invalid("facet refinement requires a cylinder and 8<=coarse<fine<=16384"));
                }
            }
        }
    }
    Ok(())
}
fn check(name: impl Into<String>, passed: bool, detail: impl Into<String>) -> RegressionCheck {
    RegressionCheck {name:name.into(), passed,detail:detail.into()}
}
fn coordinate_tolerance(points: impl IntoIterator<Item=f64>, feature_scale:f64) -> f64 {
    let largest=points.into_iter().fold(1.0_f64,|a,b|a.max(b.abs()));
    (1.0e-8 * feature_scale.abs().max(1.0e-8))
        .max(32.0*f64::EPSILON*largest).max(1.0e-10)
}
fn relative_error(a:f64,b:f64)->f64 {
    (a-b).abs()/b.abs().max(1.0e-20)
}
fn vec_error(a:[f64;3],b:[f64;3])->f64 {
    (a[0]-b[0]).hypot(a[1]-b[1]).hypot(a[2]-b[2])
}
fn coords(p:Point3)->[f64;3]{[p.x,p.y,p.z]}
fn shape_scale(body:&BrepBody)->Result<f64,GeometryError>{
    let m=body.shared()?;
    let a=m.bbox.min;let b=m.bbox.max;
    Ok((b.x-a.x).hypot(b.y-a.y).hypot(b.z-a.z))
}
fn golden_check(model:&ProbeModel, golden:&AnalyticGolden)->Result<RegressionCheck,GeometryError>{
    let (body,_)=build(model)?;
    let graph=body.shared()?;
    let volume=relative_error(graph.mass.volume,golden.volume);
    let area=relative_error(graph.mass.surface_area,golden.surface_area);
    let centroid=vec_error(coords(graph.mass.centroid),golden.centroid);
    let min=vec_error(coords(graph.bbox.min),golden.bbox_min);
    let max=vec_error(coords(graph.bbox.max),golden.bbox_max);
    let extent=shape_scale(&body)?;
    let tolerance=coordinate_tolerance(golden.centroid.into_iter().chain(golden.bbox_min)
        .chain(golden.bbox_max),extent);
    let pass=volume <= 1.0e-9 && area <= 1.0e-9 &&
        centroid <= tolerance && min <= tolerance && max <= tolerance;
    Ok(check("analytic_golden",pass,format!(
        "volume_rel={volume:.3e} area_rel={area:.3e} centroid_abs={centroid:.3e} bbox_min_abs={min:.3e} bbox_max_abs={max:.3e} position_budget={tolerance:.3e}")))
}
fn translated(input:&ProbeModel,delta:[f64;3])->Result<ProbeModel,GeometryError>{
    let add3=|origin:[f64;3]|-> [f64;3] {[origin[0]+delta[0],origin[1]+delta[1],origin[2]+delta[2]]};
    let next=match input {
        ProbeModel::Block{origin,width,height,depth} => ProbeModel::Block{
            origin:add3(*origin),width:*width,height:*height,depth:*depth},
        ProbeModel::Cylinder{origin,radius,height,facets} => ProbeModel::Cylinder{
            origin:add3(*origin),radius:*radius,height:*height,facets:*facets},
        ProbeModel::Extrusion{profile,height} if delta[1]==0.0 => ProbeModel::Extrusion{
            profile:profile.iter().map(|p|[p[0]+delta[0],p[1]+delta[2]]).collect(),height:*height},
        _=>return Err(invalid("translation relation is unsupported for this model")),
    };
    Ok(next)
}
fn scaled(input:&ProbeModel,factor:f64)->Result<ProbeModel,GeometryError>{
    let mul=|o:[f64;3]|->[f64;3]{[o[0]*factor,o[1]*factor,o[2]*factor]};
    Ok(match input {
        ProbeModel::Block{origin,width,height,depth}=>ProbeModel::Block{
            origin:mul(*origin),width:width*factor,height:height*factor,depth:depth*factor},
        ProbeModel::Cylinder{origin,radius,height,facets}=>ProbeModel::Cylinder{
            origin:mul(*origin),radius:radius*factor,height:height*factor,facets:*facets},
        ProbeModel::Extrusion{profile,height}=>ProbeModel::Extrusion{
            profile:profile.iter().map(|p|[p[0]*factor,p[1]*factor]).collect(),height:height*factor},
        _=>return Err(invalid("scaling an edit script is unsupported")),
    })
}
fn relation_check(original:&ProbeModel, relation:&MetamorphicRelation)
    ->Result<RegressionCheck,GeometryError>
{
    let name=relation.label();
    let (body,_) = build(original)?;
    let root=body.shared()?;
    match relation {
        MetamorphicRelation::Translate{offset}=>{
            let changed=build(&translated(original,*offset)?)?.0;
            let m=changed.shared()?;
            let centroid=coords(root.mass.centroid);
            let expected=[centroid[0]+offset[0],centroid[1]+offset[1],centroid[2]+offset[2]];
            let err=vec_error(coords(m.mass.centroid),expected);
            let bbox_min=coords(root.bbox.min);
            let bbox_max=coords(root.bbox.max);
            let want_min=[bbox_min[0]+offset[0],bbox_min[1]+offset[1],bbox_min[2]+offset[2]];
            let want_max=[bbox_max[0]+offset[0],bbox_max[1]+offset[1],bbox_max[2]+offset[2]];
            let bound=vec_error(coords(m.bbox.min),want_min).max(vec_error(coords(m.bbox.max),want_max));
            let scale=shape_scale(&body)?;
            let threshold=coordinate_tolerance(expected.into_iter().chain(want_min).chain(want_max),scale);
            let vol=relative_error(m.mass.volume,root.mass.volume);
            let area=relative_error(m.mass.surface_area,root.mass.surface_area);
            let pass=vol<=1e-9 && area<=1e-9 && err<=threshold && bound<=threshold;
            Ok(check(name,pass,format!("volume_rel={vol:.3e} area_rel={area:.3e} centroid_delta={err:.3e} bbox_delta={bound:.3e} tol={threshold:.3e}")))
        }
        MetamorphicRelation::Scale{factor}=>{
            let changed=build(&scaled(original,*factor)?)?.0;
            let m=changed.shared()?;
            let vol=relative_error(m.mass.volume,root.mass.volume*factor.powi(3));
            let area=relative_error(m.mass.surface_area,root.mass.surface_area*factor.powi(2));
            let scaled3=|p:Point3|->[f64;3]{[p.x*factor,p.y*factor,p.z*factor]};
            let center=scaled3(root.mass.centroid);
            let err=vec_error(coords(m.mass.centroid),center);
            let min=vec_error(coords(m.bbox.min),scaled3(root.bbox.min));
            let max=vec_error(coords(m.bbox.max),scaled3(root.bbox.max));
            let threshold=coordinate_tolerance(center.into_iter().chain(scaled3(root.bbox.min))
                .chain(scaled3(root.bbox.max)),shape_scale(&changed)?);
            let pass=vol<=1e-9&&area<=1e-9&&err<=threshold&&min<=threshold&&max<=threshold;
            Ok(check(name,pass,format!("volume_rel={vol:.3e} area_rel={area:.3e} centroid_delta={err:.3e} bbox_delta={:.3e}",min.max(max))))
        }
        MetamorphicRelation::ReverseWinding=>{
            let ProbeModel::Extrusion{profile,height}=original else {return Err(invalid("winding not applicable"))};
            let mut reverse=profile.clone();reverse.reverse();
            let changed=build(&ProbeModel::Extrusion{profile:reverse,height:*height})?.0;
            let m=changed.shared()?;
            let vol=relative_error(m.mass.volume,root.mass.volume);
            let area=relative_error(m.mass.surface_area,root.mass.surface_area);
            let err=vec_error(coords(m.mass.centroid),coords(root.mass.centroid));
            let threshold=coordinate_tolerance(coords(root.mass.centroid),shape_scale(&body)?);
            let topo_match=m.vertices.len()==root.vertices.len() &&
                m.edges.len()==root.edges.len() && m.faces.len()==root.faces.len();
            Ok(check(name,vol<=1e-9&&area<=1e-9&&err<=threshold&&topo_match,
                format!("volume_rel={vol:.3e} area_rel={area:.3e} centroid_diff={err:.3e} topology_counts_equal={topo_match}")))
        }
        MetamorphicRelation::FacetRefinement{coarse,fine}=>{
            let ProbeModel::Cylinder{radius,height,..}=original else {return Err(invalid("refinement needs cylinder"))};
            let low=body.tessellate(*coarse)?;
            let high=body.tessellate(*fine)?;
            let lo=mesh_signed_volume(&low)?;
            let hi=mesh_signed_volume(&high)?;
            let analytic=PI*radius*radius*height;
            let predicted=|n:usize|{
                let step=TAU/n as f64;
                (n as f64)*radius*radius*step.sin()*0.5*height
            };
            let rel_low=relative_error(lo,predicted(*coarse));
            let rel_high=relative_error(hi,predicted(*fine));
            let err_low=relative_error(lo,analytic);
            let err_high=relative_error(hi,analytic);
            let pass=rel_low<=1e-9&&rel_high<=1e-9&&err_high<err_low
                && root.mass.volume==body.shared()?.mass.volume;
            Ok(check(name,pass,format!("coarse_prediction_rel={rel_low:.3e} fine_prediction_rel={rel_high:.3e} coarse_facet_error={err_low:.3e} fine_facet_error={err_high:.3e}")))
        }
    }
}

fn probe_failure_type(reason:&str)->&'static str {
    if reason.starts_with("unexpected model failure:") { "unexpected_build_failure" }
    else if reason.starts_with("measurement failed:") { "measurement_failure" }
    else if reason.starts_with("invalid boundary model was accepted") { "unexpected_success" }
    else if reason.starts_with("oracle budget exceeded:") { "oracle_budget" }
    else {"other"}
}

/// Deterministically remove redundant profile vertices while keeping the same
/// *failure category*. It never weakens a test or silently promotes a failure
/// to a pass. Detailed shape minimization can be added independently later.
pub fn minimize_failing_probe(probe:&ProbeCase, budget:usize)->Option<ProbeCase>{
    let initial=evaluate_probe_case(probe.clone());
    if initial.passed {return None;}
    let failure=probe_failure_type(&initial.reason);
    let mut best=probe.clone();
    let mut changed=false;
    let mut remaining=budget.min(128);
    while remaining>0 {
        let ProbeModel::Extrusion{profile,height}=&best.input else {break;};
        if profile.len()<=3 {break;}
        let mut reduction=None;
        for index in 0..profile.len() {
            if remaining==0 {break;}
            remaining-=1;
            let mut smaller=profile.clone();smaller.remove(index);
            let mut candidate=best.clone();
            candidate.input=ProbeModel::Extrusion{profile:smaller,height:*height};
            let actual=evaluate_probe_case(candidate.clone());
            if !actual.passed && probe_failure_type(&actual.reason)==failure {
                reduction=Some(candidate);break;
            }
        }
        let Some(candidate)=reduction else {break;};
        best=candidate;changed=true;
    }
    if changed {Some(best)} else {None}
}

pub fn run_regression_corpus(json:&str)->Result<RegressionReport,GeometryError>{
    if json.len()>512_000 {return Err(invalid("regression corpus exceeds 512 KiB"));}
    let corpus:RegressionCorpus=serde_json::from_str(json)
        .map_err(|e|invalid(format!("invalid regression corpus: {e}")))?;
    if corpus.schema!=REGRESSION_SCHEMA {
        return Err(invalid(format!("unknown regression schema: {}",corpus.schema)));
    }
    if corpus.title.trim().is_empty() || corpus.title.len()>200 {
        return Err(invalid("regression title must be 1..200 chars"));
    }
    if corpus.fixtures.is_empty() || corpus.fixtures.len()>256 {
        return Err(invalid("regression corpus requires 1..256 fixtures"));
    }
    let mut unique=BTreeSet::new();
    for fixture in &corpus.fixtures {
        valid_fixture(fixture)?;
        if !unique.insert(fixture.id.as_str()) {
            return Err(invalid(format!("duplicate regression fixture id: {}",fixture.id)));
        }
    }
    let mut outcomes=Vec::with_capacity(corpus.fixtures.len());
    for fixture in corpus.fixtures {
        let probe=ProbeCase {id:fixture.id.clone(),category:"regression" ,
            expectation:fixture.expectation,input:fixture.model.clone()};
        let base=evaluate_probe_case(probe.clone());
        let mut checks=vec![check("independent_probe",base.passed,base.reason.clone())];
        if base.passed && fixture.expectation==ProbeExpectation::Accept {
            if let Some(golden)=&fixture.golden {
                match golden_check(&fixture.model,golden) {
                    Ok(item)=>checks.push(item),
                    Err(e)=>checks.push(check("analytic_golden",false,e.to_string())),
                }
            }
            for relation in &fixture.relations {
                match relation_check(&fixture.model,relation) {
                    Ok(item)=>checks.push(item),
                    Err(e)=>checks.push(check(relation.label(),false,e.to_string())),
                }
            }
        }
        let passed=checks.iter().all(|c|c.passed);
        let minimized_probe=if !base.passed {minimize_failing_probe(&probe,64)}else{None};
        outcomes.push(RegressionOutcome{fixture,passed,checks,minimized_probe});
    }
    let passed=outcomes.iter().filter(|o|o.passed).count();
    Ok(RegressionReport{schema:REGRESSION_SCHEMA.into(),title:corpus.title,
        total:outcomes.len(),passed,failed:outcomes.len()-passed,outcomes})
}

/// Repository-pinned, source-controlled fixtures. In CI the corpus runs even
/// when the seeded exploratory probes happen to choose different parameters.
pub fn run_builtin_regressions()->Result<RegressionReport,GeometryError>{
    run_regression_corpus(BUILTIN_REGRESSIONS)
}
