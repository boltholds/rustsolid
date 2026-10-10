//! Reproducible, boundary-focused geometric accuracy probes.
//!
//! This deliberately avoids Cartesian products of every modeling action.
//! Families target tolerance transitions, winding, nonconvex trims, periodic
//! seams, distant origins, and geometry-preserving editing. Measurements use
//! independent closed-form oracles and metamorphic invariants, not only the
//! kernel's own self-validation. Facet volume is *reported separately* from
//! analytic volume and compared to the correct inscribed-polygon prediction.
use crate::{BrepBody, BrepModel, CylindricalBrep, EdgeId, FaceId,
    GeometryError, GeometryTolerance, Mesh, ParameterRange, Point2, Point3, Solid};
use serde::{Serialize, Deserialize};
use std::f64::consts::{PI, TAU};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProbeExpectation { Accept, Reject }

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag="kind", rename_all="snake_case")]
pub enum ProbeModel {
    Block { origin:[f64;3], width:f64, height:f64, depth:f64 },
    Extrusion { profile:Vec<[f64;2]>, height:f64 },
    Cylinder { origin:[f64;3], radius:f64, height:f64, facets:usize },
    EditBlock { fraction:f64, split_face:bool, adjacent:bool },
}

#[derive(Debug, Clone, Serialize)]
pub struct ProbeCase {
    pub id: String,
    pub category: &'static str,
    pub expectation: ProbeExpectation,
    /// Sufficient information to reproduce failures from a seed or standalone.
    pub input: ProbeModel,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProbeMetrics {
    pub volume_relative_error: f64,
    pub area_relative_error: f64,
    pub centroid_absolute_error: f64,
    pub bbox_max_absolute_error: f64,
    /// Not the analytic kernel's error: the triangle mesh is intentionally approximate.
    pub display_mesh_volume_relative_error: f64,
    /// Difference between measured mesh error and an independent facet formula.
    pub faceting_prediction_absolute_error: f64,
    /// Consistency between 3D curve and 2D trim evaluated on its surface.
    pub max_trim_to_carrier_residual: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProbeOutcome {
    pub case: ProbeCase,
    pub passed: bool,
    pub reason: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub metrics: Option<ProbeMetrics>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProbeReport {
    pub seed: u64,
    pub total: usize,
    pub passed: usize,
    pub failed: usize,
    pub max_volume_relative_error: f64,
    pub max_area_relative_error: f64,
    pub max_mesh_volume_relative_error: f64,
    pub outcomes: Vec<ProbeOutcome>,
}
impl ProbeReport {
    pub fn all_passed(&self) -> bool { self.failed == 0 }
}

struct Splitmix64(u64);
impl Splitmix64 {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z=self.0;
        z=(z^(z>>30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z=(z^(z>>27)).wrapping_mul(0x94D049BB133111EB);
        let z=z^(z>>31);
        (z>>11) as f64 * (1.0/(1u64<<53) as f64)
    }
    fn near(&mut self,midpoint:f64,spread:f64)->f64 {
        midpoint*(1.0+spread*(2.0*self.next()-1.0))
    }
}

/// Hand-chosen boundary families, with bounded seeded perturbations of a few
/// valid near-threshold cases; no uncontrolled combinatorial enumeration.
pub fn generate_boundary_cases(seed:u64) -> Vec<ProbeCase> {
    use ProbeExpectation::{Accept as Yes, Reject as No};
    let mut random=Splitmix64(seed);
    let mut cases=Vec::new();
    let mut add=|id:&str,category:&'static str,expectation:ProbeExpectation,input:ProbeModel| {
        cases.push(ProbeCase{id:id.into(),category,expectation,input});
    };
    let origin=[0.0,0.0,0.0];
    add("block/nominal","analytic_oracle",Yes,ProbeModel::Block{origin,width:2.0,height:3.0,depth:4.0});
    add("block/scale_two","metamorphic_scale",Yes,ProbeModel::Block{origin,width:4.0,height:6.0,depth:8.0});
    add("block/large_offset","large_coordinates",Yes,ProbeModel::Block{
        origin:[1.0e12,1.0e12,-1.0e12],width:1.0,height:1.0,depth:1.0});
    add("block/sliver_valid","near_tolerance",Yes,ProbeModel::Block{
        origin,width:random.near(2.0e-8,0.15),height:1.0,depth:1.0});
    add("block/sliver_rejected","near_tolerance",No,ProbeModel::Block{
        origin,width:1.0e-12,height:1.0,depth:1.0});
    add("block/zero_height","invalid_dimension",No,ProbeModel::Block{
        origin,width:2.0,height:0.0,depth:4.0});
    add("block/negative_depth","invalid_dimension",No,ProbeModel::Block{
        origin,width:2.0,height:3.0,depth:-1.0});

    let sq=vec![[0.0,0.0],[2.0,0.0],[2.0,4.0],[0.0,4.0]];
    let mut reversed=sq.clone(); reversed.reverse(); reversed.push(reversed[0]);
    add("profile/clockwise_closed","orientation_invariance",Yes,ProbeModel::Extrusion{
        profile:reversed,height:3.0});
    let concave=vec![[0.0,0.0],[3.0,0.0],[3.0,1.0],[1.0,1.0],[1.0,3.0],[0.0,3.0]];
    add("profile/concave_L","concave_trims",Yes,ProbeModel::Extrusion{profile:concave,height:2.0});
    add("profile/near_concave_corner","near_tolerance",Yes,ProbeModel::Extrusion{
        profile:vec![[0.0,0.0],[3.0,0.0],[3.0,1.0],
            [random.near(1.0,0.04),1.0],[1.0,3.0],[0.0,3.0]],height:2.0});
    let mut collinear=sq.clone();collinear.insert(1,[1.0,0.0]);
    add("profile/collinear_control","degenerate_controls",Yes,ProbeModel::Extrusion{
        profile:collinear,height:3.0});
    add("profile/far_coordinates","large_coordinates",Yes,ProbeModel::Extrusion{
        profile:vec![[1e12,-1e12],[1e12+1.0,-1e12],
            [1e12+1.0,-1e12+1.0],[1e12,-1e12+1.0]],height:2.0});
    add("profile/self_crossing","invalid_trim",No,ProbeModel::Extrusion{
        profile:vec![[0.0,0.0],[1.0,1.0],[0.0,1.0],[1.0,0.0]],height:2.0});
    add("profile/repeated_vertex","invalid_trim",No,ProbeModel::Extrusion{
        profile:vec![[0.0,0.0],[1.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]],height:2.0});
    add("profile/sub_tolerance_height","invalid_dimension",No,ProbeModel::Extrusion{
        profile:sq,height:1e-12});

    add("cylinder/nominal","analytic_oracle",Yes,ProbeModel::Cylinder{
        origin,radius:2.0,height:5.0,facets:32});
    add("cylinder/eight_facets","mesh_refinement",Yes,ProbeModel::Cylinder{
        origin,radius:2.0,height:5.0,facets:8});
    add("cylinder/fine_facets","mesh_refinement",Yes,ProbeModel::Cylinder{
        origin,radius:2.0,height:5.0,facets:128});
    add("cylinder/far_origin","large_coordinates",Yes,ProbeModel::Cylinder{
        origin:[1e9,-1e9,1e9],radius:2.0,height:5.0,facets:32});
    add("cylinder/small_valid","near_tolerance",Yes,ProbeModel::Cylinder{
        origin,radius:random.near(4e-8,0.1),height:1.0,facets:32});
    add("cylinder/radius_below_tolerance","near_tolerance",No,ProbeModel::Cylinder{
        origin,radius:1e-12,height:1.0,facets:32});
    add("cylinder/height_below_tolerance","near_tolerance",No,ProbeModel::Cylinder{
        origin,radius:2.0,height:1e-12,facets:32});
    add("cylinder/invalid_facets","mesh_boundary",No,ProbeModel::Cylinder{
        origin,radius:2.0,height:5.0,facets:3});

    add("edit/edge_midpoint","topology_invariant",Yes,ProbeModel::EditBlock{
        fraction:0.5,split_face:false,adjacent:false});
    add("edit/edge_near_endpoint","near_tolerance",No,ProbeModel::EditBlock{
        fraction:1e-14,split_face:false,adjacent:false});
    add("edit/face_diagonal","topology_invariant",Yes,ProbeModel::EditBlock{
        fraction:0.5,split_face:true,adjacent:false});
    add("edit/adjacent_vertices","invalid_edit",No,ProbeModel::EditBlock{
        fraction:0.5,split_face:true,adjacent:true});
    cases
}

fn to3(p:[f64;3])->Point3 {Point3{x:p[0],y:p[1],z:p[2]}}
fn box_oracle(origin:[f64;3],w:f64,h:f64,d:f64)->GeometryOracle {
    GeometryOracle{volume:w*h*d,area:2.0*(w*h+w*d+h*d),
        centroid:[origin[0]+w/2.0,origin[1]+h/2.0,origin[2]+d/2.0],
        bbox_min:origin,bbox_max:[origin[0]+w,origin[1]+h,origin[2]+d],faceting:0.0}
}
#[derive(Debug, Clone, Copy)]
struct GeometryOracle {
    volume:f64, area:f64, centroid:[f64;3],bbox_min:[f64;3],bbox_max:[f64;3],faceting:f64,
}
fn prism_oracle(profile:&[[f64;2]],height:f64)->GeometryOracle {
    // Independent centroid/mass oracle calculated relative to the first vertex
    // to avoid catastrophic cancellation for small models far from origin.
    let (mut cross_sum,mut moment_x,mut moment_z,mut perimeter)=(0.0,0.0,0.0,0.0);
    let anchor=profile[0];
    for i in 0..profile.len(){
        let a=profile[i];let b=profile[(i+1)%profile.len()];
        let ux=a[0]-anchor[0];let uz=a[1]-anchor[1];
        let vx=b[0]-anchor[0];let vz=b[1]-anchor[1];
        let cross=ux*vz-vx*uz;
        cross_sum+=cross;
        moment_x+=(ux+vx)*cross;
        moment_z+=(uz+vz)*cross;
        perimeter+=(b[0]-a[0]).hypot(b[1]-a[1]);
    }
    let signed_area=cross_sum*0.5;
    let area=signed_area.abs();
    let min_x=profile.iter().map(|p|p[0]).fold(f64::INFINITY,f64::min);
    let max_x=profile.iter().map(|p|p[0]).fold(f64::NEG_INFINITY,f64::max);
    let min_z=profile.iter().map(|p|p[1]).fold(f64::INFINITY,f64::min);
    let max_z=profile.iter().map(|p|p[1]).fold(f64::NEG_INFINITY,f64::max);
    GeometryOracle { volume:area*height,area:2.0*area+perimeter*height,
        centroid:[anchor[0]+moment_x/(6.0*signed_area),height/2.0,
            anchor[1]+moment_z/(6.0*signed_area)],
        bbox_min:[min_x,0.0,min_z],bbox_max:[max_x,height,max_z],faceting:0.0 }
}
fn cylinder_oracle(origin:[f64;3],r:f64,h:f64,facets:usize)->GeometryOracle {
    let angle=TAU/facets as f64;
    GeometryOracle{volume:PI*r*r*h,area:2.0*PI*r*(r+h),
        centroid:[origin[0],origin[1]+h/2.0,origin[2]],
        bbox_min:[origin[0]-r,origin[1],origin[2]-r],
        bbox_max:[origin[0]+r,origin[1]+h,origin[2]+r],
        // An n-gon inscribed in a circle has area n*r^2*sin(2pi/n)/2.
        faceting:1.0-angle.sin()/angle}
}
fn oracle(model:&ProbeModel)->Option<GeometryOracle>{
    match model {
        ProbeModel::Block{origin,width,height,depth} if *width>0.0 && *height>0.0 && *depth>0.0 =>
            Some(box_oracle(*origin,*width,*height,*depth)),
        ProbeModel::Extrusion{profile,height} if *height>0.0 && profile.len()>=3 =>
            Some(prism_oracle(profile,*height)),
        ProbeModel::Cylinder{origin,radius,height,facets} if *radius>0.0 && *height>0.0 && *facets>=8 =>
            Some(cylinder_oracle(*origin,*radius,*height,*facets)),
        ProbeModel::EditBlock{..}=>Some(box_oracle([0.0,0.0,0.0],2.0,3.0,4.0)),
        _=>None,
    }
}
pub(crate) fn build(model:&ProbeModel)->Result<(BrepBody,usize),GeometryError>{
    match model {
        ProbeModel::Block{origin,width,height,depth}=>
            Solid::block(to3(*origin),*width,*height,*depth).map(|s|(BrepBody::from(s),32)),
        ProbeModel::Extrusion{profile,height}=>{
            let points=profile.iter().map(|p|Point2{x:p[0],z:p[1]}).collect::<Vec<_>>();
            Solid::extrude_xz(&points,*height).map(|s|(BrepBody::from(s),32))
        }
        ProbeModel::Cylinder{origin,radius,height,facets}=>
            CylindricalBrep::upright(to3(*origin),*radius,*height,GeometryTolerance::default())
                .and_then(|s| {s.tessellate(*facets)?;Ok((BrepBody::from(s),*facets))}),
        ProbeModel::EditBlock{fraction,split_face,adjacent}=>{
            let mut solid=Solid::block(to3([0.0,0.0,0.0]),2.0,3.0,4.0)?;
            let face=solid.faces[0].boundary.clone();
            solid.edit_atomic(|tx| {
                tx.split_edge(EdgeId(0),*fraction)?;
                if *split_face {
                    tx.split_face(FaceId(0),face[0],face[if *adjacent {1}else{2}])?;
                }
                Ok(())
            })?;
            Ok((BrepBody::from(solid),32))
        }
    }
}
fn dist(a:[f64;3],b:[f64;3])->f64 {
    (a[0]-b[0]).hypot(a[1]-b[1]).hypot(a[2]-b[2])
}
fn coord(p:Point3)->[f64;3] {[p.x,p.y,p.z]}
fn max_bbox_difference(a:&BrepModel,b:&GeometryOracle)->f64 {
    let d1=dist(coord(a.bbox.min),b.bbox_min);
    let d2=dist(coord(a.bbox.max),b.bbox_max);
    d1.max(d2)
}
fn relative(a:f64,b:f64)->f64 {(a-b).abs()/b.abs().max(1.0e-30)}
pub(crate) fn mesh_signed_volume(mesh:&Mesh)->Result<f64,GeometryError> {
    if mesh.vertices.is_empty() {
        return Err(GeometryError::InvalidTopology("probe mesh has no vertices".into()));
    }
    let anchor=mesh.vertices[0];
    let mut six=0.0;
    for tri in &mesh.triangles {
        if tri.iter().any(|i|*i as usize>=mesh.vertices.len()) {
            return Err(GeometryError::InvalidTopology("probe triangle outside mesh".into()));
        }
        let a=mesh.vertices[tri[0] as usize].sub(anchor);
        let b=mesh.vertices[tri[1] as usize].sub(anchor);
        let c=mesh.vertices[tri[2] as usize].sub(anchor);
        six+=a.dot(b.cross(c));
    }
    Ok(six/6.0)
}
fn trimmed_carrier_residual(graph:&BrepModel)->Result<f64,GeometryError> {
    let mut worst: f64=0.0;
    for coedge in &graph.coedges {
        let edge=graph.edges.get(coedge.edge.0 as usize)
            .ok_or_else(||GeometryError::InvalidTopology("probe coedge missing edge".into()))?;
        let geometry=graph.geometry.coedges.get(coedge.id.0 as usize)
            .ok_or_else(||GeometryError::InvalidTopology("probe coedge missing trim".into()))?;
        let edge_geometry=graph.geometry.edges.get(edge.id.0 as usize)
            .ok_or_else(||GeometryError::InvalidTopology("probe edge missing curve".into()))?;
        let face=&graph.faces[coedge.face.0 as usize];
        let surface=graph.geometry.surfaces[face.geometry.surface.0 as usize];
        let trim=graph.geometry.curves2[geometry.pcurve.0 as usize];
        let curve=graph.geometry.curves3[edge_geometry.curve.0 as usize];
        let (ParameterRange::Bounded{start:a,end:b},ParameterRange::Bounded{start:x,end:y})=
            (geometry.domain,edge_geometry.domain) else {
            return Err(GeometryError::InvalidTopology("probe needs bounded edge parameters".into()));
        };
        for q in [0.0,0.2,0.5,0.8,1.0] {
            let uv=trim.evaluate(a+(b-a)*q)?;
            let from_surface=surface.evaluate(uv.u,uv.v)?.point;
            let f=if coedge.reversed {1.0-q}else{q};
            let from_curve=curve.evaluate(x+(y-x)*f)?;
            worst=worst.max(dist(coord(from_surface),coord(from_curve)));
        }
    }
    Ok(worst)
}

fn measure(body:&BrepBody,facets:usize,expected:GeometryOracle)
    ->Result<ProbeMetrics,GeometryError>
{
    let graph=body.shared()?;
    graph.validate()?;
    let mesh=body.tessellate(facets)?;
    let volume=mesh_signed_volume(&mesh)?;
    let analytic=graph.mass.volume;
    let display_relative=relative(volume,expected.volume);
    let predicted=expected.faceting;
    Ok(ProbeMetrics{
        volume_relative_error:relative(analytic,expected.volume),
        area_relative_error:relative(graph.mass.surface_area,expected.area),
        centroid_absolute_error:dist(coord(graph.mass.centroid),expected.centroid),
        bbox_max_absolute_error:max_bbox_difference(&graph,&expected),
        display_mesh_volume_relative_error:display_relative,
        faceting_prediction_absolute_error:(display_relative-predicted).abs(),
        max_trim_to_carrier_residual:trimmed_carrier_residual(&graph)?,
    })
}
pub fn evaluate_probe_case(case:ProbeCase)->ProbeOutcome {
    let built=build(&case.input);
    match (case.expectation,built) {
        (ProbeExpectation::Reject,Err(error))=>ProbeOutcome{
            case,passed:true,reason:format!("rejected as expected: {error}"),metrics:None},
        (ProbeExpectation::Accept,Err(error))=>ProbeOutcome{
            case,passed:false,reason:format!("unexpected model failure: {error}"),metrics:None},
        (ProbeExpectation::Reject,Ok(_))=>ProbeOutcome{
            case,passed:false,reason:"invalid boundary model was accepted".into(),metrics:None},
        (ProbeExpectation::Accept,Ok((body,facets)))=>{
            let Some(oracle)=oracle(&case.input) else {
                return ProbeOutcome{case,passed:false,reason:"missing independent analytical oracle".into(),metrics:None};
            };
            match measure(&body,facets,oracle) {
                Err(error)=>ProbeOutcome{case,passed:false,reason:format!("measurement failed: {error}"),metrics:None},
                Ok(metrics)=>{
                    // Independent closed forms must match within small relative
                    // budgets. For translated 1e12 parts, centroid/bbox errors
                    // reflect the representable f64 positions, not decimal ideals.
                    let tolerance=1e-8;
                    let centroid_limit=1e-9_f64.max(oracle.volume.abs().powf(1.0/3.0)*1e-8);
                    let bbox_limit=centroid_limit;
                    let coordinate_scale=oracle.bbox_max.iter().chain(&oracle.bbox_min)
                        .fold(1.0_f64,|v,x|v.max(x.abs()));
                    let trim_limit=1e-8_f64.max(32.0*f64::EPSILON*coordinate_scale);
                    let passed=metrics.volume_relative_error<=tolerance
                        && metrics.area_relative_error<=tolerance
                        && metrics.centroid_absolute_error<=centroid_limit
                        && metrics.bbox_max_absolute_error<=bbox_limit
                        && metrics.faceting_prediction_absolute_error<=5e-8
                        && metrics.max_trim_to_carrier_residual<=trim_limit;
                    let reason=if passed {"analytic/orientation/faceting oracles agree".into()}
                        else {format!("oracle budget exceeded: {metrics:?}")};
                    ProbeOutcome{case,passed,reason,metrics:Some(metrics)}
                }
            }
        }
    }
}

/// Execute the bounded boundary corpus with explicit reproducibility via seed.
/// This is intentionally not an exhaustive CAD kernel verification.
pub fn run_boundary_probes(seed:u64) -> ProbeReport {
    let outcomes=generate_boundary_cases(seed).into_iter().map(evaluate_probe_case).collect::<Vec<_>>();
    let passed=outcomes.iter().filter(|o|o.passed).count();
    let mut max_volume_relative_error: f64=0.0;
    let mut max_area_relative_error: f64=0.0;
    let mut max_mesh_volume_relative_error: f64=0.0;
    for metrics in outcomes.iter().filter_map(|o|o.metrics.as_ref()) {
        max_volume_relative_error=max_volume_relative_error.max(metrics.volume_relative_error);
        max_area_relative_error=max_area_relative_error.max(metrics.area_relative_error);
        max_mesh_volume_relative_error=max_mesh_volume_relative_error.max(metrics.display_mesh_volume_relative_error);
    }
    ProbeReport{seed,total:outcomes.len(),passed,failed:outcomes.len()-passed,
        max_volume_relative_error,max_area_relative_error,max_mesh_volume_relative_error,outcomes}
}
