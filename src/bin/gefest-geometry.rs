//! Stateless one-request JSON stdio adapter. The "edit_solid" request runs an
//! entire batch atomically; source solids cannot be externally mutated halfway.
use rustsolid::{geometry_query, run_boundary_probes, run_builtin_regressions,
    ArenaEntityRef, ArenaHandle, ArenaOutput, ArenaStep, ArenaTransaction,
    BooleanKernel, BooleanOperation, BooleanOperand, BooleanOutcome,
    BooleanRegionClassifier, GraphRegionState, RegionHold,
    ThroughHoleBrep, HoleParameters, BlindHoleBrep, BlindHoleParameters,
    BrepCarrierEditor, CarrierHandle, CarrierOwner, CarrierId, BrepCarrierTransaction,
    BrepStructuralEditor, BrepStructuralTransaction, StructuralHoleHandle,
    HoleEdgeRole, HoleTrimRole, StructuralEdgeRef, StructuralTrimRef, StructuralFaceRef,
    BrepArena, TopologyEntity,
    CylinderSurface, Frame3, IntersectionAmbiguity, IntersectionContact, IntersectionCurve,
    ParameterRange, PlaneSurface, SurfaceIntersection, SurfaceIntersectionEngine,
    TrimmedIntersection, TrimmedIntersectionEngine, TrimContact, TrimAmbiguity,
    IntersectionGraphBuilder, GraphState, GraphNodeKind, GraphOperand, TrimPosition,
    TopologyDeltaPlanner, PlanState, PlanArcRole, PlanFaceEffect,
    PlanRequirement, PlanVertexPlacement, QueryLimits, BrepBody, BrepOrigin, CommandBatch, CommandHistory, Curve2, Curve3, CylindricalBrep, EditDelta, EditReport, EdgeId, FaceId, GeometryError, GeometryTolerance, Point2, Point3, Solid, Surface3, VertexId};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::io::{self, Read};
use std::collections::BTreeMap;

#[derive(Debug, Deserialize)]
#[serde(tag = "operation", rename_all = "snake_case")]
enum GeometryRequest {
    ExtrudeProfile { profile: Vec<[f64; 2]>, height: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    Block { origin: [f64; 3], width: f64, height: f64, depth: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    AnalyticCylinder {
        origin: [f64;3], radius: f64, height: f64,
        #[serde(default)] segments: Option<usize>,
        #[serde(default)] tolerance: Option<GeometryTolerance>,
    },
    EditSolid { source: SolidSource, edits: Vec<EditRequest> },
    InspectBrep { source: BrepInput },
    IntersectFaces { first: FaceOperand, second: FaceOperand,
        #[serde(default)] tolerance:Option<GeometryTolerance>,
        #[serde(default)] samples:Option<usize> },
    IntersectionGraph { first: FaceOperand, second: FaceOperand,
        #[serde(default)] tolerance:Option<GeometryTolerance>,
        #[serde(default)] samples:Option<usize> },
    TopologyDeltaPlan { first: FaceOperand, second: FaceOperand,
        #[serde(default)] tolerance:Option<GeometryTolerance>,
        #[serde(default)] samples:Option<usize> },
    ArenaTransaction { source:BrepInput, batches:Vec<Vec<ArenaCommandInput>>,
        #[serde(default)] undo:usize, #[serde(default)] redo:usize },
    Boolean {left:BrepInput, right:BrepInput, mode:BooleanOperation},
    CarrierTransaction { source:BrepInput, batches:Vec<Vec<CarrierCommandInput>>,
        #[serde(default)] undo:usize, #[serde(default)] redo:usize },
    StructuralTransaction {source:SolidSource,batches:Vec<Vec<StructuralCommandInput>>,
        #[serde(default)] undo:usize,#[serde(default)] redo:usize,
        #[serde(default)] display_segments:Option<usize>},
    BooleanRegionPlan { first:FaceOperand, second:FaceOperand, mode:BooleanOperation,
        #[serde(default)] tolerance:Option<GeometryTolerance> },
    GeometryQuery { source: BrepInput, query: String },
    ProbeKernel { #[serde(default)] seed: Option<u64> },
    RegressionCorpus,
    IntersectSurfaces { first: SurfaceInput, second: SurfaceInput,
        #[serde(default)] tolerance:Option<GeometryTolerance>,
        #[serde(default)] samples:Option<usize> },
    CommandHistory { source: SolidSource, #[serde(default)] feature_key: Option<String>,
        batches: Vec<CommandBatch>, #[serde(default)] undo: usize, #[serde(default)] redo: usize },
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum SolidSource {
    ExtrudeProfile { profile: Vec<[f64; 2]>, height: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
    Block { origin: [f64; 3], width: f64, height: f64, depth: f64, #[serde(default)] tolerance: Option<GeometryTolerance> },
}
impl SolidSource {
    fn build(self) -> Result<Solid, GeometryError> {
        match self {
            Self::ExtrudeProfile { profile, height, tolerance } => {
                let points: Vec<Point2> = profile.iter().map(|p| Point2 { x: p[0], z: p[1] }).collect();
                Solid::extrude_xz_with_tolerance(&points, height, tolerance.unwrap_or_default())
            }
            Self::Block { origin, width, height, depth, tolerance } => Solid::block_with_tolerance(
                Point3 { x: origin[0], y: origin[1], z: origin[2] }, width, height, depth,
                tolerance.unwrap_or_default(),
            ),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum EditRequest {
    SplitEdge { edge_id: u32, fraction: f64 },
    SplitFace { face_id: u32, start_vertex_id: u32, end_vertex_id: u32 },
}

#[derive(Serialize)]
struct Topology {
    vertices: usize,
    edges: usize,
    coedges: usize,
    loops: usize,
    shells: usize,
    faces: usize,
    euler_characteristic: isize,
}

#[derive(Serialize)]
struct Mass {
    volume: f64,
    surface_area: f64,
    centroid: [f64; 3],
}

#[derive(Serialize)]
struct SolidResponse {
    vertices: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    triangle_face_ids: Vec<u32>,
    topology: Topology,
    mass: Mass,
    bbox: [[f64; 3]; 2],
    /// Included only for the new batch editing operation; geometry.v1 fields
    /// for legacy block and extrusion remain unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    edit_report: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    history: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    analytic_brep: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    brep: Option<Value>,
}

#[derive(Serialize)]
struct Response {
    schema_version: &'static str,
    ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    solid: Option<SolidResponse>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    data: Option<Value>,
}

fn report_json(report: EditReport) -> Value {
    let changes: Vec<Value> = report.changes.into_iter().map(|event| match event {
        EditDelta::SplitEdge(s) => json!({"kind":"split_edge", "original_edge_id":s.original_edge.0,
            "inserted_vertex_id":s.inserted_vertex.0,"created_edge_id":s.created_edge.0,
            "created_coedge_ids":[s.created_coedges[0].0,s.created_coedges[1].0]}),
        EditDelta::SplitFace(s) => json!({"kind":"split_face", "original_face_id":s.original_face.0,
            "created_face_id":s.created_face.0,"created_loop_id":s.created_loop.0,
            "diagonal_edge_id":s.diagonal_edge.0,
            "created_coedge_ids":[s.created_coedges[0].0,s.created_coedges[1].0]}),
        EditDelta::KillEdgeVertex(k)=>json!({"kind":"kill_edge_vertex",
            "removed_vertex_id":k.removed_vertex.0,"removed_edge_id":k.removed_edge.0,
            "retained_edge_id":k.retained_edge.0}),
        EditDelta::KillEdgeFace(k)=>json!({"kind":"kill_edge_face",
            "removed_edge_id":k.removed_edge.0,"removed_face_id":k.removed_face.0,
            "retained_face_id":k.retained_face.0}),
    }).collect();
    let names = report.named_changes.into_iter().map(|item| json!({
        "entity": format!("{:?}",item.entity),"name":item.name.0,
        "parent_name":item.parent.0
    })).collect::<Vec<_>>();
    json!({"revision_before":report.revision_before,"revision_after":report.revision_after,
        "changes":changes,"name_changes":names,
        "journal":{"frames":report.journal.frames,
            "topology_snapshots":report.journal.topology_snapshots,
            "triangle_snapshots":report.journal.triangle_snapshots}})
}

/// Native periodic-cylinder B-rep returned alongside a display-only mesh.
/// The caller can distinguish true curved geometry through the optional
/// solid.analytic_brep metadata.
fn process_analytic_cylinder(
    origin: [f64; 3], radius: f64, height: f64,
    segments: usize, tolerance: GeometryTolerance,
) -> Response {
    let center = Point3 {x:origin[0],y:origin[1],z:origin[2]};
    let computed = CylindricalBrep::upright(center,radius,height,tolerance)
        .and_then(|body| body.tessellate(segments).map(|mesh|(body,mesh)));
    match computed {
        Err(error) => Response {schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string()),data:None},
        Ok((body,mesh)) => {
            let topo = Topology {
                vertices:body.vertices.len(),edges:body.edges.len(),
                coedges:body.coedges.len(),loops:body.loops.len(),
                shells:body.shells.len(),faces:body.faces.len(),
                euler_characteristic:body.euler_characteristic(),
            };
            let seam = body.seam_edge().0;
            let coedges=body.coedges.iter().map(|c| json!({
                "id":c.id.0,"edge_id":c.edge.0,"face_id":c.face.0,
                "loop_id":c.loop_id.0,"twin_id":c.twin.0,
                "next_id":c.next.0,"prev_id":c.prev.0,
                "reversed":c.reversed,"pcurve_id":body.coedge_geometry[c.id.0 as usize].pcurve.0,
            })).collect::<Vec<_>>();
            let analytic=json!({
                "kind":"right_circular_cylinder",
                "surface_carriers":["plane","plane","cylinder"],
                "curve_carriers":["circle","circle","line"],
                "face_surface_ids":[0,1,2],
                "edge_curve_ids":[0,1,2],
                "side_face_id":2,
                "cap_face_ids":[0,1],
                "circular_edge_ids":[0,1],
                "seam_edge_id":seam,
                "seam_coedge_ids":[4,5],
                "u_periodic":true,
                "u_range":[0.0,std::f64::consts::TAU],
                "v_range":[0.0,height],
                "coedges":coedges,
                "facet_segments":segments,
                "mass_is_analytic":true,
            });
            Response {
                schema_version:"geometry.v1",ok:true,
                solid:Some(SolidResponse {
                    vertices:mesh.vertices.iter().map(|v|v.array()).collect(),
                    faces:mesh.triangles,
                    triangle_face_ids:mesh.triangle_faces.iter().map(|f|f.0).collect(),
                    topology:topo,
                    mass:Mass {volume:body.mass.volume,
                        surface_area:body.mass.surface_area,
                        centroid:body.mass.centroid.array()},
                    bbox:[body.bbox.min.array(),body.bbox.max.array()],
                    edit_report:None,history:None,analytic_brep:Some(analytic),brep:None,
                }),
                error:None,data:None,
            }
        }
    }
}


#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum BrepInput {
    Block { origin:[f64;3], width:f64, height:f64, depth:f64,
        #[serde(default)] tolerance:Option<GeometryTolerance> },
    ExtrudeProfile { profile:Vec<[f64;2]>, height:f64,
        #[serde(default)] tolerance:Option<GeometryTolerance> },
    AnalyticCylinder { origin:[f64;3], radius:f64, height:f64,
        #[serde(default)] segments:Option<usize>,
        #[serde(default)] tolerance:Option<GeometryTolerance> },
    AnalyticThroughHole {origin:[f64;3],width:f64,height:f64,depth:f64,
        center:[f64;2],radius:f64,
        #[serde(default)] segments:Option<usize>,
        #[serde(default)] tolerance:Option<GeometryTolerance>},
    AnalyticBlindHole {origin:[f64;3],width:f64,height:f64,depth:f64,
        center:[f64;2],radius:f64,bore_depth:f64,
        #[serde(default)] segments:Option<usize>,
        #[serde(default)] tolerance:Option<GeometryTolerance>},
}
impl BrepInput {
    fn prepare(self) -> Result<(BrepBody,usize), GeometryError> {
        match self {
            Self::Block {origin,width,height,depth,tolerance} => {
                let solid=SolidSource::Block{origin,width,height,depth,tolerance}.build()?;
                Ok((BrepBody::from(solid),64))
            }
            Self::ExtrudeProfile{profile,height,tolerance} => {
                let solid=SolidSource::ExtrudeProfile{profile,height,tolerance}.build()?;
                Ok((BrepBody::from(solid),64))
            }
            Self::AnalyticCylinder{origin,radius,height,segments,tolerance} => {
                let cylinder=CylindricalBrep::upright(
                    Point3{x:origin[0],y:origin[1],z:origin[2]},
                    radius,height,tolerance.unwrap_or_default())?;
                Ok((BrepBody::from(cylinder),segments.unwrap_or(64)))
            },
            Self::AnalyticThroughHole{origin,width,height,depth,center,radius,segments,tolerance} => {
                let block=Solid::block_with_tolerance(
                    Point3{x:origin[0],y:origin[1],z:origin[2]},
                    width,height,depth,tolerance.unwrap_or_default())?;
                let body=ThroughHoleBrep::from_box(&block,
                    HoleParameters{center_x:center[0],center_z:center[1],radius})?;
                Ok((BrepBody::from(body),segments.unwrap_or(64)))
            },
            Self::AnalyticBlindHole{origin,width,height,depth,center,radius,bore_depth,segments,tolerance}=>{
                let block=Solid::block_with_tolerance(
                    Point3{x:origin[0],y:origin[1],z:origin[2]},
                    width,height,depth,tolerance.unwrap_or_default())?;
                let body=BlindHoleBrep::from_box(&block,BlindHoleParameters{
                    center_x:center[0],center_z:center[1],radius,depth:bore_depth})?;
                Ok((BrepBody::from(body),segments.unwrap_or(64)))
            }
        }
    }
}

fn inspect_brep(input:BrepInput) -> Response {
    let result=input.prepare().and_then(|(source,segments)| {
        let graph=source.shared()?;
        let mesh=source.tessellate(segments)?;
        let stats=graph.summary();
        let surfaces=graph.faces.iter().map(|f| match graph.face_surface(f.id) {
            Some(Surface3::Plane(_)) => "plane",
            Some(Surface3::Cylinder(_)) => "cylinder",
            Some(Surface3::Sphere(_)) => "sphere",
            None => "missing",
        }).collect::<Vec<_>>();
        let curves=graph.edges.iter().map(|e| match graph.edge_curve(e.id) {
            Some(Curve3::Line(_)) => "line",
            Some(Curve3::Circle(_)) => "circle",
            None => "missing",
        }).collect::<Vec<_>>();
        let trims=graph.coedges.iter().map(|c| match graph.coedge_pcurve(c.id) {
            Some(Curve2::Line(_)) => "line",
            Some(Curve2::Circle(_)) => "circle",
            None => "missing",
        }).collect::<Vec<_>>();
        let seams=graph.edges.iter().filter(|e| {
            let [a,b]=e.coedges;
            graph.coedges[a.0 as usize].face==graph.coedges[b.0 as usize].face
        }).map(|e|e.id.0).collect::<Vec<_>>();
        let closed_edges=graph.edges.iter().filter(|e|e.start==e.end)
            .map(|e|e.id.0).collect::<Vec<_>>();
        let coedges=graph.coedges.iter().map(|c| json!({
            "id":c.id.0,"edge":c.edge.0,"face":c.face.0,
            "loop":c.loop_id.0,"next_in_loop":c.next.0,
            "next_of_edge":c.twin.0,"reversed":c.reversed,
        })).collect::<Vec<_>>();
        let report=json!({
            "representation":match graph.origin {
                BrepOrigin::Polyhedral=>"polyhedral",
                BrepOrigin::AnalyticCylinder=>"analytic_cylinder",
                BrepOrigin::AnalyticThroughHole=>"analytic_through_hole",
                BrepOrigin::AnalyticMultipleThroughHoles=>"analytic_multiple_through_holes",
                BrepOrigin::AnalyticBlindHole=>"analytic_blind_hole",
            },
            "topology":{
                "vertices":stats.vertices,"edges":stats.edges,"coedges":stats.coedges,
                "loops":stats.loops,"faces":stats.faces,"shells":stats.shells,
                "euler_characteristic":stats.euler_characteristic,"genus":stats.genus,
            },
            "geometry":{
                "face_surfaces":surfaces,"edge_curves":curves,
                "coedge_trims":trims,"cylindrical_faces":stats.cylindrical_faces,
                "circular_edges":stats.circular_edges,
                "closed_edges":closed_edges,"seam_edges":seams,
            },
            "coedges":coedges,"mesh_is_display_only":true,
        });
        Ok(SolidResponse {
            vertices:mesh.vertices.iter().map(|v|v.array()).collect(),
            faces:mesh.triangles,
            triangle_face_ids:mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology{vertices:stats.vertices,edges:stats.edges,
                coedges:stats.coedges,loops:stats.loops,
                shells:stats.shells,faces:stats.faces,
                euler_characteristic:stats.euler_characteristic},
            mass:Mass{volume:graph.mass.volume,surface_area:graph.mass.surface_area,
                centroid:graph.mass.centroid.array()},
            bbox:[graph.bbox.min.array(),graph.bbox.max.array()],
            edit_report:None,history:None,analytic_brep:None,brep:Some(report),
        })
    });
    match result {
        Ok(solid)=>Response{schema_version:"geometry.v1",ok:true,solid:Some(solid),error:None,data:None},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,solid:None,error:Some(error.to_string()),data:None},
    }
}


/// Typed, read-only B-rep projection. The AST is validated before resolving
/// any fields and execution is protected by parser and output budgets.
fn evaluate_geometry_query(input: BrepInput, query: &str) -> Response {
    let data = input.prepare().map_err(|e|e.to_string()).and_then(|(body,_)| {
        let graph=body.shared().map_err(|e|e.to_string())?;
        geometry_query(&graph,query,QueryLimits::default()).map_err(|e|e.to_string())
    });
    match data {
        Ok(data)=>Response {schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(error)=>Response {schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(error),data:None},
    }
}

/// Analytic, metamorphic and near-degenerate fixtures. A failure preserves
/// the full set of reproducible case inputs in response.data.
fn evaluate_boundary_probes(seed:u64) -> Response {
    let report=run_boundary_probes(seed);
    let success=report.all_passed();
    match serde_json::to_value(&report) {
        Ok(data)=>Response {schema_version:"geometry.v1",ok:success,
            solid:None,error:if success {None}
                else {Some(format!("{} of {} boundary probes failed",report.failed,report.total))},
            data:Some(data)},
        Err(e)=>Response {schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(format!("could not serialize probe report: {e}")),data:None},
    }
}


/// Source-controlled, mathematically specified regression fixtures.
/// Unlike exploratory seed probes, these are permanently pinned in Git.
fn evaluate_regression_corpus() -> Response {
    match run_builtin_regressions() {
        Err(error)=>Response {schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string()),data:None},
        Ok(report)=>{
            let success=report.all_passed();
            match serde_json::to_value(&report) {
                Ok(data)=>Response{schema_version:"geometry.v1",ok:success,solid:None,
                    error:if success {None}
                        else {Some(format!("{} of {} pinned regressions failed",report.failed,report.total))},
                    data:Some(data)},
                Err(error)=>Response{schema_version:"geometry.v1",ok:false,solid:None,
                    error:Some(format!("cannot serialize regression report: {error}")),data:None},
            }
        }
    }
}


#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum SurfaceInput {
    Plane { origin:[f64;3], x_axis:[f64;3], y_axis:[f64;3] },
    Cylinder { origin:[f64;3], x_axis:[f64;3], y_axis:[f64;3], radius:f64 },
}
impl SurfaceInput {
    fn analytic(self,tol:GeometryTolerance)->Result<Surface3,GeometryError>{
        let vec=|p:[f64;3]|Point3{x:p[0],y:p[1],z:p[2]};
        match self {
            Self::Plane{origin,x_axis,y_axis}=>{
                let frame=Frame3::from_axes(vec(origin),vec(x_axis),vec(y_axis),tol)?;
                Ok(Surface3::Plane(PlaneSurface{frame}))
            },
            Self::Cylinder{origin,x_axis,y_axis,radius}=>{
                let frame=Frame3::from_axes(vec(origin),vec(x_axis),vec(y_axis),tol)?;
                Ok(Surface3::Cylinder(CylinderSurface::new(frame,radius,tol)?))
            },
        }
    }
}

fn evaluate_surface_intersection(
    first:SurfaceInput, second:SurfaceInput,
    tolerance:GeometryTolerance, samples:usize,
)->Response{
    let result=(||->Result<Value,GeometryError>{
        if !(3..=64).contains(&samples){
            return Err(GeometryError::InvalidEdit("intersection sample count must be 3..=64".into()));
        }
        let first=first.analytic(tolerance)?;
        let second=second.analytic(tolerance)?;
        let engine=SurfaceIntersectionEngine::new(tolerance)?;
        let intersection=engine.intersect(first,second)?;
        let data=match &intersection {
            SurfaceIntersection::Empty=>json!({"kind":"empty","curves":[]}),
            SurfaceIntersection::Coincident=>json!({"kind":"coincident","curves":[]}),
            SurfaceIntersection::Ambiguous{reason,gap}=>json!({
                "kind":"ambiguous",
                "reason":match reason {
                    IntersectionAmbiguity::NearParallel=>"near_parallel",
                    IntersectionAmbiguity::NearTangency=>"near_tangency"
                },
                "signed_gap":gap,"curves":[],"requires_refinement":true
            }),
            SurfaceIntersection::Curves{contact,curves}=>{
                let certificate=intersection.certify(first,second,tolerance)?;
                let branches=curves.iter().map(|curve| {
                    let bounded=matches!(curve.parameter_domain(),ParameterRange::Bounded{..});
                    let coords=(0..samples).map(|i|{
                        let f=i as f64/(samples-1) as f64;
                        let t=if bounded {f*std::f64::consts::TAU}else{8.0*f-4.0};
                        curve.evaluate(t).map(|p|json!({"parameter":t,"point":p.array()}))
                    }).collect::<Result<Vec<_>,_>>()?;
                    Ok(json!({
                        "kind":match curve {
                            IntersectionCurve::Line(_)=>"line",
                            IntersectionCurve::PlaneCylinderSection(_)=>"periodic_section",
                        },
                        "periodic":bounded,"samples":coords
                    }))
                }).collect::<Result<Vec<Value>,GeometryError>>()?;
                json!({
                    "kind":intersection.kind(),
                    "contact":match contact{
                        IntersectionContact::Crossing=>"crossing",
                        IntersectionContact::Tangent=>"tangent",
                    },
                    "curves":branches,
                    "certificate":{
                        "sample_count":certificate.samples,
                        "max_residual_first":certificate.max_residual_a,
                        "max_residual_second":certificate.max_residual_b,
                        "budget":certificate.budget
                    }
                })
            }
        };
        Ok(data)
    })();
    match result {
        Ok(data)=>Response {schema_version:"geometry.v1",ok:true,solid:None,error:None,data:Some(data)},
        Err(error)=>Response {schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string()),data:None},
    }
}


#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct FaceOperand {source:BrepInput,face_id:u32}

/// Bounded, read-only analytic curve clipping against two selected B-rep faces.
/// Actual face p-curves determine finite spans; display tessellation is ignored.
fn evaluate_face_intersection(
    first:FaceOperand,second:FaceOperand,
    tolerance:GeometryTolerance,samples:usize,
)->Response{
    let result=(||->Result<Value,GeometryError>{
        if !(3..=64).contains(&samples) {
            return Err(GeometryError::InvalidEdit("trimmed intersection samples must be 3..=64".into()));
        }
        let (a,_) = first.source.prepare()?;
        let (b,_) = second.source.prepare()?;
        let model_a=a.shared()?;
        let model_b=b.shared()?;
        let engine=TrimmedIntersectionEngine::new(tolerance)?;
        let trimmed=engine.intersect_faces(&model_a,FaceId(first.face_id),&model_b,FaceId(second.face_id))?;
        let data=match trimmed {
            TrimmedIntersection::Empty=>json!({"kind":"empty","branches":[],"contacts":[]}),
            TrimmedIntersection::Ambiguous{reason}=>json!({
                "kind":"ambiguous",
                "reason":match reason {
                    TrimAmbiguity::NearParallel=>"near_parallel",
                    TrimAmbiguity::NearTangency=>"near_tangency",
                    TrimAmbiguity::CoincidentSupports=>"coincident_supports",
                    TrimAmbiguity::UnresolvedBoundaryEvent=>"unresolved_boundary_event"
                },
                "branches":[],"contacts":[],"requires_refinement":true
            }),
            TrimmedIntersection::Geometry{branches,contacts,certificate}=>{
                let result_kind=if branches.is_empty() {"points"}
                    else if contacts.is_empty() {"segments"}else{"mixed"};
                let output=branches.iter().map(|branch|{
                    let coords=(0..samples).map(|i|{
                        let fraction=i as f64/(samples-1) as f64;
                        Ok(json!({"fraction":fraction,"point":branch.evaluate(fraction)?.array()}))
                    }).collect::<Result<Vec<Value>,GeometryError>>()?;
                    Ok(json!({
                        "contact":match branch.contact {
                            TrimContact::Crossing=>"crossing",
                            TrimContact::Tangent=>"tangent",
                            TrimContact::BoundaryOverlap=>"boundary_overlap"
                        },
                        "curve":match branch.carrier {
                            IntersectionCurve::Line(_)=>"line",
                            IntersectionCurve::PlaneCylinderSection(_)=>"periodic_section"
                        },
                        "parameter_interval":[branch.parameter_start,branch.parameter_end],
                        "start":branch.start.array(),"end":branch.end.array(),
                        "closed":branch.closed,"samples":coords
                    }))
                }).collect::<Result<Vec<Value>,GeometryError>>()?;
                let output_contacts=contacts.iter().map(|point|json!({
                    "point":point.point.array(),
                    "contact":match point.contact {
                        TrimContact::Crossing=>"crossing",
                        TrimContact::Tangent=>"tangent",
                        TrimContact::BoundaryOverlap=>"boundary_overlap"
                    }
                })).collect::<Vec<_>>();
                json!({
                    "kind":result_kind,
                    "branches":output,"contacts":output_contacts,
                    "sampled_support_certificate":{
                        "sample_count":certificate.samples,
                        "max_residual_first":certificate.max_residual_a,
                        "max_residual_second":certificate.max_residual_b,
                        "budget":certificate.budget
                    }
                })
            }
        };
        Ok(data)
    })();
    match result {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string()),data:None},
    }
}


/// A read-only candidate-intersection graph with per-face provenance and
/// tolerance evidence. Does NOT execute a Boolean or topology edit.
fn evaluate_intersection_graph(
    first:FaceOperand, second:FaceOperand, tolerance:GeometryTolerance,samples:usize,
)->Response{
    let result=(||->Result<Value,GeometryError>{
        if !(3..=64).contains(&samples) {
            return Err(GeometryError::InvalidEdit("intersection graph sample count must be 3..=64".into()));
        }
        let (left,_) = first.source.prepare()?;
        let (right,_) = second.source.prepare()?;
        let a=left.shared()?;
        let b=right.shared()?;
        let graph=IntersectionGraphBuilder::new(tolerance)?
            .intersect_faces(&a,FaceId(first.face_id),&b,FaceId(second.face_id))?;
        graph.validate_against(&a,FaceId(first.face_id),&b,FaceId(second.face_id))?;
        let operand=|o:GraphOperand|match o{GraphOperand::First=>"first",GraphOperand::Second=>"second"};
        let membership=|p:TrimPosition|match p{
            TrimPosition::Outside=>"outside",TrimPosition::Interior=>"interior",
            TrimPosition::Boundary=>"boundary",
        };
        let contact=|c:TrimContact|match c{
            TrimContact::Crossing=>"crossing",TrimContact::Tangent=>"tangent",
            TrimContact::BoundaryOverlap=>"boundary_overlap"
        };
        let sources=graph.source_faces.iter().map(|src|json!({
            "operand":operand(src.operand),"face_id":src.face.0,
            "surface_id":src.surface.0,
            "tolerance":{"absolute_length":src.tolerance.absolute_length,
                "relative_length":src.tolerance.relative_length,
                "angular":src.tolerance.angular}
        })).collect::<Vec<_>>();
        let nodes=graph.nodes.iter().map(|node|{
            let edges=node.boundary_sources.iter().map(|src|json!({
                "operand":operand(src.operand),"face_id":src.face.0,
                "edge_id":src.edge.0,"coedge_id":src.coedge.0,
            })).collect::<Vec<_>>();
            json!({"id":node.id.0,"position":node.point.array(),
                "kind":match node.kind{
                    GraphNodeKind::Endpoint=>"endpoint",
                    GraphNodeKind::IsolatedContact=>"isolated_contact",
                    GraphNodeKind::PeriodicAnchor=>"periodic_anchor",
                },
                "contact":node.contact.map(contact),
                "membership":[membership(node.membership[0]),membership(node.membership[1])],
                "boundary_sources":edges
            })
        }).collect::<Vec<_>>();
        let arcs=graph.arcs.iter().map(|edge|{
            let points=(0..samples).map(|i|{
                let f=i as f64/(samples-1) as f64;
                Ok(json!({"fraction":f,"point":edge.evaluate(f)?.array()}))
            }).collect::<Result<Vec<Value>,GeometryError>>()?;
            Ok(json!({"id":edge.id.0,
                "from":edge.from.0,"to":edge.to.0,
                "closed":edge.closed,
                "contact":contact(edge.contact),
                "carrier":match edge.curve {
                    IntersectionCurve::Line(_)=>"line",
                    IntersectionCurve::PlaneCylinderSection(_)=>"periodic_section"
                },
                "parameter_interval":[edge.parameter_start,edge.parameter_end],
                "samples":points
            }))
        }).collect::<Result<Vec<Value>,GeometryError>>()?;
        let status=match graph.state {
            GraphState::Empty=>"empty",GraphState::Ready=>"ready",
            GraphState::Ambiguous(_)=>"ambiguous"
        };
        let reason=match graph.state {
            GraphState::Ambiguous(TrimAmbiguity::NearParallel)=>Some("near_parallel"),
            GraphState::Ambiguous(TrimAmbiguity::NearTangency)=>Some("near_tangency"),
            GraphState::Ambiguous(TrimAmbiguity::CoincidentSupports)=>Some("coincident_supports"),
            GraphState::Ambiguous(TrimAmbiguity::UnresolvedBoundaryEvent)=>Some("unresolved_boundary_event"),
            _=>None
        };
        let certificate=graph.support_certificate.map(|c|json!({
            "sample_count":c.samples,"max_residual_first":c.max_residual_a,
            "max_residual_second":c.max_residual_b,"budget":c.budget
        }));
        Ok(json!({
            "kind":status,"reason":reason,"sources":sources,"nodes":nodes,"arcs":arcs,
            "local_tolerance":{"length":graph.tolerance.length,
                "angular":graph.tolerance.angular,
                "evaluation_budget":graph.tolerance.evaluation_budget},
            "support_certificate":certificate,
            "topology_committed":false,
            "requires_refinement":status=="ambiguous"
        }))
    })();
    match result {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(error.to_string()),data:None},
    }
}


/// Symbolic, source-bound modeling intent: NO mutation or Boolean commit.
fn evaluate_topology_delta_plan(first:FaceOperand,second:FaceOperand,
    tolerance:GeometryTolerance,samples:usize)->Response
{
    let result=(||->Result<Value,GeometryError>{
        if !(3..=64).contains(&samples) {
            return Err(GeometryError::InvalidEdit("plan sampling budget must be 3..=64".into()));
        }
        let (first_body,_) = first.source.prepare()?;
        let (second_body,_) = second.source.prepare()?;
        let a=first_body.shared()?;
        let b=second_body.shared()?;
        let fa=FaceId(first.face_id);
        let fb=FaceId(second.face_id);
        let graph=IntersectionGraphBuilder::new(tolerance)?.intersect_faces(&a,fa,&b,fb)?;
        let plan=TopologyDeltaPlanner::new().plan(&graph,&a,fa,&b,fb)?;
        plan.validate_against(&graph,&a,fa,&b,fb)?;
        let source_label=|x:GraphOperand| match x {
            GraphOperand::First=>"first",GraphOperand::Second=>"second",
        };
        let placement=|v:PlanVertexPlacement|match v {
            PlanVertexPlacement::VirtualPeriodicAnchor=>json!({"kind":"virtual_periodic_anchor"}),
            PlanVertexPlacement::ContactOnly=>json!({"kind":"contact_only"}),
            PlanVertexPlacement::ReuseVertex{vertex}=>json!({"kind":"reuse_vertex","vertex_id":vertex.0}),
            PlanVertexPlacement::SplitEdge{edge,coedge,parameter,fraction}=>json!({
                "kind":"split_edge","edge_id":edge.0,"coedge_id":coedge.0,
                "parameter":parameter,"fraction":fraction
            }),
            PlanVertexPlacement::CreateFaceVertex{face}=>json!({
                "kind":"create_face_vertex","face_id":face.0
            }),
        };
        let effect=|e:PlanFaceEffect|match e {
            PlanFaceEffect::CandidateBoundarySplit=>"candidate_boundary_split",
            PlanFaceEffect::PendingOpenTrace=>"pending_open_trace",
            PlanFaceEffect::CandidateInnerLoop=>"candidate_inner_loop",
            PlanFaceEffect::ReconcileExistingBoundary=>"reconcile_existing_boundary",
            PlanFaceEffect::RecordTangency=>"record_tangency",
        };
        let requirement=|r:PlanRequirement|match r {
            PlanRequirement::NumericalRefinement=>"numerical_refinement",
            PlanRequirement::GeneralizedBrepEditor=>"generalized_brep_editor",
            PlanRequirement::RegionClassification=>"region_classification",
            PlanRequirement::InnerLoopConstruction=>"inner_loop_construction",
            PlanRequirement::OpenTraceCompletion=>"open_trace_completion",
            PlanRequirement::CurvedEdgeSubdivision=>"curved_edge_subdivision",
            PlanRequirement::MultipleSplitsSameEdge=>"multiple_splits_same_edge",
            PlanRequirement::TangencyClassification=>"tangency_classification",
            PlanRequirement::BoundaryReconciliation=>"boundary_reconciliation",
        };
        let source_faces=plan.source_faces.iter().map(|source|json!({
            "operand":source_label(source.operand),"face_id":source.face.0,
            "surface_id":source.surface.0,
            "tolerance":{"absolute_length":source.tolerance.absolute_length,
                "relative_length":source.tolerance.relative_length,
                "angular":source.tolerance.angular}
        })).collect::<Vec<_>>();
        let nodes=plan.vertices.iter().map(|node|json!({
            "graph_node_id":node.node.0,
            "operand":source_label(node.operand),"face_id":node.face.0,
            "point":node.point.array(),"placement":placement(node.placement),
            "physical_boundary_sources":node.physical_boundary_sources.iter().map(|w|json!({
                "edge_id":w.edge.0,"coedge_id":w.coedge.0,
                "operand":source_label(w.operand),"face_id":w.face.0
            })).collect::<Vec<_>>(),
        })).collect::<Vec<_>>();
        let arcs=plan.arcs.iter().map(|arc|{
            let samples_json=(0..samples).map(|i|{
                let fraction=i as f64/(samples-1) as f64;
                arc.analytic_curve.evaluate(arc.parameter_start+
                    fraction*(arc.parameter_end-arc.parameter_start))
                .map(|point|json!({"fraction":fraction,"point":point.array()}))
            }).collect::<Result<Vec<Value>,GeometryError>>()?;
            Ok(json!({
                "graph_arc_id":arc.arc.0,"from":arc.from.0,"to":arc.to.0,
                "closed":arc.closed,
                "contact":match arc.contact {
                    TrimContact::Crossing=>"crossing",TrimContact::Tangent=>"tangent",
                    TrimContact::BoundaryOverlap=>"boundary_overlap",
                },
                "carrier_intent":match arc.role {
                    PlanArcRole::ProposeNewCarrier=>"propose_new_carrier",
                    PlanArcRole::ReconcileBoundary=>"reconcile_boundary",
                    PlanArcRole::ContactOnly=>"contact_only",
                },
                "carrier_type":match arc.analytic_curve {
                    IntersectionCurve::Line(_)=>"line",
                    IntersectionCurve::PlaneCylinderSection(_)=>"periodic_section",
                },
                "interval":[arc.parameter_start,arc.parameter_end],
                "samples":samples_json
            }))
        }).collect::<Result<Vec<Value>,GeometryError>>()?;
        let faces=plan.face_changes.iter().map(|face|json!({
            "graph_arc_id":face.arc.0,
            "operand":source_label(face.operand),
            "face_id":face.face.0,
            "effect":effect(face.effect)
        })).collect::<Vec<_>>();
        let requirements=plan.requirements.iter().copied().map(requirement).collect::<Vec<_>>();
        Ok(json!({
            "kind":"topology_delta_plan",
            "state":match plan.state{
                PlanState::NoChange=>"no_change",
                PlanState::PreviewOnly=>"preview_only",
                PlanState::Held=>"held",
            },
            "graph_state":graph.kind(),"sources":source_faces,
            "local_tolerance":{
                "length":plan.graph_tolerance.length,
                "angular":plan.graph_tolerance.angular,
                "evaluation_budget":plan.graph_tolerance.evaluation_budget
            },
            "node_bindings":nodes,"intersection_arcs":arcs,
            "face_changes":faces,"requirements":requirements,
            "split_edge_requests":plan.split_edge_requests(),
            "can_apply":plan.can_apply,"topology_committed":false
        }))
    })();
    match result {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,solid:None,error:Some(error.to_string()),data:None},
    }
}


#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ArenaOutputInput {
    /// Zero-based index of a previous command in the current JSON batch.
    step: usize,
    output: ArenaOutput,
}

#[derive(Debug, Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
enum ArenaCommandInput {
    SplitEdge {
        #[serde(default)] edge_id: Option<u32>,
        #[serde(default)] edge_ref: Option<ArenaOutputInput>,
        fraction: f64,
    },
    SplitFace {
        #[serde(default)] face_id: Option<u32>,
        #[serde(default)] face_ref: Option<ArenaOutputInput>,
        #[serde(default)] start_vertex_id: Option<u32>,
        #[serde(default)] start_vertex_ref: Option<ArenaOutputInput>,
        #[serde(default)] end_vertex_id: Option<u32>,
        #[serde(default)] end_vertex_ref: Option<ArenaOutputInput>,
    },
    KillEdgeVertex {
        #[serde(default)] vertex_id: Option<u32>,
        #[serde(default)] vertex_ref: Option<ArenaOutputInput>,
    },
    KillEdgeFace {
        #[serde(default)] edge_id: Option<u32>,
        #[serde(default)] edge_ref: Option<ArenaOutputInput>,
        #[serde(default)] removed_face_id: Option<u32>,
        #[serde(default)] removed_face_ref: Option<ArenaOutputInput>,
    },
}

// All integer ID selectors refer to entities from the pre-batch snapshot.
// A produced reference may point only to a previously staged step, not an
// arbitrary dense ID that happened to become available during this batch.
fn arena_cli_ref(
    handles: &BTreeMap<TopologyEntity, ArenaHandle>,
    steps: &[ArenaStep],
    original: Option<u32>,
    produced: Option<ArenaOutputInput>,
    kind: &str,
) -> Result<ArenaEntityRef, GeometryError> {
    match (original, produced) {
        (Some(id), None) => {
            let entity = match kind {
                "edge" => TopologyEntity::Edge(EdgeId(id)),
                "face" => TopologyEntity::Face(FaceId(id)),
                "vertex" => TopologyEntity::Vertex(VertexId(id)),
                _ => return Err(GeometryError::InvalidEdit("invalid arena selector kind".into())),
            };
            handles.get(&entity).copied().map(ArenaEntityRef::Existing)
                .ok_or_else(||GeometryError::InvalidEdit(format!(
                    "arena transaction references unknown {entity:?}")))
        }
        (None, Some(output)) => {
            let step = steps.get(output.step).copied()
                .ok_or_else(||GeometryError::InvalidEdit(
                    "arena output step must precede its use in the same batch".into()))?;
            Ok(step.output(output.output))
        }
        _ => Err(GeometryError::InvalidEdit(format!(
            "arena {kind} operand needs exactly one original ID or output ref"))),
    }
}

impl ArenaCommandInput {
    fn stage(
        self,
        tx: &mut ArenaTransaction<'_>,
        handles: &BTreeMap<TopologyEntity, ArenaHandle>,
        steps: &[ArenaStep],
    ) -> Result<ArenaStep, GeometryError> {
        let r = |id, output, kind| arena_cli_ref(handles, steps, id, output, kind);
        match self {
            Self::SplitEdge { edge_id, edge_ref, fraction } => {
                tx.split_edge_ref(r(edge_id, edge_ref, "edge")?, fraction)
            }
            Self::SplitFace { face_id, face_ref,
                start_vertex_id, start_vertex_ref, end_vertex_id, end_vertex_ref } => {
                let face = r(face_id, face_ref, "face")?;
                let start = r(start_vertex_id, start_vertex_ref, "vertex")?;
                let end = r(end_vertex_id, end_vertex_ref, "vertex")?;
                tx.split_face_ref(face, start, end)
            }
            Self::KillEdgeVertex { vertex_id, vertex_ref } => {
                tx.kill_edge_vertex_ref(r(vertex_id, vertex_ref, "vertex")?)
            }
            Self::KillEdgeFace { edge_id, edge_ref, removed_face_id, removed_face_ref } => {
                let edge = r(edge_id, edge_ref, "edge")?;
                let face = r(removed_face_id, removed_face_ref, "face")?;
                tx.kill_edge_face_ref(edge, face)
            }
        }
    }
}
fn arena_binding_json(binding:&rustsolid::ArenaBinding)->Value{
    json!({"name":binding.name.as_str(),
        "entity":format!("{:?}",binding.entity),
        "handle":{
            "arena_token":binding.handle.arena_token,
            "slot":binding.handle.slot,
            "generation":binding.handle.generation
        }
    })
}
fn evaluate_arena_transaction(
    source:BrepInput,batches:Vec<Vec<ArenaCommandInput>>,undo:usize,redo:usize,
)->Response{
    let output=(||->Result<Value,GeometryError>{
        if batches.len()>64 || undo>64 || redo>64 {
            return Err(GeometryError::InvalidEdit("arena session exceeds 64 batch/transition budget".into()));
        }
        let (source,_)=source.prepare()?;
        let mut arena=BrepArena::new(source)?;
        let mut reports=Vec::with_capacity(batches.len());
        for batch in batches {
            if batch.is_empty() || batch.len()>256 {
                return Err(GeometryError::InvalidEdit("arena batch needs 1..=256 commands".into()));
            }
            // Numeric IDs are bound to pre-batch handles; output references
            // resolve strictly against previous stages of this batch.
            let handles = arena.snapshot()?.handles;
            let mut tx = arena.begin_transaction()?;
            let mut steps = Vec::with_capacity(batch.len());
            for command in batch {
                let step = command.stage(&mut tx, &handles, &steps)?;
                steps.push(step);
            }
            let commit=tx.commit()?;
            let created=commit.created.iter().map(arena_binding_json).collect::<Vec<_>>();
            let retired=commit.retired.iter().map(arena_binding_json).collect::<Vec<_>>();
            let remapped=commit.remapped.iter().map(|m|json!({
                "name":m.name.as_str(),"before":format!("{:?}",m.before),
                "after":format!("{:?}",m.after),
                "handle":{"arena_token":m.handle.arena_token,
                    "slot":m.handle.slot,"generation":m.handle.generation}
            })).collect::<Vec<_>>();
            reports.push(json!({
                "revision_before":commit.report.revision_before,
                "revision_after":commit.report.revision_after,
                "commands":commit.report.changes.len(),
                "created":created,"retired":retired,"remapped":remapped,
                "journal_frames":commit.report.journal.frames
            }));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            let value=arena.undo()?.ok_or_else(||GeometryError::InvalidEdit(
                "arena undo exceeds committed history".into()))?;
            transitions.push(json!({"kind":"undo","revision_after":value.revision_after}));
        }
        for _ in 0..redo {
            let value=arena.redo()?.ok_or_else(||GeometryError::InvalidEdit(
                "arena redo exceeds available history".into()))?;
            transitions.push(json!({"kind":"redo","revision_after":value.revision_after}));
        }
        arena.validate()?;
        let snapshot=arena.snapshot()?;
        let stats=snapshot.model.summary();
        let mass=snapshot.model.mass;
        Ok(json!({
            "kind":"brep_arena",
            "representation":if arena.is_editable(){"polyhedral"}else{"analytic_cylinder"},
            "revision":snapshot.revision,
            "arena_token":snapshot.arena_token,
            "is_editable":arena.is_editable(),
            "live_slots":arena.live_slots(),
            "allocated_slots":arena.allocated_slots(),
            "undo_depth":arena.undo_depth(),
            "redo_depth":arena.redo_depth(),
            "mass":{"volume":mass.volume,"surface_area":mass.surface_area,
                "centroid":mass.centroid.array()},
            "topology":{"vertices":stats.vertices,"edges":stats.edges,
                "coedges":stats.coedges,"loops":stats.loops,
                "faces":stats.faces,"shells":stats.shells,
                "euler_characteristic":stats.euler_characteristic},
            "commits":reports,"transitions":transitions,
            "note":"Handle tokens are process-local and not transferable across CLI requests"
        }))
    })();
    match output {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,solid:None,error:Some(error.to_string()),data:None},
    }
}



fn evaluate_boolean_region_plan(
    first:FaceOperand,second:FaceOperand,mode:BooleanOperation,
    tolerance:GeometryTolerance
)->Response {
    let result=(||->Result<Value,GeometryError>{
        let (left,_)=first.source.prepare()?;
        let (right,_)=second.source.prepare()?;
        let fa=FaceId(first.face_id);
        let fb=FaceId(second.face_id);
        let lm=left.shared()?;
        let rm=right.shared()?;
        let graph=IntersectionGraphBuilder::new(tolerance)?.intersect_faces(&lm,fa,&rm,fb)?;
        let plan=BooleanRegionClassifier::classify_graph(
            &graph,&left,fa,&right,fb,mode)?;
        plan.validate_against(&graph,&left,fa,&right,fb)?;
        let topology=TopologyDeltaPlanner::new().plan(&graph,&lm,fa,&rm,fb)?;
        topology.validate_against(&graph,&lm,fa,&rm,fb)?;
        let outstanding=topology.requirements.iter().map(|req|match req{
            PlanRequirement::NumericalRefinement=>"numerical_refinement",
            PlanRequirement::GeneralizedBrepEditor=>"generalized_brep_editor",
            PlanRequirement::RegionClassification=>"region_classification",
            PlanRequirement::InnerLoopConstruction=>"inner_loop_construction",
            PlanRequirement::OpenTraceCompletion=>"open_trace_completion",
            PlanRequirement::CurvedEdgeSubdivision=>"curved_edge_subdivision",
            PlanRequirement::MultipleSplitsSameEdge=>"multiple_splits_same_edge",
            PlanRequirement::TangencyClassification=>"tangency_classification",
            PlanRequirement::BoundaryReconciliation=>"boundary_reconciliation",
        }).collect::<Vec<_>>();
        let position=|v|match v {
            rustsolid::VolumePosition::Inside=>"inside",
            rustsolid::VolumePosition::Outside=>"outside",
            rustsolid::VolumePosition::Boundary=>"boundary",
        };
        let hold=plan.hold.map(|r|match r{
            RegionHold::AmbiguousIntersection=>"ambiguous_intersection",
            RegionHold::Tangency=>"tangency",
            RegionHold::BoundaryOverlap=>"boundary_overlap",
            RegionHold::NearlyTangent=>"nearly_tangent",
            RegionHold::NonSeparatingTrace=>"non_separating_trace",
            RegionHold::UnresolvedProbe=>"unresolved_probe",
        });
        let samples=plan.samples.iter().map(|item|json!({
            "arc_id":item.arc.0,
            "operand":match item.operand{
                GraphOperand::First=>"first",GraphOperand::Second=>"second"},
            "face_id":item.face.0,
            "fraction":item.fraction,
            "curve_point":item.point.array(),
            "positive_probe":item.positive.array(),
            "negative_probe":item.negative.array(),
            "positive_other":position(item.positive_other),
            "negative_other":position(item.negative_other),
            "keep_positive":item.keep_positive,
            "keep_negative":item.keep_negative,
            "reverse_orientation":item.reverse_orientation,
        })).collect::<Vec<_>>();
        Ok(json!({
            "kind":"boolean_region_plan",
            "mode":mode,
            "state":match plan.state{
                GraphRegionState::NoCut=>"no_cut",
                GraphRegionState::ClassifiedPreview=>"classified_preview",
                GraphRegionState::Held=>"held",
            },
            "hold_reason":hold,
            "graph_state":graph.kind(),
            "source_faces":[fa.0,fb.0],
            "intersection_arcs":graph.arcs.len(),
            "region_samples":samples,
            "topology_intent":{
                "state":match topology.state {
                    PlanState::NoChange=>"no_change",
                    PlanState::PreviewOnly=>"preview_only",
                    PlanState::Held=>"held",
                },
                "proposed_vertices":topology.vertices.len(),
                "proposed_arcs":topology.arcs.len(),
                "face_changes":topology.face_changes.len(),
                "remaining_requirements":outstanding,
                "can_apply":false,
            },
            "topology_committed":false,
            "can_apply":false,
        }))
    })();
    match result {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(error.to_string()),data:None},
    }
}


#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
enum CarrierCommandInput {
    ReparameterizeEdge {edge_id:u32,scale:f64,offset:f64},
    ReparameterizeTrim {coedge_id:u32,scale:f64,offset:f64},
    RebaseFaceUv {face_id:u32,du:f64,dv:f64},
}
enum CarrierResolvedCommand {
    ReparameterizeEdge {handle:CarrierHandle,scale:f64,offset:f64},
    ReparameterizeTrim {handle:CarrierHandle,scale:f64,offset:f64},
    RebaseFaceUv {handle:CarrierHandle,du:f64,dv:f64},
}
impl CarrierCommandInput {
    fn resolve(self,editor:&BrepCarrierEditor)->Result<CarrierResolvedCommand,GeometryError>{
        let resolve=|owner|editor.handle(owner).ok_or_else(||
            GeometryError::InvalidEdit(format!("unknown carrier owner: {owner:?}")));
        Ok(match self {
            Self::ReparameterizeEdge{edge_id,scale,offset}=>
                CarrierResolvedCommand::ReparameterizeEdge{
                    handle:resolve(CarrierOwner::Edge(EdgeId(edge_id)))?,scale,offset},
            Self::ReparameterizeTrim{coedge_id,scale,offset}=>
                CarrierResolvedCommand::ReparameterizeTrim{
                    handle:resolve(CarrierOwner::Coedge(rustsolid::CoedgeId(coedge_id)))?,scale,offset},
            Self::RebaseFaceUv{face_id,du,dv}=>
                CarrierResolvedCommand::RebaseFaceUv{
                    handle:resolve(CarrierOwner::Face(FaceId(face_id)))?,du,dv},
        })
    }
}
impl CarrierResolvedCommand {
    fn stage(self,tx:&mut BrepCarrierTransaction<'_>)->Result<(),GeometryError>{
        match self {
            Self::ReparameterizeEdge{handle,scale,offset}=>tx.reparameterize_edge(handle,scale,offset),
            Self::ReparameterizeTrim{handle,scale,offset}=>tx.reparameterize_trim(handle,scale,offset),
            Self::RebaseFaceUv{handle,du,dv}=>tx.rebase_face_uv(handle,du,dv),
        }
    }
}

#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
enum StructuralCommandInput {
    AddThroughHole {center:[f64;2],radius:f64},
    RemoveThroughHole {slot:u32},
    ResizeThroughHole {slot:u32,center:[f64;2],radius:f64},
    ReparameterizeHoleEdge {slot:u32,role:HoleEdgeRole,scale:f64,offset:f64},
    ReparameterizeHoleTrim {slot:u32,role:HoleTrimRole,scale:f64,offset:f64},
    RebaseHoleFaceUv {slot:u32,du:f64,dv:f64},
    RebaseBaseFaceUv {face_id:u32,du:f64,dv:f64},
    ReparameterizeBaseEdge {edge_id:u32,scale:f64,offset:f64},
    ReparameterizeBaseTrim {coedge_id:u32,scale:f64,offset:f64},
}
impl StructuralCommandInput {
    fn stage(self,tx:&mut BrepStructuralTransaction<'_>)
        ->Result<(),GeometryError>
    {
        let named=|tx:&BrepStructuralTransaction<'_>,slot|tx.handle(slot)
            .ok_or_else(||GeometryError::InvalidEdit(format!(
                "unknown staged/current hole slot {slot}")));
        match self {
            Self::AddThroughHole{center,radius}=>{
                tx.add_through_hole(HoleParameters{
                    center_x:center[0],center_z:center[1],radius})?;
            }
            Self::RemoveThroughHole{slot}=>{
                let h=named(tx,slot)?;tx.remove_through_hole(h)?;
            }
            Self::ResizeThroughHole{slot,center,radius}=>{
                let h=named(tx,slot)?;
                tx.resize_through_hole(h,HoleParameters{
                    center_x:center[0],center_z:center[1],radius})?;
            }
            Self::ReparameterizeHoleEdge{slot,role,scale,offset}=>{
                let hole=named(tx,slot)?;
                tx.reparameterize_edge(StructuralEdgeRef::Hole{hole,role},scale,offset)?;
            }
            Self::ReparameterizeHoleTrim{slot,role,scale,offset}=>{
                let hole=named(tx,slot)?;
                tx.reparameterize_trim(StructuralTrimRef::Hole{hole,role},scale,offset)?;
            }
            Self::RebaseHoleFaceUv{slot,du,dv}=>{
                let hole=named(tx,slot)?;
                tx.rebase_face_uv(StructuralFaceRef::HoleSide(hole),du,dv)?;
            }
            Self::RebaseBaseFaceUv{face_id,du,dv}=>{
                tx.rebase_face_uv(StructuralFaceRef::Base(FaceId(face_id)),du,dv)?;
            }
            Self::ReparameterizeBaseEdge{edge_id,scale,offset}=>{
                tx.reparameterize_edge(StructuralEdgeRef::Base(EdgeId(edge_id)),scale,offset)?;
            }
            Self::ReparameterizeBaseTrim{coedge_id,scale,offset}=>{
                tx.reparameterize_trim(StructuralTrimRef::Base(
                    rustsolid::CoedgeId(coedge_id)),scale,offset)?;
            }
        }
        Ok(())
    }
}
fn structural_binding_json(x:rustsolid::StructuralHoleBinding)->Value{
    json!({
        "slot":x.handle.slot,"generation":x.handle.generation,
        "editor_token":x.handle.editor_token,
        "center":[x.params.center_x,x.params.center_z],"radius":x.params.radius,
        "side_face_id":x.topology.side_face.0,
        "edge_ids":[x.topology.bottom_edge.0,x.topology.top_edge.0,
            x.topology.seam_edge.0],
        "loop_ids":[x.topology.bottom_loop.0,x.topology.top_loop.0,
            x.topology.side_loop.0],
    })
}
fn evaluate_structural_transaction(
    source:SolidSource,batches:Vec<Vec<StructuralCommandInput>>,
    undo:usize,redo:usize,display_segments:Option<usize>,
)->Response{
    let result=(||->Result<(Value,Option<SolidResponse>),GeometryError>{
        if batches.len()>64 || undo>64 || redo>64 {
            return Err(GeometryError::InvalidEdit(
                "structural transaction exceeds 64 batch/transition budget".into()));
        }
        let body=source.build()?;
        let mut editor=BrepStructuralEditor::new(body)?;
        let mut commits=Vec::new();
        for batch in batches {
            if batch.is_empty() || batch.len()>256 {
                return Err(GeometryError::InvalidEdit(
                    "structural batch requires between 1 and 256 commands".into()));
            }
            let mut tx=editor.begin_transaction();
            for cmd in batch {cmd.stage(&mut tx)?;}
            let report=tx.commit()?;
            commits.push(json!({
                "revision_before":report.revision_before,
                "revision_after":report.revision_after,
                "command_count":report.command_count,
                "created":report.created.into_iter()
                    .map(structural_binding_json).collect::<Vec<_>>(),
                "retired":report.retired.into_iter()
                    .map(structural_binding_json).collect::<Vec<_>>(),
                "remapped":report.remapped.into_iter()
                    .map(|m|json!({
                        "slot":m.handle.slot,"generation":m.handle.generation,
                        "side_face_before":m.before.side_face.0,
                        "side_face_after":m.after.side_face.0,
                        "bottom_edge_before":m.before.bottom_edge.0,
                        "bottom_edge_after":m.after.bottom_edge.0,
                    })).collect::<Vec<_>>(),
            }));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            let value=editor.undo()?.ok_or_else(||GeometryError::InvalidEdit(
                "structural undo exceeds committed history".into()))?;
            transitions.push(json!({"kind":"undo","revision_after":value.revision_after}));
        }
        for _ in 0..redo {
            let value=editor.redo()?.ok_or_else(||GeometryError::InvalidEdit(
                "structural redo exceeds available history".into()))?;
            transitions.push(json!({"kind":"redo","revision_after":value.revision_after}));
        }
        editor.validate()?;
        let display=display_segments.map(|segments|editor.tessellate(segments))
            .transpose()?;
        let model=editor.model();let topology=model.summary();
        let holes=editor.bindings().into_iter()
            .map(structural_binding_json).collect::<Vec<_>>();
        let display_quality=display.as_ref().map(|m|json!({
            "segments_per_hole":m.circle_segments,
            "facet_triangles":m.mesh.triangles.len(),
            "bottom_cap_region":{
                "face_id":m.bottom_region.face.0,
                "loop_count":m.bottom_region.loop_count,
                "triangles":m.bottom_region.triangles,
                "polygonal_area":m.bottom_region.polygonal_area
            },
            "top_cap_region":{
                "face_id":m.top_region.face.0,
                "loop_count":m.top_region.loop_count,
                "triangles":m.top_region.triangles,
                "polygonal_area":m.top_region.polygonal_area
            },
            "faceted_volume":m.faceted_volume,
            "faceted_surface_area":m.faceted_surface_area,
            "exact_analytic_volume":m.analytic_volume,
            "volume_faceting_error":m.faceting_volume_error,
            "closed_oriented_manifold_certified":true
        }));
        let solid=display.map(|m|SolidResponse {
            vertices:m.mesh.vertices.iter().map(|v|v.array()).collect(),
            faces:m.mesh.triangles,
            triangle_face_ids:m.mesh.triangle_faces.iter().map(|f|f.0).collect(),
            topology:Topology {
                vertices:topology.vertices,edges:topology.edges,
                coedges:topology.coedges,loops:topology.loops,
                faces:topology.faces,shells:topology.shells,
                euler_characteristic:topology.euler_characteristic,
            },
            mass:Mass {volume:model.mass.volume,surface_area:model.mass.surface_area,
                centroid:model.mass.centroid.array()},
            bbox:[model.bbox.min.array(),model.bbox.max.array()],
            edit_report:None,history:None,
            analytic_brep:Some(json!({
                "representation":"analytic_multiple_through_holes",
                "exact_analytic_surfaces":true,
                "cylindrical_faces":topology.cylindrical_faces,
                "circular_edges":topology.circular_edges,
                "genus":topology.genus,
                "mesh_is_display_only":true
            })),
            brep:None,
        });
        let data=json!({
            "kind":"structural_brep","state":"committed",
            "representation":match model.origin{
                BrepOrigin::Polyhedral=>"polyhedral",
                BrepOrigin::AnalyticCylinder=>"analytic_cylinder",
                BrepOrigin::AnalyticThroughHole=>"analytic_through_hole",
                BrepOrigin::AnalyticMultipleThroughHoles=>"analytic_multiple_through_holes",
                BrepOrigin::AnalyticBlindHole=>"analytic_blind_hole",
            },
            "revision":editor.revision(),
            "editor_token":editor.editor_token(),
            "hole_count":editor.hole_count(),
            "holes":holes,
            "topology":{"vertices":topology.vertices,"edges":topology.edges,
                "coedges":topology.coedges,"faces":topology.faces,
                "loops":topology.loops,"shells":topology.shells,
                "euler_characteristic":topology.euler_characteristic,
                "genus":topology.genus,"cylindrical_faces":topology.cylindrical_faces,
                "circular_edges":topology.circular_edges,"seam_edges":topology.seam_edges,
                "inner_loops":model.loops.iter().filter(|l|
                    l.role==rustsolid::LoopRole::Inner).count()},
            "carrier_counts":{"curves3":model.geometry.curves3.len(),
                "curves2":model.geometry.curves2.len(),
                "surfaces":model.geometry.surfaces.len()},
            "mass":{"volume":model.mass.volume,
                "surface_area":model.mass.surface_area,
                "centroid":model.mass.centroid.array()},
            "undo_depth":editor.undo_depth(),"redo_depth":editor.redo_depth(),
            "commits":commits,"transitions":transitions,
            "topology_committed":true,
            "display_mesh_available":true,
            "display_mesh_included":solid.is_some(),
            "display_quality":display_quality,
            "generalized_topology_delta_can_apply":false,
            "note":"Independent genus-N analytic B-rep. Dense topology IDs may change on hole removal; use process-local feature handles within one transaction session"
        });
        Ok((data,solid))
    })();
    match result {
        Ok((data,solid))=>Response{schema_version:"geometry.v1",ok:true,
            solid,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(error.to_string()),data:None},
    }
}

fn evaluate_carrier_transaction(
    source:BrepInput,batches:Vec<Vec<CarrierCommandInput>>,undo:usize,redo:usize,
)->Response{
    let response=(||->Result<Value,GeometryError>{
        if batches.len()>64 || undo>64 || redo>64{
            return Err(GeometryError::InvalidEdit(
                "carrier session exceeds 64 batches/transitions".into()));
        }
        let (source,_)=source.prepare()?;
        let mut editor=BrepCarrierEditor::from_body(&source)?;
        let mut commits=Vec::new();
        for batch in batches {
            if batch.is_empty() || batch.len()>256{
                return Err(GeometryError::InvalidEdit(
                    "carrier batch needs between 1 and 256 edits".into()));
            }
            let commands=batch.into_iter()
                .map(|input|input.resolve(&editor))
                .collect::<Result<Vec<_>,_>>()?;
            let mut tx=editor.begin_transaction();
            for command in commands {command.stage(&mut tx)?;}
            let committed=tx.commit()?;
            let changed=committed.changed.iter().map(|binding|{
                let (kind,id)=match binding.handle.owner {
                    CarrierOwner::Edge(e)=>("edge",e.0),
                    CarrierOwner::Coedge(c)=>("coedge",c.0),
                    CarrierOwner::Face(f)=>("face",f.0),
                };
                let (carrier_kind,carrier_id)=match binding.id{
                    CarrierId::Curve3(c)=>("curve3",c.0),
                    CarrierId::Curve2(c)=>("curve2",c.0),
                    CarrierId::Surface(c)=>("surface",c.0),
                };
                json!({"owner_kind":kind,"owner_id":id,
                    "carrier_kind":carrier_kind,"carrier_id":carrier_id,
                    "generation":binding.handle.generation})
            }).collect::<Vec<_>>();
            commits.push(json!({"revision_before":committed.revision_before,
                "revision_after":committed.revision_after,"changed":changed}));
        }
        let mut transitions=Vec::new();
        for _ in 0..undo {
            let t=editor.undo()?.ok_or_else(||GeometryError::InvalidEdit(
                "carrier undo exceeds history".into()))?;
            transitions.push(json!({"kind":"undo","revision_after":t.revision_after}));
        }
        for _ in 0..redo {
            let t=editor.redo()?.ok_or_else(||GeometryError::InvalidEdit(
                "carrier redo exceeds history".into()))?;
            transitions.push(json!({"kind":"redo","revision_after":t.revision_after}));
        }
        editor.validate()?;
        let graph=editor.model();
        let summary=graph.summary();
        Ok(json!({
            "kind":"carrier_transaction","revision":editor.revision(),
            "editor_token":editor.editor_token(),
            "undo_depth":editor.undo_depth(),"redo_depth":editor.redo_depth(),
            "carrier_counts":{"curves3":graph.geometry.curves3.len(),
                "curves2":graph.geometry.curves2.len(),
                "surfaces":graph.geometry.surfaces.len()},
            "topology":{"vertices":summary.vertices,"edges":summary.edges,
                "coedges":summary.coedges,"faces":summary.faces,
                "loops":summary.loops,"shells":summary.shells,
                "euler_characteristic":summary.euler_characteristic,
                "genus":summary.genus},
            "mass":{"volume":graph.mass.volume,"surface_area":graph.mass.surface_area,
                "centroid":graph.mass.centroid.array()},
            "edge_domains":graph.geometry.edges.iter().map(|g|match g.domain {
                ParameterRange::Bounded{start,end}=>json!([start,end]),
                ParameterRange::Unbounded=>json!(null)
            }).collect::<Vec<_>>(),
            "trim_domains":graph.geometry.coedges.iter().map(|g|match g.domain {
                ParameterRange::Bounded{start,end}=>json!([start,end]),
                ParameterRange::Unbounded=>json!(null)
            }).collect::<Vec<_>>(),
            "commits":commits,"transitions":transitions,
            "shape_preserving":true,"topology_mutated":false,
            "note":"Carrier handles and editor_token are process-local; request does not mutate the source B-rep"
        }))
    })();
    match response {
        Ok(data)=>Response{schema_version:"geometry.v1",ok:true,
            solid:None,error:None,data:Some(data)},
        Err(error)=>Response{schema_version:"geometry.v1",ok:false,
            solid:None,error:Some(error.to_string()),data:None},
    }
}

fn evaluate_boolean(left:BrepInput,right:BrepInput,mode:BooleanOperation)->Response{
    let result=(||->Result<rustsolid::BooleanResult,GeometryError>{
        let (left,_)=left.prepare()?;
        let (right,_)=right.prepare()?;
        BooleanKernel::apply(&left,&right,mode)
            .map_err(|error|GeometryError::InvalidEdit(error.to_string()))
    })();
    match result {
        Err(error)=>Response {
            schema_version:"geometry.v1",ok:false,solid:None,
            error:Some(error.to_string()),data:None,
        },
        Ok(result)=>{
            let lineage=result.face_sources.iter().map(|f|json!({
                "result_face_id":f.result_face.0,
                "operand":match f.operand {
                    BooleanOperand::Left=>"left",BooleanOperand::Right=>"right"},
                "source_face_id":f.source_face.0,
                "reversed":f.reversed,
            })).collect::<Vec<_>>();
            match result.outcome {
                BooleanOutcome::Empty=>Response{
                    schema_version:"geometry.v1",ok:true,solid:None,error:None,
                    data:Some(json!({"kind":"boolean","mode":mode,"state":"empty",
                        "face_sources":lineage})),
                },
                BooleanOutcome::Analytic(body)=>{
                    let model=body.model();
                    let stats=model.summary();
                    let mesh=match body.tessellate(64){
                        Ok(mesh)=>mesh,
                        Err(error)=>return Response{schema_version:"geometry.v1",
                            ok:false,solid:None,data:None,error:Some(error.to_string())},
                    };
                    Response {
                        schema_version:"geometry.v1",ok:true,error:None,
                        data:Some(json!({"kind":"boolean","mode":mode,"state":"solid",
                            "representation":"analytic_cylindrical_through_hole",
                            "face_sources":lineage,"genus":stats.genus,
                            "topology_committed":true})),
                        solid:Some(SolidResponse{
                            vertices:mesh.vertices.iter().map(|p|p.array()).collect(),
                            faces:mesh.triangles,
                            triangle_face_ids:mesh.triangle_faces.iter().map(|id|id.0).collect(),
                            topology:Topology{vertices:stats.vertices,edges:stats.edges,
                                coedges:stats.coedges,loops:stats.loops,faces:stats.faces,
                                shells:stats.shells,euler_characteristic:stats.euler_characteristic},
                            mass:Mass{volume:model.mass.volume,surface_area:model.mass.surface_area,
                                centroid:model.mass.centroid.array()},
                            bbox:[model.bbox.min.array(),model.bbox.max.array()],
                            edit_report:None,history:None,
                            analytic_brep:Some(json!({"cylindrical_faces":stats.cylindrical_faces,
                                "circular_edges":stats.circular_edges,
                                "inner_loops":model.loops.iter().filter(|l|
                                    l.role==rustsolid::LoopRole::Inner).count(),
                                "genus":stats.genus,
                                "display_mesh_is_not_boundary_definition":true})),
                            brep:None,
                        }),
                    }
                }
                BooleanOutcome::BlindHole(body)=>{
                    let opening=match body.opening_trace(){
                        Ok(graph)=>graph,
                        Err(error)=>return Response{schema_version:"geometry.v1",
                            ok:false,solid:None,data:None,error:Some(error.to_string())},
                    };
                    let stats=body.model().summary();
                    let shape=body.model();
                    let display=match body.tessellate(64){
                        Ok(mesh)=>mesh,
                        Err(error)=>return Response{schema_version:"geometry.v1",
                            ok:false,solid:None,data:None,error:Some(error.to_string())},
                    };
                    Response {
                        schema_version:"geometry.v1",ok:true,error:None,
                        data:Some(json!({
                            "kind":"boolean","mode":mode,"state":"solid",
                            "representation":"analytic_cylindrical_blind_hole",
                            "face_sources":lineage,"genus":stats.genus,
                            "opening_face_id":display.opening_face.0,
                            "floor_face_id":display.floor_face.0,
                            "side_face_id":display.side_face.0,
                            "opening_inner_loops":1,
                            "source_opening_trace":{
                                "state":opening.kind(),
                                "arcs":opening.arcs.len(),
                                "closed":opening.arcs.iter().all(|arc|arc.closed),
                                "contact":match opening.arcs[0].contact {
                                    TrimContact::Crossing=>"crossing",
                                    TrimContact::Tangent=>"tangent",
                                    TrimContact::BoundaryOverlap=>"boundary_overlap",
                                }
                            },
                            "faceted_volume":display.faceted_volume,
                            "faceted_surface_area":display.faceted_surface_area,
                            "volume_faceting_error":display.faceting_volume_error,
                            "closed_oriented_display_mesh_certified":true,
                            "topology_committed":true,
                        })),
                        solid:Some(SolidResponse{
                            vertices:display.mesh.vertices.iter().map(|p|p.array()).collect(),
                            faces:display.mesh.triangles,
                            triangle_face_ids:display.mesh.triangle_faces.iter()
                                .map(|id|id.0).collect(),
                            topology:Topology{
                                vertices:stats.vertices,edges:stats.edges,
                                coedges:stats.coedges,loops:stats.loops,
                                faces:stats.faces,shells:stats.shells,
                                euler_characteristic:stats.euler_characteristic,
                            },
                            mass:Mass{volume:shape.mass.volume,
                                surface_area:shape.mass.surface_area,
                                centroid:shape.mass.centroid.array()},
                            bbox:[shape.bbox.min.array(),shape.bbox.max.array()],
                            edit_report:None,history:None,
                            analytic_brep:Some(json!({
                                "cylindrical_faces":stats.cylindrical_faces,
                                "plane_faces":stats.plane_faces,
                                "circular_edges":stats.circular_edges,
                                "inner_loops":shape.loops.iter().filter(|l|
                                    l.role==rustsolid::LoopRole::Inner).count(),
                                "floor_face_id":display.floor_face.0,
                                "genus":stats.genus,
                                "display_mesh_is_not_boundary_definition":true,
                            })),
                            brep:None,
                        })
                    }
                }
                BooleanOutcome::Solid(body)=>Response{
                    schema_version:"geometry.v1",ok:true,error:None,
                    data:Some(json!({"kind":"boolean","mode":mode,"state":"solid",
                        "representation":"exact_planar_box_arrangement",
                        "face_sources":lineage})),
                    solid:Some(SolidResponse{
                        vertices:body.mesh.vertices.iter().map(|p|p.array()).collect(),
                        faces:body.mesh.triangles,
                        triangle_face_ids:body.mesh.triangle_faces.iter().map(|id|id.0).collect(),
                        topology:Topology{
                            vertices:body.vertices.len(),edges:body.edges.len(),
                            coedges:body.coedges.len(),loops:body.loops.len(),
                            faces:body.faces.len(),shells:body.shells.len(),
                            euler_characteristic:body.vertices.len() as isize -
                                body.edges.len() as isize + body.faces.len() as isize,
                        },
                        mass:Mass {volume:body.mass.volume,surface_area:body.mass.surface_area,
                            centroid:body.mass.centroid.array()},
                        bbox:[body.bbox.min.array(),body.bbox.max.array()],
                        edit_report:None,history:None,analytic_brep:None,brep:None,
                    }),
                }
            }
        }
    }
}

fn process(input: &str) -> Response {
    let parsed: Result<GeometryRequest, _> = serde_json::from_str(input);
    let result: Result<(Solid, Option<EditReport>, Option<Value>), GeometryError> = match parsed {
        Ok(GeometryRequest::ExtrudeProfile { profile, height, tolerance }) => {
            SolidSource::ExtrudeProfile { profile, height, tolerance }.build().map(|body|(body,None,None))
        }
        Ok(GeometryRequest::Block { origin, width, height, depth, tolerance }) => {
            SolidSource::Block { origin, width, height, depth, tolerance }.build().map(|body|(body,None,None))
        }
        Ok(GeometryRequest::InspectBrep {source}) => return inspect_brep(source),
        Ok(GeometryRequest::IntersectFaces{first,second,tolerance,samples}) =>
            return evaluate_face_intersection(first,second,tolerance.unwrap_or_default(),samples.unwrap_or(9)),
        Ok(GeometryRequest::IntersectionGraph{first,second,tolerance,samples}) =>
            return evaluate_intersection_graph(first,second,tolerance.unwrap_or_default(),samples.unwrap_or(9)),
        Ok(GeometryRequest::TopologyDeltaPlan{first,second,tolerance,samples}) =>
            return evaluate_topology_delta_plan(first,second,tolerance.unwrap_or_default(),samples.unwrap_or(9)),
        Ok(GeometryRequest::ArenaTransaction{source,batches,undo,redo}) =>
            return evaluate_arena_transaction(source,batches,undo,redo),
        Ok(GeometryRequest::Boolean{left,right,mode}) =>
            return evaluate_boolean(left,right,mode),
        Ok(GeometryRequest::CarrierTransaction{source,batches,undo,redo}) =>
            return evaluate_carrier_transaction(source,batches,undo,redo),
        Ok(GeometryRequest::StructuralTransaction{source,batches,undo,redo,display_segments}) =>
            return evaluate_structural_transaction(source,batches,undo,redo,display_segments),
        Ok(GeometryRequest::BooleanRegionPlan{first,second,mode,tolerance}) =>
            return evaluate_boolean_region_plan(first,second,mode,
                tolerance.unwrap_or_default()),
        Ok(GeometryRequest::GeometryQuery {source,query}) =>
            return evaluate_geometry_query(source,&query),
        Ok(GeometryRequest::ProbeKernel {seed}) =>
            return evaluate_boundary_probes(seed.unwrap_or(13)),
        Ok(GeometryRequest::RegressionCorpus) => return evaluate_regression_corpus(),
        Ok(GeometryRequest::IntersectSurfaces {first,second,tolerance,samples}) =>
            return evaluate_surface_intersection(first,second,tolerance.unwrap_or_default(),samples.unwrap_or(9)),
        Ok(GeometryRequest::AnalyticCylinder {origin,radius,height,segments,tolerance}) => {
            return process_analytic_cylinder(origin,radius,height,
                segments.unwrap_or(64),tolerance.unwrap_or_default());
        }
        Ok(GeometryRequest::EditSolid { source, edits }) => {
            if edits.len() > 256 {
                Err(GeometryError::InvalidEdit("maximum 256 edits per request".into()))
            } else {
                source.build().and_then(|mut body| {
                    let report=body.edit_atomic(|tx| {
                        for edit in edits {
                            match edit {
                                EditRequest::SplitEdge { edge_id, fraction } => {
                                    tx.split_edge(EdgeId(edge_id),fraction)?;
                                }
                                EditRequest::SplitFace { face_id, start_vertex_id, end_vertex_id } => {
                                    tx.split_face(FaceId(face_id),VertexId(start_vertex_id),VertexId(end_vertex_id))?;
                                }
                            }
                        }
                        Ok(())
                    })?;
                    Ok((body,Some(report),None))
                })
            }
        }
        Ok(GeometryRequest::CommandHistory { source, feature_key, batches, undo, redo }) => {
            if batches.len() > 256 || undo > 256 || redo > 256 {
                Err(GeometryError::InvalidEdit("history request exceeds 256-item limit".into()))
            } else {
                source.build().and_then(|solid| {
                    let solid = match feature_key {
                        Some(key) => solid.with_feature_key(&key)?,
                        None => solid,
                    };
                    let mut session = CommandHistory::new(solid)?;
                    for batch in batches { session.execute(batch)?; }
                    for _ in 0..undo {
                        if session.undo()?.is_none() {
                            return Err(GeometryError::InvalidEdit("undo depth exceeds history".into()));
                        }
                    }
                    for _ in 0..redo {
                        if session.redo()?.is_none() {
                            return Err(GeometryError::InvalidEdit("redo depth exceeds history".into()));
                        }
                    }
                    let state = json!({
                        "revision": session.solid().revision,
                        "undo_depth": session.undo_depth(), "redo_depth": session.redo_depth(),
                        "archived_batches": session.archived_depth(),
                        "applied_batches": session.command_log(),
                    });
                    Ok((session.into_solid(), None, Some(state)))
                })
            }
        }
        Err(error) => return Response {
            schema_version: "geometry.v1", ok: false, solid: None,
            error: Some(format!("invalid JSON request: {error}")),data:None,
        },
    };
    match result {
        Ok((body,report,history_state)) => Response {
            schema_version: "geometry.v1", ok: true,
            solid: Some(SolidResponse {
                vertices: body.mesh.vertices.iter().map(|v| v.array()).collect(),
                faces: body.mesh.triangles,
                triangle_face_ids: body.mesh.triangle_faces.iter().map(|id| id.0).collect(),
                topology: Topology {
                    vertices: body.vertices.len(), edges: body.edges.len(),
                    coedges: body.coedges.len(), loops: body.loops.len(), shells: body.shells.len(),
                    faces: body.faces.len(), euler_characteristic: body.vertices.len() as isize - body.edges.len() as isize + body.faces.len() as isize,
                },
                mass: Mass { volume: body.mass.volume, surface_area: body.mass.surface_area,
                    centroid: body.mass.centroid.array() },
                bbox: [body.bbox.min.array(), body.bbox.max.array()],
                edit_report: report.map(report_json),
                history: history_state,
                analytic_brep: None,
                brep: None,
            }),
            error: None,data:None,
        },
        Err(error) => Response { schema_version: "geometry.v1", ok: false,
            solid: None, error: Some(error.to_string()), data:None },
    }
}

fn main() {
    let mut input = String::new();
    if let Err(error) = io::stdin().read_to_string(&mut input) {
        eprintln!("read failed: {error}");
        std::process::exit(1);
    }
    let response = process(&input);
    match serde_json::to_string(&response) {
        Ok(output) => println!("{output}"),
        Err(error) => { eprintln!("serialization failed: {error}"); std::process::exit(1); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_block_roundtrip() {
        let response = process(r#"{"operation":"block","origin":[1,2,3],"width":2,"height":3,"depth":4}"#);
        assert!(response.ok);
        let solid = response.solid.unwrap();
        assert_eq!(solid.topology.euler_characteristic, 2);
        assert_eq!(solid.topology.coedges, 24);
        assert_eq!(solid.topology.loops, 6);
        assert_eq!(solid.topology.shells, 1);
        assert_eq!(solid.mass.volume, 24.0);
        assert_eq!(solid.triangle_face_ids.len(), solid.faces.len());
        assert!(solid.edit_report.is_none());
    }

    #[test]
    fn malformed_request_reports_error() {
        assert!(!process(r#"{"operation":"boolean"}"#).ok);
        assert!(!process("garbage").ok);
    }

    #[test]
    fn per_request_tolerance_is_accepted_and_validated() {
        let request = r#"{"operation":"block","origin":[0,0,0],"width":0.00000001,"height":1,"depth":1,
            "tolerance":{"absolute_length":1e-12,"relative_length":1e-13}}"#;
        assert!(process(request).ok);
        let invalid = r#"{"operation":"block","origin":[0,0,0],"width":1,"height":1,"depth":1,
            "tolerance":{"absolute_length":-2}}"#;
        let response = process(invalid);
        assert!(!response.ok);
        assert!(response.error.unwrap().contains("invalid tolerance"));
    }

    #[test]
    fn json_edits_apply_atomically_and_return_provenance() {
        let response=process(r#"{"operation":"edit_solid",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "edits":[{"kind":"split_edge","edge_id":0,"fraction":0.5},
                     {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":2}]}"#);
        assert!(response.ok,"{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!(solid.topology.vertices,9);
        assert_eq!(solid.topology.faces,7);
        assert_eq!(solid.mass.volume,24.0);
        let report=solid.edit_report.unwrap();
        assert_eq!(report["revision_before"],0);
        assert_eq!(report["revision_after"],1);
        assert_eq!(report["changes"].as_array().unwrap().len(),2);
        assert_eq!(report["changes"][0]["kind"],"split_edge");
        assert_eq!(report["changes"][1]["kind"],"split_face");
    }

    #[test]
    fn json_failed_batch_does_not_return_partial_solid() {
        let response=process(r#"{"operation":"edit_solid",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "edits":[{"kind":"split_edge","edge_id":0,"fraction":0.5},
                     {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":1}]}"#);
        assert!(!response.ok);
        assert!(response.solid.is_none());
    }
    #[test]
    fn json_command_history_undo_redo() {
        let response = process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "feature_key":"part/sample",
            "batches":[
              {"commands":[{"kind":"split_edge","edge":"part/sample/edge/0","fraction":0.5}]},
              {"commands":[{"kind":"split_face","face":"part/sample/face/bottom",
                "start":"part/sample/vertex/0","end":"part/sample/vertex/2"}]}
            ],
            "undo":2,"redo":1
        }"#);
        assert!(response.ok, "{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!((solid.topology.vertices,solid.topology.faces),(9,6));
        let history=solid.history.unwrap();
        assert_eq!(history["revision"],5);
        assert_eq!(history["undo_depth"],1);
        assert_eq!(history["redo_depth"],1);
        assert_eq!(history["applied_batches"].as_array().unwrap().len(),1);
    }

    #[test]
    fn json_command_history_rejects_excess_undo() {
        let response = process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[{"commands":[{"kind":"split_edge","edge":"primitive/edge/0","fraction":0.5}]}],
            "undo":2
        }"#);
        assert!(!response.ok);
        assert!(response.solid.is_none());
    }

    #[test]
    fn json_command_history_independent_kill_and_redo() {
        let result=process(r#"{
            "operation":"command_history",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "feature_key":"part/feature",
            "batches":[
                {"commands":[{"kind":"make_edge_vertex","edge":"part/feature/edge/0","fraction":0.5}]},
                {"commands":[{"kind":"kill_edge_vertex","vertex":"part/feature/edge/0/edit-1/vertex"}]}
            ],
            "undo":1, "redo":1
        }"#);
        assert!(result.ok,"{:?}",result.error);
        let solid=result.solid.unwrap();
        assert_eq!(solid.topology.vertices,8);
        assert_eq!(solid.topology.edges,12);
        assert_eq!(solid.mass.volume,24.0);
        let history=solid.history.unwrap();
        assert_eq!(history["revision"],4);
        assert_eq!(history["applied_batches"].as_array().unwrap().len(),2);
    }

    #[test]
    fn json_analytic_cylinder_exposes_exact_brep_and_display_mesh() {
        let result=process(r#"{"operation":"analytic_cylinder",
            "origin":[0,0,0],"radius":2,"height":5,"segments":32}"#);
        assert!(result.ok,"{:?}",result.error);
        let solid=result.solid.unwrap();
        assert_eq!((solid.topology.vertices,solid.topology.edges,solid.topology.coedges,
            solid.topology.loops,solid.topology.faces,solid.topology.shells),(2,3,6,3,3,1));
        assert_eq!(solid.topology.euler_characteristic,2);
        assert_eq!(solid.faces.len(),128);
        assert_eq!(solid.triangle_face_ids.len(),solid.faces.len());
        assert!((solid.mass.volume-std::f64::consts::PI*20.0).abs()<1e-10);
        let analytic=solid.analytic_brep.unwrap();
        assert_eq!(analytic["seam_edge_id"],2);
        assert_eq!(analytic["seam_coedge_ids"],json!([4,5]));
        assert_eq!(analytic["u_periodic"],true);
        assert_eq!(analytic["coedges"].as_array().unwrap().len(),6);
    }

    #[test]
    fn json_analytic_cylinder_rejects_invalid_parameters() {
        for request in [
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":0,"height":1}"#,
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":1,"height":2,"segments":3}"#,
            r#"{"operation":"analytic_cylinder","origin":[0,0,0],"radius":1,"height":2,"tolerance":{"absolute_length":-1}}"#,
        ] {
            let response=process(request);
            assert!(!response.ok,"{request}");
            assert!(response.solid.is_none());
        }
    }

    #[test]
    fn unified_brep_query_keeps_analytic_curves_and_shared_seam() {
        let response=process(r#"{"operation":"inspect_brep",
            "source":{"kind":"analytic_cylinder",
            "origin":[0,0,0],"radius":2,"height":5,"segments":16}}"#);
        assert!(response.ok,"{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!(solid.topology.faces,3);
        assert_eq!(solid.faces.len(),64);
        let brep=solid.brep.unwrap();
        assert_eq!(brep["representation"],"analytic_cylinder");
        assert_eq!(brep["geometry"]["face_surfaces"],json!(["plane","plane","cylinder"]));
        assert_eq!(brep["geometry"]["edge_curves"],json!(["circle","circle","line"]));
        assert_eq!(brep["geometry"]["seam_edges"],json!([2]));
        assert_eq!(brep["geometry"]["closed_edges"],json!([0,1]));
        assert_eq!(brep["topology"]["euler_characteristic"],2);
        assert_eq!(brep["coedges"].as_array().unwrap().len(),6);
    }

    #[test]
    fn unified_brep_query_accepts_polyhedra_and_rejects_bad_shape() {
        let query=r#"{"operation":"inspect_brep",
            "source":{"kind":"block","origin":[1,2,3],
            "width":2,"height":3,"depth":4}}"#;
        let response=process(query);
        assert!(response.ok,"{:?}",response.error);
        let solid=response.solid.unwrap();
        assert_eq!((solid.topology.vertices,solid.topology.faces),(8,6));
        let brep=solid.brep.unwrap();
        assert_eq!(brep["representation"],"polyhedral");
        assert_eq!(brep["geometry"]["cylindrical_faces"],0);
        assert_eq!(brep["geometry"]["seam_edges"],json!([]));
        assert!(!process(r#"{"operation":"inspect_brep",
            "source":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":0,"height":5}}"#).ok);
    }

    #[test]
    fn query_operation_projects_cylinder_without_display_mesh() {
        let response=process(r#"{
            "operation":"geometry_query",
            "source":{"kind":"analytic_cylinder","origin":[0,0,0],"radius":2,"height":5},
            "query":"{ body { topology { faceCount seamEdges } faces(kind:\"cylinder\") {id surface {kind radius}} } }"
        }"#);
        assert!(response.ok,"{:?}",response.error);
        assert!(response.solid.is_none());
        let value=response.data.unwrap();
        assert_eq!(value["body"]["topology"]["seamEdges"],1);
        assert_eq!(value["body"]["faces"][0]["surface"]["kind"],"cylinder");
        assert_eq!(value["body"]["faces"][0]["surface"]["radius"],2.0);
    }

    #[test]
    fn query_rejects_unknown_field_and_keeps_legacy_outputs() {
        let request=r#"{"operation":"geometry_query",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "query":"{body {geometryMutate}}"}"#;
        let response=process(request);
        assert!(!response.ok);
        assert!(response.error.unwrap().contains("query.unknown_field"));
        let legacy=process(r#"{"operation":"block","origin":[0,0,0],
            "width":2,"height":3,"depth":4}"#);
        assert!(legacy.ok);
        assert!(legacy.data.is_none());
        assert!(legacy.solid.is_some());
    }

    #[test]
    fn probe_operation_reports_reproducible_boundary_accuracy() {
        let response=process(r#"{"operation":"probe_kernel","seed":13}"#);
        assert!(response.ok,"error={:?}; failing cases={:#?}", response.error,
            response.data.as_ref().and_then(|d|d["outcomes"].as_array())
                .map(|cases|cases.iter().filter(|c|c["passed"]==false).collect::<Vec<_>>()));
        assert!(response.solid.is_none());
        let data=response.data.unwrap();
        assert_eq!(data["seed"],13);
        assert_eq!(data["failed"],0);
        let outcomes=data["outcomes"].as_array().unwrap();
        assert!(outcomes.len()>=24);
        assert_eq!(outcomes[0]["case"]["input"]["kind"],"block");
        assert!(outcomes.iter().any(|v|v["case"]["id"]=="cylinder/eight_facets"));
    }

    #[test]
    fn pinned_first_party_regressions_are_queryable_without_external_kernels() {
        let result=process(r#"{"operation":"regression_corpus"}"#);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let report=result.data.expect("embedded regression corpus should be returned");
        assert_eq!(report["schema"],"rustsolid.regressions.v1");
        assert_eq!(report["total"],20);
        assert_eq!(report["failed"],0);
        assert_eq!(report["outcomes"].as_array().unwrap().len(),20);
        assert!(report["outcomes"].as_array().unwrap().iter()
            .any(|outcome|outcome["fixture"]["id"]=="regression/thin_face_frame_2026_10_10"));
    }

    #[test]
    fn json_intersection_returns_unfaceted_curve_and_independent_residuals(){
        let result=process(r#"{
            "operation":"intersect_surfaces",
            "first":{"kind":"plane","origin":[0,3,0],"x_axis":[1,0,0],"y_axis":[0,0,-1]},
            "second":{"kind":"cylinder","origin":[0,0,0],"x_axis":[1,0,0],"y_axis":[0,0,-1],"radius":2},
            "samples":9
        }"#);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let body=result.data.unwrap();
        assert_eq!(body["kind"],"crossing");
        assert_eq!(body["contact"],"crossing");
        assert_eq!(body["curves"][0]["periodic"],true);
        assert_eq!(body["curves"][0]["samples"].as_array().unwrap().len(),9);
        for point in body["curves"][0]["samples"].as_array().unwrap(){
            let p=point["point"].as_array().unwrap();
            let x=p[0].as_f64().unwrap();
            let y=p[1].as_f64().unwrap();
            let z=p[2].as_f64().unwrap();
            assert!((y-3.0).abs()<1e-9);
            assert!((x*x+z*z-4.0).abs()<1e-9);
        }
        assert!(body["certificate"]["max_residual_first"].as_f64().unwrap()<1e-9);
        assert!(body["certificate"]["max_residual_second"].as_f64().unwrap()<1e-9);
    }

    #[test]
    fn json_intersection_handles_tangency_uncertainty_and_limits(){
        let source=|x:f64|format!(r#"{{
            "operation":"intersect_surfaces",
            "first":{{"kind":"plane","origin":[{x},0,0],
                "x_axis":[0,1,0],"y_axis":[0,0,1]}},
            "second":{{"kind":"cylinder","origin":[0,0,0],
                "x_axis":[1,0,0],"y_axis":[0,0,-1],"radius":2}}
        }}"#);
        let tangent=process(&source(2.0));
        assert!(tangent.ok,"{:?}",tangent.error);
        assert_eq!(tangent.data.unwrap()["contact"],"tangent");
        let near=process(&source(2.0+4e-10));
        assert!(near.ok,"{:?}",near.error);
        let result=near.data.unwrap();
        assert_eq!(result["kind"],"ambiguous");
        assert_eq!(result["reason"],"near_tangency");
        assert_eq!(result["curves"],json!([]));
        let invalid=process(&source(2.0).replace("\"operation\":\"intersect_surfaces\"",
            "\"operation\":\"intersect_surfaces\",\"samples\":1000"));
        assert!(!invalid.ok);
        assert!(invalid.data.is_none());
    }

    #[test]
    fn json_bounded_circle_keeps_analytic_curve_and_residuals(){
        let response=process(r#"{
            "operation":"intersect_faces",
            "first":{"source":{"kind":"block","origin":[-3,2,-3],
                "width":6,"height":1,"depth":6},"face_id":0},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2},
            "samples":9
        }"#);
        assert!(response.ok,"{:?}",response.error);
        let data=response.data.unwrap();
        assert_eq!(data["kind"],"segments");
        assert_eq!(data["branches"].as_array().unwrap().len(),1);
        assert_eq!(data["branches"][0]["curve"],"periodic_section");
        assert_eq!(data["branches"][0]["closed"],true);
        let samples=data["branches"][0]["samples"].as_array().unwrap();
        for point in samples {
            let [x,y,z] = [point["point"][0].as_f64().unwrap(),
                point["point"][1].as_f64().unwrap(),
                point["point"][2].as_f64().unwrap()];
            assert!((y-2.0).abs()<1e-9);
            assert!((x*x+z*z-4.0).abs()<1e-9);
        }
        assert!(data["sampled_support_certificate"]["max_residual_first"].as_f64().unwrap()<1e-9);
        assert!(data["sampled_support_certificate"]["max_residual_second"].as_f64().unwrap()<1e-9);
    }
    #[test]
    fn json_trimmed_intersection_rejects_unsupported_face_and_cost(){
        let base=r#"{
            "operation":"intersect_faces",
            "first":{"source":{"kind":"block","origin":[-3,2,-3],
                "width":6,"height":1,"depth":6},"face_id":999},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2}
        }"#;
        let invalid=process(base);
        assert!(!invalid.ok);
        assert!(invalid.data.is_none());
        let oversized=base.replace("\"face_id\":999","\"face_id\":0").replace(
            "\"operation\":\"intersect_faces\"", "\"operation\":\"intersect_faces\",\"samples\":200");
        assert!(!process(&oversized).ok);
    }

    #[test]
    fn json_intersection_graph_preserves_sources_closed_loop_and_tolerance(){
        let response=process(r#"{
            "operation":"intersection_graph",
            "first":{"source":{"kind":"block","origin":[-3,2,-3],
                "width":6,"height":1,"depth":6},"face_id":0},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2},
            "samples":9
        }"#);
        assert!(response.ok,"{:?}",response.error);
        assert!(response.solid.is_none());
        let data=response.data.unwrap();
        assert_eq!(data["kind"],"ready");
        assert_eq!(data["nodes"].as_array().unwrap().len(),1);
        assert_eq!(data["arcs"].as_array().unwrap().len(),1);
        assert_eq!(data["nodes"][0]["kind"],"periodic_anchor");
        assert_eq!(data["arcs"][0]["closed"],true);
        assert_eq!(data["arcs"][0]["from"],data["arcs"][0]["to"]);
        assert_eq!(data["topology_committed"],false);
        assert_eq!(data["sources"][0]["operand"],"first");
        assert_eq!(data["sources"][1]["operand"],"second");
        assert!(data["local_tolerance"]["length"].as_f64().unwrap()>0.0);
        for position in data["arcs"][0]["samples"].as_array().unwrap(){
            let x=position["point"][0].as_f64().unwrap();
            let y=position["point"][1].as_f64().unwrap();
            let z=position["point"][2].as_f64().unwrap();
            assert!((y-2.0).abs()<1e-8);
            assert!((x*x+z*z-4.0).abs()<1e-8);
        }
    }
    #[test]
    fn json_intersection_graph_errors_and_ambiguity_are_safe(){
        let request=r#"{
            "operation":"intersection_graph",
            "first":{"source":{"kind":"block","origin":[-3,7,-3],
                "width":6,"height":1,"depth":6},"face_id":0},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2}
        }"#;
        let empty=process(request);
        assert!(empty.ok,"{:?}",empty.error);
        assert_eq!(empty.data.unwrap()["arcs"],json!([]));
        let wrong=request.replace("\"face_id\":0","\"face_id\":777");
        let rejected=process(&wrong);
        assert!(!rejected.ok);
        assert!(rejected.data.is_none());
        let oversized=request.replace("\"operation\":\"intersection_graph\"",
            "\"operation\":\"intersection_graph\",\"samples\":1000");
        assert!(!process(&oversized).ok);
    }

    #[test]
    fn json_plan_reports_symbolic_brep_intent_without_mutation(){
        let response=process(r#"{
            "operation":"topology_delta_plan",
            "first":{"source":{"kind":"block","origin":[-3,2,-3],
                "width":6,"height":1,"depth":6},"face_id":0},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2},
            "samples":7
        }"#);
        assert!(response.ok,"{:?}",response.error);
        assert!(response.solid.is_none());
        let d=response.data.unwrap();
        assert_eq!(d["state"],"preview_only");
        assert_eq!(d["can_apply"],false);
        assert_eq!(d["topology_committed"],false);
        assert_eq!(d["intersection_arcs"].as_array().unwrap().len(),1);
        assert_eq!(d["intersection_arcs"][0]["closed"],true);
        assert_eq!(d["node_bindings"][0]["placement"]["kind"],"virtual_periodic_anchor");
        assert_eq!(d["face_changes"].as_array().unwrap().len(),2);
        assert_eq!(d["face_changes"][0]["effect"],"candidate_inner_loop");
        assert!(d["requirements"].as_array().unwrap().contains(&json!("inner_loop_construction")));
        assert_eq!(d["intersection_arcs"][0]["samples"].as_array().unwrap().len(),7);
    }
    #[test]
    fn json_plan_keeps_ambiguous_intersections_on_hold(){
        let input=r#"{
            "operation":"topology_delta_plan",
            "first":{"source":{"kind":"block","origin":[2.0000000004,1,-3],
                "width":2,"height":2,"depth":6},"face_id":5},
            "second":{"source":{"kind":"analytic_cylinder",
                "origin":[0,0,0],"radius":2,"height":5},"face_id":2}
        }"#;
        let result=process(input);
        assert!(result.ok,"{:?}",result.error);
        let data=result.data.unwrap();
        assert_eq!(data["state"],"held");
        assert_eq!(data["graph_state"],"ambiguous");
        assert_eq!(data["can_apply"],false);
        assert_eq!(data["topology_committed"],false);
        assert_eq!(data["intersection_arcs"],json!([]));
        assert!(data["requirements"].as_array().unwrap().contains(&json!("numerical_refinement")));
        let invalid=input.replace("\"face_id\":5","\"face_id\":999");
        assert!(!process(&invalid).ok);
        let too_many=input.replace("\"operation\":\"topology_delta_plan\"",
            "\"operation\":\"topology_delta_plan\",\"samples\":1000");
        assert!(!process(&too_many).ok);
    }


    #[test]
    fn json_boolean_region_plan_samples_curved_analytic_cut_without_commit() {
        let response=process(r#"{"operation":"boolean_region_plan","mode":"difference",
            "first":{"source":{"kind":"block","origin":[-3,2,-3],
                "width":6,"height":1,"depth":6},"face_id":0},
            "second":{"source":{"kind":"analytic_cylinder","origin":[0,0,0],
                "radius":2,"height":5},"face_id":2}}"#);
        assert!(response.ok,"{:?}",response.error);
        assert!(response.solid.is_none());
        let data=response.data.unwrap();
        assert_eq!(data["state"],"classified_preview");
        assert_eq!(data["can_apply"],false);
        assert_eq!(data["topology_committed"],false);
        assert_eq!(data["intersection_arcs"],1);
        assert_eq!(data["region_samples"].as_array().unwrap().len(),6);
        assert_eq!(data["topology_intent"]["can_apply"],false);
        assert!(data["topology_intent"]["remaining_requirements"].as_array().unwrap()
            .contains(&json!("inner_loop_construction")));
        assert!(data["topology_intent"]["remaining_requirements"].as_array().unwrap()
            .contains(&json!("region_classification")),
            "global region completeness remains a prerequisite");
        assert!(data["region_samples"].as_array().unwrap().iter()
            .any(|v|v["operand"]=="second" && v["reverse_orientation"]==true));
    }


    #[test]
    fn json_curved_boolean_true_genus_one_and_analytic_mass(){
        let request=r#"{"operation":"boolean","mode":"difference",
          "left":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
          "right":{"kind":"analytic_cylinder","origin":[0,-1,0],"radius":1,"height":4}}"#;
        let response=process(request);
        assert!(response.ok,"{:?}",response.error);
        let info=response.data.unwrap();
        let solid=response.solid.unwrap();
        assert_eq!(info["representation"],"analytic_cylindrical_through_hole");
        assert_eq!(info["genus"],1);
        assert_eq!(info["topology_committed"],true);
        assert_eq!(solid.analytic_brep.as_ref().unwrap()["inner_loops"],2);
        assert_eq!(solid.topology.euler_characteristic,0);
        assert_eq!(solid.topology.faces,7);
        assert!((solid.mass.volume-(144.0-4.0*std::f64::consts::PI)).abs()<1e-9);
        assert_eq!(info["face_sources"].as_array().unwrap().len(),7);
    }
    #[test]
    fn json_analytic_hole_is_queryable_without_facet_driven_carriers(){
        let request=r#"{"operation":"inspect_brep","source":{
            "kind":"analytic_through_hole","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6,"center":[0,0],"radius":1,
            "segments":32}}"#;
        let response=process(request);
        assert!(response.ok,"{:?}",response.error);
        let s=response.solid.unwrap();
        let graph=s.brep.unwrap();
        assert_eq!(graph["representation"],"analytic_through_hole");
        assert_eq!(graph["topology"]["genus"],1);
        assert_eq!(graph["geometry"]["cylindrical_faces"],1);
        assert_eq!(graph["geometry"]["circular_edges"],2);
        assert_eq!(s.topology.euler_characteristic,0);
        assert!(s.vertices.len()>64);
    }
    #[test]
    fn json_boolean_rejects_blind_or_tangent_cylinders_not_fake_boolean(){
        for right in [
            r#"{"kind":"analytic_cylinder","origin":[0,-1,0],"radius":1,"height":3}"#,
            r#"{"kind":"analytic_cylinder","origin":[2.5,-1,0],"radius":0.5,"height":4}"#,
        ] {
            let request=format!(r#"{{"operation":"boolean","mode":"difference",
                "left":{{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6}},
                "right":{right}}}"#);
            let response=process(&request);
            assert!(!response.ok);
            assert!(response.solid.is_none());
            assert!(response.data.is_none());
        }
    }


    #[test]
    fn json_structural_opt_in_display_mesh_preserves_true_analytic_brep(){
        let input=r#"{"operation":"structural_transaction",
          "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
          "batches":[[
            {"kind":"add_through_hole","center":[-1.5,0],"radius":0.5},
            {"kind":"add_through_hole","center":[1.5,0],"radius":0.7}
          ]],"display_segments":48}"#;
        let result=process(input);
        assert!(result.ok,"{:?}",result.error);
        let data=result.data.unwrap();
        let solid=result.solid.unwrap();
        assert_eq!(data["hole_count"],2);
        assert_eq!(data["topology"]["genus"],2);
        assert_eq!(data["display_mesh_available"],true);
        assert_eq!(data["display_mesh_included"],true);
        assert_eq!(data["display_quality"]["segments_per_hole"],48);
        assert_eq!(data["display_quality"]["closed_oriented_manifold_certified"],true);
        assert!(data["display_quality"]["volume_faceting_error"].as_f64().unwrap()>0.0);
        assert_eq!(data["display_quality"]["bottom_cap_region"]["loop_count"],3);
        assert_eq!(solid.faces.len(),solid.triangle_face_ids.len());
        assert_eq!(solid.vertices.len(),8+2*48*2);
        assert_eq!(solid.topology.faces,8);
        assert_eq!(solid.topology.euler_characteristic,-2);
        assert_eq!(solid.analytic_brep.as_ref().unwrap()["mesh_is_display_only"],true);
        assert!((solid.mass.volume-(144.0-4.0*std::f64::consts::PI*(0.25+0.49))).abs()<1e-9);
    }
    #[test]
    fn json_structural_display_segmentation_error_does_not_return_partial_brep(){
        for segments in [0,11,513]{
            let request=format!(r#"{{"operation":"structural_transaction",
              "source":{{"kind":"block","origin":[-3,-1,-3],
                  "width":6,"height":4,"depth":6}},
              "batches":[[{{"kind":"add_through_hole",
                  "center":[0,0],"radius":1}}]],
              "display_segments":{segments}}}"#);
            let response=process(&request);
            assert!(!response.ok);
            assert!(response.solid.is_none());
            assert!(response.data.is_none());
            assert!(response.error.unwrap().contains("segments"));
        }
    }
    #[test]
    fn json_structural_brep_adds_two_true_holes_and_recomputes_genus_mass(){
        let input=r#"{"operation":"structural_transaction",
          "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
          "batches":[[
            {"kind":"add_through_hole","center":[-1.5,0],"radius":0.5},
            {"kind":"add_through_hole","center":[1.5,0],"radius":0.7},
            {"kind":"reparameterize_hole_edge","slot":1,
             "role":"bottom_circle","scale":1,"offset":0.4},
            {"kind":"rebase_base_face_uv","face_id":0,"du":0.1,"dv":-0.2}
          ]]}"#;
        let result=process(input);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let d=result.data.unwrap();
        assert_eq!(d["hole_count"],2);
        assert_eq!(d["representation"],"analytic_multiple_through_holes");
        assert_eq!(d["topology"]["genus"],2);
        assert_eq!(d["topology"]["euler_characteristic"],-2);
        assert_eq!(d["topology"]["inner_loops"],4);
        assert_eq!(d["topology"]["cylindrical_faces"],2);
        assert_eq!(d["topology_committed"],true);
        assert_eq!(d["display_mesh_available"],true);
        assert_eq!(d["display_mesh_included"],false);
        assert_eq!(d["commits"][0]["command_count"],4);
        assert!((d["mass"]["volume"].as_f64().unwrap()-
            (144.0-4.0*std::f64::consts::PI*(0.25+0.49))).abs()<1e-9);
    }
    #[test]
    fn json_structural_compaction_and_undo_redo_keep_true_topology(){
        let input=r#"{"operation":"structural_transaction",
          "source":{"kind":"block","origin":[-3,-1,-3],"width":6,"height":4,"depth":6},
          "batches":[[
            {"kind":"add_through_hole","center":[-1.5,0],"radius":0.5},
            {"kind":"add_through_hole","center":[1.5,0],"radius":0.7}
          ],[
            {"kind":"remove_through_hole","slot":0}
          ]],"undo":1,"redo":1}"#;
        let response=process(input);
        assert!(response.ok,"{:?}",response.error);
        let d=response.data.unwrap();
        assert_eq!(d["revision"],4);
        assert_eq!(d["hole_count"],1);
        assert_eq!(d["topology"]["genus"],1);
        assert_eq!(d["commits"][1]["remapped"][0]["side_face_before"],7);
        assert_eq!(d["commits"][1]["remapped"][0]["side_face_after"],6);
        assert_eq!(d["transitions"].as_array().unwrap().len(),2);
    }
    #[test]
    fn json_structural_overlap_or_invalid_carrier_cannot_expose_partial_brep(){
        for pair in [
            r#"{"kind":"add_through_hole","center":[-1.4,0],"radius":0.8}"#,
            r#"{"kind":"rebase_hole_face_uv","slot":0,"du":0.125,"dv":0}"#,
        ] {
            let request=format!(r#"{{"operation":"structural_transaction",
                "source":{{"kind":"block","origin":[-3,-1,-3],
                    "width":6,"height":4,"depth":6}},
                "batches":[[
                    {{"kind":"add_through_hole","center":[-1.5,0],"radius":0.7}},
                    {pair}
                ]]}}"#);
            let r=process(&request);
            assert!(!r.ok,"invalid topology or carrier edit should not commit");
            assert!(r.solid.is_none());
            assert!(r.data.is_none());
        }
    }
    #[test]
    fn json_canonical_carrier_batches_reparameterize_without_changing_mass(){
        let result=process(r#"{"operation":"carrier_transaction",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[[
                {"kind":"reparameterize_edge","edge_id":0,"scale":2,"offset":0.25},
                {"kind":"reparameterize_trim","coedge_id":0,"scale":1.5,"offset":0.125},
                {"kind":"rebase_face_uv","face_id":0,"du":0.5,"dv":-0.25}
            ]],"undo":1,"redo":1}"#);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let data=result.data.unwrap();
        assert_eq!(data["revision"],3);
        assert_eq!(data["shape_preserving"],true);
        assert_eq!(data["topology_mutated"],false);
        assert_eq!(data["mass"]["volume"],24.0);
        assert_eq!(data["topology"]["euler_characteristic"],2);
        assert_eq!(data["commits"][0]["changed"].as_array().unwrap()
            .iter().filter(|x|x["owner_kind"]=="face").count(),1);
    }
    #[test]
    fn json_carrier_transaction_preserves_curved_hole_and_refuses_bad_uv(){
        let source=r#"{"kind":"analytic_through_hole",
            "origin":[-3,-1,-3],"width":6,"height":4,"depth":6,
            "center":[0,0],"radius":1}"#;
        let request=format!(r#"{{"operation":"carrier_transaction",
            "source":{source},"batches":[[
                {{"kind":"reparameterize_edge","edge_id":12,"scale":1,"offset":0.5}},
                {{"kind":"rebase_face_uv","face_id":6,"du":0,"dv":0.75}}
            ]]}}"#);
        let result=process(&request);
        assert!(result.ok,"{:?}",result.error);
        let d=result.data.unwrap();
        assert_eq!(d["topology"]["genus"],1);
        assert_eq!(d["topology"]["euler_characteristic"],0);
        assert_eq!(d["shape_preserving"],true);
        let invalid=format!(r#"{{"operation":"carrier_transaction",
            "source":{source},"batches":[[
                {{"kind":"reparameterize_edge","edge_id":12,"scale":1,"offset":0.5}},
                {{"kind":"rebase_face_uv","face_id":6,"du":0.125,"dv":0}}
            ]]}}"#);
        let result=process(&invalid);
        assert!(!result.ok);
        assert!(result.solid.is_none());
        assert!(result.data.is_none());
    }

    #[test]
    fn json_blind_cylinder_difference_executes_real_floor_face_and_exact_mass(){
        let request=r#"{"operation":"boolean","mode":"difference",
            "left":{"kind":"block","origin":[-3,-1,-3],
                "width":6,"height":4,"depth":6},
            "right":{"kind":"analytic_cylinder","origin":[0,1,0],
                "radius":1,"height":2}}"#;
        let r=process(request);
        assert!(r.ok,"{:?}",r.error);
        let data=r.data.unwrap();
        let solid=r.solid.unwrap();
        assert_eq!(data["representation"],"analytic_cylindrical_blind_hole");
        assert_eq!(data["topology_committed"],true);
        assert_eq!(data["genus"],0);
        assert_eq!(data["opening_face_id"],1);
        assert_eq!(data["floor_face_id"],7);
        assert_eq!(data["side_face_id"],6);
        assert_eq!(data["opening_inner_loops"],1);
        assert_eq!(data["source_opening_trace"]["state"],"ready");
        assert_eq!(data["source_opening_trace"]["arcs"],1);
        assert_eq!(data["source_opening_trace"]["closed"],true);
        assert_eq!(data["closed_oriented_display_mesh_certified"],true);
        assert_eq!(data["face_sources"].as_array().unwrap().len(),8);
        assert_eq!(solid.analytic_brep.as_ref().unwrap()["inner_loops"],1);
        assert_eq!(solid.analytic_brep.as_ref().unwrap()["floor_face_id"],7);
        assert_eq!(solid.topology.euler_characteristic,2);
        assert_eq!(solid.topology.faces,8);
        assert_eq!(solid.faces.len(),solid.triangle_face_ids.len());
        assert!(solid.triangle_face_ids.contains(&7));
        assert!(solid.triangle_face_ids.contains(&6));
        assert!((solid.mass.volume-(144.-2.*std::f64::consts::PI)).abs()<1e-9);
        assert!((solid.mass.surface_area-(168.+4.*std::f64::consts::PI)).abs()<1e-9);
        assert!(data["volume_faceting_error"].as_f64().unwrap()>0.0);
    }
    #[test]
    fn json_blind_bore_inspection_preserves_3d_curves_and_uv_loops(){
        let r=process(r#"{"operation":"inspect_brep","source":{
            "kind":"analytic_blind_hole","origin":[-3,-1,-3],
            "width":6,"height":4,"depth":6,
            "center":[0,0],"radius":1,"bore_depth":2,"segments":32}}"#);
        assert!(r.ok,"{:?}",r.error);
        let solid=r.solid.unwrap();
        let model=solid.brep.unwrap();
        assert_eq!(model["representation"],"analytic_blind_hole");
        assert_eq!(model["topology"]["faces"],8);
        assert_eq!(model["topology"]["genus"],0);
        assert_eq!(model["topology"]["euler_characteristic"],2);
        assert_eq!(model["geometry"]["cylindrical_faces"],1);
        assert_eq!(model["geometry"]["circular_edges"],2);
        assert!((solid.mass.volume-(144.-2.*std::f64::consts::PI)).abs()<1e-9);
    }
    #[test]
    fn json_blind_bore_rejects_unmatched_top_cap_and_bottom_tangency(){
        for right in [
            r#"{"kind":"analytic_cylinder","origin":[0,1,0],"radius":1,"height":1.99999}"#,
            r#"{"kind":"analytic_cylinder","origin":[0,-1,0],"radius":1,"height":3}"#,
            r#"{"kind":"analytic_cylinder","origin":[0,1,0],"radius":3,"height":2}"#
        ]{
            let input=format!(r#"{{"operation":"boolean","mode":"difference",
                "left":{{"kind":"block","origin":[-3,-1,-3],
                    "width":6,"height":4,"depth":6}},"right":{right}}}"#);
            let r=process(&input);
            assert!(!r.ok);
            assert!(r.solid.is_none());
            assert!(r.data.is_none());
        }
    }
    #[test]
    fn json_boolean_union_and_intersection_return_exact_result_and_face_sources(){
        for (mode,expected) in [("intersection",6.0),("union",42.0),("difference",18.0)] {
            let request=format!(r#"{{"operation":"boolean","mode":"{}",
                "left":{{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4}},
                "right":{{"kind":"block","origin":[1,1,1],"width":2,"height":3,"depth":4}}}}"#,mode);
            let result=process(&request);
            assert!(result.ok,"{:?}",result.error);
            let solid=result.solid.unwrap();
            assert_eq!(solid.mass.volume,expected);
            assert_eq!(solid.topology.euler_characteristic,2);
            let info=result.data.unwrap();
            assert_eq!(info["state"],"solid");
            assert_eq!(info["representation"],"exact_planar_box_arrangement");
            assert_eq!(info["face_sources"].as_array().unwrap().len(),solid.topology.faces);
        }
    }
    #[test]
    fn json_boolean_empty_and_unsupported_result_have_distinct_contracts(){
        let empty=process(r#"{"operation":"boolean","mode":"intersection",
            "left":{"kind":"block","origin":[0,0,0],"width":1,"height":1,"depth":1},
            "right":{"kind":"block","origin":[5,0,0],"width":1,"height":1,"depth":1}}"#);
        assert!(empty.ok);
        assert!(empty.solid.is_none());
        assert_eq!(empty.data.unwrap()["state"],"empty");
        let unsupported=process(r#"{"operation":"boolean","mode":"difference",
            "left":{"kind":"block","origin":[0,0,0],"width":3,"height":3,"depth":3},
            "right":{"kind":"block","origin":[1,1,1],"width":1,"height":1,"depth":1}}"#);
        assert!(!unsupported.ok);
        assert!(unsupported.solid.is_none());
        assert!(unsupported.data.is_none());
        assert!(unsupported.error.unwrap().contains("topology not yet supported"));
    }
    #[test]
    fn json_boolean_analytic_operand_is_refused_not_faceted(){
        let unsupported=process(r#"{"operation":"boolean","mode":"difference",
            "left":{"kind":"block","origin":[0,0,0],"width":3,"height":3,"depth":3},
            "right":{"kind":"analytic_cylinder","origin":[1,1,1],"radius":1,"height":1}}"#);
        assert!(!unsupported.ok);
        assert!(unsupported.solid.is_none());
        assert!(unsupported.error.unwrap().contains("unsupported Boolean operand"));
    }
    #[test]
    fn json_arena_edits_are_one_atomic_batch_with_names_and_mass(){
        let request=r#"{"operation":"arena_transaction",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[[
                {"kind":"split_edge","edge_id":0,"fraction":0.25},
                {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":2}
            ]]}"#;
        let result=process(request);
        assert!(result.ok,"{:?}",result.error);
        assert!(result.solid.is_none());
        let value=result.data.unwrap();
        assert_eq!(value["kind"],"brep_arena");
        assert_eq!(value["revision"],1);
        assert_eq!(value["mass"]["volume"],24.0);
        assert_eq!(value["topology"]["vertices"],9);
        assert_eq!(value["topology"]["faces"],7);
        assert_eq!(value["commits"][0]["commands"],2);
        assert_eq!(value["commits"][0]["journal_frames"],2);
        assert_eq!(value["undo_depth"],1);
        assert!(value["commits"][0]["created"].as_array().unwrap().len()>=4);
    }
    #[test]
    fn json_arena_supports_intra_batch_created_edge_refs() {
        let request=r#"{"operation":"arena_transaction",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[[
                {"kind":"split_edge","edge_id":0,"fraction":0.25},
                {"kind":"split_edge","edge_ref":{"step":0,"output":"remainder_edge"},"fraction":0.5}
            ]],"undo":1,"redo":1}"#;
        let result=process(request);
        assert!(result.ok,"{:?}",result.error);
        let d=result.data.unwrap();
        assert_eq!(d["revision"],3);
        assert_eq!(d["topology"]["vertices"],10);
        assert_eq!(d["topology"]["faces"],6);
        assert_eq!(d["commits"][0]["commands"],2);
        assert_eq!(d["mass"]["volume"],24.0);
    }
    #[test]
    fn json_arena_rejects_forward_and_wrong_kind_outputs_atomically() {
        for ref_body in [
            r#"{"step":2,"output":"remainder_edge"}"#,
            r#"{"step":0,"output":"inserted_vertex"}"#,
        ] {
            let input=format!(r#"{{"operation":"arena_transaction",
                "source":{{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4}},
                "batches":[[
                    {{"kind":"split_edge","edge_id":0,"fraction":0.25}},
                    {{"kind":"split_edge","edge_ref":{},"fraction":0.5}}
                ]]}}"#,ref_body);
            let result=process(&input);
            assert!(!result.ok,"invalid output dependency must reject batch");
            assert!(result.data.is_none(),"failure cannot disclose partial geometry");
        }
    }
    #[test]
    fn json_arena_failure_does_not_expose_partial_output(){
        let result=process(r#"{"operation":"arena_transaction",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[[
                {"kind":"split_edge","edge_id":0,"fraction":0.5},
                {"kind":"split_face","face_id":0,"start_vertex_id":0,"end_vertex_id":1}
            ]]}"#);
        assert!(!result.ok);
        assert!(result.data.is_none());
        assert!(result.solid.is_none());
    }
    #[test]
    fn json_arena_undo_redo_and_curved_rejection(){
        let result=process(r#"{"operation":"arena_transaction",
            "source":{"kind":"block","origin":[0,0,0],"width":2,"height":3,"depth":4},
            "batches":[[{"kind":"split_edge","edge_id":0,"fraction":0.5}]],
            "undo":1,"redo":1}"#);
        assert!(result.ok,"{:?}",result.error);
        let value=result.data.unwrap();
        assert_eq!(value["revision"],3);
        assert_eq!(value["topology"]["vertices"],9);
        assert_eq!(value["mass"]["surface_area"],52.0);
        assert_eq!(value["transitions"].as_array().unwrap().len(),2);
        assert_eq!(value["commits"][0]["created"][0]["handle"]["arena_token"],value["arena_token"]);
        let read_only=process(r#"{"operation":"arena_transaction",
            "source":{"kind":"analytic_cylinder","origin":[0,0,0],"radius":2,"height":5},
            "batches":[]}"#);
        assert!(read_only.ok);
        let info=read_only.data.unwrap();
        assert_eq!(info["is_editable"],false);
        assert_eq!(info["topology"]["faces"],3);
        let denied=process(r#"{"operation":"arena_transaction",
            "source":{"kind":"analytic_cylinder","origin":[0,0,0],"radius":2,"height":5},
            "batches":[[{"kind":"split_edge","edge_id":0,"fraction":0.5}]]}"#);
        assert!(!denied.ok);
    }

}
