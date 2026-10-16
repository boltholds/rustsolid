//! G4b2a: complete, closed canonical B-rep candidate for two certified AABBs.
//!
//! Restricted adapter over the EXISTING exact planar box BooleanKernel.
//! Independently proves boundary coverage of all twelve original source faces,
//! vertex identity on the exact source-coordinate plane lattice, and closed
//! coedge / vertex-fan / Euler-genus incidence on the generated BrepModel.
//! This is NOT the general G1-G4a UV/cylinder Boolean assembler, an independent
//! volume oracle for curved solids, or permission for G5 atomic mutation.

use crate::{
    BooleanError, BooleanFaceSource, BooleanKernel, BooleanOperand,
    BooleanOperation, BooleanOutcome, BrepBody, BrepModel, FaceId,
    GeometryError, GeometryTolerance, GraphOperand, LoopRole, ModelUnits,
    Point3, Solid, VertexId, RadialTopologyCertificate, RadialTopologyHold,
};
use crate::radial_shell_audit::{certify_radial_shell,RadialHalfedge,RadialLoop};
use std::collections::{BTreeMap,BTreeSet};

const MAX_BOX_RESULT_FACES: usize = 256;
const MAX_BOX_RESULT_VERTICES: usize = 256;

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn coordinate(p:Point3,axis:usize)->f64{p.array()[axis]}
fn calc_volume(bounds:&AxisBox)->f64{
    (bounds.max[0]-bounds.min[0])*
    (bounds.max[1]-bounds.min[1])*
    (bounds.max[2]-bounds.min[2])
}
#[derive(Debug,Clone,Copy)]
struct AxisBox{min:[f64;3],max:[f64;3]}
impl AxisBox {
    fn from_solid(s:&Solid)->Self{
        Self{min:s.bbox.min.array(),max:s.bbox.max.array()}
    }
    fn occupies(self,grid:&[Vec<f64>;3],cell:[usize;3])->bool{
        (0..3).all(|axis|
            cell[axis]+1<grid[axis].len() &&
            grid[axis][cell[axis]]>=self.min[axis] &&
            grid[axis][cell[axis]+1]<=self.max[axis])
    }
}
fn material(a:bool,b:bool,op:BooleanOperation)->bool{
    match op{
        BooleanOperation::Union=>a||b,
        BooleanOperation::Intersection=>a&&b,
        BooleanOperation::Difference=>a&&!b,
    }
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum SourceFaceDisposition {
    /// Entire original face remains on the result boundary.
    RetainedUncut,
    /// Original face contributes no result-boundary area.
    Discarded,
    /// Exact source-coordinate tiles are partly kept and partly removed.
    Partitioned,
}
#[derive(Debug,Clone,PartialEq)]
pub struct SourceFaceCoverage {
    pub operand: GraphOperand,
    pub face: FaceId,
    pub plane_axis: u8,
    pub source_outward_positive: bool,
    pub total_tiles: usize,
    pub exposed_tiles: usize,
    pub original_area: f64,
    pub exposed_area: f64,
    pub disposition: SourceFaceDisposition,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct CandidateVertexLattice {
    /// This is a result-local canonical VertexId, never a source VertexId.
    pub result_vertex: VertexId,
    /// Exact indices of the 4 distinct input boundary planes per world axis.
    /// No proximity-based geometric snapping is involved.
    pub coordinate_index: [u8;3],
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum CanonicalBoxCandidateHold{
    UnsupportedOperand,
    NearCoincidentSourcePlanes,
    EmptyResult,
    UnsupportedResultTopology,
    MissingSourceFaceCoverage,
    IncorrectSourceFaceOrientation,
    InvalidExactVertexLattice,
    InvalidRadialTopology(RadialTopologyHold),
    InvalidAnalyticBoxMass,
    BudgetExceeded,
}
#[derive(Debug,Clone)]
pub struct CanonicalBoxCandidate {
    /// Actual owned B-rep, obtained from the independently validated
    /// planar box CSG implementation and then re-audited here.
    model:BrepModel,
    operation:BooleanOperation,
    face_sources:Vec<BooleanFaceSource>,
    source_faces:Vec<SourceFaceCoverage>,
    vertices:Vec<CandidateVertexLattice>,
    radial:RadialTopologyCertificate,
    /// Shell topology is proven for this strict planar AABB-only case.
    closed_shell_proven:bool,
    /// G5 transaction and general Boolean support are not enabled by this.
    can_apply:bool,
}
#[derive(Debug,Clone)]
pub enum CanonicalBoxCandidateResult {
    Ready(CanonicalBoxCandidate),
    Held(CanonicalBoxCandidateHold),
}
impl CanonicalBoxCandidate{
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn operation(&self)->BooleanOperation{self.operation}
    pub fn source_faces(&self)->&[SourceFaceCoverage]{&self.source_faces}
    pub fn face_provenance(&self)->&[BooleanFaceSource]{&self.face_sources}
    pub fn vertices(&self)->&[CandidateVertexLattice]{&self.vertices}
    pub fn radial_topology(&self)->&RadialTopologyCertificate{&self.radial}
    pub fn closed_shell_proven(&self)->bool{self.closed_shell_proven}
    pub fn can_apply(&self)->bool{self.can_apply}

    /// Rebuilds every original source face, exact-lattice vertex and radial
    /// incidence. Never merely accepts a saved, potentially stale preview.
    pub fn validate_against(
        &self,left:&BrepBody,right:&BrepBody,
        tolerance:GeometryTolerance,
    )->Result<(),GeometryError>{
        let CanonicalBoxCandidateResult::Ready(rebuilt)=
            Self::build(left,right,self.operation,tolerance)? else{
            return Err(invalid("G4b2 original solids no longer yield candidate"));
        };
        if self.source_faces!=rebuilt.source_faces ||
           self.face_sources!=rebuilt.face_sources ||
           self.vertices!=rebuilt.vertices ||
           self.radial!=rebuilt.radial ||
           self.closed_shell_proven!=rebuilt.closed_shell_proven ||
           self.can_apply!=rebuilt.can_apply ||
           format!("{:?}",self.model)!=format!("{:?}",rebuilt.model)
        {
            return Err(invalid("G4b2 canonical B-rep differs from source reconstruction"));
        }
        self.model.validate()
    }

    pub fn build(
        left:&BrepBody,right:&BrepBody,
        operation:BooleanOperation,tolerance:GeometryTolerance,
    )->Result<CanonicalBoxCandidateResult,GeometryError>{
        tolerance.validate()?;
        left.validate()?;right.validate()?;
        let (BrepBody::Polyhedral(a),BrepBody::Polyhedral(b))=(left,right) else{
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::UnsupportedOperand));
        };
        for s in [a,b] {
            match crate::boolean::certify_box(s) {
                Ok(())=>{},
                Err(BooleanError::InvalidGeometry(e))=>return Err(e),
                Err(_)=>return Ok(CanonicalBoxCandidateResult::Held(
                    CanonicalBoxCandidateHold::UnsupportedOperand)),
            }
        }
        let box_a=AxisBox::from_solid(a);
        let box_b=AxisBox::from_solid(b);
        let scale=(0..3).map(|i|
            box_a.max[i].max(box_b.max[i])-box_a.min[i].min(box_b.min[i])
        ).fold(1.0_f64,f64::max);
        if !scale.is_finite(){return Err(invalid("G4b2 source scale is invalid"));}
        let eps=tolerance.length_at(scale);
        // STRICT source plane separation: no ambiguous planar coincidences
        // or tolerance-snapped entities in this initial G4b2 slice.
        let grid:[Vec<f64>;3]=std::array::from_fn(|axis|{
            let mut planes=vec![
                box_a.min[axis],box_a.max[axis],
                box_b.min[axis],box_b.max[axis],
            ];
            planes.sort_by(f64::total_cmp);
            planes
        });
        if grid.iter().any(|planes|
            planes.windows(2).any(|pair|
                !pair[0].is_finite()||!pair[1].is_finite()||
                pair[1]-pair[0]<=eps
            ))
        {
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::NearCoincidentSourcePlanes));
        }
        let generated=match BooleanKernel::apply_solids(a,b,operation){
            Ok(result)=>result,
            Err(BooleanError::InvalidGeometry(e))=>return Err(e),
            Err(_)=>return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::UnsupportedResultTopology)),
        };
        let BooleanOutcome::Solid(solid)=generated.outcome else{
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::EmptyResult));
        };
        if solid.faces.len()>MAX_BOX_RESULT_FACES ||
            solid.vertices.len()>MAX_BOX_RESULT_VERTICES ||
            solid.faces.is_empty()
        {
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::BudgetExceeded));
        }
        let model=BrepModel::from_polyhedral(&solid,ModelUnits::default())?;
        let vertex_plan=match certify_vertices(&solid,&grid) {
            Some(proof)=>proof,
            None=>return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::InvalidExactVertexLattice)),
        };
        let coverage=match certify_source_faces(
            a,b,box_a,box_b,&grid,operation,
            &solid,&generated.face_sources
        )?{
            Some(v)=>v,
            None=>return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::MissingSourceFaceCoverage)),
        };
        // Box-only cross-check from original independent dimensions and
        // inclusion/exclusion, NOT a reusable curved-solid volume oracle.
        let overlap=(0..3).map(|axis|
            (box_a.max[axis].min(box_b.max[axis])-
             box_a.min[axis].max(box_b.min[axis])).max(0.0)
        ).product::<f64>();
        let expected=match operation{
            BooleanOperation::Union=>
                calc_volume(&box_a)+calc_volume(&box_b)-overlap,
            BooleanOperation::Intersection=>overlap,
            BooleanOperation::Difference=>calc_volume(&box_a)-overlap,
        };
        if !expected.is_finite()||expected<=0.0||
            (model.mass.volume-expected).abs()>expected.max(1.0)*1e-9
        {
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::InvalidAnalyticBoxMass));
        }
        let radial=match certify_model_radial(&model,eps.max(model.tolerance.absolute_length)*16.0){
            Ok(cert)=>cert,
            Err(hold)=>return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::InvalidRadialTopology(hold))),
        };
        if radial.edges!=model.edges.len() ||
            radial.vertices!=model.vertices.len()||
            radial.material_face_cells!=model.faces.len()||
            radial.shell_components!=1||
            radial.euler_characteristic!=2||
            radial.total_genus!=0
        {
            return Ok(CanonicalBoxCandidateResult::Held(
                CanonicalBoxCandidateHold::InvalidRadialTopology(
                    RadialTopologyHold::InvalidEulerOrGenus)));
        }
        Ok(CanonicalBoxCandidateResult::Ready(Self {
            model,operation,face_sources:generated.face_sources,
            source_faces:coverage,vertices:vertex_plan,radial,
            closed_shell_proven:true,can_apply:false,
        }))
    }
}

/// Exact result-vertex lattice identity: every XYZ coordinate MUST match
/// one of the source faces' original Cartesian planes bit-for-bit.
fn certify_vertices(
    solid:&Solid,grid:&[Vec<f64>;3],
)->Option<Vec<CandidateVertexLattice>>{
    let mut seen=BTreeSet::new();
    let mut result=Vec::with_capacity(solid.vertices.len());
    for (index,vertex) in solid.vertices.iter().enumerate(){
        if vertex.id.0 as usize!=index||!vertex.position.is_finite(){return None;}
        let mut indices=[0u8;3];
        for axis in 0..3 {
            let coord=coordinate(vertex.position,axis);
            let position=grid[axis].iter().position(|plane|*plane==coord)?;
            indices[axis]=u8::try_from(position).ok()?;
        }
        if !seen.insert(indices){return None;}
        result.push(CandidateVertexLattice{
            result_vertex:vertex.id,coordinate_index:indices,
        });
    }
    Some(result)
}

#[derive(Debug,Clone,Copy)]
struct SourceFacePlane{
    axis:usize,side:usize,coordinate:f64,
}
fn original_plane(solid:&Solid,face:FaceId,bounds:AxisBox)
    ->Result<SourceFacePlane,GeometryError>
{
    let orig=solid.faces.get(face.0 as usize)
        .filter(|f|f.id==face)
        .ok_or_else(||invalid("G4b2 source FaceId missing"))?;
    if orig.boundary.len()!=4{return Err(invalid("G4b2 source not a quad"));}
    let mut found=None;
    for axis in 0..3 {
        for side in 0..2{
            let plane=if side==0{bounds.min[axis]}else{bounds.max[axis]};
            if orig.boundary.iter().all(|id|
                solid.vertices.get(id.0 as usize)
                    .is_some_and(|v|coordinate(v.position,axis)==plane))
            {
                if found.is_some(){return Err(invalid("G4b2 ambiguous source face axis"));}
                found=Some(SourceFacePlane{axis,side,coordinate:plane});
            }
        }
    }
    found.ok_or_else(||invalid("G4b2 source face is not an AABB plane"))
}
#[derive(Debug,Clone,Copy)]
struct EmittedTile{
    result_face:FaceId,
    source: BooleanFaceSource,
    outward_positive:bool,
    area:f64,
}
type TileKey=(usize,u32,usize,usize);

/// Independently enumerates every complete original AABB face tile.
/// Halfspace occupancy is constant over an open grid rectangle because
/// both input boxes' exact supporting plane coordinates partition it.
fn certify_source_faces(
    a:&Solid,b:&Solid,aa:AxisBox,bb:AxisBox,
    grid:&[Vec<f64>;3],operation:BooleanOperation,
    generated:&Solid,provenance:&[BooleanFaceSource],
)->Result<Option<Vec<SourceFaceCoverage>>,GeometryError>{
    if generated.faces.len()!=provenance.len(){return Ok(None);}
    let originals=[(a,aa),(b,bb)];
    let mut emitted=BTreeMap::<TileKey,EmittedTile>::new();
    for entry in provenance {
        let (operand,source_bounds)=match entry.operand {
            BooleanOperand::Left=>(0usize,aa),
            BooleanOperand::Right=>(1usize,bb),
        };
        let source_solid=originals[operand].0;
        let plane=original_plane(source_solid,entry.source_face,source_bounds)?;
        let u=(plane.axis+1)%3;let v=(plane.axis+2)%3;
        let Some(face)=generated.faces.get(entry.result_face.0 as usize)
            .filter(|f|f.id==entry.result_face&&f.boundary.len()==4) else{
            return Ok(None);
        };
        let mut min=[f64::INFINITY;3];let mut max=[f64::NEG_INFINITY;3];
        let mut corners=BTreeSet::new();
        let mut points=Vec::new();
        for &id in &face.boundary{
            let Some(vertex)=generated.vertices.get(id.0 as usize)
                .filter(|vertex|vertex.id==id) else{return Ok(None);};
            let xyz=vertex.position.array();
            if xyz[plane.axis]!=plane.coordinate{return Ok(None);}
            for axis in 0..3{
                min[axis]=min[axis].min(xyz[axis]);
                max[axis]=max[axis].max(xyz[axis]);
            }
            points.push(vertex.position);
        }
        for point in &points{
            let xyz=point.array();
            corners.insert((xyz[u].to_bits(),xyz[v].to_bits()));
        }
        if corners.len()!=4{return Ok(None);}
        let iu=grid[u].windows(2).position(|x|x[0]==min[u]&&x[1]==max[u]);
        let iv=grid[v].windows(2).position(|x|x[0]==min[v]&&x[1]==max[v]);
        let (Some(iu),Some(iv))=(iu,iv) else{return Ok(None);};
        let du=points[1].sub(points[0]);
        let dv=points[2].sub(points[1]);
        let normal=du.cross(dv).array();
        let component=normal[plane.axis];
        if !component.is_finite()||component==0.0||
            (0..3).any(|axis|axis!=plane.axis&&normal[axis]!=0.0)
        {
            return Ok(None);
        }
        let outward_positive=component>0.0;
        let area=(max[u]-min[u])*(max[v]-min[v]);
        let key=(operand,entry.source_face.0,iu,iv);
        if !area.is_finite()||area<=0.0||
           emitted.insert(key,EmittedTile{
                result_face:entry.result_face,source:*entry,
                outward_positive,area,
           }).is_some()
        {
            return Ok(None);
        }
    }
    let mut visited=BTreeSet::new();
    let mut coverage=Vec::with_capacity(a.faces.len()+b.faces.len());
    let mut overall_area=0.0;
    for (operand,(source,bounds)) in originals.iter().enumerate(){
        for face in &source.faces{
            let plane=original_plane(source,face.id,*bounds)?;
            let coord=grid[plane.axis].iter().position(|x|*x==plane.coordinate)
                .ok_or_else(||invalid("G4b2 source plane not found in lattice"))?;
            let u=(plane.axis+1)%3;let v=(plane.axis+2)%3;
            let mut total=0usize;let mut kept=0usize;let mut area=0.0;
            let original_area=(bounds.max[u]-bounds.min[u])*
                (bounds.max[v]-bounds.min[v]);
            for iu in 0..grid[u].len()-1 {
                if grid[u][iu]<bounds.min[u]||
                   grid[u][iu+1]>bounds.max[u]{continue;}
                for iv in 0..grid[v].len()-1{
                    if grid[v][iv]<bounds.min[v]||
                       grid[v][iv+1]>bounds.max[v]{continue;}
                    total+=1;
                    let tile=(grid[u][iu+1]-grid[u][iu])*
                        (grid[v][iv+1]-grid[v][iv]);
                    let occupancy=|side:usize| ->bool {
                        if side==0&&coord==0{return false;}
                        if side==1&&coord+1==grid[plane.axis].len(){return false;}
                        let mut cell=[0usize;3];
                        cell[u]=iu;cell[v]=iv;
                        cell[plane.axis]=if side==0{coord-1}else{coord};
                        material(aa.occupies(grid,cell),bb.occupies(grid,cell),operation)
                    };
                    let minus=occupancy(0);
                    let plus=occupancy(1);
                    let boundary=minus!=plus;
                    let key=(operand,face.id.0,iu,iv);
                    let generated_tile=emitted.get(&key);
                    if boundary{
                        let Some(actual)=generated_tile else{return Ok(None);};
                        // Physical outward direction is toward EMPTY material.
                        let positive=minus&&!plus;
                        let source_positive=plane.side==1;
                        if actual.outward_positive!=positive||
                           actual.source.reversed!=(positive!=source_positive)||
                           (actual.area-tile).abs()>tile.max(1.0)*1e-12
                        {
                            return Ok(None);
                        }
                        kept+=1;area+=tile;overall_area+=tile;
                        visited.insert(key);
                    }else if generated_tile.is_some(){
                        return Ok(None);
                    }
                }
            }
            if !original_area.is_finite()||original_area<=0.0||total==0{
                return Ok(None);
            }
            let kind=if kept==0{SourceFaceDisposition::Discarded}
                else if (area-original_area).abs()<=original_area*1e-12{
                    SourceFaceDisposition::RetainedUncut
                }else{SourceFaceDisposition::Partitioned};
            coverage.push(SourceFaceCoverage{
                operand:if operand==0{GraphOperand::First}else{GraphOperand::Second},
                face:face.id,plane_axis:plane.axis as u8,
                source_outward_positive:plane.side==1,
                total_tiles:total,exposed_tiles:kept,
                original_area,exposed_area:area,disposition:kind,
            });
        }
    }
    if visited.len()!=emitted.len() ||
        (overall_area-generated.mass.surface_area).abs()>
        overall_area.max(1.0)*1e-9
    {
        return Ok(None);
    }
    Ok(Some(coverage))
}

/// A second, topology-only closure audit of the **actual** generated
/// B-rep coedges (not the pairwise G4a subdivision halfedges).
fn certify_model_radial(
    model:&BrepModel,eps:f64,
)->Result<RadialTopologyCertificate,RadialTopologyHold>{
    let mut half=Vec::<RadialHalfedge>::new();
    for coedge in &model.coedges {
        let Some(edge)=model.edges.get(coedge.edge.0 as usize) else{
            return Err(RadialTopologyHold::UnpairedEdge);
        };
        let Some(start)=model.vertices.get(coedge.start_vertex(edge).0 as usize) else{
            return Err(RadialTopologyHold::NonManifoldVertex);
        };
        let Some(end)=model.vertices.get(coedge.end_vertex(edge).0 as usize) else{
            return Err(RadialTopologyHold::NonManifoldVertex);
        };
        half.push(RadialHalfedge{
            start:start.position,end:end.position,
            face:coedge.face.0 as usize,
        });
    }
    let mut loops=Vec::<RadialLoop>::new();
    for face in &model.faces {
        for &loop_id in &face.loops {
            let Some(lp)=model.loops.get(loop_id.0 as usize)
                .filter(|lp|lp.id==loop_id&&lp.face==face.id) else{
                return Err(RadialTopologyHold::MissingLoop);
            };
            let mut ids=Vec::new();
            let mut cursor=lp.first_coedge;
            loop {
                let Some(co)=model.coedges.get(cursor.0 as usize)
                    .filter(|co|co.id==cursor&&co.loop_id==loop_id) else{
                    return Err(RadialTopologyHold::MissingLoop);
                };
                if ids.contains(&(cursor.0 as usize)) {
                    return Err(RadialTopologyHold::MissingLoop);
                }
                ids.push(cursor.0 as usize);
                cursor=co.next;
                if cursor==lp.first_coedge {break;}
                if ids.len()>model.coedges.len(){
                    return Err(RadialTopologyHold::OpenLoop);
                }
            }
            loops.push(RadialLoop{
                face:face.id.0 as usize,
                uses:ids,inner:lp.role==LoopRole::Inner,
            });
        }
    }
    let pairs=model.edges.iter().map(|edge|
        [edge.coedges[0].0 as usize,edge.coedges[1].0 as usize]
    ).collect::<Vec<_>>();
    certify_radial_shell(&half,&loops,&pairs,model.faces.len(),eps)
}
