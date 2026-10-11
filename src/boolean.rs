//! Certified first Boolean execution slice: axis-aligned box CSG.
//!
//! The result is rebuilt as an independent exact planar B-rep with a new
//! identity incarnation. No tessellated curved carriers, tolerance snapping,
//! arbitrary polyhedra, multiple shells, cavities, or Boolean commit of a
//! preview-only TopologyDeltaPlan are claimed. Unsupported cases fail closed.
use crate::{
    brep, BoundingBox, BrepBody, Face, FaceId, FaceRole, GeometryError,
    GeometryTolerance, MassProperties, Mesh, Point3, Solid, TopologyIdentity,
    Vertex, VertexId,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BooleanOperation { Union, Intersection, Difference }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BooleanOperand { Left, Right }

/// Every generated quad is traced to a planar face on one source operand.
/// The right operand contributes inward-facing source faces to A - B.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BooleanFaceSource {
    pub result_face: FaceId,
    pub operand: BooleanOperand,
    pub source_face: FaceId,
    pub reversed: bool,
}

#[derive(Debug, Clone)]
pub enum BooleanOutcome { Empty, Solid(Solid), Analytic(crate::ThroughHoleBrep) }

#[derive(Debug, Clone)]
pub struct BooleanResult {
    pub outcome: BooleanOutcome,
    pub face_sources: Vec<BooleanFaceSource>,
}
impl BooleanResult {
    pub fn solid(&self) -> Option<&Solid> {
        match &self.outcome { BooleanOutcome::Empty | BooleanOutcome::Analytic(_) => None, BooleanOutcome::Solid(s) => Some(s) }
    }
    pub fn into_solid(self) -> Option<Solid> {
        match self.outcome { BooleanOutcome::Empty | BooleanOutcome::Analytic(_) => None, BooleanOutcome::Solid(s) => Some(s) }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BooleanError {
    UnsupportedOperand(&'static str),
    AmbiguousBoundary,
    UnsupportedResultTopology(String),
    InvalidGeometry(GeometryError),
}
impl fmt::Display for BooleanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedOperand(s) => write!(f, "unsupported Boolean operand: {s}"),
            Self::AmbiguousBoundary => write!(f,
                "Boolean boundaries are distinct but closer than positional tolerance"),
            Self::UnsupportedResultTopology(s) => write!(f,
                "Boolean result needs topology not yet supported by Solid: {s}"),
            Self::InvalidGeometry(e) => write!(f, "invalid Boolean geometry: {e}"),
        }
    }
}
impl std::error::Error for BooleanError {}

pub(crate) fn certify_box(s:&Solid)->Result<(),BooleanError>{
    AxisBox::certify(s,"through-hole needs a canonical axis-aligned box")?;
    Ok(())
}
fn geometry(error: GeometryError) -> BooleanError {
    match error {
        GeometryError::InvalidTopology(detail) =>
            BooleanError::UnsupportedResultTopology(detail),
        other => BooleanError::InvalidGeometry(other),
    }
}
fn coordinate(p: Point3, axis: usize) -> f64 {
    match axis { 0 => p.x, 1 => p.y, _ => p.z }
}
fn axis_bounds(box_: &AxisBox, axis: usize) -> (f64, f64) {
    (box_.min[axis], box_.max[axis])
}
#[derive(Debug, Clone, Copy)]
struct AxisBox { min: [f64; 3], max: [f64; 3] }
impl AxisBox {
    /// Certify an exact six-quad, eight-corner AABB before treating it as an
    /// occupied halfspace arrangement. A generic polyhedron's bbox is NEVER
    /// a substitute for its actual shape.
    fn certify(s: &Solid, which: &'static str) -> Result<Self, BooleanError> {
        s.validate().map_err(BooleanError::InvalidGeometry)?;
        if s.vertices.len() != 8 || s.faces.len() != 6 ||
            s.faces.iter().any(|f| f.boundary.len() != 4) {
            return Err(BooleanError::UnsupportedOperand(which));
        }
        let bbox = Self { min: s.bbox.min.array(), max: s.bbox.max.array() };
        let mut corners = BTreeSet::new();
        for v in &s.vertices {
            let mut indices = [0u8;3];
            for axis in 0..3 {
                let coord = coordinate(v.position, axis);
                indices[axis] = if coord == bbox.min[axis] { 0 }
                    else if coord == bbox.max[axis] { 1 }
                    else { return Err(BooleanError::UnsupportedOperand(which)); };
            }
            if !corners.insert(indices) { return Err(BooleanError::UnsupportedOperand(which)); }
        }
        if corners.len() != 8 { return Err(BooleanError::UnsupportedOperand(which)); }
        let mut planes = BTreeSet::new();
        for face in &s.faces {
            let mut plane = None;
            for axis in 0..3 {
                for side in 0..2 {
                    let value = if side == 0 { bbox.min[axis] } else { bbox.max[axis] };
                    if face.boundary.iter().all(|&id| coordinate(
                        s.vertices[id.0 as usize].position, axis) == value) {
                        if plane.replace((axis, side)).is_some() {
                            return Err(BooleanError::UnsupportedOperand(which));
                        }
                    }
                }
            }
            if !planes.insert(plane.ok_or(BooleanError::UnsupportedOperand(which))?) {
                return Err(BooleanError::UnsupportedOperand(which));
            }
        }
        if planes.len()!=6 {return Err(BooleanError::UnsupportedOperand(which));}
        let dims = [bbox.max[0]-bbox.min[0],bbox.max[1]-bbox.min[1],
            bbox.max[2]-bbox.min[2]];
        if dims.iter().any(|d|!d.is_finite()||*d<=0.0) {
            return Err(BooleanError::UnsupportedOperand(which));
        }
        let exact_volume = dims[0] * dims[1] * dims[2];
        if !exact_volume.is_finite() ||
            (exact_volume - s.mass.volume).abs() > exact_volume * 1.0e-9 {
            return Err(BooleanError::UnsupportedOperand(which));
        }
        Ok(bbox)
    }
    fn contains_cell(self, axes: &[Vec<f64>;3], i:usize,j:usize,k:usize) -> bool {
        let idx = [i,j,k];
        (0..3).all(|axis| {
            axes[axis][idx[axis]] >= self.min[axis] &&
            axes[axis][idx[axis]+1] <= self.max[axis]
        })
    }
}

/// Closed-form Boolean classification on the plane arrangement of two
/// certified AABBs. Boundaries are exact source coordinates, not a voxel
/// approximation: every grid cell is a *mathematical* partition region.
pub struct BooleanKernel;
impl BooleanKernel {
    pub fn apply(left: &BrepBody, right: &BrepBody,
                 operation: BooleanOperation) -> Result<BooleanResult, BooleanError> {
        match (left,right,operation) {
            (BrepBody::Polyhedral(a),BrepBody::Cylinder(c),BooleanOperation::Difference) => {
                AxisBox::certify(a,"left must be an axis-aligned box")?;
                let analytic=crate::ThroughHoleBrep::from_operands(a,c)
                    .map_err(|_|BooleanError::UnsupportedOperand(
                        "through-hole cylinder must span the full box height and remain strictly inside its sides"))?;
                let face_sources=a.faces.iter().map(|f|BooleanFaceSource {
                    result_face:f.id,operand:BooleanOperand::Left,source_face:f.id,reversed:false,
                }).chain(std::iter::once(BooleanFaceSource {
                    result_face:FaceId(6),operand:BooleanOperand::Right,
                    source_face:FaceId(2),reversed:true,
                })).collect();
                Ok(BooleanResult {outcome:BooleanOutcome::Analytic(analytic),face_sources})
            }
            (BrepBody::Polyhedral(a),BrepBody::Polyhedral(b),operation) =>
                Self::apply_solids(a,b,operation),
            _=>Err(BooleanError::UnsupportedOperand(
                "only planar box CSG or box minus a strictly through cylinder is supported")),
        }
    }

    pub fn apply_solids(a:&Solid,b:&Solid,operation:BooleanOperation)
        -> Result<BooleanResult,BooleanError>
    {
        let left = AxisBox::certify(a,"left must be an exact axis-aligned box")?;
        let right = AxisBox::certify(b,"right must be an exact axis-aligned box")?;
        let tolerance = GeometryTolerance {
            absolute_length: a.tolerance.absolute_length.max(b.tolerance.absolute_length),
            relative_length: a.tolerance.relative_length.max(b.tolerance.relative_length),
            angular: a.tolerance.angular.max(b.tolerance.angular),
        };
        tolerance.validate().map_err(BooleanError::InvalidGeometry)?;
        let extent = (0..3).map(|axis| {
            (left.max[axis].max(right.max[axis]) -
             left.min[axis].min(right.min[axis])).abs()
        }).fold(0.0_f64,f64::max);
        if !extent.is_finite() {
            return Err(BooleanError::InvalidGeometry(GeometryError::InvalidEdit(
                "Boolean world coordinate extent overflow".into())));
        }
        let eps = tolerance.length_at(extent);
        let axes: [Vec<f64>;3] = std::array::from_fn(|axis| {
            let (a0,a1) = axis_bounds(&left,axis);
            let (b0,b1) = axis_bounds(&right,axis);
            let mut coords=vec![a0,a1,b0,b1];
            coords.sort_by(f64::total_cmp);
            coords.dedup_by(|x,y|*x==*y);
            coords
        });
        if axes.iter().any(|coords| coords.windows(2).any(|p|
            p[1]-p[0]<=eps || !(p[1]-p[0]).is_finite())) {
            return Err(BooleanError::AmbiguousBoundary);
        }
        let dims = [axes[0].len()-1, axes[1].len()-1, axes[2].len()-1];
        let (nx,ny,nz) = (dims[0],dims[1],dims[2]);
        let index = |i:usize,j:usize,k:usize| (i*ny+j)*nz+k;
        let mut inside = vec![(false,false);nx*ny*nz];
        let mut occupied = vec![false;nx*ny*nz];
        for i in 0..nx { for j in 0..ny {for k in 0..nz {
            let a = left.contains_cell(&axes,i,j,k);
            let b = right.contains_cell(&axes,i,j,k);
            inside[index(i,j,k)] = (a,b);
            occupied[index(i,j,k)] = match operation {
                BooleanOperation::Union => a || b,
                BooleanOperation::Intersection => a && b,
                BooleanOperation::Difference => a && !b,
            };
        }}}
        if !occupied.iter().any(|filled|*filled) {
            return Ok(BooleanResult {outcome:BooleanOutcome::Empty,face_sources:Vec::new()});
        }
        let mut vertices=Vec::<Vertex>::new();
        let mut by_grid=BTreeMap::<(usize,usize,usize),VertexId>::new();
        let mut faces=Vec::<Face>::new();
        let mut mesh=Mesh {vertices:Vec::new(),triangles:Vec::new(),triangle_faces:Vec::new()};
        let mut sources=Vec::new();
        let mut volume=0.0;
        let mut moments=[0.0;3];
        let mut area=0.0;
        let mut low=[f64::INFINITY;3];
        let mut high=[f64::NEG_INFINITY;3];
        let origin=[axes[0][0],axes[1][0],axes[2][0]];
        for i in 0..nx { for j in 0..ny {for k in 0..nz {
            if !occupied[index(i,j,k)] {continue;}
            let cell=[i,j,k];
            let ds=[axes[0][i+1]-axes[0][i],
                axes[1][j+1]-axes[1][j], axes[2][k+1]-axes[2][k]];
            let v=ds[0]*ds[1]*ds[2];
            if !v.is_finite() || v<=0.0 {
                return Err(BooleanError::InvalidGeometry(GeometryError::InvalidEdit(
                    "Boolean cell volume overflow or underflow".into())));
            }
            volume+=v;
            for axis in 0..3 {
                let start=axes[axis][cell[axis]];
                let end=axes[axis][cell[axis]+1];
                low[axis]=low[axis].min(start);
                high[axis]=high[axis].max(end);
                moments[axis]+=v*((start-origin[axis]) + ds[axis]*0.5);
            }
            for side in 0..6 {
                let neighbor=match side {
                    0 if i>0 => Some((i-1,j,k)),
                    1 if i+1<nx => Some((i+1,j,k)),
                    2 if j>0 => Some((i,j-1,k)),
                    3 if j+1<ny => Some((i,j+1,k)),
                    4 if k>0 => Some((i,j,k-1)),
                    5 if k+1<nz => Some((i,j,k+1)),
                    _ => None,
                };
                if neighbor.is_some_and(|(x,y,z)|occupied[index(x,y,z)]) {continue;}
                let (other_a,other_b)=neighbor.map(|(x,y,z)|inside[index(x,y,z)])
                    .unwrap_or((false,false));
                let (inside_a,inside_b)=inside[index(i,j,k)];
                // An exposed Boolean surface must coincide with a real
                // source face; pick the first source whose classification flips.
                let operand=if inside_a!=other_a {
                    BooleanOperand::Left
                } else if inside_b!=other_b {
                    BooleanOperand::Right
                } else {
                    return Err(BooleanError::UnsupportedResultTopology(
                        "boundary has no contributing source face".into()));
                };
                let coord_axis = side / 2;
                let plane_value = axes[coord_axis][cell[coord_axis]+side%2];
                let source=match operand {
                    BooleanOperand::Left=>a,
                    BooleanOperand::Right=>b,
                };
                let source_face=source.faces.iter().find(|f|f.boundary.iter().all(|id|
                    coordinate(source.vertices[id.0 as usize].position,coord_axis)==plane_value
                )).ok_or_else(||BooleanError::UnsupportedResultTopology(
                    "Boolean boundary did not match a source plane".into()))?.id;
                // Quad winding is outward for each of the six Cartesian normals.
                let q:[(usize,usize,usize);4]=match side {
                    0 => [(i,j,k),(i,j,k+1),(i,j+1,k+1),(i,j+1,k)],
                    1 => [(i+1,j,k),(i+1,j+1,k),(i+1,j+1,k+1),(i+1,j,k+1)],
                    2 => [(i,j,k),(i+1,j,k),(i+1,j,k+1),(i,j,k+1)],
                    3 => [(i,j+1,k),(i,j+1,k+1),(i+1,j+1,k+1),(i+1,j+1,k)],
                    4 => [(i,j,k),(i,j+1,k),(i+1,j+1,k),(i+1,j,k)],
                    _ => [(i,j,k+1),(i+1,j,k+1),(i+1,j+1,k+1),(i,j+1,k+1)],
                };
                let mut ring=Vec::with_capacity(4);
                for corner in q {
                    let id=if let Some(&id)=by_grid.get(&corner){id}else{
                        let id=VertexId(u32::try_from(vertices.len()).map_err(|_|
                            BooleanError::UnsupportedResultTopology("too many Boolean vertices".into()))?);
                        let p=Point3 {x:axes[0][corner.0],y:axes[1][corner.1],
                            z:axes[2][corner.2]};
                        vertices.push(Vertex{id,position:p});
                        by_grid.insert(corner,id);
                        id
                    };
                    ring.push(id);
                }
                let fid=FaceId(u32::try_from(faces.len()).map_err(|_|
                    BooleanError::UnsupportedResultTopology("too many Boolean faces".into()))?);
                faces.push(Face{id:fid,role:FaceRole::Side(fid.0),
                    boundary:ring.clone(),loops:Vec::new()});
                mesh.triangles.push([ring[0].0,ring[1].0,ring[2].0]);
                mesh.triangles.push([ring[0].0,ring[2].0,ring[3].0]);
                mesh.triangle_faces.extend([fid,fid]);
                area+=ds[(coord_axis+1)%3]*ds[(coord_axis+2)%3];
                sources.push(BooleanFaceSource {
                    result_face:fid,operand,source_face,
                    reversed:operation==BooleanOperation::Difference &&
                        operand==BooleanOperand::Right,
                });
            }
        }}}
        if !volume.is_finite() || !area.is_finite() {
            return Err(BooleanError::InvalidGeometry(GeometryError::InvalidEdit(
                "Boolean accumulated mass overflow".into())));
        }
        mesh.vertices=vertices.iter().map(|v|v.position).collect();
        let (edges,coedges,loops,shells)=brep::from_faces(&mut faces).map_err(geometry)?;
        let centroid=Point3 {x:origin[0]+moments[0]/volume,
            y:origin[1]+moments[1]/volume,
            z:origin[2]+moments[2]/volume};
        let mut candidate=Solid {
            vertices, edges, coedges, loops, faces, shells,
            identities:TopologyIdentity::empty(), tolerance, revision:0, mesh,
            bbox:BoundingBox {min:Point3{x:low[0],y:low[1],z:low[2]},
                max:Point3{x:high[0],y:high[1],z:high[2]}},
            mass:MassProperties{volume,surface_area:area,centroid},
        };
        let mut identity=TopologyIdentity::empty();
        identity.initialize(&candidate).map_err(geometry)?;
        candidate.identities=identity;
        candidate.validate().map_err(geometry)?;
        // The result is an independent body. Its names are local to the
        // generated arrangement, not falsely preserved from either input.
        let candidate=candidate.with_feature_key("boolean").map_err(geometry)?;
        Ok(BooleanResult{outcome:BooleanOutcome::Solid(candidate),face_sources:sources})
    }
}
