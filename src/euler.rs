//! Independent, validated inverse Euler edits for the current planar B-rep.
//!
//! `kill_edge_vertex` collapses a valence-two collinear split vertex;
//! `kill_edge_face` removes an interior edge between two coplanar faces.
//! The new body is prepared separately and installed only after complete
//! topology, mesh and identity validation. This first kill path deliberately
//! rebuilds contiguous topology indices and mesh triangles; it is not an
//! in-place O(1) general Euler operator, and remapped runtime handles expire.
use crate::{brep, edit, CoedgeId, EdgeId, Face, FaceId, GeometryError, LoopId, Mesh,
            Solid, TopologyEntity, Vertex, VertexId};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct KillEdgeVertex {
    pub removed_vertex: VertexId,
    pub removed_edge: EdgeId,
    pub retained_edge: EdgeId,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct KillEdgeFace {
    pub removed_edge: EdgeId,
    pub removed_face: FaceId,
    pub retained_face: FaceId,
}

fn invalid(msg:impl Into<String>)->GeometryError {GeometryError::InvalidEdit(msg.into())}

/// Maps original vertex/face IDs to the compacted body IDs; `None` means
/// the entity was retired. Pairs can be used to remap selections by name.
pub(crate) struct RebuildMappings {
    pub vertices: Vec<Option<VertexId>>,
    pub faces: Vec<Option<FaceId>>,
    pub coedge_face: Vec<Option<FaceId>>,
    pub merged_edge: Option<(EdgeId, [VertexId;2])>,
    pub retired_edge: EdgeId,
}

fn distance(a:crate::Point3,b:crate::Point3)->f64{a.sub(b).length_squared().sqrt()}

fn rebuild(original:&Solid, mut vertices:Vec<Vertex>, mut faces:Vec<Face>, map:RebuildMappings)
    ->Result<Solid,GeometryError>
{
    for (i,vertex) in vertices.iter_mut().enumerate(){
        vertex.id=VertexId(u32::try_from(i).map_err(|_|invalid("too many vertices"))?);
    }
    for (i,face) in faces.iter_mut().enumerate(){
        face.id=FaceId(u32::try_from(i).map_err(|_|invalid("too many faces"))?);
        face.loops.clear();
    }
    let (edges,coedges,loops,shells)=brep::from_faces(&mut faces)?;
    let mut candidate=original.clone();
    candidate.vertices=vertices;
    candidate.faces=faces;
    candidate.edges=edges;
    candidate.coedges=coedges;
    candidate.loops=loops;
    candidate.shells=shells;
    candidate.mesh=Mesh {vertices:candidate.vertices.iter().map(|v|v.position).collect(),
        triangles:Vec::new(),triangle_faces:Vec::new()};
    for face in &candidate.faces {
        let triangles=edit::triangulate_ring(&candidate,&face.boundary)?;
        for tri in triangles {
            candidate.mesh.triangles.push(tri);
            candidate.mesh.triangle_faces.push(face.id);
        }
    }
    candidate.identities=original.identities.rebind_after_kill(original,&candidate,&map)?;
    candidate.validate()?;
    Ok(candidate)
}

impl Solid {
    /// Inverse of the *local split* MEV: collapses a valence-two collinear
    /// vertex. The vertex may have been introduced many commands ago. No
    /// dependency on the Undo stack or on it being the most recent slot.
    pub fn kill_edge_vertex(&mut self,vertex:VertexId)->Result<KillEdgeVertex,GeometryError>{
        self.validate()?;
        let old=self.vertices.get(vertex.0 as usize)
            .ok_or_else(||invalid("vertex not found"))?;
        let incident:Vec<_>=self.edges.iter().filter(|edge|edge.start==vertex||edge.end==vertex)
            .collect();
        if incident.len()!=2 {return Err(invalid("KEV requires a valence-two vertex"));}
        let left=if incident[0].start==vertex {incident[0].end}else{incident[0].start};
        let right=if incident[1].start==vertex {incident[1].end}else{incident[1].start};
        if left==right {return Err(invalid("KEV needs two distinct neighboring vertices"));}
        let a=self.vertices[left.0 as usize].position;
        let b=self.vertices[right.0 as usize].position;
        let p=old.position;
        let d=b.sub(a);
        let length=distance(a,b);
        if !length.is_finite() || length<=self.tolerance.length_at(length) {
            return Err(invalid("KEV merged edge is too short"));
        }
        let fraction=p.sub(a).dot(d)/d.length_squared();
        let projection=a.add(d.scale(fraction));
        if fraction<=0.0||fraction>=1.0 || distance(p,projection)>self.tolerance.length_at(length) {
            return Err(invalid("KEV vertex is not inside the straight merged segment"));
        }
        let adj0: BTreeSet<_>=incident[0].coedges.iter().map(|id|self.coedges[id.0 as usize].face).collect();
        let adj1: BTreeSet<_>=incident[1].coedges.iter().map(|id|self.coedges[id.0 as usize].face).collect();
        if adj0!=adj1 || adj0.len()!=2 {
            return Err(invalid("KEV edges must share exactly two incident faces"));
        }
        let mut touched=0;
        for face in &self.faces {
            if face.boundary.contains(&vertex) {
                if !adj0.contains(&face.id) {return Err(invalid("KEV vertex has unexpected face incidence"));}
                touched+=1;
            }
        }
        if touched!=2 {return Err(invalid("KEV expects exactly two affected face loops"));}
        let keep=incident[0].id.min(incident[1].id);
        let drop=incident[0].id.max(incident[1].id);
        let mut vertex_map=vec![None;self.vertices.len()];
        let vertices=self.vertices.iter().filter(|v|v.id!=vertex).enumerate().map(|(new,v)|{
            let id=VertexId(new as u32);
            vertex_map[v.id.0 as usize]=Some(id);
            Vertex{id,position:v.position}
        }).collect::<Vec<_>>();
        let mut faces=self.faces.clone();
        for face in &mut faces {
            if face.boundary.contains(&vertex){face.boundary.retain(|&v|v!=vertex);}
            if face.boundary.len()<3 {return Err(invalid("KEV would collapse a face"));}
            for v in &mut face.boundary { *v=vertex_map[v.0 as usize].ok_or_else(||invalid("missing vertex remap"))?; }
        }
        let e1=vertex_map[left.0 as usize].ok_or_else(||invalid("missing left endpoint"))?;
        let e2=vertex_map[right.0 as usize].ok_or_else(||invalid("missing right endpoint"))?;
        let map=RebuildMappings{
            vertices:vertex_map,
            faces:self.faces.iter().map(|f|Some(f.id)).collect(),
            coedge_face:self.faces.iter().map(|f|Some(f.id)).collect(),
            merged_edge:Some((keep,[e1,e2])), retired_edge:drop,
        };
        let mut candidate=rebuild(self,vertices,faces,map)?;
        candidate.revision=self.revision.checked_add(1).ok_or_else(||invalid("revision overflow"))?;
        *self=candidate;
        Ok(KillEdgeVertex{removed_vertex:vertex,removed_edge:drop,retained_edge:keep})
    }

    /// Inverse MEF for an interior edge separating coplanar, consistently
    /// oriented planar faces. Allows deleting an earlier split diagonal while
    /// unrelated later faces/edges are present elsewhere in the shell.
    pub fn kill_edge_face(&mut self,edge:EdgeId,removed:FaceId)->Result<KillEdgeFace,GeometryError>{
        self.validate()?;
        let e=self.edges.get(edge.0 as usize).ok_or_else(||invalid("edge not found"))?;
        let use_a=&self.coedges[e.coedges[0].0 as usize];
        let use_b=&self.coedges[e.coedges[1].0 as usize];
        let retained=if use_a.face==removed{use_b.face}else if use_b.face==removed{use_a.face}
            else {return Err(invalid("KEF removed face is not incident to edge"));};
        let (outer,inner)=(&self.faces[retained.0 as usize],&self.faces[removed.0 as usize]);
        if outer.role!=inner.role {return Err(invalid("KEF requires equal face surface provenance"));}
        // Each shared edge of the two faces is a potential conflicting join.
        let common=self.edges.iter().filter(|e|{
            let a=self.coedges[e.coedges[0].0 as usize].face;
            let b=self.coedges[e.coedges[1].0 as usize].face;
            (a==retained&&b==removed)||(a==removed&&b==retained)
        }).count();
        if common!=1 {return Err(invalid("KEF requires exactly one shared edge"));}
        let na=crate::face_support_plane(self,retained)?;
        let nb=crate::face_support_plane(self,removed)?;
        let parallel=na.normal.cross(nb.normal).length_squared().sqrt();
        if parallel>self.tolerance.angular.sin() || na.normal.dot(nb.normal)<=0.0 ||
            na.signed_distance(nb.origin)?.abs()>self.tolerance.absolute_length {
            return Err(invalid("KEF faces do not belong to the same oriented plane"));
        }
        // Pick the directed edge in the retained loop, and its reverse in the
        // removed loop. The merged ring follows their non-shared paths.
        fn shared(ring:&[VertexId],a:VertexId,b:VertexId)->Option<usize>{
            (0..ring.len()).find(|&i|{
                let u=ring[i];let v=ring[(i+1)%ring.len()];
                (u==a&&v==b)||(u==b&&v==a)
            })
        }
        let i=shared(&outer.boundary,e.start,e.end).ok_or_else(||invalid("edge absent in retained face"))?;
        let j=shared(&inner.boundary,e.start,e.end).ok_or_else(||invalid("edge absent in removed face"))?;
        let (u,v)=(outer.boundary[i],outer.boundary[(i+1)%outer.boundary.len()]);
        if inner.boundary[j]!=v || inner.boundary[(j+1)%inner.boundary.len()]!=u {
            return Err(invalid("KEF shared edge directions are inconsistent"));
        }
        // Path after the retained shared edge: v -> ... -> u.
        fn path(ring:&[VertexId],start:usize,end:usize)->Vec<VertexId>{
            let mut result=vec![ring[start]];
            let mut cur=start;
            while cur!=end {cur=(cur+1)%ring.len();result.push(ring[cur]);}
            result
        }
        let outside=path(&outer.boundary,(i+1)%outer.boundary.len(),i);
        // Removed face contributes u -> ... -> v, excluding its shared edge.
        let other=path(&inner.boundary,(j+1)%inner.boundary.len(),j);
        let mut joined=outside;
        joined.extend_from_slice(&other[1..other.len()-1]);
        if joined.len()<3||joined.iter().copied().collect::<BTreeSet<_>>().len()!=joined.len() {
            return Err(invalid("KEF merged face is degenerate or self-touching"));
        }
        let mut face_map=vec![None;self.faces.len()];
        let mut faces=Vec::with_capacity(self.faces.len()-1);
        for face in &self.faces {
            if face.id==removed {continue;}
            let next_id=FaceId(faces.len() as u32);
            face_map[face.id.0 as usize]=Some(next_id);
            let mut current=face.clone();current.id=next_id;
            if face.id==retained { current.boundary=joined.clone(); }
            faces.push(current);
        }
        let kept_new=face_map[retained.0 as usize].ok_or_else(||invalid("missing retained face"))?;
        let mut coedge_map=face_map.clone();coedge_map[removed.0 as usize]=Some(kept_new);
        let map=RebuildMappings{
            vertices:self.vertices.iter().map(|v|Some(v.id)).collect(),
            faces:face_map,coedge_face:coedge_map,merged_edge:None,retired_edge:edge,
        };
        let mut candidate=rebuild(self,self.vertices.clone(),faces,map)?;
        candidate.revision=self.revision.checked_add(1).ok_or_else(||invalid("revision overflow"))?;
        *self=candidate;
        Ok(KillEdgeFace{removed_edge:edge,removed_face:removed,retained_face:retained})
    }
}
