//! Combinatorial G4b radial-shell preflight.
//!
//! Input halfedges must already be certified against source B-rep geometry.
//! Vertex equivalence is generated ONLY by opposite radial coedges and
//! consecutive source-face loop uses, never by proximity-based welding.
//! Geometric collisions between DISTINCT vertex fans are Held. This audit
//! cannot replace construction/validation of a canonical BrepModel or an
//! independent volume oracle.

use crate::Point3;
use std::collections::{BTreeMap, BTreeSet};

const MAX_RADIAL_HALFEDGES: usize = 4096;

#[derive(Debug,Clone,PartialEq,Eq)]
pub struct RadialTopologyCertificate {
    pub vertices: usize,
    pub edges: usize,
    pub material_face_cells: usize,
    pub hole_cycles: usize,
    pub shell_components: usize,
    pub euler_characteristic: isize,
    pub total_genus: usize,
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum RadialTopologyHold {
    BudgetExceeded,
    MissingLoop,
    OpenLoop,
    UnpairedEdge,
    NonManifoldEdge,
    NonManifoldVertex,
    CoincidentUncertifiedVertices,
    InvalidEulerOrGenus,
}

#[derive(Debug,Clone,Copy)]
pub(crate) struct RadialHalfedge {
    pub start: Point3,
    pub end: Point3,
    pub face: usize,
}
#[derive(Debug,Clone)]
pub(crate) struct RadialLoop {
    pub face: usize,
    pub uses: Vec<usize>,
    pub inner: bool,
}
#[derive(Debug)]
struct SetUnion {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl SetUnion {
    fn new(n:usize)->Self{
        Self{parent:(0..n).collect(),rank:vec![0;n]}
    }
    fn root(&mut self,x:usize)->usize {
        if self.parent[x]!=x {
            let parent=self.parent[x];
            self.parent[x]=self.root(parent);
        }
        self.parent[x]
    }
    fn join(&mut self,a:usize,b:usize){
        let mut x=self.root(a);let mut y=self.root(b);
        if x==y{return;}
        if self.rank[x]<self.rank[y]{std::mem::swap(&mut x,&mut y);}
        self.parent[y]=x;
        if self.rank[x]==self.rank[y]{self.rank[x]+=1;}
    }
}
fn close(a:Point3,b:Point3,eps:f64)->bool {
    a.is_finite()&&b.is_finite()&&
    (a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)<=eps
}

/// This topology-only input is private to the RustSolid kernel. Public
/// callers must pass through source-revalidated G4b patch matching.
pub(crate) fn certify_radial_shell(
    uses:&[RadialHalfedge],
    loops:&[RadialLoop],
    radial_pairs:&[[usize;2]],
    face_count:usize,
    eps:f64,
)->Result<RadialTopologyCertificate,RadialTopologyHold>{
    if uses.is_empty()||uses.len()>MAX_RADIAL_HALFEDGES ||
        uses.len()!=radial_pairs.len()*2 || face_count==0 ||
        !eps.is_finite()||eps<=0.0
    {
        return Err(RadialTopologyHold::BudgetExceeded);
    }
    let mut next=vec![usize::MAX;uses.len()];
    let mut twin=vec![usize::MAX;uses.len()];
    let mut faces=SetUnion::new(face_count);
    let mut vertices=SetUnion::new(uses.len()*2);
    let mut outer=vec![0usize;face_count];
    let mut holes=vec![0usize;face_count];

    for lp in loops {
        if lp.face>=face_count||lp.uses.is_empty() {
            return Err(RadialTopologyHold::MissingLoop);
        }
        if lp.inner {holes[lp.face]+=1;}else{outer[lp.face]+=1;}
        for (index,&id) in lp.uses.iter().enumerate() {
            let j=lp.uses[(index+1)%lp.uses.len()];
            let (Some(a),Some(b))=(uses.get(id),uses.get(j)) else{
                return Err(RadialTopologyHold::MissingLoop);
            };
            if a.face!=lp.face||b.face!=lp.face||next[id]!=usize::MAX{
                return Err(RadialTopologyHold::MissingLoop);
            }
            if !close(a.end,b.start,eps){
                return Err(RadialTopologyHold::OpenLoop);
            }
            next[id]=j;
            vertices.join(2*id+1,2*j);
        }
    }
    if outer.iter().any(|n|*n!=1)||next.contains(&usize::MAX){
        return Err(RadialTopologyHold::MissingLoop);
    }

    for pair in radial_pairs {
        let [u,v]=*pair;
        let (Some(a),Some(b))=(uses.get(u),uses.get(v)) else {
            return Err(RadialTopologyHold::UnpairedEdge);
        };
        if u==v||twin[u]!=usize::MAX||twin[v]!=usize::MAX||
            a.face==b.face
        {
            // A pair between two distinct cells on one original face is
            // not a true physical radial twin of a new Boolean cut.
            return Err(RadialTopologyHold::NonManifoldEdge);
        }
        if !close(a.start,b.end,eps)||!close(a.end,b.start,eps){
            return Err(RadialTopologyHold::NonManifoldEdge);
        }
        twin[u]=v;twin[v]=u;
        faces.join(a.face,b.face);
        vertices.join(2*u,2*v+1);
        vertices.join(2*u+1,2*v);
    }
    if twin.contains(&usize::MAX){
        return Err(RadialTopologyHold::UnpairedEdge);
    }

    // Equivalence is purely topological, independently checked in XYZ.
    let mut representatives=BTreeMap::<usize,Point3>::new();
    let mut outgoing=BTreeMap::<usize,Vec<usize>>::new();
    for i in 0..uses.len() {
        for (is_end,p) in [(false,uses[i].start),(true,uses[i].end)] {
            let root=vertices.root(2*i+usize::from(is_end));
            if !p.is_finite(){return Err(RadialTopologyHold::NonManifoldVertex);}
            if let Some(existing)=representatives.get(&root){
                if !close(*existing,p,eps){
                    return Err(RadialTopologyHold::NonManifoldVertex);
                }
            }else{representatives.insert(root,p);}
            if !is_end {outgoing.entry(root).or_default().push(i);}
        }
    }
    // A pinched surface may have disjoint incidence fans at identical XYZ.
    // Without a physical source identity proof, it must not become one
    // candidate vertex merely because the coordinates coincide.
    let entries=representatives.iter().collect::<Vec<_>>();
    for (i,(root,a)) in entries.iter().enumerate(){
        for (other,b) in entries.iter().skip(i+1) {
            if *root!=*other && close(**a,**b,eps) {
                return Err(RadialTopologyHold::CoincidentUncertifiedVertices);
            }
        }
    }
    for (root,arms) in &outgoing {
        let mut visited=BTreeSet::new();
        let start=arms[0];
        let mut at=start;
        loop {
            if !visited.insert(at){
                if at!=start {return Err(RadialTopologyHold::NonManifoldVertex);}
                break;
            }
            // Rotate around this corner with the opposite radial use,
            // then advance in the adjacent face's directed boundary.
            at=next[twin[at]];
            if vertices.root(2*at)!=*root ||
                visited.len()>arms.len()
            {
                return Err(RadialTopologyHold::NonManifoldVertex);
            }
        }
        if visited.len()!=arms.len(){
            return Err(RadialTopologyHold::NonManifoldVertex);
        }
    }

    let mut components=BTreeMap::<usize,(BTreeSet<usize>,BTreeSet<usize>,BTreeSet<usize>,usize)>::new();
    for face in 0..face_count {
        let root=faces.root(face);
        let entry=components.entry(root).or_default();
        entry.2.insert(face);
        entry.3+=holes[face];
    }
    for (i,halfedge) in uses.iter().enumerate(){
        let root=faces.root(halfedge.face);
        components.get_mut(&root).unwrap().0.insert(vertices.root(2*i));
    }
    for (edge,[u,_]) in radial_pairs.iter().enumerate(){
        let root=faces.root(uses[*u].face);
        components.get_mut(&root).unwrap().1.insert(edge);
    }
    let mut genus_total=0usize;
    let mut euler_total=0isize;
    for (_, (vs,es,fs,holes)) in &components {
        let chi=vs.len() as isize-es.len() as isize+
            fs.len() as isize-*holes as isize;
        let deficit=2-chi;
        if deficit<0||deficit%2!=0 {
            return Err(RadialTopologyHold::InvalidEulerOrGenus);
        }
        genus_total+=deficit as usize/2;
        euler_total+=chi;
    }
    Ok(RadialTopologyCertificate {
        vertices:representatives.len(),
        edges:radial_pairs.len(),
        material_face_cells:face_count,
        hole_cycles:holes.iter().sum(),
        shell_components:components.len(),
        euler_characteristic:euler_total,
        total_genus:genus_total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
    fn tetra(offset:f64)->(Vec<RadialHalfedge>,Vec<RadialLoop>,Vec<[usize;2]>) {
        let v=[
            p(offset,0.,0.),p(offset+1.,0.,0.),
            p(offset,1.,0.),p(offset,0.,1.),
        ];
        // Consistent outward orientation on the four source triangular faces.
        let faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]];
        let mut half=Vec::new();
        let mut loops=Vec::new();
        let mut pairs=BTreeMap::<(usize,usize),Vec<usize>>::new();
        for (face_idx,vertices) in faces.iter().enumerate() {
            let mut uses=Vec::new();
            for j in 0..3 {
                let a=vertices[j];let b=vertices[(j+1)%3];
                let idx=half.len();
                half.push(RadialHalfedge{start:v[a],end:v[b],face:face_idx});
                uses.push(idx);
                pairs.entry((a.min(b),a.max(b))).or_default().push(idx);
            }
            loops.push(RadialLoop{face:face_idx,uses,inner:false});
        }
        let radial=pairs.into_values().map(|u|[u[0],u[1]]).collect();
        (half,loops,radial)
    }
    #[test]
    fn tetrahedron_has_one_closed_oriented_genus_zero_shell(){
        let (half,loops,edges)=tetra(0.);
        let cert=certify_radial_shell(&half,&loops,&edges,4,1e-9).unwrap();
        assert_eq!(cert.vertices,4);
        assert_eq!(cert.edges,6);
        assert_eq!(cert.material_face_cells,4);
        assert_eq!(cert.euler_characteristic,2);
        assert_eq!(cert.shell_components,1);
        assert_eq!(cert.total_genus,0);
    }
    #[test]
    fn duplicate_nonmanifold_oriented_edge_is_held(){
        let (half,loops,mut edges)=tetra(0.);
        edges[0][1]=edges[0][0];
        assert_eq!(certify_radial_shell(&half,&loops,&edges,4,1e-9),
            Err(RadialTopologyHold::NonManifoldEdge));
    }
    #[test]
    fn face_with_open_loop_cannot_claim_watertight_shell(){
        let (mut half,loops,edges)=tetra(0.);
        half[0].end.x+=0.5;
        assert_eq!(certify_radial_shell(&half,&loops,&edges,4,1e-9),
            Err(RadialTopologyHold::OpenLoop));
    }
    #[test]
    fn vertex_only_contact_between_two_shells_cannot_weld_by_xyz(){
        let (mut a,mut al,mut ae)=tetra(0.);
        let (b,bl,be)=tetra(0.); // Same coordinates, DIFFERENT topology.
        let edge_shift=a.len();
        let face_shift=al.len();
        for mut h in b {h.face+=face_shift;a.push(h);}
        for mut lp in bl {lp.face+=face_shift;
            for u in &mut lp.uses {*u+=edge_shift;}
            al.push(lp);
        }
        for [u,v] in be {ae.push([u+edge_shift,v+edge_shift]);}
        assert_eq!(certify_radial_shell(&a,&al,&ae,8,1e-9),
            Err(RadialTopologyHold::CoincidentUncertifiedVertices));
    }
    #[test]
    fn disjoint_multiple_components_keep_independent_euler_genus(){
        let (mut a,mut al,mut ae)=tetra(0.);
        let (b,bl,be)=tetra(3.);
        let edge_shift=a.len();let face_shift=al.len();
        for mut h in b {h.face+=face_shift;a.push(h);}
        for mut lp in bl {lp.face+=face_shift;
            for u in &mut lp.uses {*u+=edge_shift;}
            al.push(lp);
        }
        for [u,v] in be {ae.push([u+edge_shift,v+edge_shift]);}
        let cert=certify_radial_shell(&a,&al,&ae,8,1e-9).unwrap();
        assert_eq!((cert.vertices,cert.edges,cert.material_face_cells),
            (8,12,8));
        assert_eq!((cert.shell_components,cert.euler_characteristic,cert.total_genus),
            (2,4,0));
    }
}
