//! Face-owned watertight display for a source-preserving oblique cylinder
//! imprint. One canonical sample ring is used for the physical elliptical
//! cut edge, shared by both cylindrical faces; no seam duplication.
use std::collections::BTreeMap;
use std::f64::consts::TAU;
use crate::{ClosedCylinderImprint, Curve3, Curve2, EdgeId, FaceId,
    GeometryError, Mesh, Point3, Surface3};

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn norm(a:Point3)->f64{a.x.hypot(a.y).hypot(a.z)}
fn distance(a:Point3,b:Point3)->f64{norm(a.sub(b))}
fn face_index(f:FaceId)->Result<usize,GeometryError>{
    if f.0<4 {Ok(f.0 as usize)}
    else{Err(invalid("imprinted mesh has unknown face owner"))}
}
const MIN_SEGMENTS:usize=12;
const MAX_SEGMENTS:usize=4096;

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ImprintMeshCertificate {
    pub source_revision:u64,
    pub angular_segments:usize,
    pub triangle_count:usize,
    pub vertex_count:usize,
    pub facet_volume:f64,
    pub inscribed_prism_volume:f64,
    pub max_analytic_edge_residual:f64,
    pub face_area:[f64;4],
    pub mesh_euler:isize,
}
#[derive(Debug,Clone)]
pub struct CertifiedImprintDisplay {
    pub mesh:Mesh,
    pub certificate:ImprintMeshCertificate,
}

impl ClosedCylinderImprint {
    /// This display is defined only for a source cylinder with one completed
    /// imprint. It does not change the B-rep or historical revision.
    pub fn tessellate(&self,segments:usize)->Result<CertifiedImprintDisplay,GeometryError>{
        if !(MIN_SEGMENTS..=MAX_SEGMENTS).contains(&segments) {
            return Err(invalid("imprint tessellation requires 12..=4096 segments"));
        }
        let model=self.model();
        model.validate()?;
        self.source_cylinder().validate()?;
        if model.vertices.len()!=3 || model.edges.len()!=5 ||
            model.coedges.len()!=10 || model.faces.len()!=4 ||
            model.loops.len()!=4 || model.shells.len()!=1 {
            return Err(invalid("imprint display requires a source-certified closed side split"));
        }
        let mut vertices=Vec::with_capacity(3*segments+2);
        for eid in [EdgeId(0),EdgeId(4),EdgeId(1)] {
            let curve=model.edge_curve(eid).ok_or_else(||invalid("missing physical ring carrier"))?;
            for i in 0..segments {
                let angle=TAU*(i as f64)/(segments as f64);
                vertices.push(curve.evaluate(angle)?);
            }
        }
        let bottom_center=vertices.len() as u32;
        vertices.push(self.source_cylinder().frame.origin);
        let top_center=vertices.len() as u32;
        vertices.push(self.source_cylinder().frame.origin.add(
            self.source_cylinder().frame.normal.scale(self.source_cylinder().height)));
        let mut triangles=Vec::with_capacity(6*segments);
        let mut triangle_faces=Vec::with_capacity(6*segments);
        let n=segments as u32;
        for i in 0..n {
            let j=(i+1)%n;
            let bottom_i=i;let bottom_j=j;
            let ellipse_i=n+i;let ellipse_j=n+j;
            let top_i=2*n+i;let top_j=2*n+j;
            for (tri,face) in [
                ([bottom_center,bottom_j,bottom_i],FaceId(0)),
                ([top_center,top_i,top_j],FaceId(1)),
                ([bottom_i,bottom_j,ellipse_j],FaceId(2)),
                ([bottom_i,ellipse_j,ellipse_i],FaceId(2)),
                ([ellipse_i,ellipse_j,top_j],FaceId(3)),
                ([ellipse_i,top_j,top_i],FaceId(3)),
            ] {
                triangles.push(tri);
                triangle_faces.push(face);
            }
        }
        let mesh=Mesh{vertices,triangles,triangle_faces};
        let certificate=certify_mesh(self,&mesh,segments)?;
        Ok(CertifiedImprintDisplay{mesh,certificate})
    }

    /// Recompute independent mesh checks, including face ownership, manifold
    /// incidence, edge sample provenance and signed inscribed-prism mass.
    pub fn validate_display(&self,display:&CertifiedImprintDisplay)
        ->Result<(),GeometryError>
    {
        if display.certificate.source_revision!=self.revision() {
            return Err(invalid("imprint display belongs to a stale revision"));
        }
        let found=certify_mesh(self,&display.mesh,display.certificate.angular_segments)?;
        let declared=display.certificate;
        let scale=self.source_cylinder().mass.volume.max(1.0);
        if found.angular_segments!=declared.angular_segments ||
            found.triangle_count!=declared.triangle_count ||
            found.vertex_count!=declared.vertex_count ||
            found.mesh_euler!=declared.mesh_euler ||
            (found.facet_volume-declared.facet_volume).abs()>1e-10*scale ||
            (found.inscribed_prism_volume-declared.inscribed_prism_volume).abs()>1e-10*scale ||
            found.face_area.iter().zip(declared.face_area.iter()).any(|(a,b)|
                (a-b).abs()>1e-10*scale) ||
            (found.max_analytic_edge_residual-declared.max_analytic_edge_residual).abs()>
                self.source_cylinder().tolerance.absolute_length {
            return Err(invalid("imprint mesh certificate is inconsistent with triangles"));
        }
        Ok(())
    }
}

fn certify_mesh(owner:&ClosedCylinderImprint,mesh:&Mesh,n:usize)
    ->Result<ImprintMeshCertificate,GeometryError>
{
    let model=owner.model();
    let source=owner.source_cylinder();
    model.validate()?;
    if !(MIN_SEGMENTS..=MAX_SEGMENTS).contains(&n) ||
       mesh.vertices.len()!=3*n+2 || mesh.triangles.len()!=6*n ||
       mesh.triangle_faces.len()!=mesh.triangles.len() ||
       mesh.vertices.iter().any(|p|!p.is_finite()) {
        return Err(invalid("imprint facet shape, budget or vertex data invalid"));
    }
    if model.mass!=source.mass || model.bbox!=source.bbox {
        return Err(invalid("imprint display changed analytic source mass"));
    }
    let source_extent=source.radius.max(source.height);
    let local_epsilon=source.tolerance.length_at(source_extent);
    let world_max=model.bbox.max.x.abs().max(model.bbox.max.y.abs())
        .max(model.bbox.max.z.abs())
        .max(model.bbox.min.x.abs()).max(model.bbox.min.y.abs())
        .max(model.bbox.min.z.abs()).max(1.0);
    let eps=local_epsilon.max(64.0*f64::EPSILON*world_max);
    let side=model.face_surface(FaceId(2))
        .ok_or_else(||invalid("missing imprinted lower support"))?;
    if side!=model.face_surface(FaceId(3))
        .ok_or_else(||invalid("missing imprinted upper support"))? {
        return Err(invalid("imprinted side patches have incompatible surfaces"));
    }
    if !matches!(side,Surface3::Cylinder(_)) ||
       !matches!(model.edge_curve(EdgeId(4)),Some(Curve3::Harmonic(_))) ||
       !matches!(model.coedge_pcurve(crate::CoedgeId(8)),
           Some(Curve2::CylinderSection(_))) ||
       !matches!(model.coedge_pcurve(crate::CoedgeId(9)),
           Some(Curve2::CylinderSection(_))) {
        return Err(invalid("imprint display lacks canonical analytic ellipse owners"));
    }
    let mut maximal_residual=0.0f64;
    for (eid,ring,vertex) in [
        (EdgeId(0),0usize,Some(crate::VertexId(0))),
        (EdgeId(4),1usize,Some(crate::VertexId(2))),
        (EdgeId(1),2usize,Some(crate::VertexId(1))),
    ] {
        let curve=model.edge_curve(eid)
            .ok_or_else(||invalid("mesh ring is missing its source edge"))?;
        for i in 0..n {
            let param=TAU*i as f64/n as f64;
            let actual=mesh.vertices[ring*n+i];
            let expected=curve.evaluate(param)?;
            let gap=distance(actual,expected);
            maximal_residual=maximal_residual.max(gap);
            if gap>eps{return Err(invalid("mesh ring is detached from physical analytic edge"));}
            if ring==1 {
                let Surface3::Cylinder(_) = side else{unreachable!()};
                let [c0,c1]=[crate::CoedgeId(8),crate::CoedgeId(9)];
                for (coedge,t) in [(c0,param),(c1,TAU-param)] {
                    let g=model.geometry.coedges[coedge.0 as usize];
                    let trim=model.geometry.curves2.get(g.pcurve.0 as usize)
                        .ok_or_else(||invalid("missing ellipse pcurve"))?;
                    let q=trim.evaluate(t)?;
                    let on_side=side.evaluate(q.u,q.v)?.point;
                    let residual=distance(actual,on_side);
                    maximal_residual=maximal_residual.max(residual);
                    if residual>eps {
                        return Err(invalid("mesh ellipse does not match owned UV/3D carriers"));
                    }
                }
            }
        }
        let origin=model.vertices[vertex.unwrap().0 as usize].position;
        if distance(mesh.vertices[ring*n],origin)>eps {
            return Err(invalid("mesh seam ring does not reuse physical B-rep vertex"));
        }
    }
    let bc=source.frame.origin;
    let tc=bc.add(source.frame.normal.scale(source.height));
    if distance(mesh.vertices[3*n],bc)>eps ||
       distance(mesh.vertices[3*n+1],tc)>eps {
        return Err(invalid("cap mesh center differs from analytic cylinder center"));
    }

    let mut incidences=BTreeMap::<(u32,u32),Vec<(FaceId,bool)>>::new();
    let mut area=[0.0f64;4];
    let mut triangle_counts=[0usize;4];
    let mut volume_six=0.0;
    for (tri,face) in mesh.triangles.iter().zip(mesh.triangle_faces.iter()) {
        let face_idx=face_index(*face)?;
        if tri.iter().any(|x|(*x as usize)>=mesh.vertices.len()) ||
            tri[0]==tri[1] || tri[1]==tri[2] || tri[2]==tri[0] {
            return Err(invalid("imprint mesh has malformed triangles"));
        }
        let indices:Vec<usize>=tri.iter().map(|i|*i as usize).collect();
        let allowed=match face_idx {
            0=>indices.iter().filter(|&&j|j<n).count()==2 &&
                indices.contains(&(3*n)),
            1=>indices.iter().filter(|&&j|j>=2*n && j<3*n).count()==2 &&
                indices.contains(&(3*n+1)),
            2=>indices.iter().all(|&j|j<2*n) &&
                indices.iter().any(|&j|j<n) &&
                indices.iter().any(|&j|j>=n),
            3=>indices.iter().all(|&j|j>=n && j<3*n) &&
                indices.iter().any(|&j|j<2*n) &&
                indices.iter().any(|&j|j>=2*n),
            _=>unreachable!(),
        };
        if !allowed {return Err(invalid("triangle owner disagrees with analytic face interval"));}
        let a=mesh.vertices[indices[0]];
        let b=mesh.vertices[indices[1]];
        let c=mesh.vertices[indices[2]];
        let cross=b.sub(a).cross(c.sub(a));
        let twice=norm(cross);
        if !twice.is_finite() || twice<=0.0 {
            return Err(invalid("degenerate imprint facet"));
        }
        let orientation=if face_idx==0 {-cross.dot(source.frame.normal)}
            else if face_idx==1 {cross.dot(source.frame.normal)}
            else {
                let midpoint=a.add(b).add(c).scale(1.0/3.0);
                let axial=midpoint.sub(source.frame.origin).dot(source.frame.normal);
                let radial=midpoint.sub(source.frame.origin)
                    .sub(source.frame.normal.scale(axial));
                cross.dot(radial)
            };
        if orientation<=0.0 {
            return Err(invalid("imprint facet is inconsistent with outward material normal"));
        }
        area[face_idx]+=twice*0.5;
        triangle_counts[face_idx]+=1;
        volume_six+=a.sub(bc).dot(b.sub(bc).cross(c.sub(bc)));
        for (i,j) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])] {
            incidences.entry((i.min(j),i.max(j))).or_default().push((*face,i<j));
        }
    }
    if triangle_counts!=[n,n,2*n,2*n] {
        return Err(invalid("mesh does not cover all four source B-rep faces"));
    }
    if incidences.values().any(|owners| owners.len()!=2 || owners[0].1==owners[1].1) {
        return Err(invalid("imprint mesh is not closed and consistently oriented"));
    }
    let e=incidences.len() as isize;
    let chi=mesh.vertices.len() as isize-e+mesh.triangles.len() as isize;
    if chi!=model.euler_characteristic() || chi!=2 {
        return Err(invalid("imprint display Euler characteristic differs from B-rep"));
    }
    // Every segment of the elliptical cut ring must be seen once by each
    // of the two new side faces, in opposite orientation. None is a seam.
    for i in 0..n {
        let j=(i+1)%n;
        let ei=(n+i) as u32;let ej=(n+j) as u32;
        let owners=incidences.get(&(ei.min(ej),ei.max(ej)))
            .ok_or_else(||invalid("physical ellipse mesh segment is missing"))?;
        if owners.len()!=2 || owners[0].0==owners[1].0 ||
           owners[0].1==owners[1].1 ||
           !owners.iter().any(|x|x.0==FaceId(2)) ||
           !owners.iter().any(|x|x.0==FaceId(3)) {
            return Err(invalid("physical ellipse facet edge lacks opposite face owners"));
        }
    }
    let facet_volume=volume_six/6.0;
    let inscribed_prism_volume=0.5*(n as f64)*source.radius.powi(2)*
        (TAU/n as f64).sin()*source.height;
    let bound=(inscribed_prism_volume.abs()*1e-8).max(eps*source_extent.powi(2)*32.0);
    if !facet_volume.is_finite() || facet_volume<=0.0 ||
       (facet_volume-inscribed_prism_volume).abs()>bound ||
       facet_volume>source.mass.volume*(1.0+1e-9) {
        return Err(invalid("mesh signed volume is inconsistent with source/inscribed oracle"));
    }
    Ok(ImprintMeshCertificate{
        source_revision:owner.revision(),
        angular_segments:n,
        triangle_count:mesh.triangles.len(),
        vertex_count:mesh.vertices.len(),
        facet_volume,inscribed_prism_volume,max_analytic_edge_residual:maximal_residual,
        face_area:area,mesh_euler:chi,
    })
}
