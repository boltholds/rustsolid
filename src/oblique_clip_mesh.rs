//! Independent face-owned display and polygonal volume oracle for the
//! restricted oblique halfspace cut. Certified triangles never replace the
//! analytic B-rep or its closed-form mass/centroid.
use std::collections::BTreeMap;
use std::f64::consts::TAU;
use crate::{Curve3,Curve2,EdgeId,FaceId,CoedgeId,GeometryError,
    Mesh,ObliqueHalfspaceEditor,ObliqueHalfspaceSide,Point3,Surface3,
    CylindricalBrep,BrepModel};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn magnitude(x:Point3)->f64{x.x.hypot(x.y).hypot(x.z)}
fn delta(x:Point3,y:Point3)->f64{magnitude(x.sub(y))}

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ObliqueClipMeshCertificate {
    pub revision:u64,
    pub angular_segments:usize,
    pub facet_volume:f64,
    pub source_polygon_volume:f64,
    pub triangle_faces:[usize;3],
    pub area_by_face:[f64;3],
    pub mesh_euler:isize,
    pub max_edge_residual:f64,
}
#[derive(Debug,Clone)]
pub struct CertifiedObliqueClipDisplay {
    pub mesh:Mesh,
    pub certificate:ObliqueClipMeshCertificate,
}
impl ObliqueHalfspaceEditor {
    /// Source-dependent watertight facets of one committed lower/upper clip.
    /// Reject uncommitted history states rather than mislabeling four original
    /// imprinted faces as a three-face clipped solid.
    pub fn tessellate(&self,n:usize)->Result<CertifiedObliqueClipDisplay,GeometryError>{
        if !(12..=4096).contains(&n){
            return Err(invalid("oblique clipped mesh requires 12..=4096 segments"));
        }
        self.model().validate()?;
        let keep=self.applied().ok_or_else(||invalid("no committed oblique halfspace cut"))?;
        make_display(self.source(),self.model(),keep,self.revision(),n)
    }

    /// This validator rejects stale, tampered, non-manifold or incorrectly
    /// face-owned meshes without mutating the committed B-rep.
    pub fn validate_display(&self,display:&CertifiedObliqueClipDisplay)
        ->Result<(),GeometryError>
    {
        if display.certificate.revision!=self.revision() {
            return Err(invalid("clipped mesh belongs to an obsolete source revision"));
        }
        let keep=self.applied().ok_or_else(||invalid("not clipped"))?;
        let actual=certify_clip(self.source(),self.model(),keep,self.revision(),
            &display.mesh,display.certificate.angular_segments)?;
        let expected=display.certificate;
        let volume_eps=(self.source().mass.volume*1e-10).max(1e-10);
        if actual.angular_segments!=expected.angular_segments ||
           actual.triangle_faces!=expected.triangle_faces ||
           actual.mesh_euler!=expected.mesh_euler ||
           (actual.facet_volume-expected.facet_volume).abs()>volume_eps ||
           (actual.source_polygon_volume-expected.source_polygon_volume).abs()>volume_eps ||
           actual.area_by_face.iter().zip(expected.area_by_face.iter())
               .any(|(a,b)|(a-b).abs()>volume_eps) ||
           (actual.max_edge_residual-expected.max_edge_residual).abs()>1e-9 {
            return Err(invalid("clipped mesh certificate is inconsistent with evidence"));
        }
        Ok(())
    }
}

/// The candidate is triangulated and certified *before* publication. This
/// is called by the bounded Boolean transaction and stays independent of
/// the editor's mutation/history bookkeeping.
pub(crate) fn preflight_clip_display(
    source:&CylindricalBrep,model:&BrepModel,keep:ObliqueHalfspaceSide,
)->Result<(),GeometryError>{
    make_display(source,model,keep,0,64)?;
    Ok(())
}

fn make_display(
    source:&CylindricalBrep,model:&BrepModel,
    keep:ObliqueHalfspaceSide,revision:u64,n:usize,
)->Result<CertifiedObliqueClipDisplay,GeometryError>{
    model.validate()?;
    if !(12..=4096).contains(&n) {
        return Err(invalid("clipped mesh segment budget is invalid"));
    }
        if model.faces.len()!=3 || model.edges.len()!=3 ||
            model.coedges.len()!=6 || model.loops.len()!=3 ||
            model.vertices.len()!=2 {
            return Err(invalid("clipped body lacks expected manifold incidence"));
        }
        let mut vertices=Vec::with_capacity(2*n+2);
        for id in [EdgeId(0),EdgeId(1)]{
            let curve=model.edge_curve(id).ok_or_else(||invalid("missing clipped physical edge"))?;
            for i in 0..n {
                vertices.push(curve.evaluate(TAU*i as f64/n as f64)?);
            }
        }
        let ellipse_edge=if keep==ObliqueHalfspaceSide::Negative {EdgeId(1)}else{EdgeId(0)};
        let Curve3::Harmonic(ellipse)=model.edge_curve(ellipse_edge)
            .ok_or_else(||invalid("clipped plane has no ellipse"))?
            else{return Err(invalid("clipped boundary is not harmonic ellipse"));};
        let base=source.frame.origin;
        let top=base.add(source.frame.normal.scale(source.height));
        let [cap0,cap1]=if keep==ObliqueHalfspaceSide::Negative{
            [base,ellipse.center]
        }else{[ellipse.center,top]};
        let center0=vertices.len() as u32;vertices.push(cap0);
        let center1=vertices.len() as u32;vertices.push(cap1);
        let mut triangles=Vec::with_capacity(4*n);
        let mut owners=Vec::with_capacity(4*n);
        for i in 0..n {
            let j=(i+1)%n;
            let (a,b,c,d)=(i as u32,j as u32,(n+i) as u32,(n+j) as u32);
            for (tri,face) in [
                ([center0,b,a],FaceId(0)),
                ([center1,c,d],FaceId(1)),
                ([a,b,d],FaceId(2)),
                ([a,d,c],FaceId(2)),
            ] {
                triangles.push(tri);
                owners.push(face);
            }
        }
        let mesh=Mesh{vertices,triangles,triangle_faces:owners};
        let cert=certify_clip(source,model,keep,revision,&mesh,n)?;
        Ok(CertifiedObliqueClipDisplay{mesh,certificate:cert})
    }



fn certify_clip(
    source:&CylindricalBrep,model:&BrepModel,keep:ObliqueHalfspaceSide,
    revision:u64,mesh:&Mesh,n:usize,
)->Result<ObliqueClipMeshCertificate,GeometryError>
{
    model.validate()?;
    if !(12..=4096).contains(&n) || mesh.vertices.len()!=2*n+2 ||
       mesh.triangles.len()!=4*n ||
       mesh.triangle_faces.len()!=4*n ||
       mesh.vertices.iter().any(|p|!p.is_finite()) {
        return Err(invalid("clipped display arrays or segment budget are invalid"));
    }
    let global_scale=source.bbox.max.x.abs().max(source.bbox.min.x.abs())
        .max(source.bbox.max.y.abs()).max(source.bbox.min.y.abs())
        .max(source.bbox.max.z.abs()).max(source.bbox.min.z.abs()).max(1.0);
    let tolerance=source.tolerance.length_at(source.radius.max(source.height))
        .max(64.0*f64::EPSILON*global_scale);
    let mut residual=0.0f64;
    let ellipse_ring=if keep==ObliqueHalfspaceSide::Negative{1usize}else{0};
    let ellipse_edge=EdgeId(ellipse_ring as u32);
    let ellipse=model.edge_curve(ellipse_edge)
        .ok_or_else(||invalid("missing clipped edge"))?;
    if !matches!(ellipse,Curve3::Harmonic(_)){
        return Err(invalid("clipped ellipsoid ring is not harmonic"));
    }
    for ring in 0..2 {
        let curve=model.edge_curve(EdgeId(ring as u32))
            .ok_or_else(||invalid("missing analytical clipped ring"))?;
        for i in 0..n {
            let param=TAU*i as f64/n as f64;
            let actual=mesh.vertices[ring*n+i];
            let gap=delta(actual,curve.evaluate(param)?);
            residual=residual.max(gap);
            if gap>tolerance {
                return Err(invalid("clipped mesh vertices deviate from analytical physical ring"));
            }
            if ring==ellipse_ring {
                // Each physical ellipse sample reconstructs on BOTH incident
                // surface charts, with opposite coedge parameter directions.
                let edge=&model.edges[ellipse_edge.0 as usize];
                for co_id in edge.coedges {
                    let co=&model.coedges[co_id.0 as usize];
                    let pcurve=model.coedge_pcurve(co_id)
                        .ok_or_else(||invalid("missing clipped elliptical coedge pcurve"))?;
                    if !matches!(pcurve,Curve2::Harmonic(_) | Curve2::CylinderSection(_)) {
                        return Err(invalid("unsupported elliptical coedge parameter curve"));
                    }
                    let t=if co.reversed {TAU-param}else{param};
                    let uv=pcurve.evaluate(t)?;
                    let surface=model.face_surface(co.face)
                        .ok_or_else(||invalid("missing clipped ellipse support surface"))?;
                    let reconstructed=surface.evaluate(uv.u,uv.v)?.point;
                    let error=delta(actual,reconstructed);
                    residual=residual.max(error);
                    if error>tolerance{
                        return Err(invalid("clipped ellipse face UV disagrees with 3D carrier"));
                    }
                }
            }
        }
        if delta(mesh.vertices[ring*n],model.vertices[ring].position)>tolerance{
            return Err(invalid("clipped ring seam vertex does not match physical vertex"));
        }
    }
    let (base,cut)=(source.frame.origin,model.edge_curve(ellipse_edge)
        .ok_or_else(||invalid("missing elliptical carrier"))?);
    let Curve3::Harmonic(ellipse)=cut else{unreachable!()};
    let top=base.add(source.frame.normal.scale(source.height));
    let (center0,center1)=if keep==ObliqueHalfspaceSide::Negative{
        (base,ellipse.center)
    }else{(ellipse.center,top)};
    if delta(center0,mesh.vertices[2*n])>tolerance ||
       delta(center1,mesh.vertices[2*n+1])>tolerance {
        return Err(invalid("clipped cap centers depart from source plane support"));
    }
    let mut incidences=BTreeMap::<(u32,u32),Vec<(FaceId,bool)>>::new();
    let mut counts=[0usize;3];
    let mut areas=[0.0f64;3];
    let mut total_six=0.0;
    for (triangle,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        if face.0>2 || triangle.iter().any(|v|*v as usize>=mesh.vertices.len()) ||
           triangle[0]==triangle[1] || triangle[0]==triangle[2] ||
           triangle[1]==triangle[2] {
            return Err(invalid("clipped facet geometry has invalid owner or vertex"));
        }
        let f=face.0 as usize;
        let pts=triangle.map(|index|index as usize);
        let valid=match f {
            0=>pts.contains(&(2*n)) && pts.iter().filter(|&&j|j<n).count()==2,
            1=>pts.contains(&(2*n+1)) &&
                pts.iter().filter(|&&j|j>=n && j<2*n).count()==2,
            2=>pts.iter().all(|&j|j<2*n) &&
                pts.iter().any(|&j|j<n) &&
                pts.iter().any(|&j|j>=n),
            _=>unreachable!(),
        };
        if !valid{return Err(invalid("triangle crosses an unrelated clipped face interval"));}
        let a=mesh.vertices[pts[0]];
        let b=mesh.vertices[pts[1]];
        let c=mesh.vertices[pts[2]];
        let cross=b.sub(a).cross(c.sub(a));
        let magnitude=magnitude(cross);
        if !magnitude.is_finite() || magnitude<=0.0 {
            return Err(invalid("clipped facet has nonpositive area"));
        }
        let oriented=if f<=1 {
            let Surface3::Plane(p)=model.face_surface(*face)
                .ok_or_else(||invalid("missing clipped plane cap"))?
                else{return Err(invalid("clipped facet cap surface is not planar"));};
            cross.dot(p.frame.normal)
        }else{
            let midpoint=a.add(b).add(c).scale(1.0/3.0);
            let q=midpoint.sub(base);
            let radial=q.sub(source.frame.normal.scale(q.dot(source.frame.normal)));
            cross.dot(radial)
        };
        if oriented<=0.0{
            return Err(invalid("clipped display normal disagrees with exterior material sense"));
        }
        counts[f]+=1;
        areas[f]+=magnitude*0.5;
        total_six+=a.sub(base).dot(b.sub(base).cross(c.sub(base)));
        for (u,v) in [(triangle[0],triangle[1]),
            (triangle[1],triangle[2]),(triangle[2],triangle[0])] {
            incidences.entry((u.min(v),u.max(v))).or_default().push((*face,u<v));
        }
    }
    if counts!=[n,n,2*n] {
        return Err(invalid("clipped display does not cover all three B-rep faces"));
    }
    if incidences.values().any(|v|v.len()!=2||v[0].1==v[1].1){
        return Err(invalid("clipped shell faceting is open or inconsistently wound"));
    }
    let mesh_euler=mesh.vertices.len() as isize-
        incidences.len() as isize+mesh.triangles.len() as isize;
    if mesh_euler!=2 || model.euler_characteristic()!=2 {
        return Err(invalid("clipped display Euler invariant is not two"));
    }
    // The elliptical cut is one physical ring with opposite face owners:
    // (side, elliptic cap), independent of the cylinder chart seam.
    for i in 0..n{
        let j=(i+1)%n;
        let a=(ellipse_ring*n+i) as u32;
        let b=(ellipse_ring*n+j) as u32;
        let uses=incidences.get(&(a.min(b),a.max(b)))
            .ok_or_else(||invalid("clipped ellipse edge is not meshed"))?;
        let plane=if keep==ObliqueHalfspaceSide::Negative{FaceId(1)}
            else{FaceId(0)};
        if uses.len()!=2 || uses[0].0==uses[1].0 ||
            !uses.iter().any(|(f,_)|*f==FaceId(2)) ||
            !uses.iter().any(|(f,_)|*f==plane) {
            return Err(invalid("elliptic cut mesh ownership is inconsistent"));
        }
    }
    let mean=ellipse.center.sub(base).dot(source.frame.normal);
    let axial_height=if keep==ObliqueHalfspaceSide::Negative{
        mean
    }else{source.height-mean};
    let prism=0.5*(n as f64)*source.radius.powi(2)*
        (TAU/n as f64).sin()*axial_height;
    let mesh_volume=total_six/6.0;
    let volume_budget=(prism.abs()*1e-8).max(
        tolerance*source.radius.max(source.height).powi(2)*16.0);
    if !mesh_volume.is_finite() || mesh_volume<=0.0 ||
       (mesh_volume-prism).abs()>volume_budget ||
       mesh_volume>model.mass.volume*(1.0+1e-9) {
        return Err(invalid("clipped triangle volume differs from independent polygon oracle"));
    }
    Ok(ObliqueClipMeshCertificate{
        revision,
        angular_segments:n,
        facet_volume:mesh_volume,
        source_polygon_volume:prism,
        triangle_faces:counts,area_by_face:areas,
        mesh_euler,max_edge_residual:residual,
    })
}
