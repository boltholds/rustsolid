//! A bounded, independently certified display mesher for connected analytic
//! B-reps with planar cap regions containing multiple circular inner trims.
//!
//! The 2D triangulator is an implementation detail. Actual geometry remains
//! in canonical 3D circle edges, cylindrical surfaces and coedge-specific UV
//! loops. The result MUST pass area, ownership, orientation, manifold and
//! closed-volume certificates before a single triangle is exposed.
use crate::{
    BrepModel, BrepStructuralEditor, Curve3, FaceId, GeometryError,
    Mesh, Point3, Solid,
};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

const MIN_SEGMENTS:usize=12;
const MAX_SEGMENTS:usize=512;
const MAX_TOTAL_CIRCULAR_VERTICES:usize=16384;
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn magnitude(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn xy_area2(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn polygon_area(points:&[[f64;2]])->f64{
    let p0=points[0];
    (0..points.len()).map(|i|xy_area2(p0,points[i],
        points[(i+1)%points.len()])).sum::<f64>()*0.5
}
fn polygon_perimeter(points:&[[f64;2]])->f64{
    (0..points.len()).map(|i|{
        let a=points[i];let b=points[(i+1)%points.len()];
        (a[0]-b[0]).hypot(a[1]-b[1])
    }).sum()
}
fn inside_ring(point:[f64;2],points:&[[f64;2]])->bool{
    let mut odd=false;
    for k in 0..points.len(){
        let a=points[k];let b=points[(k+1)%points.len()];
        if (a[1]>point[1])!=(b[1]>point[1]) &&
            point[0] < a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]) {
            odd=!odd;
        }
    }
    odd
}

/// A planar face with a single outer polygon and zero or more interior
/// circular loops, projected into a stable local XZ chart.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct FaceRegionCertificate {
    pub face:FaceId,
    pub loop_count:usize,
    pub triangles:usize,
    /// 2D tessellated material area, not the analytic circular domain area.
    pub polygonal_area:f64,
}

/// Independently calculated properties of the *display* polyhedron. Neither
/// is written into BrepModel.mass, which stays mathematically analytic.
#[derive(Debug,Clone)]
pub struct CertifiedDisplayMesh {
    pub mesh:Mesh,
    pub bottom_region:FaceRegionCertificate,
    pub top_region:FaceRegionCertificate,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
    pub circle_segments:usize,
}

#[derive(Debug,Clone,Copy,Default)]
pub struct FaceRegionMesher;
impl FaceRegionMesher {
    /// The source editor is independently checked against its canonical
    /// analytic features; no caller may pass an arbitrary, unchecked polygon
    /// masquerading as a multi-loop face.
    pub fn tessellate(editor:&BrepStructuralEditor,segments:usize)
        ->Result<CertifiedDisplayMesh,GeometryError>
    {
        if !(MIN_SEGMENTS..=MAX_SEGMENTS).contains(&segments) {
            return Err(invalid("multi-loop display segments must be 12..=512"));
        }
        editor.validate()?;
        let source=editor.source_box();
        let model=editor.model();
        let bindings=editor.bindings();
        if bindings.len().saturating_mul(segments)>MAX_TOTAL_CIRCULAR_VERTICES {
            return Err(invalid("multi-loop face sampling exceeds 16,384 circular vertices"));
        }
        if source.vertices.len()!=8 || source.faces.len()!=6 ||
            source.faces[0].boundary.len()!=4 ||
            model.vertices.len()!=8+2*bindings.len() ||
            model.faces.len()!=6+bindings.len() ||
            model.faces[0].loops.len()!=1+bindings.len() ||
            model.faces[1].loops.len()!=1+bindings.len(){
            return Err(invalid("noncanonical cap ownership or multi-loop topology"));
        }
        let bounds=source.bbox;
        let height=bounds.max.y-bounds.min.y;
        let scale=magnitude(bounds.max.sub(bounds.min));
        let coord=[bounds.min,bounds.max].into_iter()
            .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
            .fold(1.0_f64,f64::max);
        let epsilon=source.tolerance.length_at(scale)
            .max(32.0*f64::EPSILON*coord);
        let mut mesh=Mesh{
            vertices:model.vertices[..8].iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        // x/z coordinates are anchored to the source bounding box to avoid
        // catastrophic cancellation when a small part is far from origin.
        let xy=|p:Point3|[p.x-bounds.min.x,p.z-bounds.min.z];
        let mut bottom_indices=Vec::<u32>::new();
        let mut top_indices=Vec::<u32>::new();
        let mut contour_xy=Vec::<[f64;2]>::new();
        let mut inner_contours=Vec::<Vec<[f64;2]>>::new();
        for &id in &source.faces[0].boundary {
            let bottom=id.0 as usize;
            let top=bottom+4;
            if bottom>=4 || top>=model.vertices.len(){
                return Err(invalid("source cap corner vertex indices not canonical"));
            }
            let lower=model.vertices[bottom].position;
            let upper=model.vertices[top].position;
            if (lower.x-upper.x).hypot(lower.z-upper.z)>epsilon ||
               (lower.y-bounds.min.y).abs()>epsilon ||
               (upper.y-bounds.max.y).abs()>epsilon {
                return Err(invalid("box cap vertices do not share vertical corners"));
            }
            bottom_indices.push(bottom as u32);
            top_indices.push(top as u32);
            contour_xy.push(xy(lower));
        }
        if polygon_area(&contour_xy)<=0.0{
            return Err(invalid("box outer cap has nonpositive XZ winding"));
        }
        let outer_area=polygon_area(&contour_xy);
        let mut hole_offsets=Vec::new();
        let mut projected_hole_area=0.0;
        let mut faceted_hole_perimeter=0.0;
        for binding in &bindings {
            let ids=binding.topology;
            let lower_edge=model.geometry.edges.get(ids.bottom_edge.0 as usize)
                .ok_or_else(||invalid("missing bottom analytic circle binding"))?;
            let upper_edge=model.geometry.edges.get(ids.top_edge.0 as usize)
                .ok_or_else(||invalid("missing top analytic circle binding"))?;
            let Curve3::Circle(lower)=model.edge_curve(ids.bottom_edge)
                .ok_or_else(||invalid("missing bottom analytic circle carrier"))?
                else{return Err(invalid("inner lower edge is not analytic circle"))};
            let Curve3::Circle(upper)=model.edge_curve(ids.top_edge)
                .ok_or_else(||invalid("missing top analytic circle carrier"))?
                else{return Err(invalid("inner upper edge is not analytic circle"))};
            let (start_a,end_a)=match lower_edge.domain {
                crate::ParameterRange::Bounded{start,end}=> (start,end),
                _=>return Err(invalid("bottom circular edge must be periodic bounded"))
            };
            let (start_b,end_b)=match upper_edge.domain {
                crate::ParameterRange::Bounded{start,end}=>(start,end),
                _=>return Err(invalid("top circular edge must be periodic bounded"))
            };
            if (end_a-start_a-TAU).abs()>source.tolerance.angular ||
               (end_b-start_b-TAU).abs()>source.tolerance.angular {
                return Err(invalid("circular B-rep edge parameter interval is not a complete turn"));
            }
            hole_offsets.push(contour_xy.len());
            let mut ring=Vec::with_capacity(segments);
            for j in 0..segments {
                let fraction=j as f64/segments as f64;
                let down=lower.evaluate(start_a+(end_a-start_a)*fraction)?;
                let up=upper.evaluate(start_b+(end_b-start_b)*fraction)?;
                if (down.x-up.x).hypot(down.z-up.z)>epsilon ||
                   (down.y-bounds.min.y).abs()>epsilon ||
                   (up.y-bounds.max.y).abs()>epsilon ||
                   ((down.x-binding.params.center_x)
                        .hypot(down.z-binding.params.center_z)-binding.params.radius).abs()
                        >epsilon {
                    return Err(invalid("top/bottom circular B-rep carriers do not align"));
                }
                let new_index=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("display mesh vertex index overflow"))?;
                mesh.vertices.push(down);
                mesh.vertices.push(up);
                bottom_indices.push(new_index);
                top_indices.push(new_index+1);
                let point=xy(down);
                ring.push(point);
                contour_xy.push(point);
            }
            let area=polygon_area(&ring);
            if !area.is_finite() || area>=-epsilon*epsilon {
                return Err(invalid("sampled inner cap curve must have clockwise XZ winding"));
            }
            projected_hole_area-=area;
            faceted_hole_perimeter+=polygon_perimeter(&ring);
            inner_contours.push(ring);
        }
        // Earcut is a triangle *proposal*; its output must match the full
        // topology-derived planar region area, stay outside all excluded
        // polygonal holes, and produce a closed orientable 3D boundary.
        let mut coords=Vec::with_capacity(contour_xy.len()*2);
        for p in &contour_xy { coords.extend_from_slice(p); }
        let proposed=earcutr::earcut(&coords,&hole_offsets,2)
            .map_err(|e|invalid(format!("planar multi-loop earcut failed: {e:?}")))?;
        if proposed.is_empty() || proposed.len()%3!=0 {
            return Err(invalid("multi-loop triangulation returned no valid triangles"));
        }
        let mut bottom_area=0.0;
        let count=proposed.len()/3;
        for triple in proposed.chunks_exact(3) {
            let mut local=[triple[0],triple[1],triple[2]];
            if local.iter().any(|&i|i>=contour_xy.len()) ||
                local[0]==local[1] || local[1]==local[2] ||
                local[2]==local[0] {
                return Err(invalid("planar triangulation returned degenerate indices"));
            }
            let mut area2=xy_area2(contour_xy[local[0]],
                contour_xy[local[1]],contour_xy[local[2]]);
            if !area2.is_finite() || area2.abs()<=source.tolerance.area_at(scale)*0.1 {
                return Err(invalid("planar triangulation has degenerate sub-tolerance face"));
            }
            if area2<0.0 {local.swap(1,2);area2=-area2;}
            let centroid=[
                (contour_xy[local[0]][0]+contour_xy[local[1]][0]+
                    contour_xy[local[2]][0])/3.0,
                (contour_xy[local[0]][1]+contour_xy[local[1]][1]+
                    contour_xy[local[2]][1])/3.0,
            ];
            if !inside_ring(centroid,&contour_xy[..4]) ||
                inner_contours.iter().any(|ring|inside_ring(centroid,ring)) {
                return Err(invalid("planar triangle covers an excluded region"));
            }
            bottom_area+=area2/2.0;
            let down=local.map(|i|bottom_indices[i]);
            let mut top=local.map(|i|top_indices[i]);
            top.swap(1,2);
            mesh.triangles.extend([down,top]);
            mesh.triangle_faces.extend([FaceId(0),FaceId(1)]);
        }
        let area_expect=outer_area-projected_hole_area;
        let area_budget=source.tolerance.area_at(scale)*count as f64+
            epsilon*scale*32.0;
        if !bottom_area.is_finite() || bottom_area<=0.0 ||
            (bottom_area-area_expect).abs()>area_budget {
            return Err(invalid(format!(
                "planar face region area is inconsistent: got {bottom_area}, expected {area_expect}")));
        }
        // Side wall triangles are original source-box facets; reuse both
        // their exact indices and source B-rep face lineage.
        for (tri,owner) in source.mesh.triangles.iter()
            .zip(&source.mesh.triangle_faces){
            if owner.0>=2 {
                mesh.triangles.push(*tri);
                mesh.triangle_faces.push(*owner);
            }
        }
        // The analytical side is sampled only for display; both bottom/top
        // caps and the cylinder wall share the very same ring vertex IDs.
        let outer_vertex_count=4usize;
        for (ordinal,binding) in bindings.iter().enumerate(){
            let offset=outer_vertex_count+ordinal*segments;
            for j in 0..segments {
                let k=(j+1)%segments;
                let b0=bottom_indices[offset+j];
                let t0=top_indices[offset+j];
                let b1=bottom_indices[offset+k];
                let t1=top_indices[offset+k];
                mesh.triangles.extend([
                    [b0,t0,t1],
                    [b0,t1,b1],
                ]);
                mesh.triangle_faces.extend([
                    binding.topology.side_face,binding.topology.side_face]);
            }
        }
        let (faceted_volume,faceted_surface_area)=
            certify_closed_facets(source,model,&bindings,&mesh,epsilon)?;
        // The cap removes the sampled *polygon* at each circular hole.
        // That analytical polygonal prism volume is an independent oracle
        // for the mesh's signed tetrahedra.
        let expected_volume=source.mass.volume-projected_hole_area*height;
        let expected_area=source.mass.surface_area-2.0*projected_hole_area+
            faceted_hole_perimeter*height;
        let volume_budget=source.mass.volume.abs()*1e-8+
            epsilon*scale*scale*16.0;
        let surface_budget=source.mass.surface_area.abs()*1e-8+
            epsilon*scale*64.0;
        if !expected_volume.is_finite() || expected_volume<=0.0 ||
           !expected_area.is_finite() || expected_area<=0.0 ||
           (faceted_volume-expected_volume).abs()>volume_budget ||
           (faceted_surface_area-expected_area).abs()>surface_budget {
            return Err(invalid(format!("facet area or signed volume disagree with independent polygonal oracle: V={faceted_volume}, expected={expected_volume}; S={faceted_surface_area}, expected={expected_area}")));
        }
        let analytic_volume=model.mass.volume;
        if faceted_volume+volume_budget<analytic_volume{
            return Err(invalid("cylindrical display mesh underestimates material volume"));
        }
        Ok(CertifiedDisplayMesh{
            mesh,
            bottom_region:FaceRegionCertificate{face:FaceId(0),
                loop_count:bindings.len()+1,triangles:count,
                polygonal_area:bottom_area},
            top_region:FaceRegionCertificate{face:FaceId(1),
                loop_count:bindings.len()+1,triangles:count,
                polygonal_area:bottom_area},
            faceted_volume,faceted_surface_area,analytic_volume,
            faceting_volume_error:faceted_volume-analytic_volume,
            circle_segments:segments,
        })
    }
}
fn certify_closed_facets(source:&Solid,model:&BrepModel,
    holes:&[crate::StructuralHoleBinding],mesh:&Mesh,eps:f64)
    ->Result<(f64,f64),GeometryError>
{
    if mesh.vertices.iter().any(|p|!p.is_finite()) ||
        mesh.triangles.len()!=mesh.triangle_faces.len() {
        return Err(invalid("invalid display mesh vertices or triangle ownership"));
    }
    let anchor=source.bbox.min;
    let mut six_volume=0.0;
    let mut area=0.0;
    let mut pair_directions=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let extent=magnitude(source.bbox.max.sub(source.bbox.min));
    for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        if owner.0 as usize>=model.faces.len() ||
            tri.iter().any(|&i|i as usize>=mesh.vertices.len()) {
            return Err(invalid("triangle references invalid face or vertex"));
        }
        let [a,b,c]=*tri;
        if a==b || b==c || a==c{return Err(invalid("mesh has repeated triangle indices"));}
        let v0=mesh.vertices[a as usize];
        let v1=mesh.vertices[b as usize];
        let v2=mesh.vertices[c as usize];
        let normal=v1.sub(v0).cross(v2.sub(v0));
        let norm=magnitude(normal);
        if !norm.is_finite() || norm<=source.tolerance.area_at(extent)*0.1 {
            return Err(invalid("zero-area display triangle"));
        }
        let centroid=v0.add(v1).add(v2).scale(1.0/3.0);
        match model.face_surface(*owner) {
            Some(crate::Surface3::Plane(plane))=>{
                let dot=normal.dot(plane.frame.normal);
                if dot<=0.0 ||
                    [v0,v1,v2].iter().any(|v|
                        (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps) {
                    return Err(invalid("display planar triangle has incorrect face plane or winding"));
                }
            }
            Some(crate::Surface3::Cylinder(_))=>{
                let Some(hole)=holes.iter().find(|h|h.topology.side_face==*owner)
                    else{return Err(invalid("unknown cylindrical display face"))};
                let inward=Point3{x:hole.params.center_x-centroid.x,
                    y:0.0,z:hole.params.center_z-centroid.z};
                if normal.dot(inward)<=0.0 ||
                    [v0,v1,v2].iter().any(|p|{
                        let r=(p.x-hole.params.center_x).hypot(p.z-hole.params.center_z);
                        (r-hole.params.radius).abs()>eps
                    }){
                    return Err(invalid("cylindrical display facet has incorrect analytic boundary or inward winding"));
                }
            }
            _=>return Err(invalid("unsupported display face surface carrier")),
        }
        for (u,v) in [(a,b),(b,c),(c,a)]{
            pair_directions.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
        six_volume+=v0.sub(anchor).dot(v1.sub(anchor).cross(v2.sub(anchor)));
        area+=norm*0.5;
    }
    if pair_directions.values().any(|uses|uses.len()!=2 || uses[0]==uses[1]){
        return Err(invalid("multi-loop display mesh is open or inconsistently oriented"));
    }
    let volume=six_volume/6.0;
    if !volume.is_finite() || volume<=0.0 || !area.is_finite(){
        return Err(invalid("closed display mesh mass is invalid"));
    }
    Ok((volume,area))
}
