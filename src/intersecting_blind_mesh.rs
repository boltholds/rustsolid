//! Closed, source-face-owned viewport mesh for the v0.27 two-arc crossed
//! cylindrical cut. This is only a *derived representation*: analytical
//! curves, UV loops, face senses and exact mass stay in the canonical B-rep.
//! The independent facet certificate checks polygon-domain areas, both
//! cylindrical patches, manifold edge uses, signed volume and surface area.
use crate::{
    BrepModel, Curve3, FaceId, GeometryError, IntersectingBlindBrep,
    Mesh, Point3, Solid, Surface3,
};
use std::collections::BTreeMap;

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(s.into())
}
fn length(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn area2(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn polygon_area(points:&[[f64;2]])->f64{
    let anchor=points[0];
    (0..points.len()).map(|i|area2(anchor,points[i],
        points[(i+1)%points.len()])).sum::<f64>()*0.5
}
fn perimeter(points:&[[f64;2]])->f64{
    (0..points.len()).map(|i|{
        let a=points[i];let b=points[(i+1)%points.len()];
        (a[0]-b[0]).hypot(a[1]-b[1])
    }).sum()
}
fn inside(point:[f64;2],polygon:&[[f64;2]])->bool{
    let mut result=false;
    for i in 0..polygon.len(){
        let a=polygon[i];let b=polygon[(i+1)%polygon.len()];
        if (a[1]>point[1])!=(b[1]>point[1]) &&
            point[0]<a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]){
            result=!result;
        }
    }
    result
}
#[derive(Debug,Clone)]
pub struct CrossedBlindDisplay {
    pub mesh:Mesh,
    pub segments_per_arc:usize,
    pub opening_triangles:usize,
    pub floor_triangles:usize,
    pub faceted_opening_area:f64,
    pub faceted_floor_area:f64,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
}
pub struct CrossedBlindMesher;
impl CrossedBlindMesher{
    pub fn tessellate(body:&IntersectingBlindBrep,segments:usize)
        ->Result<CrossedBlindDisplay,GeometryError>
    {
        if !(12..=512).contains(&segments){
            return Err(invalid("intersecting blind arcs require 12..=512 display segments per arc"));
        }
        body.validate()?;
        let b=body.source_box();
        let m=body.model();
        let region=body.region();
        let holes=body.parameters();
        if m.faces.len()!=9 || m.vertices.len()!=12 ||
           m.faces[1].loops.len()!=2 || m.faces[8].loops.len()!=1 {
            return Err(invalid("display input is not certified crossing-disk B-rep"));
        }
        let bbox=b.bbox;
        let h=region.cut_depth;
        let top=bbox.max.y;
        let floor=top-h;
        let span=length(bbox.max.sub(bbox.min));
        let absolute=[bbox.min,bbox.max].into_iter()
            .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
            .fold(1.0_f64,f64::max);
        let eps=b.tolerance.length_at(span)
            .max(32.0*f64::EPSILON*absolute);
        let mut mesh=Mesh{
            vertices:b.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        let mut outer=Vec::with_capacity(4);
        let mut upper_ids=Vec::with_capacity(4);
        for &vid in &b.faces[0].boundary{
            let id=vid.0 as usize;
            if id>=4{return Err(invalid("display box requires canonical bottom corner indices"));}
            let lower=b.vertices[id].position;
            let upper=b.vertices[id+4].position;
            if (lower.x-upper.x).hypot(lower.z-upper.z)>eps ||
                (lower.y-bbox.min.y).abs()>eps ||
                (upper.y-top).abs()>eps {
                return Err(invalid("box source cap vertices are inconsistent"));
            }
            outer.push([lower.x-bbox.min.x,lower.z-bbox.min.z]);
            upper_ids.push((id+4) as u32);
        }
        if polygon_area(&outer)<=0.0 {
            return Err(invalid("outer planar XZ cap must have positive winding"));
        }
        let ring_count=segments*2;
        let mut contour=Vec::with_capacity(ring_count);
        let mut top_ring=Vec::with_capacity(ring_count);
        let mut floor_ring=Vec::with_capacity(ring_count);
        let mut owner=Vec::with_capacity(ring_count);
        for (i,arc) in [region.arc_a,region.arc_b].into_iter().enumerate(){
            let lower_curve=m.edge_curve(crate::EdgeId(12+2*i as u32))
                .ok_or_else(||invalid("missing analytic floor circular arc"))?;
            let upper_curve=m.edge_curve(crate::EdgeId(13+2*i as u32))
                .ok_or_else(||invalid("missing analytic top circular arc"))?;
            if !matches!(lower_curve,Curve3::Circle(_)) ||
                !matches!(upper_curve,Curve3::Circle(_)) {
                return Err(invalid("crossed bore mesh can sample only exact circular carriers"));
            }
            for j in 0..segments{
                let u=arc.parameter_start+
                    arc.angular_span()*j as f64/segments as f64;
                let (low,high)=if j==0{
                    // The two circle portions share one exact crossing
                    // vertex on each cap, not two nearly equal copies.
                    let crossing=if i==0 {0}else{1};
                    (m.vertices[8+crossing].position,
                     m.vertices[10+crossing].position)
                }else{(lower_curve.evaluate(u)?,upper_curve.evaluate(u)?)};
                if !low.is_finite()||!high.is_finite() ||
                    (high.y-top).abs()>eps || (low.y-floor).abs()>eps ||
                    (low.x-high.x).hypot(low.z-high.z)>eps ||
                    ((low.x-holes[i].center_x)
                        .hypot(low.z-holes[i].center_z)-holes[i].radius).abs()>eps {
                    return Err(invalid("crossed bore displayed arc leaves source circle or height"));
                }
                let vid=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("display vertex index exhausted"))?;
                mesh.vertices.extend([low,high]);
                floor_ring.push(vid);
                top_ring.push(vid+1);
                contour.push([low.x-bbox.min.x,low.z-bbox.min.z]);
                owner.push(FaceId(6+i as u32));
            }
        }
        let signed=polygon_area(&contour);
        if !signed.is_finite() || signed>=0.0 {
            return Err(invalid("split union opening must be one simple clockwise XZ loop"));
        }
        let opening=-signed;
        let relative_tol=b.tolerance.area_at(span)*(ring_count as f64+8.0)
            +eps*span*32.0;
        if opening<=relative_tol || opening>=polygon_area(&outer)-relative_tol {
            return Err(invalid("sampled union aperture area is degenerate"));
        }
        // Triangulate the source top face with exactly one compound
        // inner-loop derived from both analytic circle arcs.
        let mut top_points=outer.clone();
        top_points.extend_from_slice(&contour);
        let mut coordinates=Vec::with_capacity(top_points.len()*2);
        for v in &top_points{coordinates.extend_from_slice(v);}
        let top_triangles=earcutr::earcut(&coordinates,&[4usize],2)
            .map_err(|e|invalid(format!("top union-loop triangulation failed: {e:?}")))?;
        if top_triangles.is_empty()||top_triangles.len()%3!=0 {
            return Err(invalid("compound top opening has no planar triangles"));
        }
        let mut opening_area=0.0;
        for tri in top_triangles.chunks_exact(3){
            let mut t=[tri[0],tri[1],tri[2]];
            if t.iter().any(|&v|v>=top_points.len())||
                t[0]==t[1]||t[1]==t[2]||t[2]==t[0] {
                return Err(invalid("invalid top planar triangle indices"));
            }
            let mut signed=area2(top_points[t[0]],top_points[t[1]],top_points[t[2]]);
            if !signed.is_finite()||signed.abs()<=b.tolerance.area_at(span)*0.1 {
                return Err(invalid("top union-loop triangle is degenerate"));
            }
            if signed>0.0{t.swap(1,2);signed=-signed;}
            let center=[
                (top_points[t[0]][0]+top_points[t[1]][0]+top_points[t[2]][0])/3.0,
                (top_points[t[0]][1]+top_points[t[1]][1]+top_points[t[2]][1])/3.0,
            ];
            if !inside(center,&outer)||inside(center,&contour){
                return Err(invalid("top triangular region violates compound hole exclusion"));
            }
            opening_area+=-signed/2.0;
            let f=|id:usize|if id<4 {upper_ids[id]}else{top_ring[id-4]};
            mesh.triangles.push(t.map(f));
            mesh.triangle_faces.push(FaceId(1));
        }
        if (opening_area-(polygon_area(&outer)-opening)).abs()>relative_tol {
            return Err(invalid("top compound annular material area mismatches independent region area"));
        }
        // The floor is an independent upward-facing *single disk union*
        // with one outer loop made of both exposed circular arcs.
        let mut floor_coords=Vec::with_capacity(contour.len()*2);
        for point in &contour{floor_coords.extend_from_slice(point);}
        let floor_tris=earcutr::earcut(&floor_coords,&[],2)
            .map_err(|e|invalid(format!("crossed floor triangulation failed: {e:?}")))?;
        if floor_tris.is_empty()||floor_tris.len()%3!=0 {
            return Err(invalid("crossed opening floor has no planar triangles"));
        }
        let mut floor_area=0.0;
        for tri in floor_tris.chunks_exact(3){
            let mut t=[tri[0],tri[1],tri[2]];
            if t.iter().any(|&x|x>=contour.len()) ||
                t[0]==t[1]||t[1]==t[2]||t[2]==t[0] {
                return Err(invalid("floor disk union triangulation indices invalid"));
            }
            let mut signed=area2(contour[t[0]],contour[t[1]],contour[t[2]]);
            if !signed.is_finite()||signed.abs()<=b.tolerance.area_at(span)*0.1 {
                return Err(invalid("crossed floor has degenerate triangular patch"));
            }
            if signed>0.0{t.swap(1,2);signed=-signed;}
            let center=[
                (contour[t[0]][0]+contour[t[1]][0]+contour[t[2]][0])/3.0,
                (contour[t[0]][1]+contour[t[1]][1]+contour[t[2]][1])/3.0,
            ];
            if !inside(center,&contour){
                return Err(invalid("floor facet escaped union of source circles"));
            }
            floor_area+=-signed/2.0;
            mesh.triangles.push(t.map(|id|floor_ring[id]));
            mesh.triangle_faces.push(FaceId(8));
        }
        if (floor_area-opening).abs()>relative_tol {
            return Err(invalid("compound floor area diverges from top circle union"));
        }
        // Reuse the existing original bottom and four side planes.
        for (triangle,face) in b.mesh.triangles.iter()
            .zip(&b.mesh.triangle_faces){
            if face.0!=1 {
                mesh.triangles.push(*triangle);
                mesh.triangle_faces.push(*face);
            }
        }
        // The same physical ring vertices close both cylinders at their
        // two shared crossing-line edges: no duplicate seam coordinates.
        for i in 0..ring_count{
            let next=(i+1)%ring_count;
            let current_owner=owner[i];
            // Both split cylinder patches contain exactly `segments`
            // adjacent quads. The final quad of one patch ends on the
            // *same exact ring vertex* as the first of the other.
            mesh.triangles.extend([
                [floor_ring[i],top_ring[i],top_ring[next]],
                [floor_ring[i],top_ring[next],floor_ring[next]],
            ]);
            mesh.triangle_faces.extend([current_owner,current_owner]);
        }
        let (faceted_volume,faceted_surface_area)=
            certify_mesh(b,m,holes,&mesh,eps)?;
        let expected_volume=b.mass.volume-opening*h;
        let expected_surface=b.mass.surface_area+perimeter(&contour)*h;
        let volume_slack=b.mass.volume.abs()*1.0e-8+eps*span*span*16.0;
        let area_slack=b.mass.surface_area.abs()*1.0e-8+eps*span*64.0;
        if (faceted_volume-expected_volume).abs()>volume_slack ||
            (faceted_surface_area-expected_surface).abs()>area_slack ||
            faceted_volume+volume_slack<m.mass.volume {
            return Err(invalid(format!(
                "crossed blind display volume/area fails independent polygon oracle: volume={faceted_volume}, expected={expected_volume}, surface={faceted_surface_area}, expected={expected_surface}")));
        }
        Ok(CrossedBlindDisplay{
            mesh,segments_per_arc:segments,
            opening_triangles:top_triangles.len()/3,
            floor_triangles:floor_tris.len()/3,
            faceted_opening_area:opening_area,
            faceted_floor_area:floor_area,
            faceted_volume,faceted_surface_area,
            analytic_volume:m.mass.volume,
            faceting_volume_error:faceted_volume-m.mass.volume,
        })
    }
}
fn certify_mesh(box_:&Solid,model:&BrepModel,
    holes:[crate::BlindHoleParameters;2],mesh:&Mesh,eps:f64)
    ->Result<(f64,f64),GeometryError>
{
    if mesh.triangles.len()!=mesh.triangle_faces.len() ||
        mesh.vertices.iter().any(|v|!v.is_finite()){
        return Err(invalid("crossed blind display has invalid face ownership or coordinates"));
    }
    let anchor=box_.bbox.min;
    let mut use_map=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut six_volume=0.0;
    let mut surface_area=0.0;
    for (tri,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        let [a,b,c]=*tri;
        if a==b||b==c||c==a||
            tri.iter().any(|&id|id as usize>=mesh.vertices.len()) ||
            face.0 as usize>=model.faces.len() {
            return Err(invalid("triangle references unknown vertex or analytical face"));
        }
        let [v0,v1,v2]=tri.map(|id|mesh.vertices[id as usize]);
        let normal=v1.sub(v0).cross(v2.sub(v0));
        let length=length(normal);
        if !length.is_finite()||length<=0.0 {
            return Err(invalid("degenerate crossed circular display triangle"));
        }
        let center=v0.add(v1).add(v2).scale(1.0/3.0);
        match face.0{
            6|7=>{
                let h=holes[(face.0-6) as usize];
                let inward=Point3{x:h.center_x-center.x,
                    y:0.0,z:h.center_z-center.z};
                if model.face_sense(*face)!=Some(crate::FaceSense::Reversed) ||
                    normal.dot(inward)<=0.0 ||
                    [v0,v1,v2].iter().any(|v|
                        ((v.x-h.center_x).hypot(v.z-h.center_z)-h.radius).abs()>eps) {
                    return Err(invalid("cylindrical crossed bore facet has invalid carrier or inward face normal"));
                }
            }
            8=>{
                if model.face_sense(*face)!=Some(crate::FaceSense::Reversed) ||
                   normal.y<=0.0 ||
                   [v0,v1,v2].iter().any(|v|
                       (v.y-(box_.bbox.max.y-(box_.bbox.max.y-center.y))).abs()>eps) {
                    // Check plane independently against canonical floor
                    // support below; the comparison above guards deviation
                    // across three vertices without trusting source normals.
                    return Err(invalid("compound blind floor has non-upward normal or warped triangle"));
                }
                let Some(Surface3::Plane(floor))=model.face_surface(*face)else{
                    return Err(invalid("missing true analytical floor support"));
                };
                if [v0,v1,v2].iter().any(|v|
                    (v.sub(floor.frame.origin).dot(floor.frame.normal)).abs()>eps){
                    return Err(invalid("compound floor triangle is off its analytical plane"));
                }
            }
            _=>{
                let Some(Surface3::Plane(plane))=model.face_surface(*face)else {
                    return Err(invalid("non-planar old box owner"));
                };
                if normal.dot(plane.frame.normal)<=0.0 ||
                   [v0,v1,v2].iter().any(|v|
                       (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps){
                    return Err(invalid("original box facet has invalid support normal"));
                }
            }
        }
        for (from,to) in [(a,b),(b,c),(c,a)]{
            use_map.entry((from.min(to),from.max(to)))
                .or_default().push(from<to);
        }
        six_volume+=v0.sub(anchor).dot(v1.sub(anchor).cross(v2.sub(anchor)));
        surface_area+=length*0.5;
    }
    if use_map.values().any(|incidences|
        incidences.len()!=2 || incidences[0]==incidences[1]){
        return Err(invalid("crossed blind display mesh has open/nonorientable facet edges"));
    }
    let volume=six_volume/6.0;
    if !volume.is_finite()||volume<=0.0||
        !surface_area.is_finite()||surface_area<=0.0{
        return Err(invalid("crossed blind display facet integrals are invalid"));
    }
    Ok((volume,surface_area))
}
