//! Independent watertight DISPLAY meshing of an analytical heterogeneous
//! B-rep containing through-bores and disjoint blind regions, where a blind
//! region may be a single circle or two crossing analytical circle arcs.

//! The analytical B-rep and mass are never defined by, or changed by,
//! the resultant polygonal mesh.
use crate::{
    BlindRegionKind,MixedThroughBlindBrep,Curve3,FaceId,GeometryError,
    Mesh,Point3,Solid,Surface3,
};
use std::collections::BTreeMap;
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn mag(v:Point3)->f64{v.x.hypot(v.y).hypot(v.z)}
fn area2(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn signed_area(points:&[[f64;2]])->f64{
    let p=points[0];
    (0..points.len()).map(|i|area2(p,points[i],
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
#[derive(Debug,Clone,PartialEq)]
pub struct MixedHoleFloorRegion{
    pub component_index:usize,
    pub floor_face:FaceId,
    pub loop_vertices:usize,
    pub floor_triangles:usize,
    pub polygonal_area:f64,
    pub polygonal_perimeter:f64,
}
#[derive(Debug,Clone)]
pub struct MixedHoleDisplay{
    pub mesh:Mesh,
    pub segments_per_circle:usize,
    pub top_face_triangles:usize,
    pub bottom_face_triangles:usize,
    pub top_face_loop_count:usize,
    pub bottom_face_loop_count:usize,
    pub floor_regions:Vec<MixedHoleFloorRegion>,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
}
#[derive(Debug,Clone)]
struct Ring{
    points:Vec<[f64;2]>,
    floor_vertices:Vec<u32>,
    top_vertices:Vec<u32>,
    source_face_per_segment:Vec<FaceId>,
    floor_face:Option<FaceId>,
    depth:f64,
    through:bool,
}
pub struct MixedHoleMesher;
impl MixedHoleMesher {
    pub fn tessellate(body:&MixedThroughBlindBrep,segments:usize)
        ->Result<MixedHoleDisplay,GeometryError>
    {
        if !(12..=512).contains(&segments){
            return Err(invalid("compound blind display requires 12..=512 segments per source circle"));
        }
        body.validate()?;
        let source=body.source_box();
        let m=body.model();
        let regions=body.components();
        let holes=body.blind_holes();
        let through=body.through_holes();
        let total_circles=through.len()+
            regions.iter().map(|c|c.circle_count()).sum::<usize>();
        if total_circles.saturating_mul(segments)>8192{
            return Err(invalid("compound blind display work budget exceeds 8192 curve samples"));
        }
        if m.faces[1].loops.len()!=regions.len()+through.len()+1 ||
           m.faces[0].loops.len()!=through.len()+1 ||
           source.vertices.len()!=8 {
            return Err(invalid("compound top face loops do not match certified source components"));
        }
        let bbox=source.bbox;
        let span=mag(bbox.max.sub(bbox.min));
        let absolute=[bbox.min,bbox.max].into_iter()
            .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
            .fold(1.0_f64,f64::max);
        let eps=source.tolerance.length_at(span)
            .max(32.0*f64::EPSILON*absolute);
        let mut mesh=Mesh{vertices:source.vertices.iter()
            .map(|v|v.position).collect(),triangles:Vec::new(),
            triangle_faces:Vec::new()};
        let mut outer=Vec::with_capacity(4);
        let mut outer_top=Vec::with_capacity(4);
        for &id in &source.faces[0].boundary{
            let v=id.0 as usize;
            if v>=4{return Err(invalid("source box is not a canonical four-corner bottom cap"));}
            let p=source.vertices[v].position;
            let q=source.vertices[v+4].position;
            if (p.x-q.x).hypot(p.z-q.z)>eps ||
                (p.y-bbox.min.y).abs()>eps ||
                (q.y-bbox.max.y).abs()>eps {
                return Err(invalid("source vertical box corners disagree"));
            }
            outer.push([p.x-bbox.min.x,p.z-bbox.min.z]);
            outer_top.push((v+4) as u32);
        }
        if signed_area(&outer)<=0.0{
            return Err(invalid("outer box cap is not counterclockwise in local XZ coordinates"));
        }
        let mut rings=Vec::with_capacity(regions.len()+through.len());
        let mut contours=outer.clone();
        let mut top_indices=outer_top;
        let mut hole_starts=Vec::<usize>::new();
        let mut polygonal_removed=0.0;
        // The lower planar box cap has holes ONLY for the full-height
        // bores; blind floors never punch through to this plane.
        let mut bottom_contours=outer.clone();
        let mut bottom_indices=source.faces[0].boundary.iter()
            .map(|id|id.0).collect::<Vec<u32>>();
        let mut bottom_hole_starts=Vec::<usize>::new();
        let mut bottom_removed=0.0;
        let full_height=bbox.max.y-bbox.min.y;
        for (ti,cut) in through.iter().enumerate() {
            let lower_edge=crate::EdgeId(12+3*ti as u32);
            let upper_edge=crate::EdgeId(13+3*ti as u32);
            let side=FaceId(6+ti as u32);
            let lower=m.edge_curve(lower_edge)
                .ok_or_else(||invalid("missing analytical through-hole lower circle"))?;
            let upper=m.edge_curve(upper_edge)
                .ok_or_else(||invalid("missing analytical through-hole upper circle"))?;
            if !matches!(lower,Curve3::Circle(_))||
               !matches!(upper,Curve3::Circle(_)) ||
               m.face_sense(side)!=Some(crate::FaceSense::Reversed) {
                return Err(invalid("through hole display requires actual inward analytic cylinder and circles"));
            }
            let domain=m.geometry.edges[lower_edge.0 as usize].domain;
            if m.geometry.edges[upper_edge.0 as usize].domain!=domain{
                return Err(invalid("through circle caps have different parameter intervals"));
            }
            let crate::ParameterRange::Bounded{start,end}=domain else{
                return Err(invalid("through-bore 3D circle is unbounded"));
            };
            hole_starts.push(contours.len());
            bottom_hole_starts.push(bottom_contours.len());
            let mut ring=Ring{
                points:Vec::with_capacity(segments),
                floor_vertices:Vec::with_capacity(segments),
                top_vertices:Vec::with_capacity(segments),
                source_face_per_segment:Vec::with_capacity(segments),
                floor_face:None,depth:full_height,through:true,
            };
            for j in 0..segments {
                let u=start+(end-start)*j as f64/segments as f64;
                let down=lower.evaluate(u)?;let up=upper.evaluate(u)?;
                if !down.is_finite()||!up.is_finite() ||
                   (down.y-bbox.min.y).abs()>eps ||
                   (up.y-bbox.max.y).abs()>eps ||
                   (down.x-up.x).hypot(down.z-up.z)>eps ||
                   ((down.x-cut.center_x).hypot(down.z-cut.center_z)-cut.radius).abs()>eps {
                    return Err(invalid("through-circle display sampling leaves its physical source"));
                }
                let index=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("mixed display vertex budget exhausted"))?;
                mesh.vertices.extend([down,up]);
                ring.floor_vertices.push(index);
                ring.top_vertices.push(index+1);
                ring.source_face_per_segment.push(side);
                let point=[down.x-bbox.min.x,down.z-bbox.min.z];
                ring.points.push(point);
                contours.push(point);top_indices.push(index+1);
                bottom_contours.push(point);bottom_indices.push(index);
            }
            let a=signed_area(&ring.points);
            if !a.is_finite() || a>=0.0 {
                return Err(invalid("sampled through circular inner contour is not clockwise"));
            }
            polygonal_removed-=a;
            bottom_removed-=a;
            rings.push(ring);
        }
        for (ci,c) in regions.iter().enumerate(){
            hole_starts.push(contours.len());
            let mut ring=Ring{
                points:Vec::with_capacity(c.circle_count()*segments),
                floor_vertices:Vec::with_capacity(c.circle_count()*segments),
                top_vertices:Vec::with_capacity(c.circle_count()*segments),
                source_face_per_segment:Vec::with_capacity(c.circle_count()*segments),
                floor_face:Some(c.floor_face),depth:c.depth,through:false,
            };
            for local_arc in 0..c.circle_count(){
                let low_edge=crate::EdgeId(c.first_edge.0+2*local_arc as u32);
                let high_edge=crate::EdgeId(c.first_edge.0+2*local_arc as u32+1);
                let low=m.edge_curve(low_edge)
                    .ok_or_else(||invalid("missing source floor circle carrier"))?;
                let high=m.edge_curve(high_edge)
                    .ok_or_else(||invalid("missing source top circle carrier"))?;
                if !matches!(low,Curve3::Circle(_)) || !matches!(high,Curve3::Circle(_)){
                    return Err(invalid("source aperture edge is not an analytic circular arc"));
                }
                let domain=m.geometry.edges[low_edge.0 as usize].domain;
                let crate::ParameterRange::Bounded{start,end}=domain else{
                    return Err(invalid("source circle must have bounded parameter interval"));
                };
                if m.geometry.edges[high_edge.0 as usize].domain!=domain{
                    return Err(invalid("top/floor circle arcs have inconsistent parameter extents"));
                }
                for j in 0..segments{
                    let u=start+(end-start)*j as f64/segments as f64;
                    let (lower,upper)=if c.kind==BlindRegionKind::TwoCrossingCircles && j==0{
                        let bottom_id=m.edges[low_edge.0 as usize].start.0 as usize;
                        let top_id=m.edges[high_edge.0 as usize].start.0 as usize;
                        (m.vertices[bottom_id].position,m.vertices[top_id].position)
                    }else{(low.evaluate(u)?,high.evaluate(u)?)};
                    if !lower.is_finite()||!upper.is_finite() ||
                        (lower.y-(bbox.max.y-c.depth)).abs()>eps ||
                        (upper.y-bbox.max.y).abs()>eps ||
                        (lower.x-upper.x).hypot(lower.z-upper.z)>eps {
                        return Err(invalid("analytic floor and top arc carriers are not vertical"));
                    }
                    let original=holes[c.source_indices[local_arc]];
                    if ((lower.x-original.center_x).hypot(
                        lower.z-original.center_z)-original.radius).abs()>eps{
                        return Err(invalid("display arc leaves its original source circle"));
                    }
                    let new_id=u32::try_from(mesh.vertices.len())
                        .map_err(|_|invalid("compound mesh vertex budget exhausted"))?;
                    mesh.vertices.extend([lower,upper]);
                    ring.floor_vertices.push(new_id);
                    ring.top_vertices.push(new_id+1);
                    ring.source_face_per_segment.push(c.cylinder_faces[local_arc]);
                    let point=[lower.x-bbox.min.x,lower.z-bbox.min.z];
                    ring.points.push(point);
                    contours.push(point);
                    top_indices.push(new_id+1);
                }
            }
            let signed=signed_area(&ring.points);
            if !signed.is_finite()||signed>=0.0{
                return Err(invalid("compound source region is not a clockwise aperture loop"));
            }
            polygonal_removed-=signed;
            rings.push(ring);
            let _=ci;
        }
        let mut coords=Vec::with_capacity(contours.len()*2);
        for point in &contours{coords.extend_from_slice(point);}
        let top_triangles=earcutr::earcut(&coords,&hole_starts,2)
            .map_err(|e|invalid(format!("compound top planar triangulation failed: {e:?}")))?;
        if top_triangles.is_empty() || top_triangles.len()%3!=0{
            return Err(invalid("compound top face returned invalid triangles"));
        }
        let mut cap_area=0.0;
        for t in top_triangles.chunks_exact(3){
            let mut ids=[t[0],t[1],t[2]];
            if ids.iter().any(|&i|i>=contours.len()) ||
                ids[0]==ids[1]||ids[1]==ids[2]||ids[2]==ids[0]{
                return Err(invalid("compound top triangle has invalid index"));
            }
            let mut a=area2(contours[ids[0]],contours[ids[1]],contours[ids[2]]);
            if !a.is_finite()||a.abs()<=source.tolerance.area_at(span)*0.1{
                return Err(invalid("compound cap has degenerate projected triangle"));
            }
            if a>0.0{ids.swap(1,2);a=-a;}
            let center=[
                (contours[ids[0]][0]+contours[ids[1]][0]+contours[ids[2]][0])/3.0,
                (contours[ids[0]][1]+contours[ids[1]][1]+contours[ids[2]][1])/3.0
            ];
            if !inside(center,&outer)||
                rings.iter().any(|ring|inside(center,&ring.points)){
                return Err(invalid("compound cap triangle violates an excluded analytic opening"));
            }
            cap_area+=-a/2.0;
            mesh.triangles.push(ids.map(|i|top_indices[i]));
            mesh.triangle_faces.push(FaceId(1));
        }
        let area_slack=source.tolerance.area_at(span)*
            (top_triangles.len()+contours.len()) as f64+eps*span*32.0;
        if (cap_area-(signed_area(&outer)-polygonal_removed)).abs()>area_slack{
            return Err(invalid("compound cap projected area disagrees with polygon-with-holes domain"));
        }
        // The bottom cap needs a separate annular triangulation with only
        // the through holes; use the same physically shared lower ring IDs
        // as the cylinder wall facets. No blind circle appears here.
        let mut bottom_coords=Vec::with_capacity(bottom_contours.len()*2);
        for point in &bottom_contours{bottom_coords.extend_from_slice(point);}
        let bottom_tris=earcutr::earcut(&bottom_coords,&bottom_hole_starts,2)
            .map_err(|e|invalid(format!("mixed bottom planar triangulation failed: {e:?}")))?;
        if bottom_tris.is_empty()||bottom_tris.len()%3!=0 {
            return Err(invalid("mixed bottom cap has no valid triangles"));
        }
        let mut lower_cap_area=0.0;
        for t in bottom_tris.chunks_exact(3) {
            let mut ids=[t[0],t[1],t[2]];
            if ids.iter().any(|&j|j>=bottom_contours.len()) ||
                ids[0]==ids[1]||ids[1]==ids[2]||ids[2]==ids[0] {
                return Err(invalid("mixed bottom cap has invalid vertex IDs"));
            }
            let mut a=area2(bottom_contours[ids[0]],bottom_contours[ids[1]],
                bottom_contours[ids[2]]);
            if !a.is_finite()||a.abs()<=source.tolerance.area_at(span)*0.1 {
                return Err(invalid("mixed bottom cap has sub-tolerance triangle"));
            }
            if a<0.0{ids.swap(1,2);a=-a;}
            let center=[
                (bottom_contours[ids[0]][0]+bottom_contours[ids[1]][0]+bottom_contours[ids[2]][0])/3.0,
                (bottom_contours[ids[0]][1]+bottom_contours[ids[1]][1]+bottom_contours[ids[2]][1])/3.0,
            ];
            if !inside(center,&outer) ||
                rings.iter().filter(|r|r.through).any(|r|inside(center,&r.points)) {
                return Err(invalid("mixed bottom planar triangle covers a through bore void"));
            }
            lower_cap_area+=a/2.0;
            mesh.triangles.push(ids.map(|j|bottom_indices[j]));
            mesh.triangle_faces.push(FaceId(0));
        }
        let area_slack=source.tolerance.area_at(span)*
            (top_triangles.len()+bottom_tris.len()+contours.len()) as f64+
            eps*span*32.0;
        if (lower_cap_area-(signed_area(&outer)-bottom_removed)).abs()>area_slack {
            return Err(invalid("mixed bottom cap area diverges from through-opening region"));
        }
        let mut floor_regions=Vec::with_capacity(regions.len());
        let mut expected_surface=source.mass.surface_area;
        let mut expected_volume=source.mass.volume;
        for ring in rings.iter().filter(|r|r.through) {
            let polygonal_area=-signed_area(&ring.points);
            let perimeter=perimeter(&ring.points);
            expected_volume-=polygonal_area*full_height;
            expected_surface+=perimeter*full_height-2.0*polygonal_area;
        }
        for (i,ring) in rings.iter().enumerate(){
            if ring.through {continue;}
            let floor_face=ring.floor_face
                .ok_or_else(||invalid("blind floor is missing a source face"))?;
            let mut coords=Vec::with_capacity(ring.points.len()*2);
            for point in &ring.points{coords.extend_from_slice(point);}
            let triangles=earcutr::earcut(&coords,&[],2)
                .map_err(|e|invalid(format!("component floor triangulation failed: {e:?}")))?;
            if triangles.is_empty()||triangles.len()%3!=0{
                return Err(invalid("source region floor has no valid planar triangles"));
            }
            let mut floor_area=0.0;
            for t in triangles.chunks_exact(3){
                let mut ids=[t[0],t[1],t[2]];
                if ids.iter().any(|&j|j>=ring.points.len())||
                    ids[0]==ids[1]||ids[1]==ids[2]||ids[0]==ids[2]{
                    return Err(invalid("invalid component floor triangle index"));
                }
                let mut a=area2(ring.points[ids[0]],ring.points[ids[1]],ring.points[ids[2]]);
                if !a.is_finite()||a.abs()<=source.tolerance.area_at(span)*0.1{
                    return Err(invalid("degenerate component floor triangle"));
                }
                if a>0.0{ids.swap(1,2);a=-a;}
                floor_area+=-a/2.0;
                mesh.triangles.push(ids.map(|j|ring.floor_vertices[j]));
                mesh.triangle_faces.push(floor_face);
            }
            let expected=-signed_area(&ring.points);
            if (floor_area-expected).abs()>area_slack{
                return Err(invalid("component floor area diverges from source aperture domain"));
            }
            expected_volume-=floor_area*ring.depth;
            let p=perimeter(&ring.points);
            expected_surface+=p*ring.depth;
            floor_regions.push(MixedHoleFloorRegion{
                component_index:i-through.len(),floor_face,
                loop_vertices:ring.points.len(),
                floor_triangles:triangles.len()/3,
                polygonal_area:floor_area,polygonal_perimeter:p,
            });
        }
        // Preserve unaffected original bottom and four box walls.
        for (tri,face) in source.mesh.triangles.iter()
            .zip(&source.mesh.triangle_faces){
            if face.0>=2 {
                mesh.triangles.push(*tri);
                mesh.triangle_faces.push(*face);
            }
        }
        for ring in &rings {
            for j in 0..ring.points.len(){
                let next=(j+1)%ring.points.len();
                let f0=ring.floor_vertices[j];
                let t0=ring.top_vertices[j];
                let f1=ring.floor_vertices[next];
                let t1=ring.top_vertices[next];
                mesh.triangles.extend([[f0,t0,t1],[f0,t1,f1]]);
                mesh.triangle_faces.extend([
                    ring.source_face_per_segment[j],
                    ring.source_face_per_segment[j],
                ]);
            }
        }
        let (faceted_volume,faceted_surface_area)=
            validate_display(source,m,through,holes,body.components(),&mesh,eps)?;
        let volume_slack=source.mass.volume.abs()*1e-8+eps*span*span*16.0;
        let surface_slack=source.mass.surface_area.abs()*1e-8+eps*span*64.0;
        if !expected_volume.is_finite()||!expected_surface.is_finite()||
            (faceted_volume-expected_volume).abs()>volume_slack ||
            (faceted_surface_area-expected_surface).abs()>surface_slack ||
            faceted_volume+volume_slack<m.mass.volume {
            return Err(invalid("compound display facets violate independently integrated polygonal volume/area"));
        }
        Ok(MixedHoleDisplay{
            mesh,segments_per_circle:segments,
            top_face_triangles:top_triangles.len()/3,
            bottom_face_triangles:bottom_tris.len()/3,
            top_face_loop_count:regions.len()+through.len()+1,
            bottom_face_loop_count:through.len()+1,
            floor_regions,faceted_volume,faceted_surface_area,
            analytic_volume:m.mass.volume,
            faceting_volume_error:faceted_volume-m.mass.volume
        })
    }
}
fn validate_display(source:&Solid,model:&crate::BrepModel,
    through:&[crate::HoleParameters],
    holes:&[crate::BlindHoleParameters],
    components:&[crate::BlindRegionComponent],mesh:&Mesh,eps:f64)
    ->Result<(f64,f64),GeometryError>
{
    if mesh.vertices.iter().any(|p|!p.is_finite())||
        mesh.triangles.len()!=mesh.triangle_faces.len(){
        return Err(invalid("compound mesh missing face provenance or finite vertices"));
    }
    let anchor=source.bbox.min;
    let mut edge_use=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut six_volume=0.0;
    let mut area=0.0;
    for (tri,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        let [a,b,c]=*tri;
        if a==b||b==c||a==c||
            tri.iter().any(|&id|id as usize>=mesh.vertices.len())||
            face.0 as usize>=model.faces.len() {
            return Err(invalid("compound mesh triangle has invalid vertex or source face"));
        }
        let [v0,v1,v2]=tri.map(|id|mesh.vertices[id as usize]);
        let n=v1.sub(v0).cross(v2.sub(v0));
        let nlen=mag(n);
        if !nlen.is_finite()||nlen<=0.0{
            return Err(invalid("compound mesh contains zero-area triangle"));
        }
        let centroid=v0.add(v1).add(v2).scale(1.0/3.0);
        if face.0>=6 && (face.0-6) < through.len() as u32 {
            let src=through[(face.0-6) as usize];
            let inward=Point3{x:src.center_x-centroid.x,
                y:0.0,z:src.center_z-centroid.z};
            if n.dot(inward)<=0.0 ||
                model.face_sense(*face)!=Some(crate::FaceSense::Reversed) ||
                [v0,v1,v2].iter().any(|v|
                    ((v.x-src.center_x).hypot(v.z-src.center_z)-src.radius).abs()>eps) {
                return Err(invalid("mixed through-wall display facet is outside the cylinder or inverted"));
            }
        }else if face.0>=6 {
            let component=components.iter().find(|comp|
                comp.floor_face==*face||comp.cylinder_faces.contains(face))
                .ok_or_else(||invalid("compound triangle belongs to unknown cut component"))?;
            if component.floor_face==*face {
                if n.y<=0.0 || model.face_sense(*face)!=Some(crate::FaceSense::Reversed){
                    return Err(invalid("compound blind floor faces away from removed material"));
                }
                let Some(Surface3::Plane(plane))=model.face_surface(*face)else{
                    return Err(invalid("missing analytical floor support"));
                };
                if [v0,v1,v2].iter().any(|v|
                    (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps){
                    return Err(invalid("compound floor triangle leaves its supporting plane"));
                }
            }else{
                let local=component.cylinder_faces.iter().position(|f|f==face)
                    .ok_or_else(||invalid("unknown component cylindrical patch"))?;
                let src=holes[component.source_indices[local]];
                let inward=Point3{x:src.center_x-centroid.x,
                    y:0.0,z:src.center_z-centroid.z};
                if n.dot(inward)<=0.0||
                    model.face_sense(*face)!=Some(crate::FaceSense::Reversed)||
                    [v0,v1,v2].iter().any(|v|
                        ((v.x-src.center_x).hypot(v.z-src.center_z)-src.radius).abs()>eps){
                    return Err(invalid("compound cylindrical facet is outside source carrier or inverted"));
                }
            }
        }else{
            let Some(Surface3::Plane(plane))=model.face_surface(*face)else{
                return Err(invalid("unmodified box face is no longer planar"));
            };
            if n.dot(plane.frame.normal)<=0.0 ||
                [v0,v1,v2].iter().any(|v|
                    (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps){
                return Err(invalid("original box facet is no longer outward on its plane"));
            }
        }
        for (from,to) in [(a,b),(b,c),(c,a)]{
            edge_use.entry((from.min(to),from.max(to)))
                .or_default().push(from<to);
        }
        six_volume+=v0.sub(anchor).dot(v1.sub(anchor).cross(v2.sub(anchor)));
        area+=nlen*0.5;
    }
    if edge_use.values().any(|v|v.len()!=2 || v[0]==v[1]){
        return Err(invalid("compound display mesh has open or inconsistent edge incidence"));
    }
    let volume=six_volume/6.0;
    if !volume.is_finite()||volume<=0.0||!area.is_finite()||area<=0.0{
        return Err(invalid("compound display mesh has invalid signed volume or surface area"));
    }
    Ok((volume,area))
}
