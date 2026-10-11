//! Derived, independently certified display mesh for one *mixed* cap
//! arrangement: top openings include every cut, bottom openings only cuts
//! that penetrate the entire box; each blind cut adds a real circular floor.
//! These triangles never define the canonical B-rep or its exact mass.
use crate::{
    AnalyticCut, BrepCompositeEditor, BrepModel, CompositeCutBinding,
    Curve3, FaceId, FaceSense, GeometryError, Mesh, Point3, Surface3,
};
use std::collections::BTreeMap;
use std::f64::consts::{PI,TAU};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn len(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn twice(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn signed_area(ring:&[[f64;2]])->f64 {
    let p=ring[0];
    (0..ring.len()).map(|i|twice(p,ring[i],ring[(i+1)%ring.len()]))
        .sum::<f64>()*0.5
}
fn ring_perimeter(ring:&[[f64;2]])->f64{
    (0..ring.len()).map(|i|{
        let a=ring[i];let b=ring[(i+1)%ring.len()];
        (a[0]-b[0]).hypot(a[1]-b[1])
    }).sum()
}
fn in_polygon(q:[f64;2],ring:&[[f64;2]])->bool{
    let mut inside=false;
    for i in 0..ring.len(){
        let a=ring[i];let b=ring[(i+1)%ring.len()];
        if (a[1]>q[1])!=(b[1]>q[1]) &&
           q[0] < a[0]+(q[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]){
            inside=!inside;
        }
    }
    inside
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CompositeCapCertificate {
    pub face:FaceId,
    pub loops:usize,
    pub triangles:usize,
    pub polygonal_area:f64,
}
#[derive(Debug,Clone)]
pub struct CompositeDisplay {
    pub mesh:Mesh,
    pub bottom:CompositeCapCertificate,
    pub top:CompositeCapCertificate,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
    pub segments:usize,
}
struct CircleRing{
    binding:CompositeCutBinding,
    coords:Vec<[f64;2]>,
    lower:Vec<u32>,
    upper:Vec<u32>,
    area:f64,
    perimeter:f64,
}
fn cap(mesh:&mut Mesh,face:FaceId,outer:&[[f64;2];4],
       outer_indices:&[u32;4],rings:&[&CircleRing],top:bool,
       tolerance:f64,area_tolerance:f64)
    ->Result<CompositeCapCertificate,GeometryError>
{
    let mut points=outer.to_vec();
    let mut indices=outer_indices.to_vec();
    let mut ring_start=Vec::new();
    let mut total_holes=0.0;
    for r in rings{
        ring_start.push(points.len());
        points.extend_from_slice(&r.coords);
        indices.extend_from_slice(if top{&r.upper}else{&r.lower});
        total_holes+=r.area;
    }
    let mut packed=Vec::with_capacity(2*points.len());
    for p in &points {packed.extend_from_slice(p);}
    let triangles=earcutr::earcut(&packed,&ring_start,2)
        .map_err(|e|invalid(format!("composite cap triangulation failed: {e:?}")))?;
    if triangles.len()%3!=0 || triangles.is_empty(){
        return Err(invalid("composite cap triangulator returned malformed triangles"));
    }
    let mut area=0.0;
    for triangle in triangles.chunks_exact(3){
        let mut p=[triangle[0],triangle[1],triangle[2]];
        if p.iter().any(|&i|i>=points.len()) ||
            p[0]==p[1] ||p[1]==p[2] ||p[2]==p[0] {
            return Err(invalid("composite cap has invalid or duplicate vertex index"));
        }
        let mut area2=twice(points[p[0]],points[p[1]],points[p[2]]);
        if !area2.is_finite() || area2.abs()<=area_tolerance*0.1 {
            return Err(invalid("composite cap has degenerate projected triangle"));
        }
        if area2<0.0{p.swap(1,2);area2=-area2;}
        let q=[(points[p[0]][0]+points[p[1]][0]+points[p[2]][0])/3.0,
            (points[p[0]][1]+points[p[1]][1]+points[p[2]][1])/3.0];
        if !in_polygon(q,outer)||
            rings.iter().any(|r|in_polygon(q,&r.coords)){
            return Err(invalid("composite cap facet crosses removed analytic opening"));
        }
        let mut vertex=p.map(|i|indices[i]);
        if top {vertex.swap(1,2);}
        mesh.triangles.push(vertex);
        mesh.triangle_faces.push(face);
        area+=area2/2.0;
    }
    let expected=signed_area(outer)-total_holes;
    if !area.is_finite() || area<=0.0 ||
        (area-expected).abs()>tolerance {
        return Err(invalid(format!(
            "composite cap area {area} disagrees with source contour area {expected}"
        )));
    }
    Ok(CompositeCapCertificate{face,loops:1+rings.len(),
        triangles:triangles.len()/3,polygonal_area:area})
}
#[derive(Debug,Clone,Copy,Default)]
pub struct CompositeRegionMesher;
impl CompositeRegionMesher{
    pub fn tessellate(editor:&BrepCompositeEditor,segments:usize)
        ->Result<CompositeDisplay,GeometryError>
    {
        if !(12..=512).contains(&segments){
            return Err(invalid("composite display requires 12..=512 circumference segments"));
        }
        editor.validate()?;
        let cuts=editor.bindings();
        if cuts.len().saturating_mul(segments)>16384 {
            return Err(invalid("composite circle sample budget exceeded 16384 vertices"));
        }
        let source=editor.source_box();
        let model=editor.model();
        let bounds=source.bbox;
        let scale=len(bounds.max.sub(bounds.min));
        let largest=[bounds.min,bounds.max].into_iter()
            .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
            .fold(1.0_f64,f64::max);
        let epsilon=source.tolerance.length_at(scale)
            .max(32.0*f64::EPSILON*largest);
        let area_budget=source.tolerance.area_at(scale)*
            (2+segments.saturating_mul(cuts.len())) as f64+
            epsilon*scale*32.0;
        let mut mesh=Mesh{
            vertices:source.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        let mut outer=Vec::new();
        let mut lower_ids=Vec::new();
        let mut upper_ids=Vec::new();
        for vertex in &source.faces[0].boundary{
            let index=vertex.0 as usize;
            if index>=4 {return Err(invalid("composite box lower cap is not canonical"));}
            let lower=source.vertices[index].position;
            let upper=source.vertices[index+4].position;
            if (lower.x-upper.x).hypot(lower.z-upper.z)>epsilon ||
               (lower.y-bounds.min.y).abs()>epsilon ||
               (upper.y-bounds.max.y).abs()>epsilon{
                return Err(invalid("composite box top and bottom corners disagree"));
            }
            outer.push([lower.x-bounds.min.x,lower.z-bounds.min.z]);
            lower_ids.push(vertex.0);
            upper_ids.push(vertex.0+4);
        }
        let outer:[[_;2];4]=outer.try_into()
            .map_err(|_|invalid("composite cap must have exactly 4 corners"))?;
        let lower_ids:[u32;4]=lower_ids.try_into().unwrap();
        let upper_ids:[u32;4]=upper_ids.try_into().unwrap();
        if signed_area(&outer)<=0.0{
            return Err(invalid("composite box outer cap must have positive XZ winding"));
        }
        let mut rings=Vec::new();
        let mut polygon_removed=0.0;
        let mut side_area=0.0;
        let mut exact_depth_volume=0.0;
        let mut floor_centers=Vec::<(FaceId,[f64;2],f64)>::new();
        for binding in cuts{
            let lower_edge=model.geometry.edges.get(binding.topology.lower_edge.0 as usize)
                .ok_or_else(||invalid("missing composite lower edge geometry"))?;
            let upper_edge=model.geometry.edges.get(binding.topology.upper_edge.0 as usize)
                .ok_or_else(||invalid("missing composite upper edge geometry"))?;
            let Curve3::Circle(lower)=model.edge_curve(binding.topology.lower_edge)
                .ok_or_else(||invalid("missing composite lower 3D curve"))?
                else {return Err(invalid("composite lower cut edge is not circle"))};
            let Curve3::Circle(upper)=model.edge_curve(binding.topology.upper_edge)
                .ok_or_else(||invalid("missing composite upper 3D curve"))?
                else {return Err(invalid("composite upper cut edge is not circle"))};
            let range=|r:crate::ParameterRange|{
                match r {
                    crate::ParameterRange::Bounded{start,end}
                        if (end-start-TAU).abs()<=source.tolerance.angular=>Ok((start,end)),
                    _=>Err(invalid("composite circle domain must cover exactly one turn")),
                }
            };
            let (la,lb)=range(lower_edge.domain)?;
            let (ua,ub)=range(upper_edge.domain)?;
            let depth=binding.cut.depth(source);
            let bottom=bounds.max.y-depth;
            let center=binding.cut.center();
            let mut coords=Vec::with_capacity(segments);
            let mut low=Vec::with_capacity(segments);
            let mut high=Vec::with_capacity(segments);
            for i in 0..segments{
                let frac=i as f64/segments as f64;
                let a=lower.evaluate(la+(lb-la)*frac)?;
                let b=upper.evaluate(ua+(ub-ua)*frac)?;
                if (a.x-b.x).hypot(a.z-b.z)>epsilon ||
                   (a.y-bottom).abs()>epsilon ||
                   (b.y-bounds.max.y).abs()>epsilon ||
                   ((b.x-center[0]).hypot(b.z-center[1])-binding.cut.radius()).abs()>epsilon{
                    return Err(invalid("composite cut ring carrier does not lie on source cylinder"));
                }
                let idx=id(mesh.vertices.len())?;
                mesh.vertices.extend([a,b]);
                low.push(idx);high.push(idx+1);
                coords.push([a.x-bounds.min.x,a.z-bounds.min.z]);
            }
            let signed=signed_area(&coords);
            if !signed.is_finite() || signed>=0.0 {
                return Err(invalid("composite circle winding is not clockwise in XZ chart"));
            }
            let area=-signed;
            let perim=ring_perimeter(&coords);
            polygon_removed+=area*depth;
            side_area+=perim*depth;
            if binding.cut.is_through() {
                side_area-=2.0*area;
            }
            exact_depth_volume+=PI*binding.cut.radius().powi(2)*depth;
            if let Some(floor)=binding.topology.floor_face {
                floor_centers.push((floor,center,bottom));
            }
            rings.push(CircleRing{binding,coords,lower:low,upper:high,
                area,perimeter:perim});
        }
        // Top has EVERY opening. Bottom excludes BLIND bores.
        let top_refs=rings.iter().collect::<Vec<_>>();
        let bottom_refs=rings.iter().filter(|r|r.binding.cut.is_through())
            .collect::<Vec<_>>();
        let bottom_cert=cap(&mut mesh,FaceId(0),&outer,&lower_ids,
            &bottom_refs,false,area_budget,
            source.tolerance.area_at(scale))?;
        let top_cert=cap(&mut mesh,FaceId(1),&outer,&upper_ids,
            &top_refs,true,area_budget,
            source.tolerance.area_at(scale))?;
        // Original four box vertical face triangles are unchanged.
        for (tri,face) in source.mesh.triangles.iter().zip(&source.mesh.triangle_faces){
            if face.0>=2 {
                mesh.triangles.push(*tri);
                mesh.triangle_faces.push(*face);
            }
        }
        // All material-facing cylindrical normals point inward to the void.
        // The cap and wall share IDENTICAL ring vertex indices.
        for ring in &rings{
            for j in 0..segments{
                let k=(j+1)%segments;
                let a=ring.lower[j];let b=ring.upper[j];
                let c=ring.upper[k];let d=ring.lower[k];
                mesh.triangles.extend([[a,b,c],[a,c,d]]);
                mesh.triangle_faces.extend([
                    ring.binding.topology.side_face,
                    ring.binding.topology.side_face,
                ]);
            }
            if let Some(floor)=ring.binding.topology.floor_face {
                let [cx,cz]=ring.binding.cut.center();
                let y=bounds.max.y-ring.binding.cut.depth(source);
                let center=id(mesh.vertices.len())?;
                mesh.vertices.push(Point3{x:cx,y,z:cz});
                for j in 0..segments {
                    mesh.triangles.push([
                        center,ring.lower[j],ring.lower[(j+1)%segments]
                    ]);
                    mesh.triangle_faces.push(floor);
                }
            }
        }
        let (facet_volume,facet_area)=certify_facets(
            source,model,&rings,&mesh,epsilon)?;
        let expected_volume=source.mass.volume-polygon_removed;
        let expected_area=source.mass.surface_area+side_area;
        let volume_budget=source.mass.volume.abs()*1e-8+
            epsilon*scale*scale*16.0;
        let area_budget=source.mass.surface_area.abs()*1e-8+
            epsilon*scale*64.0;
        if !expected_volume.is_finite()||expected_volume<=0.0||
           !expected_area.is_finite()||expected_area<=0.0||
           (facet_volume-expected_volume).abs()>volume_budget ||
           (facet_area-expected_area).abs()>area_budget {
            return Err(invalid(format!(
                "composite facet integrals diverge from polygon-source oracle: volume {facet_volume} vs {expected_volume}, area {facet_area} vs {expected_area}")));
        }
        if (model.mass.volume-(source.mass.volume-exact_depth_volume)).abs()
            >volume_budget {
            return Err(invalid("composite analytic mass mismatches independent cylinder sum"));
        }
        if facet_volume+volume_budget<model.mass.volume{
            return Err(invalid("composite display underestimates exact material mass"));
        }
        Ok(CompositeDisplay{
            mesh,bottom:bottom_cert,top:top_cert,
            faceted_volume:facet_volume,faceted_surface_area:facet_area,
            analytic_volume:model.mass.volume,
            faceting_volume_error:facet_volume-model.mass.volume,
            segments,
        })
    }
}
fn certify_facets(source:&crate::Solid,model:&BrepModel,rings:&[CircleRing],
                  mesh:&Mesh,eps:f64)->Result<(f64,f64),GeometryError>
{
    if mesh.triangles.len()!=mesh.triangle_faces.len()||
       mesh.vertices.iter().any(|v|!v.is_finite()){
        return Err(invalid("composite display mesh lacks valid facet lineage"));
    }
    let mut edges=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let anchor=source.bbox.min;
    let mut signed_six=0.0;
    let mut area=0.0;
    for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        let [a,b,c]=*tri;
        if a==b||b==c||c==a ||
           [a,b,c].iter().any(|&i|i as usize>=mesh.vertices.len()){
            return Err(invalid("composite display degenerate or invalid facet indices"));
        }
        let [p,q,r]=tri.map(|i|mesh.vertices[i as usize]);
        let normal=q.sub(p).cross(r.sub(p));
        let norm=len(normal);
        if !norm.is_finite()||norm<=0.0{
            return Err(invalid("composite display has degenerate 3D facet"));
        }
        if model.face_sense(*owner).is_none(){
            return Err(invalid("composite display facet has no owner face"));
        }
        if owner.0>=6 {
            let ring=rings.iter().find(|ring|
                ring.binding.topology.side_face==*owner ||
                ring.binding.topology.floor_face==Some(*owner))
                .ok_or_else(||invalid("composite facet owned by an unknown analytic cut"))?;
            if model.face_sense(*owner)!=Some(FaceSense::Reversed){
                return Err(invalid("composite cavity face must reverse analytic support"));
            }
            if ring.binding.topology.side_face==*owner{
                let cx=ring.binding.cut.center()[0];
                let cz=ring.binding.cut.center()[1];
                let centroid=p.add(q).add(r).scale(1.0/3.0);
                let inward=Point3{x:cx-centroid.x,y:0.0,z:cz-centroid.z};
                if normal.dot(inward)<=0.0 ||
                   [p,q,r].iter().any(|v|{
                       ((v.x-cx).hypot(v.z-cz)-ring.binding.cut.radius()).abs()>eps
                   }){
                    return Err(invalid("composite cylindrical facet not inward or off analytic radius"));
                }
            }else{
                let floor=source.bbox.max.y-ring.binding.cut.depth(source);
                if normal.y<=0.0 ||
                   [p,q,r].iter().any(|v|(v.y-floor).abs()>eps){
                    return Err(invalid("composite blind disk facet must face +Y on its floor plane"));
                }
            }
        }else{
            let Some(Surface3::Plane(plane))=model.face_surface(*owner)else{
                return Err(invalid("composite outer facet does not bind plane"));
            };
            if model.face_sense(*owner)!=Some(FaceSense::Forward) ||
                normal.dot(plane.frame.normal)<=0.0 ||
               [p,q,r].iter().any(|v|
                   (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps){
                return Err(invalid("composite outer facet is not an outward box plane"));
            }
        }
        for (u,v) in [(a,b),(b,c),(c,a)]{
            edges.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
        signed_six+=p.sub(anchor).dot(q.sub(anchor).cross(r.sub(anchor)));
        area+=norm*0.5;
    }
    if edges.values().any(|uses|uses.len()!=2 || uses[0]==uses[1]){
        return Err(invalid("composite display mesh fails closed two-use directed-edge certificate"));
    }
    let volume=signed_six/6.0;
    if !volume.is_finite()||volume<=0.0||
       !area.is_finite()||area<=0.0{
        return Err(invalid("composite display has invalid signed volume or area"));
    }
    Ok((volume,area))
}
