//! Independently checked watertight display mesh for the strictly scoped
//! one-through/one-blind transverse-crossing analytic B-rep.
//! The source B-rep never uses these facets as its geometry or mass oracle.
use std::collections::BTreeMap;
use std::f64::consts::{PI,TAU};

use crate::{
    Curve3,FaceId,GeometryError,Mesh,MixedCrossingBrep,Point3,
};

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
fn psub(a:Point3,b:Point3)->Point3{
    Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}
}
fn cross(a:Point3,b:Point3)->Point3{
    Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}
}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn mag(a:Point3)->f64{a.x.hypot(a.y).hypot(a.z)}
fn dist(a:Point3,b:Point3)->f64{mag(psub(a,b))}
fn shoelace(mesh:&Mesh,loop_:&[u32])->f64{
    let mut sum=0.0;
    let origin=mesh.vertices[loop_[0] as usize];
    for i in 0..loop_.len(){
        let a=psub(mesh.vertices[loop_[i] as usize],origin);
        let b=psub(mesh.vertices[loop_[(i+1)%loop_.len()] as usize],origin);
        sum+=a.x*b.z-a.z*b.x;
    }
    sum/2.0
}
fn cycle(a:&[u32],b:&[u32])->Vec<u32>{
    let mut out=a[..a.len()-1].to_vec();
    out.extend_from_slice(&b[..b.len()-1]);
    out
}
fn sample_arc(
    mesh:&mut Mesh,curve:Curve3,u0:f64,u1:f64,start:u32,end:u32,
    segments:usize,
)->Result<Vec<u32>,GeometryError>{
    let mut ids=Vec::with_capacity(segments+1);
    ids.push(start);
    for i in 1..segments{
        let parameter=u0+(u1-u0)*i as f64/segments as f64;
        let point=curve.evaluate(parameter)?;
        let index=u32::try_from(mesh.vertices.len())
            .map_err(|_|invalid("display vertices exceed u32 budget"))?;
        mesh.vertices.push(point);
        ids.push(index);
    }
    ids.push(end);
    Ok(ids)
}
fn add_triangle(mesh:&mut Mesh,mut ids:[u32;3],face:FaceId,
                desired_y:i8,tolerance:f64)->Result<(),GeometryError>{
    let [a,b,c]=ids.map(|i|mesh.vertices[i as usize]);
    let n=cross(psub(b,a),psub(c,a));
    if !n.y.is_finite()||n.y.abs()<=tolerance {
        return Err(invalid("degenerate or nonplanar cap triangle"));
    }
    if (n.y>0.0)!=(desired_y>0){ids.swap(1,2);}
    mesh.triangles.push(ids);
    mesh.triangle_faces.push(face);
    Ok(())
}
fn triangulate_cap(
    mesh:&mut Mesh,outer:&[u32],inner:Option<&[u32]>,
    face:FaceId,desired_y:i8,area_eps:f64,
)->Result<f64,GeometryError>{
    let mut indices=outer.to_vec();
    let mut holes=Vec::new();
    if let Some(loop_)=inner{
        holes.push(indices.len());
        indices.extend_from_slice(loop_);
    }
    let mut flat=Vec::with_capacity(indices.len()*2);
    for id in &indices {
        let p=mesh.vertices[*id as usize];
        flat.extend_from_slice(&[p.x,p.z]);
    }
    let triangles=earcutr::earcut(&flat,&holes,2)
        .map_err(|e|invalid(format!("planar cap triangulation failed: {e:?}")))?;
    if triangles.is_empty()||triangles.len()%3!=0{
        return Err(invalid("planar cap returned zero or incomplete triangles"));
    }
    let mut area=0.0;
    for tri in triangles.chunks_exact(3){
        let ids=[indices[tri[0]],indices[tri[1]],indices[tri[2]]];
        if ids[0]==ids[1]||ids[1]==ids[2]||ids[2]==ids[0]{
            return Err(invalid("planar triangulation reused a vertex"));
        }
        let [a,b,c]=ids.map(|i|mesh.vertices[i as usize]);
        let signed=(b.x-a.x)*(c.z-a.z)-(b.z-a.z)*(c.x-a.x);
        if !signed.is_finite()||signed.abs()<=area_eps {
            return Err(invalid("planar cap triangle area below tolerance"));
        }
        area+=signed.abs()/2.0;
        add_triangle(mesh,ids,face,desired_y,area_eps)?;
    }
    let expected=(shoelace(mesh,outer).abs()-
        inner.map(|r|shoelace(mesh,r).abs()).unwrap_or(0.0)).abs();
    let allowable=area_eps*(triangles.len()+indices.len()) as f64*32.0;
    if (area-expected).abs()>allowable{
        return Err(invalid("triangulated cap area differs from polygonal source contour"));
    }
    Ok(area)
}
fn wall(mesh:&mut Mesh,lower:&[u32],upper:&[u32],face:FaceId,
        tol:f64)->Result<f64,GeometryError>{
    if lower.len()<2||lower.len()!=upper.len(){
        return Err(invalid("wall curves require equal polygonization"));
    }
    let mut perimeter=0.0;
    for i in 0..lower.len()-1{
        let a=lower[i];let b=lower[i+1];
        let c=upper[i+1];let d=upper[i];
        if a==b||c==d{return Err(invalid("degenerate cylinder segment"));}
        let bottom=mesh.vertices[a as usize];let other=mesh.vertices[b as usize];
        let height=dist(bottom,mesh.vertices[d as usize]);
        let chord=dist(bottom,other);
        if chord<=tol || height<=tol{return Err(invalid("wall below feature tolerance"));}
        perimeter+=chord;
        // Increasing world XZ angle crossed with +Y gives inward radial
        // normal: this is the physical material side of a bored void.
        mesh.triangles.push([a,b,c]);
        mesh.triangle_faces.push(face);
        mesh.triangles.push([a,c,d]);
        mesh.triangle_faces.push(face);
    }
    Ok(perimeter)
}

/// Certified display only: committed analytic volume and faceted volume are
/// independent, and their difference is explicitly reported.
#[derive(Debug,Clone)]
pub struct MixedCrossingDisplay {
    pub mesh:Mesh,
    pub facet_volume:f64,
    pub facet_area:f64,
    pub analytic_volume:f64,
    pub facet_volume_error:f64,
    pub expected_polygonal_volume:f64,
    pub top_cap_polygon_area:f64,
    pub bottom_cap_polygon_area:f64,
    pub floor_polygon_area:f64,
    pub triangle_edges:usize,
    pub genus:usize,
    pub circle_segments_per_arc:usize,
}
pub struct MixedCrossingMesher;
impl MixedCrossingMesher {
    pub fn tessellate(body:&MixedCrossingBrep,segments:usize)
        ->Result<MixedCrossingDisplay,GeometryError>
    {
        if !(12..=512).contains(&segments){
            return Err(invalid("crossing display needs 12..=512 segments per analytic arc"));
        }
        body.validate()?;
        let model=body.model();
        let source=body.source_box();
        let arc=body.plan().footprint.exposed_arcs
            .ok_or_else(||invalid("missing source arc partition"))?;
        let points=body.plan().footprint.crossing_points
            .ok_or_else(||invalid("missing source crossings"))?;
        let mut mesh=Mesh{
            vertices:source.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        // First-order physical junctions at bottom, floor, top: shared by
        // both circle families and every face using the same crossing.
        let mut junctions=[[0u32;2];3];
        for (layer,y) in [
            source.bbox.min.y,body.plan().splitting_height_y,source.bbox.max.y
        ].into_iter().enumerate(){
            for (index,which) in [arc[0].start_crossing,arc[0].end_crossing]
                .into_iter().enumerate(){
                let p=points[which];
                junctions[layer][index]=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("crossing display vertex budget exceeded"))?;
                mesh.vertices.push(Point3{x:p.x,y,z:p.z});
            }
        }
        let [a0,a1]=[arc[0].start_angle,arc[0].start_angle+arc[0].sweep_radians];
        let [c0,c1]=[a1,arc[0].start_angle+TAU];
        let [b0,b1]=[arc[1].start_angle,arc[1].start_angle+arc[1].sweep_radians];
        // Analytic curves were constructed from the same source angles.
        // Downward/upward donor rings sample matching angular fractions,
        // so actual physical mesh edges are shared, not merely coincident.
        let edge=|id:u32|->Result<Curve3,GeometryError>{
            model.edge_curve(crate::EdgeId(id))
                .ok_or_else(||invalid("crossing display edge lacks carrier"))
        };
        let bottom_to=sample_arc(&mut mesh,edge(12)?,a0,a1,junctions[0][0],junctions[0][1],segments)?;
        let bottom_ti=sample_arc(&mut mesh,edge(13)?,c0,c1,junctions[0][1],junctions[0][0],segments)?;
        let floor_to=sample_arc(&mut mesh,edge(14)?,a0,a1,junctions[1][0],junctions[1][1],segments)?;
        let floor_ti=sample_arc(&mut mesh,edge(15)?,c0,c1,junctions[1][1],junctions[1][0],segments)?;
        let floor_bo=sample_arc(&mut mesh,edge(16)?,b0,b1,junctions[1][1],junctions[1][0],segments)?;
        let top_to=sample_arc(&mut mesh,edge(17)?,a0,a1,junctions[2][0],junctions[2][1],segments)?;
        let top_bo=sample_arc(&mut mesh,edge(18)?,b0,b1,junctions[2][1],junctions[2][0],segments)?;
        let top_aperture=cycle(&top_to,&top_bo);
        let lower_aperture=cycle(&bottom_to,&bottom_ti);
        let mut reversed_b=floor_bo.clone();reversed_b.reverse();
        let floor_crescent=cycle(&reversed_b,&floor_ti);
        let bbox=source.bbox;
        let extent=(bbox.max.x-bbox.min.x).max(bbox.max.y-bbox.min.y)
            .max(bbox.max.z-bbox.min.z);
        let area_eps=source.tolerance.area_at(extent)*0.01;
        let eps=source.tolerance.length_at(extent);
        let bottom_outer=source.faces[0].boundary.iter().map(|v|v.0).collect::<Vec<_>>();
        let top_outer=source.faces[1].boundary.iter().map(|v|v.0).collect::<Vec<_>>();
        let bottom_area=triangulate_cap(&mut mesh,&bottom_outer,Some(&lower_aperture),
            FaceId(0),-1,area_eps)?;
        let top_area=triangulate_cap(&mut mesh,&top_outer,Some(&top_aperture),
            FaceId(1),1,area_eps)?;
        let floor_area=triangulate_cap(&mut mesh,&floor_crescent,None,
            FaceId(10),1,area_eps)?;
        let mut facet_side_area=0.0;
        for (tri,face) in source.mesh.triangles.iter().zip(&source.mesh.triangle_faces){
            if face.0<2 {continue;}
            mesh.triangles.push(*tri);
            mesh.triangle_faces.push(*face);
            let [a,b,c]=tri.map(|i|mesh.vertices[i as usize]);
            facet_side_area+=mag(cross(psub(b,a),psub(c,a)))/2.0;
        }
        let perimeter_lower_to=wall(&mut mesh,&bottom_to,&floor_to,FaceId(6),eps)?;
        let perimeter_lower_ti=wall(&mut mesh,&bottom_ti,&floor_ti,FaceId(7),eps)?;
        let perimeter_upper_to=wall(&mut mesh,&floor_to,&top_to,FaceId(8),eps)?;
        let perimeter_upper_bo=wall(&mut mesh,&floor_bo,&top_bo,FaceId(9),eps)?;
        let hl=body.plan().splitting_height_y-bbox.min.y;
        let hu=bbox.max.y-body.plan().splitting_height_y;
        let target_area=top_area+bottom_area+floor_area+facet_side_area+
            hl*(perimeter_lower_to+perimeter_lower_ti)+
            hu*(perimeter_upper_to+perimeter_upper_bo);
        // Every triangle edge must have EXACTLY two opposite oriented
        // uses. This is a stronger mesh closure check than Euler alone.
        let mut usage=BTreeMap::<(u32,u32),Vec<bool>>::new();
        let mut signed_volume=0.0;
        let mut facet_area=0.0;
        let c0=source.mass.centroid;
        for tri in &mesh.triangles{
            let [a,b,c]=tri.map(|i|mesh.vertices[i as usize]);
            let va=psub(a,c0);let vb=psub(b,c0);let vc=psub(c,c0);
            signed_volume+=dot(va,cross(vb,vc))/6.0;
            facet_area+=mag(cross(psub(b,a),psub(c,a)))/2.0;
            for (x,y) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]{
                if x==y{return Err(invalid("display has degenerate edge"));}
                usage.entry((x.min(y),x.max(y))).or_default().push(x<y);
            }
        }
        if usage.values().any(|directions|
            directions.len()!=2||directions[0]==directions[1]){
            return Err(invalid("crossing display is not a closed oriented watertight mesh"));
        }
        let expected_polygon_volume=hl*bottom_area+hu*top_area;
        let facet_volume_error=signed_volume-model.mass.volume;
        let slack=(source.tolerance.length_at(extent)*extent*extent)
            .max(1e-8*source.mass.volume.abs());
        if !signed_volume.is_finite()||signed_volume<=0.0||
            (signed_volume-expected_polygon_volume).abs()>slack||
            (facet_area-target_area).abs()>slack||
            facet_volume_error< -slack {
            return Err(invalid("independent faceted volume/area certificate failed"));
        }
        let chi=mesh.vertices.len() as isize -usage.len() as isize+
            mesh.triangles.len() as isize;
        if chi!=0 {return Err(invalid("crossing display genus-one Euler relation failed"));}
        if mesh.triangles.len()!=mesh.triangle_faces.len(){
            return Err(invalid("facet ownership count differs from mesh triangles"));
        }
        if mesh.triangle_faces.iter().any(|id|id.0>=11){
            return Err(invalid("facet references nonexistent source face"));
        }
        Ok(MixedCrossingDisplay{
            mesh,facet_volume:signed_volume,facet_area,
            analytic_volume:model.mass.volume,facet_volume_error,
            expected_polygonal_volume:expected_polygon_volume,
            top_cap_polygon_area:top_area,bottom_cap_polygon_area:bottom_area,
            floor_polygon_area:floor_area,triangle_edges:usage.len(),
            genus:1,circle_segments_per_arc:segments,
        })
    }
}
