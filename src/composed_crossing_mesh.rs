//! Faceted DISPLAY of a source-certified multi-cluster analytical B-rep.
//! Each independently certified donor supplies real curved side/floor facets.
//! The original block top and bottom are triangulated ONCE, with all true
//! source aperture loops and globally shared donor boundary vertices.
use std::collections::BTreeMap;
use crate::{
    BlindHoleBrep,ComposedCrossingBrep,FaceId,GeometryError,Mesh,
    MixedCrossingBrep,Point3,
};
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn point_sub(a:Point3,b:Point3)->Point3{
    Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}
}
fn cross(a:Point3,b:Point3)->Point3{
    Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}
}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn norm(a:Point3)->f64{a.x.hypot(a.y).hypot(a.z)}
fn signed_area(mesh:&Mesh,points:&[u32])->Result<f64,GeometryError>{
    if points.len()<3{return Err(invalid("cap loop too small"));}
    let first=mesh.vertices[points[0] as usize];
    let mut sum=0.0;
    for k in 0..points.len(){
        let a=point_sub(mesh.vertices[points[k] as usize],first);
        let b=point_sub(mesh.vertices[points[(k+1)%points.len()] as usize],first);
        sum+=a.x*b.z-a.z*b.x;
    }
    Ok(sum/2.0)
}

/// Extract an actual single hole boundary from a donor's certified planar
/// cap facets, not a fresh or independently sampled approximation. This
/// enforces the same physical vertex IDs across donor curved walls and cap.
fn donor_cap_ring(mesh:&Mesh,face:FaceId)->Result<Vec<u32>,GeometryError>{
    let mut edges=BTreeMap::<(u32,u32),(usize,(u32,u32))>::new();
    for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        if *owner!=face{continue;}
        for (a,b) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]{
            let entry=edges.entry((a.min(b),a.max(b))).or_insert((0,(a,b)));
            entry.0+=1;
        }
    }
    let mut next=BTreeMap::new();
    let mut incoming=BTreeMap::<u32,usize>::new();
    for (_, (count,(a,b))) in edges {
        if count>2{return Err(invalid("donor cap edge has nonmanifold incidence"));}
        if count==1 && a>=8 && b>=8{
            if next.insert(a,b).is_some(){return Err(invalid("cap hole branches"));}
            *incoming.entry(b).or_default()+=1;
        }
    }
    if next.len()<12 || incoming.len()!=next.len() ||
        incoming.values().any(|&n|n!=1){
        return Err(invalid("donor has no single complete sampled inner boundary"));
    }
    let start=*next.keys().next().ok_or_else(||invalid("missing cap cycle"))?;
    let mut ring=Vec::with_capacity(next.len());
    let mut at=start;
    while ring.len()<next.len(){
        ring.push(at);
        at=*next.get(&at).ok_or_else(||invalid("missing cap ring successor"))?;
        if at==start{break;}
    }
    if at!=start || ring.len()!=next.len(){
        return Err(invalid("donor cap boundary cycle is not closed"));
    }
    Ok(ring)
}
fn remap_donor(
    target:&mut Mesh,donor:&Mesh,face_offset:u32,has_bottom:bool
)->Result<(Vec<u32>,Option<Vec<u32>>),GeometryError>{
    if donor.vertices.len()<8 ||
        target.vertices.len()<8 ||
        donor.vertices[..8]!=target.vertices[..8]{
        return Err(invalid("donor display does not share exact source box vertices"));
    }
    let mut vmap=(0u32..8).collect::<Vec<_>>();
    for p in donor.vertices.iter().skip(8) {
        let id=u32::try_from(target.vertices.len())
            .map_err(|_|invalid("display vertex budget overflow"))?;
        target.vertices.push(*p);
        vmap.push(id);
    }
    let top=donor_cap_ring(donor,FaceId(1))?;
    let bottom=if has_bottom {
        Some(donor_cap_ring(donor,FaceId(0))?)
    }else{None};
    let translate=|ids:Vec<u32>| -> Result<Vec<u32>,GeometryError>{
        ids.into_iter().map(|i|vmap.get(i as usize).copied()
            .ok_or_else(||invalid("donor loop contains missing vertex")))
            .collect()
    };
    for (tri,owner) in donor.triangles.iter().zip(&donor.triangle_faces){
        if owner.0<6 {continue;}
        let remapped=tri.map(|i|vmap.get(i as usize).copied());
        let [Some(a),Some(b),Some(c)]=remapped else{
            return Err(invalid("donor has invalid triangle vertex"));
        };
        let new_face=owner.0.checked_sub(6).and_then(|x|face_offset.checked_add(x))
            .ok_or_else(||invalid("donor face remap overflow"))?;
        target.triangles.push([a,b,c]);
        target.triangle_faces.push(FaceId(new_face));
    }
    Ok((translate(top)?,bottom.map(translate).transpose()?))
}
fn triangulate_cap(
    mesh:&mut Mesh,outer:&[u32],holes:&[Vec<u32>],
    owner:FaceId,sign_y:f64,eps_area:f64
)->Result<f64,GeometryError>{
    let mut indices=outer.to_vec();
    let mut hole_starts=Vec::new();
    for hole in holes{
        if hole.len()<12{return Err(invalid("composed cap has malformed hole")); }
        hole_starts.push(indices.len());
        indices.extend_from_slice(hole);
    }
    let p0=mesh.vertices[indices[0] as usize];
    let mut coords=Vec::with_capacity(indices.len()*2);
    for &id in &indices{
        let p=mesh.vertices[id as usize];
        coords.extend_from_slice(&[p.x-p0.x,p.z-p0.z]);
    }
    let triangles=earcutr::earcut(&coords,&hole_starts,2)
        .map_err(|e|invalid(format!("multi-opening cap triangulation failed: {e:?}")))?;
    if triangles.is_empty()||triangles.len()%3!=0{
        return Err(invalid("cap tessellator returned no triangles"));
    }
    let mut area=0.0;
    for t in triangles.chunks_exact(3){
        let mut tri=[indices[t[0]],indices[t[1]],indices[t[2]]];
        if tri[0]==tri[1]||tri[1]==tri[2]||tri[2]==tri[0]{
            return Err(invalid("degenerate cap facet"));
        }
        let [a,b,c]=tri.map(|i|mesh.vertices[i as usize]);
        let n=cross(point_sub(b,a),point_sub(c,a));
        if !n.y.is_finite()||n.y.abs()<=eps_area {
            return Err(invalid("new cap triangle below area budget"));
        }
        area+=n.y.abs()/2.0;
        if (n.y>0.0)!=(sign_y>0.0){tri.swap(1,2);}
        mesh.triangles.push(tri);
        mesh.triangle_faces.push(owner);
    }
    let expected=signed_area(mesh,outer)?.abs()-
        holes.iter().map(|h|signed_area(mesh,h).map(f64::abs))
            .collect::<Result<Vec<_>,_>>()?.iter().sum::<f64>();
    if !expected.is_finite()||expected<=0.0 ||
        (area-expected).abs()>eps_area*(indices.len()+triangles.len()) as f64*64.0{
        return Err(invalid("triangulated cap does not match source cut-loop areas"));
    }
    Ok(area)
}
#[derive(Debug,Clone)]
pub struct ComposedCrossingDisplay{
    pub mesh:Mesh,
    pub facet_volume:f64,
    pub facet_surface_area:f64,
    pub analytic_volume:f64,
    pub volume_faceting_error:f64,
    pub expected_polygonal_volume:f64,
    pub triangle_edges:usize,
    pub genus:usize,
    pub crossings:usize,
    pub isolated_blind:usize,
    pub segments_per_source_arc:usize,
}
pub struct ComposedCrossingMesher;
impl ComposedCrossingMesher{
    pub fn tessellate(body:&ComposedCrossingBrep,segments:usize)
        ->Result<ComposedCrossingDisplay,GeometryError>
    {
        if !(12..=256).contains(&segments){
            return Err(invalid("composed display needs 12..=256 arc segments"));
        }
        body.validate()?;
        let source=body.source_box();
        let params=body.parameters();
        let expected_faces=body.model().faces.len();
        let mut mesh=Mesh{vertices:source.vertices.iter().map(|v|v.position)
            .collect(),triangles:Vec::new(),triangle_faces:Vec::new()};
        let mut top_holes=Vec::new();
        let mut bottom_holes=Vec::new();
        let mut face_offset=6u32;
        let mut volume_removed=0.0;
        let mut additional_area=0.0;
        let v0=source.mass.volume;
        for pair in &params.pairs{
            let donor=MixedCrossingBrep::from_box(source,pair.through,pair.blind)?;
            let display=donor.tessellate(segments)?;
            let (top,bottom)=remap_donor(&mut mesh,&display.mesh,face_offset,true)?;
            top_holes.push(top);
            bottom_holes.push(bottom.ok_or_else(||invalid("lost crossing bottom aperture"))?);
            face_offset+=5;
            volume_removed+=v0-display.facet_volume;
            additional_area+=display.facet_area-source.mass.surface_area;
        }
        for blind in &params.isolated_blind{
            let donor=BlindHoleBrep::from_box(source,*blind)?;
            let display=donor.tessellate(segments)?;
            let (top,bottom)=remap_donor(&mut mesh,&display.mesh,face_offset,false)?;
            if bottom.is_some(){return Err(invalid("isolated blind unexpectedly crosses bottom"));}
            top_holes.push(top);
            face_offset+=2;
            volume_removed+=v0-display.faceted_volume;
            additional_area+=display.faceted_surface_area-source.mass.surface_area;
        }
        if face_offset as usize!=expected_faces{
            return Err(invalid("source face mapping does not match shell topology"));
        }
        let base=source.faces[0].boundary.iter().map(|x|x.0).collect::<Vec<_>>();
        let top=source.faces[1].boundary.iter().map(|x|x.0).collect::<Vec<_>>();
        let span=source.bbox.max.sub(source.bbox.min);
        let extent=span.x.max(span.y).max(span.z);
        let eps_area=source.tolerance.area_at(extent)*0.02;
        triangulate_cap(&mut mesh,&base,&bottom_holes,FaceId(0),-1.0,eps_area)?;
        triangulate_cap(&mut mesh,&top,&top_holes,FaceId(1),1.0,eps_area)?;
        for (tri,face) in source.mesh.triangles.iter().zip(&source.mesh.triangle_faces){
            if face.0<2{continue;}
            mesh.triangles.push(*tri);
            mesh.triangle_faces.push(*face);
        }
        let mut adjacency=BTreeMap::<(u32,u32),Vec<bool>>::new();
        let center=source.mass.centroid;
        let mut volume=0.0;
        let mut area=0.0;
        for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
            if owner.0 as usize>=expected_faces{return Err(invalid("unknown analytic face owner"));}
            let [a,b,c]=tri.map(|x|mesh.vertices[x as usize]);
            let a=point_sub(a,center);
            let b=point_sub(b,center);
            let c=point_sub(c,center);
            let n=cross(b,c);
            volume+=dot(a,n)/6.0;
            area+=norm(cross(point_sub(b,a),point_sub(c,a)))/2.0;
            for (i,j) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]{
                if i==j{return Err(invalid("display has a zero-length vertex edge"));}
                adjacency.entry((i.min(j),i.max(j))).or_default().push(i<j);
            }
        }
        if adjacency.values().any(|v|v.len()!=2||v[0]==v[1]){
            return Err(invalid("multi-opening display is not a closed oriented manifold"));
        }
        let euler=mesh.vertices.len() as isize-
            adjacency.len() as isize+mesh.triangles.len() as isize;
        let target_euler=2-2*params.pairs.len() as isize;
        if euler!=target_euler{
            return Err(invalid("multi-opening display Euler characteristic disagrees with shell"));
        }
        let expected_polygonal_volume=v0-volume_removed;
        let expected_area=source.mass.surface_area+additional_area;
        let roundoff=v0.abs()*1e-8+
            source.tolerance.length_at(extent)*extent*extent*64.0;
        if !volume.is_finite()||!area.is_finite()||volume<=0.0||
           (volume-expected_polygonal_volume).abs()>roundoff||
           (area-expected_area).abs()>roundoff*16.0||
           volume<body.model().mass.volume-roundoff {
            return Err(invalid("multi-component polygonized volume/area certificate failed"));
        }
        if mesh.triangles.len()!=mesh.triangle_faces.len(){
            return Err(invalid("mesh provenance length differs from facets"));
        }
        Ok(ComposedCrossingDisplay{
            mesh,facet_volume:volume,facet_surface_area:area,
            analytic_volume:body.model().mass.volume,
            volume_faceting_error:volume-body.model().mass.volume,
            expected_polygonal_volume,triangle_edges:adjacency.len(),
            genus:params.pairs.len(),crossings:params.pairs.len(),
            isolated_blind:params.isolated_blind.len(),
            segments_per_source_arc:segments
        })
    }
}
