//! Strictly derived viewport triangulation for the v0.35 connected 3-cylinder
//! analytic B-rep. Shared physical edges own ONE sample polyline globally.
//! The analytic B-rep, volume and centroid are never defined by triangles.
use std::collections::BTreeMap;
use crate::{
    BrepModel,ConnectedTripleBrep,Curve3,EdgeId,FaceId,GeometryError,
    Mesh,Point3,Surface3,
};
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn sub(a:Point3,b:Point3)->Point3{
    Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}
}
fn cross(a:Point3,b:Point3)->Point3{
    Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}
}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn norm(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn signed_loop(mesh:&Mesh,poly:&[u32])->f64{
    let first=mesh.vertices[poly[0] as usize];
    (0..poly.len()).map(|i|{
        let a=sub(mesh.vertices[poly[i] as usize],first);
        let b=sub(mesh.vertices[poly[(i+1)%poly.len()] as usize],first);
        (a.x*b.z-a.z*b.x)*0.5
    }).sum()
}
fn build_edge_samples(m:&BrepModel,mesh:&mut Mesh,segments:usize)
    ->Result<Vec<Vec<u32>>,GeometryError>
{
    let mut result=Vec::with_capacity(m.edges.len());
    let span=m.bbox.max.sub(m.bbox.min);
    let eps=m.tolerance.length_at(span.x.max(span.y).max(span.z));
    for e in &m.edges{
        let c=m.edge_curve(e.id).ok_or_else(||invalid("missing 3D edge carrier"))?;
        let crate::ParameterRange::Bounded{start,end}=m.geometry.edges[e.id.0 as usize].domain
            else{return Err(invalid("display needs bounded physical edge"))};
        let n=if matches!(c,Curve3::Circle(_)){segments}else{1};
        let mut sampled=Vec::with_capacity(n+1);
        sampled.push(e.start.0);
        for i in 1..n{
            let t=start+(end-start)*i as f64/n as f64;
            let p=c.evaluate(t)?;
            if !p.is_finite(){return Err(invalid("nonfinite analytic edge display sample"));}
            let id=u32::try_from(mesh.vertices.len())
                .map_err(|_|invalid("display vertex index capacity exceeded"))?;
            mesh.vertices.push(p);sampled.push(id);
        }
        sampled.push(e.end.0);
        for (t,id) in [(start,e.start.0),(end,e.end.0)]{
            let expected=c.evaluate(t)?;
            let actual=mesh.vertices[id as usize];
            if norm(sub(expected,actual))>eps*16.0{
                return Err(invalid("source carrier endpoint does not match physical B-rep vertex"));
            }
        }
        result.push(sampled);
    }
    Ok(result)
}
fn loop_ring(m:&BrepModel,loopid:crate::LoopId,polylines:&[Vec<u32>])
    ->Result<Vec<u32>,GeometryError>
{
    let lp=&m.loops[loopid.0 as usize];
    let mut cursor=lp.first_coedge;
    let mut ring=Vec::new();
    loop {
        let c=&m.coedges[cursor.0 as usize];
        let mut edge=polylines[c.edge.0 as usize].clone();
        if c.reversed{edge.reverse();}
        ring.extend_from_slice(&edge[..edge.len()-1]);
        cursor=c.next;
        if cursor==lp.first_coedge{break;}
        if ring.len()>32768{return Err(invalid("unbounded display face loop"));}
    }
    if ring.len()<3||ring.windows(2).any(|w|w[0]==w[1]){
        return Err(invalid("display face loop has degenerate sample sequence"));
    }
    Ok(ring)
}
fn tessellate_planar(m:&BrepModel,mesh:&mut Mesh,face:FaceId,
    polylines:&[Vec<u32>],eps:f64)->Result<f64,GeometryError>
{
    let f=&m.faces[face.0 as usize];
    let mut rings=Vec::new();
    for &lp in &f.loops{rings.push(loop_ring(m,lp,polylines)?);}
    let mut flat=Vec::new();
    let mut indices=Vec::new();
    let mut hole_starts=Vec::new();
    let origin=mesh.vertices[rings[0][0] as usize];
    for (i,ring) in rings.iter().enumerate(){
        if i>0{hole_starts.push(indices.len());}
        for &v in ring{
            let p=mesh.vertices[v as usize];
            if (p.y-origin.y).abs()>eps*8.0{
                return Err(invalid("planar face loop departs from its supporting floor height"));
            }
            flat.extend_from_slice(&[p.x-origin.x,p.z-origin.z]);
            indices.push(v);
        }
    }
    let triangles=earcutr::earcut(&flat,&hole_starts,2)
        .map_err(|e|invalid(format!("analytic planar face tessellation failed: {e:?}")))?;
    if triangles.is_empty()||triangles.len()%3!=0 {
        return Err(invalid("analytic planar face triangulator returned no complete triangles"));
    }
    let desired_y=if face.0==0{-1.0}else{1.0};
    let mut area=0.0;
    for tri in triangles.chunks_exact(3){
        let mut ids=[indices[tri[0]],indices[tri[1]],indices[tri[2]]];
        let [a,b,c]=ids.map(|i|mesh.vertices[i as usize]);
        let normal=cross(sub(b,a),sub(c,a));
        if !normal.y.is_finite()||normal.y.abs()<=eps*eps*0.001 {
            return Err(invalid("planar triangle below geometric resolution"));
        }
        area+=normal.y.abs()*0.5;
        if normal.y.signum()!=desired_y {ids.swap(1,2);}
        mesh.triangles.push(ids);
        mesh.triangle_faces.push(face);
    }
    let expected=signed_loop(mesh,&rings[0]).abs()-
        rings.iter().skip(1).map(|r|signed_loop(mesh,r).abs()).sum::<f64>();
    let budget=eps*(m.bbox.max.x-m.bbox.min.x).max(m.bbox.max.z-m.bbox.min.z)*
        (triangles.len()+indices.len()) as f64*24.0;
    if !expected.is_finite()||expected<=0.0||
        (expected-area).abs()>budget{
        return Err(invalid("planar facet area disagrees with exact source contour polygon"));
    }
    Ok(area)
}
fn tessellate_cylinder(m:&BrepModel,mesh:&mut Mesh,face:FaceId,
    polylines:&[Vec<u32>])->Result<f64,GeometryError>
{
    let f=&m.faces[face.0 as usize];
    if f.loops.len()!=1{return Err(invalid("cylinder patch cannot contain secondary loop"));}
    let mut cursor=m.loops[f.loops[0].0 as usize].first_coedge;
    let mut circles=Vec::<EdgeId>::new();
    let mut sides=0usize;
    loop{
        let c=&m.coedges[cursor.0 as usize];
        match m.edge_curve(c.edge){
            Some(Curve3::Circle(_))=>circles.push(c.edge),
            Some(Curve3::Line(_))=>sides+=1,
            _=>return Err(invalid("cylinder patch has no analytic edge geometry"))
        }
        cursor=c.next;
        if cursor==m.loops[f.loops[0].0 as usize].first_coedge{break;}
    }
    if circles.len()!=2||sides!=2{
        return Err(invalid("cylinder patch is not a closed 4-arc trim"));
    }
    let first=&polylines[circles[0].0 as usize];
    let second=&polylines[circles[1].0 as usize];
    if first.len()!=second.len(){return Err(invalid("cylinder wall sampling intervals differ"));}
    let (lower,upper)=if mesh.vertices[first[0] as usize].y <
        mesh.vertices[second[0] as usize].y {
        (first,second)
    }else{(second,first)};
    let mut area=0.0;
    for i in 0..lower.len()-1{
        let (a,b,c,d)=(lower[i],lower[i+1],upper[i+1],upper[i]);
        for tri in [[a,b,c],[a,c,d]]{
            let [u,v,w]=tri.map(|id|mesh.vertices[id as usize]);
            let n=cross(sub(v,u),sub(w,u));
            if !norm(n).is_finite()||norm(n)==0.0{
                return Err(invalid("degenerate analytic cylinder wall facet"));
            }
            area+=norm(n)*0.5;
            mesh.triangles.push(tri);
            mesh.triangle_faces.push(face);
        }
    }
    Ok(area)
}
#[derive(Debug,Clone)]
pub struct ConnectedTripleDisplay {
    pub mesh:Mesh,
    pub analytic_volume:f64,
    pub facet_volume:f64,
    pub facet_area:f64,
    pub faceting_volume_error:f64,
    pub independently_polygonized_volume:f64,
    pub triangle_edges:usize,
    pub genus:usize,
    pub floor_faces:usize,
    pub segments_per_arc:usize,
}
pub struct ConnectedTripleMesher;
impl ConnectedTripleMesher {
    pub fn tessellate(body:&ConnectedTripleBrep,segments:usize)
        ->Result<ConnectedTripleDisplay,GeometryError>
    {
        if !(12..=128).contains(&segments){
            return Err(invalid("connected triple display requires 12..=128 samples per arc"));
        }
        body.validate()?;
        let m=body.model();let source=body.source_box();
        let mut mesh=Mesh{vertices:m.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new()};
        let polylines=build_edge_samples(m,&mut mesh,segments)?;
        let extent=m.bbox.max.sub(m.bbox.min);
        let eps=m.tolerance.length_at(extent.x.max(extent.y).max(extent.z));
        // Planar top/bottom and both genuine blind floors share their
        // physical edge polylines with adjacent cylinder face polygons.
        for f in &m.faces{
            match m.face_surface(f.id) {
                Some(Surface3::Plane(_)) if f.id.0==0||f.id.0==1||f.id.0>=6 => {
                    tessellate_planar(m,&mut mesh,f.id,&polylines,eps)?;
                }
                Some(Surface3::Cylinder(_))=>{
                    tessellate_cylinder(m,&mut mesh,f.id,&polylines)?;
                }
                Some(Surface3::Plane(_))=>{},
                _=>return Err(invalid("unsupported support surface in triple mesh")),
            }
        }
        // The original box exterior walls keep original face/vertex
        // ownership and need to appear exactly once.
        for (tri,face) in source.mesh.triangles.iter().zip(&source.mesh.triangle_faces) {
            if (2..=5).contains(&face.0){
                mesh.triangles.push(*tri);mesh.triangle_faces.push(*face);
            }
        }
        let mut adjacency=BTreeMap::<(u32,u32),Vec<bool>>::new();
        let mut facet_volume=0.0;
        let mut facet_area=0.0;
        let center=source.mass.centroid;
        for (tri,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
            if face.0 as usize>=m.faces.len(){
                return Err(invalid("mesh triangle references a non-existent B-rep face"));
            }
            let [a,b,c]=tri.map(|i|mesh.vertices[i as usize]);
            let (u,v,w)=(sub(a,center),sub(b,center),sub(c,center));
            facet_volume+=dot(u,cross(v,w))/6.0;
            facet_area+=norm(cross(sub(b,a),sub(c,a)))*0.5;
            for (i,j) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]{
                if i==j{return Err(invalid("display has identical edge endpoints"));}
                adjacency.entry((i.min(j),i.max(j))).or_default().push(i<j);
            }
        }
        if mesh.triangles.is_empty()||mesh.triangles.len()!=mesh.triangle_faces.len()||
            adjacency.values().any(|dir|dir.len()!=2||dir[0]==dir[1]){
            return Err(invalid("triple display is not watertight and oppositely oriented"));
        }
        let euler=mesh.vertices.len() as isize-
            adjacency.len() as isize+mesh.triangles.len() as isize;
        if euler!=0 {
            return Err(invalid("faceted connected through tunnel has wrong Euler genus"));
        }
        // Independent polygonal horizontal slices of exactly the sampled
        // exposed cylinder boundaries (not the analytic Green integrals).
        let mut polygon_removed=0.0;
        for slab in &body.plan().slabs{
            let mut planar_area=0.0;
            for f in m.faces.iter().skip(6){
                let Some(Surface3::Cylinder(_))=m.face_surface(f.id) else{continue};
                let mut circles=Vec::new();
                let mut cursor=m.loops[f.loops[0].0 as usize].first_coedge;
                loop {
                    let coedge=&m.coedges[cursor.0 as usize];
                    if let Some(Curve3::Circle(_))=m.edge_curve(coedge.edge){
                        circles.push(coedge.edge);
                    }
                    cursor=coedge.next;
                    if cursor==m.loops[f.loops[0].0 as usize].first_coedge{break;}
                }
                if circles.len()!=2{return Err(invalid("source cylinder patch lost floor/top arcs"));}
                let lo=if mesh.vertices[polylines[circles[0].0 as usize][0] as usize].y <
                    mesh.vertices[polylines[circles[1].0 as usize][0] as usize].y{
                    circles[0]
                }else{circles[1]};
                let hi=if lo==circles[0]{circles[1]}else{circles[0]};
                let lower_y=mesh.vertices[polylines[lo.0 as usize][0] as usize].y;
                let upper_y=mesh.vertices[polylines[hi.0 as usize][0] as usize].y;
                if (lower_y-slab.y_lower).abs()>eps*4.0 ||
                    (upper_y-slab.y_upper).abs()>eps*4.0{
                    continue;
                }
                let points=&polylines[lo.0 as usize];
                let anchor=source.mass.centroid;
                for v in points.windows(2){
                    let a=sub(mesh.vertices[v[0] as usize],anchor);
                    let b=sub(mesh.vertices[v[1] as usize],anchor);
                    planar_area+=(a.x*b.z-a.z*b.x)*0.5;
                }
            }
            if planar_area<=0.0{return Err(invalid("polygonized cutter union area is nonpositive"));}
            polygon_removed+=planar_area*(slab.y_upper-slab.y_lower);
        }
        let polygon_volume=source.mass.volume-polygon_removed;
        let slack=source.mass.volume.abs()*1e-8+
            eps*extent.x.max(extent.y).max(extent.z).powi(2)*64.0;
        if !facet_volume.is_finite()||facet_volume<=0.0||
            !facet_area.is_finite()||facet_area<=0.0 ||
            (facet_volume-polygon_volume).abs()>slack ||
            facet_volume<m.mass.volume-slack{
            return Err(invalid("triple faceted mass disagrees with independent polygon slices"));
        }
        let floors=m.faces.iter().skip(6).filter(|f|
            matches!(m.face_surface(f.id),Some(Surface3::Plane(_)))).count();
        if floors!=2 {
            return Err(invalid("triple display missing one of two material-facing floors"));
        }
        Ok(ConnectedTripleDisplay{
            mesh,analytic_volume:m.mass.volume,facet_volume,facet_area,
            faceting_volume_error:facet_volume-m.mass.volume,
            independently_polygonized_volume:polygon_volume,
            triangle_edges:adjacency.len(),genus:1,
            floor_faces:floors,segments_per_arc:segments
        })
    }
}
