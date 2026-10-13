//! Derived, watertight triangulation of an analytically certified genus-one
//! crossing through/blind body. Polygonization NEVER defines the B-rep or
//! its mathematical mass. Face IDs trace every facet to the canonical model.
use std::collections::BTreeMap;
use crate::{Curve3,FaceId,GeometryError,Mesh,Point3,ThroughBlindCrossingBrep};

fn invalid(msg:impl Into<String>)->GeometryError{
    GeometryError::InvalidTopology(msg.into())
}
fn norm(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn signed2(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn polygon_area(p:&[[f64;2]])->f64 {
    if p.is_empty(){return 0.0;}
    let anchor=p[0];
    (0..p.len()).map(|i|signed2(anchor,p[i],p[(i+1)%p.len()]))
        .sum::<f64>()*0.5
}
fn perimeter(p:&[[f64;2]])->f64{
    (0..p.len()).map(|i|{
        let a=p[i];let b=p[(i+1)%p.len()];
        (a[0]-b[0]).hypot(a[1]-b[1])
    }).sum()
}
fn inside(point:[f64;2],polygon:&[[f64;2]])->bool{
    let mut parity=false;
    for i in 0..polygon.len(){
        let a=polygon[i];let b=polygon[(i+1)%polygon.len()];
        if (a[1]>point[1])!=(b[1]>point[1]) &&
            point[0]<a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]){
            parity=!parity;
        }
    }
    parity
}
fn xy(p:Point3,origin:Point3)->[f64;2]{
    [p.x-origin.x,p.z-origin.z]
}
fn append_cap(mesh:&mut Mesh,outer:&[u32],inner:&[u32],
              face:FaceId,outward_positive_y:bool,
              origin:Point3,area_budget:f64)->Result<usize,GeometryError>{
    let points=outer.iter().chain(inner).map(|i|xy(mesh.vertices[*i as usize],origin))
        .collect::<Vec<_>>();
    let n=outer.len();
    let mut flat=Vec::with_capacity(points.len()*2);
    for p in &points {flat.extend_from_slice(p);}
    let triangles=earcutr::earcut(&flat,&[n],2)
        .map_err(|e|invalid(format!("cap polygon-with-hole triangulation failed: {e:?}")))?;
    if triangles.len()%3!=0||triangles.is_empty(){
        return Err(invalid("analytical cap produced no display triangles"));
    }
    let outer_area=polygon_area(&points[..n]).abs();
    let hole_area=polygon_area(&points[n..]).abs();
    if hole_area<=area_budget || outer_area<=hole_area+area_budget {
        return Err(invalid("invalid triangulated cap opening area"));
    }
    let mut area=0.0;
    for t in triangles.chunks_exact(3) {
        let mut tri=[t[0],t[1],t[2]];
        if tri.iter().any(|&id|id>=points.len()) {
            return Err(invalid("cap triangulation index escaped ring"));
        }
        let signed=signed2(points[tri[0]],points[tri[1]],points[tri[2]]);
        if !signed.is_finite()||signed.abs()<=area_budget*0.001 {
            return Err(invalid("degenerate cap triangle"));
        }
        let centroid=[
            (points[tri[0]][0]+points[tri[1]][0]+points[tri[2]][0])/3.0,
            (points[tri[0]][1]+points[tri[1]][1]+points[tri[2]][1])/3.0,
        ];
        if !inside(centroid,&points[..n]) || inside(centroid,&points[n..]) {
            return Err(invalid("cap facet violates physical material domain"));
        }
        area+=signed.abs()*0.5;
        if (signed<0.0)!=outward_positive_y {tri.swap(1,2);}
        let id=|local:usize| if local<n{outer[local]}else{inner[local-n]};
        mesh.triangles.push(tri.map(id));
        mesh.triangle_faces.push(face);
    }
    let tolerance=area_budget*(triangles.len() as f64+4.0);
    if (area-(outer_area-hole_area)).abs()>tolerance {
        return Err(invalid("triangulated cap area disagrees with original polygonal region"));
    }
    Ok(triangles.len()/3)
}
fn append_floor(mesh:&mut Mesh,ring:&[u32],origin:Point3,
                area_budget:f64)->Result<usize,GeometryError>{
    let points=ring.iter().map(|i|xy(mesh.vertices[*i as usize],origin)).collect::<Vec<_>>();
    let mut flat=Vec::with_capacity(2*points.len());
    for p in &points{flat.extend_from_slice(p);}
    let triangles=earcutr::earcut(&flat,&[],2)
        .map_err(|e|invalid(format!("crescent floor triangulation failed: {e:?}")))?;
    if triangles.len()%3!=0 || triangles.is_empty(){
        return Err(invalid("crescent floor triangulation is empty"));
    }
    let mut area=0.0;
    for tri in triangles.chunks_exact(3){
        let mut t=[tri[0],tri[1],tri[2]];
        if t.iter().any(|&v|v>=ring.len()){
            return Err(invalid("crescent triangle index escaped ring"));
        }
        let signed=signed2(points[t[0]],points[t[1]],points[t[2]]);
        if !signed.is_finite()||signed.abs()<=area_budget*0.001{
            return Err(invalid("degenerate crescent triangle"));
        }
        let centroid=[
            (points[t[0]][0]+points[t[1]][0]+points[t[2]][0])/3.0,
            (points[t[0]][1]+points[t[1]][1]+points[t[2]][1])/3.0,
        ];
        if !inside(centroid,&points){
            return Err(invalid("crescent facet is outside source planar region"));
        }
        area+=signed.abs()*0.5;
        if signed>0.0{t.swap(1,2);}
        mesh.triangles.push(t.map(|k|ring[k]));
        mesh.triangle_faces.push(FaceId(8));
    }
    let expected=polygon_area(&points).abs();
    if expected<=area_budget || (area-expected).abs()>
        area_budget*(triangles.len() as f64+4.0){
        return Err(invalid("crescent planar area certificate failed"));
    }
    Ok(triangles.len()/3)
}
fn wall(mesh:&mut Mesh,low:&[u32],high:&[u32],owners:&[FaceId])
    ->Result<(),GeometryError>{
    let n=low.len();
    if n<3||high.len()!=n||owners.len()!=n{
        return Err(invalid("cylindrical display rings have inconsistent sample counts"));
    }
    for j in 0..n{
        let k=(j+1)%n;
        mesh.triangles.extend([
            [low[j],high[j],high[k]],[low[j],high[k],low[k]]
        ]);
        mesh.triangle_faces.extend([owners[j],owners[j]]);
    }
    Ok(())
}
fn certify(body:&ThroughBlindCrossingBrep,mesh:&Mesh,
           eps:f64)->Result<(f64,f64),GeometryError>{
    let m=body.model();
    let source=body.source_box();
    let (through,blind)=body.parameters();
    let floor_y=source.bbox.max.y-blind.depth;
    if mesh.triangles.len()!=mesh.triangle_faces.len(){
        return Err(invalid("display face ownership count diverges from facets"));
    }
    let mut edges=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let origin=source.bbox.min;
    let mut six_volume=0.0;
    let mut area=0.0;
    for (triangle,face) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        let [a,b,c]=*triangle;
        if a==b||b==c||c==a||face.0 as usize>=m.faces.len()||
            triangle.iter().any(|&x|x as usize>=mesh.vertices.len()){
            return Err(invalid("display triangle has invalid indices or face owner"));
        }
        let [x,y,z]=triangle.map(|i|mesh.vertices[i as usize]);
        let normal=y.sub(x).cross(z.sub(x));
        let magnitude=norm(normal);
        if !magnitude.is_finite()||magnitude<=0.0 {
            return Err(invalid("display triangle has degenerate physical area"));
        }
        let center=x.add(y).add(z).scale(1.0/3.0);
        match face.0 {
            6|7|9=>{
                let params=if face.0==7 { (blind.center_x,blind.center_z,blind.radius) }
                    else {(through.center_x,through.center_z,through.radius)};
                let inward=Point3{x:params.0-center.x,y:0.0,z:params.1-center.z};
                if normal.dot(inward)<=0.0 ||
                    [x,y,z].iter().any(|p|
                        ((p.x-params.0).hypot(p.z-params.1)-params.2).abs()>eps){
                    return Err(invalid("cylinder facet is not inward-facing on analytic source"));
                }
            }
            0|1|8=>{
                let height=if face.0==0{origin.y}
                    else if face.0==1{source.bbox.max.y}
                    else{floor_y};
                let desired=if face.0==0 {-1.0}else{1.0};
                if normal.y*desired<=0.0 || [x,y,z].iter().any(|p|(p.y-height).abs()>eps){
                    return Err(invalid("planar cap/floor is off support or reversed"));
                }
            }
            _=>{
                let Some(crate::Surface3::Plane(surface))=m.face_surface(*face)else{
                    return Err(invalid("original box wall lacks analytic plane"));
                };
                if normal.dot(surface.frame.normal)<=0.0 ||
                   [x,y,z].iter().any(|p|
                      (p.sub(surface.frame.origin).dot(surface.frame.normal)).abs()>eps){
                    return Err(invalid("source box wall mesh is off support or reversed"));
                }
            }
        }
        for (u,v) in [(a,b),(b,c),(c,a)]{
            edges.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
        six_volume+=x.sub(origin).dot(y.sub(origin).cross(z.sub(origin)));
        area+=0.5*magnitude;
    }
    if edges.values().any(|v|v.len()!=2||v[0]==v[1]){
        return Err(invalid("derived triangular shell is not watertight and oriented"));
    }
    let volume=six_volume/6.0;
    if !volume.is_finite()||volume<=0.0||!area.is_finite()||area<=0.0{
        return Err(invalid("display signed volume or area invalid"));
    }
    Ok((volume,area))
}

#[derive(Debug,Clone)]
pub struct ThroughBlindCrossingDisplay{
    pub mesh:Mesh,
    pub segments_per_arc:usize,
    pub top_triangles:usize,
    pub bottom_triangles:usize,
    pub crescent_triangles:usize,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
}
pub struct ThroughBlindCrossingMesher;
impl ThroughBlindCrossingMesher {
    pub fn tessellate(body:&ThroughBlindCrossingBrep,segments:usize)
        ->Result<ThroughBlindCrossingDisplay,GeometryError>
    {
        if !(12..=512).contains(&segments){
            return Err(invalid("through/blind display needs 12..=512 segments per source arc"));
        }
        body.validate()?;
        let source=body.source_box();
        let m=body.model();
        let (through,blind)=body.parameters();
        if m.vertices.len()!=13||m.edges.len()!=21||m.faces.len()!=10{
            return Err(invalid("display source is not certified through/blind crossing"));
        }
        let bbox=source.bbox;
        let extent=norm(bbox.max.sub(bbox.min));
        let coord=[bbox.min,bbox.max].into_iter()
            .flat_map(|p|[p.x.abs(),p.y.abs(),p.z.abs()])
            .fold(1.0_f64,f64::max);
        let eps=source.tolerance.length_at(extent).max(32.0*f64::EPSILON*coord);
        let area_budget=source.tolerance.area_at(extent).max(eps*extent*8.0);
        let mut mesh=Mesh{
            vertices:source.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        let mut outer_bottom=Vec::with_capacity(4);
        let mut outer_top=Vec::with_capacity(4);
        for &id in &source.faces[0].boundary {
            let n=id.0 as usize;
            if n>=4{return Err(invalid("source bottom face has unexpected vertex order"));}
            outer_bottom.push(n as u32);
            outer_top.push((n+4) as u32);
        }
        let outer=outer_bottom.iter().map(|v|xy(mesh.vertices[*v as usize],bbox.min))
            .collect::<Vec<_>>();
        if polygon_area(&outer)<=0.0{
            return Err(invalid("source outer cap must have positive XZ winding"));
        }
        // Sample the two certified exterior arcs of the true union opening.
        // Lower and upper rings share the same crossing vertices with both
        // cylinder sides. The through-only lower wall later reuses the exact
        // floor ring vertices, including the 2 crossing points.
        let mut floor_union=Vec::with_capacity(2*segments);
        let mut top_union=Vec::with_capacity(2*segments);
        let mut upper_owner=Vec::with_capacity(2*segments);
        for arc_index in 0..2 {
            let bottom_edge=crate::EdgeId(12+2*arc_index);
            let top_edge=crate::EdgeId(13+2*arc_index);
            let Curve3::Circle(lower)=m.edge_curve(bottom_edge)
                .ok_or_else(||invalid("missing lower exterior circle arc"))?
                else{return Err(invalid("lower arc lost analytical circle carrier"))};
            let Curve3::Circle(upper)=m.edge_curve(top_edge)
                .ok_or_else(||invalid("missing upper exterior circle arc"))?
                else{return Err(invalid("upper arc lost analytical circle carrier"))};
            let crate::ParameterRange::Bounded{start,end}=m.geometry.edges[bottom_edge.0 as usize].domain
                else{return Err(invalid("upper arc must have bounded circle domain"))};
            if m.geometry.edges[top_edge.0 as usize].domain!=m.geometry.edges[bottom_edge.0 as usize].domain {
                return Err(invalid("upper/floor source arc domains diverge"));
            }
            for j in 0..segments{
                let u=start+(end-start)*j as f64/segments as f64;
                let (f,t)=if j==0 {
                    let lf=m.edges[bottom_edge.0 as usize].start.0 as usize;
                    let lt=m.edges[top_edge.0 as usize].start.0 as usize;
                    (m.vertices[lf].position,m.vertices[lt].position)
                }else{(lower.evaluate(u)?,upper.evaluate(u)?)};
                if !f.is_finite()||!t.is_finite()||
                    (t.y-bbox.max.y).abs()>eps ||
                    (f.y-(bbox.max.y-blind.depth)).abs()>eps ||
                    (t.x-f.x).hypot(t.z-f.z)>eps {
                    return Err(invalid("arc sample drifts away from shared vertical carrier"));
                }
                let v=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("display vertex budget exceeded"))?;
                mesh.vertices.extend([f,t]);
                floor_union.push(v);top_union.push(v+1);
                upper_owner.push(FaceId(6+arc_index as u32));
            }
        }
        let union_coords=top_union.iter().map(|id|xy(mesh.vertices[*id as usize],bbox.min))
            .collect::<Vec<_>>();
        if polygon_area(&union_coords)>=-area_budget {
            return Err(invalid("certified union opening must be clockwise in XZ"));
        }

        // Full lower through circumference: exterior arc from the union
        // contour, plus the hidden complementary through arc at blind floor.
        let mut floor_through=floor_union[..segments].to_vec();
        let mut bottom_through=Vec::with_capacity(2*segments);
        let Curve3::Circle(bottom)=m.edge_curve(crate::EdgeId(19))
            .ok_or_else(||invalid("bottom through circle unavailable"))?
            else{return Err(invalid("bottom through edge must be truly circular"))};
        let crate::ParameterRange::Bounded{start:inside_start,end:inside_end}=
            m.geometry.edges[18].domain else {
                return Err(invalid("inner crescent circle needs bounded arc"));
            };
        let crate::ParameterRange::Bounded{start:outside_start,end:outside_end}=
            m.geometry.edges[12].domain else {
                return Err(invalid("exterior through circle needs bounded arc"));
            };
        for j in 0..segments {
            let u=outside_start+(outside_end-outside_start)*j as f64/segments as f64;
            let p=bottom.evaluate(u)?;
            bottom_through.push(u32::try_from(mesh.vertices.len())
                .map_err(|_|invalid("bottom circle mesh index overflow"))?);
            mesh.vertices.push(p);
        }
        let Curve3::Circle(inside)=m.edge_curve(crate::EdgeId(18))
            .ok_or_else(||invalid("crescent inner floor arc unavailable"))?
            else{return Err(invalid("crescent floor inner edge must be circular"))};
        for j in 0..segments {
            let u=inside_start+(inside_end-inside_start)*j as f64/segments as f64;
            let floor_id=if j==0 {floor_union[segments]} else {
                let p=inside.evaluate(u)?;
                let id=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("floor crescent index overflow"))?;
                mesh.vertices.push(p);id
            };
            floor_through.push(floor_id);
            bottom_through.push(u32::try_from(mesh.vertices.len())
                .map_err(|_|invalid("bottom through ring index overflow"))?);
            mesh.vertices.push(bottom.evaluate(u)?);
        }
        let full_coords=floor_through.iter().map(|id|xy(mesh.vertices[*id as usize],bbox.min))
            .collect::<Vec<_>>();
        if polygon_area(&full_coords)>=-area_budget{
            return Err(invalid("full through bottom circle must sample one clockwise contour"));
        }

        // Real crescent floor: exterior blind arc (floor_union second half),
        // followed by the *reversed complementary through arc*. Reuse all
        // physical vertices along both upper cylinder patches and T lower wall.
        let mut crescent=floor_union[segments..].to_vec();
        crescent.push(floor_union[0]);
        for j in (1..segments).rev(){
            crescent.push(floor_through[segments+j]);
        }
        let crescent_coords=crescent.iter().map(|id|xy(mesh.vertices[*id as usize],bbox.min))
            .collect::<Vec<_>>();
        let crescent_area=polygon_area(&crescent_coords).abs();
        let through_area=polygon_area(&full_coords).abs();
        let union_area=polygon_area(&union_coords).abs();
        if crescent_area<=area_budget||
            (through_area+crescent_area-union_area).abs()>
            area_budget*(segments as f64*10.0){
            return Err(invalid("polygonal through/blind union area is inconsistent"));
        }
        let top_triangles=append_cap(&mut mesh,&outer_top,&top_union,
            FaceId(1),true,bbox.min,area_budget)?;
        let bottom_triangles=append_cap(&mut mesh,&outer_bottom,&bottom_through,
            FaceId(0),false,bbox.min,area_budget)?;
        let crescent_triangles=append_floor(&mut mesh,&crescent,bbox.min,area_budget)?;
        // Preserve box wall facets with their own preexisting analytic face IDs.
        for (triangle,face) in source.mesh.triangles.iter().zip(&source.mesh.triangle_faces){
            if face.0>=2 {
                mesh.triangles.push(*triangle);mesh.triangle_faces.push(*face);
            }
        }
        wall(&mut mesh,&floor_union,&top_union,&upper_owner)?;
        wall(&mut mesh,&bottom_through,&floor_through,
            &vec![FaceId(9);floor_through.len()])?;
        let (faceted_volume,faceted_surface_area)=certify(body,&mesh,eps)?;
        let box_height=bbox.max.y-bbox.min.y;
        let lower_height=box_height-blind.depth;
        let exposed_perimeter=perimeter(&union_coords);
        let polygon_through_perimeter=perimeter(&full_coords);
        let expected_volume=source.mass.volume-through_area*box_height-
            crescent_area*blind.depth;
        let expected_surface=source.mass.surface_area-2.0*through_area+
            polygon_through_perimeter*lower_height+
            exposed_perimeter*blind.depth;
        let volume_slack=source.mass.volume.abs()*1e-8+
            eps*extent*extent*128.0;
        let area_slack=source.mass.surface_area.abs()*1e-8+
            eps*extent*128.0;
        if (faceted_volume-expected_volume).abs()>volume_slack ||
            (faceted_surface_area-expected_surface).abs()>area_slack {
            return Err(invalid(format!(
                "closed mesh integral differs from independent polygon solid: V={faceted_volume}, expected={expected_volume}; S={faceted_surface_area}, expected={expected_surface}"
            )));
        }
        Ok(ThroughBlindCrossingDisplay{
            mesh,segments_per_arc:segments,top_triangles,bottom_triangles,
            crescent_triangles,faceted_volume,faceted_surface_area,
            analytic_volume:m.mass.volume,
            faceting_volume_error:faceted_volume-m.mass.volume,
        })
    }
}
