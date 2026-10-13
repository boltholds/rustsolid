//! Independently checked derived watertight display of a source-certified
//! connected T-A-B triple-cylinder ANALYTICAL B-rep.
//!
//! All horizontal arc samples reuse the exact same IDs on top/bottom caps,
//! depth floors and adjoining cylindrical patches. No stacked independent
//! donor meshes, duplicate caps or faceted CAD ownership are allowed.
use std::collections::{BTreeMap,BTreeSet};
use crate::{
    ConnectedTripleBrep,Curve3,EdgeId,FaceId,GeometryError,Mesh,
    Point3,TripleContour,TriplePhysicalEdgeKey,TripleStitchPlan,
};
use crate::mixed_crossing_mesh::{shoelace,triangulate_cap,wall};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn sub(a:Point3,b:Point3)->Point3{
    Point3{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}
}
fn cross(a:Point3,b:Point3)->Point3{
    Point3{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}
}
fn dot(a:Point3,b:Point3)->f64{a.x*b.x+a.y*b.y+a.z*b.z}
fn mag(a:Point3)->f64{a.x.hypot(a.y).hypot(a.z)}
fn source_circles(plan:&TripleStitchPlan)->[(f64,f64,f64);3]{
    let t=plan.input.through;
    let a=plan.input.deeper_blind;
    let b=plan.input.shallower_blind;
    [(t.center_x,t.center_z,t.radius),
     (a.center_x,a.center_z,a.radius),
     (b.center_x,b.center_z,b.radius)]
}
fn sampled_contour(contour:&TripleContour,level:usize,
    sample:&BTreeMap<(usize,usize),Vec<u32>>)
    ->Result<Vec<u32>,GeometryError>{
    let mut output=Vec::new();
    let mut starting=None::<u32>;
    let mut previous=None::<u32>;
    for part in &contour.parts {
        let reference=sample.get(&(part.arc,level))
            .ok_or_else(||invalid("3D cap or floor references missing shared arc sample"))?;
        let mut curve=reference.clone();
        if part.reversed{curve.reverse();}
        let begin=*curve.first().ok_or_else(||invalid("empty shared analytic arc"))?;
        let end=*curve.last().ok_or_else(||invalid("empty shared arc endpoint"))?;
        if let Some(last)=previous {
            if begin!=last{return Err(invalid("cap/floor arc endpoint chain is not continuous"));}
        }else {starting=Some(begin);}
        output.extend_from_slice(&curve[..curve.len()-1]);
        previous=Some(end);
    }
    if previous!=starting||output.len()<3{
        return Err(invalid("planar polygonal contour fails source-shared closure"));
    }
    if output.windows(2).any(|w|w[0]==w[1]){
        return Err(invalid("source contour contains duplicate consecutive vertices"));
    }
    Ok(output)
}
#[derive(Debug,Clone)]
pub struct ConnectedTripleDisplay {
    pub mesh:Mesh,
    pub segments_per_source_arc:usize,
    pub analytical_volume:f64,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub faceting_volume_error:f64,
    pub predicted_polygonized_volume:f64,
    pub predicted_polygonized_area:f64,
    pub floor_polygon_areas:[f64;2],
    pub upper_cap_polygon_area:f64,
    pub bottom_cap_polygon_area:f64,
    pub triangle_edges:usize,
    pub genus:usize,
}
pub struct ConnectedTripleMesher;
impl ConnectedTripleMesher {
    pub fn tessellate(body:&ConnectedTripleBrep,segments:usize)
        ->Result<ConnectedTripleDisplay,GeometryError>{
        if !(12..=256).contains(&segments) {
            return Err(invalid("connected triple display needs 12..=256 samples per source arc"));
        }
        body.validate()?;
        let m=body.model();
        let p=body.stitch();
        let source=body.source_box();
        if p.slabs.len()!=3||p.floors.len()!=2||p.walls.len()!=10 ||
            m.vertices.len()!=22||m.edges.len()!=38||m.faces.len()!=18{
            return Err(invalid("display requires certified connected triple source topology"));
        }
        let mut mesh=Mesh{
            vertices:m.vertices.iter().map(|v|v.position).collect(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        let mut mapping=BTreeMap::<(usize,usize),Vec<u32>>::new();
        let circles=source_circles(p);
        let span=source.bbox.max.sub(source.bbox.min);
        let extent=span.x.max(span.y).max(span.z);
        let eps=source.tolerance.length_at(extent);
        let area_eps=source.tolerance.area_at(extent)*0.01;
        for (i,edge) in p.shared_physical_edges.iter().enumerate(){
            let TriplePhysicalEdgeKey::CircleArc{arc,level}=edge.key else {continue};
            let entity=m.edges.get(12+i)
                .ok_or_else(||invalid("missing physical analytical circle edge"))?;
            let curve=m.edge_curve(entity.id)
                .ok_or_else(||invalid("source 3D circular carrier absent"))?;
            if !matches!(curve,Curve3::Circle(_)){
                return Err(invalid("planned circular arc has noncircular B-rep carrier"));
            }
            let atom=&p.atomic_arcs[arc];
            let (cx,cz,r)=circles[atom.source];
            let mut indices=Vec::with_capacity(segments+1);
            indices.push(entity.start.0);
            for j in 1..segments {
                let u=atom.start_angle+(atom.end_angle-atom.start_angle)*
                    (j as f64/segments as f64);
                let v=curve.evaluate(u)?;
                if !v.is_finite()||
                    ( (v.x-cx).hypot(v.z-cz)-r ).abs()>eps*8.0 ||
                    (v.y-p.levels[level]).abs()>eps*8.0 {
                    return Err(invalid("sampled source circular curve has invalid radius or height"));
                }
                let new=u32::try_from(mesh.vertices.len())
                    .map_err(|_|invalid("too many derived facet vertices"))?;
                mesh.vertices.push(v);indices.push(new);
            }
            indices.push(entity.end.0);
            if mapping.insert((arc,level),indices).is_some(){
                return Err(invalid("duplicate physical arc/height edge ID"));
            }
        }
        let outer_bottom=source.faces[0].boundary.iter()
            .map(|v|v.0).collect::<Vec<_>>();
        let outer_top=source.faces[1].boundary.iter()
            .map(|v|v.0).collect::<Vec<_>>();
        let bottom=sampled_contour(p.bottom_cap.as_ref()
            .ok_or_else(||invalid("missing bottom cap inner contour"))?,0,&mapping)?;
        let top=sampled_contour(p.top_cap.as_ref()
            .ok_or_else(||invalid("missing top cap inner contour"))?,3,&mapping)?;
        let bottom_area=triangulate_cap(
            &mut mesh,&outer_bottom,Some(&bottom),FaceId(0),-1,area_eps)?;
        let top_area=triangulate_cap(
            &mut mesh,&outer_top,Some(&top),FaceId(1),1,area_eps)?;
        let mut floor_areas=[0.0;2];
        for (i,floor) in p.floors.iter().enumerate(){
            let poly=sampled_contour(&floor.contour,i+1,&mapping)?;
            floor_areas[i]=triangulate_cap(
                &mut mesh,&poly,None,FaceId(16+i as u32),1,area_eps)?;
        }
        for (tri,face) in source.mesh.triangles.iter()
            .zip(&source.mesh.triangle_faces){
            if face.0<2{continue;}
            mesh.triangles.push(*tri);
            mesh.triangle_faces.push(*face);
        }
        let mut polygonal_wall_area=0.0;
        for (idx,w) in p.walls.iter().enumerate(){
            let lower=mapping.get(&(w.arc,w.slab))
                .ok_or_else(||invalid("source wall bottom edge lost"))?;
            let upper=mapping.get(&(w.arc,w.slab+1))
                .ok_or_else(||invalid("source wall upper edge lost"))?;
            let perimeter=wall(&mut mesh,lower,upper,
                FaceId(6+idx as u32),eps)?;
            polygonal_wall_area+=perimeter*(w.y_upper-w.y_lower);
        }
        let mut polygon_areas=Vec::new();
        for (i,layer) in p.slabs.iter().enumerate(){
            let polygon=sampled_contour(
                &layer.exposed_contour,i,&mapping)?;
            let area=shoelace(&mesh,&polygon)?.abs();
            if !area.is_finite()||area<=0.0{
                return Err(invalid("source slab polygonized union is nonphysical"));
            }
            polygon_areas.push(area);
        }
        let expected_faceted_volume=source.mass.volume-
            p.slabs.iter().zip(&polygon_areas)
                .map(|(layer,area)|area*(layer.y_upper-layer.y_lower))
                .sum::<f64>();
        let predicted_faceted_area=source.mass.surface_area-
            polygon_areas[0]-polygon_areas[2]+
            floor_areas.iter().sum::<f64>()+polygonal_wall_area;
        let mut signed_volume=0.0;
        let mut area=0.0;
        let mut edges=BTreeMap::<(u32,u32),Vec<bool>>::new();
        let origin=source.mass.centroid;
        for (triangle,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
            let [a,b,c]=*triangle;
            if a==b||b==c||c==a||owner.0 as usize>=m.faces.len() {
                return Err(invalid("degenerate triangle or invalid source face owner"));
            }
            let pts=[a,b,c].map(|v|mesh.vertices.get(v as usize).copied());
            let [Some(x),Some(y),Some(z)]=pts else{
                return Err(invalid("triangle references nonexistent shared vertex"));
            };
            let va=sub(x,origin);
            let vb=sub(y,origin);
            let vc=sub(z,origin);
            let n=cross(sub(y,x),sub(z,x));
            let facet_area=mag(n)/2.0;
            if !facet_area.is_finite()||facet_area<=area_eps{
                return Err(invalid("degenerate finite-area source facet"));
            }
            area+=facet_area;
            signed_volume+=dot(va,cross(vb,vc))/6.0;
            if owner.0>=16 && n.y<=0.0 {
                return Err(invalid("triple blind floor display has incorrect +Y normal"));
            }
            if owner.0>=6 && owner.0<16 {
                let circle=circles[p.walls[owner.0 as usize-6].source];
                let c=Point3{x:(x.x+y.x+z.x)/3.0,
                    y:(x.y+y.y+z.y)/3.0,
                    z:(x.z+y.z+z.z)/3.0};
                let toward=Point3{x:circle.0-c.x,y:0.0,z:circle.1-c.z};
                if dot(n,toward)<=0.0{
                    return Err(invalid("display cylinder patch is not inward material-facing"));
                }
            }
            for (a,b) in [(a,b),(b,c),(c,a)]{
                edges.entry((a.min(b),a.max(b))).or_default().push(a<b);
            }
        }
        if edges.values().any(|v|v.len()!=2||v[0]==v[1]){
            return Err(invalid("triple display is not a closed oriented watertight mesh"));
        }
        let chi=mesh.vertices.len() as isize-edges.len() as isize+
            mesh.triangles.len() as isize;
        if chi!=0{
            return Err(invalid("triangle shell Euler/genus differs from analytic source"));
        }
        if mesh.triangles.len()!=mesh.triangle_faces.len(){
            return Err(invalid("faceted source owner array has incorrect length"));
        }
        let slack=source.mass.volume.abs()*1e-8+
            eps*extent*extent*64.0;
        if !signed_volume.is_finite()||signed_volume<=0.0||
            !area.is_finite()||
            (signed_volume-expected_faceted_volume).abs()>slack||
            (area-predicted_faceted_area).abs()>slack*16.0 ||
            signed_volume<body.model().mass.volume-slack{
            return Err(invalid("triple faceted volume/area not independently certified"));
        }
        Ok(ConnectedTripleDisplay{
            mesh,segments_per_source_arc:segments,
            analytical_volume:m.mass.volume,
            faceted_volume:signed_volume,
            faceted_surface_area:area,
            faceting_volume_error:signed_volume-m.mass.volume,
            predicted_polygonized_volume:expected_faceted_volume,
            predicted_polygonized_area:predicted_faceted_area,
            floor_polygon_areas:floor_areas,
            upper_cap_polygon_area:top_area,
            bottom_cap_polygon_area:bottom_area,
            triangle_edges:edges.len(),genus:1,
        })
    }
}
