//! Exact-in-form, source-certified analytic B-rep for a **connected chain**
//! of one upright through cylinder and two top-entry blind cylinders.
//! The deep blind crosses the through hole; the shallow blind crosses the
//! deep blind, but remains strictly separated from the through footprint.
//!
//! v0.34's global source circle-union contours provide planar classifications
//! for three depth slabs. Real physical horizontal edges and vertical edges
//! share source vertices across each slab. Each exposed source sector owns
//! one four-trim analytic cylinder face; every depth transition contributes
//! its true planar floor crescent.
//! No generic TopologyDeltaPlan or arbitrary triple-junction Boolean is enabled.
use std::collections::BTreeMap;
use std::f64::consts::{PI,TAU};
use crate::{
    BlindHoleParameters,BrepFace,BrepModel,BrepOrigin,BrepSource,
    CircleUnionArc,ConnectedCutTransitionPlanner,Curve3,CylinderSurface,
    EdgeId,FaceGeometry,FaceId,FaceSense,Frame3,GeometryError,
    HoleParameters,LoopId,LoopRole,MassProperties,Mesh,ModelUnits,
    MultiCutPlanState,MultiCutRegionPlan,MultiCutRegionPlanner,
    ParameterDomain2,ParameterRange,PlaneSurface,Point3,Solid,Surface3,
    VertexId,
};
use crate::mixed_crossing_brep::{
    add_circle,add_surface,add_vertex,add_vertical,append_loop,bounded,twins,
};

const MAX_HISTORY:usize=128;
fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct TripleCutParameters {
    pub through:HoleParameters,
    pub deep_blind:BlindHoleParameters,
    pub shallow_blind:BlindHoleParameters,
}
#[derive(Debug,Clone,Copy)]
struct Disk{x:f64,z:f64,r:f64}
impl Disk {
    fn point(self,theta:f64,y:f64)->Point3{
        Point3{x:self.x+self.r*theta.cos(),y,z:self.z+self.r*theta.sin()}
    }
    fn angle(self,p:Point3)->f64{
        (p.z-self.z).atan2(p.x-self.x).rem_euclid(TAU)
    }
    fn gap(self,p:Point3)->f64{
        (p.x-self.x).hypot(p.z-self.z)-self.r
    }
}
fn source_disks(p:TripleCutParameters)->[Disk;3]{
    [Disk{x:p.through.center_x,z:p.through.center_z,r:p.through.radius},
     Disk{x:p.deep_blind.center_x,z:p.deep_blind.center_z,r:p.deep_blind.radius},
     Disk{x:p.shallow_blind.center_x,z:p.shallow_blind.center_z,r:p.shallow_blind.radius}]
}
fn classify_pair(a:Disk,b:Disk,eps:f64,expected_crossing:bool)
    ->Result<(),GeometryError>
{
    let d=(a.x-b.x).hypot(a.z-b.z);
    if !d.is_finite()||!a.r.is_finite()||!b.r.is_finite() ||
        (d-(a.r+b.r)).abs()<=eps*32.0 ||
        (d-(a.r-b.r).abs()).abs()<=eps*32.0 {
        return Err(invalid("unresolved tangency or nonfinite source pair"));
    }
    if expected_crossing {
        if d<= (a.r-b.r).abs()+eps*32.0 ||
            d>=a.r+b.r-eps*32.0 {
            return Err(invalid("source chain lacks strict pairwise transverse crossing"));
        }
    }else if d<=a.r+b.r+eps*32.0{
        return Err(invalid("outer source cutters must remain strictly separated"));
    }
    Ok(())
}
fn angular_cuts(plan:&MultiCutRegionPlan,disk:Disk,index:usize,eps:f64)
    ->Result<Vec<f64>,GeometryError>
{
    let mut angles=vec![0.0,PI,TAU];
    for slab in &plan.slabs{
        for event in &slab.vertices {
            if event.is_virtual_seam || (event.disk_a!=index && event.disk_b!=index) {
                continue;
            }
            let mut theta=disk.angle(event.point);
            if (theta-TAU).abs()<f64::EPSILON*32.0{theta=0.0;}
            angles.push(theta);
        }
    }
    angles.sort_by(f64::total_cmp);
    let angular_eps=(eps/disk.r).max(32.0*f64::EPSILON)*16.0;
    angles.dedup_by(|a,b|(*b-*a).abs()<=angular_eps);
    if angles[0]!=0.0 || (*angles.last().unwrap()-TAU).abs()>angular_eps ||
        angles.windows(2).any(|v|v[1]-v[0]<=angular_eps){
        return Err(invalid("source event angles lack resolved unique circular sectors"));
    }
    Ok(angles)
}
fn material_exposure(disk:usize,theta:f64,active:&[usize],
    all:&[Disk;3],eps:f64)->Result<bool,GeometryError>
{
    if !active.contains(&disk){return Ok(false);}
    let p=all[disk].point(theta,0.0);
    for &other in active {
        if other==disk{continue;}
        let gap=all[other].gap(p);
        if !gap.is_finite()||gap.abs()<=eps*8.0{
            return Err(invalid("boundary sample too close to another cylinder"));
        }
        if gap<0.0{return Ok(false);}
    }
    Ok(true)
}
fn area_moment(arc:&CircleUnionArc,disk:Disk,anchor:Point3)->(f64,f64){
    // Independent analytic first moments of the disk-union area in XZ,
    // using 1/2 ∮ x^2 dz and -1/2 ∮ z^2 dx with relative coordinates.
    let a=arc.start_angle;
    let b=arc.end_angle;
    let delta=b-a;
    let (sa,sb,ca,cb)=(a.sin(),b.sin(),a.cos(),b.cos());
    let sin2= (2.0*b).sin()-(2.0*a).sin();
    let cos_sq=delta*0.5+sin2*0.25;
    let sin_sq=delta*0.5-sin2*0.25;
    let cos1=sb-sa;
    let sin1=ca-cb;
    let cos3=cos1-(sb.powi(3)-sa.powi(3))/3.0;
    let sin3=sin1+(cb.powi(3)-ca.powi(3))/3.0;
    let x=disk.x-anchor.x;let z=disk.z-anchor.z;
    let mx=0.5*disk.r*(x*x*cos1+
        2.0*x*disk.r*cos_sq+disk.r*disk.r*cos3);
    let mz=0.5*disk.r*(z*z*sin1+
        2.0*z*disk.r*sin_sq+disk.r*disk.r*sin3);
    (mx,mz)
}
fn source_mass(source:&Solid,plan:&MultiCutRegionPlan,disks:&[Disk;3])
    ->Result<MassProperties,GeometryError>
{
    let reference=source.mass.centroid;
    let mut removed_x=0.0;let mut removed_y=0.0;let mut removed_z=0.0;
    let mut removed=0.0;
    for slab in &plan.slabs{
        let h=slab.y_upper-slab.y_lower;
        let (mut mx,mut mz)=(0.0,0.0);
        for lp in &slab.boundary_loops{
            for arc in &lp.arcs{
                let (x,z)=area_moment(arc,disks[arc.disk_index],reference);
                mx+=x;mz+=z;
            }
        }
        removed+=slab.union_area*h;
        removed_x+=mx*h;removed_z+=mz*h;
        removed_y+=slab.union_area*h*
            ((slab.y_lower+slab.y_upper)*0.5-reference.y);
    }
    let v=source.mass.volume-removed;
    if !v.is_finite()||v<=0.0{return Err(invalid("invalid global analytic material volume"));}
    let expected=plan.predicted_remaining_volume
        .ok_or_else(||invalid("missing independent multilayer volume"))?;
    if (v-expected).abs()>source.tolerance.length_at(1.0)*128.0 {
        return Err(invalid("independent depth-layer mass does not match source integral"));
    }
    Ok(MassProperties{
        volume:v,surface_area:0.0,
        centroid:Point3{x:reference.x-removed_x/v,
            y:reference.y-removed_y/v,z:reference.z-removed_z/v},
    })
}
fn face_id(model:&BrepModel)->Result<FaceId,GeometryError>{
    Ok(FaceId(u32::try_from(model.faces.len())
        .map_err(|_|invalid("face index overflow"))?))
}
fn vertex_at(model:&mut BrepModel,point:Point3,eps:f64)->Result<VertexId,GeometryError>{
    // Real pair intersections on different circles at the same layer must
    // resolve to ONE physical topological vertex. Distinct near-coincident
    // support vertices would violate the source angular-resolution gate.
    for v in model.vertices.iter().skip(8){
        let p=v.position;
        let separation=(p.x-point.x).hypot(p.y-point.y).hypot(p.z-point.z);
        if separation<=eps*8.0{return Ok(v.id);}
    }
    add_vertex(model,point)
}
fn ordered_parts(model:&BrepModel,parts:&[(EdgeId,bool)])
    ->Result<Vec<(EdgeId,bool)>,GeometryError>
{
    if parts.len()<2 {return Err(invalid("planar face needs at least two real edges"));}
    let mut start=BTreeMap::<u32,usize>::new();
    let mut incoming=BTreeMap::<u32,usize>::new();
    for (i,&(id,reversed)) in parts.iter().enumerate(){
        let e=&model.edges[id.0 as usize];
        let v0=if reversed{e.end.0}else{e.start.0};
        let v1=if reversed{e.start.0}else{e.end.0};
        if start.insert(v0,i).is_some(){return Err(invalid("face boundary branches at source vertex"));}
        *incoming.entry(v1).or_default()+=1;
    }
    if start.len()!=parts.len()||incoming.len()!=parts.len()||
       incoming.values().any(|&n|n!=1){
        return Err(invalid("face boundary does not have one predecessor and successor"));
    }
    let mut ordered=Vec::with_capacity(parts.len());
    let mut seen=vec![false;parts.len()];
    let first=*start.keys().next().unwrap();
    let mut at=first;
    loop {
        let i=*start.get(&at).ok_or_else(||invalid("open analytic face boundary"))?;
        if seen[i]{return Err(invalid("prematurely repeated face edge"));}
        seen[i]=true;
        let (id,reversed)=parts[i];
        ordered.push((id,reversed));
        let e=&model.edges[id.0 as usize];
        at=if reversed{e.start.0}else{e.end.0};
        if at==first{break;}
    }
    if ordered.len()!=parts.len(){return Err(invalid("face boundary contains multiple disconnected cycles"));}
    Ok(ordered)
}

#[derive(Debug,Clone)]
pub struct ConnectedTripleBrep {
    source:Solid,
    parameters:TripleCutParameters,
    model:BrepModel,
    plan:MultiCutRegionPlan,
    revision:u64,
    undo:Vec<TripleCutParameters>,
    redo:Vec<TripleCutParameters>,
}
impl ConnectedTripleBrep {
    pub fn from_box(source:&Solid,parameters:TripleCutParameters)
        ->Result<Self,GeometryError>{
        let (model,plan)=Self::build(source,parameters)?;
        Ok(Self{source:source.clone(),parameters,model,plan,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn parameters(&self)->TripleCutParameters{self.parameters}
    pub fn plan(&self)->&MultiCutRegionPlan{&self.plan}
    pub(crate) fn source_box(&self)->&Solid{&self.source}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::ConnectedTripleDisplay,GeometryError>{
        crate::ConnectedTripleMesher::tessellate(self,segments)
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (expected,plan)=Self::build(&self.source,self.parameters)?;
        let a=&self.model;let b=&expected;
        if a.summary()!=b.summary() || a.origin!=b.origin ||
            a.mass!=b.mass || a.bbox!=b.bbox ||
            a.tolerance!=b.tolerance || a.units!=b.units ||
            a.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>()!=
                b.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>() ||
            a.edges.iter().zip(&b.edges).any(|(x,y)|
                x.id!=y.id||x.start!=y.start||x.end!=y.end||x.coedges!=y.coedges) ||
            a.coedges.iter().zip(&b.coedges).any(|(x,y)|
                x.id!=y.id||x.edge!=y.edge||x.face!=y.face||
                x.loop_id!=y.loop_id||x.reversed!=y.reversed||
                x.next!=y.next||x.prev!=y.prev||x.twin!=y.twin) ||
            a.loops.iter().zip(&b.loops).any(|(x,y)|
                x.id!=y.id||x.face!=y.face||x.role!=y.role||
                x.first_coedge!=y.first_coedge) ||
            a.faces.iter().zip(&b.faces).any(|(x,y)|
                x.id!=y.id||x.loops!=y.loops||x.sense!=y.sense||
                x.geometry!=y.geometry) ||
            a.shells.iter().zip(&b.shells).any(|(x,y)|
                x.id!=y.id||x.closed!=y.closed||x.faces!=y.faces) ||
            a.geometry.edges!=b.geometry.edges ||
            a.geometry.coedges!=b.geometry.coedges ||
            a.geometry.curves2!=b.geometry.curves2 ||
            a.geometry.curves3!=b.geometry.curves3 ||
            a.geometry.surfaces!=b.geometry.surfaces ||
            self.plan!=plan{
            return Err(invalid("triple B-rep diverges from independently rebuilt source"));
        }
        Ok(())
    }
    fn build(source:&Solid,p:TripleCutParameters)
        ->Result<(BrepModel,MultiCutRegionPlan),GeometryError>
    {
        crate::boolean::certify_box(source)
            .map_err(|e|invalid(format!("triple source is not a six-face box: {e}")))?;
        let through=[p.through];
        let blind=[p.deep_blind,p.shallow_blind];
        let plan=MultiCutRegionPlanner::plan(source,&through,&blind)?;
        if plan.state!=MultiCutPlanState::Ready ||
            plan.slabs.len()!=3 ||
            plan.slabs[0].active_disks!=[0] ||
            plan.slabs[1].active_disks!=[0,1] ||
            plan.slabs[2].active_disks!=[0,1,2] ||
            plan.slabs.iter().any(|s|s.connected_union_components!=1||s.planar_holes!=0) {
            return Err(invalid("triple chain requires three connected, hole-free source slabs"));
        }
        plan.validate_against(source,&through,&blind)?;
        let transitions=ConnectedCutTransitionPlanner::plan(source,&through,&blind)?;
        transitions.validate_against(source,&through,&blind)?;
        if transitions.transitions.len()!=2{
            return Err(invalid("connected triple has incorrect number of physical floors"));
        }
        let disks=source_disks(p);
        let eps=plan.length_budget;
        classify_pair(disks[0],disks[1],eps,true)?;
        classify_pair(disks[1],disks[2],eps,true)?;
        classify_pair(disks[0],disks[2],eps,false)?;
        if !(p.deep_blind.depth>p.shallow_blind.depth+eps*32.0){
            return Err(invalid("first interacting blind must reach strictly deeper than second"));
        }
        let cuts=(0..3).map(|i|angular_cuts(&plan,disks[i],i,eps))
            .collect::<Result<Vec<_>,_>>()?;
        let ys=[plan.slabs[0].y_lower,plan.slabs[0].y_upper,
                plan.slabs[1].y_upper,plan.slabs[2].y_upper];
        let mut exposed=Vec::new();
        for slab in &plan.slabs{
            let mut active=Vec::new();
            for disk in 0..3{
                let mut flags=Vec::new();
                for ab in cuts[disk].windows(2){
                    let mid=ab[0]+(ab[1]-ab[0])*0.5;
                    flags.push(material_exposure(disk,mid,&slab.active_disks,&disks,eps)?);
                }
                active.push(flags);
            }
            exposed.push(active);
        }
        let mut model=BrepModel::from_polyhedral(source,ModelUnits::default())?;
        if model.vertices.len()!=8 || model.edges.len()!=12 ||
            model.faces.len()!=6 || model.coedges.len()!=24 ||
            model.loops.len()!=6 || model.shells.len()!=1 {
            return Err(invalid("triple source has noncanonical box graph"));
        }
        let x=Point3{x:1.0,y:0.0,z:0.0};
        let z=Point3{x:0.0,y:0.0,z:1.0};
        let frames=(0..3).map(|i|Frame3::from_axes(
            Point3{x:disks[i].x,y:ys[3],z:disks[i].z},x,z,model.tolerance
        )).collect::<Result<Vec<_>,_>>()?;
        let mut horizontal=BTreeMap::<(usize,usize,usize),EdgeId>::new();
        for layer in 0..ys.len(){
            for disk in 0..3{
                for seg in 0..cuts[disk].len()-1{
                    let belongs=(layer>0 && exposed[layer-1][disk][seg]) ||
                        (layer<plan.slabs.len() && exposed[layer][disk][seg]);
                    if !belongs{continue;}
                    let v0=vertex_at(&mut model,disks[disk].point(cuts[disk][seg],ys[layer]),eps)?;
                    let v1=vertex_at(&mut model,disks[disk].point(cuts[disk][seg+1],ys[layer]),eps)?;
                    if v0==v1{return Err(invalid("source arc endpoints merged within tolerance"));}
                    let edge=add_circle(&mut model,frames[disk],disks[disk].r,ys[layer],
                        (cuts[disk][seg],cuts[disk][seg+1]),v0,v1)?;
                    horizontal.insert((layer,disk,seg),edge);
                }
            }
        }
        let mut vertical=BTreeMap::<(usize,u32,u32),EdgeId>::new();
        let mut cylinder_area=0.0;
        for slab in 0..plan.slabs.len(){
            for disk in 0..3{
                for seg in 0..cuts[disk].len()-1 {
                    if !exposed[slab][disk][seg]{continue;}
                    let lower=*horizontal.get(&(slab,disk,seg))
                        .ok_or_else(||invalid("missing lower physical arc"))?;
                    let upper=*horizontal.get(&(slab+1,disk,seg))
                        .ok_or_else(||invalid("missing upper physical arc"))?;
                    let (a,b)=(cuts[disk][seg],cuts[disk][seg+1]);
                    let (vl0,vl1)=(model.edges[lower.0 as usize].start,
                        model.edges[lower.0 as usize].end);
                    let (vu0,vu1)=(model.edges[upper.0 as usize].start,
                        model.edges[upper.0 as usize].end);
                    let v_edge=|model:&mut BrepModel,cache:&mut BTreeMap<(usize,u32,u32),EdgeId>,
                        vlow:VertexId,vhigh:VertexId|->Result<EdgeId,GeometryError>{
                        let key=(slab,vlow.0,vhigh.0);
                        if let Some(&id)=cache.get(&key){return Ok(id);}
                        let e=add_vertical(model,vlow,vhigh)?;
                        cache.insert(key,e);
                        Ok(e)
                    };
                    let start=v_edge(&mut model,&mut vertical,vl0,vu0)?;
                    let end=v_edge(&mut model,&mut vertical,vl1,vu1)?;
                    if start==end{return Err(invalid("coincident vertical seam endpoints"));}
                    let face=face_id(&model)?;
                    let sid=add_surface(&mut model,Surface3::Cylinder(
                        CylinderSurface::new(frames[disk],disks[disk].r,model.tolerance)?))?;
                    let uv=ParameterDomain2{
                        u:bounded(a,b)?,
                        v:bounded(ys[3]-ys[slab+1],ys[3]-ys[slab])?,
                        u_periodic:false,v_periodic:false
                    };
                    model.faces.push(BrepFace{id:face,loops:Vec::new(),
                        sense:FaceSense::Reversed,geometry:FaceGeometry{
                            surface:sid,bounds:uv,outer_loop:LoopId(u32::MAX),
                            tolerance:model.tolerance
                        }});
                    model.shells[0].faces.push(face);
                    append_loop(&mut model,face,LoopRole::Outer,
                        &[(upper,true),(start,true),(lower,false),(end,false)])?;
                    cylinder_area+=disks[disk].r*(b-a)*(ys[slab+1]-ys[slab]);
                }
            }
        }
        let bottom=horizontal.iter().filter(|((layer,disk,seg),_)|
            *layer==0&&exposed[0][*disk][*seg])
            .map(|(_,id)|(*id,true)).collect::<Vec<_>>();
        let top=horizontal.iter().filter(|((layer,disk,seg),_)|
            *layer==3&&exposed[2][*disk][*seg])
            .map(|(_,id)|(*id,false)).collect::<Vec<_>>();
        append_loop(&mut model,FaceId(0),LoopRole::Inner,
            &ordered_parts(&model,&bottom)?)?;
        append_loop(&mut model,FaceId(1),LoopRole::Inner,
            &ordered_parts(&model,&top)?)?;
        let mut floor_area=0.0;
        for layer in 1..3{
            let mut paths=Vec::new();
            for disk in 0..3{
                for seg in 0..cuts[disk].len()-1{
                    let below=exposed[layer-1][disk][seg];
                    let above=exposed[layer][disk][seg];
                    if below==above{continue;}
                    let edge=*horizontal.get(&(layer,disk,seg))
                        .ok_or_else(||invalid("floor requires a missing physical edge"))?;
                    // Floor uses the opposite orientation to its bounding
                    // cylinder wall at the same physical elevation.
                    paths.push((edge,above));
                }
            }
            let ordered=ordered_parts(&model,&paths)?;
            let plane=PlaneSurface{frame:Frame3::from_axes(
                Point3{x:source.bbox.min.x,y:ys[layer],z:source.bbox.min.z},
                x,z,model.tolerance)?};
            let surface=add_surface(&mut model,Surface3::Plane(plane))?;
            let face=face_id(&model)?;
            model.faces.push(BrepFace{id:face,loops:Vec::new(),
                sense:FaceSense::Reversed,geometry:FaceGeometry{
                    surface,bounds:ParameterDomain2{
                        u:ParameterRange::Unbounded,v:ParameterRange::Unbounded,
                        u_periodic:false,v_periodic:false},
                    outer_loop:LoopId(u32::MAX),tolerance:model.tolerance
                }});
            model.shells[0].faces.push(face);
            append_loop(&mut model,face,LoopRole::Outer,&ordered)?;
            let area=plan.slabs[layer].union_area-plan.slabs[layer-1].union_area;
            if !area.is_finite()||area<=eps*eps*64.0{
                return Err(invalid("blind floor has unresolved material area"));
            }
            floor_area+=area;
        }
        twins(&mut model,12)?;
        let mut mass=source_mass(source,&plan,&disks)?;
        mass.surface_area=source.mass.surface_area-
            plan.slabs[0].union_area-plan.slabs[2].union_area+
            floor_area+cylinder_area;
        if !mass.surface_area.is_finite()||mass.surface_area<=0.0{
            return Err(invalid("triple analytic surface area is invalid"));
        }
        model.mass=mass;
        model.origin=BrepOrigin::AnalyticCrossingThroughBlind;
        model.validate()?;
        let stats=model.summary();
        if stats.shells!=1||stats.genus!=1||stats.euler_characteristic!=0||
            model.faces[0].loops.len()!=2||
            model.faces[1].loops.len()!=2||
            model.faces.iter().filter(|f|f.id.0>=6 &&
                matches!(model.face_surface(f.id),Some(Surface3::Plane(_)))).count()!=2 ||
            model.edges.len()-12!=(model.coedges.len()-24)/2 {
            return Err(invalid("connected triple physical B-rep does not close as genus one"));
        }
        Ok((model,plan))
    }
    pub fn begin_edit(&mut self)->ConnectedTripleTransaction<'_>{
        let staged=self.parameters;
        ConnectedTripleTransaction{owner:self,staged,poisoned:false,changes:0}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let Some(next)= (if redo{&self.redo}else{&self.undo}).last()else{return Ok(false)};
        let (model,plan)=Self::build(&self.source,*next)?;
        let revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("connected triple source revision overflow"))?;
        let before=self.parameters;
        let candidate=*next;
        if redo{self.redo.pop();self.undo.push(before);}
        else{self.undo.pop();self.redo.push(before);}
        self.parameters=candidate;self.model=model;self.plan=plan;
        self.revision=revision;
        Ok(true)
    }
}
impl BrepSource for ConnectedTripleBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct ConnectedTripleTransaction<'a>{
    owner:&'a mut ConnectedTripleBrep,
    staged:TripleCutParameters,
    poisoned:bool,
    changes:usize,
}
impl ConnectedTripleTransaction<'_>{
    pub fn replace_through(&mut self,value:HoleParameters)->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("triple transaction poisoned"));}
        if !value.radius.is_finite()||value.radius<=0.0||
           !value.center_x.is_finite()||!value.center_z.is_finite(){
            self.poisoned=true;return Err(invalid("invalid staged through cutter"));
        }
        self.staged.through=value;self.changes+=1;Ok(())
    }
    pub fn replace_blind(&mut self,index:usize,value:BlindHoleParameters)
        ->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("triple transaction poisoned"));}
        if index>1||!value.radius.is_finite()||value.radius<=0.0||
           !value.depth.is_finite()||value.depth<=0.0||
           !value.center_x.is_finite()||!value.center_z.is_finite(){
            self.poisoned=true;return Err(invalid("invalid staged blind cutter/index"));
        }
        if index==0{self.staged.deep_blind=value;}
        else{self.staged.shallow_blind=value;}
        self.changes+=1;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned||self.changes==0{return Err(invalid("poisoned/empty triple transaction"));}
        let (model,plan)=ConnectedTripleBrep::build(&self.owner.source,self.staged)?;
        let revision=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("connected triple revision overflow"))?;
        if self.owner.undo.len()>=MAX_HISTORY{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters);
        self.owner.redo.clear();
        self.owner.model=model;self.owner.plan=plan;
        self.owner.parameters=self.staged;self.owner.revision=revision;
        Ok(revision)
    }
}
