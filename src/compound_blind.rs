//! Canonical composition of multiple independent analytical blind-bore
//! regions. Each connected aperture component is either one full circular
//! blind bore or exactly two partially overlapping, equal-depth bores.
//! This is a bounded executable Boolean arrangement, not a general CSG
//! executor. The source plane/cylinder witnesses remain analytical and
//! each candidate rebuild validates the whole manifold B-rep atomically.
use crate::{
    BlindHoleBrep,BlindHoleParameters,BrepFace,BrepGeometry,BrepModel,BrepOrigin,
    CoedgeId,Curve2Id,Curve3Id,EdgeId,FaceId,GeometryError,LoopId,LoopRole,
    MassProperties,ModelUnits,Point3,Solid,Surface3,SurfaceId,VertexId,
    IntersectingBlindBrep,
};
use std::f64::consts::{PI,TAU};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn idx(i:usize)->Result<u32,GeometryError>{
    u32::try_from(i).map_err(|_|invalid("blind region graph exceeds topological ID capacity"))
}
fn map_id(id:u32,base:u32,first:u32)->Result<u32,GeometryError>{
    if id<base {Ok(id)} else{
        first.checked_add(id-base).ok_or_else(||invalid("analytical region ID overflow"))
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum BlindRegionKind { SingleCircle, TwoCrossingCircles }
#[derive(Debug,Clone,PartialEq)]
pub struct BlindRegionComponent {
    pub kind:BlindRegionKind,
    pub source_indices:Vec<usize>,
    pub opening_loop:LoopId,
    pub floor_face:FaceId,
    pub cylinder_faces:Vec<FaceId>,
    pub first_edge:EdgeId,
    pub region_area:f64,
    pub exposed_perimeter:f64,
    pub depth:f64,
}
impl BlindRegionComponent {
    pub fn circle_count(&self)->usize{self.source_indices.len()}
    pub fn triangle_ring_count(&self,segments_per_circle:usize)->usize{
        self.circle_count()*segments_per_circle
    }
}

/// Geometry-only compositor. The published B-rep, including UV trim
/// ownership and floor/source lineage, is independent of any display mesh.
#[derive(Debug,Clone,Copy,Default)]
pub struct BlindRegionComposer;
impl BlindRegionComposer {
    pub fn build(source:&Solid,holes:&[BlindHoleParameters])
        ->Result<(BrepModel,Vec<BlindRegionComponent>),GeometryError>
    {
        crate::boolean::certify_box(source).map_err(|e|invalid(e.to_string()))?;
        let plan=plan_groups(source,holes)?;
        let mut model=BrepModel::from_polyhedral(source,ModelUnits::default())?;
        let mut components=Vec::with_capacity(plan.len());
        let mut removed_volume=0.0;
        let mut removed_relative_moment=Point3{x:0.0,y:0.0,z:0.0};
        let mut exposed_area=source.mass.surface_area;
        let old_c=source.mass.centroid;
        for group in plan {
            let (donor,area,perimeter,depth,relative_area_moment)=
                if group.len()==1 {
                    let p=holes[group[0]];
                    let shape=BlindHoleBrep::from_box(source,p)?;
                    shape.opening_trace()?;
                    shape.validate()?;
                    let area=PI*p.radius*p.radius;
                    (shape.model().clone(),area,TAU*p.radius,p.depth,
                        [(p.center_x-old_c.x)*area,
                         (p.center_z-old_c.z)*area])
                }else{
                    let parameters=[holes[group[0]],holes[group[1]]];
                    let shape=IntersectingBlindBrep::from_box(source,parameters)?;
                    shape.validate()?;
                    let region=shape.region();
                    (shape.model().clone(),
                        region.disk_union_area,region.exposed_perimeter,
                        region.cut_depth,region.relative_area_moment)
                };
            let kind=if group.len()==1 {
                BlindRegionKind::SingleCircle
            }else{BlindRegionKind::TwoCrossingCircles};
            let bound=graft_region(&mut model,&donor,&group,kind,
                area,perimeter,depth)?;
            components.push(bound);
            let dv=area*depth;
            removed_volume+=dv;
            removed_relative_moment=removed_relative_moment.add(Point3{
                x:relative_area_moment[0]*depth,
                y:dv*(source.bbox.max.y-depth*0.5-old_c.y),
                z:relative_area_moment[1]*depth,
            });
            exposed_area+=perimeter*depth;
        }
        let remaining=source.mass.volume-removed_volume;
        if !remaining.is_finite()||remaining<=0.0 ||
            !exposed_area.is_finite()||exposed_area<=0.0 {
            return Err(invalid("combined blind aperture has nonpositive analytic volume/area"));
        }
        model.mass=MassProperties{
            volume:remaining,surface_area:exposed_area,
            centroid:old_c.sub(removed_relative_moment.scale(1.0/remaining))
        };
        model.origin=BrepOrigin::AnalyticCompoundBlindRegions;
        model.validate()?;
        let expected_v=8+components.iter().map(|c|
            if c.kind==BlindRegionKind::SingleCircle{2}else{4}).sum::<usize>();
        let expected_e=12+components.iter().map(|c|
            if c.kind==BlindRegionKind::SingleCircle{3}else{6}).sum::<usize>();
        let expected_f=6+components.iter().map(|c|
            if c.kind==BlindRegionKind::SingleCircle{2}else{3}).sum::<usize>();
        let summary=model.summary();
        if summary.vertices!=expected_v || summary.edges!=expected_e ||
            summary.coedges!=2*expected_e ||
            summary.faces!=expected_f ||
            summary.loops!=6+components.iter().map(|c|
                if c.kind==BlindRegionKind::SingleCircle{3}else{4}).sum::<usize>() ||
            summary.shells!=1 || summary.euler_characteristic!=2 ||
            summary.genus!=0 || model.faces[1].loops.len()!=1+components.len() ||
            model.loops.iter().filter(|l|l.role==LoopRole::Inner).count()!=components.len() ||
            model.faces[0].loops.len()!=1 {
            return Err(invalid("compound blind B-rep has inconsistent shell genus or region ownership"));
        }
        for c in &components {
            if model.face_sense(c.floor_face)!=Some(crate::FaceSense::Reversed)
                || c.cylinder_faces.iter().any(|&f|
                    model.face_sense(f)!=Some(crate::FaceSense::Reversed)) {
                return Err(invalid("compound analytical cut floor/wall material sense mismatch"));
            }
        }
        Ok((model,components))
    }
}
fn plan_groups(source:&Solid,holes:&[BlindHoleParameters])
    ->Result<Vec<Vec<usize>>,GeometryError>
{
    if holes.is_empty()||holes.len()>16 {
        return Err(invalid("compound blind builder accepts 1..=16 top-entry cylinders"));
    }
    let b=source.bbox;
    let scale=(b.max.x-b.min.x).max(b.max.y-b.min.y).max(b.max.z-b.min.z);
    let epsilon=source.tolerance.length_at(scale);
    // Single-cut source certificates cannot be waived by a group that
    // incidentally produces a plausible circle-circle crossing.
    for &h in holes {
        let shape=BlindHoleBrep::from_box(source,h)?;
        shape.opening_trace()?;
    }
    let mut partner=vec![None;holes.len()];
    for i in 0..holes.len(){
        for j in i+1..holes.len(){
            let a=holes[i];let c=holes[j];
            let d=(a.center_x-c.center_x).hypot(a.center_z-c.center_z);
            let radius_sum=a.radius+c.radius;
            let radius_diff=(a.radius-c.radius).abs();
            if !d.is_finite() {
                return Err(invalid("nonfinite distance between blind circle centers"));
            }
            if d>radius_sum+epsilon {
                continue;
            }
            if d>=radius_sum-epsilon {
                return Err(invalid("tangent or unresolved near-tangent independent apertures"));
            }
            if d<=radius_diff+epsilon {
                return Err(invalid("contained, coincident or internally tangent cutting circles"));
            }
            if partner[i].is_some()||partner[j].is_some() {
                return Err(invalid("compound regions with three-way or chained circle crossings are not implemented"));
            }
            let slack=16.0*f64::EPSILON*
                [a.depth.abs(),c.depth.abs(),b.max.y.abs()].into_iter()
                .fold(1.0_f64,f64::max);
            if (a.depth-c.depth).abs()>slack {
                return Err(invalid("two overlapping bores must share exactly one planar floor level"));
            }
            partner[i]=Some(j);
            partner[j]=Some(i);
        }
    }
    let mut visited=vec![false;holes.len()];
    let mut groups=Vec::new();
    for i in 0..holes.len(){
        if visited[i]{continue;}
        visited[i]=true;
        if let Some(j)=partner[i]{
            if visited[j]||i>=j||partner[j]!=Some(i){
                return Err(invalid("non-bijective overlapping circle component assignment"));
            }
            visited[j]=true;
            groups.push(vec![i,j]);
        }else{groups.push(vec![i]);}
    }
    Ok(groups)
}
fn graft_region(model:&mut BrepModel,donor:&BrepModel,
    source_indices:&[usize],kind:BlindRegionKind,
    area:f64,perimeter:f64,depth:f64
)->Result<BlindRegionComponent,GeometryError>
{
    donor.validate()?;
    let base_v=8usize;let base_e=12usize;let base_c=24usize;
    let base_l=6usize;let base_f=6usize;
    if donor.vertices.len()<=base_v||donor.edges.len()<=base_e ||
       donor.coedges.len()<=base_c||donor.faces.len()<=base_f ||
       donor.loops.len()<=base_l ||
       donor.faces[0].loops.len()!=1 || donor.faces[1].loops.len()!=2 {
        return Err(invalid("donor is not an independently source-validated blind aperture"));
    }
    let target_v=idx(model.vertices.len())?;
    let target_e=idx(model.edges.len())?;
    let target_c=idx(model.coedges.len())?;
    let target_l=idx(model.loops.len())?;
    let target_f=idx(model.faces.len())?;
    let target_cr3=idx(model.geometry.curves3.len())?;
    let target_cr2=idx(model.geometry.curves2.len())?;
    let target_s=idx(model.geometry.surfaces.len())?;
    let base_cr3=idx(base_e)?;let base_cr2=idx(base_c)?;
    let base_s=idx(base_f)?;
    if donor.geometry.curves3.len()<base_e ||
       donor.geometry.curves2.len()<base_c ||
       donor.geometry.surfaces.len()<base_f {
        return Err(invalid("donor analytic carrier catalog is truncated"));
    }
    // Verify shared uncut block catalog instead of blindly attaching
    // locally parameterized edges to a different source box.
    if donor.geometry.curves3[..base_e]!=model.geometry.curves3[..base_e] ||
       donor.geometry.curves2[..base_c]!=model.geometry.curves2[..base_c] ||
       donor.geometry.surfaces[..base_f]!=model.geometry.surfaces[..base_f] {
        return Err(invalid("donor box analytical carrier prefix differs from shared source"));
    }
    model.geometry.curves3.extend_from_slice(&donor.geometry.curves3[base_e..]);
    model.geometry.curves2.extend_from_slice(&donor.geometry.curves2[base_c..]);
    model.geometry.surfaces.extend_from_slice(&donor.geometry.surfaces[base_f..]);
    for old in &donor.vertices[base_v..]{
        let mut next=old.clone();
        next.id=VertexId(map_id(old.id.0,base_v as u32,target_v)?);
        model.vertices.push(next);
    }
    for (old,g) in donor.edges[base_e..].iter()
        .zip(&donor.geometry.edges[base_e..]){
        let mut next=old.clone();
        next.id=EdgeId(map_id(old.id.0,base_e as u32,target_e)?);
        next.start=VertexId(map_id(old.start.0,base_v as u32,target_v)?);
        next.end=VertexId(map_id(old.end.0,base_v as u32,target_v)?);
        next.coedges=[
            CoedgeId(map_id(old.coedges[0].0,base_c as u32,target_c)?),
            CoedgeId(map_id(old.coedges[1].0,base_c as u32,target_c)?),
        ];
        model.edges.push(next);
        let mut geometry=*g;
        geometry.curve=Curve3Id(map_id(g.curve.0,base_cr3,target_cr3)?);
        model.geometry.edges.push(geometry);
    }
    for (old,g) in donor.coedges[base_c..].iter()
        .zip(&donor.geometry.coedges[base_c..]){
        let mut next=old.clone();
        next.id=CoedgeId(map_id(old.id.0,base_c as u32,target_c)?);
        next.edge=EdgeId(map_id(old.edge.0,base_e as u32,target_e)?);
        next.face=FaceId(map_id(old.face.0,base_f as u32,target_f)?);
        next.loop_id=LoopId(map_id(old.loop_id.0,base_l as u32,target_l)?);
        next.next=CoedgeId(map_id(old.next.0,base_c as u32,target_c)?);
        next.prev=CoedgeId(map_id(old.prev.0,base_c as u32,target_c)?);
        next.twin=CoedgeId(map_id(old.twin.0,base_c as u32,target_c)?);
        model.coedges.push(next);
        let mut geometry=*g;
        geometry.pcurve=Curve2Id(map_id(g.pcurve.0,base_cr2,target_cr2)?);
        geometry.face=FaceId(map_id(g.face.0,base_f as u32,target_f)?);
        geometry.loop_id=LoopId(map_id(g.loop_id.0,base_l as u32,target_l)?);
        model.geometry.coedges.push(geometry);
    }
    for old in &donor.loops[base_l..]{
        let mut next=old.clone();
        next.id=LoopId(map_id(old.id.0,base_l as u32,target_l)?);
        next.face=FaceId(map_id(old.face.0,base_f as u32,target_f)?);
        next.first_coedge=CoedgeId(map_id(
            old.first_coedge.0,base_c as u32,target_c)?);
        model.loops.push(next);
    }
    for old in &donor.faces[base_f..]{
        let mut next=old.clone();
        next.id=FaceId(map_id(old.id.0,base_f as u32,target_f)?);
        next.loops=old.loops.iter().map(|id|
            map_id(id.0,base_l as u32,target_l).map(LoopId))
            .collect::<Result<Vec<_>,_>>()?;
        next.geometry.surface=SurfaceId(map_id(
            old.geometry.surface.0,base_s,target_s)?);
        next.geometry.outer_loop=LoopId(map_id(
            old.geometry.outer_loop.0,base_l as u32,target_l)?);
        model.faces.push(next);
        model.shells[0].faces.push(FaceId(map_id(
            old.id.0,base_f as u32,target_f)?));
    }
    let opening=LoopId(map_id(donor.faces[1].loops[1].0,
        base_l as u32,target_l)?);
    model.faces[1].loops.push(opening);
    let floor_local=donor.faces.last().ok_or_else(||invalid("missing donor floor"))?;
    if !matches!(donor.face_surface(floor_local.id),Some(Surface3::Plane(_))) ||
       floor_local.sense!=crate::FaceSense::Reversed {
        return Err(invalid("donor does not terminate in a real upward-facing floor"));
    }
    let floor=FaceId(map_id(floor_local.id.0,base_f as u32,target_f)?);
    let cylinders=donor.faces[base_f..].iter().filter(|f|
        matches!(donor.face_surface(f.id),Some(Surface3::Cylinder(_))))
        .map(|f|map_id(f.id.0,base_f as u32,target_f).map(FaceId))
        .collect::<Result<Vec<_>,_>>()?;
    let expected_cylinders=if kind==BlindRegionKind::SingleCircle{1}else{2};
    if cylinders.len()!=expected_cylinders {
        return Err(invalid("donor crossed region lost a cylindrical face"));
    }
    Ok(BlindRegionComponent{
        kind,source_indices:source_indices.to_vec(),
        opening_loop:opening,floor_face:floor,
        cylinder_faces:cylinders,first_edge:EdgeId(target_e),
        region_area:area,exposed_perimeter:perimeter,depth,
    })
}

#[derive(Debug,Clone)]
pub struct CompoundBlindBrep {
    source:Solid,
    holes:Vec<BlindHoleParameters>,
    components:Vec<BlindRegionComponent>,
    model:BrepModel,
    revision:u64,
    undo:Vec<Vec<BlindHoleParameters>>,
    redo:Vec<Vec<BlindHoleParameters>>,
}
impl CompoundBlindBrep {
    pub fn from_box(source:&Solid,holes:Vec<BlindHoleParameters>)
        ->Result<Self,GeometryError>
    {
        let (model,components)=BlindRegionComposer::build(source,&holes)?;
        Ok(Self{source:source.clone(),holes,components,model,revision:0,
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn holes(&self)->&[BlindHoleParameters]{&self.holes}
    pub fn components(&self)->&[BlindRegionComponent]{&self.components}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub(crate) fn source_box(&self)->&Solid{&self.source}
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (independent,components)=BlindRegionComposer::build(&self.source,&self.holes)?;
        if self.components!=components ||
            self.model.origin!=independent.origin ||
            self.model.summary()!=independent.summary() ||
            self.model.mass!=independent.mass ||
            self.model.bbox!=independent.bbox ||
            self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>()!=
                independent.vertices.iter().map(|v|v.position).collect::<Vec<_>>() ||
            self.model.geometry.curves3!=independent.geometry.curves3 ||
            self.model.geometry.curves2!=independent.geometry.curves2 ||
            self.model.geometry.surfaces!=independent.geometry.surfaces ||
            self.model.tolerance!=independent.tolerance ||
            self.model.units!=independent.units ||
            self.model.edges.iter().zip(&independent.edges).any(|(a,b)|
                a.id!=b.id || a.start!=b.start || a.end!=b.end ||
                a.coedges!=b.coedges) ||
            self.model.coedges.iter().zip(&independent.coedges).any(|(a,b)|
                a.id!=b.id || a.edge!=b.edge || a.face!=b.face ||
                a.loop_id!=b.loop_id || a.reversed!=b.reversed ||
                a.next!=b.next || a.prev!=b.prev || a.twin!=b.twin) ||
            self.model.loops.iter().zip(&independent.loops).any(|(a,b)|
                a.id!=b.id || a.face!=b.face || a.role!=b.role ||
                a.first_coedge!=b.first_coedge) ||
            self.model.faces.iter().zip(&independent.faces).any(|(a,b)|
                a.id!=b.id || a.loops!=b.loops ||
                a.sense!=b.sense || a.geometry!=b.geometry) ||
            self.model.shells.iter().zip(&independent.shells).any(|(a,b)|
                a.id!=b.id || a.faces!=b.faces || a.closed!=b.closed) ||
            self.model.geometry.edges!=independent.geometry.edges ||
            self.model.geometry.coedges!=independent.geometry.coedges {
            return Err(invalid("compound blind body diverges from source-rebuilt analytical regions"));
        }
        Ok(())
    }
    pub fn begin_edit(&mut self)->CompoundBlindTransaction<'_>{
        let next=self.holes.clone();
        CompoundBlindTransaction{owner:self,staged:next,poisoned:false,changes:0}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let next=if redo{&self.redo}else{&self.undo};
        let Some(params)=next.last() else{return Ok(false)};
        let (model,components)=BlindRegionComposer::build(&self.source,params)?;
        let nextrev=self.revision.checked_add(1)
            .ok_or_else(||invalid("compound analytical revision overflow"))?;
        let restored=params.clone();
        if redo{self.redo.pop();self.undo.push(self.holes.clone());}
        else{self.undo.pop();self.redo.push(self.holes.clone());}
        self.model=model;self.components=components;self.holes=restored;
        self.revision=nextrev;
        Ok(true)
    }
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::CompoundBlindDisplay,GeometryError>{
        crate::CompoundBlindMesher::tessellate(self,segments)
    }
}
impl crate::BrepSource for CompoundBlindBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<crate::Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct CompoundBlindTransaction<'a>{
    owner:&'a mut CompoundBlindBrep,
    staged:Vec<BlindHoleParameters>,
    poisoned:bool,
    changes:usize,
}
impl CompoundBlindTransaction<'_>{
    pub fn replace(&mut self,index:usize,value:BlindHoleParameters)
        ->Result<(),GeometryError>
    {
        if self.poisoned{return Err(invalid("compound blind transaction is poisoned"));}
        if index>=self.staged.len()||self.changes>=256 ||
            !value.center_x.is_finite()||!value.center_z.is_finite()||
            !value.radius.is_finite()||!value.depth.is_finite()||
            value.radius<=0.0||value.depth<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid compound source index, parameter or edit budget"));
        }
        self.staged[index]=value;self.changes+=1;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned compound blind transaction"));}
        if self.changes==0{return Err(invalid("empty compound blind transaction"));}
        let (candidate,components)=BlindRegionComposer::build(&self.owner.source,&self.staged)?;
        let next=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("compound blind revision overflow"))?;
        if self.owner.undo.len()>=128{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.holes.clone());
        self.owner.redo.clear();
        self.owner.model=candidate;self.owner.components=components;
        self.owner.holes=self.staged;self.owner.revision=next;
        Ok(next)
    }
}
