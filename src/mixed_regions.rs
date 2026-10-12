//! Heterogeneous analytic aperture composition: one closed connected
//! axis-aligned box shell with both genuine through-cylinder channels
//! and independent blind-cut regions, each independently source certified.
//!
//! A blind region is either one complete circular bore or two partially
//! overlapping equal-depth bores. Through footprints must be strictly
//! separated from every blind footprint. This first slice does not claim
//! crossings between through/finite-depth cuts, arbitrary depth steps,
//! universal polygonal CSG or a general TopologyDeltaPlan executor.
//! ALL region boundaries are actual 3D curves, 2D coedge trims and surfaces.
use crate::{
    BlindHoleBrep,BlindHoleParameters,BlindRegionComponent,BlindRegionKind,
    BrepModel,BrepOrigin,BrepSource,CrossedDiskRegion,FaceId,FaceSense,
    GeometryError,HoleParameters,LoopRole,MassProperties,Mesh,
    Point3,Solid,Surface3,ThroughHoleBrep,
};
use crate::compound_blind::{graft_region,plan_groups};
use crate::multi_through_hole::build_multiple;
use std::f64::consts::{PI,TAU};

const MAX_TOTAL_CUTS:usize=16;
const MAX_HISTORY:usize=128;
fn invalid(msg:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(msg.into())
}
#[derive(Debug,Clone,PartialEq)]
pub struct MixedThroughBlindParameters {
    pub through:Vec<HoleParameters>,
    pub blind:Vec<BlindHoleParameters>,
}
/// A source-proven, mathematical (not faceted) mixed-region B-rep compositor.
#[derive(Debug,Clone,Copy,Default)]
pub struct MixedRegionComposer;
impl MixedRegionComposer {
    pub fn build(source:&Solid,parameters:&MixedThroughBlindParameters)
        ->Result<(BrepModel,Vec<BlindRegionComponent>),GeometryError>
    {
        crate::boolean::certify_box(source).map_err(|e|invalid(e.to_string()))?;
        let through=&parameters.through;
        let blind=&parameters.blind;
        if through.is_empty() || blind.is_empty() ||
            through.len()+blind.len()>MAX_TOTAL_CUTS {
            return Err(invalid("heterogeneous cut requires at least one through and one blind bore, with at most 16 total cylinders"));
        }
        // This independent builder is also responsible for geometry and
        // exact mass checks of the *entire* through-cut family.
        let mut model=build_multiple(source,through)?;
        let b=source.bbox;
        let dx=b.max.x-b.min.x;
        let dy=b.max.y-b.min.y;
        let dz=b.max.z-b.min.z;
        let eps=source.tolerance.length_at(dx.max(dy).max(dz));
        for t in through {
            for h in blind {
                let gap=(t.center_x-h.center_x)
                    .hypot(t.center_z-h.center_z)-(t.radius+h.radius);
                if !gap.is_finite()||gap<=eps {
                    return Err(invalid("through and blind circle footprints must be strictly disjoint with proven clearance"));
                }
            }
        }
        // Shared source-plane proof for every blind circle and every
        // possible two-circle crossing group, including strict rejection
        // of third-way crossings, containment, tangency and depth steps.
        let groups=plan_groups(source,blind)?;
        let old=source.mass.centroid;
        let mut removed=0.0;
        let mut relative=Point3{x:0.0,y:0.0,z:0.0};
        let mut surface_area=source.mass.surface_area;
        for t in through {
            let area=PI*t.radius*t.radius;
            let volume=area*dy;
            removed+=volume;
            relative=relative.add(Point3{
                x:(t.center_x-old.x)*volume,
                y:0.0,
                z:(t.center_z-old.z)*volume,
            });
            surface_area+=TAU*t.radius*dy-2.0*area;
        }
        let mut components=Vec::with_capacity(groups.len());
        for group in groups {
            let (donor,area,perimeter,depth,moment,kind)=if group.len()==1{
                let h=blind[group[0]];
                let body=BlindHoleBrep::from_box(source,h)?;
                body.opening_trace()?;
                body.validate()?;
                let area=PI*h.radius*h.radius;
                (body.model().clone(),area,TAU*h.radius,h.depth,
                    [(h.center_x-old.x)*area,(h.center_z-old.z)*area],
                    BlindRegionKind::SingleCircle)
            }else{
                let h=[blind[group[0]],blind[group[1]]];
                let body=crate::IntersectingBlindBrep::from_box(source,h)?;
                body.validate()?;
                let r=body.region();
                // Recheck the mathematical overlap independently of the
                // mixed graft and of facet triangulation.
                let independently=CrossedDiskRegion::from_box(source,h)?;
                if (r.disk_union_area-independently.disk_union_area).abs()>eps*eps {
                    return Err(invalid("overlapping blind region area diverged from source proof"));
                }
                (body.model().clone(),r.disk_union_area,
                    r.exposed_perimeter,r.cut_depth,r.relative_area_moment,
                    BlindRegionKind::TwoCrossingCircles)
            };
            let descriptor=graft_region(&mut model,&donor,
                &group,kind,area,perimeter,depth)?;
            let delta_volume=area*depth;
            removed+=delta_volume;
            relative=relative.add(Point3{
                x:moment[0]*depth,
                y:delta_volume*(b.max.y-depth/2.0-old.y),
                z:moment[1]*depth,
            });
            surface_area+=perimeter*depth;
            components.push(descriptor);
        }
        let volume=source.mass.volume-removed;
        if !volume.is_finite()||volume<=0.0||
           !surface_area.is_finite()||surface_area<=0.0||
           !relative.is_finite() {
            return Err(invalid("mixed region analytical mass is invalid"));
        }
        model.mass=MassProperties{
            volume,surface_area,
            centroid:old.sub(relative.scale(1.0/volume)),
        };
        model.origin=BrepOrigin::AnalyticMixedThroughBlindRegions;
        model.validate()?;
        // The **genus is the number of through bores**. Adding a blind
        // floor does not create another tunnel through the body.
        let t=through.len();
        let singles=components.iter().filter(|c|
            c.kind==BlindRegionKind::SingleCircle).count();
        let crossed=components.len()-singles;
        let stats=model.summary();
        if stats.vertices!=8+2*t+2*singles+4*crossed ||
           stats.edges!=12+3*t+3*singles+6*crossed ||
           stats.coedges!=2*stats.edges ||
           stats.faces!=6+t+2*singles+3*crossed ||
           stats.loops!=6+3*t+3*singles+4*crossed ||
           stats.genus!=t ||
           stats.euler_characteristic!=2-2*t as isize ||
           stats.shells!=1 ||
           stats.circular_edges!=2*t+2*singles+4*crossed ||
           stats.cylindrical_faces!=t+singles+2*crossed ||
           model.faces[0].loops.len()!=1+t ||
           model.faces[1].loops.len()!=1+t+components.len() ||
           model.loops.iter().filter(|l|l.role==LoopRole::Inner).count()!=
                2*t+components.len() {
            return Err(invalid("heterogeneous analytical trim ownership, Euler characteristic or genus mismatch"));
        }
        for f in &model.faces {
            let expected=if f.id.0>=6 {
                FaceSense::Reversed
            }else{FaceSense::Forward};
            if f.sense!=expected {
                return Err(invalid("heterogeneous material shell face orientation disagrees with source"));
            }
        }
        for c in &components{
            if model.face_sense(c.floor_face)!=Some(FaceSense::Reversed) ||
               !matches!(model.face_surface(c.floor_face),Some(Surface3::Plane(_))) ||
               c.cylinder_faces.iter().any(|&id|
                   model.face_sense(id)!=Some(FaceSense::Reversed) ||
                   !matches!(model.face_surface(id),Some(Surface3::Cylinder(_)))) {
                return Err(invalid("heterogeneous blind floor/wall analytical carrier was not preserved"));
            }
        }
        for id in 6..6+t{
            let face=FaceId(id as u32);
            if model.face_sense(face)!=Some(FaceSense::Reversed) ||
               !matches!(model.face_surface(face),Some(Surface3::Cylinder(_))) {
                return Err(invalid("through-bore analytical cylindrical carrier lost source identity"));
            }
        }
        Ok((model,components))
    }
}
/// Strong-atomic parametric owner for a single connected genus-N mixed B-rep.
/// The model is always rebuilt from mathematical source features, not from
/// a faceted Boolean result or by mutating a shared donor in-place.
#[derive(Debug,Clone)]
pub struct MixedThroughBlindBrep {
    source:Solid,
    parameters:MixedThroughBlindParameters,
    components:Vec<BlindRegionComponent>,
    model:BrepModel,
    revision:u64,
    undo:Vec<MixedThroughBlindParameters>,
    redo:Vec<MixedThroughBlindParameters>,
}
impl MixedThroughBlindBrep {
    pub fn from_box(source:&Solid,params:MixedThroughBlindParameters)
        ->Result<Self,GeometryError>
    {
        let (model,components)=MixedRegionComposer::build(source,&params)?;
        Ok(Self{source:source.clone(),parameters:params,components,model,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn source_box(&self)->&Solid{&self.source}
    pub fn parameters(&self)->&MixedThroughBlindParameters{&self.parameters}
    pub fn through_holes(&self)->&[HoleParameters]{&self.parameters.through}
    pub fn blind_holes(&self)->&[BlindHoleParameters]{&self.parameters.blind}
    pub fn components(&self)->&[BlindRegionComponent]{&self.components}
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (next,components)=MixedRegionComposer::build(
            &self.source,&self.parameters)?;
        if next.summary()!=self.model.summary() ||
           next.origin!=self.model.origin ||
           next.mass!=self.model.mass ||
           next.bbox!=self.model.bbox ||
           components!=self.components ||
           next.geometry.curves3!=self.model.geometry.curves3 ||
           next.geometry.curves2!=self.model.geometry.curves2 ||
           next.geometry.surfaces!=self.model.geometry.surfaces ||
           next.vertices.iter().map(|v|v.position).collect::<Vec<_>>()!=
                self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>() {
            return Err(invalid("mixed B-rep differs from source-proven independent canonical rebuild"));
        }
        Ok(())
    }
    pub fn begin_edit(&mut self)->MixedThroughBlindTransaction<'_>{
        let staged=self.parameters.clone();
        MixedThroughBlindTransaction{owner:self,staged,poisoned:false}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let desired=if redo{&self.redo}else{&self.undo};
        let Some(target)=desired.last()else{return Ok(false)};
        let next=target.clone();
        let (model,components)=MixedRegionComposer::build(&self.source,&next)?;
        let revision=self.revision.checked_add(1)
            .ok_or_else(||invalid("mixed-region revision overflow"))?;
        let old=self.parameters.clone();
        if redo{self.redo.pop();self.undo.push(old);}
        else{self.undo.pop();self.redo.push(old);}
        self.parameters=next;self.model=model;self.components=components;
        self.revision=revision;
        Ok(true)
    }
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::MixedHoleDisplay,GeometryError>
    {
        crate::MixedHoleMesher::tessellate(self,segments)
    }
}
pub struct MixedThroughBlindTransaction<'a>{
    owner:&'a mut MixedThroughBlindBrep,
    staged:MixedThroughBlindParameters,
    poisoned:bool,
}
impl MixedThroughBlindTransaction<'_>{
    pub fn replace_through(&mut self,index:usize,next:HoleParameters)
        ->Result<(),GeometryError>
    {
        if self.poisoned{return Err(invalid("mixed-region transaction is poisoned"));}
        if index>=self.staged.through.len() ||
           !next.center_x.is_finite()||!next.center_z.is_finite()||
           !next.radius.is_finite()||next.radius<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid through-bore replacement selector/parameters"));
        }
        self.staged.through[index]=next;Ok(())
    }
    pub fn replace_blind(&mut self,index:usize,next:BlindHoleParameters)
        ->Result<(),GeometryError>
    {
        if self.poisoned{return Err(invalid("mixed-region transaction is poisoned"));}
        if index>=self.staged.blind.len() ||
           !next.center_x.is_finite()||!next.center_z.is_finite()||
           !next.radius.is_finite()||next.radius<=0.0 ||
           !next.depth.is_finite()||next.depth<=0.0 {
            self.poisoned=true;
            return Err(invalid("invalid blind-bore replacement selector/parameters"));
        }
        self.staged.blind[index]=next;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned heterogeneous B-rep transaction"));}
        let (model,components)=MixedRegionComposer::build(
            &self.owner.source,&self.staged)?;
        let revision=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("heterogeneous edit revision overflow"))?;
        if self.owner.undo.len()>=MAX_HISTORY{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters.clone());
        self.owner.redo.clear();
        self.owner.parameters=self.staged;
        self.owner.model=model;self.owner.components=components;
        self.owner.revision=revision;
        Ok(revision)
    }
}
impl BrepSource for MixedThroughBlindBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
