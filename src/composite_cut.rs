//! Canonical, strongly atomic composition of disjoint analytic cylindrical
//! subtractive features. A feature can be a through-hole or a top-entry
//! blind bore. Each cut is an independently certified existing RustSolid
//! producer; this module transplants its *real* topology/UV/3D carriers into
//! a single B-rep with shared original box faces, then recomputes exact mass.
//!
//! This is a bounded compositional executor, not arbitrary Boolean graph
//! arrangement or general TopologyDeltaPlan execution.
use crate::{
    BlindHoleBrep, BlindHoleParameters, BrepModel, BrepOrigin, BrepSource,
    Curve2Id, Curve3Id, EdgeId, FaceId, GeometryError, HoleParameters,
    HistoryDirection, LoopId, MassProperties, Mesh, ModelUnits, Point3, Solid,
    SurfaceId, ThroughHoleBrep, VertexId, CoedgeId,
};
use std::sync::atomic::{AtomicU64,Ordering};
use std::f64::consts::PI;

static NEXT_COMPOSITE_EDITOR:AtomicU64=AtomicU64::new(1);
static NEXT_COMPOSITE_GENERATION:AtomicU64=AtomicU64::new(1);
const CUT_LIMIT:usize=32;
const UNDO_LIMIT:usize=128;
const BATCH_LIMIT:usize=256;

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
fn fresh(counter:&AtomicU64,what:&'static str)->Result<u64,GeometryError>{
    let value=counter.fetch_add(1,Ordering::Relaxed);
    if value==0||value==u64::MAX{
        return Err(invalid(format!("{what} token overflow")));
    }
    Ok(value)
}
fn id(len:usize)->Result<u32,GeometryError>{
    u32::try_from(len).map_err(|_|invalid("composite B-rep dense identity overflow"))
}
fn mapped(original:u32,base:u32,start:u32)->u32{
    if original<base {original} else {start+(original-base)}
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub enum AnalyticCut {
    Through(HoleParameters),
    Blind(BlindHoleParameters),
}
impl AnalyticCut {
    pub fn radius(self)->f64{
        match self{Self::Through(h)=>h.radius,Self::Blind(h)=>h.radius}
    }
    pub fn center(self)->[f64;2]{
        match self{
            Self::Through(h)=>[h.center_x,h.center_z],
            Self::Blind(h)=>[h.center_x,h.center_z],
        }
    }
    pub fn depth(self,block:&Solid)->f64{
        match self {
            Self::Through(_)=>block.bbox.max.y-block.bbox.min.y,
            Self::Blind(h)=>h.depth,
        }
    }
    pub fn is_through(self)->bool{matches!(self,Self::Through(_))}
    pub fn kind(self)->&'static str{
        match self{Self::Through(_)=>"through",Self::Blind(_)=>"blind"}
    }
}
#[derive(Debug,Clone,Copy,PartialEq,Eq,Hash)]
pub struct CompositeCutHandle{
    pub editor_token:u64,
    pub slot:u32,
    pub generation:u64,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum CompositeHandleError{ForeignEditor,StaleFeature}
impl std::fmt::Display for CompositeHandleError{
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{
        match self {
            Self::ForeignEditor=>write!(f,"cut handle belongs to another composite editor"),
            Self::StaleFeature=>write!(f,"cut feature has been retired or reincarnated"),
        }
    }
}
impl std::error::Error for CompositeHandleError{}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct CompositeTopology{
    pub side_face:FaceId,
    pub floor_face:Option<FaceId>,
    pub lower_edge:EdgeId,
    pub upper_edge:EdgeId,
    pub seam_edge:EdgeId,
    pub lower_loop:LoopId,
    pub upper_loop:LoopId,
    pub side_loop:LoopId,
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CompositeCutBinding{
    pub handle:CompositeCutHandle,
    pub cut:AnalyticCut,
    pub topology:CompositeTopology,
}
#[derive(Debug,Clone)]
pub struct CompositeCommit{
    pub revision_before:u64,
    pub revision_after:u64,
    pub command_count:usize,
    pub created:Vec<CompositeCutBinding>,
    pub retired:Vec<CompositeCutBinding>,
    pub remapped:Vec<(CompositeCutHandle,CompositeTopology,CompositeTopology)>,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct CompositeTransition{
    pub direction:HistoryDirection,
    pub revision_before:u64,
    pub revision_after:u64,
}
#[derive(Debug,Clone)]
pub struct CompositeSnapshot{
    pub editor_token:u64,
    pub revision:u64,
    pub model:BrepModel,
    pub cuts:Vec<CompositeCutBinding>,
}
#[derive(Debug,Clone,Copy)]
struct Slot{
    generation:u64,
    current:Option<AnalyticCut>,
}
#[derive(Debug,Clone)]
struct State{
    slots:Vec<Slot>,
    model:BrepModel,
}
fn active(slots:&[Slot])->Vec<AnalyticCut>{
    slots.iter().filter_map(|s|s.current).collect()
}
fn check_set(source:&Solid,cuts:&[AnalyticCut])->Result<(),GeometryError>{
    crate::boolean::certify_box(source).map_err(|e|invalid(e.to_string()))?;
    if cuts.len()>CUT_LIMIT{
        return Err(invalid("more than 32 disjoint composite analytic cuts"));
    }
    let extent=source.bbox.max.sub(source.bbox.min);
    let epsilon=source.tolerance.length_at(extent.x.hypot(extent.y).hypot(extent.z));
    for (i,cut) in cuts.iter().copied().enumerate(){
        match cut {
            AnalyticCut::Through(p)=>{ThroughHoleBrep::from_box(source,p)?;},
            AnalyticCut::Blind(p)=>{BlindHoleBrep::from_box(source,p)?;},
        }
        for other in &cuts[..i] {
            let a=cut.center();let b=other.center();
            if (a[0]-b[0]).hypot(a[1]-b[1])
                <=cut.radius()+other.radius()+epsilon {
                return Err(invalid("composite cuts overlap or touch on their common opening face"));
            }
        }
    }
    Ok(())
}
/// Transplant the independently validated single-cut B-rep as a *genuine*
/// set of vertices, edges, coedges, loops, surfaces and directed UV trims.
/// The original 6 box faces/edges and shared carrier catalog prefix are
/// retained. Each fragment's cap inner loops are attached to the same
/// canonical box face; its extra side/floor faces get new dense IDs.
fn attach(model:&mut BrepModel,source:&BrepModel,fragment:&BrepModel)
    ->Result<CompositeTopology,GeometryError>
{
    fragment.validate()?;
    let source_v=id(source.vertices.len())?;
    let source_e=id(source.edges.len())?;
    let source_c=id(source.coedges.len())?;
    let source_l=id(source.loops.len())?;
    let source_f=id(source.faces.len())?;
    let start_v=id(model.vertices.len())?;
    let start_e=id(model.edges.len())?;
    let start_c=id(model.coedges.len())?;
    let start_l=id(model.loops.len())?;
    let start_f=id(model.faces.len())?;
    let source_g3=id(source.geometry.curves3.len())?;
    let source_g2=id(source.geometry.curves2.len())?;
    let source_gs=id(source.geometry.surfaces.len())?;
    let start_g3=id(model.geometry.curves3.len())?;
    let start_g2=id(model.geometry.curves2.len())?;
    let start_gs=id(model.geometry.surfaces.len())?;
    if source_v!=8 || source_e!=12 || source_c!=24 ||
       source_l!=6 || source_f!=6 ||
       fragment.vertices.len()!=source.vertices.len()+2 ||
       fragment.edges.len()!=source.edges.len()+3 ||
       fragment.coedges.len()!=source.coedges.len()+6 ||
       fragment.loops.len()!=source.loops.len()+3 ||
       !(fragment.faces.len()==7 || fragment.faces.len()==8){
        return Err(invalid("unsupported single-cut fragment topology for compositional assembly"));
    }
    let mv=|x:VertexId|VertexId(mapped(x.0,source_v,start_v));
    let me=|x:EdgeId|EdgeId(mapped(x.0,source_e,start_e));
    let mc=|x:CoedgeId|CoedgeId(mapped(x.0,source_c,start_c));
    let ml=|x:LoopId|LoopId(mapped(x.0,source_l,start_l));
    let mf=|x:FaceId|FaceId(mapped(x.0,source_f,start_f));
    let g3=|x:Curve3Id|Curve3Id(mapped(x.0,source_g3,start_g3));
    let g2=|x:Curve2Id|Curve2Id(mapped(x.0,source_g2,start_g2));
    let gs=|x:SurfaceId|SurfaceId(mapped(x.0,source_gs,start_gs));
    model.geometry.curves3.extend_from_slice(
        &fragment.geometry.curves3[source_g3 as usize..]);
    model.geometry.curves2.extend_from_slice(
        &fragment.geometry.curves2[source_g2 as usize..]);
    model.geometry.surfaces.extend_from_slice(
        &fragment.geometry.surfaces[source_gs as usize..]);
    for vertex in &fragment.vertices[source_v as usize..] {
        let mut clone=vertex.clone();
        clone.id=mv(vertex.id);
        model.vertices.push(clone);
    }
    for edge in &fragment.edges[source_e as usize..]{
        let mut clone=edge.clone();
        clone.id=me(edge.id);
        clone.start=mv(edge.start);clone.end=mv(edge.end);
        clone.coedges=edge.coedges.map(mc);
        model.edges.push(clone);
    }
    for coedge in &fragment.coedges[source_c as usize..] {
        let mut clone=coedge.clone();
        clone.id=mc(coedge.id);clone.edge=me(coedge.edge);
        clone.face=mf(coedge.face);clone.loop_id=ml(coedge.loop_id);
        clone.next=mc(coedge.next);clone.prev=mc(coedge.prev);
        clone.twin=mc(coedge.twin);
        model.coedges.push(clone);
    }
    for lp in &fragment.loops[source_l as usize..]{
        let mut clone=lp.clone();
        clone.id=ml(lp.id);clone.face=mf(lp.face);
        clone.first_coedge=mc(lp.first_coedge);
        model.loops.push(clone);
    }
    for cap in 0..2 {
        for inner in fragment.faces[cap].loops.iter().skip(1){
            model.faces[cap].loops.push(ml(*inner));
        }
    }
    for face in &fragment.faces[source_f as usize..] {
        let mut clone=face.clone();
        clone.id=mf(face.id);
        clone.loops=face.loops.iter().copied().map(ml).collect();
        clone.geometry.outer_loop=ml(face.geometry.outer_loop);
        clone.geometry.surface=gs(face.geometry.surface);
        model.shells[0].faces.push(clone.id);
        model.faces.push(clone);
    }
    for g in &fragment.geometry.edges[source_e as usize..]{
        let mut clone=*g;clone.curve=g3(g.curve);
        model.geometry.edges.push(clone);
    }
    for g in &fragment.geometry.coedges[source_c as usize..]{
        let mut clone=*g;
        clone.pcurve=g2(g.pcurve);
        clone.face=mf(g.face);clone.loop_id=ml(g.loop_id);
        model.geometry.coedges.push(clone);
    }
    let blind=fragment.faces.len()==8;
    Ok(CompositeTopology{
        side_face:mf(FaceId(6)),
        floor_face:blind.then(||mf(FaceId(7))),
        lower_edge:me(EdgeId(12)),upper_edge:me(EdgeId(13)),
        seam_edge:me(EdgeId(14)),
        lower_loop:ml(LoopId(6)),
        upper_loop:ml(LoopId(7)),
        side_loop:ml(LoopId(8)),
    })
}
fn construct(source:&Solid,cuts:&[AnalyticCut])
    ->Result<(BrepModel,Vec<CompositeTopology>),GeometryError>
{
    check_set(source,cuts)?;
    let baseline=BrepModel::from_polyhedral(source,ModelUnits::default())?;
    let mut result=baseline.clone();
    let mut topology=Vec::new();
    let original=source.mass.centroid;
    let mut removed=0.0;
    let mut correction=Point3{x:0.,y:0.,z:0.};
    let mut area=source.mass.surface_area;
    for &cut in cuts {
        let fragment=match cut {
            AnalyticCut::Through(p)=>ThroughHoleBrep::from_box(source,p)?.model().clone(),
            AnalyticCut::Blind(p)=>{
                let blind=BlindHoleBrep::from_box(source,p)?;
                // A direct parametric blind feature must pass the same
                // source-face closed IntersectionGraph certificate as a
                // cutter promoted through BooleanKernel::apply_sequence.
                blind.opening_trace()?;
                blind.model().clone()
            },
        };
        let t=attach(&mut result,&baseline,&fragment)?;
        topology.push(t);
        let depth=cut.depth(source);
        let r=cut.radius();
        let dv=PI*r*r*depth;
        if !dv.is_finite(){return Err(invalid("nonfinite cut volume"));}
        removed+=dv;
        let centroid=cut.center();
        let center_y=if cut.is_through(){original.y}
            else {source.bbox.max.y-depth/2.0};
        correction=correction.add(Point3{
            x:(original.x-centroid[0])*dv,
            y:(original.y-center_y)*dv,
            z:(original.z-centroid[1])*dv,
        });
        area+=if cut.is_through(){
            2.0*PI*r*depth-2.0*PI*r*r
        }else{2.0*PI*r*depth};
    }
    let volume=source.mass.volume-removed;
    if !volume.is_finite() || volume<=0.0 ||
        !area.is_finite() || area<=0.0 {
        return Err(invalid("composite cut material mass or surface area invalid"));
    }
    result.mass=MassProperties{
        volume,surface_area:area,
        centroid:original.add(correction.scale(1.0/volume)),
    };
    result.origin=BrepOrigin::AnalyticCompositeCuts;
    result.validate()?;
    let through_count=cuts.iter().filter(|c|c.is_through()).count();
    let expected_euler=2-2*through_count as isize;
    if result.euler_characteristic()!=expected_euler ||
        result.summary().genus!=through_count {
        return Err(invalid("compositional cuts violate Euler/genus contract"));
    }
    Ok((result,topology))
}

fn bindings(token:u64,slots:&[Slot])->Vec<CompositeCutBinding>{
    // Each cut adds exactly 2 vertices, 3 edges and 6 coedges, but the
    // number of new faces depends on the cut kind. Recompute from kinds.
    let mut edge=12u32;let mut face=6u32;let mut lp=6u32;
    let mut items=Vec::new();
    for (i,s) in slots.iter().enumerate(){
        if let Some(cut)=s.current {
            items.push(CompositeCutBinding{
                handle:CompositeCutHandle{editor_token:token,
                    slot:i as u32,generation:s.generation},
                cut,
                topology:CompositeTopology{
                    side_face:FaceId(face),
                    floor_face:matches!(cut,AnalyticCut::Blind(_))
                        .then(||FaceId(face+1)),
                    lower_edge:EdgeId(edge),upper_edge:EdgeId(edge+1),
                    seam_edge:EdgeId(edge+2),
                    lower_loop:LoopId(lp),upper_loop:LoopId(lp+1),
                    side_loop:LoopId(lp+2),
                },
            });
            face+=if cut.is_through(){1}else{2};
            edge+=3;lp+=3;
        }
    }
    items
}
/// Independent source-owned editor. Real surfaces and loops are rebuilt
/// off to the side; neither stale handles nor partial candidates are exposed.
#[derive(Debug)]
pub struct BrepCompositeEditor{
    token:u64,
    revision:u64,
    source:Solid,
    state:State,
    undo:Vec<State>,
    redo:Vec<State>,
}
/// Cloning an editor creates a **new** feature-handle namespace. Two
/// independently editable branches must never accept each other's
/// process-local handles, even if their dense topology is initially equal.
impl Clone for BrepCompositeEditor {
    fn clone(&self)->Self{
        let token=fresh(&NEXT_COMPOSITE_EDITOR,"cloned composite editor")
            .expect("composite editor token capacity exhausted");
        let mut state=self.state.clone();
        for slot in &mut state.slots {
            slot.generation=fresh(&NEXT_COMPOSITE_GENERATION,
                "cloned composite feature generation")
                .expect("composite handle generation exhausted");
        }
        Self{token,revision:self.revision,source:self.source.clone(),
            state,undo:self.undo.clone(),redo:self.redo.clone()}
    }
}
impl BrepCompositeEditor{
    pub fn new(source:Solid)->Result<Self,GeometryError>{
        crate::boolean::certify_box(&source).map_err(|e|invalid(e.to_string()))?;
        let token=fresh(&NEXT_COMPOSITE_EDITOR,"composite editor")?;
        let model=construct(&source,&[])?.0;
        Ok(Self{token,revision:0,source,
            state:State{slots:Vec::new(),model},undo:Vec::new(),redo:Vec::new()})
    }
    pub fn from_cuts(source:Solid,cuts:&[AnalyticCut])->Result<Self,GeometryError>{
        let mut out=Self::new(source)?;
        let mut slots=Vec::with_capacity(cuts.len());
        for &cut in cuts{
            slots.push(Slot{generation:fresh(&NEXT_COMPOSITE_GENERATION,
                "composite cut generation")?,current:Some(cut)});
        }
        out.state=State{slots,model:construct(&out.source,cuts)?.0};
        Ok(out)
    }
    pub fn source_box(&self)->&Solid{&self.source}
    pub fn model(&self)->&BrepModel{&self.state.model}
    pub fn revision(&self)->u64{self.revision}
    pub fn editor_token(&self)->u64{self.token}
    pub fn cut_count(&self)->usize{active(&self.state.slots).len()}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn bindings(&self)->Vec<CompositeCutBinding>{
        bindings(self.token,&self.state.slots)
    }
    pub fn handle(&self,slot:u32)->Option<CompositeCutHandle>{
        self.state.slots.get(slot as usize).filter(|s|s.current.is_some())
            .map(|s|CompositeCutHandle{editor_token:self.token,
                slot,generation:s.generation})
    }
    pub fn resolve(&self,h:CompositeCutHandle)
        ->Result<CompositeCutBinding,CompositeHandleError>
    {
        if h.editor_token!=self.token{
            return Err(CompositeHandleError::ForeignEditor);
        }
        self.bindings().into_iter().find(|b|b.handle==h)
            .ok_or(CompositeHandleError::StaleFeature)
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.state.model.validate()?;
        let (expected,_)=construct(&self.source,&active(&self.state.slots))?;
        if self.state.model.summary()!=expected.summary() ||
           self.state.model.mass!=expected.mass ||
           self.state.model.bbox!=expected.bbox ||
           self.state.model.geometry.curves3!=expected.geometry.curves3 ||
           self.state.model.geometry.curves2!=expected.geometry.curves2 ||
           self.state.model.geometry.surfaces!=expected.geometry.surfaces ||
           self.state.model.faces.iter().zip(&expected.faces)
               .any(|(a,b)|a.id!=b.id || a.sense!=b.sense || a.loops!=b.loops) {
            return Err(invalid("composite edit diverged from independent canonical analytic rebuild"));
        }
        Ok(())
    }
    pub fn snapshot(&self)->CompositeSnapshot{
        CompositeSnapshot{editor_token:self.token,revision:self.revision,
            model:self.state.model.clone(),cuts:self.bindings()}
    }
    pub fn check_snapshot(&self,s:&CompositeSnapshot)->Result<(),GeometryError>{
        if s.editor_token!=self.token {
            return Err(invalid("snapshot belongs to another composite editor"));
        }
        if s.revision!=self.revision{
            return Err(invalid("stale composite revision"));
        }
        Ok(())
    }
    pub fn begin_transaction(&mut self)->CompositeCutTransaction<'_>{
        let slots=self.state.slots.clone();
        CompositeCutTransaction{owner:self,slots,failed:false,commands:0}
    }
    pub fn undo(&mut self)->Result<Option<CompositeTransition>,GeometryError>{
        self.transition(false)
    }
    pub fn redo(&mut self)->Result<Option<CompositeTransition>,GeometryError>{
        self.transition(true)
    }
    fn transition(&mut self,redo:bool)->Result<Option<CompositeTransition>,GeometryError>{
        let Some(target)=(if redo {self.redo.last()}else{self.undo.last()}) else {
            return Ok(None);
        };
        let mut candidate=target.clone();
        // A resurrected incarnation always obtains a fresh generation;
        // equal sparse slots alone never confer identity on replacement.
        for (i,slot) in candidate.slots.iter_mut().enumerate(){
            let current=self.state.slots.get(i);
            if slot.current.is_some() {
                if !current.is_some_and(|old|
                    old.current.is_some() && old.generation==slot.generation) {
                    slot.generation=fresh(&NEXT_COMPOSITE_GENERATION,
                        "restored composite generation")?;
                }
            }
        }
        let expected=construct(&self.source,&active(&candidate.slots))?.0;
        candidate.model=expected;
        let before=self.revision;
        let after=before.checked_add(1)
            .ok_or_else(||invalid("composite revision overflow"))?;
        let old=std::mem::replace(&mut self.state,candidate);
        if redo {self.redo.pop();self.undo.push(old);}
        else {self.undo.pop();self.redo.push(old);}
        self.revision=after;
        Ok(Some(CompositeTransition{direction:if redo{
            HistoryDirection::Redo}else{HistoryDirection::Undo},
            revision_before:before,revision_after:after}))
    }
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::CompositeDisplay,GeometryError>{
        crate::CompositeRegionMesher::tessellate(self,segments)
    }
}
impl BrepSource for BrepCompositeEditor{
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model().clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
/// Compound feature and topology mutation, with one revision, one
/// Undo step, and no observable partial result.
pub struct CompositeCutTransaction<'a>{
    owner:&'a mut BrepCompositeEditor,
    slots:Vec<Slot>,
    failed:bool,
    commands:usize,
}
impl CompositeCutTransaction<'_>{
    fn guard(&mut self)->Result<(),GeometryError>{
        if self.failed {return Err(invalid("composite transaction poisoned"));}
        if self.commands>=BATCH_LIMIT{
            self.failed=true;
            return Err(invalid("composite batch exceeds 256 commands"));
        }
        Ok(())
    }
    fn checked(&mut self,h:CompositeCutHandle)->Result<usize,GeometryError>{
        if h.editor_token!=self.owner.token {
            self.failed=true;return Err(invalid("foreign composite cut handle"));
        }
        let index=h.slot as usize;
        if !self.slots.get(index).is_some_and(|slot|
            slot.current.is_some()&&slot.generation==h.generation){
            self.failed=true;return Err(invalid("stale or retired composite cut handle"));
        }
        Ok(index)
    }
    pub fn handle(&self,slot:u32)->Option<CompositeCutHandle>{
        if self.failed{return None;}
        self.slots.get(slot as usize).filter(|s|s.current.is_some())
            .map(|s|CompositeCutHandle{editor_token:self.owner.token,slot,
                generation:s.generation})
    }
    pub fn add(&mut self,cut:AnalyticCut)->Result<CompositeCutHandle,GeometryError>{
        self.guard()?;
        let index=self.slots.iter().position(|s|s.current.is_none())
            .unwrap_or(self.slots.len());
        let generation=fresh(&NEXT_COMPOSITE_GENERATION,"composite generation")?;
        if index==self.slots.len(){
            self.slots.push(Slot{current:Some(cut),generation});
        }else{
            self.slots[index]=Slot{current:Some(cut),generation};
        }
        if let Err(e)=check_set(&self.owner.source,&active(&self.slots)){
            self.failed=true;return Err(e);
        }
        self.commands+=1;
        Ok(CompositeCutHandle{editor_token:self.owner.token,slot:index as u32,generation})
    }
    pub fn remove(&mut self,h:CompositeCutHandle)->Result<(),GeometryError>{
        self.guard()?;
        let i=self.checked(h)?;
        self.slots[i].current=None;
        self.commands+=1;Ok(())
    }
    pub fn replace(&mut self,h:CompositeCutHandle,new_cut:AnalyticCut)
        ->Result<(),GeometryError>
    {
        self.guard()?;
        let i=self.checked(h)?;
        self.slots[i].current=Some(new_cut);
        if let Err(e)=check_set(&self.owner.source,&active(&self.slots)){
            self.failed=true;return Err(e);
        }
        self.commands+=1;Ok(())
    }
    pub fn staged_count(&self)->usize{self.commands}
    pub fn commit(self)->Result<CompositeCommit,GeometryError>{
        if self.failed{return Err(invalid("cannot commit poisoned composite transaction"));}
        if self.commands==0{return Err(invalid("empty composite edit batch"));}
        let cuts=active(&self.slots);
        let (model,_)=construct(&self.owner.source,&cuts)?;
        let before=self.owner.revision;
        let after=before.checked_add(1)
            .ok_or_else(||invalid("composite revision overflow"))?;
        let old=bindings(self.owner.token,&self.owner.state.slots);
        let new=bindings(self.owner.token,&self.slots);
        let created=new.iter().filter(|b|!old.iter().any(|o|o.handle==b.handle))
            .copied().collect();
        let retired=old.iter().filter(|b|!new.iter().any(|n|n.handle==b.handle))
            .copied().collect();
        let mut remapped=Vec::new();
        for n in &new {
            if let Some(o)=old.iter().find(|old|old.handle==n.handle){
                if o.topology!=n.topology {
                    remapped.push((n.handle,o.topology,n.topology));
                }
            }
        }
        let state=State{slots:self.slots,model};
        let previous=std::mem::replace(&mut self.owner.state,state);
        if self.owner.undo.len()>=UNDO_LIMIT {self.owner.undo.remove(0);}
        self.owner.undo.push(previous);
        self.owner.redo.clear();
        self.owner.revision=after;
        Ok(CompositeCommit{revision_before:before,revision_after:after,
            command_count:self.commands,created,retired,remapped})
    }
}
