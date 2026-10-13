//! Narrow executable composition of separated through/blind crossing clusters.
//! Each cluster has one certified v0.31 transverse crossing. Isolated blind
//! disks and other crossing clusters must have provable positive clearance.
//! Shared top/bottom caps, coedges, and carrier references form ONE B-rep shell.
use crate::{
    BlindHoleBrep,BlindHoleParameters,BrepModel,BrepOrigin,BrepSource,
    CoedgeId,Curve2Id,Curve3Id,EdgeId,FaceId,GeometryError,HoleParameters,
    LoopId,LoopRole,MassProperties,Mesh,MixedCrossingBrep,
    MultiCutPlanState,MultiCutRegionPlan,MultiCutRegionPlanner,Point3,
    Solid,SurfaceId,VertexId,
};
const MAX_HISTORY:usize=128;
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn id(n:usize)->Result<u32,GeometryError>{
    u32::try_from(n).map_err(|_|invalid("composed shell entity ID overflow"))
}
fn mv(old:u32,base:u32,target:u32)->Result<u32,GeometryError>{
    if old<base {Ok(old)}else{target.checked_add(old-base)
        .ok_or_else(||invalid("grafted entity index overflow"))}
}
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CertifiedCrossingPair {
    pub through:HoleParameters,
    pub blind:BlindHoleParameters,
}
#[derive(Debug,Clone,PartialEq)]
pub struct ComposedCrossingParameters {
    pub pairs:Vec<CertifiedCrossingPair>,
    pub isolated_blind:Vec<BlindHoleParameters>,
}
fn certify_group_clearance(p:&ComposedCrossingParameters,eps:f64)
    ->Result<(),GeometryError>{
    let mut groups=p.pairs.iter().map(|x|vec![
        (x.through.center_x,x.through.center_z,x.through.radius),
        (x.blind.center_x,x.blind.center_z,x.blind.radius),
    ]).collect::<Vec<_>>();
    groups.extend(p.isolated_blind.iter().map(|b|
        vec![(b.center_x,b.center_z,b.radius)]));
    for i in 0..groups.len(){
        for j in i+1..groups.len(){
            for &(ax,az,ar) in &groups[i]{
                for &(bx,bz,br) in &groups[j]{
                    let gap=(ax-bx).hypot(az-bz)-ar-br;
                    if !gap.is_finite()||gap<=eps*32.0{
                        return Err(invalid("inter-cluster source circles are not strictly separated"));
                    }
                }
            }
        }
    }
    Ok(())
}
/// The first 8 vertices, 12 edges, 24 coedges, six faces/loops and their
/// analytical carrier catalogs belong to a single shared original box.
/// Copy only certified donor additions and remap all references atomically.
fn graft(m:&mut BrepModel,d:&BrepModel,through:bool)->Result<(),GeometryError>{
    d.validate()?;
    const V:u32=8;const E:u32=12;const C:u32=24;
    const L:u32=6;const F:u32=6;
    if d.vertices.len()<=V as usize||d.edges.len()<=E as usize||
        d.coedges.len()<=C as usize||d.faces.len()<=F as usize||
        d.loops.len()<=L as usize||d.shells.len()!=1||
        d.faces[0].loops.len()!=if through{2}else{1}||
        d.faces[1].loops.len()!=2 {
        return Err(invalid("donor not a certified box-based aperture"));
    }
    if d.vertices[..V as usize].iter().map(|v|v.position).collect::<Vec<_>>()!=
            m.vertices[..V as usize].iter().map(|v|v.position).collect::<Vec<_>>()||
        d.geometry.curves3[..E as usize]!=m.geometry.curves3[..E as usize]||
        d.geometry.curves2[..C as usize]!=m.geometry.curves2[..C as usize]||
        d.geometry.surfaces[..F as usize]!=m.geometry.surfaces[..F as usize]||
        d.tolerance!=m.tolerance || d.units!=m.units{
        return Err(invalid("donor source geometry does not match box prefix"));
    }
    let (v,e,c,l,f)=(id(m.vertices.len())?,id(m.edges.len())?,
        id(m.coedges.len())?,id(m.loops.len())?,id(m.faces.len())?);
    let (c3,c2,s)=(id(m.geometry.curves3.len())?,
        id(m.geometry.curves2.len())?,id(m.geometry.surfaces.len())?);
    m.geometry.curves3.extend_from_slice(&d.geometry.curves3[E as usize..]);
    m.geometry.curves2.extend_from_slice(&d.geometry.curves2[C as usize..]);
    m.geometry.surfaces.extend_from_slice(&d.geometry.surfaces[F as usize..]);
    for a in &d.vertices[V as usize..]{
        let mut n=a.clone();n.id=VertexId(mv(a.id.0,V,v)?);m.vertices.push(n);
    }
    for (a,g) in d.edges[E as usize..].iter()
        .zip(&d.geometry.edges[E as usize..]){
        let mut n=a.clone();
        n.id=EdgeId(mv(a.id.0,E,e)?);
        n.start=VertexId(mv(a.start.0,V,v)?);
        n.end=VertexId(mv(a.end.0,V,v)?);
        n.coedges=[CoedgeId(mv(a.coedges[0].0,C,c)?),
            CoedgeId(mv(a.coedges[1].0,C,c)?)];
        m.edges.push(n);
        let mut x=*g;x.curve=Curve3Id(mv(g.curve.0,E,c3)?);
        m.geometry.edges.push(x);
    }
    for (a,g) in d.coedges[C as usize..].iter()
        .zip(&d.geometry.coedges[C as usize..]){
        let mut n=a.clone();
        n.id=CoedgeId(mv(a.id.0,C,c)?);
        n.edge=EdgeId(mv(a.edge.0,E,e)?);
        n.face=FaceId(mv(a.face.0,F,f)?);
        n.loop_id=LoopId(mv(a.loop_id.0,L,l)?);
        n.next=CoedgeId(mv(a.next.0,C,c)?);
        n.prev=CoedgeId(mv(a.prev.0,C,c)?);
        n.twin=CoedgeId(mv(a.twin.0,C,c)?);
        m.coedges.push(n);
        let mut x=*g;
        x.pcurve=Curve2Id(mv(g.pcurve.0,C,c2)?);
        x.face=FaceId(mv(g.face.0,F,f)?);
        x.loop_id=LoopId(mv(g.loop_id.0,L,l)?);
        m.geometry.coedges.push(x);
    }
    for a in &d.loops[L as usize..]{
        let mut n=a.clone();
        n.id=LoopId(mv(a.id.0,L,l)?);
        n.face=FaceId(mv(a.face.0,F,f)?);
        n.first_coedge=CoedgeId(mv(a.first_coedge.0,C,c)?);
        m.loops.push(n);
    }
    for a in &d.faces[F as usize..]{
        let mut n=a.clone();
        n.id=FaceId(mv(a.id.0,F,f)?);
        n.loops=a.loops.iter().map(|x|mv(x.0,L,l).map(LoopId))
            .collect::<Result<Vec<_>,_>>()?;
        n.geometry.surface=SurfaceId(mv(a.geometry.surface.0,F,s)?);
        n.geometry.outer_loop=LoopId(mv(a.geometry.outer_loop.0,L,l)?);
        m.shells[0].faces.push(n.id);
        m.faces.push(n);
    }
    if through{
        let lp=LoopId(mv(d.faces[0].loops[1].0,L,l)?);
        if m.loops[lp.0 as usize].face!=FaceId(0)||
            m.loops[lp.0 as usize].role!=LoopRole::Inner{
            return Err(invalid("crossing bottom cap donor lost real inner loop"));
        }
        m.faces[0].loops.push(lp);
    }
    let lp=LoopId(mv(d.faces[1].loops[1].0,L,l)?);
    if m.loops[lp.0 as usize].face!=FaceId(1)||
        m.loops[lp.0 as usize].role!=LoopRole::Inner{
        return Err(invalid("top cap donor has invalid inner loop"));
    }
    m.faces[1].loops.push(lp);
    Ok(())
}
#[derive(Debug,Clone)]
pub struct ComposedCrossingBrep {
    source:Solid,
    parameters:ComposedCrossingParameters,
    model:BrepModel,
    slab_plan:MultiCutRegionPlan,
    revision:u64,
    undo:Vec<ComposedCrossingParameters>,
    redo:Vec<ComposedCrossingParameters>,
}
impl ComposedCrossingBrep{
    pub fn from_box(source:&Solid,p:ComposedCrossingParameters)
        ->Result<Self,GeometryError>{
        let (model,plan)=Self::build(source,&p)?;
        Ok(Self{source:source.clone(),parameters:p,model,slab_plan:plan,
            revision:0,undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn parameters(&self)->&ComposedCrossingParameters{&self.parameters}
    pub fn slab_plan(&self)->&MultiCutRegionPlan{&self.slab_plan}
    pub(crate) fn source_box(&self)->&Solid{&self.source}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn tessellate(&self,n:usize)->Result<crate::ComposedCrossingDisplay,GeometryError>{
        crate::ComposedCrossingMesher::tessellate(self,n)
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (m,p)=Self::build(&self.source,&self.parameters)?;
        let matching=self.model.summary()==m.summary()&&
            self.model.origin==m.origin&&self.model.mass==m.mass&&
            self.model.bbox==m.bbox&&self.model.tolerance==m.tolerance&&
            self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>()==
                m.vertices.iter().map(|v|v.position).collect::<Vec<_>>()&&
            self.model.edges.iter().zip(&m.edges).all(|(a,b)|
                a.id==b.id&&a.start==b.start&&a.end==b.end&&a.coedges==b.coedges)&&
            self.model.coedges.iter().zip(&m.coedges).all(|(a,b)|
                a.id==b.id&&a.edge==b.edge&&a.face==b.face&&
                a.loop_id==b.loop_id&&a.reversed==b.reversed&&
                a.next==b.next&&a.prev==b.prev&&a.twin==b.twin)&&
            self.model.loops.iter().zip(&m.loops).all(|(a,b)|
                a.id==b.id&&a.face==b.face&&a.role==b.role&&
                a.first_coedge==b.first_coedge)&&
            self.model.faces.iter().zip(&m.faces).all(|(a,b)|
                a.id==b.id&&a.loops==b.loops&&a.sense==b.sense&&
                a.geometry==b.geometry)&&
            self.model.shells.iter().zip(&m.shells).all(|(a,b)|
                a.id==b.id&&a.faces==b.faces&&a.closed==b.closed)&&
            self.model.geometry.curves3==m.geometry.curves3&&
            self.model.geometry.curves2==m.geometry.curves2&&
            self.model.geometry.surfaces==m.geometry.surfaces&&
            self.model.geometry.edges==m.geometry.edges&&
            self.model.geometry.coedges==m.geometry.coedges&&
            self.slab_plan==p;
        if !matching{return Err(invalid("composed shell differs from source replay"));}
        Ok(())
    }
    fn build(source:&Solid,p:&ComposedCrossingParameters)
        ->Result<(BrepModel,MultiCutRegionPlan),GeometryError>{
        crate::boolean::certify_box(source)
            .map_err(|e|invalid(format!("source is not a certified box: {e}")))?;
        let n=2*p.pairs.len()+p.isolated_blind.len();
        if p.pairs.is_empty()||p.pairs.len()>4||n>12{
            return Err(invalid("requires 1..=4 crossing pairs and <=12 cutters"));
        }
        let b=source.bbox;let span=b.max.sub(b.min);
        let eps=source.tolerance.length_at(span.x.max(span.y).max(span.z));
        certify_group_clearance(p,eps)?;
        let through=p.pairs.iter().map(|q|q.through).collect::<Vec<_>>();
        let blind=p.pairs.iter().map(|q|q.blind)
            .chain(p.isolated_blind.iter().copied()).collect::<Vec<_>>();
        let plan=MultiCutRegionPlanner::plan(source,&through,&blind)?;
        if plan.state!=MultiCutPlanState::Ready {
            return Err(invalid("multi-depth source arrangement is held/unresolved"));
        }
        plan.validate_against(source,&through,&blind)?;
        let mut model=BrepModel::from_polyhedral(source,crate::ModelUnits::default())?;
        let v0=source.mass.volume;
        let c0=source.mass.centroid;
        let mut removed=0.0;
        let mut extra_area=0.0;
        let mut relative=Point3{x:0.0,y:0.0,z:0.0};
        for pair in &p.pairs{
            let donor=MixedCrossingBrep::from_box(source,pair.through,pair.blind)?;
            donor.validate()?;
            graft(&mut model,donor.model(),true)?;
            let m=donor.model().mass;
            removed+=v0-m.volume;
            extra_area+=m.surface_area-source.mass.surface_area;
            relative=relative.add(m.centroid.sub(c0).scale(m.volume));
        }
        for hole in &p.isolated_blind{
            let donor=BlindHoleBrep::from_box(source,*hole)?;
            donor.validate()?;
            graft(&mut model,donor.model(),false)?;
            let m=donor.model().mass;
            removed+=v0-m.volume;
            extra_area+=m.surface_area-source.mass.surface_area;
            relative=relative.add(m.centroid.sub(c0).scale(m.volume));
        }
        let volume=v0-removed;
        let independent=plan.predicted_remaining_volume
            .ok_or_else(||invalid("missing independent slab volume"))?;
        let budget=eps*span.x.max(span.y).max(span.z).powi(2)*256.0+
            v0.abs()*2048.0*f64::EPSILON;
        if !volume.is_finite()||volume<=0.0||
            (volume-independent).abs()>budget{
            return Err(invalid("cross-layer volume differs from independent slab integral"));
        }
        model.mass=MassProperties{
            volume,surface_area:source.mass.surface_area+extra_area,
            centroid:c0.add(relative.scale(1.0/volume))
        };
        model.origin=BrepOrigin::AnalyticMixedThroughBlindRegions;
        model.validate()?;
        let genus=p.pairs.len();
        let s=model.summary();
        if s.shells!=1||s.genus!=genus||
           s.euler_characteristic!=2-2*genus as isize||
           model.faces[0].loops.len()!=1+genus||
           model.faces[1].loops.len()!=1+genus+p.isolated_blind.len()||
           model.loops.iter().filter(|lp|lp.role==LoopRole::Inner).count()!=
               2*genus+p.isolated_blind.len()||
           model.faces.len()!=6+5*genus+2*p.isolated_blind.len(){
            return Err(invalid("joined crossing shell topology/genus mismatch"));
        }
        Ok((model,plan))
    }
    pub fn begin_edit(&mut self)->ComposedCrossingTransaction<'_>{
        let staged=self.parameters.clone();
        ComposedCrossingTransaction{owner:self,staged,failed:false,changes:0}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let history=if redo{&self.redo}else{&self.undo};
        let Some(candidate)=history.last()else{return Ok(false)};
        let staged=candidate.clone();
        let (model,plan)=Self::build(&self.source,&staged)?;
        let rev=self.revision.checked_add(1)
            .ok_or_else(||invalid("source revision overflow"))?;
        if redo{self.redo.pop();self.undo.push(self.parameters.clone());}
        else{self.undo.pop();self.redo.push(self.parameters.clone());}
        self.parameters=staged;self.model=model;self.slab_plan=plan;self.revision=rev;
        Ok(true)
    }
}
impl BrepSource for ComposedCrossingBrep{
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,n:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(n)?.mesh)
    }
}
pub struct ComposedCrossingTransaction<'a>{
    owner:&'a mut ComposedCrossingBrep,
    staged:ComposedCrossingParameters,
    failed:bool,
    changes:usize,
}
impl ComposedCrossingTransaction<'_>{
    pub fn replace_pair(&mut self,index:usize,pair:CertifiedCrossingPair)
        ->Result<(),GeometryError>{
        if self.failed{return Err(invalid("transaction is poisoned"));}
        if index>=self.staged.pairs.len()||self.changes>=64||
           !pair.through.center_x.is_finite()||
           !pair.through.center_z.is_finite()||
           !pair.through.radius.is_finite()||
           !pair.blind.center_x.is_finite()||
           !pair.blind.center_z.is_finite()||
           !pair.blind.radius.is_finite()||
           !pair.blind.depth.is_finite()||
           pair.through.radius<=0.0||pair.blind.radius<=0.0||
           pair.blind.depth<=0.0{
            self.failed=true;return Err(invalid("invalid pair edit/index"));
        }
        self.staged.pairs[index]=pair;self.changes+=1;Ok(())
    }
    pub fn replace_isolated_blind(&mut self,index:usize,hole:BlindHoleParameters)
        ->Result<(),GeometryError>{
        if self.failed{return Err(invalid("transaction is poisoned"));}
        if index>=self.staged.isolated_blind.len()||self.changes>=64||
           !hole.center_x.is_finite()||
           !hole.center_z.is_finite()||
           !hole.radius.is_finite()||
           !hole.depth.is_finite()||
           hole.radius<=0.0||hole.depth<=0.0{
            self.failed=true;return Err(invalid("invalid blind edit/index"));
        }
        self.staged.isolated_blind[index]=hole;self.changes+=1;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.failed{return Err(invalid("cannot commit poisoned source transaction"));}
        if self.changes==0{return Err(invalid("empty source edit batch"));}
        let (model,plan)=ComposedCrossingBrep::build(
            &self.owner.source,&self.staged)?;
        let revision=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("source revision overflow"))?;
        if self.owner.undo.len()>=MAX_HISTORY{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters.clone());
        self.owner.redo.clear();
        self.owner.parameters=self.staged;
        self.owner.model=model;
        self.owner.slab_plan=plan;
        self.owner.revision=revision;
        Ok(revision)
    }
}
