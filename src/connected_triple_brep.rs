//! Exactly one source-certified connected triple T-A-B chain as a real
//! cylindrical/planar canonical B-rep, built from v0.34 global 3D edge uses.
//!
//! T is full height; A is a deeper blind circle crossing T; B is a shallower
//! blind circle crossing A but not T. Only strict transverse intersections
//! with independently separated T/B are executable. A three-way intersection,
//! tangency or coalescing floor levels is rejected before topology mutation.
//!
//! The B-rep is analytical: 3D Circle3 and Line3 carriers, true plane/cylinder
//! UV trims, two-use reciprocal coedges, closed loops and correct material
//! sense. It is NOT a general curved Boolean operator or faceted CAD model.
use std::collections::{BTreeMap,BTreeSet};
use std::f64::consts::TAU;
use crate::{
    BrepFace,BrepModel,BrepOrigin,BrepSource,Circle3,Curve3,CylinderSurface,
    EdgeId,FaceGeometry,FaceId,FaceSense,Frame3,GeometryError,
    LoopId,LoopRole,MassProperties,Mesh,ModelUnits,ParameterDomain2,
    ParameterRange,PlaneSurface,Point3,Solid,Surface3,VertexId,
    ConnectedTripleStitchPlanner,TripleChainParameters,TripleStitchPlan,
    TripleStitchState,TriplePhysicalEdgeKey,TripleContour,TripleWallPatch,
    TripleArc,
};
use crate::mixed_crossing_brep::{
    add_circle,add_surface,add_vertex,add_vertical,append_loop,twins,
};

const HISTORY_LIMIT:usize=128;
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn bounded(a:f64,b:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(a,b)
}
fn sum_compensated(values:impl IntoIterator<Item=f64>)->f64{
    let (mut out,mut err)=(0.0,0.0);
    for v in values {
        let adjusted=v-err;
        let tmp=out+adjusted;
        err=(tmp-out)-adjusted;
        out=tmp;
    }
    out
}
#[derive(Debug,Clone,Copy)]
struct Circle{
    x:f64,z:f64,r:f64
}
fn circle_sources(params:TripleChainParameters)->[Circle;3]{
    [
        Circle{x:params.through.center_x,z:params.through.center_z,
            r:params.through.radius},
        Circle{x:params.deeper_blind.center_x,z:params.deeper_blind.center_z,
            r:params.deeper_blind.radius},
        Circle{x:params.shallower_blind.center_x,z:params.shallower_blind.center_z,
            r:params.shallower_blind.radius},
    ]
}
/// Green's line integral of the planar FIRST MOMENT about an independent
/// coordinate origin, derived from 1/2 integral x^2 dz and -1/2 integral
/// z^2 dx. This does not use display facets or mass from donors.
fn first_moments(circle:Circle,a:f64,b:f64,origin:(f64,f64))->(f64,f64){
    let x=circle.x-origin.0;
    let z=circle.z-origin.1;
    let r=circle.r;
    let delta=b-a;
    let sin1=b.sin()-a.sin();
    let sin2=(2.0*b).sin()-(2.0*a).sin();
    let cos1=a.cos()-b.cos();
    let cos3=b.cos().powi(3)-a.cos().powi(3);
    let sin3=b.sin().powi(3)-a.sin().powi(3);
    let cos_sq=0.5*delta+0.25*sin2;
    let sin_sq=0.5*delta-0.25*sin2;
    let cos_cub=sin1-sin3/3.0;
    let sin_cub=cos1+cos3/3.0;
    let mx=0.5*r*(x*x*sin1+2.0*x*r*cos_sq+r*r*cos_cub);
    let mz=0.5*r*(z*z*cos1+2.0*z*r*sin_sq+r*r*sin_cub);
    (mx,mz)
}
fn orient_world_axes(top:f64,x:f64,z:f64,
                     tolerance:crate::GeometryTolerance)->Result<Frame3,GeometryError>{
    Frame3::from_axes(
        Point3{x,y:top,z},
        Point3{x:1.0,y:0.0,z:0.0},
        Point3{x:0.0,y:0.0,z:1.0},tolerance
    )
}
/// Planar source contours are oriented CCW for Green area. For the two
/// cap inner loops, top preserves CCW and bottom reverses. For an internal
/// blind floor (material normal +Y), the physical B-rep loop reverses the
/// CCW difference contour, like the single-pair v0.31 analytical donor.
fn contour_parts(contour:&TripleContour,level:usize,
                 reverse:bool,edges:&BTreeMap<TriplePhysicalEdgeKey,EdgeId>)
    ->Result<Vec<(EdgeId,bool)>,GeometryError>{
    let ordered=if reverse{
        contour.parts.iter().rev().collect::<Vec<_>>()
    }else{
        contour.parts.iter().collect::<Vec<_>>()
    };
    ordered.into_iter().map(|p|{
        let id=*edges.get(&TriplePhysicalEdgeKey::CircleArc{
            arc:p.arc,level
        }).ok_or_else(||invalid("source contour references absent physical edge"))?;
        Ok((id,p.reversed!=reverse))
    }).collect()
}
fn patch_parts(wall:&TripleWallPatch,atoms:&[TripleArc],
               edges:&BTreeMap<TriplePhysicalEdgeKey,EdgeId>)
    ->Result<[(EdgeId,bool);4],GeometryError>{
    let atom=&atoms[wall.arc];
    let get=|key:TriplePhysicalEdgeKey|edges.get(&key).copied()
        .ok_or_else(||invalid("missing source wall edge"));
    Ok([
        (get(TriplePhysicalEdgeKey::CircleArc{
            arc:wall.arc,level:wall.slab+1})?,true),
        (get(TriplePhysicalEdgeKey::VerticalTrace{
            junction:atom.start_junction,slab:wall.slab})?,true),
        (get(TriplePhysicalEdgeKey::CircleArc{
            arc:wall.arc,level:wall.slab})?,false),
        (get(TriplePhysicalEdgeKey::VerticalTrace{
            junction:atom.end_junction,slab:wall.slab})?,false),
    ])
}
fn analytic_mass(source:&Solid,p:&TripleStitchPlan)
    ->Result<MassProperties,GeometryError>{
    let remaining=p.predicted_remaining_volume
        .ok_or_else(||invalid("missing ready volume witness"))?;
    let surface=p.analytic_surface_area
        .ok_or_else(||invalid("missing ready area witness"))?;
    let center=source.mass.centroid;
    let disks=circle_sources(p.input);
    let origin=(center.x,center.z);
    let mut moments_x=Vec::new();
    let mut moments_y=Vec::new();
    let mut moments_z=Vec::new();
    for slab in &p.slabs{
        let h=slab.y_upper-slab.y_lower;
        let mut contributions_x=Vec::new();
        let mut contributions_z=Vec::new();
        for edge in &slab.exposed_contour.parts{
            let atom=&p.atomic_arcs[edge.arc];
            let circle=disks[atom.source];
            let (mx,mz)=first_moments(circle,
                atom.start_angle,atom.end_angle,origin);
            contributions_x.push(if edge.reversed{-mx}else{mx});
            contributions_z.push(if edge.reversed{-mz}else{mz});
        }
        moments_x.push(sum_compensated(contributions_x)*h);
        moments_z.push(sum_compensated(contributions_z)*h);
        moments_y.push(slab.union_area*h*
            ((slab.y_lower+slab.y_upper)*0.5-center.y));
    }
    let delta_x=sum_compensated(moments_x);
    let delta_y=sum_compensated(moments_y);
    let delta_z=sum_compensated(moments_z);
    let centroid=Point3{
        x:center.x-delta_x/remaining,
        y:center.y-delta_y/remaining,
        z:center.z-delta_z/remaining,
    };
    if !centroid.is_finite()||!surface.is_finite()||
        !remaining.is_finite()||remaining<=0.0||surface<=0.0{
        return Err(invalid("triple source mass integration is nonfinite"));
    }
    Ok(MassProperties{volume:remaining,
        surface_area:surface,centroid})
}
#[derive(Debug,Clone)]
pub struct ConnectedTripleBrep {
    source:Solid,
    inputs:TripleChainParameters,
    model:BrepModel,
    stitch:TripleStitchPlan,
    revision:u64,
    undo:Vec<TripleChainParameters>,
    redo:Vec<TripleChainParameters>,
}
impl ConnectedTripleBrep {
    pub fn from_box(source:&Solid,inputs:TripleChainParameters)
        ->Result<Self,GeometryError>{
        let (model,stitch)=Self::build(source,inputs)?;
        Ok(Self{source:source.clone(),inputs,model,stitch,revision:0,
            undo:Vec::new(),redo:Vec::new()})
    }
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn stitch(&self)->&TripleStitchPlan{&self.stitch}
    pub fn inputs(&self)->TripleChainParameters{self.inputs}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub(crate) fn source_box(&self)->&Solid{&self.source}
    pub fn tessellate(&self,segments:usize)
        ->Result<crate::ConnectedTripleDisplay,GeometryError>{
        crate::ConnectedTripleMesher::tessellate(self,segments)
    }
    pub fn validate(&self)->Result<(),GeometryError>{
        self.model.validate()?;
        let (expected,stitch)=Self::build(&self.source,self.inputs)?;
        let same=self.stitch==stitch&&
            self.model.summary()==expected.summary()&&
            self.model.origin==expected.origin&&
            self.model.mass==expected.mass&&
            self.model.bbox==expected.bbox&&
            self.model.tolerance==expected.tolerance&&
            self.model.geometry.curves3==expected.geometry.curves3&&
            self.model.geometry.curves2==expected.geometry.curves2&&
            self.model.geometry.surfaces==expected.geometry.surfaces&&
            self.model.geometry.edges==expected.geometry.edges&&
            self.model.geometry.coedges==expected.geometry.coedges&&
            self.model.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>()==
                expected.vertices.iter().map(|v|(v.id,v.position)).collect::<Vec<_>>()&&
            self.model.edges.iter().zip(&expected.edges).all(|(a,b)|
                a.id==b.id&&a.start==b.start&&a.end==b.end&&a.coedges==b.coedges)&&
            self.model.coedges.iter().zip(&expected.coedges).all(|(a,b)|
                a.id==b.id&&a.edge==b.edge&&a.face==b.face&&
                a.loop_id==b.loop_id&&a.reversed==b.reversed&&
                a.next==b.next&&a.prev==b.prev&&a.twin==b.twin)&&
            self.model.loops.iter().zip(&expected.loops).all(|(a,b)|
                a.id==b.id&&a.face==b.face&&a.role==b.role&&
                a.first_coedge==b.first_coedge)&&
            self.model.faces.iter().zip(&expected.faces).all(|(a,b)|
                a.id==b.id&&a.loops==b.loops&&a.sense==b.sense&&
                a.geometry==b.geometry)&&
            self.model.shells.iter().zip(&expected.shells).all(|(a,b)|
                a.id==b.id&&a.faces==b.faces&&a.closed==b.closed);
        if !same{return Err(invalid("triple B-rep topology diverges from source replay"));}
        Ok(())
    }
    fn build(source:&Solid,input:TripleChainParameters)
        ->Result<(BrepModel,TripleStitchPlan),GeometryError>{
        let p=ConnectedTripleStitchPlanner::plan(source,input)?;
        if p.state!=TripleStitchState::Ready ||
            p.predicted_euler!=Some(0)||p.predicted_genus!=Some(1){
            return Err(invalid("cannot commit unvalidated/held triple halfedge plan"));
        }
        p.validate_against(source,input)?;
        let mut model=BrepModel::from_polyhedral(source,ModelUnits::default())?;
        if model.vertices.len()!=8||model.edges.len()!=12||
            model.coedges.len()!=24||model.loops.len()!=6||
            model.faces.len()!=6||model.shells.len()!=1{
            return Err(invalid("connected triple requires canonical source box"));
        }
        let inputs=circle_sources(input);
        let tol=source.tolerance;
        let top=source.bbox.max.y;
        let mut frames=Vec::new();
        for circle in inputs {
            frames.push(orient_world_axes(top,circle.x,circle.z,tol)?);
        }
        let mut used_vertices=BTreeSet::new();
        for e in &p.shared_physical_edges{
            match e.key{
                TriplePhysicalEdgeKey::CircleArc{arc,level}=>{
                    let atom=&p.atomic_arcs[arc];
                    used_vertices.insert((atom.start_junction,level));
                    used_vertices.insert((atom.end_junction,level));
                },
                TriplePhysicalEdgeKey::VerticalTrace{junction,slab}=>{
                    used_vertices.insert((junction,slab));
                    used_vertices.insert((junction,slab+1));
                }
            }
        }
        let mut vertices=BTreeMap::new();
        for &(junction,level) in &used_vertices{
            let point=p.junctions[junction].position_xz;
            let v=add_vertex(&mut model,Point3{
                x:point[0],z:point[1],y:p.levels[level]
            })?;
            vertices.insert((junction,level),v);
        }
        let vert=|junction:usize,level:usize|
            ->Result<VertexId,GeometryError>{
            vertices.get(&(junction,level)).copied().ok_or_else(||
                invalid("physical source junction is missing at height level"))
        };
        let mut edges=BTreeMap::<TriplePhysicalEdgeKey,EdgeId>::new();
        for e in &p.shared_physical_edges{
            let created=match e.key{
                TriplePhysicalEdgeKey::CircleArc{arc,level}=>{
                    let atom=&p.atomic_arcs[arc];
                    let disk=inputs[atom.source];
                    add_circle(&mut model,frames[atom.source],disk.r,
                        p.levels[level],
                        (atom.start_angle,atom.end_angle),
                        vert(atom.start_junction,level)?,
                        vert(atom.end_junction,level)?)?
                },
                TriplePhysicalEdgeKey::VerticalTrace{junction,slab}=>{
                    add_vertical(&mut model,
                        vert(junction,slab)?,vert(junction,slab+1)?)?
                }
            };
            edges.insert(e.key,created);
        }
        let mut side_face_ids=Vec::new();
        for w in &p.walls{
            let disk=inputs[w.source];
            let support=add_surface(&mut model,Surface3::Cylinder(
                CylinderSurface::new(frames[w.source],disk.r,tol)?))?;
            let face=FaceId(u32::try_from(model.faces.len())
                .map_err(|_|invalid("face ID overflow"))?);
            model.faces.push(BrepFace{
                id:face,loops:Vec::new(),sense:FaceSense::Reversed,
                geometry:FaceGeometry{
                    surface:support,outer_loop:LoopId(u32::MAX),tolerance:tol,
                    bounds:ParameterDomain2{
                        u:bounded(w.uv_bounds[0],w.uv_bounds[1])?,
                        v:bounded(w.uv_bounds[2],w.uv_bounds[3])?,
                        u_periodic:false,v_periodic:false,
                    }
                }
            });
            model.shells[0].faces.push(face);
            side_face_ids.push(face);
        }
        let mut floor_face_ids=Vec::new();
        for floor in &p.floors{
            let support=add_surface(&mut model,Surface3::Plane(
                PlaneSurface{frame:orient_world_axes(floor.y,0.0,0.0,tol)?}))?;
            let face=FaceId(u32::try_from(model.faces.len())
                .map_err(|_|invalid("face ID overflow"))?);
            model.faces.push(BrepFace{
                id:face,loops:Vec::new(),sense:FaceSense::Reversed,
                geometry:FaceGeometry{
                    surface:support,outer_loop:LoopId(u32::MAX),tolerance:tol,
                    bounds:ParameterDomain2{
                        u:ParameterRange::Unbounded,
                        v:ParameterRange::Unbounded,
                        u_periodic:false,v_periodic:false,
                    }
                }
            });
            model.shells[0].faces.push(face);
            floor_face_ids.push(face);
        }
        // Caps have one TRUE inner contour each; the physical through hole
        // persists across the bottom and top. Do not fabricate separate
        // top cap faces for individual source circles.
        let bottom=p.bottom_cap.as_ref().ok_or_else(||
            invalid("source bottom cap has no chain contour"))?;
        let top_cap=p.top_cap.as_ref().ok_or_else(||
            invalid("source upper cap has no chain contour"))?;
        append_loop(&mut model,FaceId(0),LoopRole::Inner,
            &contour_parts(bottom,0,true,&edges)?)?;
        append_loop(&mut model,FaceId(1),LoopRole::Inner,
            &contour_parts(top_cap,3,false,&edges)?)?;
        for (floor,face) in p.floors.iter().zip(&floor_face_ids) {
            let level=if floor.y==p.levels[1]{1}else{2};
            append_loop(&mut model,*face,LoopRole::Outer,
                &contour_parts(&floor.contour,level,true,&edges)?)?;
        }
        for (w,face) in p.walls.iter().zip(&side_face_ids){
            let parts=patch_parts(w,&p.atomic_arcs,&edges)?;
            append_loop(&mut model,*face,LoopRole::Outer,&parts)?;
        }
        twins(&mut model,12)?;
        model.mass=analytic_mass(source,&p)?;
        // A closed shell of genus one has chi=0. All physically shared
        // horizontal curves, vertical meeting edges and UV p-curves are
        // checked by the canonical BrepModel validator, not a copied donor.
        model.origin=BrepOrigin::AnalyticMixedThroughBlindRegions;
        model.validate()?;
        let summary=model.summary();
        if summary.vertices!=p.predicted_shell_vertices.unwrap()||
            summary.edges!=p.predicted_shell_edges.unwrap()||
            summary.faces!=p.predicted_shell_faces.unwrap()||
            summary.coedges!=2*summary.edges||
            summary.euler_characteristic!=0||summary.genus!=1||
            summary.shells!=1||model.faces[0].loops.len()!=2||
            model.faces[1].loops.len()!=2||
            model.loops.iter().filter(|lp|lp.role==LoopRole::Inner).count()!=2 ||
            summary.circular_edges!=p.shared_physical_edges.iter()
                .filter(|e|matches!(e.key,
                    TriplePhysicalEdgeKey::CircleArc{..})).count() ||
            summary.cylindrical_faces!=p.walls.len(){
            return Err(invalid("global triple source topology disagrees with edge-use certificate"));
        }
        Ok((model,p))
    }
    pub fn begin_edit(&mut self)->ConnectedTripleTransaction<'_>{
        let staged=self.inputs;
        ConnectedTripleTransaction{owner:self,staged,poisoned:false,changes:0}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{self.transition(false)}
    pub fn redo(&mut self)->Result<bool,GeometryError>{self.transition(true)}
    fn transition(&mut self,redo:bool)->Result<bool,GeometryError>{
        let history=if redo{&self.redo}else{&self.undo};
        let Some(&next)=history.last()else{return Ok(false)};
        let (model,stitch)=Self::build(&self.source,next)?;
        let rev=self.revision.checked_add(1)
            .ok_or_else(||invalid("triple chain revision overflow"))?;
        let before=self.inputs;
        if redo{self.redo.pop();self.undo.push(before);}
        else{self.undo.pop();self.redo.push(before);}
        self.inputs=next;self.model=model;self.stitch=stitch;self.revision=rev;
        Ok(true)
    }
}
impl BrepSource for ConnectedTripleBrep{
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct ConnectedTripleTransaction<'a> {
    owner:&'a mut ConnectedTripleBrep,
    staged:TripleChainParameters,
    poisoned:bool,changes:usize,
}
impl ConnectedTripleTransaction<'_>{
    pub fn replace_sources(&mut self,next:TripleChainParameters)
        ->Result<(),GeometryError>{
        if self.poisoned{return Err(invalid("triple chain edit is poisoned"));}
        let numbers=[
            next.through.center_x,next.through.center_z,next.through.radius,
            next.deeper_blind.center_x,next.deeper_blind.center_z,
            next.deeper_blind.radius,next.deeper_blind.depth,
            next.shallower_blind.center_x,next.shallower_blind.center_z,
            next.shallower_blind.radius,next.shallower_blind.depth
        ];
        if numbers.iter().any(|x|!x.is_finite())||
            next.through.radius<=0.0||next.deeper_blind.radius<=0.0||
            next.shallower_blind.radius<=0.0||
            next.deeper_blind.depth<=0.0||next.shallower_blind.depth<=0.0{
            self.poisoned=true;
            return Err(invalid("invalid source update parameters"));
        }
        if self.changes>=64{
            self.poisoned=true;return Err(invalid("triple source edit budget exhausted"));
        }
        self.staged=next;self.changes+=1;Ok(())
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.poisoned{return Err(invalid("cannot commit poisoned triple source edit"));}
        if self.changes==0{return Err(invalid("empty triple source transaction"));}
        let (model,stitch)=ConnectedTripleBrep::build(&self.owner.source,self.staged)?;
        let rev=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("triple chain revision overflow"))?;
        if self.owner.undo.len()>=HISTORY_LIMIT{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.inputs);
        self.owner.redo.clear();
        self.owner.inputs=self.staged;self.owner.model=model;
        self.owner.stitch=stitch;self.owner.revision=rev;
        Ok(rev)
    }
}
