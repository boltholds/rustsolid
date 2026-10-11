//! Canonical reconstruction for a bounded family of true analytic genus-N
//! solids: a certified box with pairwise-disjoint, fully penetrating circular
//! cylinders. Insertion creates real vertices, edges, coedges, cap inner
//! loops, periodic cylindrical faces and 2D/3D carriers, never faceted CSG.
//! Every candidate is validated as a complete B-rep before publication.
use crate::{
    BrepFace, BrepModel, BrepOrigin, Circle2, Circle3, Coedge, CoedgeGeometry,
    CoedgeId, Curve2, Curve2Id, Curve3, Curve3Id, CylinderSurface, Edge,
    EdgeGeometry, EdgeId, FaceGeometry, FaceId, Frame3, GeometryError,
    HoleParameters, Line2, Line3, Loop, LoopId, LoopRole, MassProperties,
    ModelUnits, ParameterDomain2, ParameterRange, Point2Param, Point3,
    Solid, Surface3, SurfaceId, ThroughHoleBrep, Vertex, VertexId,
};
use std::f64::consts::{PI,TAU};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn uv(u:f64,v:f64)->Point2Param{Point2Param{u,v}}
fn uv_line(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,delta:uv(b.u-a.u,b.v-a.v)})
}
fn next_id(len:usize)->Result<u32,GeometryError>{
    u32::try_from(len).map_err(|_|invalid("B-rep topology or carrier capacity exhausted"))
}

/// All holes must lie strictly inside the four side planes. Zero spacing,
/// touching/tangent cylinders, invalid or sub-tolerance features are held.
pub(crate) fn check_hole_set(block:&Solid,holes:&[HoleParameters])
    ->Result<(),GeometryError>
{
    crate::boolean::certify_box(block).map_err(|e|invalid(e.to_string()))?;
    if holes.len()>32{return Err(invalid("maximum 32 independent through-hole features"));}
    let b=block.bbox;
    let dx=b.max.x-b.min.x;
    let dy=b.max.y-b.min.y;
    let dz=b.max.z-b.min.z;
    let eps=block.tolerance.length_at(dx.max(dy).max(dz));
    for (i,h) in holes.iter().enumerate(){
        if !h.center_x.is_finite() || !h.center_z.is_finite() ||
           !h.radius.is_finite() || h.radius<=eps {
            return Err(invalid("invalid or sub-tolerance hole parameters"));
        }
        if h.center_x-h.radius<=b.min.x+eps || h.center_x+h.radius>=b.max.x-eps ||
           h.center_z-h.radius<=b.min.z+eps || h.center_z+h.radius>=b.max.z-eps {
            return Err(invalid("hole must remain strictly inside all four box side planes"));
        }
        for other in &holes[..i] {
            let separation=(h.center_x-other.center_x)
                .hypot(h.center_z-other.center_z);
            if !separation.is_finite() || separation<=h.radius+other.radius+eps{
                return Err(invalid("holes intersect, touch, or lack provable clearance"));
            }
        }
    }
    Ok(())
}

/// Index origin of a hole within the densely reconstructed topology.
/// The first six box faces, eight box vertices and twelve box edges retain
/// their topological IDs. A feature's slot ID is *not* its dense order.
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub(crate) struct HoleTopology {
    pub bottom_edge:EdgeId,
    pub top_edge:EdgeId,
    pub seam_edge:EdgeId,
    pub bottom_trim:CoedgeId,
    pub top_trim:CoedgeId,
    pub bottom_side_trim:CoedgeId,
    pub top_side_trim:CoedgeId,
    pub seam_positive:CoedgeId,
    pub seam_negative:CoedgeId,
    pub side_face:FaceId,
    pub bottom_loop:LoopId,
    pub top_loop:LoopId,
    pub side_loop:LoopId,
}
pub(crate) fn topology_for(ordinal:usize)->Result<HoleTopology,GeometryError>{
    let edge=next_id(12+ordinal*3)?;
    let coedge=next_id(24+ordinal*6)?;
    let face=next_id(6+ordinal)?;
    let lp=next_id(6+ordinal*3)?;
    Ok(HoleTopology{
        bottom_edge:EdgeId(edge),top_edge:EdgeId(edge+1),
        seam_edge:EdgeId(edge+2),
        bottom_trim:CoedgeId(coedge),top_trim:CoedgeId(coedge+1),
        bottom_side_trim:CoedgeId(coedge+2),top_side_trim:CoedgeId(coedge+3),
        seam_positive:CoedgeId(coedge+4),seam_negative:CoedgeId(coedge+5),
        side_face:FaceId(face),
        bottom_loop:LoopId(lp),top_loop:LoopId(lp+1),side_loop:LoopId(lp+2),
    })
}
/// A verified whole-body candidate, not a partial topology delta. Every
/// newly added hole changes the genus and analytic mass of the result.
pub(crate) fn build_multiple(block:&Solid,holes:&[HoleParameters])
    ->Result<BrepModel,GeometryError>
{
    check_hole_set(block,holes)?;
    if holes.is_empty(){return BrepModel::from_polyhedral(block,ModelUnits::default());}
    let first=ThroughHoleBrep::from_box(block,holes[0])?;
    let mut model=first.model().clone();
    for (ordinal,&params) in holes.iter().enumerate().skip(1){
        append_hole(&mut model,block,params,ordinal)?;
    }
    let height=block.bbox.max.y-block.bbox.min.y;
    let original=block.mass.centroid;
    let mut removed=0.0;
    let mut correction=Point3{x:0.0,y:0.0,z:0.0};
    let mut exposed_area=block.mass.surface_area;
    for h in holes {
        let dv=PI*h.radius*h.radius*height;
        removed+=dv;
        exposed_area+=TAU*h.radius*height-2.0*PI*h.radius*h.radius;
        correction=correction.add(Point3{
            x:(original.x-h.center_x)*dv,y:0.0,
            z:(original.z-h.center_z)*dv,
        });
    }
    let volume=block.mass.volume-removed;
    if !volume.is_finite() || volume<=0.0 ||
       !exposed_area.is_finite() || exposed_area<=0.0 {
        return Err(invalid("multi-hole analytic mass is invalid"));
    }
    model.mass=MassProperties{
        volume,surface_area:exposed_area,
        centroid:original.add(correction.scale(1.0/volume)),
    };
    model.origin=if holes.len()==1{
        BrepOrigin::AnalyticThroughHole
    }else{BrepOrigin::AnalyticMultipleThroughHoles};
    model.validate()?;
    if model.summary().genus!=holes.len(){
        return Err(invalid("multi-hole genus does not match independent feature count"));
    }
    Ok(model)
}
fn append_hole(model:&mut BrepModel,block:&Solid,params:HoleParameters,
               ordinal:usize)->Result<(),GeometryError>
{
    let ids=topology_for(ordinal)?;
    if model.vertices.len()!=8+2*ordinal || model.edges.len()!=12+3*ordinal ||
       model.coedges.len()!=24+6*ordinal || model.faces.len()!=6+ordinal ||
       model.loops.len()!=6+3*ordinal {
        return Err(invalid("existing B-rep does not match canonical hole insertion offsets"));
    }
    let tol=block.tolerance;
    let height=block.bbox.max.y-block.bbox.min.y;
    let origin=Point3{x:params.center_x,y:block.bbox.min.y,z:params.center_z};
    let frame=Frame3::from_axes(origin,
        Point3{x:1.0,y:0.0,z:0.0},
        Point3{x:0.0,y:0.0,z:-1.0},tol)?;
    let top_origin=Point3{y:block.bbox.max.y,..origin};
    let top_frame=Frame3{origin:top_origin,..frame};
    let bottom_seam=Point3{x:params.center_x+params.radius,
        y:block.bbox.min.y,z:params.center_z};
    let upper_seam=Point3{y:block.bbox.max.y,..bottom_seam};
    let vid=next_id(model.vertices.len())?;
    model.vertices.push(Vertex{id:VertexId(vid),position:bottom_seam});
    model.vertices.push(Vertex{id:VertexId(vid+1),position:upper_seam});
    let cg=&mut model.geometry;
    let curve_lower=Curve3Id(next_id(cg.curves3.len())?);
    cg.curves3.push(Curve3::Circle(Circle3::new(frame,params.radius,tol)?));
    let curve_upper=Curve3Id(next_id(cg.curves3.len())?);
    cg.curves3.push(Curve3::Circle(Circle3::new(top_frame,params.radius,tol)?));
    let curve_seam=Curve3Id(next_id(cg.curves3.len())?);
    cg.curves3.push(Curve3::Line(Line3{origin:bottom_seam,direction:frame.normal}));
    let eb=ids.bottom_edge;let et=ids.top_edge;let es=ids.seam_edge;
    let cb=ids.bottom_trim;let ct=ids.top_trim;
    let cbs=ids.bottom_side_trim;let cts=ids.top_side_trim;
    let csp=ids.seam_positive;let csn=ids.seam_negative;
    let side=ids.side_face;
    model.edges.extend([
        Edge{id:eb,start:VertexId(vid),end:VertexId(vid),coedges:[cbs,cb]},
        Edge{id:et,start:VertexId(vid+1),end:VertexId(vid+1),coedges:[cts,ct]},
        Edge{id:es,start:VertexId(vid),end:VertexId(vid+1),coedges:[csp,csn]},
    ]);
    let full=ParameterRange::Bounded{start:0.0,end:TAU};
    let axial=ParameterRange::Bounded{start:0.0,end:height};
    cg.edges.extend([
        EdgeGeometry{curve:curve_lower,domain:full,tolerance:tol},
        EdgeGeometry{curve:curve_upper,domain:full,tolerance:tol},
        EdgeGeometry{curve:curve_seam,domain:axial,tolerance:tol},
    ]);
    let support=SurfaceId(next_id(cg.surfaces.len())?);
    cg.surfaces.push(Surface3::Cylinder(CylinderSurface::new(
        frame,params.radius,tol)?));
    let bottom_plane=match cg.surfaces[model.faces[0].geometry.surface.0 as usize]{
        Surface3::Plane(p)=>p,_=>return Err(invalid("bottom support is not planar"))
    };
    let top_plane=match cg.surfaces[model.faces[1].geometry.surface.0 as usize]{
        Surface3::Plane(p)=>p,_=>return Err(invalid("top support is not planar"))
    };
    let p_lower=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(Curve2::Circle(Circle2::new(
        bottom_plane.frame.project(origin)?,params.radius,tol)?));
    let p_upper=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(Curve2::Circle(Circle2::new(
        top_plane.frame.project(top_origin)?,params.radius,tol)?));
    let p_bottom=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(uv_line(uv(0.0,0.0),uv(TAU,0.0)));
    let p_top=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(uv_line(uv(TAU,height),uv(0.0,height)));
    let p_seam1=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(uv_line(uv(TAU,0.0),uv(TAU,height)));
    let p_seam2=Curve2Id(next_id(cg.curves2.len())?);
    cg.curves2.push(uv_line(uv(0.0,height),uv(0.0,0.0)));
    let lb=ids.bottom_loop;let lt=ids.top_loop;let ls=ids.side_loop;
    model.coedges.extend([
        Coedge{id:cb,edge:eb,face:FaceId(0),loop_id:lb,reversed:true,
            next:cb,prev:cb,twin:cbs},
        Coedge{id:ct,edge:et,face:FaceId(1),loop_id:lt,reversed:false,
            next:ct,prev:ct,twin:cts},
        Coedge{id:cbs,edge:eb,face:side,loop_id:ls,reversed:false,
            next:csp,prev:csn,twin:cb},
        Coedge{id:cts,edge:et,face:side,loop_id:ls,reversed:true,
            next:csn,prev:csp,twin:ct},
        Coedge{id:csp,edge:es,face:side,loop_id:ls,reversed:false,
            next:cts,prev:cbs,twin:csn},
        Coedge{id:csn,edge:es,face:side,loop_id:ls,reversed:true,
            next:cbs,prev:cts,twin:csp},
    ]);
    let unit=ParameterRange::Bounded{start:0.0,end:1.0};
    cg.coedges.extend([
        CoedgeGeometry{pcurve:p_lower,domain:full,face:FaceId(0),loop_id:lb},
        CoedgeGeometry{pcurve:p_upper,domain:full,face:FaceId(1),loop_id:lt},
        CoedgeGeometry{pcurve:p_bottom,domain:unit,face:side,loop_id:ls},
        CoedgeGeometry{pcurve:p_top,domain:unit,face:side,loop_id:ls},
        CoedgeGeometry{pcurve:p_seam1,domain:unit,face:side,loop_id:ls},
        CoedgeGeometry{pcurve:p_seam2,domain:unit,face:side,loop_id:ls},
    ]);
    model.loops.extend([
        Loop{id:lb,face:FaceId(0),role:LoopRole::Inner,first_coedge:cb},
        Loop{id:lt,face:FaceId(1),role:LoopRole::Inner,first_coedge:ct},
        Loop{id:ls,face:side,role:LoopRole::Outer,first_coedge:cbs},
    ]);
    model.faces[0].loops.push(lb);
    model.faces[1].loops.push(lt);
    model.faces.push(BrepFace{id:side,loops:vec![ls],geometry:FaceGeometry{
        surface:support,outer_loop:ls,tolerance:tol,bounds:ParameterDomain2{
            u:full,v:axial,u_periodic:true,v_periodic:false,
        },
    }});
    model.shells[0].faces.push(side);
    Ok(())
}
