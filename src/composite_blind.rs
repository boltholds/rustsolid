//! Exact bounded multi-blind-bore composition, *not* polygonal CSG.
//! Each candidate is rebuilt using independently validated single-blind
//! analytic B-reps, then copied into a shared face/loop/carrier graph.
//! Independent closed circular source intersections are checked per bore.
//! Intersecting and near-tangent bore regions fail closed.
use crate::{
    BlindHoleBrep, BlindHoleParameters, BrepModel, BrepOrigin, CoedgeId,
    Curve2Id, Curve3Id, EdgeId, FaceId, GeometryError, LoopId,
    MassProperties, Point3, Solid, SurfaceId, VertexId,
};
use std::f64::consts::PI;
fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}

fn validate_features(block:&Solid,holes:&[BlindHoleParameters])->Result<(),GeometryError>{
    crate::boolean::certify_box(block).map_err(|e|invalid(e.to_string()))?;
    if holes.is_empty()||holes.len()>16 {
        return Err(invalid("multi-blind composition requires 1..=16 closed bores"));
    }
    let b=block.bbox;
    let scale=(b.max.x-b.min.x).max(b.max.y-b.min.y).max(b.max.z-b.min.z);
    let eps=block.tolerance.length_at(scale);
    for (i,h) in holes.iter().enumerate() {
        for old in &holes[..i] {
            let distance=(h.center_x-old.center_x).hypot(h.center_z-old.center_z);
            if !distance.is_finite() || distance<=h.radius+old.radius+eps {
                return Err(invalid("blind bore footprints overlap, touch or lack certified clearance"));
            }
        }
    }
    Ok(())
}

/// Candidate creator for a family of independent top-entry blind cylinders.
/// No display tessellation is used and the returned graph is owned.
pub struct CompositeBlindBrep;
impl CompositeBlindBrep {
    pub fn build(block:&Solid,holes:&[BlindHoleParameters])->Result<BrepModel,GeometryError>{
        validate_features(block,holes)?;
        let mut donors=Vec::with_capacity(holes.len());
        for &h in holes {
            let donor=BlindHoleBrep::from_box(block,h)?;
            donor.opening_trace()?;
            donor.validate()?;
            donors.push(donor);
        }
        let mut model=donors[0].model().clone();
        for (index,donor) in donors.iter().enumerate().skip(1){
            graft(&mut model,donor.model(),index)?;
        }
        let mut removed=0.0;
        let mut moment=Point3{x:0.,y:0.,z:0.};
        let mut area=block.mass.surface_area;
        for h in holes {
            let dv=PI*h.radius*h.radius*h.depth;
            let center=Point3{x:h.center_x,y:block.bbox.max.y-h.depth/2.0,
                z:h.center_z};
            removed+=dv;
            moment=moment.add(block.mass.centroid.sub(center).scale(dv));
            area+=2.0*PI*h.radius*h.depth;
        }
        let volume=block.mass.volume-removed;
        if !volume.is_finite()||volume<=0.0||!area.is_finite() {
            return Err(invalid("nonpositive composite blind-hole mass"));
        }
        model.mass=MassProperties{
            volume,surface_area:area,
            centroid:block.mass.centroid.add(moment.scale(1.0/volume)),
        };
        model.origin=BrepOrigin::AnalyticMultipleBlindHoles;
        model.validate()?;
        let summary=model.summary();
        let n=holes.len();
        if summary.vertices!=8+2*n || summary.edges!=12+3*n ||
           summary.coedges!=24+6*n || summary.faces!=6+2*n ||
           summary.loops!=6+3*n || summary.genus!=0 ||
           summary.euler_characteristic!=2 ||
           model.faces[1].loops.len()!=1+n ||
           model.loops.iter().filter(|lp|lp.role==crate::LoopRole::Inner).count()!=n {
            return Err(invalid("multi-blind B-rep invariant failed: missing floor/trim or genus mismatch"));
        }
        Ok(model)
    }
}
fn graft(m:&mut BrepModel,d:&BrepModel,i:usize)->Result<(),GeometryError>{
    let vertex=u32::try_from(8+2*i).map_err(|_|invalid("vertex index exhausted"))?;
    let edge=u32::try_from(12+3*i).map_err(|_|invalid("edge index exhausted"))?;
    let coedge=u32::try_from(24+6*i).map_err(|_|invalid("coedge index exhausted"))?;
    let lp=u32::try_from(6+3*i).map_err(|_|invalid("loop index exhausted"))?;
    let side=u32::try_from(6+2*i).map_err(|_|invalid("face index exhausted"))?;
    if m.vertices.len()!=vertex as usize||m.edges.len()!=edge as usize||
       m.coedges.len()!=coedge as usize||m.loops.len()!=lp as usize||
       m.faces.len()!=side as usize {
        return Err(invalid("composite graft has inconsistent dense topology offsets"));
    }
    let curve3=u32::try_from(m.geometry.curves3.len())
        .map_err(|_|invalid("curve3 catalog exhausted"))?;
    let curve2=u32::try_from(m.geometry.curves2.len())
        .map_err(|_|invalid("curve2 catalog exhausted"))?;
    let surf=u32::try_from(m.geometry.surfaces.len())
        .map_err(|_|invalid("surface catalog exhausted"))?;
    let cv3base=d.geometry.edges[12].curve.0;
    let cv2base=d.geometry.coedges[24].pcurve.0;
    let sbase=d.faces[6].geometry.surface.0;
    for v in &d.vertices[8..10]{
        let mut x=v.clone();x.id=VertexId(vertex+v.id.0-8);m.vertices.push(x);
    }
    for c in &d.geometry.curves3[cv3base as usize..cv3base as usize+3] {
        m.geometry.curves3.push(*c);
    }
    for c in &d.geometry.curves2[cv2base as usize..cv2base as usize+6]{
        m.geometry.curves2.push(*c);
    }
    for c in &d.geometry.surfaces[sbase as usize..sbase as usize+2]{
        m.geometry.surfaces.push(*c);
    }
    for original in &d.edges[12..15]{
        let mut e=original.clone();e.id=EdgeId(edge+original.id.0-12);
        e.start=VertexId(vertex+original.start.0-8);
        e.end=VertexId(vertex+original.end.0-8);
        e.coedges=e.coedges.map(|c|CoedgeId(coedge+c.0-24));
        m.edges.push(e);
    }
    for g in &d.geometry.edges[12..15]{
        let mut e=*g;e.curve=Curve3Id(curve3+g.curve.0-cv3base);
        m.geometry.edges.push(e);
    }
    for original in &d.coedges[24..30]{
        let mut c=original.clone();c.id=CoedgeId(coedge+c.id.0-24);
        c.edge=EdgeId(edge+c.edge.0-12);
        c.face=if c.face.0==6{FaceId(side)}
            else if c.face.0==7{FaceId(side+1)}
            else {c.face};
        c.loop_id=LoopId(lp+c.loop_id.0-6);
        c.next=CoedgeId(coedge+c.next.0-24);
        c.prev=CoedgeId(coedge+c.prev.0-24);
        c.twin=CoedgeId(coedge+c.twin.0-24);
        m.coedges.push(c);
    }
    for original in &d.geometry.coedges[24..30]{
        let mut c=*original;c.pcurve=Curve2Id(curve2+c.pcurve.0-cv2base);
        c.face=if c.face.0==6{FaceId(side)}
            else if c.face.0==7{FaceId(side+1)}else{c.face};
        c.loop_id=LoopId(lp+c.loop_id.0-6);
        m.geometry.coedges.push(c);
    }
    for original in &d.loops[6..9] {
        let mut l=original.clone();
        l.id=LoopId(lp+l.id.0-6);
        l.face=if l.face.0==6{FaceId(side)}
            else if l.face.0==7{FaceId(side+1)}else{l.face};
        l.first_coedge=CoedgeId(coedge+l.first_coedge.0-24);
        m.loops.push(l);
    }
    for original in &d.faces[6..8] {
        let mut f=original.clone();
        f.id=FaceId(side+f.id.0-6);
        f.loops=f.loops.iter().map(|x|LoopId(lp+x.0-6)).collect();
        f.geometry.outer_loop=LoopId(lp+f.geometry.outer_loop.0-6);
        f.geometry.surface=SurfaceId(surf+f.geometry.surface.0-sbase);
        m.faces.push(f);
    }
    m.faces[1].loops.push(LoopId(lp+1));
    m.shells[0].faces.extend([FaceId(side),FaceId(side+1)]);
    Ok(())
}
