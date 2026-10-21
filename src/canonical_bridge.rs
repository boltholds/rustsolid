//! Temporary production migration boundary from existing G7/G8 analytic B-reps.
//! Non-authoritative: no volume or Boolean acceptance is implied by this conversion.
//! Delete when G7/G8 assembly emits BrepDraft directly.
use crate as old;
use rustsolid_geometry::ParameterRange;
use rustsolid_topology as topo;
fn invalid(message:impl Into<String>)->old::GeometryError{
    old::GeometryError::InvalidTopology(message.into())
}
fn edge_vertices_extent(model:&old::BrepModel,edge:&old::Edge)->Result<f64,old::GeometryError>{
    let a=model.vertices.get(edge.start.0 as usize).ok_or_else(||invalid("missing start vertex"))?.position;
    let b=model.vertices.get(edge.end.0 as usize).ok_or_else(||invalid("missing end vertex"))?.position;
    let d=a.sub(b);
    Ok(d.x.hypot(d.y).hypot(d.z))
}
fn interval(range:ParameterRange)->Result<(f64,f64),old::GeometryError>{
    match range {
        ParameterRange::Bounded{start,end} if start.is_finite() && end.is_finite() && end>start=>
            Ok((start,end)),
        _=>Err(old::GeometryError::InvalidTopology("source requires bounded edge/trim parameter".into())),
    }
}

fn conversion_for(
    model:&old::BrepModel, coedge:&old::Coedge,
    edge:&old::Edge, new:&topo::BrepEdge,
)->Result<topo::AffineParameterMap,old::GeometryError>{
    let id=edge.id.0 as usize;
    let cid=coedge.id.0 as usize;
    let (t0,t1)=interval(model.geometry.edges[id].domain)?;
    let (u0,u1)=interval(model.geometry.coedges[cid].domain)?;
    let curve=model.edge_curve(edge.id).ok_or_else(||invalid("missing 3D curve"))?;
    let chart=model.coedge_pcurve(coedge.id).ok_or_else(||invalid("missing UV trim"))?;
    let surface=model.face_surface(coedge.face).ok_or_else(||invalid("missing support surface"))?;
    let slope=(u1-u0)/(t1-t0);
    let tries=[
        topo::AffineParameterMap{scale:slope,offset:u0-slope*t0},
        topo::AffineParameterMap{scale:-slope,offset:u1+slope*t0},
    ];
    let mut winner=None;
    for map in tries {
        let mut peak=0.0_f64;
        for fraction in [0.0,0.25,0.5,0.75,1.0] {
            let t=t0+(t1-t0)*fraction;
            let xyz=curve.evaluate(t)?;
            let mapped=map.map(t).map_err(|e|invalid(format!("{e}")))?;
            let uv=chart.evaluate(mapped)?;
            let back=surface.evaluate(uv.u,uv.v)?.point;
            let delta=xyz.sub(back);
            let residual=delta.x.hypot(delta.y).hypot(delta.z);
            if !residual.is_finite(){return Err(invalid("nonfinite source UV residual"));}
            peak=peak.max(residual);
        }
        if winner.map_or(true,|(_,best)|peak<best){winner=Some((map,peak));}
    }
    let (map,peak)=winner.ok_or_else(||invalid("no curve mapping"))?;
    let local_extent=edge_vertices_extent(model,edge)?;
    if peak>model.tolerance.length_at(local_extent.max(1.0))*32.0 {
        return Err(invalid(format!("source coedge {:?} does not match 3D curve",coedge.id)));
    }
    if new.id!=edge.id{return Err(invalid("incorrect physical edge mapping"));}
    Ok(map)
}

pub(crate) fn from_existing(model:&old::BrepModel)->Result<topo::BrepModel,old::GeometryError> {
    model.validate()?;
    let vertices=model.vertices.iter().map(|v|topo::BrepVertex{
        id:v.id,point:v.position,
    }).collect();
    let edges=model.edges.iter().map(|e|{
        topo::BrepEdge{
            id:e.id,vertices:[e.start,e.end],
            curve:topo::Curve3Id(model.geometry.edges[e.id.0 as usize].curve.0),
            domain:model.geometry.edges[e.id.0 as usize].domain,
            uses:e.coedges,
            kind:if model.coedges[e.coedges[0].0 as usize].face==
                    model.coedges[e.coedges[1].0 as usize].face {
                topo::EdgeKind::PeriodicSeam
            }else{topo::EdgeKind::Physical},
        }
    }).collect::<Vec<_>>();
    let coedges=model.coedges.iter().map(|c|{
        let edge=&model.edges[c.edge.0 as usize];
        let map=conversion_for(model,c,edge,&edges[c.edge.0 as usize])?;
        Ok::<_,old::GeometryError>(topo::BrepCoedge{
            id:c.id,edge:c.edge,face:c.face,loop_id:c.loop_id,
            pcurve:topo::Curve2Id(model.geometry.coedges[c.id.0 as usize].pcurve.0),
            edge_to_pcurve:map,
            reversed:c.reversed,next:c.next,prev:c.prev,twin:c.twin,
        })
    }).collect::<Result<Vec<_>,_>>()?;
    let loops=model.loops.iter().map(|l|topo::BrepLoop{
        id:l.id,face:l.face,role:l.role,first:l.first_coedge,
    }).collect();
    let faces=model.faces.iter().map(|f|topo::BrepFace{
        id:f.id,surface:topo::SurfaceId(f.geometry.surface.0),loops:f.loops.clone(),
        sense:match f.sense{old::FaceSense::Forward=>topo::FaceSense::Forward,
                           old::FaceSense::Reversed=>topo::FaceSense::Reversed},
    }).collect();
    let shells=model.shells.iter().map(|s|topo::BrepShell{
        id:s.id,faces:s.faces.clone(),sense:topo::ShellSense::Outward,
    }).collect();
    let carriers=topo::GeometryCatalog{
        curves3:model.geometry.curves3.clone(),
        curves2:model.geometry.curves2.clone(),
        surfaces:model.geometry.surfaces.clone(),
    };
    topo::BrepDraft{
        incarnation:topo::LocalBrepIncarnation::new(17)
            .map_err(|e|invalid(format!("{e}")))?,
        units:model.units,tolerance:model.tolerance,
        vertices,edges,coedges,loops,faces,shells,carriers
    }.check().map_err(|e|invalid(format!("{e}")))
}

