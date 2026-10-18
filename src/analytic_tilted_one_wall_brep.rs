//! G8c: source-owned analytic B-rep for a tilted cylinder intersecting only
//! the MIN-X face of a canonical box. Physical cut curves are EXACTLY the
//! four independently certified G1/G8a source curves, not fitted from mesh.
//! This staged shell cannot be committed: actual flux is checked separately.

use std::collections::{BTreeMap,BTreeSet};
use std::f64::consts::TAU;
use crate::{
    BooleanOperation,BrepBody,BrepFace,BrepModel,BrepOrigin,Coedge,
    CoedgeGeometry,CoedgeId,Curve2,Curve2Id,Curve3,Curve3Id,
    Edge,EdgeGeometry,EdgeId,FaceGeometry,FaceId,FaceSense,GeometryError,
    GeometryTolerance,GlobalArcRef,GlobalEventId,GraphOperand,
    HarmonicCurve2,IndependentTiltedVolumeOracle,IndependentTiltedVolumeResult,
    Line2,Line3,Loop,LoopId,LoopRole,MassProperties,ParameterDomain2,
    ParameterRange,Point2Param,Point3,RadialTopologyCertificate,
    Shell,ShellId,SourceCertifiedCutArc,Surface3,SurfaceId,
    TiltedOneWallSourcePlan,TiltedOneWallSourceResult,Vertex,VertexId,
};
use crate::analytic_face_patch_partition::{
    PatchBuilder,find_source_plane,certify_original_material_edge,
    radial_certificate,
};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidTopology(s.into())}
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{
    let d=a.sub(b);d.x.hypot(d.y).hypot(d.z)
}
fn range(start:f64,end:f64)->Result<ParameterRange,GeometryError>{
    ParameterRange::bounded(start,end)
}
fn ends(span:ParameterRange)->Result<(f64,f64),GeometryError>{
    match span{
        ParameterRange::Bounded{start,end} if start.is_finite()&&end.is_finite()&&start<end=>
            Ok((start,end)),
        _=>Err(invalid("G8c requires bounded analytic cut carrier")),
    }
}
fn vertex_pair(a:usize,b:usize,c:usize,d:usize)->bool{
    (a==c&&b==d)||(a==d&&b==c)
}
fn uvline(a:Point2Param,b:Point2Param)->Curve2{
    Curve2::Line(Line2{origin:a,
        delta:Point2Param{u:b.u-a.u,v:b.v-a.v}})
}
fn reflect(source:Curve2,span:ParameterRange)->Result<Curve2,GeometryError>{
    let (lo,hi)=ends(span)?;
    let offset=lo+hi;
    Ok(match source{
        Curve2::Line(l)=>Curve2::Line(Line2{
            origin:Point2Param{
                u:l.origin.u+l.delta.u*offset,
                v:l.origin.v+l.delta.v*offset,
            },
            delta:Point2Param{u:-l.delta.u,v:-l.delta.v},
        }),
        Curve2::CylinderSection(mut s)=>{
            s.angle_offset+=s.direction*offset;
            s.direction= -s.direction;
            Curve2::CylinderSection(s)
        },
        Curve2::Harmonic(h)=>{
            let (sin,cos)=offset.sin_cos();
            Curve2::Harmonic(HarmonicCurve2{
                center:h.center,
                cosine:Point2Param{
                    u:h.cosine.u*cos+h.sine.u*sin,
                    v:h.cosine.v*cos+h.sine.v*sin,
                },
                sine:Point2Param{
                    u:h.cosine.u*sin-h.sine.u*cos,
                    v:h.cosine.v*sin-h.sine.v*cos,
                },
            })
        },
        _=>return Err(invalid("G8c unsupported source cut pcurve")),
    })
}
fn bounds_for_curve(
    curve:Curve2,span:ParameterRange,
)->Result<([f64;2],[f64;2]),GeometryError>{
    let (lo,hi)=ends(span)?;
    let a=curve.evaluate(lo)?;
    let b=curve.evaluate(hi)?;
    let mut u=[a.u.min(b.u),a.u.max(b.u)];
    let mut v=[a.v.min(b.v),a.v.max(b.v)];
    match curve{
        Curve2::Harmonic(c)=>{
            let ru=c.cosine.u.hypot(c.sine.u);
            let rv=c.cosine.v.hypot(c.sine.v);
            u=[c.center.u-ru,c.center.u+ru];
            v=[c.center.v-rv,c.center.v+rv];
        },
        Curve2::CylinderSection(c)=>{
            let rv=c.axial_cosine.hypot(c.axial_sine);
            v=[c.axial_center-rv,c.axial_center+rv];
        },
        Curve2::Line(_)=>{},
        _=>return Err(invalid("G8c unsupported planar or cylinder pcurve")),
    }
    if ![u[0],u[1],v[0],v[1]].iter().all(|x|x.is_finite()){
        return Err(invalid("G8c nonfinite trim bound"));
    }
    Ok((u,v))
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum TiltedShellHold{
    UnsupportedOperand,InvalidSourceVolume,MissingOriginalMaterialEdge,
    MissingPhysicalG1Arc,DuplicatePhysicalG1Arc,UncertainPhysicalEvent,
    MissingSourceFace,SourceCurveMismatch,InvalidUvChart,
    InvalidMaterialSide,InvalidRadialTopology,UnexpectedShellStructure,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct TiltedCutBinding{
    pub source_arc:GlobalArcRef,
    pub result_edge:EdgeId,
    pub source_endpoints:[GlobalEventId;2],
    pub result_vertices:[VertexId;2],
    pub first_face:FaceId,
    pub second_face:FaceId,
}
#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub struct TiltedFaceProvenance {
    pub result_face:FaceId,
    pub original_face:FaceId,
    pub operand:GraphOperand,
    pub reversed_for_difference:bool,
}
#[derive(Debug,Clone)]
pub struct TiltedOneWallShellPlan{
    model:BrepModel,
    /// Bottom, top, Zmin, Xmax, Zmax, two split MIN-X faces, cylinder.
    pub faces:[FaceId;8],
    pub edges:[EdgeId;18],
    pub source_faces:[TiltedFaceProvenance;8],
    pub physical_cuts:[TiltedCutBinding;4],
    pub radial:RadialTopologyCertificate,
    pub source_volume:IndependentTiltedVolumeOracle,
    pub closed_shell_proven:bool,
    pub boundary_flux_certified:bool,
    pub can_apply:bool,
}
#[derive(Debug,Clone)]
pub enum TiltedOneWallShellResult{
    Ready(TiltedOneWallShellPlan),
    Held(TiltedShellHold),
}
impl TiltedOneWallShellPlan{
    pub fn brep(&self)->&BrepModel{&self.model}
    pub fn validate_against(
        &self,first:&BrepBody,second:&BrepBody,
        tol:GeometryTolerance,
    )->Result<(),GeometryError>{
        let TiltedOneWallSourceResult::Ready(src)=
            TiltedOneWallSourcePlan::build(first,second,
                BooleanOperation::Difference,tol)?
        else{return Err(invalid("G8c sources changed from certified family"));};
        self.source_volume.validate_against(first,second,&src,tol)?;
        let TiltedOneWallShellResult::Ready(rebuilt)=
            Self::build(first,second,&src,&self.source_volume,tol)?
        else{return Err(invalid("G8c source-derived shell no longer closes"));};
        if self.faces!=rebuilt.faces||self.edges!=rebuilt.edges
            ||self.source_faces!=rebuilt.source_faces
            ||self.physical_cuts!=rebuilt.physical_cuts
            ||self.radial!=rebuilt.radial
            ||self.closed_shell_proven!=rebuilt.closed_shell_proven
            ||self.boundary_flux_certified!=rebuilt.boundary_flux_certified
            ||self.can_apply!=rebuilt.can_apply
            ||format!("{:?}",self.model)!=format!("{:?}",rebuilt.model){
            return Err(invalid("G8c physical B-rep differs from rebuilt sources"));
        }
        self.model.validate()
    }
    pub fn build(
        first:&BrepBody,second:&BrepBody,
        source:&TiltedOneWallSourcePlan,
        volume:&IndependentTiltedVolumeOracle,tol:GeometryTolerance,
    )->Result<TiltedOneWallShellResult,GeometryError>{
        use TiltedShellHold as Hold;
        tol.validate()?;
        source.validate_against(first,second,tol)?;
        volume.validate_against(first,second,source,tol)?;
        if !volume.source_bracket_established||volume.result_shell_certified
            ||volume.can_apply ||source.closed_shell_proven||source.can_apply{
            return Ok(TiltedOneWallShellResult::Held(Hold::InvalidSourceVolume));
        }
        let (BrepBody::Polyhedral(block),BrepBody::Cylinder(cyl))=
            (first,second) else{
            return Ok(TiltedOneWallShellResult::Held(Hold::UnsupportedOperand));
        };
        let a=first.shared()?;let b=second.shared()?;
        let bbox=block.bbox;
        let scale=bbox.max.sub(bbox.min).length_squared().sqrt()
            .max(cyl.radius).max(cyl.height).max(1.);
        let eps=tol.length_at(scale)*128.;
        let faces_src=[
            source.box_bottom,source.box_top,
            find_source_plane(&a,2,bbox.min.z,-1.,eps)
                .ok_or_else(||invalid("G8c missing source Zmin face"))?,
            find_source_plane(&a,0,bbox.max.x,1.,eps)
                .ok_or_else(||invalid("G8c missing source Xmax face"))?,
            find_source_plane(&a,2,bbox.max.z,1.,eps)
                .ok_or_else(||invalid("G8c missing source Zmax face"))?,
            source.box_min_x,source.box_min_x,source.cutter_side,
        ];
        let bottom=source.bottom.wall_points;
        let top=source.top.wall_points;
        let (bup,bdown)=if bottom[0].z>bottom[1].z
            {(bottom[0],bottom[1])}else{(bottom[1],bottom[0])};
        let (tup,tdown)=if top[0].z>top[1].z
            {(top[0],top[1])}else{(top[1],top[0])};
        let xmin=bbox.min.x;let xmax=bbox.max.x;
        let ymin=bbox.min.y;let ymax=bbox.max.y;
        let zmin=bbox.min.z;let zmax=bbox.max.z;
        let xyz=[
            p(xmin,ymin,zmin),p(xmax,ymin,zmin),
            p(xmax,ymin,zmax),p(xmin,ymin,zmax),
            bup,bdown,
            p(xmin,ymax,zmin),p(xmax,ymax,zmin),
            p(xmax,ymax,zmax),p(xmin,ymax,zmax),
            tup,tdown,
        ];
        let mut builder=PatchBuilder::new(bbox,tol,a.units);
        builder.model.origin=BrepOrigin::AnalyticTiltedOneWallCylinderDifference;
        for (i,&position) in xyz.iter().enumerate(){
            builder.model.vertices.push(Vertex{
                id:VertexId(i as u32),position,
            });
        }
        let mut pairs=[
            (0,1),(1,2),(2,3),(3,4),(4,5),(5,0),
            (6,7),(7,8),(8,9),(9,10),(10,11),(11,6),
            (0,6),(1,7),(2,8),(3,9),(4,10),(5,11),
        ];
        // Link every new curved edge directly to its G1 source arc.
        // No endpoint welding by proximity: only the four G8a physical
        // source event IDs are allowed to designate result vertices.
        let mut bindings=Vec::new();
        let mut arcs=[None::<SourceCertifiedCutArc>;18];
        let mut events=BTreeMap::<GlobalEventId,VertexId>::new();
        let cut_indices=[4usize,10usize,16usize,17usize];
        let source_vertices=[[4usize,5usize],[10,11],[4,10],[5,11]];
        for source_arc in &source.source_cuts.arcs{
            let Some(ids)=source_arc.physical_endpoints else {
                return Ok(TiltedOneWallShellResult::Held(Hold::UncertainPhysicalEvent));
            };
            let (lo,hi)=ends(source_arc.parameter_domain)?;
            let endpoints=[
                source_arc.carrier.evaluate(lo)?,
                source_arc.carrier.evaluate(hi)?,
            ];
            let mut candidate=None;
            for (i,indices) in source_vertices.iter().enumerate(){
                let a_id=indices[0];let b_id=indices[1];
                if distance(endpoints[0],xyz[a_id])<eps &&
                   distance(endpoints[1],xyz[b_id])<eps {
                    if candidate.replace((cut_indices[i],a_id,b_id)).is_some(){
                        return Ok(TiltedOneWallShellResult::Held(Hold::DuplicatePhysicalG1Arc));
                    }
                }else if distance(endpoints[1],xyz[a_id])<eps &&
                          distance(endpoints[0],xyz[b_id])<eps {
                    if candidate.replace((cut_indices[i],b_id,a_id)).is_some(){
                        return Ok(TiltedOneWallShellResult::Held(Hold::DuplicatePhysicalG1Arc));
                    }
                }
            }
            let Some((index,start,end))=candidate else{
                return Ok(TiltedOneWallShellResult::Held(Hold::MissingPhysicalG1Arc));
            };
            if arcs[index].replace(*source_arc).is_some(){
                return Ok(TiltedOneWallShellResult::Held(Hold::DuplicatePhysicalG1Arc));
            }
            pairs[index]=(start,end);
            for (event,vertex) in [(ids[0],VertexId(start as u32)),
                (ids[1],VertexId(end as u32))]{
                if let Some(previous)=events.insert(event,vertex){
                    if previous!=vertex{
                        return Ok(TiltedOneWallShellResult::Held(Hold::UncertainPhysicalEvent));
                    }
                }
            }
            bindings.push(TiltedCutBinding{
                source_arc:source_arc.source,
                result_edge:EdgeId(index as u32),
                source_endpoints:ids,
                result_vertices:[VertexId(start as u32),VertexId(end as u32)],
                first_face:source_arc.first_face,second_face:source_arc.second_face,
            });
        }
        if bindings.len()!=4 ||cut_indices.iter().any(|&i|arcs[i].is_none())
            ||events.len()!=4{
            return Ok(TiltedOneWallShellResult::Held(Hold::MissingPhysicalG1Arc));
        }
        let mut result_edges=[EdgeId(0);18];
        for (index,&(start,end)) in pairs.iter().enumerate(){
            let id=EdgeId(builder.model.edges.len() as u32);
            if index!=id.0 as usize{
                return Err(invalid("G8c physical edges have noncanonical ordering"));
            }
            let (curve,span)=if let Some(cut)=arcs[index]{
                let (lo,hi)=ends(cut.parameter_domain)?;
                if distance(cut.carrier.evaluate(lo)?,xyz[start])>eps
                    ||distance(cut.carrier.evaluate(hi)?,xyz[end])>eps{
                    return Ok(TiltedOneWallShellResult::Held(Hold::SourceCurveMismatch));
                }
                (cut.carrier,cut.parameter_domain)
            }else{
                let (line,span)=Line3::between(xyz[start],xyz[end],tol)?;
                (Curve3::Line(line),span)
            };
            let curve_id=Curve3Id(builder.model.geometry.curves3.len() as u32);
            builder.model.geometry.curves3.push(curve);
            builder.model.geometry.edges.push(EdgeGeometry{
                curve:curve_id,domain:span,tolerance:tol,
            });
            builder.model.edges.push(Edge{
                id,start:VertexId(start as u32),end:VertexId(end as u32),
                coedges:[CoedgeId(u32::MAX);2],
            });
            result_edges[index]=id;
        }
        let originals=[
            (0,[faces_src[0],faces_src[2]]),
            (1,[faces_src[0],faces_src[3]]),
            (2,[faces_src[0],faces_src[4]]),
            (3,[faces_src[0],faces_src[5]]),
            (5,[faces_src[0],faces_src[5]]),
            (6,[faces_src[1],faces_src[2]]),
            (7,[faces_src[1],faces_src[3]]),
            (8,[faces_src[1],faces_src[4]]),
            (9,[faces_src[1],faces_src[5]]),
            (11,[faces_src[1],faces_src[5]]),
            (12,[faces_src[2],faces_src[5]]),
            (13,[faces_src[2],faces_src[3]]),
            (14,[faces_src[4],faces_src[3]]),
            (15,[faces_src[4],faces_src[5]]),
        ];
        let mut used=BTreeSet::<EdgeId>::new();
        for (idx,source_faces) in originals{
            if !certify_original_material_edge(&a,
                &builder.model.edges[idx],&builder.model.vertices,
                source_faces,eps,&mut used){
                return Ok(TiltedOneWallShellResult::Held(
                    Hold::MissingOriginalMaterialEdge));
            }
        }
        let face_specs=[
            (vec![0,1,2,3,4,5],p(0.,-1.,0.)),
            (vec![6,7,8,9,10,11],p(0.,1.,0.)),
            (vec![0,1,7,6],p(0.,0.,-1.)),
            (vec![1,2,8,7],p(1.,0.,0.)),
            (vec![2,3,9,8],p(0.,0.,1.)),
            (vec![3,4,10,9],p(-1.,0.,0.)),
            (vec![5,0,6,11],p(-1.,0.,0.)),
            (vec![4,5,11,10],p(-1.,0.,0.)),
        ];
        let mut sources=Vec::<TiltedFaceProvenance>::new();
        for (face_index,(initial,expected_normal)) in
            face_specs.into_iter().enumerate()
        {
            let first_operand=face_index!=7;
            let source_surface=if first_operand{
                a.face_surface(faces_src[face_index])
            }else{
                b.face_surface(faces_src[face_index])
            }.ok_or_else(||invalid("G8c source support surface disappeared"))?;
            if matches!(source_surface,Surface3::Cylinder(_))!=(!first_operand){
                return Ok(TiltedOneWallShellResult::Held(Hold::MissingSourceFace));
            }
            let mut local=initial;
            let q0=xyz[local[0]];
            let q1=xyz[local[1]];
            let q2=xyz[local[2]];
            let sign=q1.sub(q0).cross(q2.sub(q1)).dot(expected_normal);
            if !sign.is_finite()||sign.abs()<=eps*eps{
                return Ok(TiltedOneWallShellResult::Held(Hold::InvalidMaterialSide));
            }
            if sign<0. {local.reverse();}
            let fid=FaceId(builder.model.faces.len() as u32);
            let lid=LoopId(builder.model.loops.len() as u32);
            let sid=SurfaceId(builder.model.geometry.surfaces.len() as u32);
            builder.model.geometry.surfaces.push(source_surface);
            // Bounds will be recomputed from analytic pcurves *after*
            // unwinding the cylinder UV chart.
            builder.model.faces.push(BrepFace{
                id:fid,loops:vec![lid],
                sense:if first_operand{FaceSense::Forward}else{FaceSense::Reversed},
                geometry:FaceGeometry{
                    surface:sid,bounds:ParameterDomain2{
                        u:range(-1.,1.)?,v:range(-1.,1.)?,
                        u_periodic:false,v_periodic:false,
                    },
                    outer_loop:lid,tolerance:tol,
                },
            });
            let first_co=CoedgeId(builder.model.coedges.len() as u32);
            builder.model.loops.push(Loop{
                id:lid,face:fid,role:LoopRole::Outer,first_coedge:first_co,
            });
            let count=local.len();
            for i in 0..count{
                let from_id=local[i];
                let to_id=local[(i+1)%count];
                let Some((index,reversed))=pairs.iter().enumerate()
                    .find_map(|(idx,&(x,y))|{
                        if from_id==x&&to_id==y{Some((idx,false))}
                        else if from_id==y&&to_id==x{Some((idx,true))}
                        else{None}
                    }) else{return Err(invalid("G8c loop has no physical source edge"));};
                let co_id=CoedgeId(builder.model.coedges.len() as u32);
                builder.model.coedges.push(Coedge{
                    id:co_id,edge:result_edges[index],face:fid,loop_id:lid,
                    reversed,
                    next:CoedgeId(first_co.0+((i+1)%count) as u32),
                    prev:CoedgeId(first_co.0+((i+count-1)%count) as u32),
                    twin:CoedgeId(u32::MAX),
                });
                let (trim,trim_range)=if let Some(cut)=arcs[index] {
                    let source_curve=if first_operand{
                        if cut.first_face!=faces_src[face_index]{
                            return Ok(TiltedOneWallShellResult::Held(
                                Hold::SourceCurveMismatch));
                        }
                        cut.face_pcurves[0]
                    }else{
                        if cut.second_face!=faces_src[face_index]{
                            return Ok(TiltedOneWallShellResult::Held(
                                Hold::SourceCurveMismatch));
                        }
                        cut.face_pcurves[1]
                    };
                    (if reversed{reflect(source_curve,cut.parameter_domain)?}
                      else{source_curve},cut.parameter_domain)
                }else{
                    let Surface3::Plane(plane)=source_surface else{
                        return Ok(TiltedOneWallShellResult::Held(
                            Hold::SourceCurveMismatch));
                    };
                    let aa=plane.frame.project(xyz[from_id])?;
                    let bb=plane.frame.project(xyz[to_id])?;
                    (uvline(aa,bb),range(0.,1.)?)
                };
                let trim_id=Curve2Id(builder.model.geometry.curves2.len() as u32);
                builder.model.geometry.curves2.push(trim);
                builder.model.geometry.coedges.push(CoedgeGeometry{
                    pcurve:trim_id,domain:trim_range,face:fid,loop_id:lid,
                });
            }
            sources.push(TiltedFaceProvenance{
                result_face:fid,original_face:faces_src[face_index],
                operand:if first_operand{GraphOperand::First}
                    else{GraphOperand::Second},
                reversed_for_difference:!first_operand,
            });
        }
        // Every physical edge occurs in exactly two opposed coedges.
        for edge in &mut builder.model.edges{
            let co=builder.model.coedges.iter().filter(|co|co.edge==edge.id)
                .map(|co|co.id).collect::<Vec<_>>();
            if co.len()!=2 ||
                builder.model.coedges[co[0].0 as usize].reversed==
                builder.model.coedges[co[1].0 as usize].reversed{
                return Ok(TiltedOneWallShellResult::Held(Hold::InvalidRadialTopology));
            }
            builder.model.coedges[co[0].0 as usize].twin=co[1];
            builder.model.coedges[co[1].0 as usize].twin=co[0];
            edge.coedges=[co[0],co[1]];
        }
        builder.model.shells.push(Shell{
            id:ShellId(0),faces:(0..8).map(FaceId).collect(),closed:true,
        });
        // Normalize analytic cylinder seam angles to ONE unwrapped chart.
        // A permissible 2π shift changes no 3D geometry or source provenance.
        let cyl_loop=builder.model.loops[7].first_coedge;
        let mut co_ids=Vec::new();
        let mut cursor=cyl_loop;
        loop {
            co_ids.push(cursor);
            cursor=builder.model.coedges[cursor.0 as usize].next;
            if cursor==cyl_loop {break;}
            if co_ids.len()>4{return Err(invalid("G8c cut wall loop not quadrilateral"));}
        }
        if co_ids.len()!=4{
            return Ok(TiltedOneWallShellResult::Held(Hold::InvalidUvChart));
        }
        let read_uv=|model:&BrepModel,id:CoedgeId,last:bool|->Result<Point2Param,GeometryError>{
            let co=&model.geometry.coedges[id.0 as usize];
            let (lo,hi)=ends(co.domain)?;
            model.coedge_pcurve(id).ok_or_else(||invalid("missing cylinder UV"))?
                .evaluate(if last{hi}else{lo})
        };
        let first_uv=read_uv(&builder.model,co_ids[0],false)?;
        let mut previous=read_uv(&builder.model,co_ids[0],true)?;
        for &cid in co_ids.iter().skip(1){
            let before=read_uv(&builder.model,cid,false)?;
            let turns=((previous.u-before.u)/TAU).round();
            if !turns.is_finite()||turns.abs()>3. {
                return Ok(TiltedOneWallShellResult::Held(Hold::InvalidUvChart));
            }
            let shift=turns*TAU;
            let geom=builder.model.geometry.coedges[cid.0 as usize];
            let curve=&mut builder.model.geometry.curves2[geom.pcurve.0 as usize];
            match curve {
                Curve2::Line(line)=>line.origin.u+=shift,
                Curve2::CylinderSection(section)=>section.angle_offset+=shift,
                _=>return Ok(TiltedOneWallShellResult::Held(Hold::InvalidUvChart)),
            }
            let after=read_uv(&builder.model,cid,false)?;
            if (after.u-previous.u).abs()>eps||
                (after.v-previous.v).abs()>eps{
                return Ok(TiltedOneWallShellResult::Held(Hold::InvalidUvChart));
            }
            previous=read_uv(&builder.model,cid,true)?;
        }
        if (previous.u-first_uv.u).abs()>eps ||
            (previous.v-first_uv.v).abs()>eps {
            return Ok(TiltedOneWallShellResult::Held(Hold::InvalidUvChart));
        }
        // Bound every genuine analytic trim, including interior extrema
        // of elliptic edges (a vertex-only bound would be incorrect).
        for i in 0..builder.model.faces.len(){
            let mut min_u=f64::INFINITY;let mut max_u=f64::NEG_INFINITY;
            let mut min_v=f64::INFINITY;let mut max_v=f64::NEG_INFINITY;
            for co in builder.model.coedges.iter().filter(|c|c.face.0 as usize==i){
                let cg=builder.model.geometry.coedges[co.id.0 as usize];
                let curve=builder.model.geometry.curves2[cg.pcurve.0 as usize];
                let (u,v)=bounds_for_curve(curve,cg.domain)?;
                min_u=min_u.min(u[0]);max_u=max_u.max(u[1]);
                min_v=min_v.min(v[0]);max_v=max_v.max(v[1]);
            }
            let padding=eps.max(1e-8);
            builder.model.faces[i].geometry.bounds=ParameterDomain2{
                u:range(min_u-padding,max_u+padding)?,
                v:range(min_v-padding,max_v+padding)?,
                u_periodic:false,v_periodic:false,
            };
        }
        // Staged read-only BrepModel requires mass fields. Volume is the
        // source interval midpoint; surface area and centroid are explicitly
        // placeholders until an independently certified mass-properties pass.
        // No editor, history or persisted result may consume these proxies.
        builder.model.mass=MassProperties{
            volume:(volume.retained.lower+volume.retained.upper)*0.5,
            surface_area:block.mass.surface_area,
            centroid:block.mass.centroid,
        };
        let model=builder.finish()?;
        let radial=radial_certificate(&model,eps)
            .map_err(|e|invalid(format!("G8c radial audit: {e:?}")))?;
        if model.shells.len()!=1 ||model.vertices.len()!=12||
            model.edges.len()!=18||model.coedges.len()!=36||
            model.faces.len()!=8||model.loops.len()!=8||
            model.euler_characteristic()!=2||model.summary().genus!=0||
            radial.shell_components!=1||radial.total_genus!=0||
            radial.euler_characteristic!=2{
            return Ok(TiltedOneWallShellResult::Held(Hold::UnexpectedShellStructure));
        }
        Ok(TiltedOneWallShellResult::Ready(Self{
            model,faces:(0..8).map(FaceId).collect::<Vec<_>>().try_into()
                .map_err(|_|invalid("G8c missing eight faces"))?,
            edges:result_edges,
            source_faces:sources.try_into()
                .map_err(|_|invalid("G8c missing source face bindings"))?,
            physical_cuts:bindings.try_into()
                .map_err(|_|invalid("G8c missing physical cut bindings"))?,
            radial,source_volume:volume.clone(),
            closed_shell_proven:true,boundary_flux_certified:false,
            can_apply:false,
        }))
    }
}
