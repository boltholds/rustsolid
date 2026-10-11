//! A real, bounded genus-zero analytic Boolean: a vertical blind bore
//! entering the *top* plane of a certified axis-aligned box.
//!
//! The top face receives a circular inner loop; the bore side retains a
//! periodic cylindrical support and seam; a *new* planar circular floor
//! face closes the cavity. No source triangle ever defines these carriers.
//! Only cutters strictly inside the side walls and ending strictly above
//! the block bottom are admitted. Their top cap must exactly meet the box
//! top. Region classification and mass are analytic, not inferred from
//! the viewport triangulation.
use crate::{
    BrepBody, BrepFace, BrepModel, BrepOrigin, BrepSource, Circle2, Curve2,
    Curve3, EdgeId, FaceGeometry, FaceId, Frame3, GeometryError, LoopRole,
    MassProperties, Mesh, ParameterRange, PlaneSurface, Point2Param, Point3,
    Solid, Surface3, SurfaceId, ThroughHoleBrep, GeometryTolerance,
    GraphState, TrimContact, IntersectionGraph, IntersectionGraphBuilder,
};
use std::collections::BTreeMap;
use std::f64::consts::{PI,TAU};

fn invalid(s:impl Into<String>)->GeometryError{GeometryError::InvalidEdit(s.into())}
fn length(p:Point3)->f64{p.x.hypot(p.y).hypot(p.z)}
fn sign_area(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn signed_polygon_area(points:&[[f64;2]])->f64 {
    let first=points[0];
    (0..points.len()).map(|i|sign_area(first,points[i],
        points[(i+1)%points.len()])).sum::<f64>()*0.5
}
fn perimeter(points:&[[f64;2]])->f64 {
    (0..points.len()).map(|i|{
        let a=points[i];let b=points[(i+1)%points.len()];
        (a[0]-b[0]).hypot(a[1]-b[1])
    }).sum()
}
fn probe_inside(point:[f64;2],ring:&[[f64;2]])->bool{
    let mut inside=false;
    for i in 0..ring.len() {
        let a=ring[i];let b=ring[(i+1)%ring.len()];
        if (a[1]>point[1])!=(b[1]>point[1]) &&
            point[0]<a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]) {
            inside=!inside;
        }
    }
    inside
}

/// Source-bound analytic opening trace, not an unverified circle obtained
/// from a cutter bounding box. The upper block face (left FaceId 1) and the
/// real periodic cutter side (right FaceId 2) must intersect in exactly one
/// closed, transverse, full-circle trace whose samples lie on the cap plane
/// and the cylinder's mathematical radius. Near tangencies, fragmented
/// branches and coincident support ambiguities are refused.
fn certify_opening_trace(block:&Solid,cutter:&crate::CylindricalBrep)
    ->Result<IntersectionGraph,GeometryError>
{
    let first=BrepModel::from_polyhedral(block,crate::ModelUnits::default())?;
    let second=BrepModel::from_cylinder(cutter)?;
    let tolerance=GeometryTolerance {
        absolute_length:block.tolerance.absolute_length
            .max(cutter.tolerance.absolute_length),
        relative_length:block.tolerance.relative_length
            .max(cutter.tolerance.relative_length),
        angular:block.tolerance.angular.max(cutter.tolerance.angular),
    };
    let graph=IntersectionGraphBuilder::new(tolerance)?
        .intersect_faces(&first,FaceId(1),&second,FaceId(2))?;
    graph.validate_against(&first,FaceId(1),&second,FaceId(2))?;
    if graph.state!=GraphState::Ready||graph.arcs.len()!=1||
       !graph.arcs[0].closed ||
       graph.arcs[0].contact==TrimContact::Tangent {
        return Err(invalid("blind opening has no unique certified closed transverse circle"));
    }
    let b=block.bbox;
    let scale=length(b.max.sub(b.min)).max(cutter.radius*2.0);
    let floating=16.0*f64::EPSILON*
        [b.max.y,cutter.frame.origin.y].into_iter()
        .map(f64::abs).fold(1.0_f64,f64::max);
    let eps=tolerance.length_at(scale).max(floating);
    for fraction in [0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0] {
        let sample=graph.arcs[0].evaluate(fraction)?;
        if !sample.is_finite() ||
            (sample.y-b.max.y).abs()>eps ||
            ((sample.x-cutter.frame.origin.x).hypot(
                sample.z-cutter.frame.origin.z)-cutter.radius).abs()>eps {
            return Err(invalid("blind opening intersection trace leaves physical cap or cylinder"));
        }
    }
    Ok(graph)
}

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct BlindHoleParameters {
    pub center_x:f64,
    pub center_z:f64,
    pub radius:f64,
    /// Positive depth measured downward from the box upper (+Y) face.
    pub depth:f64,
}

#[derive(Debug,Clone)]
pub struct BlindHoleBrep {
    template:Solid,
    parameters:BlindHoleParameters,
    model:BrepModel,
    revision:u64,
    undo:Vec<BlindHoleParameters>,
    redo:Vec<BlindHoleParameters>,
}

/// Signed facet integrals are independent of the true analytic B-rep mass.
#[derive(Debug,Clone)]
pub struct BlindHoleDisplay {
    pub mesh:Mesh,
    pub faceted_volume:f64,
    pub faceted_surface_area:f64,
    pub analytic_volume:f64,
    pub faceting_volume_error:f64,
    pub opening_face:FaceId,
    pub floor_face:FaceId,
    pub side_face:FaceId,
    pub opening_triangles:usize,
    pub circle_segments:usize,
}

impl BlindHoleBrep {
    pub fn from_box(block:&Solid,parameters:BlindHoleParameters)
        ->Result<Self,GeometryError>
    {
        crate::boolean::certify_box(block).map_err(|e|invalid(e.to_string()))?;
        let model=Self::build(block,parameters)?;
        Ok(Self{template:block.clone(),parameters,model,revision:0,
            undo:Vec::new(),redo:Vec::new()})
    }
    /// A positive-Y cylinder is accepted iff its upper planar cap is on
    /// the box upper face, its lower cap is strictly inside the material,
    /// and the full circle remains strictly inside the four side planes.
    pub(crate) fn from_operands(block:&Solid,cutter:&crate::CylindricalBrep)
        ->Result<Self,GeometryError>
    {
        crate::boolean::certify_box(block).map_err(|e|invalid(e.to_string()))?;
        cutter.validate()?;
        let axis=cutter.frame.normal;
        let angular=block.tolerance.angular.max(cutter.tolerance.angular);
        if axis.x.abs()>angular || (axis.y-1.0).abs()>angular ||
            axis.z.abs()>angular {
            return Err(invalid("blind bore requires an upright +Y cylinder"));
        }
        let b=block.bbox;
        let extent=length(b.max.sub(b.min)).max(cutter.radius*2.0);
        let eps=block.tolerance.length_at(extent)
            .max(cutter.tolerance.length_at(extent));
        let upper=cutter.frame.origin.y+cutter.height;
        let representable=[upper,b.max.y,cutter.frame.origin.y,cutter.height]
            .into_iter().map(f64::abs).fold(1.0_f64,f64::max);
        let plane_budget=16.0*f64::EPSILON*representable;
        if !upper.is_finite() ||
            (upper-b.max.y).abs()>plane_budget ||
            cutter.frame.origin.y<=b.min.y+eps ||
            cutter.frame.origin.y>=b.max.y-eps {
            return Err(invalid("blind bore requires upper cap on box top and lower cap strictly inside"));
        }
        certify_opening_trace(block,cutter)?;
        let p=BlindHoleParameters{
            center_x:cutter.frame.origin.x,center_z:cutter.frame.origin.z,
            radius:cutter.radius,depth:cutter.height,
        };
        Self::from_box(block,p)
    }
    pub fn parameters(&self)->BlindHoleParameters{self.parameters}
    /// Reconstruct and validate the exact analytic source-face intersection
    /// defining the closed opening. The graph is read-only and source-bound.
    pub fn opening_trace(&self)->Result<IntersectionGraph,GeometryError>{
        let b=self.template.bbox;
        let p=self.parameters;
        let cutter=crate::CylindricalBrep::upright(
            Point3{x:p.center_x,y:b.max.y-p.depth,z:p.center_z},
            p.radius,p.depth,self.template.tolerance)?;
        certify_opening_trace(&self.template,&cutter)
    }
    pub fn template_solid(&self)->&Solid{&self.template}
    pub fn model(&self)->&BrepModel{&self.model}
    pub fn revision(&self)->u64{self.revision}
    pub fn undo_depth(&self)->usize{self.undo.len()}
    pub fn redo_depth(&self)->usize{self.redo.len()}
    pub fn validate(&self)->Result<(),GeometryError>{
        let expected=Self::build(&self.template,self.parameters)?;
        if self.model.origin!=expected.origin ||
            self.model.summary()!=expected.summary() ||
            self.model.vertices.iter().map(|v|v.position).collect::<Vec<_>>()!=
                expected.vertices.iter().map(|v|v.position).collect::<Vec<_>>() ||
            self.model.geometry.curves3!=expected.geometry.curves3 ||
            self.model.geometry.curves2!=expected.geometry.curves2 ||
            self.model.geometry.surfaces!=expected.geometry.surfaces ||
            self.model.mass!=expected.mass || self.model.bbox!=expected.bbox {
            return Err(invalid("blind hole no longer matches its independently reconstructed analytic source"));
        }
        self.model.validate()
    }
    pub fn begin_edit(&mut self)->BlindHoleTransaction<'_>{
        let staged=self.parameters;
        BlindHoleTransaction{owner:self,staged,failed:false}
    }
    pub fn undo(&mut self)->Result<bool,GeometryError>{
        let Some(parameters)=self.undo.last().copied()else{return Ok(false)};
        let candidate=Self::build(&self.template,parameters)?;
        let next=self.revision.checked_add(1)
            .ok_or_else(||invalid("blind-hole revision overflow"))?;
        self.undo.pop();
        self.redo.push(self.parameters);
        self.parameters=parameters;self.model=candidate;self.revision=next;
        Ok(true)
    }
    pub fn redo(&mut self)->Result<bool,GeometryError>{
        let Some(parameters)=self.redo.last().copied()else{return Ok(false)};
        let candidate=Self::build(&self.template,parameters)?;
        let next=self.revision.checked_add(1)
            .ok_or_else(||invalid("blind-hole revision overflow"))?;
        self.redo.pop();
        self.undo.push(self.parameters);
        self.parameters=parameters;self.model=candidate;self.revision=next;
        Ok(true)
    }

    fn build(block:&Solid,p:BlindHoleParameters)->Result<BrepModel,GeometryError>{
        let b=block.bbox;
        let height=b.max.y-b.min.y;
        let local=length(b.max.sub(b.min));
        let eps=block.tolerance.length_at(local);
        if !p.center_x.is_finite()||!p.center_z.is_finite()||
            !p.radius.is_finite()||!p.depth.is_finite() ||
            p.depth<=eps||p.depth>=height-eps {
            return Err(invalid("blind bore depth must lie strictly within the box height"));
        }
        let through=ThroughHoleBrep::from_box(block,
            crate::HoleParameters{center_x:p.center_x,center_z:p.center_z,
                radius:p.radius})?;
        let mut m=through.model().clone();
        let floor_y=b.max.y-p.depth;
        let delta=floor_y-b.min.y;
        let floor_origin=Point3{x:p.center_x,y:floor_y,z:p.center_z};
        let top_origin=Point3{y:b.max.y,..floor_origin};
        // Rebuild the *lower* circular cut edge and the seam from their
        // existing exact carriers. The top circular edge and all box
        // geometry are left unchanged.
        m.vertices[8].position.y=floor_y;
        let lower_id=m.geometry.edges[12].curve.0 as usize;
        let seam_id=m.geometry.edges[14].curve.0 as usize;
        let Curve3::Circle(mut bottom_circle)=m.geometry.curves3[lower_id]else{
            return Err(invalid("canonical lower edge is no longer a true circle"));
        };
        bottom_circle.frame.origin.y=floor_y;
        m.geometry.curves3[lower_id]=Curve3::Circle(bottom_circle);
        m.geometry.curves3[seam_id]=Curve3::Line(crate::Line3{
            origin:m.vertices[8].position,direction:bottom_circle.frame.normal,
        });
        m.geometry.edges[14].domain=ParameterRange::Bounded{
            start:0.0,end:p.depth};
        // The cylindrical side is physically bounded by the new planar
        // floor and the original top face. Its two seam pcurves are kept
        // distinct in the periodic chart.
        m.faces[6].geometry.bounds.v=ParameterRange::Bounded{
            start:0.0,end:p.depth};
        let side_id=m.faces[6].geometry.surface.0 as usize;
        let Surface3::Cylinder(mut side_carrier)=m.geometry.surfaces[side_id]else {
            return Err(invalid("canonical bore side lost its cylindrical support"));
        };
        // Reanchor cylindrical v=0 to the new floor plane. Merely shortening
        // the v interval without moving the support would leave the lower
        // circle geometrically disconnected from its 3D/UV side coedges.
        side_carrier.frame.origin.y=floor_y;
        m.geometry.surfaces[side_id]=Surface3::Cylinder(side_carrier);
        for (coedge_id,from,to) in [
            (27usize,Point2Param{u:TAU,v:p.depth},Point2Param{u:0.0,v:p.depth}),
            (28usize,Point2Param{u:TAU,v:0.0},Point2Param{u:TAU,v:p.depth}),
            (29usize,Point2Param{u:0.0,v:p.depth},Point2Param{u:0.0,v:0.0}),
        ]{
            let idx=m.geometry.coedges[coedge_id].pcurve.0 as usize;
            m.geometry.curves2[idx]=Curve2::Line(crate::Line2{
                origin:from,delta:Point2Param{u:to.u-from.u,v:to.v-from.v}
            });
        }
        // In the current B-rep chart convention a floor circular pcurve
        // borrows the box bottom-plane (-Y) parameterization. Its distinct
        // physical *shell* orientation (+Y) is preserved by the certified
        // display mesher; the common graph does not yet have an explicit
        // face-sense bit separate from support-surface parameterization.
        let lower_plane=match m.face_surface(FaceId(0)){
            Some(Surface3::Plane(plane))=>plane,
            _=>return Err(invalid("source box bottom plane missing")),
        };
        let floor_frame=Frame3 {
            origin:Point3{y:lower_plane.frame.origin.y+delta,..lower_plane.frame.origin},
            ..lower_plane.frame
        };
        let floor_surface=SurfaceId(u32::try_from(m.geometry.surfaces.len())
            .map_err(|_|invalid("floor surface catalog exhausted"))?);
        m.geometry.surfaces.push(Surface3::Plane(PlaneSurface{frame:floor_frame}));
        // Transplant the now-internal lower inner loop onto a new disk face.
        // It has exactly one coedge whose actual 3D circle is shared with
        // the cylinder side's lower boundary.
        m.faces[0].loops.retain(|&lp|lp.0!=6);
        m.loops[6].face=FaceId(7);
        m.loops[6].role=LoopRole::Outer;
        m.coedges[24].face=FaceId(7);
        m.geometry.coedges[24].face=FaceId(7);
        let mut floor_geometry=m.faces[0].geometry.clone();
        floor_geometry.surface=floor_surface;
        floor_geometry.outer_loop=crate::LoopId(6);
        // Circle2 only stores a center and radius; translation of the
        // floor chart along its normal preserves the original UV center.
        let floor_uv=m.geometry.coedges[24].pcurve.0 as usize;
        let Curve2::Circle(circle_uv)=m.geometry.curves2[floor_uv]else{
            return Err(invalid("floor UV pcurve must remain circular"));
        };
        m.geometry.curves2[floor_uv]=Curve2::Circle(Circle2::new(
            floor_frame.project(floor_origin)?,circle_uv.radius,block.tolerance)?);
        m.faces.push(BrepFace{id:FaceId(7),loops:vec![crate::LoopId(6)],
            geometry:floor_geometry});
        m.shells[0].faces.push(FaceId(7));
        let removed=PI*p.radius*p.radius*p.depth;
        let volume=block.mass.volume-removed;
        if !volume.is_finite()||volume<=0.0 {
            return Err(invalid("blind-hole volume must remain positive"));
        }
        let original=block.mass.centroid;
        let bore_center=Point3{x:p.center_x,
            y:(floor_origin.y+top_origin.y)/2.0,z:p.center_z};
        let shift=original.sub(bore_center).scale(removed/volume);
        m.mass=MassProperties{
            volume,
            surface_area:block.mass.surface_area+TAU*p.radius*p.depth,
            centroid:original.add(shift),
        };
        m.origin=BrepOrigin::AnalyticBlindHole;
        m.validate()?;
        if m.summary().genus!=0 || m.euler_characteristic()!=2 ||
            m.faces[0].loops.len()!=1 || m.faces[1].loops.len()!=2 ||
            m.loops.iter().filter(|lp|lp.role==LoopRole::Inner).count()!=1 {
            return Err(invalid("blind bore should produce one closed genus-zero shell with one top inner loop"));
        }
        Ok(m)
    }

    /// A purely derived display tessellation; triangle count does not
    /// influence exact analytic carrier geometry or computed mass.
    pub fn tessellate(&self,segments:usize)->Result<BlindHoleDisplay,GeometryError>{
        if !(12..=512).contains(&segments) {
            return Err(invalid("blind-hole display segments must be 12..=512"));
        }
        self.validate()?;
        let p=self.parameters;
        let b=self.template.bbox;
        let bottom=b.min.y;let top=b.max.y;let floor=top-p.depth;
        let scale=length(b.max.sub(b.min));
        let coord=[b.min,b.max].into_iter().flat_map(|v|
            [v.x.abs(),v.y.abs(),v.z.abs()]).fold(1.0f64,f64::max);
        let eps=self.template.tolerance.length_at(scale)
            .max(32.0*f64::EPSILON*coord);
        let mut mesh=Mesh {
            vertices:self.template.mesh.vertices.clone(),
            triangles:Vec::new(),triangle_faces:Vec::new(),
        };
        // Retain the source box vertices in their original canonical order,
        // not the per-triangle display mesh's arbitrary arrangement.
        mesh.vertices=self.template.vertices.iter().map(|v|v.position).collect();
        if mesh.vertices.len()!=8 {
            return Err(invalid("blind-hole display requires eight box vertices"));
        }
        let mut outer=Vec::new();let mut upper_ids=Vec::new();
        for &id in &self.template.faces[0].boundary {
            let index=id.0 as usize;
            if index>=4 {return Err(invalid("noncanonical bottom-cap vertex"));}
            let v=mesh.vertices[index];
            outer.push([v.x-b.min.x,v.z-b.min.z]);
            upper_ids.push(u32::try_from(index+4).unwrap());
        }
        if outer.len()!=4 || signed_polygon_area(&outer)<=0.0 {
            return Err(invalid("box outer boundary has incorrect XZ winding"));
        }
        let top_circle=self.model.edge_curve(EdgeId(13))
            .ok_or_else(||invalid("missing upper 3D circular carrier"))?;
        let floor_circle=self.model.edge_curve(EdgeId(12))
            .ok_or_else(||invalid("missing floor 3D circular carrier"))?;
        let mut ring=Vec::with_capacity(segments);
        let mut floor_indices=Vec::with_capacity(segments);
        let mut top_indices=Vec::with_capacity(segments);
        let mut coords=outer.iter().flat_map(|v|*v).collect::<Vec<f64>>();
        for i in 0..segments {
            let t=TAU*i as f64/segments as f64;
            let down=floor_circle.evaluate(t)?;
            let up=top_circle.evaluate(t)?;
            if (down.x-up.x).hypot(down.z-up.z)>eps ||
                (down.y-floor).abs()>eps||(up.y-top).abs()>eps{
                return Err(invalid("circular blind bore top/floor carriers disagree"));
            }
            let index=u32::try_from(mesh.vertices.len())
                .map_err(|_|invalid("blind-hole display vertex budget exceeded"))?;
            mesh.vertices.extend([down,up]);
            floor_indices.push(index);top_indices.push(index+1);
            let point=[up.x-b.min.x,up.z-b.min.z];
            ring.push(point);
            coords.extend_from_slice(&point);
        }
        if signed_polygon_area(&ring)>=0.0 {
            return Err(invalid("upper circular trim does not follow clockwise XZ contour"));
        }
        let triangle_ids=earcutr::earcut(&coords,&[4usize],2)
            .map_err(|e|invalid(format!("blind cap triangulation failed: {e:?}")))?;
        if triangle_ids.is_empty()||triangle_ids.len()%3!=0{
            return Err(invalid("blind cap triangulation has invalid length"));
        }
        let mut cap_area=0.0;
        for triangle in triangle_ids.chunks_exact(3){
            let mut a=triangle[0];let mut c=triangle[1];let mut d=triangle[2];
            if [a,c,d].iter().any(|&id|id>=4+segments)||
                a==c||c==d||a==d {
                return Err(invalid("blind opening triangle indexes invalid contour"));
            }
            let get=|i:usize|if i<4{outer[i]}else{ring[i-4]};
            let mut signed=sign_area(get(a),get(c),get(d));
            if !signed.is_finite() ||
                signed.abs()<=self.template.tolerance.area_at(scale)*0.1{
                return Err(invalid("blind upper cap contains degenerate triangle"));
            }
            if signed<0.0{std::mem::swap(&mut c,&mut d);signed=-signed;}
            let center=[(get(a)[0]+get(c)[0]+get(d)[0])/3.0,
                (get(a)[1]+get(c)[1]+get(d)[1])/3.0];
            if !probe_inside(center,&outer)||probe_inside(center,&ring){
                return Err(invalid("opening cap triangle crosses excluded circular domain"));
            }
            cap_area+=signed*0.5;
            let to_id=|i:usize|if i<4{upper_ids[i]}else{top_indices[i-4]};
            mesh.triangles.push([to_id(a),to_id(d),to_id(c)]);
            mesh.triangle_faces.push(FaceId(1));
        }
        let ring_area=-signed_polygon_area(&ring);
        let area_expected=signed_polygon_area(&outer)-ring_area;
        let budget=self.template.tolerance.area_at(scale)*(segments+8) as f64+
            eps*scale*32.0;
        if !cap_area.is_finite()||cap_area<=0.0 ||
            (cap_area-area_expected).abs()>budget {
            return Err(invalid("blind upper opening triangulated region area mismatch"));
        }
        // The box's original bottom and four vertical side faces are
        // unchanged. Each retains its actual face ownership.
        for (triangle,face) in self.template.mesh.triangles.iter()
            .zip(&self.template.mesh.triangle_faces){
            if face.0!=1 {
                mesh.triangles.push(*triangle);
                mesh.triangle_faces.push(*face);
            }
        }
        let center_idx=u32::try_from(mesh.vertices.len())
            .map_err(|_|invalid("floor disk center vertex budget exceeded"))?;
        mesh.vertices.push(Point3{x:p.center_x,y:floor,z:p.center_z});
        for j in 0..segments {
            let k=(j+1)%segments;
            let f0=floor_indices[j];let f1=floor_indices[k];
            let t0=top_indices[j];let t1=top_indices[k];
            mesh.triangles.extend([
                [center_idx,f0,f1], // positive Y material-to-void floor
                [f0,t0,t1],[f0,t1,f1], // inward cylinder side
            ]);
            mesh.triangle_faces.extend([FaceId(7),FaceId(6),FaceId(6)]);
        }
        let (faceted_volume,faceted_surface_area)=
            certify_blind_display(&self.template,&self.model,&mesh,p,eps)?;
        let area_floor=ring_area;
        let expected_faceted_volume=self.template.mass.volume-area_floor*p.depth;
        let expected_faceted_surface=self.template.mass.surface_area+
            perimeter(&ring)*p.depth;
        let vol_eps=self.template.mass.volume*1.0e-8+eps*scale*scale*16.0;
        let area_eps=self.template.mass.surface_area*1.0e-8+eps*scale*64.0;
        if (faceted_volume-expected_faceted_volume).abs()>vol_eps ||
           (faceted_surface_area-expected_faceted_surface).abs()>area_eps ||
           faceted_volume+vol_eps<self.model.mass.volume {
            return Err(invalid("blind-bore display area or volume disagrees with polygonal analytic oracle"));
        }
        Ok(BlindHoleDisplay{mesh,faceted_volume,faceted_surface_area,
            analytic_volume:self.model.mass.volume,
            faceting_volume_error:faceted_volume-self.model.mass.volume,
            opening_face:FaceId(1),floor_face:FaceId(7),side_face:FaceId(6),
            opening_triangles:triangle_ids.len()/3,circle_segments:segments})
    }
}
fn certify_blind_display(box_:&Solid,model:&BrepModel,mesh:&Mesh,
    p:BlindHoleParameters,eps:f64)->Result<(f64,f64),GeometryError>{
    if mesh.triangles.len()!=mesh.triangle_faces.len()||
        mesh.vertices.iter().any(|v|!v.is_finite()) {
        return Err(invalid("blind-hole display missing per-triangle face provenance"));
    }
    let anchor=box_.bbox.min;
    let mut edge_usage=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut six_volume=0.0;
    let mut area=0.0;
    for (tri,owner) in mesh.triangles.iter().zip(&mesh.triangle_faces){
        let [a,b,c]=*tri;
        if a==b||b==c||a==c ||
            [a,b,c].iter().any(|&id|id as usize>=mesh.vertices.len()) ||
            owner.0 as usize>=model.faces.len(){
            return Err(invalid("blind display triangle has invalid vertices or parent face"));
        }
        let [v0,v1,v2]=tri.map(|i|mesh.vertices[i as usize]);
        let normal=v1.sub(v0).cross(v2.sub(v0));
        let norm=length(normal);
        if !norm.is_finite() || norm<=0.0 {
            return Err(invalid("blind display contains zero-area facet"));
        }
        let centroid=v0.add(v1).add(v2).scale(1.0/3.0);
        match *owner {
            FaceId(6)=>{
                if normal.dot(Point3{x:p.center_x-centroid.x,y:0.0,
                    z:p.center_z-centroid.z})<=0.0 ||
                   [v0,v1,v2].iter().any(|v|{
                       (v.x-p.center_x).hypot(v.z-p.center_z).sub_abs(p.radius)>eps
                   }){
                    return Err(invalid("blind cylindrical facet has incorrect inward orientation or radius"));
                }
            },
            FaceId(7)=>{
                if normal.y<=0.0 ||
                    [v0,v1,v2].iter().any(|v|
                        (v.y-(box_.bbox.max.y-p.depth)).abs()>eps){
                    return Err(invalid("blind circular floor is not oriented upward"));
                }
            },
            _=>{
                let surface=model.face_surface(*owner)
                    .ok_or_else(||invalid("display references unsupported planar source"))?;
                let Surface3::Plane(plane)=surface else{
                    return Err(invalid("non-planar outer box triangle"));
                };
                if normal.dot(plane.frame.normal)<=0.0 ||
                   [v0,v1,v2].iter().any(|v|
                        (v.sub(plane.frame.origin).dot(plane.frame.normal)).abs()>eps){
                    return Err(invalid("box face triangle lost source plane or outward winding"));
                }
            }
        }
        for (from,to) in [(a,b),(b,c),(c,a)]{
            edge_usage.entry((from.min(to),from.max(to))).or_default()
                .push(from<to);
        }
        six_volume+=v0.sub(anchor).dot(v1.sub(anchor).cross(v2.sub(anchor)));
        area+=norm*0.5;
    }
    if edge_usage.values().any(|uses|uses.len()!=2 || uses[0]==uses[1]){
        return Err(invalid("blind display mesh fails closed two-use manifold edge certificate"));
    }
    let volume=six_volume/6.0;
    if !volume.is_finite()||volume<=0.0||!area.is_finite()||area<=0.0{
        return Err(invalid("blind display signed volume or area invalid"));
    }
    Ok((volume,area))
}
// Kept local to this constrained analytic floor classifier.
trait AbsoluteDifference {fn sub_abs(self,rhs:f64)->f64;}
impl AbsoluteDifference for f64 {
    fn sub_abs(self,rhs:f64)->f64{(self-rhs).abs()}
}

impl BrepSource for BlindHoleBrep {
    fn brep_model(&self)->Result<BrepModel,GeometryError>{
        self.validate()?;
        Ok(self.model.clone())
    }
    fn display_mesh(&self,segments:usize)->Result<Mesh,GeometryError>{
        Ok(self.tessellate(segments)?.mesh)
    }
}
pub struct BlindHoleTransaction<'a>{
    owner:&'a mut BlindHoleBrep,
    staged:BlindHoleParameters,
    failed:bool,
}
impl BlindHoleTransaction<'_> {
    fn set(&mut self,change:impl FnOnce(&mut BlindHoleParameters))
        ->Result<(),GeometryError>{
        if self.failed{return Err(invalid("blind-hole transaction is poisoned"));}
        change(&mut self.staged);
        let p=self.staged;
        if !p.center_x.is_finite()||!p.center_z.is_finite()||
           !p.radius.is_finite()||p.radius<=0.0||
           !p.depth.is_finite()||p.depth<=0.0 {
            self.failed=true;
            return Err(invalid("invalid staged blind-hole parameters"));
        }
        Ok(())
    }
    pub fn set_radius(&mut self,r:f64)->Result<(),GeometryError>{
        self.set(|p|p.radius=r)
    }
    pub fn set_center(&mut self,x:f64,z:f64)->Result<(),GeometryError>{
        self.set(|p|{p.center_x=x;p.center_z=z;})
    }
    pub fn set_depth(&mut self,depth:f64)->Result<(),GeometryError>{
        self.set(|p|p.depth=depth)
    }
    pub fn commit(self)->Result<u64,GeometryError>{
        if self.failed{return Err(invalid("cannot commit poisoned blind-hole transaction"));}
        let candidate=BlindHoleBrep::build(&self.owner.template,self.staged)?;
        let next=self.owner.revision.checked_add(1)
            .ok_or_else(||invalid("blind-hole revision overflow"))?;
        if self.owner.undo.len()>=128{self.owner.undo.remove(0);}
        self.owner.undo.push(self.owner.parameters);
        self.owner.redo.clear();
        self.owner.model=candidate;
        self.owner.parameters=self.staged;
        self.owner.revision=next;
        Ok(next)
    }
}
