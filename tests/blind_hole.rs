//! Source-controlled mathematical and topological witnesses for a real
//! blind circular bore. The result contains a planar floor that is NOT a
//! mesh approximation of a hole or an extra phantom bottom-cap opening.
use rustsolid::{
    BlindHoleBrep,BlindHoleParameters,BooleanKernel,BooleanOperation,
    BooleanOutcome,BooleanOperand,BooleanRegionClassifier,VolumePosition,
    BrepArena,BrepBody,CylindricalBrep,FaceId,EdgeId,GeometryTolerance,
    LoopRole,Point3,Solid,Surface3,Curve3,HoleParameters,ThroughHoleBrep,
    TopologyEntity,
};
use std::f64::consts::PI;
use std::collections::BTreeMap;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn par(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn blind()->BlindHoleBrep{BlindHoleBrep::from_box(&box_(),par(0.,0.,1.,2.)).unwrap()}
fn close(a:f64,b:f64){assert!((a-b).abs()<1e-9,"{a} != {b} (err {})",(a-b).abs());}
fn assert_brep(b:&BlindHoleBrep){
    b.validate().unwrap();
    let m=b.model();m.validate().unwrap();
    let s=m.summary();
    assert_eq!(s.vertices,10);
    assert_eq!(s.edges,15);
    assert_eq!(s.coedges,30);
    assert_eq!(s.faces,8);
    assert_eq!(s.loops,9);
    assert_eq!(s.shells,1);
    assert_eq!(s.euler_characteristic,2);
    assert_eq!(s.genus,0);
    assert_eq!(s.circular_edges,2);
    assert_eq!(s.cylindrical_faces,1);
    assert_eq!(s.plane_faces,7);
    assert_eq!(s.seam_edges,1);
    assert_eq!(s.closed_edges,2);
    assert_eq!(m.faces[0].loops.len(),1);
    assert_eq!(m.faces[1].loops.len(),2);
    assert_eq!(m.faces[7].loops.len(),1);
    assert_eq!(m.loops.iter().filter(|lp|lp.role==LoopRole::Inner).count(),1);
    assert_eq!(m.loops[6].role,LoopRole::Outer);
    assert_eq!(m.loops[6].face,FaceId(7));
    assert_eq!(m.coedges[24].face,FaceId(7));
    assert_eq!(m.coedges[26].face,FaceId(6));
    assert_eq!(m.coedges[26].twin,rustsolid::CoedgeId(24));
    assert!(matches!(m.face_surface(FaceId(7)),Some(Surface3::Plane(_))));
    assert!(matches!(m.face_surface(FaceId(6)),Some(Surface3::Cylinder(_))));
    assert!(matches!(m.edge_curve(EdgeId(12)),Some(Curve3::Circle(_))));
    assert!(matches!(m.edge_curve(EdgeId(13)),Some(Curve3::Circle(_))));
}
fn prove_display(b:&BlindHoleBrep,n:usize)->rustsolid::BlindHoleDisplay{
    let m=b.tessellate(n).unwrap();
    let mut use_edges=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut counts=BTreeMap::<u32,usize>::new();
    let mut signed_six=0.0;
    let anchor=box_().bbox.min;
    for (tri,face) in m.mesh.triangles.iter().zip(&m.mesh.triangle_faces){
        *counts.entry(face.0).or_insert(0)+=1;
        let [a,b,c]=*tri;
        for (u,v) in [(a,b),(b,c),(c,a)] {
            use_edges.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
        let v0=m.mesh.vertices[a as usize];
        let v1=m.mesh.vertices[b as usize];
        let v2=m.mesh.vertices[c as usize];
        let u=p(v0.x-anchor.x,v0.y-anchor.y,v0.z-anchor.z);
        let v=p(v1.x-anchor.x,v1.y-anchor.y,v1.z-anchor.z);
        let w=p(v2.x-anchor.x,v2.y-anchor.y,v2.z-anchor.z);
        signed_six+=u.x*(v.y*w.z-v.z*w.y)
            +u.y*(v.z*w.x-v.x*w.z)
            +u.z*(v.x*w.y-v.y*w.x);
    }
    assert!(use_edges.values().all(|v|v.len()==2 && v[0]!=v[1]));
    assert_eq!(counts[&0],2);
    assert_eq!(counts[&1],n+4);
    for face in 2..6 {assert_eq!(counts[&face],2);}
    assert_eq!(counts[&6],2*n);
    assert_eq!(counts[&7],n);
    assert_eq!(m.mesh.vertices.len(),9+2*n);
    assert_eq!(m.mesh.triangles.len(),4*n+14);
    assert_eq!(m.mesh.triangles.len(),m.mesh.triangle_faces.len());
    close(signed_six/6.0,m.faceted_volume);
    assert!(m.faceting_volume_error>0.0);
    m
}
#[test]
fn blind_bore_has_new_planar_disk_and_genus_zero_not_a_second_through_cap(){
    let b=blind();assert_brep(&b);
    close(b.model().mass.volume,144.-2.*PI);
    close(b.model().mass.surface_area,168.+4.*PI);
    close(b.model().mass.centroid.x,0.);
    close(b.model().mass.centroid.y,(144.-4.*PI)/(144.-2.*PI));
    close(b.model().mass.centroid.z,0.);
}
#[test]
fn top_entry_boolean_builds_real_blind_hole_with_floor_provenance(){
    let block=BrepBody::from(box_());
    let cutter=BrepBody::from(CylindricalBrep::upright(
        p(0.,1.,0.),1.,2.,GeometryTolerance::default()).unwrap());
    let result=BooleanKernel::apply(&block,&cutter,BooleanOperation::Difference).unwrap();
    assert_eq!(result.face_sources.len(),8);
    let lower=result.face_sources.iter().find(|x|x.result_face==FaceId(7)).unwrap();
    assert_eq!(lower.operand,BooleanOperand::Right);
    assert_eq!(lower.source_face,FaceId(0));
    assert!(lower.reversed);
    let side=result.face_sources.iter().find(|x|x.result_face==FaceId(6)).unwrap();
    assert_eq!(side.source_face,FaceId(2));
    assert!(side.reversed);
    match result.outcome {
        BooleanOutcome::BlindHole(body)=>assert_brep(&body),
        _=>panic!("expected independent exact blind-hole outcome"),
    }
    for op in [BooleanOperation::Union,BooleanOperation::Intersection]{
        assert!(BooleanKernel::apply(&block,&cutter,op).is_err());
    }
}
#[test]
fn old_through_hole_and_exact_planar_box_booleans_are_unchanged(){
    let cube=BrepBody::from(box_());
    let through=BrepBody::from(CylindricalBrep::upright(
        p(0.,-1.,0.),1.,4.,GeometryTolerance::default()).unwrap());
    let outcome=BooleanKernel::apply(&cube,&through,BooleanOperation::Difference).unwrap();
    assert!(matches!(outcome.outcome,BooleanOutcome::Analytic(_)));
    let box2=BrepBody::from(Solid::block(p(-2.,0.,-2.),2.,2.,2.).unwrap());
    let inter=BooleanKernel::apply(&cube,&box2,BooleanOperation::Intersection).unwrap();
    assert!(matches!(inter.outcome,BooleanOutcome::Solid(_)));
}
#[test]
fn invalid_partial_depth_and_contact_geometries_fail_closed(){
    for q in [
        par(0.,0.,1.,0.),
        par(0.,0.,1.,4.),
        par(0.,0.,1.,4.000001),
        par(0.,0.,1.,1e-12),
        par(2.5,0.,0.5,1.),
        par(0.,0.,3.,1.),
        par(f64::NAN,0.,1.,1.),
        par(0.,f64::INFINITY,1.,1.),
        par(0.,0.,-1.,1.),
    ]{assert!(BlindHoleBrep::from_box(&box_(),q).is_err(),"invalid {q:?}");}
    for (origin_y,height) in [(-1.,3.),(0.,2.),(1.,1.),
        (1.,3.),(-1.0+1e-11,4.0-1e-11)] {
        let cutter=CylindricalBrep::upright(p(0.,origin_y,0.),1.,height,
            GeometryTolerance::default()).unwrap();
        assert!(BooleanKernel::apply(&BrepBody::from(box_()),
            &BrepBody::from(cutter),BooleanOperation::Difference).is_err(),
            "mismatched cylinder cap or bottom clearance must not cut");
    }
}
#[test]
fn cylinder_void_and_real_floor_have_correct_three_way_membership(){
    let b=BrepBody::from(blind());
    let samples=[
        (p(0.,2.,0.),VolumePosition::Outside),
        (p(0.,3.,0.),VolumePosition::Outside),
        (p(0.,1.5,0.),VolumePosition::Outside),
        (p(0.,1.,0.),VolumePosition::Boundary),
        (p(0.,0.5,0.),VolumePosition::Inside),
        (p(1.,2.,0.),VolumePosition::Boundary),
        (p(2.,2.,0.),VolumePosition::Inside),
        (p(0.,-1.,0.),VolumePosition::Boundary),
        (p(3.,0.,0.),VolumePosition::Boundary),
        (p(4.,0.,0.),VolumePosition::Outside),
        (p(0.,3.001,0.),VolumePosition::Outside)
    ];
    for (point,expected) in samples {
        assert_eq!(BooleanRegionClassifier::classify_point(&b,point).unwrap(),
            expected,"incorrect material at {point:?}");
    }
}
#[test]
fn certified_blind_display_has_closed_opposite_oriented_triangle_pairs(){
    let b=blind();assert_brep(&b);
    let display=prove_display(&b,48);
    close(display.analytic_volume,144.-2.*PI);
    assert_eq!((display.opening_face,display.floor_face,display.side_face),
        (FaceId(1),FaceId(7),FaceId(6)));
    assert_eq!(display.opening_triangles,52);
}
#[test]
fn increasing_blind_circle_segments_converges_to_analytic_mass(){
    let b=blind();let exact=b.model().mass;
    let a=prove_display(&b,12);
    let c=prove_display(&b,48);
    let d=prove_display(&b,192);
    assert!(a.faceting_volume_error>c.faceting_volume_error);
    assert!(c.faceting_volume_error>d.faceting_volume_error);
    assert!(d.faceting_volume_error<c.faceting_volume_error*0.08);
    close(b.model().mass.volume,exact.volume);
    close(b.model().mass.surface_area,exact.surface_area);
}
#[test]
fn off_center_floor_removal_changes_centroid_at_correct_depth(){
    let p=par(0.75,-0.5,0.7,1.25);
    let b=BlindHoleBrep::from_box(&box_(),p).unwrap();assert_brep(&b);
    let removed=PI*0.7*0.7*1.25;
    let mass=144.-removed;
    close(b.model().mass.volume,mass);
    close(b.model().mass.centroid.x,-0.75*removed/mass);
    close(b.model().mass.centroid.z,0.5*removed/mass);
    let cut_y=3.-1.25/2.;
    close(b.model().mass.centroid.y,1.+(1.-cut_y)*removed/mass);
}
#[test]
fn radius_center_and_depth_edits_are_one_strong_atomic_history_step(){
    let mut b=blind();
    let before=b.model().mass.volume;
    let mut tx=b.begin_edit();
    tx.set_radius(0.8).unwrap();
    tx.set_depth(2.5).unwrap();
    tx.set_center(0.2,-0.2).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    assert_brep(&b);
    close(b.model().mass.volume,144.-PI*0.8*0.8*2.5);
    assert_eq!(b.undo_depth(),1);
    assert!(b.undo().unwrap());
    close(b.model().mass.volume,before);
    assert_eq!(b.revision(),2);
    assert!(b.redo().unwrap());
    assert_eq!(b.revision(),3);
    close(b.model().mass.volume,144.-PI*0.8*0.8*2.5);
    assert_brep(&b);
}
#[test]
fn aborted_stage_and_unsupported_depth_leave_revision_and_mass_untouched(){
    let mut b=blind();let original=b.model().mass;
    {
        let mut tx=b.begin_edit();
        tx.set_depth(1.5).unwrap();
        assert!(tx.set_radius(f64::NAN).is_err());
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=b.begin_edit();
        tx.set_depth(5.).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(b.revision(),0);
    assert_eq!(b.undo_depth(),0);
    assert_eq!(b.model().mass,original);
    b.validate().unwrap();
}
#[test]
fn arena_updates_depth_atomically_and_reuses_live_semantic_handles(){
    let mut arena=BrepArena::from_blind_hole(blind()).unwrap();
    let snapshot=arena.snapshot().unwrap();
    let floor=arena.handle(TopologyEntity::Face(FaceId(7))).unwrap();
    let cylinder=arena.handle(TopologyEntity::Face(FaceId(6))).unwrap();
    let old=arena.model().unwrap().mass.volume;
    assert_eq!(arena.update_blind_hole(par(0.3,0.,1.25,1.5)).unwrap(),1);
    arena.validate().unwrap();
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(7))),Some(floor));
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(6))),Some(cylinder));
    assert!(arena.check_snapshot(&snapshot).is_err());
    assert_eq!(arena.undo_depth(),1);
    assert!(arena.update_blind_hole(par(2.9,0.,0.5,1.5)).is_err());
    assert_eq!(arena.revision(),1);
    assert!(arena.undo().unwrap().is_some());
    close(arena.model().unwrap().mass.volume,old);
    assert!(arena.redo().unwrap().is_some());
    assert_eq!(arena.revision(),3);
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(7))),Some(floor));
}
#[test]
fn shape_stays_queryable_via_shared_analytic_brep_and_detached_carrier_editor(){
    let input=BrepBody::from(blind());
    let graph=input.shared().unwrap();graph.validate().unwrap();
    assert_eq!(graph.summary().genus,0);
    let mut carrier=rustsolid::BrepCarrierEditor::from_body(&input).unwrap();
    let edge=carrier.handle(rustsolid::CarrierOwner::Edge(EdgeId(12))).unwrap();
    let mut tx=carrier.begin_transaction();
    tx.reparameterize_edge(edge,1.,0.25).unwrap();
    tx.commit().unwrap();
    carrier.validate().unwrap();
    close(carrier.model().mass.volume,graph.mass.volume);
    let values=rustsolid::geometry_query(&graph,
        "query {body {kind topology {faceCount genus circularEdges} mass {volume}}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(values["body"]["kind"],"analytic_blind_hole");
    assert_eq!(values["body"]["topology"]["faceCount"],8);
}
#[test]
fn localized_world_translation_does_not_change_exact_mass(){
    let shift=1e8;
    let block=Solid::block(p(shift-3.,shift-1.,shift-3.),6.,4.,6.).unwrap();
    let b=BlindHoleBrep::from_box(&block,par(shift+0.4,shift-0.2,0.8,2.))
        .unwrap();
    b.validate().unwrap();
    close(b.model().mass.volume,144.-PI*0.8*0.8*2.);
    let display=b.tessellate(48).unwrap();
    assert!(display.faceting_volume_error>0.);
    assert!((display.faceted_volume-display.analytic_volume).abs()<0.02);
}
#[test]
fn direct_cylinder_with_partial_depth_reversed_axis_does_not_silently_flip(){
    let box_=BrepBody::from(box_());
    let frame=rustsolid::Frame3::from_axes(
        p(0.,1.,0.),p(1.,0.,0.),p(0.,0.,1.),GeometryTolerance::default()).unwrap();
    let cutter=CylindricalBrep::new(frame,1.,2.,GeometryTolerance::default()).unwrap();
    assert!(BooleanKernel::apply(&box_,&BrepBody::from(cutter),
        BooleanOperation::Difference).is_err());
}

#[test]
fn source_intersection_graph_retains_closed_circular_opening_trace(){
    let block=BrepBody::from(box_());
    let cutter=BrepBody::from(CylindricalBrep::upright(
        p(0.,1.,0.),1.,2.,GeometryTolerance::default()).unwrap());
    let b=block.shared().unwrap();
    let c=cutter.shared().unwrap();
    let graph=rustsolid::IntersectionGraphBuilder::new(GeometryTolerance::default())
        .unwrap().intersect_faces(&b,FaceId(1),&c,FaceId(2)).unwrap();
    assert_eq!(graph.kind(),"ready","expected circular source-bound cut trace");
    assert_eq!(graph.arcs.len(),1);
    assert!(graph.arcs[0].closed);
    graph.validate_against(&b,FaceId(1),&c,FaceId(2)).unwrap();
    for fraction in [0.0,0.125,0.3,0.5,0.875,1.0] {
        let q=graph.arcs[0].evaluate(fraction).unwrap();
        close(q.y,3.);
        close(q.x.hypot(q.z),1.);
    }
}

#[test]
fn constructed_blind_result_exposes_source_bound_analytic_opening_graph(){
    let shape=blind();
    let graph=shape.opening_trace().unwrap();
    assert_eq!(graph.state,rustsolid::GraphState::Ready);
    assert_eq!(graph.arcs.len(),1);
    assert!(graph.arcs[0].closed);
    assert_ne!(graph.arcs[0].contact,rustsolid::TrimContact::Tangent);
    for fraction in [0.0,0.125,0.5,0.875]{
        let point=graph.arcs[0].evaluate(fraction).unwrap();
        close(point.y,3.);
        close(point.x.hypot(point.z),1.);
    }
}

#[test]
fn physical_face_sense_distinguishes_analytic_support_from_material_normal(){
    let b=blind();let m=b.model();
    assert_eq!(m.face_sense(FaceId(6)),Some(rustsolid::FaceSense::Reversed));
    assert_eq!(m.face_sense(FaceId(7)),Some(rustsolid::FaceSense::Reversed));
    assert_eq!(m.face_sense(FaceId(1)),Some(rustsolid::FaceSense::Forward));
    let cylinder=m.oriented_face_normal(FaceId(6),0.0,1.0).unwrap();
    let floor=m.oriented_face_normal(FaceId(7),0.0,0.0).unwrap();
    let opening=m.oriented_face_normal(FaceId(1),0.0,0.0).unwrap();
    close(cylinder.x,-1.0);
    close(cylinder.y,0.0);
    close(floor.y,1.0);
    close(opening.y,1.0);
    let q=rustsolid::geometry_query(m,
        "query {body {faces(id:7) { id sense surface { kind } }}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(q["body"]["faces"][0]["sense"],"reversed");
    assert_eq!(q["body"]["faces"][0]["surface"]["kind"],"plane");
}
#[test]
fn existing_through_hole_and_standalone_cylinder_have_correct_material_face_sense(){
    let through=ThroughHoleBrep::from_box(&box_(),HoleParameters{
        center_x:0.,center_z:0.,radius:1.}).unwrap();
    assert_eq!(through.model().face_sense(FaceId(6)),
        Some(rustsolid::FaceSense::Reversed));
    let inward=through.model().oriented_face_normal(FaceId(6),0.0,1.0).unwrap();
    close(inward.x,-1.);
    let stand=CylindricalBrep::upright(p(0.,0.,0.),1.,2.,
        GeometryTolerance::default()).unwrap();
    let canonical=BrepBody::from(stand).shared().unwrap();
    assert_eq!(canonical.face_sense(FaceId(2)),
        Some(rustsolid::FaceSense::Forward));
    let outward=canonical.oriented_face_normal(FaceId(2),0.0,1.0).unwrap();
    close(outward.x,1.);
}

#[test]
fn canonical_brep_validator_rejects_tampered_blind_floor_and_side_orientation(){
    let b=blind();
    let mut changed=b.model().clone();
    changed.faces[7].sense=rustsolid::FaceSense::Forward;
    assert!(changed.validate().is_err());
    let mut changed=b.model().clone();
    changed.faces[6].sense=rustsolid::FaceSense::Forward;
    assert!(changed.validate().is_err());
    b.validate().unwrap();
}
#[test]
fn carrier_uv_edits_preserve_material_face_sense_on_blind_floor(){
    let body=BrepBody::from(blind());
    let mut editor=rustsolid::BrepCarrierEditor::from_body(&body).unwrap();
    let original=editor.snapshot();
    let floor=editor.handle(rustsolid::CarrierOwner::Face(FaceId(7))).unwrap();
    let side=editor.handle(rustsolid::CarrierOwner::Face(FaceId(6))).unwrap();
    let mut tx=editor.begin_transaction();
    tx.rebase_face_uv(floor,0.15,-0.3).unwrap();
    tx.rebase_face_uv(side,0.0,0.2).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();
    assert_eq!(editor.model().face_sense(FaceId(7)),
        Some(rustsolid::FaceSense::Reversed));
    assert_eq!(editor.model().face_sense(FaceId(6)),
        Some(rustsolid::FaceSense::Reversed));
    close(editor.model().oriented_face_normal(FaceId(7),0.0,0.0).unwrap().y,1.0);
    close(editor.model().oriented_face_normal(FaceId(6),0.0,0.5).unwrap().x,-1.0);
    assert_eq!(editor.model().mass,original.model.mass);
    editor.undo().unwrap().unwrap();
    assert_eq!(editor.model().face_sense(FaceId(7)),
        Some(rustsolid::FaceSense::Reversed));
}
