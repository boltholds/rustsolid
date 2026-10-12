//! Executable intersection topology witnesses: two overlapping but
//! nonnested equal-depth blind cylinders are joined along true crossing
//! vertical edges with two exposed analytical circle arcs per cap.
use rustsolid::{
    BlindHoleParameters, CrossedDiskRegion, Curve3, EdgeId, FaceId,
    FaceSense, IntersectingBlindBrep, LoopRole, Point3, Solid, Surface3,
    CoedgeId, ModelUnits,
};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn hole(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn shape()->IntersectingBlindBrep{
    IntersectingBlindBrep::from_box(&block(),
        [hole(-0.7,0.0,1.0,2.0),hole(0.7,0.0,1.0,2.0)]).unwrap()
}
fn eq(a:f64,b:f64,eps:f64){
    assert!((a-b).abs()<=eps,"{a} vs {b}, difference {}",(a-b).abs());
}
fn certify(s:&IntersectingBlindBrep){
    s.validate().unwrap();
    let b=s.model();b.validate().unwrap();
    let summary=b.summary();
    assert_eq!((summary.vertices,summary.edges,summary.coedges),
        (12,18,36));
    assert_eq!((summary.loops,summary.faces,summary.shells),(10,9,1));
    assert_eq!(summary.euler_characteristic,2);
    assert_eq!(summary.genus,0);
    assert_eq!(summary.circular_edges,4);
    assert_eq!(summary.cylindrical_faces,2);
    assert_eq!(summary.seam_edges,0,"crossing is a pair of real shared edges, not a periodic seam");
    assert_eq!(b.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),1);
    assert_eq!(b.faces[0].loops.len(),1);
    assert_eq!(b.faces[1].loops.len(),2);
    assert_eq!(b.faces[8].loops.len(),1);
    for id in 6..=8{
        assert_eq!(b.face_sense(FaceId(id)),Some(FaceSense::Reversed));
    }
    assert!(matches!(b.face_surface(FaceId(8)),Some(Surface3::Plane(_))));
    assert!(matches!(b.face_surface(FaceId(6)),Some(Surface3::Cylinder(_))));
    assert!(matches!(b.face_surface(FaceId(7)),Some(Surface3::Cylinder(_))));
    for id in 12..16{
        assert!(matches!(b.edge_curve(EdgeId(id)),Some(Curve3::Circle(_))));
        let g=b.geometry.edges[id as usize];
        let rustsolid::ParameterRange::Bounded{start,end}=g.domain
            else{panic!("split circle arc should have a bounded domain")};
        assert!(end>start);
        assert!(end-start<2.*PI);
    }
    for id in 16..18{
        assert!(matches!(b.edge_curve(EdgeId(id)),Some(Curve3::Line(_))));
    }
    assert_eq!(b.coedges[32].twin,CoedgeId(33));
    assert_eq!(b.coedges[34].twin,CoedgeId(35));
    assert_eq!(b.loops[6].face,FaceId(1));
    assert_eq!(b.loops[7].face,FaceId(8));
}
#[test]
fn true_crossing_reconstructs_one_connected_two_arc_opening_and_floor(){
    let s=shape();
    certify(&s);
    let r=s.region();
    assert_eq!(r.arc_a.end_crossing,r.arc_b.start_crossing);
    assert_eq!(r.arc_b.end_crossing,r.arc_a.start_crossing);
    assert!(r.arc_a.arc_length()>0.0);
    assert!(r.arc_b.arc_length()>0.0);
    let b=block();
    eq(s.model().mass.volume,b.mass.volume-r.disk_union_area*2.0,1e-9);
    eq(s.model().mass.surface_area,
        b.mass.surface_area+r.exposed_perimeter*2.0,1e-9);
    eq(s.model().mass.centroid.y,
        1.0-(r.disk_union_area*2.0)*(2.0-1.0)/
            (144.0-r.disk_union_area*2.0),1e-9);
}
#[test]
fn circle_crossing_graph_has_two_real_vertical_3d_edges_and_shared_vertices(){
    let s=shape();
    let m=s.model();let region=s.region();
    let floor=3.0-2.0;
    for i in 0..2{
        let q=region.crossings_xz[if i==0{region.arc_a.start_crossing}else{region.arc_a.end_crossing}];
        let bottom=m.vertices[8+i].position;
        let top=m.vertices[10+i].position;
        eq(bottom.x,q[0],1e-9);eq(bottom.z,q[1],1e-9);
        eq(bottom.y,floor,1e-9);eq(top.y,3.0,1e-9);
        let seam=m.edge_curve(EdgeId((16+i) as u32)).unwrap();
        let g=m.geometry.edges[16+i];
        let rustsolid::ParameterRange::Bounded{start,end}=g.domain else{panic!()};
        let p0=seam.evaluate(start).unwrap();
        let p1=seam.evaluate(end).unwrap();
        eq(p0.x,bottom.x,1e-9);eq(p0.z,bottom.z,1e-9);
        eq(p1.x,top.x,1e-9);eq(p1.y,top.y,1e-9);
    }
}
#[test]
fn partially_overlapped_equal_circles_have_exact_classical_lens_area(){
    let r=shape().region().clone();
    let c1:f64=1.0;
    let d:f64=1.4;
    let x=d/2.0;
    let overlap=2.0*(c1*c1*(x/c1).acos()-x*(c1*c1-x*x).sqrt());
    eq(r.disk_union_area,2.0*PI-overlap,1e-12);
}
#[test]
fn unequal_radii_and_rotated_centers_produce_exact_genus_zero_shell(){
    for holes in [
        [hole(-0.75,0.,1.4,1.5),hole(0.8,0.,0.8,1.5)],
        [hole(-0.5,-0.5,1.25,2.5),hole(0.9,0.65,0.95,2.5)],
        [hole(0.,-0.55,1.1,1.0),hole(0.,0.75,1.0,1.0)],
    ]{
        let s=IntersectingBlindBrep::from_box(&block(),holes).unwrap();
        certify(&s);
        assert!(s.region().disk_union_area>PI*holes[0].radius*holes[0].radius);
        assert!(s.region().exposed_perimeter>0.0);
        let m=s.model();
        eq(m.mass.volume,144.0-s.region().disk_union_area*holes[0].depth,1e-9);
    }
}
#[test]
fn swapping_operands_preserves_boolean_volume_area_and_centroid(){
    let a=hole(-0.5,-0.35,1.25,2.5);
    let b=hole(0.8,0.5,1.0,2.5);
    let first=IntersectingBlindBrep::from_box(&block(),[a,b]).unwrap();
    let second=IntersectingBlindBrep::from_box(&block(),[b,a]).unwrap();
    certify(&first);certify(&second);
    eq(first.model().mass.volume,second.model().mass.volume,1e-9);
    eq(first.model().mass.surface_area,second.model().mass.surface_area,1e-9);
    eq(first.model().mass.centroid.x,second.model().mass.centroid.x,1e-9);
    eq(first.model().mass.centroid.y,second.model().mass.centroid.y,1e-9);
    eq(first.model().mass.centroid.z,second.model().mass.centroid.z,1e-9);
}
#[test]
fn equal_depth_and_strict_transverse_crossings_are_required(){
    let a=hole(-0.7,0.,1.,2.);
    for b in [
        hole(0.7,0.,1.,2.1),
        hole(0.7,0.,1.,2.+1e-8),
        hole(1.3,0.,1.,2.), // nearly tangent outside
        hole(-0.7,0.,1.,2.), // coincident supports
        hole(-0.6,0.,0.3,2.), // contained circle
        hole(2.5,0.,1.,2.), // separated
        hole(0.7,0.,3.0,2.), // outside box
    ]{
        assert!(IntersectingBlindBrep::from_box(&block(),[a,b]).is_err(),
            "invalid pair {a:?}, {b:?}");
    }
}
#[test]
fn failed_late_crossing_region_commit_does_not_publish_any_geometry(){
    let mut s=shape();
    let original=s.model().mass;
    let p=s.parameters();
    let mut tx=s.begin_edit();
    tx.replace(0,hole(-0.8,0.,1.1,2.)).unwrap();
    tx.replace(1,hole(0.8,0.,1.0,1.5)).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(s.revision(),0);
    assert_eq!(s.undo_depth(),0);
    assert_eq!(s.model().mass,original);
    assert_eq!(s.parameters(),p);
    certify(&s);
}
#[test]
fn poisoned_stage_must_abort_even_after_a_valid_previous_stage(){
    let mut s=shape();
    let mut tx=s.begin_edit();
    tx.replace(0,hole(-0.8,0.,1.1,2.)).unwrap();
    assert!(tx.replace(2,hole(0.,0.,1.,2.)).is_err());
    assert!(tx.commit().is_err());
    certify(&s);
    assert_eq!(s.revision(),0);
    assert_eq!(s.undo_depth(),0);
}
#[test]
fn complete_analytic_region_rebuild_is_one_undo_step(){
    let mut s=shape();
    let before=s.model().mass;
    let mut tx=s.begin_edit();
    tx.replace(0,hole(-0.85,0.05,1.05,2.25)).unwrap();
    tx.replace(1,hole(0.7,-0.05,0.95,2.25)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    certify(&s);
    assert_eq!(s.undo_depth(),1);
    assert_ne!(before.volume,s.model().mass.volume);
    assert!(s.undo().unwrap());
    assert_eq!(s.revision(),2);
    eq(s.model().mass.volume,before.volume,1e-9);
    certify(&s);
    assert!(s.redo().unwrap());
    assert_eq!(s.revision(),3);
    certify(&s);
}
#[test]
fn canonical_query_can_inspect_split_circle_arc_carriers_and_face_senses(){
    let s=shape();
    let json=rustsolid::geometry_query(s.model(),
        "query {body {kind topology {faceCount edgeCount genus circularEdges} faces(id:8) {id sense surface {kind}} }}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(json["body"]["kind"],"analytic_intersecting_blind_holes");
    assert_eq!(json["body"]["topology"]["faceCount"],9);
    assert_eq!(json["body"]["topology"]["edgeCount"],18);
    assert_eq!(json["body"]["topology"]["genus"],0);
    assert_eq!(json["body"]["faces"][0]["sense"],"reversed");
}
#[test]
fn translation_keeps_local_disk_mass_and_relative_centroid_consistent(){
    let b=block();
    let s=shape();
    let d=1e4;
    let translated=Solid::block(p(d-3.,d-1.,d-3.),6.,4.,6.).unwrap();
    let shape2=IntersectingBlindBrep::from_box(&translated,[
        hole(d-0.7,d,1.,2.),hole(d+0.7,d,1.,2.)
    ]).unwrap();
    certify(&shape2);
    eq(s.model().mass.volume,shape2.model().mass.volume,1e-8);
    eq(s.model().mass.centroid.x+d,shape2.model().mass.centroid.x,1e-8);
    eq(s.model().mass.centroid.y+d,shape2.model().mass.centroid.y,1e-8);
    eq(s.model().mass.centroid.z+d,shape2.model().mass.centroid.z,1e-8);
    eq(b.mass.volume,translated.mass.volume,1e-9);
}

#[test]
fn common_brep_body_and_region_membership_cover_union_void_and_crossed_floor(){
    use rustsolid::{BrepBody,BooleanRegionClassifier,VolumePosition};
    let body=BrepBody::from(shape());
    body.validate().unwrap();
    let model=body.shared().unwrap();
    assert_eq!(model.summary().faces,9);
    let mesh=body.tessellate(32).unwrap();
    assert!(!mesh.triangles.is_empty());
    for (point,expected) in [
        (p(-0.7,2.,0.),VolumePosition::Outside),
        (p(0.7,2.,0.),VolumePosition::Outside),
        (p(0.,2.,0.),VolumePosition::Outside),
        (p(0.,1.,0.),VolumePosition::Boundary),
        (p(0.,0.5,0.),VolumePosition::Inside),
        (p(-1.7,2.,0.),VolumePosition::Boundary),
        (p(2.5,2.,0.),VolumePosition::Inside),
        (p(0.,3.,0.),VolumePosition::Outside),
        (p(0.,-1.,0.),VolumePosition::Boundary),
        (p(4.,0.,0.),VolumePosition::Outside),
    ] {
        assert_eq!(BooleanRegionClassifier::classify_point(&body,point).unwrap(),
            expected,"incorrect union-bore material at {point:?}");
    }
}
#[test]
fn arena_reassembles_crossing_cylinders_without_losing_semantic_handles(){
    use rustsolid::{BrepArena,TopologyEntity};
    let mut arena=BrepArena::from_intersecting_blind(shape()).unwrap();
    let snapshot=arena.snapshot().unwrap();
    let floor=arena.handle(TopologyEntity::Face(FaceId(8))).unwrap();
    let side=arena.handle(TopologyEntity::Face(FaceId(6))).unwrap();
    let initial=arena.model().unwrap().mass.volume;
    let next=[
        hole(-0.8,0.,1.1,2.5),
        hole(0.7,0.,0.95,2.5),
    ];
    assert_eq!(arena.update_intersecting_blind(next).unwrap(),1);
    arena.validate().unwrap();
    assert!(arena.check_snapshot(&snapshot).is_err());
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(8))),Some(floor));
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(6))),Some(side));
    assert_eq!(arena.undo_depth(),1);
    assert_ne!(arena.model().unwrap().mass.volume,initial);
    assert!(arena.update_intersecting_blind([
        hole(-0.8,0.,1.1,2.5),
        hole(0.7,0.,0.95,2.4),
    ]).is_err());
    assert_eq!(arena.revision(),1);
    assert!(arena.undo().unwrap().is_some());
    eq(arena.model().unwrap().mass.volume,initial,1e-9);
    assert!(arena.redo().unwrap().is_some());
    assert_eq!(arena.revision(),3);
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(8))),Some(floor));
    arena.validate().unwrap();
}
#[test]
fn detached_analytic_carrier_editor_preserves_crossed_uv_trim_provenance(){
    use rustsolid::{BrepBody,BrepCarrierEditor,CarrierOwner,CoedgeId};
    let input=BrepBody::from(shape());
    let mut editor=BrepCarrierEditor::from_body(&input).unwrap();
    let model=editor.model().clone();
    let curve=editor.handle(CarrierOwner::Edge(EdgeId(13))).unwrap();
    let side=editor.handle(CarrierOwner::Face(FaceId(6))).unwrap();
    let mut tx=editor.begin_transaction();
    tx.reparameterize_edge(curve,1.0,0.3).unwrap();
    tx.rebase_face_uv(side,0.0,0.15).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();
    eq(editor.model().mass.volume,model.mass.volume,1e-9);
    assert_eq!(editor.model().coedges[33].twin,CoedgeId(32));
    assert_eq!(editor.model().face_sense(FaceId(6)),
        Some(rustsolid::FaceSense::Reversed));
}

#[test]
fn independent_simpson_cross_section_oracle_confirms_unequal_disk_area_and_centroid(){
    // Integrate union cross-section heights directly, without evaluating
    // circular-segment formulae used by the analytical kernel.
    let a=hole(-0.65,-0.55,1.25,2.3);
    let b=hole(0.8,0.5,0.95,2.3);
    let shape=IntersectingBlindBrep::from_box(&block(),[a,b]).unwrap();
    let r=shape.region();
    let dx=b.center_x-a.center_x;
    let dz=b.center_z-a.center_z;
    let d=dx.hypot(dz);let ux=dx/d;let uz=dz/d;
    let n=60000usize;
    let lower=-a.radius;
    let upper=d+b.radius;
    let step=(upper-lower)/n as f64;
    let mut area=0.0;
    let mut axial_moment=0.0;
    for i in 0..=n {
        let x=lower+i as f64*step;
        let h1=(a.radius*a.radius-x*x).max(0.0).sqrt();
        let h2=(b.radius*b.radius-(x-d)*(x-d)).max(0.0).sqrt();
        let vertical=2.0*h1.max(h2);
        let weight=if i==0 || i==n {1.0}
            else if i%2==0 {2.0}else{4.0};
        area+=weight*vertical;
        axial_moment+=weight*vertical*x;
    }
    area*=step/3.0;
    axial_moment*=step/3.0;
    let absolute_axial=area*(a.center_x*ux+a.center_z*uz)+axial_moment;
    let analytical_axial=r.relative_area_moment[0]*ux+
        r.relative_area_moment[1]*uz;
    eq(area,r.disk_union_area,8e-5);
    eq(absolute_axial,analytical_axial,8e-5);
    let px=-uz;let pz=ux;
    let perpendicular=r.relative_area_moment[0]*px+
        r.relative_area_moment[1]*pz;
    eq(perpendicular,area*(a.center_x*px+a.center_z*pz),8e-5);
}
