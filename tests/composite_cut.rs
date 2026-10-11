//! Independent shape, topology, volume, mixed planar-face region and
//! feature lifecycle checks for compositional blind/through analytical cuts.
use rustsolid::{
    AnalyticCut, BlindHoleParameters, BrepCompositeEditor, BrepBody,
    BrepArena, BooleanRegionClassifier, CompositeHandleError, FaceId,
    FaceSense, HoleParameters, Point3, Solid, TopologyEntity, VolumePosition,
    GeometryTolerance, HistoryDirection, Curve3, Surface3, 
};
use std::collections::BTreeMap;
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn through(x:f64,z:f64,r:f64)->AnalyticCut{
    AnalyticCut::Through(HoleParameters{center_x:x,center_z:z,radius:r})
}
fn blind(x:f64,z:f64,r:f64,d:f64)->AnalyticCut{
    AnalyticCut::Blind(BlindHoleParameters{
        center_x:x,center_z:z,radius:r,depth:d
    })
}
fn mix()->Vec<AnalyticCut>{vec![through(-1.5,0.,0.6),blind(1.4,0.3,0.75,1.5)]}
fn approx(a:f64,b:f64,t:f64){
    assert!((a-b).abs()<=t,"{a} versus {b}, err={} budget {t}",(a-b).abs());
}
fn model(cuts:&[AnalyticCut])->BrepCompositeEditor{
    BrepCompositeEditor::from_cuts(box_(),cuts).unwrap()
}
fn prove(e:&BrepCompositeEditor,through:usize,blind:usize){
    e.validate().unwrap();
    let m=e.model();
    let s=m.summary();let n=through+blind;
    assert_eq!(s.vertices,8+2*n);
    assert_eq!(s.edges,12+3*n);
    assert_eq!(s.coedges,24+6*n);
    assert_eq!(s.loops,6+3*n);
    assert_eq!(s.faces,6+through+2*blind);
    assert_eq!(s.shells,1);
    assert_eq!(s.circular_edges,2*n);
    assert_eq!(s.cylindrical_faces,n);
    assert_eq!(s.seam_edges,n);
    assert_eq!(s.euler_characteristic,2-2*through as isize);
    assert_eq!(s.genus,through);
    assert_eq!(m.loops.iter().filter(|lp|lp.role==rustsolid::LoopRole::Inner)
        .count(),2*through+blind);
    assert_eq!(m.faces[0].loops.len(),1+through);
    assert_eq!(m.faces[1].loops.len(),1+n);
    for face in 6..m.faces.len(){
        assert_eq!(m.face_sense(FaceId(face as u32)),Some(FaceSense::Reversed));
    }
    let mut removed=0.0;
    let mut surface=168.0;
    let mut moment=p(0.,0.,0.);
    for binding in e.bindings(){
        let c=binding.cut;let [x,z]=c.center();
        let d=c.depth(e.source_box());
        let volume=PI*c.radius().powi(2)*d;
        removed+=volume;
        surface+=2.0*PI*c.radius()*d;
        if c.is_through(){surface-=2.0*PI*c.radius().powi(2);}
        let y=if c.is_through(){1.0}else{3.-d/2.};
        moment=moment.add(p(-volume*x,volume*(1.-y),-volume*z));
    }
    let v=144.-removed;
    approx(m.mass.volume,v,1e-9);
    approx(m.mass.surface_area,surface,1e-9);
    approx(m.mass.centroid.x,moment.x/v,1e-9);
    approx(m.mass.centroid.y,1.+moment.y/v,1e-9);
    approx(m.mass.centroid.z,moment.z/v,1e-9);
}
fn prove_mesh(e:&BrepCompositeEditor,segments:usize){
    let result=e.tessellate(segments).unwrap();
    let n=e.cut_count();
    let cuts=e.bindings();
    let t=cuts.iter().filter(|b|b.cut.is_through()).count();
    let b=n-t;
    assert_eq!(result.mesh.vertices.len(),8+2*n*segments+b);
    let bottom=2+t*(segments+2);
    let top=2+n*(segments+2);
    assert_eq!(result.bottom.loops,1+t);
    assert_eq!(result.top.loops,1+n);
    assert_eq!(result.bottom.triangles,bottom);
    assert_eq!(result.top.triangles,top);
    assert_eq!(result.mesh.triangles.len(),
        bottom+top+8+2*n*segments+b*segments);
    assert_eq!(result.mesh.triangle_faces.len(),result.mesh.triangles.len());
    let mut owners=BTreeMap::<u32,usize>::new();
    let mut uses=BTreeMap::<(u32,u32),Vec<bool>>::new();
    let mut six_volume=0.0;
    let a=e.source_box().bbox.min;
    for (tri,face) in result.mesh.triangles.iter()
        .zip(&result.mesh.triangle_faces){
        let [i,j,k]=*tri;
        *owners.entry(face.0).or_default()+=1;
        for (x,y) in [(i,j),(j,k),(k,i)]{
            uses.entry((x.min(y),x.max(y))).or_default().push(x<y);
        }
        let [v0,v1,v2]=tri.map(|u|result.mesh.vertices[u as usize]);
        six_volume+=v0.sub(a).dot(v1.sub(a).cross(v2.sub(a)));
    }
    assert!(uses.values().all(|v|v.len()==2&&v[0]!=v[1]));
    assert_eq!(owners[&0],bottom);
    assert_eq!(owners[&1],top);
    for face in 2..6{assert_eq!(owners[&face],2);}
    for cut in cuts{
        assert_eq!(owners[&cut.topology.side_face.0],2*segments);
        if let Some(floor)=cut.topology.floor_face{
            assert_eq!(owners[&floor.0],segments);
        }
    }
    approx(six_volume/6.,result.faceted_volume,1e-8);
    approx(result.analytic_volume,e.model().mass.volume,1e-9);
    assert!(result.faceting_volume_error>=0.0);
}
#[test]
fn one_through_and_one_blind_share_caps_with_genuine_disk_floor(){
    let e=model(&mix());
    prove(&e,1,1);
    let x=e.bindings();
    assert_eq!(x[0].topology.side_face,FaceId(6));
    assert!(x[0].topology.floor_face.is_none());
    assert_eq!(x[1].topology.side_face,FaceId(7));
    assert_eq!(x[1].topology.floor_face,Some(FaceId(8)));
    assert!(matches!(e.model().face_surface(FaceId(8)),Some(Surface3::Plane(_))));
    assert!(matches!(e.model().edge_curve(x[1].topology.lower_edge),
        Some(Curve3::Circle(_))));
    assert_eq!(e.model().loops[x[1].topology.lower_loop.0 as usize].face,FaceId(8));
}
#[test]
fn analytic_exact_mass_is_independent_of_cut_insertion_order(){
    let a=model(&mix());
    let mut other=mix();other.reverse();
    let b=model(&other);
    prove(&a,1,1);prove(&b,1,1);
    approx(a.model().mass.volume,b.model().mass.volume,1e-9);
    approx(a.model().mass.surface_area,b.model().mass.surface_area,1e-9);
    approx(a.model().mass.centroid.x,b.model().mass.centroid.x,1e-9);
    approx(a.model().mass.centroid.y,b.model().mass.centroid.y,1e-9);
}
#[test]
fn zero_and_homogeneous_kinds_have_expected_euler_characteristics(){
    prove(&model(&[]),0,0);
    prove(&model(&[through(0.,0.,0.75)]),1,0);
    prove(&model(&[blind(0.,0.,0.75,1.5)]),0,1);
    prove(&model(&[through(-1.5,0.,0.5),through(1.5,0.,0.5)]),2,0);
    prove(&model(&[blind(-1.5,0.,0.5,1.5),blind(1.5,0.,0.5,2.25)]),0,2);
}
#[test]
fn both_mixed_cap_regions_and_true_blind_floor_are_watertight(){
    let e=model(&mix());prove_mesh(&e,48);
}
#[test]
fn homogeneous_through_and_blind_meshing_remains_closed(){
    prove_mesh(&model(&[through(0.,0.,1.)]),32);
    prove_mesh(&model(&[blind(0.,0.,1.,2.)]),32);
    prove_mesh(&model(&[]),32);
}
#[test]
fn mixed_faceting_error_converges_with_resolution_without_changing_analytic_mass(){
    let e=model(&mix());let original=e.snapshot();
    let a=e.tessellate(12).unwrap();
    let b=e.tessellate(48).unwrap();
    let c=e.tessellate(192).unwrap();
    assert!(a.faceting_volume_error>b.faceting_volume_error);
    assert!(b.faceting_volume_error>c.faceting_volume_error);
    assert!(c.faceting_volume_error<b.faceting_volume_error*0.08);
    for m in [&a,&b,&c]{
        approx(m.analytic_volume,original.model.mass.volume,1e-9);
    }
    e.check_snapshot(&original).unwrap();
}
#[test]
fn committed_mutation_and_undo_redo_change_genus_but_keep_live_handle(){
    let mut e=model(&mix());
    let first=e.handle(0).unwrap();let second=e.handle(1).unwrap();
    let mut tx=e.begin_transaction();
    tx.replace(second,through(1.4,0.3,0.75)).unwrap();
    let commit=tx.commit().unwrap();
    assert_eq!(commit.command_count,1);
    assert_eq!(e.resolve(first).unwrap().handle,first);
    assert_eq!(e.resolve(second).unwrap().handle,second);
    prove(&e,2,0);
    let transition=e.undo().unwrap().unwrap();
    assert_eq!(transition.direction,HistoryDirection::Undo);
    prove(&e,1,1);
    e.redo().unwrap().unwrap();
    prove(&e,2,0);
}
#[test]
fn deleting_earlier_cut_compacts_dense_topology_but_not_surviving_feature_handle(){
    let mut e=model(&mix());
    let early=e.handle(0).unwrap();let later=e.handle(1).unwrap();
    let prev=e.resolve(later).unwrap().topology;
    let mut tx=e.begin_transaction();
    tx.remove(early).unwrap();
    let result=tx.commit().unwrap();
    assert_eq!(result.retired.len(),1);
    assert_eq!(result.remapped.len(),1);
    assert_eq!(result.remapped[0].0,later);
    assert_eq!(result.remapped[0].1,prev);
    assert_eq!(e.resolve(later).unwrap().topology.side_face,FaceId(6));
    assert_eq!(e.resolve(early),Err(CompositeHandleError::StaleFeature));
    prove(&e,0,1);prove_mesh(&e,24);
}
#[test]
fn abandoned_staging_and_remove_replace_same_slot_do_not_allow_aba(){
    let mut e=model(&mix());
    let early=e.handle(0).unwrap();
    let replacement={
        let mut tx=e.begin_transaction();
        tx.remove(early).unwrap();
        let n=tx.add(blind(-1.5,0.,0.6,2.)).unwrap();
        assert_eq!(n.slot,early.slot);
        assert_ne!(n.generation,early.generation);
        tx.commit().unwrap();n
    };
    assert_eq!(e.resolve(early),Err(CompositeHandleError::StaleFeature));
    e.undo().unwrap().unwrap();
    let restored=e.handle(early.slot).unwrap();
    assert_ne!(restored,early);
    assert_ne!(restored,replacement);
    e.redo().unwrap().unwrap();
    assert_eq!(e.resolve(restored),Err(CompositeHandleError::StaleFeature));
    assert_eq!(e.resolve(replacement),Err(CompositeHandleError::StaleFeature));

}
#[test]
fn foreign_and_abandoned_cut_handles_fail_closed(){
    let mut a=BrepCompositeEditor::new(box_()).unwrap();
    let mut b=BrepCompositeEditor::new(box_()).unwrap();
    let abandoned={
        let mut t=a.begin_transaction();
        t.add(blind(0.,0.,1.,2.)).unwrap()
    };
    let actual={
        let mut t=a.begin_transaction();
        let id=t.add(blind(0.,0.,1.,2.)).unwrap();
        t.commit().unwrap();id
    };
    assert_eq!(a.resolve(abandoned),Err(CompositeHandleError::StaleFeature));
    let foreign={
        let mut t=b.begin_transaction();
        let h=t.add(blind(0.,0.,1.,2.)).unwrap();
        t.commit().unwrap();h
    };
    assert_eq!(a.resolve(foreign),Err(CompositeHandleError::ForeignEditor));
    let snapshot=a.snapshot();
    let mut t=a.begin_transaction();
    assert!(t.remove(foreign).is_err());
    assert!(t.commit().is_err());
    a.check_snapshot(&snapshot).unwrap();
}
#[test]
fn overlapping_or_tangent_cuts_are_rejected_without_publishing_partial_body(){
    let mut e=BrepCompositeEditor::new(box_()).unwrap();
    let before=e.snapshot();
    let mut tx=e.begin_transaction();
    tx.add(through(-1.,0.,0.8)).unwrap();
    assert!(tx.add(blind(-0.25,0.,0.8,1.5)).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(e.revision(),0);
    assert_eq!(e.undo_depth(),0);
    e.check_snapshot(&before).unwrap();
    prove(&e,0,0);
}
#[test]
fn mixed_cut_analytic_volume_material_classifier_knows_each_floor(){
    let body=BrepBody::from(model(&mix()));
    let samples=[
        (p(-1.5,1.,0.),VolumePosition::Outside),
        (p(1.4,2.5,0.3),VolumePosition::Outside),
        (p(1.4,1.5,0.3),VolumePosition::Boundary),
        (p(1.4,1.,0.3),VolumePosition::Inside),
        (p(2.8,2.,0.),VolumePosition::Inside),
        (p(-3.,0.,0.),VolumePosition::Boundary),
        (p(3.1,0.,0.),VolumePosition::Outside),
    ];
    for (point,expected) in samples {
        assert_eq!(BooleanRegionClassifier::classify_point(&body,point).unwrap(),
            expected,"at {point:?}");
    }
}
#[test]
fn read_only_brep_query_and_arena_detachment_preserve_source(){
    let e=model(&mix());let body=BrepBody::from(e.clone());
    let graph=body.shared().unwrap();
    assert_eq!(graph.summary().genus,1);
    let q=rustsolid::geometry_query(&graph,
        "query {body {kind topology {faceCount genus circularEdges} mass {volume}}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(q["body"]["kind"],"analytic_composite_cuts");
    assert_eq!(q["body"]["topology"]["faceCount"],9);
    assert_eq!(q["body"]["topology"]["genus"],1);
    let arena=BrepArena::from_composite(e.clone()).unwrap();
    let snapshot=arena.snapshot().unwrap();
    let face=arena.handle(TopologyEntity::Face(FaceId(8))).unwrap();
    let detached=arena.detached_composite_editor().unwrap();
    prove(&detached,1,1);
    arena.check_snapshot(&snapshot).unwrap();
    assert_eq!(arena.handle(TopologyEntity::Face(FaceId(8))),Some(face));
}
#[test]
fn invalid_segmentation_is_rejected_and_unmodified_revision_survives(){
    let e=model(&mix());
    let snap=e.snapshot();
    for n in [0,5,11,513]{
        assert!(e.tessellate(n).is_err());
    }
    e.check_snapshot(&snap).unwrap();
}
#[test]
fn far_world_coordinates_preserve_analytic_mass_and_mesh_closure(){
    let offset=1e8;
    let source=Solid::block(p(offset-3.,offset-1.,offset-3.),
        6.,4.,6.).unwrap();
    let editor=BrepCompositeEditor::from_cuts(source,&[
        through(offset-1.5,offset,0.6),
        blind(offset+1.4,offset+0.3,0.75,1.5),
    ]).unwrap();
    editor.validate().unwrap();
    prove_mesh(&editor,48);
    approx(editor.model().mass.volume,
        144.-PI*(0.6*0.6*4.+0.75*0.75*1.5),1e-9);
}
#[test]
fn arbitrary_cutting_axis_and_partial_contour_crossings_still_not_accepted(){
    let e=BrepCompositeEditor::from_cuts(box_(),&[blind(0.,0.,1.,2.)]).unwrap();
    assert_eq!(e.model().summary().genus,0);
    let mut tx=e.clone();let mut stage=tx.begin_transaction();
    assert!(stage.add(blind(0.,0.,0.8,1.)).is_err());
    assert!(stage.commit().is_err());
}
