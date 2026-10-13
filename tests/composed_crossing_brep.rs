//! First-party topology, depth-slab and atomicity checks for v0.33.
use rustsolid::{
    BlindHoleBrep,BlindHoleParameters,CertifiedCrossingPair,
    ComposedCrossingBrep,ComposedCrossingParameters,
    HoleParameters,LoopRole,MixedCrossingBrep,MultiCutPlanState,
    Point3,Solid,Surface3,FaceSense,
};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-5.,-1.,-4.),10.,4.,8.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn pair(x:f64,z:f64,r:f64,d:f64)->CertifiedCrossingPair{
    CertifiedCrossingPair{through:t(x,z,r),blind:b(x+r,z,r,d)}
}
fn params(n:usize,isolated:bool)->ComposedCrossingParameters{
    let pairs=[pair(-2.,0.,0.45,1.5),pair(1.5,0.,0.5,2.5)];
    ComposedCrossingParameters{
        pairs:pairs[..n].to_vec(),
        isolated_blind:if isolated{vec![b(0.,2.,0.3,1.0)]}else{vec![]},
    }
}
fn close(a:f64,b:f64,eps:f64){
    assert!((a-b).abs()<=eps,"value {a:.15}, expected {b:.15}, residual {}",(a-b).abs())
}
#[test]
fn two_source_certified_crossing_clusters_share_one_closed_genus_two_shell(){
    let source=box_();
    let parameters=params(2,false);
    let body=ComposedCrossingBrep::from_box(&source,parameters.clone()).unwrap();
    body.validate().unwrap();
    let m=body.model();let st=m.summary();
    assert_eq!((st.vertices,st.edges,st.coedges,st.loops,st.faces,st.shells),
        (20,34,68,20,16,1));
    assert_eq!((st.euler_characteristic,st.genus),(-2,2));
    assert_eq!(m.faces[0].loops.len(),3);
    assert_eq!(m.faces[1].loops.len(),3);
    assert_eq!(m.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),4);
    assert!(m.shells[0].closed);
    for f in &m.faces[6..] {
        assert_eq!(f.sense,FaceSense::Reversed);
        assert!(matches!(m.face_surface(f.id),Some(Surface3::Cylinder(_)|Surface3::Plane(_))));
    }
    for edge in &m.edges[12..]{
        let a=&m.coedges[edge.coedges[0].0 as usize];
        let b=&m.coedges[edge.coedges[1].0 as usize];
        assert_ne!(a.reversed,b.reversed);
        assert_eq!(a.twin,b.id);assert_eq!(b.twin,a.id);
    }
    let direct=parameters.pairs.iter().map(|pair|
        MixedCrossingBrep::from_box(&source,pair.through,pair.blind).unwrap())
        .collect::<Vec<_>>();
    let expected=source.mass.volume-
        direct.iter().map(|d|source.mass.volume-d.model().mass.volume).sum::<f64>();
    close(m.mass.volume,expected,1e-10);
    close(m.mass.volume,body.slab_plan().predicted_remaining_volume.unwrap(),1e-10);
    assert_eq!(body.slab_plan().state,MultiCutPlanState::Ready);
    // Two distinct blind depths partition full height into THREE slabs.
    assert_eq!(body.slab_plan().slabs.len(),3);
}

#[test]
fn one_crossing_plus_distinct_depth_isolated_blind_produces_one_shared_shell(){
    let source=box_();let parameters=params(1,true);
    let body=ComposedCrossingBrep::from_box(&source,parameters.clone()).unwrap();
    body.validate().unwrap();
    let m=body.model();let st=m.summary();
    assert_eq!((st.vertices,st.edges,st.coedges,st.loops,st.faces),(16,26,52,16,13));
    assert_eq!((st.euler_characteristic,st.genus),(0,1));
    assert_eq!(m.faces[0].loops.len(),2);
    assert_eq!(m.faces[1].loops.len(),3);
    assert_eq!(m.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),3);
    assert_eq!(body.slab_plan().slabs.len(),3);
    let crossed=MixedCrossingBrep::from_box(&source,
        parameters.pairs[0].through,parameters.pairs[0].blind).unwrap();
    let solo=BlindHoleBrep::from_box(&source,parameters.isolated_blind[0]).unwrap();
    let v=source.mass.volume-
        (source.mass.volume-crossed.model().mass.volume)-
        (source.mass.volume-solo.model().mass.volume);
    close(m.mass.volume,v,1e-10);
    close(m.mass.volume,body.slab_plan().predicted_remaining_volume.unwrap(),1e-10);
}

#[test]
fn B_rep_mesh_is_watertight_with_two_global_apertures_and_source_provenance(){
    let body=ComposedCrossingBrep::from_box(&box_(),params(2,false)).unwrap();
    let d=body.tessellate(32).unwrap();
    assert_eq!(d.genus,2);
    assert_eq!(d.crossings,2);
    assert_eq!(d.isolated_blind,0);
    assert!(d.mesh.vertices.len()>150);
    assert_eq!(d.mesh.triangles.len(),d.mesh.triangle_faces.len());
    assert!(d.mesh.triangle_faces.iter().any(|x|x.0==0));
    assert!(d.mesh.triangle_faces.iter().any(|x|x.0==1));
    assert!(d.mesh.triangle_faces.iter().any(|x|x.0==15));
    assert_eq!(d.mesh.vertices.len() as isize-d.triangle_edges as isize+
        d.mesh.triangles.len() as isize,-2);
    close(d.analytic_volume,body.model().mass.volume,1e-10);
    close(d.facet_volume,d.expected_polygonal_volume,1e-7);
    assert!(d.volume_faceting_error>=-1e-7);
}

#[test]
fn B_rep_mesh_supports_crossed_pair_plus_unconnected_blind_floor(){
    let body=ComposedCrossingBrep::from_box(&box_(),params(1,true)).unwrap();
    let d=body.tessellate(32).unwrap();
    assert_eq!((d.genus,d.crossings,d.isolated_blind),(1,1,1));
    assert!(d.mesh.triangle_faces.iter().any(|x|x.0==12));
    close(d.facet_volume,d.expected_polygonal_volume,1e-7);
}

#[test]
fn analytic_volume_is_unchanged_by_meshing_and_facet_error_decreases(){
    let body=ComposedCrossingBrep::from_box(&box_(),params(2,false)).unwrap();
    let initial=body.model().mass;
    let low=body.tessellate(16).unwrap();
    let mid=body.tessellate(32).unwrap();
    let high=body.tessellate(64).unwrap();
    assert_eq!(body.model().mass,initial);
    assert!(low.volume_faceting_error>=mid.volume_faceting_error-1e-7);
    assert!(mid.volume_faceting_error>=high.volume_faceting_error-1e-7);
    assert!(low.mesh.triangles.len()<high.mesh.triangles.len());
}

#[test]
fn atomic_source_edits_and_undo_redo_reconstruct_whole_shell(){
    let mut body=ComposedCrossingBrep::from_box(&box_(),params(2,true)).unwrap();
    let before=body.model().mass;
    assert_eq!(body.revision(),0);
    let mut tx=body.begin_edit();
    tx.replace_pair(0,pair(-2.1,0.1,0.46,1.6)).unwrap();
    tx.replace_isolated_blind(0,b(0.,2.,0.35,1.2)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    let edited=body.model().mass;
    assert_ne!(edited.volume,before.volume);
    body.validate().unwrap();
    assert!(body.undo().unwrap());
    assert_eq!(body.revision(),2);
    close(body.model().mass.volume,before.volume,1e-11);
    assert!(body.redo().unwrap());
    assert_eq!(body.revision(),3);
    close(body.model().mass.volume,edited.volume,1e-11);
    assert_eq!(body.tessellate(24).unwrap().genus,2);
}

#[test]
fn near_touch_and_connected_three_cutter_interactions_fail_closed(){
    let source=box_();
    let mut parameters=params(1,true);
    // Extra blind disk touches the crossing's through footprint.
    parameters.isolated_blind[0]=b(-2.,0.85,0.4,1.0);
    assert!(ComposedCrossingBrep::from_box(&source,parameters).is_err());
    let mut p=params(2,false);
    p.pairs[1]=pair(-1.1,0.,0.5,2.5);
    assert!(ComposedCrossingBrep::from_box(&source,p).is_err());
    let mut p=params(2,true);
    p.isolated_blind[0].radius=f64::NAN;
    assert!(ComposedCrossingBrep::from_box(&source,p).is_err());
}
#[test]
fn strong_atomic_failure_preserves_mass_history_and_all_carriers(){
    let mut body=ComposedCrossingBrep::from_box(&box_(),params(2,true)).unwrap();
    let before=body.model().mass;
    let mut tx=body.begin_edit();
    tx.replace_pair(0,pair(-2.,0.,0.48,1.4)).unwrap();
    tx.replace_isolated_blind(0,b(-2.,0.0,0.4,1.0)).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(body.revision(),0);
    assert_eq!(body.model().mass,before);
    assert_eq!((body.undo_depth(),body.redo_depth()),(0,0));
    body.validate().unwrap();
    let mut tx=body.begin_edit();
    assert!(tx.replace_pair(99,pair(-2.,0.,0.4,1.0)).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(body.revision(),0);
}
