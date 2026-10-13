//! Derived facet shell is checked independently of the committed analytic
//! triple-cylinder B-rep. Every physical source arc has globally shared
//! vertices across top/bottom caps, side walls and both blind floors.
use rustsolid::{
    BlindHoleParameters,ConnectedTripleBrep,HoleParameters,Point3,
    Solid,TripleChainParameters,
};
use std::collections::BTreeMap;
fn pos(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn source()->Solid{Solid::block(pos(-4.,-2.,-3.),8.,5.,6.).unwrap()}
fn initial()->TripleChainParameters{
    TripleChainParameters{
        through:HoleParameters{center_x:-1.15,center_z:0.0,radius:0.8},
        deeper_blind:BlindHoleParameters{
            center_x:0.0,center_z:0.0,radius:0.8,depth:2.8},
        shallower_blind:BlindHoleParameters{
            center_x:1.15,center_z:0.0,radius:0.8,depth:1.2}
    }
}
fn body()->ConnectedTripleBrep{
    ConnectedTripleBrep::from_box(&source(),initial()).unwrap()
}
fn close(a:f64,e:f64,eps:f64){
    assert!((a-e).abs()<=eps,"value {a:.14}, expected {e:.14}, error {:.3e}",(a-e).abs());
}
#[test]
fn one_connected_triple_mesh_has_opposite_twins_and_four_cap_floor_owners(){
    let b=body();
    let result=b.tessellate(32).unwrap();
    let m=&result.mesh;
    assert_eq!(result.genus,1);
    assert_eq!(result.segments_per_source_arc,32);
    assert!(result.faceted_volume>0.0);
    assert!(result.faceted_surface_area>0.0);
    assert_eq!(m.triangles.len(),m.triangle_faces.len());
    let mut face_count=[0usize;18];
    for id in &m.triangle_faces{
        face_count[id.0 as usize]+=1;
    }
    assert!(face_count.iter().all(|c|*c>0));
    assert!(face_count[0]>0&&face_count[1]>0);
    assert!(face_count[16]>0&&face_count[17]>0);
    let mut uses=BTreeMap::<(u32,u32),Vec<bool>>::new();
    for t in &m.triangles{
        for (a,b) in [(t[0],t[1]),(t[1],t[2]),(t[2],t[0])] {
            uses.entry((a.min(b),a.max(b))).or_default().push(a<b);
        }
    }
    assert!(uses.values().all(|x|x.len()==2&&x[0]!=x[1]));
    assert_eq!(m.vertices.len() as isize-uses.len() as isize+
        m.triangles.len() as isize,0);
    assert_eq!(result.triangle_edges,uses.len());
    close(result.faceted_volume,result.predicted_polygonized_volume,2e-6);
    close(result.faceted_surface_area,result.predicted_polygonized_area,2e-6);
    close(result.analytical_volume,b.model().mass.volume,1e-11);
    assert!(result.faceting_volume_error>=-2e-6);
}
#[test]
fn triple_mesh_refinement_preserves_source_mass_and_reduces_faceting_error(){
    let b=body();
    let before=b.model().mass;
    let mut errors=Vec::new();
    let mut triangles=Vec::new();
    for n in [12,24,48,96]{
        let m=b.tessellate(n).unwrap();
        close(m.analytical_volume,before.volume,1e-11);
        close(m.faceted_volume,m.predicted_polygonized_volume,2e-6);
        assert!(m.faceting_volume_error.is_finite());
        errors.push(m.faceting_volume_error.abs());
        triangles.push(m.mesh.triangles.len());
    }
    assert_eq!(b.model().mass,before);
    for pair in errors.windows(2){
        assert!(pair[1]<pair[0],
            "mesh refinement did not improve analytic error: {:?}",errors);
    }
    for pair in triangles.windows(2){assert!(pair[1]>pair[0]);}
    assert!(errors[3]<0.002);
}
#[test]
fn bad_display_budget_and_replay_transactions_leave_source_untouched(){
    let mut b=body();
    let first=b.tessellate(48).unwrap();
    for n in [0,1,7,11,257,10000]{
        assert!(b.tessellate(n).is_err());
        assert_eq!(b.revision(),0);
    }
    let mut next=initial();
    next.shallower_blind.center_x=1.2;
    next.shallower_blind.depth=1.4;
    b.begin_edit().replace_sources(next).unwrap();
    // Dropping a staged transaction never publishes a mutation.
    assert_eq!(b.revision(),0);
    let mut tx=b.begin_edit();
    tx.replace_sources(next).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    b.validate().unwrap();
    assert_ne!(b.tessellate(48).unwrap().analytical_volume,
        first.analytical_volume);
    assert!(b.undo().unwrap());
    close(b.tessellate(48).unwrap().faceted_volume,first.faceted_volume,1e-8);
    assert!(b.redo().unwrap());
    b.validate().unwrap();
    assert_eq!(b.tessellate(24).unwrap().genus,1);
}
#[test]
fn rotated_and_translated_triple_mesh_preserves_independent_genus(){
    let source=source();
    let mut rotated=initial();
    rotated.through.center_x=0.0;
    rotated.through.center_z=-1.15;
    rotated.deeper_blind.center_x=0.0;
    rotated.deeper_blind.center_z=0.0;
    rotated.shallower_blind.center_x=0.0;
    rotated.shallower_blind.center_z=1.15;
    let a=ConnectedTripleBrep::from_box(&source,rotated).unwrap();
    let mesh=a.tessellate(32).unwrap();
    assert_eq!(mesh.genus,1);
    close(a.model().mass.volume,body().model().mass.volume,1e-8);
    close(mesh.analytical_volume,a.model().mass.volume,1e-11);
}
