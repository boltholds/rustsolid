//! Independent watertight facet certificates for *one* genus-N source
//! shell containing full-height tunnels, finite-depth blind floors and
//! optional two-arc crossing blind components. Mass remains analytical.
use rustsolid::{
    BlindHoleParameters,HoleParameters,MixedThroughBlindBrep,
    MixedThroughBlindParameters,Mesh,Point3,Solid,
};
use std::collections::BTreeMap;
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn through(x:f64,z:f64,r:f64)->HoleParameters{
    HoleParameters{center_x:x,center_z:z,radius:r}
}
fn blind(x:f64,z:f64,r:f64,depth:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth}
}
fn setup(t:Vec<HoleParameters>,b:Vec<BlindHoleParameters>)->MixedThroughBlindBrep{
    MixedThroughBlindBrep::from_box(&block(),MixedThroughBlindParameters{
        through:t,blind:b}).unwrap()
}
fn simple()->MixedThroughBlindBrep{
    setup(vec![through(-1.7,0.,0.5)],vec![blind(1.5,0.,0.75,2.)])
}
fn approx(a:f64,b:f64,eps:f64){
    assert!((a-b).abs()<=eps,"got {a}, expected {b}, diff {}",(a-b).abs());
}
fn check_mesh(mesh:&Mesh,expected_vertices:usize,expected_faces:usize){
    assert_eq!(mesh.vertices.len(),expected_vertices);
    assert_eq!(mesh.triangles.len(),expected_faces);
    assert_eq!(mesh.triangle_faces.len(),mesh.triangles.len());
    let mut pairs=BTreeMap::<(u32,u32),Vec<bool>>::new();
    for t in &mesh.triangles{
        let [a,b,c]=*t;
        for (u,v) in [(a,b),(b,c),(c,a)]{
            pairs.entry((u.min(v),u.max(v))).or_default().push(u<v);
        }
    }
    assert!(pairs.values().all(|dirs|dirs.len()==2 && dirs[0]!=dirs[1]),
        "every facet edge must have exactly two oppositely directed uses");
}
#[test]
fn one_through_one_blind_has_two_different_cap_loop_counts_and_one_real_floor(){
    let body=simple();
    for segments in [12,48]{
        let display=body.tessellate(segments).unwrap();
        assert_eq!(display.top_face_loop_count,3);
        assert_eq!(display.bottom_face_loop_count,2);
        assert_eq!(display.top_face_triangles,2*segments+6);
        assert_eq!(display.bottom_face_triangles,segments+4);
        assert_eq!(display.floor_regions.len(),1);
        assert_eq!(display.floor_regions[0].floor_triangles,segments-2);
        assert_eq!(display.floor_regions[0].polygonal_area>0.,true);
        check_mesh(&display.mesh,8+4*segments,8*segments+16);
        assert_eq!(body.model().summary().genus,1);
        approx(display.analytic_volume,body.model().mass.volume,1e-9);
        assert!(display.faceting_volume_error>0.0);
    }
}
#[test]
fn two_through_plus_one_crossed_pair_and_blind_single_has_true_genus_two_mesh(){
    let body=setup(
        vec![through(-1.8,0.,0.4),through(-0.4,-2.,0.35)],
        vec![
            blind(1.1,0.,0.5,1.8),blind(2.0,0.,0.5,1.8),
            blind(0.,1.8,0.3,2.6),
        ]
    );
    for segments in [12,48]{
        let display=body.tessellate(segments).unwrap();
        assert_eq!(display.top_face_loop_count,5);
        assert_eq!(display.bottom_face_loop_count,3);
        assert_eq!(display.floor_regions.len(),2);
        assert!(display.floor_regions.iter().any(|region|
            region.loop_vertices==2*segments));
        assert!(display.floor_regions.iter().any(|region|
            region.loop_vertices==segments));
        check_mesh(&display.mesh,8+10*segments,20*segments+20);
        assert_eq!(body.model().summary().genus,2);
        assert_eq!(body.model().summary().euler_characteristic,-2);
    }
}
#[test]
fn facet_volume_converges_to_exact_analytic_mass_without_changing_revision(){
    let body=simple();
    let mass=body.model().mass;
    let a=body.tessellate(12).unwrap();
    let b=body.tessellate(48).unwrap();
    let c=body.tessellate(192).unwrap();
    assert!(a.faceting_volume_error>b.faceting_volume_error);
    assert!(b.faceting_volume_error>c.faceting_volume_error);
    assert!(c.faceting_volume_error<b.faceting_volume_error*0.08);
    assert_eq!(body.model().mass,mass);
    assert_eq!(body.revision(),0);
}
#[test]
fn optional_display_fidelity_never_replaces_exact_surface_area(){
    let body=simple();
    let through_source=0.5;
    let blind_source=0.75;
    let true_area=168.-2.*PI*through_source*through_source+
        2.*PI*through_source*4.+2.*PI*blind_source*2.;
    approx(body.model().mass.surface_area,true_area,1e-9);
    for n in [12,48,192]{
        let display=body.tessellate(n).unwrap();
        assert!(display.faceted_surface_area<true_area);
        approx(body.model().mass.surface_area,true_area,1e-9);
    }
}
#[test]
fn entire_mesh_is_derived_after_edit_undo_redo_from_current_analytic_carriers(){
    let mut body=simple();
    let original=body.tessellate(32).unwrap();
    {
        let mut tx=body.begin_edit();
        tx.replace_through(0,through(-1.8,0.15,0.45)).unwrap();
        tx.replace_blind(0,blind(1.6,-0.2,0.8,2.3)).unwrap();
        tx.commit().unwrap();
    }
    let updated=body.tessellate(32).unwrap();
    check_mesh(&updated.mesh,8+4*32,8*32+16);
    assert_ne!(original.analytic_volume,updated.analytic_volume);
    assert!(body.undo().unwrap());
    let undone=body.tessellate(32).unwrap();
    approx(original.analytic_volume,undone.analytic_volume,1e-9);
    assert!(body.redo().unwrap());
    let redone=body.tessellate(32).unwrap();
    approx(updated.analytic_volume,redone.analytic_volume,1e-9);
}
#[test]
fn source_permutation_gives_identical_facet_integrals(){
    let throughs=vec![through(-1.6,0.,0.4),through(-1.6,1.4,0.4)];
    let blinds=vec![
        blind(1.0,0.,0.55,1.8),blind(1.95,0.,0.55,1.8),
        blind(0.5,1.9,0.3,2.)
    ];
    let mut rt=throughs.clone();rt.reverse();
    let mut rb=blinds.clone();rb.reverse();
    let left=setup(throughs,blinds);
    let right=setup(rt,rb);
    let a=left.tessellate(32).unwrap();
    let b=right.tessellate(32).unwrap();
    approx(a.faceted_volume,b.faceted_volume,1e-8);
    approx(a.faceted_surface_area,b.faceted_surface_area,1e-8);
    approx(a.analytic_volume,b.analytic_volume,1e-9);
}
#[test]
fn invalid_segment_budget_refuses_without_partial_result_or_state_mutation(){
    let body=simple();
    let before=body.model().mass;
    for n in [0,3,11,513,9999]{
        assert!(body.tessellate(n).is_err());
    }
    assert_eq!(body.model().mass,before);
    assert_eq!(body.revision(),0);
}
#[test]
fn translated_input_beyond_origin_uses_box_local_polygon_charts(){
    let shift=1e4;
    let source=Solid::block(
        p(shift-3.,shift-1.,shift-3.),6.,4.,6.).unwrap();
    let body=MixedThroughBlindBrep::from_box(&source,MixedThroughBlindParameters{
        through:vec![through(shift-1.7,shift,0.5)],
        blind:vec![blind(shift+1.5,shift,0.75,2.)]
    }).unwrap();
    let mesh=body.tessellate(48).unwrap();
    check_mesh(&mesh.mesh,8+4*48,8*48+16);
    approx(mesh.analytic_volume,
        144.-PI*0.5*0.5*4.-PI*0.75*0.75*2.,1e-8);
}
#[test]
fn repeated_display_calls_are_deterministic_for_identical_analytic_sources(){
    let body=simple();
    let a=body.tessellate(24).unwrap();
    let b=body.tessellate(24).unwrap();
    assert_eq!(a.mesh.vertices,b.mesh.vertices);
    assert_eq!(a.mesh.triangles,b.mesh.triangles);
    assert_eq!(a.mesh.triangle_faces,b.mesh.triangle_faces);
    assert_eq!(a.faceted_volume,b.faceted_volume);
    assert_eq!(a.top_face_triangles,b.top_face_triangles);
    assert_eq!(a.bottom_face_triangles,b.bottom_face_triangles);
}
