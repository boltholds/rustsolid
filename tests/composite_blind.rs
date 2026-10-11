use rustsolid::{CompositeBlindBrep,BlindHoleParameters,Point3,Solid,GeometryTolerance,BrepModel,FaceId,Surface3,FaceSense};
use std::f64::consts::PI;
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn h(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn eq(a:f64,b:f64){assert!((a-b).abs()<1e-9,"{a} vs {b}");}
fn check(model:&BrepModel,n:usize){
    model.validate().unwrap();
    let s=model.summary();
    assert_eq!(s.vertices,8+2*n);
    assert_eq!(s.edges,12+3*n);
    assert_eq!(s.coedges,24+6*n);
    assert_eq!(s.faces,6+2*n);
    assert_eq!(s.loops,6+3*n);
    assert_eq!(s.genus,0);
    assert_eq!(s.euler_characteristic,2);
    assert_eq!(s.shells,1);
    assert_eq!(s.circular_edges,2*n);
    assert_eq!(s.cylindrical_faces,n);
    assert_eq!(model.faces[1].loops.len(),1+n);
    assert_eq!(model.faces[0].loops.len(),1);
    for i in 0..n{
        let side=FaceId(6+2*i as u32);
        let floor=FaceId(7+2*i as u32);
        assert_eq!(model.face_sense(side),Some(FaceSense::Reversed));
        assert_eq!(model.face_sense(floor),Some(FaceSense::Reversed));
        assert!(matches!(model.face_surface(side),Some(Surface3::Cylinder(_))));
        assert!(matches!(model.face_surface(floor),Some(Surface3::Plane(_))));
    }
}
#[test]
fn two_blind_bores_have_distinct_real_floors_and_true_analytic_mass(){
    let holes=[h(-1.5,0.,0.6,1.25),h(1.5,0.,0.8,2.4)];
    let m=CompositeBlindBrep::build(&box_(),&holes).unwrap();
    check(&m,2);
    let removed=PI*(0.6*0.6*1.25+0.8*0.8*2.4);
    eq(m.mass.volume,144.-removed);
    eq(m.mass.surface_area,168.+2.*PI*(0.6*1.25+0.8*2.4));
    let cy=holes.iter().map(|hole|{
        let v=PI*hole.radius*hole.radius*hole.depth;
        (1.-(3.-hole.depth/2.))*v
    }).sum::<f64>();
    eq(m.mass.centroid.y,1.+cy/(144.-removed));
}
#[test]
fn three_disjoint_bores_preserve_topology_and_euler_two(){
    let holes=[h(-1.6,-0.5,0.45,1.1),h(1.6,-0.5,0.45,2.),
        h(0.,1.5,0.5,3.2)];
    let m=CompositeBlindBrep::build(&box_(),&holes).unwrap();
    check(&m,3);
    assert_eq!(m.loops.iter().filter(|x|x.role==rustsolid::LoopRole::Inner).count(),3);
}
#[test]
fn insertion_permutation_preserves_analytic_volume_area_centroid(){
    let holes=[h(-1.4,0.,0.6,1.),h(1.4,0.,0.75,2.5)];
    let a=CompositeBlindBrep::build(&box_(),&holes).unwrap();
    let b=CompositeBlindBrep::build(&box_(),&[holes[1],holes[0]]).unwrap();
    check(&a,2);check(&b,2);
    eq(a.mass.volume,b.mass.volume);
    eq(a.mass.surface_area,b.mass.surface_area);
    eq(a.mass.centroid.x,b.mass.centroid.x);
    eq(a.mass.centroid.y,b.mass.centroid.y);
}
#[test]
fn overlapping_and_tangent_circular_openings_refuse_any_result(){
    for x in [0.0,0.5,0.5+1e-11]{
        let holes=[h(-1.,0.,0.75,1.),h(x,0.,0.75,2.)];
        assert!(CompositeBlindBrep::build(&box_(),&holes).is_err());
    }
}
#[test]
fn invalid_individual_cuts_and_empty_set_are_rejected(){
    assert!(CompositeBlindBrep::build(&box_(),&[]).is_err());
    for params in [
        h(0.,0.,1.,0.),h(0.,0.,1.,4.),h(0.,0.,3.,1.),
        h(f64::NAN,0.,1.,1.)
    ]{
        assert!(CompositeBlindBrep::build(&box_(),&[params]).is_err());
    }
}
#[test]
fn composite_must_remain_queryable_by_canonical_brep_query(){
    let m=CompositeBlindBrep::build(&box_(),&[
        h(-1.5,0.,0.6,2.),h(1.5,0.,0.6,2.)
    ]).unwrap();
    let result=rustsolid::geometry_query(&m,
        "{body {kind topology {faceCount genus circularEdges} mass {volume} faces(id:9) {id sense surface {kind}}}}",
        rustsolid::QueryLimits::default()).unwrap();
    assert_eq!(result["body"]["kind"],"analytic_multiple_blind_holes");
    assert_eq!(result["body"]["topology"]["genus"],0);
    assert_eq!(result["body"]["topology"]["faceCount"],10);
    assert_eq!(result["body"]["faces"][0]["sense"],"reversed");
}
#[test]
fn translated_world_keeps_local_mass_moment_accuracy(){
    let s=1e4;
    let source=Solid::block(p(s-3.,s-1.,s-3.),6.,4.,6.).unwrap();
    let holes=[h(s-1.5,s,0.6,1.2),h(s+1.5,s,0.6,2.4)];
    let m=CompositeBlindBrep::build(&source,&holes).unwrap();
    check(&m,2);
    eq(m.mass.volume,144.-PI*0.36*(1.2+2.4));
}
