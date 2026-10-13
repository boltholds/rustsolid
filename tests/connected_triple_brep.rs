//! Independent numeric and topology acceptance for the dedicated v0.35
//! one-through, two-blind connected crossing chain.
//! The numerical oracle integrates source disk strips, not Green arc moments.
use std::collections::BTreeMap;
use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters,ConnectedTripleBrep,TripleCutParameters,
    HoleParameters,LoopRole,FaceId,FaceSense,Point3,Solid,Surface3,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn box_()->Solid{Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{
    center_x:x,center_z:z,radius:r
}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{BlindHoleParameters{
    center_x:x,center_z:z,radius:r,depth:d
}}
fn fixture()->TripleCutParameters{TripleCutParameters{
    through:t(0.,0.,0.75),
    deep_blind:b(0.7,0.,0.75,2.),
    shallow_blind:b(1.4,0.25,0.6,1.),
}}
fn close(a:f64,e:f64,eps:f64){
    assert!((a-e).abs()<eps,
        "expected {e:.13}, got {a:.13}, error {:.5e}",(a-e).abs());
}
fn disks(x:TripleCutParameters)->[(f64,f64,f64);3]{
    [(x.through.center_x,x.through.center_z,x.through.radius),
     (x.deep_blind.center_x,x.deep_blind.center_z,x.deep_blind.radius),
     (x.shallow_blind.center_x,x.shallow_blind.center_z,x.shallow_blind.radius)]
}
/// Numerical union slice z-intervals and z-moment.
fn width_zmoment(x:f64,source:&[(f64,f64,f64)])->(f64,f64){
    let mut spans=source.iter().filter_map(|&(cx,cz,r)|{
        let h2=r*r-(x-cx)*(x-cx);
        if h2<=0.0{None}else{
            let h=h2.sqrt();Some((cz-h,cz+h))
        }
    }).collect::<Vec<_>>();
    spans.sort_by(|a,b|a.0.total_cmp(&b.0));
    if spans.is_empty(){return (0.0,0.0);}
    let mut lo=spans[0].0;
    let mut hi=spans[0].1;
    let(mut area,mut mz)=(0.0,0.0);
    for &(a,b) in spans.iter().skip(1){
        if a>hi {
            area+=hi-lo;
            mz+=0.5*(hi*hi-lo*lo);
            lo=a;hi=b;
        }else{hi=hi.max(b);}
    }
    area+=hi-lo;
    mz+=0.5*(hi*hi-lo*lo);
    (area,mz)
}
/// Cosine-smoothed composite Simpson quadrature of union of intervals.
fn numeric_area_first_moments(disks:&[(f64,f64,f64)])->(f64,f64,f64){
    let left=disks.iter().map(|&(x,_,r)|x-r).fold(f64::INFINITY,f64::min);
    let right=disks.iter().map(|&(x,_,r)|x+r).fold(f64::NEG_INFINITY,f64::max);
    let n=16384usize;
    let h=PI/n as f64;
    let (mut a,mut mx,mut mz)=(0.0,0.0,0.0);
    for i in 0..=n{
        let theta=i as f64*h;
        let x=left+0.5*(right-left)*(1.0-theta.cos());
        let jac=0.5*(right-left)*theta.sin();
        let (width,z_moment)=width_zmoment(x,disks);
        let w=if i==0||i==n{1.0}else if i%2==0{2.0}else{4.0};
        a+=w*width*jac;
        mx+=w*x*width*jac;
        mz+=w*z_moment*jac;
    }
    (a*h/3.0,mx*h/3.0,mz*h/3.0)
}
#[test]
fn connected_chain_builds_real_shared_coedges_and_two_upward_floors(){
    let body=ConnectedTripleBrep::from_box(&box_(),fixture()).unwrap();
    body.validate().unwrap();
    let m=body.model();
    let s=m.summary();
    assert_eq!(s.shells,1);
    assert_eq!(s.euler_characteristic,0);
    assert_eq!(s.genus,1);
    assert_eq!(m.faces[0].loops.len(),2);
    assert_eq!(m.faces[1].loops.len(),2);
    assert_eq!(m.loops.iter().filter(|l|l.role==LoopRole::Inner).count(),2);
    let floors=m.faces.iter().skip(6).filter(|f|
        matches!(m.face_surface(f.id),Some(Surface3::Plane(_))))
        .collect::<Vec<_>>();
    assert_eq!(floors.len(),2);
    for f in floors {
        assert_eq!(f.sense,FaceSense::Reversed);
        assert_eq!(f.loops.len(),1);
        let n=m.oriented_face_normal(f.id,0.0,0.0).unwrap();
        close(n.y,1.0,1e-12);
    }
    assert_eq!(m.coedges.len(),2*m.edges.len());
    assert!(s.faces>11);
    for edge in &m.edges[12..]{
        let [a,b]=edge.coedges.map(|id|&m.coedges[id.0 as usize]);
        assert_eq!(a.twin,b.id);
        assert_eq!(b.twin,a.id);
        assert_ne!(a.reversed,b.reversed);
    }
    assert_eq!(body.plan().slabs.len(),3);
    close(body.model().mass.volume,
        body.plan().predicted_remaining_volume.unwrap(),1e-9);
}

#[test]
fn independent_numeric_strips_prove_triple_analytic_volume_and_centroid(){
    let source=box_();let pars=fixture();
    let body=ConnectedTripleBrep::from_box(&source,pars).unwrap();
    let m=body.model();
    let disks=disks(pars);
    let mut removed=0.0;
    let mut mx=0.0;
    let mut my=0.0;
    let mut mz=0.0;
    let c0=source.mass.centroid;
    for slab in &body.plan().slabs {
        let active=slab.active_disks.iter().map(|&idx|disks[idx])
            .collect::<Vec<_>>();
        let (area,ax,az)=numeric_area_first_moments(&active);
        close(slab.union_area,area,3e-5);
        let h=slab.y_upper-slab.y_lower;
        removed+=area*h;
        mx+=(ax-c0.x*area)*h;
        mz+=(az-c0.z*area)*h;
        my+=area*h*((slab.y_lower+slab.y_upper)*0.5-c0.y);
    }
    let volume=source.mass.volume-removed;
    close(m.mass.volume,volume,5e-5);
    close(m.mass.centroid.x,c0.x-mx/volume,5e-6);
    close(m.mass.centroid.y,c0.y-my/volume,5e-6);
    close(m.mass.centroid.z,c0.z-mz/volume,5e-6);
}
#[test]
fn derived_triangle_mesh_is_oriented_closed_with_analytic_volume_independent(){
    let body=ConnectedTripleBrep::from_box(&box_(),fixture()).unwrap();
    let mut prev=f64::INFINITY;
    for segments in [12,24,48]{
        let display=body.tessellate(segments).unwrap();
        assert_eq!(display.genus,1);
        assert_eq!(display.floor_faces,2);
        assert!(display.mesh.triangles.len()>100);
        assert_eq!(display.mesh.triangles.len(),display.mesh.triangle_faces.len());
        let mut edge_use=BTreeMap::<(u32,u32),Vec<bool>>::new();
        for triangle in &display.mesh.triangles{
            for (a,b) in [(triangle[0],triangle[1]),
                          (triangle[1],triangle[2]),
                          (triangle[2],triangle[0])] {
                edge_use.entry((a.min(b),a.max(b))).or_default().push(a<b);
            }
        }
        assert!(edge_use.values().all(|v|v.len()==2&&v[0]!=v[1]));
        assert_eq!(display.triangle_edges,edge_use.len());
        close(display.facet_volume,display.independently_polygonized_volume,1e-7);
        close(display.analytic_volume,body.model().mass.volume,1e-12);
        assert!(display.faceting_volume_error>=-1e-7);
        assert!(display.faceting_volume_error.abs()<prev);
        prev=display.faceting_volume_error.abs();
    }
    assert!(prev<0.03);
    assert!(body.model().summary().genus==1);
}
#[test]
fn parametric_undo_redo_publishes_entire_new_source_or_nothing(){
    let mut body=ConnectedTripleBrep::from_box(&box_(),fixture()).unwrap();
    let before=body.model().mass;
    let mut tx=body.begin_edit();
    tx.replace_blind(1,b(1.43,0.28,0.58,1.1)).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    body.validate().unwrap();
    assert_ne!(body.model().mass.volume,before.volume);
    assert!(body.undo().unwrap());
    close(body.model().mass.volume,before.volume,1e-10);
    assert!(body.redo().unwrap());
    body.validate().unwrap();
    assert_eq!(body.revision(),3);
    let current=body.model().mass;
    let mut tx=body.begin_edit();
    tx.replace_through(t(0.,0.,0.78)).unwrap();
    tx.replace_blind(1,b(0.,0.,0.6,1.0)).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(body.revision(),3);
    assert_eq!(body.model().mass,current);
    body.validate().unwrap();
}
#[test]
fn unsupported_triple_junctions_tangencies_and_reversed_depths_fail_closed(){
    let box_=box_();
    let original=fixture();
    for p in [
        TripleCutParameters{
            shallow_blind:b(1.4,0.25,0.6,3.),..original
        },
        TripleCutParameters{
            shallow_blind:b(0.2,0.,0.6,1.),..original
        },
        TripleCutParameters{
            shallow_blind:b(2.05,0.,0.6,1.),..original
        },
        TripleCutParameters{
            deep_blind:b(0.7,0.,0.75,4.0),..original
        },
    ]{
        assert!(ConnectedTripleBrep::from_box(&box_,p).is_err());
    }
}
#[test]
fn failed_display_budget_does_not_edit_source_or_history(){
    let body=ConnectedTripleBrep::from_box(&box_(),fixture()).unwrap();
    let initial=(body.revision(),body.model().mass,body.undo_depth());
    for n in [0,4,9,129,10000]{
        assert!(body.tessellate(n).is_err());
        assert_eq!((body.revision(),body.model().mass,body.undo_depth()),initial);
    }
}
