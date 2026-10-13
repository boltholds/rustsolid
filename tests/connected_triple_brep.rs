//! First-party source-verified connected triple cylinder B-rep checks.
//! All topology, curves, UV trims, mass and history are independently
//! rebuilt from mathematical source parameters; no vendor kernel oracle.
use std::f64::consts::PI;
use rustsolid::{
    BlindHoleParameters,BrepSource,ConnectedTripleBrep,
    HoleParameters,Point3,Solid,TripleChainParameters,
    FaceId,FaceSense,Surface3,LoopRole,Curve3,
};
fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->Solid{Solid::block(p(-4.,-2.,-3.),8.,5.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn fixture()->TripleChainParameters{
    TripleChainParameters{
        through:t(-1.15,0.,0.8),
        deeper_blind:b(0.,0.,0.8,2.8),
        shallower_blind:b(1.15,0.,0.8,1.2)
    }
}
fn build()->ConnectedTripleBrep{
    ConnectedTripleBrep::from_box(&block(),fixture()).unwrap()
}
fn close(a:f64,e:f64,eps:f64){
    assert!((a-e).abs()<eps,"got {a:.14}, expected {e:.14}, err {:.4e}",(a-e).abs())
}
#[test]
fn true_triple_brep_is_one_connected_closed_genus_one_analytic_shell(){
    let body=build();
    body.validate().unwrap();
    let m=body.model();m.validate().unwrap();
    let s=m.summary();
    assert_eq!((s.vertices,s.edges,s.coedges,s.faces,s.loops,s.shells),
        (22,38,76,18,20,1));
    assert_eq!((s.euler_characteristic,s.genus),(0,1));
    assert_eq!((s.circular_edges,s.cylindrical_faces),(16,10));
    assert!(m.shells[0].closed);
    assert_eq!(m.faces[0].loops.len(),2);
    assert_eq!(m.faces[1].loops.len(),2);
    assert_eq!(m.loops.iter().filter(|x|x.role==LoopRole::Inner).count(),2);
    for f in &m.faces[6..]{
        assert_eq!(f.sense,FaceSense::Reversed);
    }
    for f in &m.faces[6..16]{
        assert!(matches!(m.face_surface(f.id),Some(Surface3::Cylinder(_))));
    }
    for f in &m.faces[16..]{
        assert!(matches!(m.face_surface(f.id),Some(Surface3::Plane(_))));
        let n=m.oriented_face_normal(f.id,0.0,0.0).unwrap();
        close(n.y,1.0,1e-12);
    }
    assert_eq!(body.brep_model().unwrap().summary(),s);
}
#[test]
fn every_new_physical_circle_or_vertical_edge_has_opposed_twins(){
    let body=build();let m=body.model();
    for edge in &m.edges[12..]{
        let [a,b]=edge.coedges;
        let a=&m.coedges[a.0 as usize];
        let b=&m.coedges[b.0 as usize];
        assert_eq!(a.twin,b.id);
        assert_eq!(b.twin,a.id);
        assert_ne!(a.reversed,b.reversed);
        assert_ne!(a.face,b.face);
        assert_eq!(m.next_in_loop(a.prev).unwrap(),a.id);
        assert_eq!(m.next_in_loop(b.prev).unwrap(),b.id);
        assert!(matches!(m.edge_curve(edge.id),
            Some(Curve3::Circle(_)|Curve3::Line(_))));
        assert!(m.coedge_pcurve(a.id).is_some());
        assert!(m.coedge_pcurve(b.id).is_some());
    }
    let floor0=body.stitch().floors[0].y;
    let floor1=body.stitch().floors[1].y;
    close(floor0,0.2,1e-12);
    close(floor1,1.8,1e-12);
    assert!(floor0<floor1);
}
#[test]
fn source_first_moments_match_an_independent_strip_integration_oracle(){
    let source=block();
    let b=build();
    let m=b.model();
    let input=fixture();
    let disks=[
        (input.through.center_x,input.through.center_z,input.through.radius),
        (input.deeper_blind.center_x,input.deeper_blind.center_z,input.deeper_blind.radius),
        (input.shallower_blind.center_x,input.shallower_blind.center_z,input.shallower_blind.radius),
    ];
    let c=source.mass.centroid;
    let mut removed=0.0;
    let mut moment_x=0.0;
    let mut moment_y=0.0;
    let mut moment_z=0.0;
    let n=32768usize;
    for (slab,layer) in b.stitch().slabs.iter().enumerate(){
        let active=&disks[..slab+1];
        let left=active.iter().map(|d|d.0-d.2).fold(f64::INFINITY,f64::min);
        let right=active.iter().map(|d|d.0+d.2).fold(f64::NEG_INFINITY,f64::max);
        let delta=PI/n as f64;
        let (mut a,mut mx,mut mz)=(0.0,0.0,0.0);
        for j in 0..=n{
            let theta=j as f64*delta;
            let x=left+0.5*(right-left)*(1.0-theta.cos());
            let jac=(right-left)*0.5*theta.sin();
            let mut intervals=active.iter().filter_map(|&(cx,cz,r)|{
                let h2=r*r-(x-cx).powi(2);
                if h2>0.0{
                    let h=h2.sqrt();
                    Some((cz-h,cz+h))
                }else{None}
            }).collect::<Vec<_>>();
            intervals.sort_by(|a,b|a.0.total_cmp(&b.0));
            let mut length=0.0;
            let mut local_z=0.0;
            if let Some(&(mut lo,mut hi))=intervals.first(){
                for &(start,end) in intervals.iter().skip(1){
                    if start>hi{
                        length+=hi-lo;
                        local_z+=((hi-c.z).powi(2)-(lo-c.z).powi(2))*0.5;
                        lo=start;hi=end;
                    }else{hi=hi.max(end);}
                }
                length+=hi-lo;
                local_z+=((hi-c.z).powi(2)-(lo-c.z).powi(2))*0.5;
            }
            let w=if j==0||j==n{1.0}else if j%2==0{2.0}else{4.0};
            a+=w*length*jac;
            mx+=w*(x-c.x)*length*jac;
            mz+=w*local_z*jac;
        }
        a*=delta/3.0;mx*=delta/3.0;mz*=delta/3.0;
        close(a,layer.union_area,6e-5);
        let h=layer.y_upper-layer.y_lower;
        removed+=a*h;
        moment_x+=mx*h;
        moment_z+=mz*h;
        moment_y+=a*h*((layer.y_upper+layer.y_lower)*0.5-c.y);
    }
    let v=source.mass.volume-removed;
    close(m.mass.volume,v,4e-4);
    close(m.mass.centroid.x,c.x-moment_x/v,1e-5);
    close(m.mass.centroid.y,c.y-moment_y/v,1e-5);
    close(m.mass.centroid.z,c.z-moment_z/v,1e-5);
}
#[test]
fn analytic_area_and_mass_do_not_require_a_mesh_as_the_truth_oracle(){
    let body=build();
    let p=body.stitch();
    close(body.model().mass.volume,
        p.predicted_remaining_volume.unwrap(),1e-12);
    close(body.model().mass.surface_area,
        p.analytic_surface_area.unwrap(),1e-12);
    let mesh=body.tessellate(32).unwrap();
    assert_eq!(mesh.genus,1);
    assert_eq!(mesh.mesh.triangles.len(),mesh.mesh.triangle_faces.len());
    close(mesh.analytical_volume,body.model().mass.volume,1e-11);
    assert!(mesh.faceted_volume>=mesh.analytical_volume-1e-7);
    assert_eq!(body.display_mesh(32).unwrap().triangles.len(),
        mesh.mesh.triangles.len());
}
#[test]
fn parameter_edits_are_strong_atomic_with_source_replay_and_undo_redo(){
    let mut body=build();let original=body.model().mass;
    let mut changed=fixture();
    changed.through.center_x=-1.2;
    changed.shallower_blind.depth=1.3;
    let mut tx=body.begin_edit();
    tx.replace_sources(changed).unwrap();
    assert_eq!(tx.commit().unwrap(),1);
    body.validate().unwrap();
    assert_ne!(body.model().mass.volume,original.volume);
    assert!(body.undo().unwrap());
    close(body.model().mass.volume,original.volume,1e-12);
    assert!(body.redo().unwrap());
    assert_eq!(body.revision(),3);
    body.validate().unwrap();
    let before=body.model().mass;
    let mut held=changed;
    held.shallower_blind.center_x=-0.4; // third now touches through
    let mut tx=body.begin_edit();
    tx.replace_sources(held).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(body.revision(),3);
    assert_eq!(body.model().mass,before);
    let mut tx=body.begin_edit();
    let mut broken=changed;
    broken.deeper_blind.radius=f64::NAN;
    assert!(tx.replace_sources(broken).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(body.model().mass,before);
}
#[test]
fn near_tangencies_and_nonbox_rejected_without_partial_commits(){
    let source=block();
    let mut p=fixture();
    p.shallower_blind.center_x=1.6;
    assert!(ConnectedTripleBrep::from_box(&source,p).is_err());
    let mut p=fixture();
    p.shallower_blind.center_x=-0.4;
    assert!(ConnectedTripleBrep::from_box(&source,p).is_err());
    let mut p=fixture();
    p.shallower_blind.depth=p.deeper_blind.depth;
    assert!(ConnectedTripleBrep::from_box(&source,p).is_err());
    let poly=Solid::extrude_xz(&[
        rustsolid::Point2{x:-4.,z:-3.},rustsolid::Point2{x:4.,z:-3.},
        rustsolid::Point2{x:3.,z:3.},rustsolid::Point2{x:-4.,z:3.}
    ],5.).unwrap();
    assert!(ConnectedTripleBrep::from_box(&poly,fixture()).is_err());
}
