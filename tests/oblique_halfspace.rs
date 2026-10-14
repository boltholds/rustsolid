//! Closed, source-analytic halfspace cut contracts. This is not a general
//! Boolean between two arbitrary closed bodies: the infinite cutting plane
//! and material-side selection are part of this intentionally bounded gate.
use std::f64::consts::PI;
use rustsolid::{
    BrepModel,ClosedCylinderImprint,Curve3,CylindricalBrep,FaceId,Frame3,
    GeometryTolerance,GraphState,ImprintResult,
    IntersectionGraphBuilder,ModelUnits,ObliqueClipOutcome,
    ObliqueCutHold,ObliqueHalfspaceEditor,ObliqueHalfspaceSide,
    Point3,Solid,Surface3,TopologyDeltaPlanner,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn distance(a:Point3,b:Point3)->f64{(a.x-b.x).hypot(a.y-b.y).hypot(a.z-b.z)}
fn close(a:f64,b:f64,eps:f64){assert!((a-b).abs()<eps,"{a} != {b} (eps={eps})");}
fn tol()->GeometryTolerance{GeometryTolerance::default()}

fn imprinted(tilt:f64,plane_height:f64)->ClosedCylinderImprint{
    let (s,c)=tilt.sin_cos();
    let frame=Frame3::from_axes(p(0.,0.,0.),p(c,s,0.),p(0.,0.,-1.),tol()).unwrap();
    let source=CylindricalBrep::new(frame,2.0,5.0,tol()).unwrap();
    let side=BrepModel::from_cylinder(&source).unwrap();
    let plane_body=Solid::block(p(-9.,plane_height,-9.),18.,1.,18.).unwrap();
    let plane=BrepModel::from_polyhedral(&plane_body,ModelUnits::default()).unwrap();
    let face=plane.faces.iter().find_map(|f|{
        let Surface3::Plane(x)=plane.face_surface(f.id)? else{return None};
        ((x.frame.origin.y-plane_height).abs()<1e-9 &&
         x.frame.normal.y.abs()>0.9).then_some(f.id)
    }).unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&plane,face,&side,FaceId(2)).unwrap();
    assert_eq!(graph.state,GraphState::Ready);
    let plan=TopologyDeltaPlanner::new().plan(&graph,&plane,face,&side,FaceId(2)).unwrap();
    let mut imprint=ClosedCylinderImprint::new(source).unwrap();
    assert!(matches!(imprint.imprint(&graph,&plan,&plane,&side).unwrap(),
        ImprintResult::Committed(_)));
    imprint
}

#[test]
fn classifies_entire_imprinted_faces_by_analytic_halfspace_inequalities(){
    let imprint=imprinted(0.3,2.5);
    let mut editor=ObliqueHalfspaceEditor::new(&imprint).unwrap();
    let region=editor.classify().unwrap();
    assert_eq!(region.source_revision,imprint.revision());
    assert_eq!(region.face_sides,[
        ObliqueHalfspaceSide::Negative,
        ObliqueHalfspaceSide::Positive,
        ObliqueHalfspaceSide::Negative,
        ObliqueHalfspaceSide::Positive,
    ]);
    let m=2.5/0.3_f64.cos();
    close(region.axial_mean,m,1e-9);
    close(region.axial_amplitude,2.0*0.3_f64.tan(),1e-9);
    assert_eq!(editor.applied(),None);
    assert_eq!(editor.revision(),0);
    let source_before=format!("{:?}",imprint.model());
    assert!(matches!(editor.clip(ObliqueHalfspaceSide::Negative).unwrap(),
        ObliqueClipOutcome::Committed(_)));
    assert_eq!(source_before,format!("{:?}",imprint.model()),
        "bounded region decisions cannot mutate the original imprint");
}

#[test]
fn lower_and_upper_closed_halfspace_solids_partition_original_volume_and_centroid(){
    for (tilt,height) in [(0.0,2.5),(0.1,2.2),(0.3,2.5),(0.48,2.0)]{
        let imprint=imprinted(tilt,height);
        let mut lower=ObliqueHalfspaceEditor::new(&imprint).unwrap();
        let mut upper=ObliqueHalfspaceEditor::new(&imprint).unwrap();
        let a=lower.clip(ObliqueHalfspaceSide::Negative).unwrap();
        let b=upper.clip(ObliqueHalfspaceSide::Positive).unwrap();
        assert!(matches!(a,ObliqueClipOutcome::Committed(_)));
        assert!(matches!(b,ObliqueClipOutcome::Committed(_)));
        let lo=lower.model();let hi=upper.model();
        lo.validate().unwrap();hi.validate().unwrap();
        for model in [lo,hi]{
            assert_eq!((model.vertices.len(),model.edges.len(),
                model.coedges.len(),model.faces.len(),model.loops.len()),(2,3,6,3,3));
            assert_eq!(model.euler_characteristic(),2);
            assert_eq!(model.shells.len(),1);
            assert!(model.mass.surface_area>0.);
            assert!(model.mass.volume>0.);
            assert!(model.mass.centroid.is_finite());
        }
        let source=&imprint.model().mass;
        let m=height/tilt.cos();
        close(lo.mass.volume,PI*4.0*m,1e-7);
        close(hi.mass.volume,PI*4.0*(5.0-m),1e-7);
        close(lo.mass.volume+hi.mass.volume,source.volume,1e-7);
        let center=p(
            (lo.mass.centroid.x*lo.mass.volume+
             hi.mass.centroid.x*hi.mass.volume)/source.volume,
            (lo.mass.centroid.y*lo.mass.volume+
             hi.mass.centroid.y*hi.mass.volume)/source.volume,
            (lo.mass.centroid.z*lo.mass.volume+
             hi.mass.centroid.z*hi.mass.volume)/source.volume);
        assert!(distance(center,source.centroid)<1e-8,
            "two volume-region moments must reconstruct the source centroid");
        // Both retained bodies share a true cut plane; the new cap
        // is finite elliptical geometry with a matching 3D edge.
        let cut_lo=lo.face_surface(FaceId(1)).unwrap();
        let cut_hi=hi.face_surface(FaceId(0)).unwrap();
        let (Surface3::Plane(a),Surface3::Plane(b))=(cut_lo,cut_hi)
            else{panic!("missing analytic oblique cut plane");};
        close(a.frame.normal.x+b.frame.normal.x,0.0,1e-8);
        close(a.frame.normal.y+b.frame.normal.y,0.0,1e-8);
        close(a.frame.normal.z+b.frame.normal.z,0.0,1e-8);
        assert!(lo.mass.volume<source.volume && hi.mass.volume<source.volume);
    }
}

#[test]
fn clipped_brep_atomic_undo_redo_and_alternate_side_after_undo(){
    let imprint=imprinted(0.3,2.5);
    let mut editor=ObliqueHalfspaceEditor::new(&imprint).unwrap();
    let original=format!("{:?}",editor.model());
    let result=editor.clip(ObliqueHalfspaceSide::Negative).unwrap();
    let ObliqueClipOutcome::Committed(commit)=result else{panic!("{result:?}");};
    assert_eq!((commit.revision_before,commit.revision_after),(0,1));
    let cut=format!("{:?}",editor.model());
    assert_ne!(original,cut);
    assert_eq!(editor.clip(ObliqueHalfspaceSide::Positive).unwrap(),
        ObliqueClipOutcome::Held(ObliqueCutHold::AlreadyClipped));
    assert_eq!(editor.revision(),1);
    assert!(editor.undo().unwrap());
    assert_eq!(editor.revision(),2);
    assert_eq!(format!("{:?}",editor.model()),original);
    assert!(editor.redo().unwrap());
    assert_eq!(editor.revision(),3);
    assert_eq!(format!("{:?}",editor.model()),cut);
    assert!(editor.undo().unwrap());
    assert_eq!(editor.revision(),4);
    assert!(matches!(editor.clip(ObliqueHalfspaceSide::Positive).unwrap(),
        ObliqueClipOutcome::Committed(_)));
    assert_eq!(editor.revision(),5);
    assert_eq!(editor.applied(),Some(ObliqueHalfspaceSide::Positive));
    assert!(!editor.redo().unwrap());
}

#[test]
fn cut_cannot_be_initialized_from_unimprinted_source(){
    let source=CylindricalBrep::upright(p(0.,0.,0.),2.,5.,tol()).unwrap();
    let imprint=ClosedCylinderImprint::new(source).unwrap();
    assert!(ObliqueHalfspaceEditor::new(&imprint).is_err());
}

#[test]
fn source_integral_independently_confirms_clipped_mass_and_first_moments(){
    // Midpoint polar quadrature of a linear-height field over the true disk.
    // This is algorithmically distinct from the closed-form moments in the
    // production clip and independently exercises the center of mass.
    let imprint=imprinted(0.41,2.3);
    let source=imprint.source_cylinder();
    let Curve3::Harmonic(ellipse)=imprint.model().edge_curve(rustsolid::EdgeId(4)).unwrap()
        else{panic!("missing original analytical section");};
    let axis=source.frame.normal;
    let mean=(ellipse.center.x-source.frame.origin.x)*axis.x+
        (ellipse.center.y-source.frame.origin.y)*axis.y+
        (ellipse.center.z-source.frame.origin.z)*axis.z;
    let cosine=ellipse.cosine.x*axis.x+ellipse.cosine.y*axis.y+
        ellipse.cosine.z*axis.z;
    let sine=ellipse.sine.x*axis.x+ellipse.sine.y*axis.y+
        ellipse.sine.z*axis.z;
    let (nr,na)=(128usize,256usize);
    let radial_step=source.radius/nr as f64;
    let angular_step=2.0*PI/na as f64;
    let mut totals=[[0.0f64;4];2];
    for i in 0..nr{
        let radius=(i as f64+0.5)*radial_step;
        for j in 0..na {
            let theta=(j as f64+0.5)*angular_step;
            let (sn,cs)=theta.sin_cos();
            let surface_height=mean+radius/source.radius*(cosine*cs+sine*sn);
            let dx=radius*cs;let dy=radius*sn;
            let measure=radius*radial_step*angular_step;
            for (side,total) in totals.iter_mut().enumerate() {
                let segment_height=if side==0{surface_height}
                    else{source.height-surface_height};
                let z_moment=if side==0{
                    0.5*surface_height*surface_height
                }else{
                    0.5*(source.height*source.height-surface_height*surface_height)
                };
                total[0]+=segment_height*measure;
                total[1]+=dx*segment_height*measure;
                total[2]+=dy*segment_height*measure;
                total[3]+=z_moment*measure;
            }
        }
    }
    for (i,side) in [ObliqueHalfspaceSide::Negative,
                     ObliqueHalfspaceSide::Positive].into_iter().enumerate() {
        let mut editor=ObliqueHalfspaceEditor::new(&imprint).unwrap();
        assert!(matches!(editor.clip(side).unwrap(),ObliqueClipOutcome::Committed(_)));
        let mass=editor.model().mass;
        let terms=totals[i];
        close(mass.volume,terms[0],1e-6);
        let weighted_local=source.frame.origin;
        let expected=p(
            weighted_local.x+source.frame.x.x*terms[1]/terms[0]+
                source.frame.y.x*terms[2]/terms[0]+axis.x*terms[3]/terms[0],
            weighted_local.y+source.frame.x.y*terms[1]/terms[0]+
                source.frame.y.y*terms[2]/terms[0]+axis.y*terms[3]/terms[0],
            weighted_local.z+source.frame.x.z*terms[1]/terms[0]+
                source.frame.y.z*terms[2]/terms[0]+axis.z*terms[3]/terms[0],
        );
        assert!(distance(expected,mass.centroid)<5e-5,
            "independent radial integration disproves analytic centroid");
    }
}
