//! First-party carrier editing invariants for polyhedral, periodic and
//! genus-one analytic B-reps. No faceted curved Boolean oracle.
use rustsolid::{
    BrepBody,BrepCarrierEditor,CarrierHandleError,CarrierId,CarrierOwner,
    Curve2,Curve3,EdgeId,FaceId,CoedgeId,HistoryDirection,
    HoleParameters,Point2Param,Point3,ParameterRange,PlaneSurface,
    Solid,CylindricalBrep,ThroughHoleBrep,GeometryTolerance,Surface3,
};
use std::f64::consts::TAU;

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn block()->BrepBody{
    BrepBody::from(Solid::block(p(0.,0.,0.),2.,3.,4.).unwrap())
}
fn cylinder()->BrepBody{
    BrepBody::from(CylindricalBrep::upright(
        p(0.,0.,0.),2.,5.,GeometryTolerance::default()).unwrap())
}
fn through_hole()->BrepBody {
    let block=Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(&block,
        HoleParameters{center_x:0.,center_z:0.,radius:1.}).unwrap())
}
fn editor(body:&BrepBody)->BrepCarrierEditor{
    BrepCarrierEditor::from_body(body).unwrap()
}
fn edge(e:&BrepCarrierEditor,index:u32)->rustsolid::CarrierHandle{
    e.handle(CarrierOwner::Edge(EdgeId(index))).unwrap()
}
fn trim(e:&BrepCarrierEditor,index:u32)->rustsolid::CarrierHandle{
    e.handle(CarrierOwner::Coedge(CoedgeId(index))).unwrap()
}
fn face(e:&BrepCarrierEditor,index:u32)->rustsolid::CarrierHandle{
    e.handle(CarrierOwner::Face(FaceId(index))).unwrap()
}
fn assert_mass_equal(a:&rustsolid::BrepModel,b:&rustsolid::BrepModel){
    assert_eq!(a.mass,b.mass);
    assert_eq!(a.bbox,b.bbox);
    assert_eq!(a.summary(),b.summary());
}

#[test]
fn compound_planar_edge_trim_and_uv_rebasing_is_one_atomic_step(){
    let mut e=editor(&block());
    let original=e.snapshot();
    let owner_edge=edge(&e,0);
    let owner_trim=trim(&e,0);
    let owner_face=face(&e,0);
    let old_curve=e.model().edge_curve(EdgeId(0)).unwrap();
    let old_trim=e.model().coedge_pcurve(CoedgeId(0)).unwrap();
    let old_surface=e.model().face_surface(FaceId(0)).unwrap();
    assert!(matches!(old_curve,Curve3::Line(_)));
    assert!(matches!(old_trim,Curve2::Line(_)));
    assert!(matches!(old_surface,Surface3::Plane(_)));
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(owner_edge,2.0,0.25).unwrap();
    tx.reparameterize_trim(owner_trim,3.0,0.125).unwrap();
    tx.rebase_face_uv(owner_face,0.5,-0.25).unwrap();
    assert!(tx.staged_len()>=3);
    let commit=tx.commit().unwrap();
    assert_eq!((commit.revision_before,commit.revision_after),(0,1));
    assert!(commit.changed.iter().any(|x|x.handle==owner_edge));
    assert!(commit.changed.iter().any(|x|x.handle==owner_trim));
    assert!(commit.changed.iter().any(|x|x.handle==owner_face));
    e.validate().unwrap();
    assert_mass_equal(&original.model,e.model());
    assert_eq!(e.resolve(owner_edge).unwrap().handle,owner_edge);
    assert_eq!(e.resolve(owner_trim).unwrap().handle,owner_trim);
    assert_eq!(e.resolve(owner_face).unwrap().handle,owner_face);
    assert!(e.check_snapshot(&original).is_err());
    assert_eq!((e.undo_depth(),e.redo_depth()),(1,0));
    let undone=e.undo().unwrap().unwrap();
    assert_eq!(undone.direction,HistoryDirection::Undo);
    assert_eq!((undone.revision_before,undone.revision_after),(1,2));
    assert_eq!(e.resolve(owner_trim).unwrap().handle,owner_trim);
    assert_eq!(e.model().face_surface(FaceId(0)),Some(old_surface));
    assert_eq!(e.model().coedge_pcurve(CoedgeId(0)),Some(old_trim));
    let redone=e.redo().unwrap().unwrap();
    assert_eq!(redone.direction,HistoryDirection::Redo);
    assert_eq!((redone.revision_before,redone.revision_after),(2,3));
    e.validate().unwrap();
    assert_mass_equal(&original.model,e.model());
}

#[test]
fn periodic_cylinder_seam_trim_reparameterizes_without_moving_any_point(){
    let mut e=editor(&cylinder());
    let first=e.snapshot();
    let circle=edge(&e,0);
    let circular_cap_trim=trim(&e,1);
    let seam=trim(&e,4);
    let side=face(&e,2);
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(circle,1.0,0.7).unwrap();
    tx.reparameterize_trim(circular_cap_trim,1.0,TAU).unwrap();
    tx.reparameterize_trim(seam,2.0,0.25).unwrap();
    tx.rebase_face_uv(side,0.0,0.125).unwrap();
    tx.commit().unwrap();
    e.validate().unwrap();
    assert_mass_equal(&first.model,e.model());
    let stats=e.model().summary();
    assert_eq!(stats.cylindrical_faces,1);
    assert_eq!(stats.circular_edges,2);
    assert_eq!(stats.seam_edges,1);
    assert_eq!(stats.genus,0);
    let a=e.model().geometry.coedges[4].pcurve;
    let b=e.model().geometry.coedges[5].pcurve;
    assert_ne!(a,b,"opposite periodic seam coedges retain separate UV trims");
}

#[test]
fn through_hole_genus_one_and_inner_circle_trims_survive_edit(){
    let mut e=editor(&through_hole());
    let first=e.snapshot();
    let bottom_circle=edge(&e,12);
    let bottom_trim=trim(&e,24);
    let cylinder_side=face(&e,6);
    let bottom_face=face(&e,0);
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(bottom_circle,1.0,0.125).unwrap();
    tx.reparameterize_trim(bottom_trim,1.0,-TAU).unwrap();
    tx.rebase_face_uv(bottom_face,0.125,-0.33).unwrap();
    tx.rebase_face_uv(cylinder_side,0.0,0.75).unwrap();
    tx.commit().unwrap();
    e.validate().unwrap();
    assert_mass_equal(&first.model,e.model());
    let stats=e.model().summary();
    assert_eq!((stats.euler_characteristic,stats.genus),(0,1));
    assert_eq!(e.model().loops.iter().filter(|l|
        l.role==rustsolid::LoopRole::Inner).count(),2);
    assert_eq!(e.model().coedges[28].twin,CoedgeId(29));
    assert_eq!(e.resolve(bottom_trim).unwrap().handle,bottom_trim);
}

#[test]
fn invalid_staging_poison_must_abort_even_when_err_is_ignored(){
    let mut e=editor(&block());
    let original=e.snapshot();
    let good=edge(&e,0);
    let wrong=face(&e,0);
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(good,2.0,0.0).unwrap();
    assert!(tx.reparameterize_edge(wrong,1.0,0.0).is_err());
    assert!(tx.reparameterize_edge(good,1.0,0.0).is_err());
    assert!(tx.commit().is_err());
    e.check_snapshot(&original).unwrap();
    assert_eq!(e.revision(),0);
    assert_eq!((e.undo_depth(),e.redo_depth()),(0,0));
    e.validate().unwrap();
    assert_mass_equal(&original.model,e.model());
    let first_trim=trim(&e,0);
    let second_trim=trim(&e,1);
    let mut tx=e.begin_transaction();
    tx.reparameterize_trim(first_trim,1.0,0.0).unwrap();
    // invalid scale while a prior stage is prepared cannot be ignored.
    assert!(tx.reparameterize_trim(second_trim,f64::NAN,0.0).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(e.revision(),0);
}

#[test]
fn foreign_forged_and_stale_editor_snapshots_are_rejected(){
    let mut a=editor(&block());
    let b=editor(&block());
    let foreign=edge(&b,0);
    assert_eq!(a.resolve(foreign),Err(CarrierHandleError::ForeignEditor));
    let mut forged=edge(&a,0);
    forged.generation=forged.generation.wrapping_add(1);
    assert_eq!(a.resolve(forged),Err(CarrierHandleError::StaleCarrier));
    let original=a.snapshot();
    {
        let mut tx=a.begin_transaction();
        assert!(tx.reparameterize_edge(foreign,1.0,0.0).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(a.revision(),0);
    let stable=edge(&a,0);
    let mut tx=a.begin_transaction();
    tx.reparameterize_edge(stable,1.5,0.1).unwrap();
    tx.commit().unwrap();
    assert!(a.check_snapshot(&original).is_err());
    assert_eq!(a.resolve(stable).unwrap().handle,stable);
}

#[test]
fn physically_deforming_edges_or_surfaces_fail_closed(){
    let mut e=editor(&block());
    let original=e.snapshot();
    // Retarget the edge 0 carrier away from its actual endpoints.
    let edge0=edge(&e,0);
    let Curve3::Line(line)=e.model().edge_curve(EdgeId(0)).unwrap()
        else {panic!("expected line")};
    let wrong=Curve3::Line(rustsolid::Line3{
        origin:p(line.origin.x+0.1,line.origin.y,line.origin.z),
        direction:line.direction,
    });
    {
        let old=e.model().geometry.edges[0];
        let mut tx=e.begin_transaction();
        tx.replace_edge(edge0,wrong,old.domain).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(e.revision(),0);
    // Plane moved along normal without matching geometry is never allowed.
    let owner=face(&e,0);
    let old=e.model().face_surface(FaceId(0)).unwrap();
    let Surface3::Plane(p0)=old else {panic!("expected plane")};
    {
        let mut tx=e.begin_transaction();
        tx.replace_surface(owner,Surface3::Plane(PlaneSurface{
            frame:rustsolid::Frame3{origin:p(p0.frame.origin.x+p0.frame.normal.x*0.1,
                p0.frame.origin.y+p0.frame.normal.y*0.1,
                p0.frame.origin.z+p0.frame.normal.z*0.1),..p0.frame}
        })).unwrap();
        assert!(tx.commit().is_err());
    }
    assert_eq!(e.revision(),0);
    e.check_snapshot(&original).unwrap();
    assert_mass_equal(&original.model,e.model());
}
#[test]
fn wrong_circle_trim_phase_and_periodic_u_chart_changes_poison(){
    let mut e=editor(&cylinder());
    let cap=trim(&e,1);
    let cylinder_face=face(&e,2);
    {
        let mut tx=e.begin_transaction();
        assert!(tx.reparameterize_trim(cap,1.0,TAU/2.0).is_err());
        assert!(tx.commit().is_err());
    }
    {
        let mut tx=e.begin_transaction();
        assert!(tx.rebase_face_uv(cylinder_face,0.2,0.0).is_err());
        assert!(tx.commit().is_err());
    }
    assert_eq!(e.revision(),0);
    e.validate().unwrap();
}
#[test]
fn no_orphan_carriers_after_compaction_and_handles_survive_dense_reindex(){
    let mut model=block().shared().unwrap();
    let original=model.edge_curve(EdgeId(0)).unwrap();
    // Give one edge a noncanonical dense curve position with an unused entry,
    // while keeping all of its physical carrier geometry valid.
    model.geometry.curves3.insert(0,original);
    for edge in &mut model.geometry.edges {edge.curve.0+=1;}
    model.validate().unwrap();
    let mut e=BrepCarrierEditor::new(model).unwrap();
    let handle=edge(&e,0);
    let before=e.resolve(handle).unwrap();
    assert_eq!(before.id,CarrierId::Curve3(rustsolid::Curve3Id(1)));
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(handle,2.0,0.1).unwrap();
    tx.commit().unwrap();
    let after=e.resolve(handle).unwrap();
    assert_eq!(after.id,CarrierId::Curve3(rustsolid::Curve3Id(0)));
    assert_eq!(after.handle,before.handle);
    assert_eq!(e.model().geometry.curves3.len(),e.model().edges.len());
    e.validate().unwrap();
}
#[test]
fn committed_edit_after_undo_discards_redo_branch(){
    let mut e=editor(&block());
    let h=edge(&e,0);
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(h,2.0,0.1).unwrap();
    tx.commit().unwrap();
    e.undo().unwrap().unwrap();
    assert_eq!(e.redo_depth(),1);
    let mut tx=e.begin_transaction();
    tx.reparameterize_edge(h,0.5,0.25).unwrap();
    tx.commit().unwrap();
    assert_eq!(e.redo_depth(),0);
    assert!(e.redo().unwrap().is_none());
    e.validate().unwrap();
}

#[test]
fn direct_curve_replacement_and_trim_rebinding_are_source_equivalent(){
    let mut e=editor(&block());
    let original=e.snapshot();
    let h=edge(&e,0);
    let th=trim(&e,0);
    let Curve3::Line(line)=e.model().edge_curve(EdgeId(0)).unwrap()
        else {panic!("block edge is linear")};
    let Curve2::Line(pcurve)=e.model().coedge_pcurve(CoedgeId(0)).unwrap()
        else {panic!("block coedge is linear")};
    let old_edge=e.model().geometry.edges[0].domain;
    let old_trim=e.model().geometry.coedges[0].domain;
    let ParameterRange::Bounded{start:a,end:b}=old_edge
        else {panic!("edge must have bounded parameters")};
    let ParameterRange::Bounded{start:c,end:d}=old_trim
        else {panic!("trim must have bounded parameters")};
    let edge_offset=0.3;
    let edge_scale=1.75;
    let trim_offset=0.125;
    let trim_scale=1.25;
    let shifted=Curve3::Line(rustsolid::Line3{
        origin:p(line.origin.x+line.direction.x*edge_offset,
            line.origin.y+line.direction.y*edge_offset,
            line.origin.z+line.direction.z*edge_offset),
        direction:p(line.direction.x*edge_scale,
            line.direction.y*edge_scale,line.direction.z*edge_scale),
    });
    let shifted_uv=Curve2::Line(rustsolid::Line2{
        origin:Point2Param{u:pcurve.origin.u+pcurve.delta.u*trim_offset,
            v:pcurve.origin.v+pcurve.delta.v*trim_offset},
        delta:Point2Param{u:pcurve.delta.u*trim_scale,
            v:pcurve.delta.v*trim_scale},
    });
    let mut tx=e.begin_transaction();
    tx.replace_edge(h,shifted,
        ParameterRange::bounded((a-edge_offset)/edge_scale,
            (b-edge_offset)/edge_scale).unwrap()).unwrap();
    tx.replace_trim(th,shifted_uv,
        ParameterRange::bounded((c-trim_offset)/trim_scale,
            (d-trim_offset)/trim_scale).unwrap()).unwrap();
    tx.commit().unwrap();
    e.validate().unwrap();
    assert_mass_equal(&original.model,e.model());
}

#[test]
fn shared_periodic_pcurve_is_copy_on_write_per_coedge_owner(){
    let mut model=cylinder().shared().unwrap();
    // Both circular cap UV curves have identical geometry, but the two
    // coedges must preserve distinct ownership if one is reparameterized.
    let old_a=model.geometry.coedges[1].pcurve;
    let old_b=model.geometry.coedges[3].pcurve;
    assert_eq!(model.geometry.curves2[old_a.0 as usize],
        model.geometry.curves2[old_b.0 as usize]);
    model.geometry.coedges[3].pcurve=old_a;
    model.validate().unwrap();
    let mut e=BrepCarrierEditor::new(model).unwrap();
    let a=trim(&e,1);
    let b=trim(&e,3);
    assert_eq!(e.resolve(a).unwrap().id,e.resolve(b).unwrap().id);
    let before_b=e.model().coedge_pcurve(CoedgeId(3));
    let before_domain=e.model().geometry.coedges[3].domain;
    let mut tx=e.begin_transaction();
    tx.reparameterize_trim(a,1.0,TAU).unwrap();
    tx.commit().unwrap();
    e.validate().unwrap();
    assert_ne!(e.resolve(a).unwrap().id,e.resolve(b).unwrap().id);
    assert_eq!(e.model().coedge_pcurve(CoedgeId(3)),before_b);
    assert_eq!(e.model().geometry.coedges[3].domain,before_domain);
    assert_eq!(e.resolve(a).unwrap().handle,a);
    assert_eq!(e.resolve(b).unwrap().handle,b);
}

#[test]
fn editing_the_same_owner_repeatedly_still_respects_transaction_budget(){
    let mut e=editor(&block());
    let snap=e.snapshot();
    let handle=edge(&e,0);
    let mut tx=e.begin_transaction();
    for _ in 0..256 {
        tx.reparameterize_edge(handle,1.0,0.0).unwrap();
    }
    assert_eq!(tx.staged_len(),256);
    assert!(tx.reparameterize_edge(handle,1.0,0.0).is_err());
    assert!(tx.commit().is_err());
    assert_eq!(e.revision(),0);
    e.check_snapshot(&snap).unwrap();
    assert_mass_equal(&snap.model,e.model());
}
#[test]
fn arena_can_export_an_independent_canonical_editor_without_mutation(){
    let solid=Solid::block(p(-3.,-1.,-3.),6.,4.,6.).unwrap();
    let original=ThroughHoleBrep::from_box(&solid,HoleParameters{
        center_x:0.,center_z:0.,radius:1.}).unwrap();
    let arena=rustsolid::BrepArena::from_through_hole(original).unwrap();
    let frozen=arena.snapshot().unwrap();
    let mut detached=arena.detached_carrier_editor().unwrap();
    let h=edge(&detached,12);
    let mut tx=detached.begin_transaction();
    tx.reparameterize_edge(h,1.0,0.375).unwrap();
    tx.commit().unwrap();
    detached.validate().unwrap();
    assert_eq!(detached.model().summary().genus,1);
    assert_eq!(detached.revision(),1);
    assert_eq!(arena.revision(),0);
    arena.check_snapshot(&frozen).unwrap();
    assert_eq!(arena.model().unwrap().geometry.curves3,
        frozen.model.geometry.curves3);
}

#[test]
fn nonorthonormal_circle_frame_is_rejected_even_when_all_points_match(){
    let mut e=editor(&cylinder());
    let snap=e.snapshot();
    let handle=edge(&e,0);
    let Curve3::Circle(circle)=e.model().edge_curve(EdgeId(0)).unwrap()
        else {panic!("expected circular carrier")};
    // A dishonest parameterization scales the frame basis by 2 and divides
    // the radius by 2: every sampled point on the curve is still identical.
    // The underlying Frame3 must nevertheless remain orthonormal.
    let invalid=Curve3::Circle(rustsolid::Circle3{
        frame:rustsolid::Frame3{
            x:p(circle.frame.x.x*2.0,circle.frame.x.y*2.0,
                circle.frame.x.z*2.0),
            y:p(circle.frame.y.x*2.0,circle.frame.y.y*2.0,
                circle.frame.y.z*2.0),
            ..circle.frame
        },
        radius:circle.radius*0.5,
    });
    let old_domain=e.model().geometry.edges[0].domain;
    let mut tx=e.begin_transaction();
    tx.replace_edge(handle,invalid,old_domain).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(e.revision(),0);
    e.check_snapshot(&snap).unwrap();
    e.validate().unwrap();
}

#[test]
fn periodic_trim_outside_its_uv_chart_is_not_a_valid_reparameterization(){
    let mut e=editor(&cylinder());
    let baseline=e.snapshot();
    let source=e.model().coedge_pcurve(CoedgeId(0)).unwrap();
    let Curve2::Line(line)=source else {panic!("periodic cylindrical side has line pcurve")};
    // u+2*pi represents the same 3D cylinder point but crosses the
    // declared [0,2*pi] UV chart. Silent seam relabeling is forbidden.
    let shifted=Curve2::Line(rustsolid::Line2{
        origin:Point2Param{u:line.origin.u+TAU,v:line.origin.v},
        delta:line.delta,
    });
    let domain=e.model().geometry.coedges[0].domain;
    let handle=trim(&e,0);
    let mut tx=e.begin_transaction();
    tx.replace_trim(handle,shifted,domain).unwrap();
    assert!(tx.commit().is_err());
    assert_eq!(e.revision(),0);
    e.check_snapshot(&baseline).unwrap();
    e.validate().unwrap();
}
#[test]
fn committed_cylinder_uv_rebase_keeps_trimmed_intersection_graph_usable(){
    let other=Solid::block(p(-3.,2.,-3.),6.,1.,6.).unwrap();
    let other_model=other.geometry_store(rustsolid::ModelUnits::default()).unwrap();
    let _=other_model; // confirm source geometry adapters remain compatible
    let first=block().shared().unwrap();
    let mut editor=editor(&cylinder());
    let cylindrical_side=face(&editor,2);
    let mut tx=editor.begin_transaction();
    tx.rebase_face_uv(cylindrical_side,0.0,0.25).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();
    let flat=rustsolid::BrepModel::from_polyhedral(
        &other,rustsolid::ModelUnits::default()).unwrap();
    let graph=rustsolid::IntersectionGraphBuilder::new(
        GeometryTolerance::default()).unwrap()
        .intersect_faces(&flat,FaceId(0),editor.model(),FaceId(2)).unwrap();
    assert_eq!(graph.arcs.len(),1);
    assert!(graph.arcs[0].closed);
    assert_eq!(editor.model().summary().genus,0);
    let _=first;
}
