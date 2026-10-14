//! Independent 2D area and topological-incidence oracles for a real B-rep
//! planar cap with exact circular inner loops and physical crossing arcs.
use std::f64::consts::PI;
use rustsolid::{
    BrepBody, BrepStructuralEditor, FaceId, GeometryTolerance, GraphOperand, GraphState,
    HoleParameters, IntersectionGraph, IntersectionGraphBuilder,
    Point3, Solid, Surface3, ThroughHoleBrep, UvArrangementHold,
    UvBoundaryProvenance, UvCycleKind,
    UvFaceArrangement, UvFaceRegionPlan,
};

fn p(x:f64,y:f64,z:f64)->Point3{Point3{x,y,z}}
fn tol()->GeometryTolerance{GeometryTolerance::default()}
fn block(origin:[f64;3],size:[f64;3])->BrepBody{
    BrepBody::from(Solid::block(
        p(origin[0],origin[1],origin[2]),size[0],size[1],size[2]).unwrap())
}
fn ring(radius:f64,at_x:f64)->BrepBody{
    let cube=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(&cube,
        HoleParameters{center_x:at_x,center_z:0.,radius}).unwrap())
}
fn vertical_face(x:f64)->(BrepBody,FaceId){
    let b=block([x,4.,-4.],[2.,2.,8.]);
    let model=b.shared().unwrap();
    let fid=model.faces.iter().find_map(|face|{
        let Surface3::Plane(s)=model.face_surface(face.id)? else{return None};
        ((s.frame.origin.x-x).abs()<1e-9 &&
         s.frame.normal.x.abs()>0.9).then_some(face.id)
    }).expect("missing finite planar cutter");
    (b,fid)
}
fn arrange(left:&BrepBody,other:&BrepBody,fb:FaceId)
    ->(IntersectionGraph,UvFaceArrangement){
    let a=left.shared().unwrap();
    let b=other.shared().unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&a,FaceId(1),&b,fb).unwrap();
    assert_eq!(graph.state,GraphState::Ready);
    let result=UvFaceRegionPlan::build(
        &a,&b,&graph,GraphOperand::First,tol()).unwrap();
    (graph,result)
}
fn close(x:f64,y:f64){assert!((x-y).abs()<3e-7,"{x} vs {y}");}

#[test]
fn line_bisects_planar_annulus_into_two_complete_connected_uv_cells(){
    let model=ring(1.,0.);
    let (cut,fb)=vertical_face(0.);
    let (graph,result)=arrange(&model,&cut,fb);
    let UvFaceArrangement::Ready(plan)=result else{
        panic!("full annular UV region arrangement failed: {result:?}");
    };
    assert_eq!(graph.arcs.len(),2);
    assert_eq!(plan.cut_twins.len(),2);
    assert!(!plan.can_apply);
    close(plan.original_face_area,36.0-PI);
    let mut areas=plan.cells.iter().map(|x|x.area).collect::<Vec<_>>();
    areas.sort_by(f64::total_cmp);
    close(areas[0],18.0-PI/2.0);
    close(areas[1],18.0-PI/2.0);
    for side in &plan.cells{
        assert_eq!(side.cycles.len(),1);
        assert_eq!(side.cycles[0].kind,UvCycleKind::Outer);
        close(side.cycles[0].signed_area,side.area);
        // Every hole portion retains the native B-rep circle, and
        // both material cells contain one use of each new cut edge.
        let parts=side.cycles.iter().flat_map(|c|c.spans.iter())
            .collect::<Vec<_>>();
        assert_eq!(parts.iter().filter(|p|p.circular.is_some()).count(),1);
        let mut ids=parts.iter().filter_map(|s|match s.provenance{
            UvBoundaryProvenance::Cut{arc}=>Some(arc.0),_=>None
        }).collect::<Vec<_>>();
        ids.sort_unstable();
        assert_eq!(ids.len(),2);
        assert_ne!(ids[0],ids[1]);
        assert!(parts.iter().any(|part|match part.provenance{
            UvBoundaryProvenance::Source{coedge,..}=>
                model.shared().unwrap().coedges.get(coedge.0 as usize).is_some(),
            _=>false
        }));
    }
    for pair in &plan.cut_twins{
        close(pair.positive_start.u,pair.negative_end.u);
        close(pair.positive_start.v,pair.negative_end.v);
        close(pair.positive_end.u,pair.negative_start.u);
        close(pair.positive_end.v,pair.negative_start.v);
    }
    assert_ne!(plan.cells[0].side,plan.cells[1].side);
    plan.validate_against(&model.shared().unwrap(),&cut.shared().unwrap(),
        &graph,tol()).unwrap();
    let mut altered=plan.clone();
    altered.cells[0].area+=0.2;
    assert!(altered.validate_against(&model.shared().unwrap(),
        &cut.shared().unwrap(),&graph,tol()).is_err());
}

#[test]
fn nonintersected_circular_hole_stays_a_real_negative_area_inner_cycle(){
    let left=ring(0.5,1.5);
    let (cut,fb)=vertical_face(0.0);
    let (graph,result)=arrange(&left,&cut,fb);
    let UvFaceArrangement::Ready(plan)=result else{panic!("{result:?}");};
    assert_eq!(graph.arcs.len(),1);
    assert_eq!(plan.cut_twins.len(),1);
    close(plan.original_face_area,36.0-PI*0.25);
    let side_hole=plan.cells.iter().find(|x|
        x.cycles.iter().any(|y|y.kind==UvCycleKind::Inner)).unwrap();
    let clear=plan.cells.iter().find(|x|x.side!=side_hole.side).unwrap();
    assert_eq!(side_hole.cycles.len(),2);
    assert_eq!(side_hole.cycles.iter().filter(|c|c.kind==UvCycleKind::Inner).count(),1);
    let circle=side_hole.cycles.iter().find(|c|c.kind==UvCycleKind::Inner).unwrap();
    close(circle.signed_area,-PI*0.25);
    assert_eq!(circle.spans.len(),1,
        "whole native circle hole must not gain a fake chart seam vertex");
    close(side_hole.area,18.0-PI*0.25);
    close(clear.area,18.0);
}

#[test]
fn rectangle_without_holes_is_a_certified_two_cell_uv_arrangement(){
    let left=block([-3.,0.,-3.],[6.,5.,6.]);
    let (cut,fb)=vertical_face(0.);
    let (graph,result)=arrange(&left,&cut,fb);
    let UvFaceArrangement::Ready(plan)=result else{panic!("{result:?}");};
    assert_eq!(graph.arcs.len(),1);
    close(plan.original_face_area,36.);
    close(plan.cells[0].area,18.);
    close(plan.cells[1].area,18.);
    assert!(plan.cells.iter().all(|x|x.cycles.len()==1));
}

#[test]
fn tangent_hole_causes_no_false_face_cells_or_cut_twins(){
    let left=ring(1.,0.);
    let (cut,fb)=vertical_face(1.0+1e-10);
    let a=left.shared().unwrap();let b=cut.shared().unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&a,FaceId(1),&b,fb).unwrap();
    let result=UvFaceRegionPlan::build(
        &a,&b,&graph,GraphOperand::First,tol()).unwrap();
    assert!(matches!(result,UvFaceArrangement::Held(_)),"{result:?}");
}

#[test]
fn unsupported_nonconvex_outer_profile_remains_held(){
    let points=[
        rustsolid::Point2{x:-3.,z:-3.},rustsolid::Point2{x:3.,z:-3.},
        rustsolid::Point2{x:3.,z:3.},rustsolid::Point2{x:0.5,z:1.},
        rustsolid::Point2{x:-3.,z:3.},
    ];
    let solid=Solid::extrude_xz(&points,5.).unwrap();
    let left=BrepBody::from(solid);
    let (cut,fb)=vertical_face(0.);
    let a=left.shared().unwrap();let b=cut.shared().unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&a,FaceId(1),&b,fb).unwrap();
    let result=UvFaceRegionPlan::build(&a,&b,&graph,
        GraphOperand::First,tol()).unwrap();
    assert!(matches!(result,
        UvFaceArrangement::Held(UvArrangementHold::NonConvexOuterLoop)),
        "unexpected nonconvex arrangement: {result:?}");
}

#[test]
fn three_disjoint_source_holes_with_one_intersected_have_exact_two_cell_areas(){
    let source=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    let mut editor=BrepStructuralEditor::new(source).unwrap();
    let mut tx=editor.begin_transaction();
    tx.add_through_hole(HoleParameters{center_x:0.,center_z:-1.,radius:0.5}).unwrap();
    tx.add_through_hole(HoleParameters{center_x:-1.5,center_z:1.,radius:0.4}).unwrap();
    tx.add_through_hole(HoleParameters{center_x:1.5,center_z:1.,radius:0.6}).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();
    let model=editor.model().clone();
    assert_eq!(model.faces[1].loops.len(),4);
    let (cutter,face)=vertical_face(0.);
    let other=cutter.shared().unwrap();
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&model,FaceId(1),&other,face).unwrap();
    assert_eq!(graph.state,GraphState::Ready);
    assert_eq!(graph.arcs.len(),2);
    let outcome=UvFaceRegionPlan::build(&model,&other,&graph,
        GraphOperand::First,tol()).unwrap();
    let UvFaceArrangement::Ready(plan)=outcome else{panic!("{outcome:?}");};
    assert_eq!(plan.cut_twins.len(),2);
    close(plan.original_face_area,36.0-PI*(0.25+0.16+0.36));
    assert!(plan.cells.iter().all(|c|c.cycles.len()==2));
    let mut areas=plan.cells.iter().map(|c|c.area).collect::<Vec<_>>();
    areas.sort_by(f64::total_cmp);
    close(areas[0],18.0-PI*(0.125+0.36));
    close(areas[1],18.0-PI*(0.125+0.16));
    assert!(plan.cells.iter().all(|cell|
        cell.cycles.iter().filter(|c|c.kind==UvCycleKind::Inner).count()==1));
    plan.validate_against(&model,&other,&graph,tol()).unwrap();
    // The preview may not damage source history/mass or attempt a commit.
    assert!(!plan.can_apply);
    assert_eq!(editor.revision(),1);
    editor.model().validate().unwrap();
}
