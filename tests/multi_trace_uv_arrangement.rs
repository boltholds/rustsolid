//! RustSolid G2: certified selected-pair multi-trace UV regions with true
//! circular source holes. Tests do not imply a global solid Boolean.
use rustsolid::{
    BrepBody, BrepModel, FaceId, FacePairGraph, FacePairId, GeometryTolerance,
    GraphOperand, GraphState, HoleParameters, IntersectionGraphBuilder,
    MultiTraceUvHold, MultiTraceUvPlan, MultiTraceUvResult,
    MultiUvBoundaryProvenance, Point3, Solid, Surface3, ThroughHoleBrep,
    UvCycleKind,
};

fn p(x:f64,y:f64,z:f64)->Point3 {Point3{x,y,z}}
fn tol()->GeometryTolerance {GeometryTolerance::default()}
fn annular_cap(radius:f64)->BrepModel {
    let body=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    BrepBody::from(ThroughHoleBrep::from_box(
        &body,HoleParameters{center_x:0.,center_z:0.,radius},
    ).unwrap()).shared().unwrap()
}
fn cutter_block(x0:f64,x1:f64,z0:f64,z1:f64)->BrepModel {
    BrepBody::from(Solid::block(
        p(x0,4.,z0),x1-x0,2.,z1-z0,
    ).unwrap()).shared().unwrap()
}
fn plane_x(model:&BrepModel,x:f64)->FaceId {
    model.faces.iter().find_map(|face|{
        let Surface3::Plane(plane)=model.face_surface(face.id)? else{return None};
        ((plane.frame.origin.x-x).abs()<1e-8 &&
            plane.frame.normal.x.abs()>0.9).then_some(face.id)
    }).expect("real vertical cutter face")
}
fn selected_pairs(a:&BrepModel,b:&BrepModel,xs:&[f64])->Vec<FacePairGraph>{
    xs.iter().enumerate().map(|(i,&x)|{
        let fb=plane_x(b,x);
        let graph=IntersectionGraphBuilder::new(tol()).unwrap()
            .intersect_faces(a,FaceId(1),b,fb).unwrap();
        assert_eq!(graph.state,GraphState::Ready,
            "expected supported transverse intersection for x={x}: {graph:?}");
        FacePairGraph{
            id:FacePairId(i as u32),
            first_face:FaceId(1),second_face:fb,graph,
        }
    }).collect()
}
fn ready(result:MultiTraceUvResult)->MultiTraceUvPlan {
    match result{
        MultiTraceUvResult::Ready(plan)=>plan,
        MultiTraceUvResult::Held(reason)=>panic!("unexpected G2 Held: {reason:?}"),
    }
}
fn approx(a:f64,b:f64){assert!((a-b).abs()<5e-7,"{a} vs {b}");}

#[test]
fn two_parallel_traces_produce_three_analytic_cells_and_one_bisected_hole(){
    let a=annular_cap(1.);
    let b=cutter_block(0.,1.7,-4.,4.);
    let pairs=selected_pairs(&a,&b,&[0.,1.7]);
    assert_eq!(pairs[0].graph.arcs.len(),2);
    assert_eq!(pairs[1].graph.arcs.len(),1);
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap();
    let plan=ready(result);
    assert_eq!(plan.cells.len(),3);
    assert_eq!(plan.cuts.len(),3);
    assert!(!plan.can_apply);
    approx(plan.original_face_area,36.-std::f64::consts::PI);
    approx(plan.cells.iter().map(|cell|cell.area).sum(),
        plan.original_face_area);
    let mut areas=plan.cells.iter().map(|cell|cell.area).collect::<Vec<_>>();
    areas.sort_by(f64::total_cmp);
    let mut expected=vec![18.-std::f64::consts::PI*0.5,
        1.7*6.-std::f64::consts::PI*0.5,(3.-1.7)*6.];
    expected.sort_by(f64::total_cmp);
    for (actual,expected) in areas.into_iter().zip(expected){approx(actual,expected);}
    let circles=plan.cells.iter().flat_map(|cell|&cell.cycles)
        .flat_map(|cycle|&cycle.spans)
        .filter(|span|span.circular.is_some()).count();
    assert!(circles>=2,"hole boundary must remain analytic circular arcs");
    plan.validate_against(&a,&b,&pairs,tol()).unwrap();
    let mut tampered=plan.clone();
    tampered.cells[0].area+=0.1;
    assert!(tampered.validate_against(&a,&b,&pairs,tol()).is_err());
}

#[test]
fn untouched_circular_hole_is_an_inner_cycle_of_middle_cell(){
    let a=annular_cap(0.7);
    let b=cutter_block(-1.5,1.5,-4.,4.);
    let pairs=selected_pairs(&a,&b,&[-1.5,1.5]);
    let plan=ready(MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap());
    assert_eq!(plan.cells.len(),3);
    let cycles=plan.cells.iter().flat_map(|cell|&cell.cycles)
        .filter(|cycle|cycle.kind==UvCycleKind::Inner).collect::<Vec<_>>();
    assert_eq!(cycles.len(),1);
    assert!(cycles[0].spans.iter().all(|span|span.circular.is_some()));
    approx(cycles[0].signed_area,-std::f64::consts::PI*0.7*0.7);
    approx(plan.cells.iter().map(|cell|cell.area).sum(),
        36.-std::f64::consts::PI*0.7*0.7);
}

#[test]
fn incomplete_finite_cutter_face_must_not_extend_its_line(){
    let a=annular_cap(1.);
    // Source x=0 intersects only the lower part of the cap and ends at
    // the opposite cutter face trim; it cannot be promoted to infinite cut.
    let b=cutter_block(0.,1.7,-4.,-0.2);
    let pairs=selected_pairs(&a,&b,&[0.]);
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap();
    assert!(matches!(result,MultiTraceUvResult::Held(
        MultiTraceUvHold::IncompletePhysicalTrace)), "{result:?}");
}

#[test]
fn coincident_selected_traces_never_generate_duplicate_cut_edges(){
    let a=annular_cap(0.7);
    let b=cutter_block(-1.5,1.5,-4.,4.);
    let mut pairs=selected_pairs(&a,&b,&[-1.5]);
    let mut clone=pairs[0].clone();
    clone.id=FacePairId(1);
    pairs.push(clone);
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap();
    assert!(matches!(result,MultiTraceUvResult::Held(_)), "{result:?}");
}

#[test]
fn every_cut_piece_keeps_a_pair_scoped_arc_reference_and_two_uses(){
    let a=annular_cap(1.);
    let b=cutter_block(0.,1.7,-4.,4.);
    let pairs=selected_pairs(&a,&b,&[0.,1.7]);
    let plan=ready(MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap());
    for cut in &plan.cuts {
        assert!(pairs.iter().any(|pair|pair.id==cut.source.pair &&
            pair.graph.arcs.iter().any(|arc|arc.id==cut.source.local)));
        let uses=plan.cells.iter().flat_map(|cell|&cell.cycles)
            .flat_map(|cycle|&cycle.spans)
            .filter(|span|matches!(span.provenance,
                MultiUvBoundaryProvenance::Cut{segment,source}
                if segment==cut.segment && source==cut.source)).count();
        assert_eq!(uses,2);
    }
}

#[test]
fn tangent_hole_is_held_upstream_not_reconstructed_as_boolean(){
    let a=annular_cap(1.);
    let b=cutter_block(1.,1.7,-4.,4.);
    let fb=plane_x(&b,1.);
    let graph=IntersectionGraphBuilder::new(tol()).unwrap()
        .intersect_faces(&a,FaceId(1),&b,fb).unwrap();
    // At EXACT tangency the upstream graph may retain two line branches
    // meeting on the physical circle. G2 must still reject the vertex event.
    assert!(matches!(graph.state,GraphState::Ready|GraphState::Ambiguous(_)),
        "{graph:?}");
    let pairs=[FacePairGraph {
        id:FacePairId(0),first_face:FaceId(1),second_face:fb,graph,
    }];
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap();
    assert!(matches!(result,
        MultiTraceUvResult::Held(MultiTraceUvHold::TangencyOrVertexEvent |
            MultiTraceUvHold::UpstreamAmbiguity)), "{result:?}");
}

#[test]
fn three_real_holes_and_two_traces_preserve_all_analytic_boundaries(){
    use rustsolid::BrepStructuralEditor;
    let source=Solid::block(p(-3.,0.,-3.),6.,5.,6.).unwrap();
    let mut editor=BrepStructuralEditor::new(source).unwrap();
    let mut tx=editor.begin_transaction();
    tx.add_through_hole(HoleParameters {
        center_x:0.,center_z:-1.,radius:0.5,
    }).unwrap();
    tx.add_through_hole(HoleParameters {
        center_x:-1.5,center_z:1.,radius:0.4,
    }).unwrap();
    tx.add_through_hole(HoleParameters {
        center_x:1.5,center_z:1.,radius:0.6,
    }).unwrap();
    tx.commit().unwrap();
    editor.validate().unwrap();

    let a=editor.model().clone();
    let b=cutter_block(0.,1.7,-4.,4.);
    let pairs=selected_pairs(&a,&b,&[0.,1.7]);
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol()).unwrap();
    let plan=ready(result);
    assert_eq!(plan.cells.len(),3);
    assert_eq!(plan.cuts.len(),4);
    assert_eq!(plan.cells.iter().flat_map(|c|&c.cycles)
        .filter(|c|c.kind==UvCycleKind::Inner).count(),1);
    approx(plan.original_face_area,36.0-
        std::f64::consts::PI*(0.25+0.16+0.36));
    approx(plan.cells.iter().map(|c|c.area).sum::<f64>(),
        plan.original_face_area);
    assert!(!plan.can_apply);
    plan.validate_against(&a,&b,&pairs,tol()).unwrap();
    assert_eq!(editor.revision(),1);
    editor.model().validate().unwrap();
}

#[test]
fn fabricated_empty_pair_is_not_a_valid_g2_source_certificate(){
    let a=annular_cap(1.);
    let b=cutter_block(0.,1.7,-4.,4.);
    let mut pairs=selected_pairs(&a,&b,&[0.,1.7]);
    pairs[1].graph.state=GraphState::Empty;
    pairs[1].graph.nodes.clear();
    pairs[1].graph.arcs.clear();
    pairs[1].graph.support_certificate=None;
    let result=MultiTraceUvPlan::build_selected(
        &a,&b,&pairs,GraphOperand::First,FaceId(1),tol());
    assert!(result.is_err(),
        "empty-source forgery must be rejected by source recomputation");
}
