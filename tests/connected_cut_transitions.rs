use rustsolid::{
    BlindHoleParameters,ConnectedCutTransitionPlanner,HoleParameters,
    MultiCutPlanState,Point3,Solid,TransitionArcKind,
};
fn block()->Solid{Solid::block(Point3{x:-3.,y:-1.,z:-3.},6.,4.,6.).unwrap()}
fn t(x:f64,z:f64,r:f64)->HoleParameters{HoleParameters{center_x:x,center_z:z,radius:r}}
fn b(x:f64,z:f64,r:f64,d:f64)->BlindHoleParameters{
    BlindHoleParameters{center_x:x,center_z:z,radius:r,depth:d}
}
fn inputs()->(Vec<HoleParameters>,Vec<BlindHoleParameters>){
    (vec![t(0.,0.,0.75)],
     vec![b(0.7,0.,0.75,2.),b(1.4,0.25,0.6,1.)])
}
#[test]
fn connected_three_source_arcs_have_cross_layer_material_boundaries(){
    let (through,blind)=inputs();
    let plan=ConnectedCutTransitionPlanner::plan(&block(),&through,&blind).unwrap();
    assert_eq!(plan.source_plan.state,MultiCutPlanState::Ready);
    assert_eq!(plan.transitions.len(),2);
    assert!(!plan.can_apply);
    assert_eq!(plan.transitions[0].activated_blind,vec![1]);
    assert_eq!(plan.transitions[1].activated_blind,vec![2]);
    assert!(plan.transitions.iter().all(|step|step.floor_material_area>0.));
    assert!(plan.transitions.iter().all(|step|step.arcs.iter().any(|a|
        a.kind==TransitionArcKind::Continuous)));
    assert!(plan.transitions.iter().all(|step|step.arcs.iter().any(|a|
        a.kind==TransitionArcKind::Emerging)));
    for step in &plan.transitions{
        for arc in &step.arcs{
            assert!(arc.theta_end>arc.theta_start);
            let (x,z,r)=if arc.disk_index==0 {
                (through[0].center_x,through[0].center_z,through[0].radius)
            }else{
                let q=blind[arc.disk_index-1];
                (q.center_x,q.center_z,q.radius)
            };
            for v in [arc.start,arc.end]{
                assert!((v.y-step.height_y).abs()<1e-12);
                assert!(((v.x-x).hypot(v.z-z)-r).abs()<1e-9);
            }
        }
    }
    assert_eq!(plan.predicted_volume,
        plan.source_plan.predicted_remaining_volume.unwrap());
    plan.validate_against(&block(),&through,&blind).unwrap();
}
#[test]
fn transitions_reject_tampered_circle_supports_and_depth_witnesses(){
    let (through,blind)=inputs();
    let plan=ConnectedCutTransitionPlanner::plan(&block(),&through,&blind).unwrap();
    let mut modified=plan.clone();
    modified.transitions[1].arcs[0].start.x+=0.05;
    assert!(modified.validate_against(&block(),&through,&blind).is_err());
    let mut modified=plan.clone();
    modified.transitions[1].floor_material_area+=0.01;
    assert!(modified.validate_against(&block(),&through,&blind).is_err());
    let mut modified=plan.clone();
    modified.can_apply=true;
    assert!(modified.validate_against(&block(),&through,&blind).is_err());
    assert!(plan.validate_against(&block(),&through,
        &[blind[0],b(1.4,0.25,0.6,1.3)]).is_err());
}
#[test]
fn tangency_and_triple_junction_refuse_to_generate_executable_geometry(){
    let source=block();
    let bad=ConnectedCutTransitionPlanner::plan(&source,
        &[t(0.,0.,0.7)],
        &[b(1.4,0.,0.7,2.),b(0.3,0.4,0.5,1.)]);
    assert!(bad.is_err());
    let plain=ConnectedCutTransitionPlanner::plan(&source,
        &[t(-1.3,0.,0.5)],
        &[b(-0.8,0.,0.5,2.),b(1.5,0.,0.4,1.)]);
    assert!(plain.is_err(),
        "disconnected three-source clusters cannot satisfy connected triple gate");
}
