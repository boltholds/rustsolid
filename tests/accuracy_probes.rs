use rustsolid::{generate_boundary_cases, run_boundary_probes, ProbeExpectation, ProbeModel};
use std::collections::BTreeSet;

#[test]
fn boundary_suite_passes_at_multiple_reproducible_seeds() {
    for seed in [0,1,42,2026,123456789u64] {
        let report=run_boundary_probes(seed);
        let failed=report.outcomes.iter().filter(|o|!o.passed)
            .map(|o|format!("{}: {}",o.case.id,o.reason)).collect::<Vec<_>>();
        assert!(report.all_passed(),"seed={seed} failed: {failed:#?}");
        assert!(report.total>=24 && report.total<=40);
        assert_eq!(report.passed,report.total);
        assert_eq!(report.failed,0);
    }
}

#[test]
fn probes_target_boundary_families_instead_of_cartesian_product() {
    let cases=generate_boundary_cases(42);
    let ids: BTreeSet<_>=cases.iter().map(|c|c.id.as_str()).collect();
    assert_eq!(cases.len(),ids.len());
    assert!(cases.iter().any(|c|c.expectation==ProbeExpectation::Accept));
    assert!(cases.iter().any(|c|c.expectation==ProbeExpectation::Reject));
    for category in ["analytic_oracle","near_tolerance","large_coordinates","invalid_trim",
        "orientation_invariance","mesh_refinement","topology_invariant"] {
        assert!(cases.iter().any(|c|c.category==category),"missing boundary family {category}");
    }
}

#[test]
fn seed_and_case_specs_are_serializable_for_reproduction() {
    let first=generate_boundary_cases(13);
    let again=generate_boundary_cases(13);
    assert_eq!(serde_json::to_value(&first).unwrap(),serde_json::to_value(&again).unwrap());
    let changed=generate_boundary_cases(14);
    assert_ne!(serde_json::to_value(&first).unwrap(),serde_json::to_value(&changed).unwrap());
    let encoded=serde_json::to_value(&first).unwrap();
    assert_eq!(encoded[0]["id"],"block/nominal");
    assert_eq!(encoded[0]["input"]["kind"],"block");
    assert!(first.iter().any(|c|matches!(c.input, ProbeModel::Cylinder{..})));
}

#[test]
fn cylinder_faceting_is_reported_separately_from_kernel_accuracy() {
    let report=run_boundary_probes(13);
    let coarse=report.outcomes.iter().find(|o|o.case.id=="cylinder/eight_facets")
        .unwrap().metrics.as_ref().unwrap();
    let fine=report.outcomes.iter().find(|o|o.case.id=="cylinder/fine_facets")
        .unwrap().metrics.as_ref().unwrap();
    assert!(coarse.display_mesh_volume_relative_error>fine.display_mesh_volume_relative_error*30.0);
    assert!(coarse.display_mesh_volume_relative_error>0.09);
    assert!(fine.display_mesh_volume_relative_error<0.001);
    assert!(coarse.volume_relative_error<1e-12);
    assert!(fine.volume_relative_error<1e-12);
    assert!(coarse.faceting_prediction_absolute_error<1e-10);
    assert!(fine.faceting_prediction_absolute_error<1e-10);
}

#[test]
fn analytic_oracles_detect_scale_and_translation_invariance() {
    let report=run_boundary_probes(77);
    let find=|id:&str|report.outcomes.iter().find(|r|r.case.id==id)
        .unwrap().metrics.as_ref().unwrap();
    for id in ["block/nominal","block/scale_two","block/large_offset",
        "profile/clockwise_closed","profile/far_coordinates", "edit/edge_midpoint",
        "edit/face_diagonal","cylinder/far_origin"] {
        let m=find(id);
        assert!(m.volume_relative_error<1e-9,"{id}");
        assert!(m.area_relative_error<1e-9,"{id}");
    }
}
