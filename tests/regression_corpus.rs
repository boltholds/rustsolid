use rustsolid::{
    generate_boundary_cases, minimize_failing_probe, run_builtin_regressions,
    run_regression_corpus, ProbeCase, ProbeExpectation, ProbeModel,
    BUILTIN_REGRESSIONS,
};

#[test]
fn independently_specified_regression_fixtures_pass_in_ci() {
    let report=run_builtin_regressions().expect("fixture schema must be valid");
    let failures=report.outcomes.iter().filter(|o|!o.passed).map(|o| {
        format!("{}: {}",o.fixture.id,o.checks.iter().filter(|c|!c.passed)
            .map(|c| format!("{}: {}",c.name,c.detail)).collect::<Vec<_>>().join("; "))
    }).collect::<Vec<_>>();
    assert!(report.all_passed(),"independent regressions: {failures:#?}");
    assert_eq!(report.total,20);
    assert!(report.outcomes.iter().any(|o|o.fixture.id.contains("thin_face_frame")));
    assert!(report.outcomes.iter().any(|o|o.fixture.id.contains("cylinder_periodic_seam")));
    assert!(report.outcomes.iter().any(|o|o.checks.iter().any(|c|c.name=="uniform_scale")));
    assert!(report.outcomes.iter().any(|o|o.checks.iter().any(|c|c.name=="facet_refinement")));
}

#[test]
fn fixture_results_are_deterministic_and_serializable() {
    let first=run_builtin_regressions().unwrap();
    let second=run_builtin_regressions().unwrap();
    let a=serde_json::to_string(&first).unwrap();
    let b=serde_json::to_string(&second).unwrap();
    assert_eq!(a,b);
    assert!(a.contains("known/block_2_3_4"));
    assert!(a.contains("analytic_golden"));
}

#[test]
fn forged_golden_reference_is_detected_without_any_external_kernel() {
    let mut value:serde_json::Value=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    value["fixtures"][0]["golden"]["volume"]=serde_json::json!(25.0);
    let report=run_regression_corpus(&value.to_string()).unwrap();
    assert_eq!(report.failed,1);
    let fail=report.outcomes.iter().find(|o|!o.passed).unwrap();
    assert_eq!(fail.fixture.id,"known/block_2_3_4");
    assert!(fail.checks.iter().any(|c|c.name=="analytic_golden"&&!c.passed));
}

#[test]
fn invalid_fixture_and_duplicate_identifiers_are_rejected() {
    let mut value:serde_json::Value=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    value["schema"]=serde_json::json!("rustsolid.regressions.v999");
    assert!(run_regression_corpus(&value.to_string()).is_err());
    let mut value:serde_json::Value=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    let cloned=value["fixtures"][0].clone();
    value["fixtures"].as_array_mut().unwrap().push(cloned);
    assert!(run_regression_corpus(&value.to_string()).is_err());
    let mut value:serde_json::Value=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    value["fixtures"][0]["relations"][0]["offset"]=serde_json::json!([0.0,"bad",1.0]);
    assert!(run_regression_corpus(&value.to_string()).is_err());
    let mut value:serde_json::Value=serde_json::from_str(BUILTIN_REGRESSIONS).unwrap();
    value["fixtures"][0]["relations"][1]["factor"]=serde_json::json!(0.0);
    assert!(run_regression_corpus(&value.to_string()).is_err());
}

#[test]
fn shrinker_preserves_failure_class_while_reducing_invalid_profile() {
    let bad=ProbeCase {
        id:"counterexample/minimize".into(),category:"invalid_trim",
        expectation:ProbeExpectation::Accept, // deliberately incorrect oracle
        input:ProbeModel::Extrusion {
            profile:vec![[0.0,0.0],[1.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]],
            height:2.0,
        },
    };
    let reduced=minimize_failing_probe(&bad,64).expect("5 vertices should reduce while still failing");
    let ProbeModel::Extrusion{profile,..}=reduced.input else{panic!("lost model kind")};
    assert!(profile.len()<5);
    assert!(profile.len()>=3);
    assert_eq!(reduced.id,bad.id);
}

#[test]
fn exploratory_seed_sweep_is_reproducible_but_not_the_golden_corpus() {
    let seeds=[0_u64,1,7,13,42,99,2026,987654321];
    for seed in seeds {
        let first=generate_boundary_cases(seed);
        let second=generate_boundary_cases(seed);
        assert_eq!(serde_json::to_value(&first).unwrap(),serde_json::to_value(&second).unwrap());
        assert!(first.len()>=24);
    }
}
