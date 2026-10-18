//! Source-certified G7 geometry.v1 project handoff between independent processes.
use std::io::Write;
use std::process::{Command,Stdio};
use serde_json::{json,Value};

fn call(request:Value)->Value{
    let mut child=Command::new(env!("CARGO_BIN_EXE_gefest-geometry"))
        .stdin(Stdio::piped()).stdout(Stdio::piped())
        .stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap()
        .write_all(request.to_string().as_bytes()).unwrap();
    let output=child.wait_with_output().unwrap();
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    serde_json::from_slice(&output.stdout).unwrap()
}
fn one_wall()->Value{
    json!({
        "operation":"g7_project","feature_key":"real_one_wall",
        "left":{"kind":"block","origin":[-1,1,-3],
            "width":4,"height":3,"depth":6},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":1.5,"height":5},
        "actions":[{"kind":"commit"}],
    })
}
#[test]
fn g7_cli_schema3_handoff_rebuilds_real_brep_and_historical_undo_redo(){
    let created=call(one_wall());
    assert_eq!(created["ok"],true,"{created:?}");
    assert!(created["solid"].is_null());
    let info=&created["data"];
    assert_eq!(info["kind"],"g7_one_wall_boolean_document");
    assert_eq!(info["archive_schema_version"],3);
    assert_eq!(info["document_revision"],1);
    assert_eq!(info["topology"]["shells"],1);
    assert_eq!(info["topology"]["vertices"],12);
    assert_eq!(info["topology"]["edges"],18);
    assert_eq!(info["topology"]["coedges"],36);
    assert_eq!(info["topology"]["faces"],8);
    assert_eq!(info["topology"]["euler_characteristic"],2);
    assert_eq!(info["semantic_names"].as_array().unwrap().len(),84);
    assert_eq!(info["feature_graph"]["output_bindings"].as_array().unwrap().len(),84);
    assert_eq!(info["undo_depth"],1);
    assert!(info["volume"].as_f64().unwrap()>0.);
    let undone=call(json!({
        "operation":"g7_project","archive":info["archive"],
        "actions":[{"kind":"undo"},{"kind":"reopen"}],
    }));
    assert_eq!(undone["ok"],true,"{undone:?}");
    assert_eq!(undone["data"]["has_result"],false);
    assert_eq!(undone["data"]["redo_depth"],1);
    assert_eq!(undone["data"]["semantic_names"].as_array().unwrap().len(),0);
    assert_eq!(undone["data"]["feature_graph"]["output_bindings"].as_array().unwrap().len(),84);
    let redone=call(json!({
        "operation":"g7_project","archive":undone["data"]["archive"],
        "actions":[{"kind":"redo"}],
    }));
    assert_eq!(redone["ok"],true,"{redone:?}");
    assert_eq!(redone["data"]["document_revision"],3);
    assert_eq!(redone["data"]["topology"]["shells"],1);
    assert_eq!(redone["data"]["semantic_names"],info["semantic_names"]);
}
#[test]
fn schema1_schema2_and_mutated_source_volume_or_dag_fail_closed(){
    let original=call(one_wall());
    assert_eq!(original["ok"],true);
    let archive=original["data"]["archive"].clone();
    for version in [1,2]{
        let mut wrong=archive.clone();
        wrong["schema_version"]=json!(version);
        let denied=call(json!({
            "operation":"g7_project","archive":wrong,"actions":[],
        }));
        assert_eq!(denied["ok"],false);
        assert!(denied["data"].is_null());
        assert!(denied["solid"].is_null());
    }
    let mut changed=archive.clone();
    changed["current"]["result"]["source_lower"]=json!(1.);
    let denied=call(json!({
        "operation":"g7_project","archive":changed,"actions":[],
    }));
    assert_eq!(denied["ok"],false);
    assert!(denied["data"].is_null());
    let mut changed=archive.clone();
    changed["graph"]["nodes"][2]["left"]=json!("source_missing");
    let denied=call(json!({
        "operation":"g7_project","archive":changed,"actions":[],
    }));
    assert_eq!(denied["ok"],false);
    assert!(denied["data"].is_null());
}
#[test]
fn g7_rejects_two_wall_source_without_exposing_partial_analytic_shell(){
    let denied=call(json!({
        "operation":"g7_project","feature_key":"wrong_shape",
        "left":{"kind":"block","origin":[-1.2,1,-3],
            "width":2.4,"height":3,"depth":6},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":1.5,"height":5},
        "actions":[{"kind":"commit"}],
    }));
    assert_eq!(denied["ok"],false);
    assert!(denied["solid"].is_null());
    assert!(denied["data"].is_null());
}
#[test]
fn legacy_g5_and_g6_project_schemas_are_not_regressed_by_g7(){
    let old=call(json!({
        "operation":"g5_project","feature_key":"old_four",
        "left":{"kind":"block","origin":[-1.5,1,-1.5],
            "width":3,"height":3,"depth":3},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":2,"height":5},
        "actions":[{"kind":"commit"}],
    }));
    assert_eq!(old["ok"],true,"{old:?}");
    assert_eq!(old["data"]["archive_schema_version"],1);
    assert_eq!(old["data"]["topology"]["shells"],4);
    assert_eq!(old["data"]["semantic_names"].as_array().unwrap().len(),177);
    let two=call(json!({
        "operation":"g6_project","feature_key":"old_two",
        "left":{"kind":"block","origin":[-1.2,1,-3],
            "width":2.4,"height":3,"depth":6},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":1.5,"height":5},
        "actions":[{"kind":"commit"}],
    }));
    assert_eq!(two["ok"],true,"{two:?}");
    assert_eq!(two["data"]["archive_schema_version"],2);
    assert_eq!(two["data"]["topology"]["shells"],2);
    assert_eq!(two["data"]["semantic_names"].as_array().unwrap().len(),115);
}
