//! G5b geometry.v1 API, separate-process archive handoff.
use std::io::Write;
use std::process::{Command,Stdio};
use serde_json::{json,Value};

fn call(request:Value)->Value{
    let mut child=Command::new(env!("CARGO_BIN_EXE_gefest-geometry"))
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
        .spawn().unwrap();
    child.stdin.take().unwrap()
        .write_all(request.to_string().as_bytes()).unwrap();
    let output=child.wait_with_output().unwrap();
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    serde_json::from_slice(&output.stdout).unwrap()
}
fn input()->Value{
    json!({
        "operation":"g5_project",
        "feature_key":"gefest_cell_A",
        "left":{"kind":"block","origin":[-1.5,1.0,-1.5],
            "width":3.0,"height":3.0,"depth":3.0},
        "right":{"kind":"analytic_cylinder","origin":[0.,0.,0.],
            "radius":2.0,"height":5.0},
        "actions":[{"kind":"commit"}]
    })
}
#[test]
fn geom_v1_reopen_in_new_process_keeps_provenance_and_real_undo_redo(){
    let original=call(input());
    assert_eq!(original["ok"],true,"{original:?}");
    assert!(original["solid"].is_null());
    let data=&original["data"];
    assert_eq!(data["archive_schema_version"],1);
    assert_eq!(data["document_revision"],1);
    assert_eq!(data["topology"]["shells"],4);
    assert_eq!(data["topology"]["euler_characteristic"],8);
    assert_eq!(data["semantic_names"].as_array().unwrap().len(),177);
    assert_eq!(data["feature_graph"]["nodes"].as_array().unwrap().len(),3);
    assert_eq!(data["feature_graph"]["output_bindings"].as_array().unwrap().len(),177);
    assert_eq!(data["undo_depth"],1);
    let archive=data["archive"].clone();
    let undone=call(json!({
        "operation":"g5_project","archive":archive,
        "actions":[{"kind":"undo"},{"kind":"reopen"}]
    }));
    assert_eq!(undone["ok"],true,"{undone:?}");
    assert_eq!(undone["data"]["has_result"],false);
    assert_eq!(undone["data"]["redo_depth"],1);
    assert_eq!(undone["data"]["feature_graph"]["output_bindings"]
        .as_array().unwrap().len(),177);
    let redone=call(json!({
        "operation":"g5_project","archive":undone["data"]["archive"],
        "actions":[{"kind":"redo"}]
    }));
    assert_eq!(redone["ok"],true,"{redone:?}");
    assert_eq!(redone["data"]["document_revision"],3);
    assert_eq!(redone["data"]["has_result"],true);
    assert_eq!(redone["data"]["topology"]["shells"],4);
    assert_eq!(redone["data"]["semantic_names"],data["semantic_names"]);
    assert_eq!(redone["data"]["redo_depth"],0);
}
#[test]
fn geom_v1_detects_corrupt_provenance_without_partial_geometry_leak(){
    let original=call(input());
    assert_eq!(original["ok"],true);
    let mut archive=original["data"]["archive"].clone();
    archive["graph"]["nodes"][2]["left"]=json!("source_missing");
    let response=call(json!({
        "operation":"g5_project","archive":archive,
        "actions":[{"kind":"redo"}]
    }));
    assert_eq!(response["ok"],false);
    assert!(response["data"].is_null());
    assert!(response["solid"].is_null());
    assert!(response["error"].as_str().unwrap().contains("G5b"));
}
#[test]
fn geom_v1_rejects_edit_command_failure_and_forged_archive_merge(){
    let original=call(input());
    assert_eq!(original["ok"],true);
    let denied=call(json!({
        "operation":"g5_project",
        "archive":original["data"]["archive"],
        "actions":[{"kind":"commit"}]
    }));
    assert_eq!(denied["ok"],false);
    assert!(denied["data"].is_null());
    let denied=call(json!({
        "operation":"g5_project",
        "archive":original["data"]["archive"],
        "feature_key":"another_key",
        "actions":[]
    }));
    assert_eq!(denied["ok"],false);
    assert!(denied["data"].is_null());
}
