use std::io::Write;
use std::process::{Command,Stdio};
use serde_json::{json,Value};
fn run(request:Value)->Value {
    let mut child=Command::new(env!("CARGO_BIN_EXE_gefest-geometry"))
        .stdin(Stdio::piped()).stdout(Stdio::piped())
        .stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap()
        .write_all(request.to_string().as_bytes()).unwrap();
    let out=child.wait_with_output().unwrap();
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    serde_json::from_slice(&out.stdout).unwrap()
}
fn initial()->Value{
    json!({
        "operation":"g6_project","feature_key":"two_wall_cut",
        "left":{"kind":"block","origin":[-1.2,1,-3],
            "width":2.4,"height":3,"depth":6},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":1.5,"height":5},
        "actions":[{"kind":"commit"}],
    })
}
#[test]
fn native_g6_schema2_opens_in_another_process_with_two_shells_and_history(){
    let original=run(initial());
    assert_eq!(original["ok"],true,"{original:?}");
    assert!(original["solid"].is_null());
    let result=&original["data"];
    assert_eq!(result["kind"],"g6_two_wall_boolean_document");
    assert_eq!(result["archive_schema_version"],2);
    assert_eq!(result["document_revision"],1);
    assert_eq!(result["topology"]["shells"],2);
    assert_eq!(result["topology"]["vertices"],16);
    assert_eq!(result["topology"]["coedges"],48);
    assert_eq!(result["topology"]["euler_characteristic"],4);
    assert_eq!(result["semantic_names"].as_array().unwrap().len(),115);
    assert_eq!(result["feature_graph"]["output_bindings"].as_array().unwrap().len(),115);
    assert_eq!(result["undo_depth"],1);
    let undid=run(json!({
        "operation":"g6_project","archive":result["archive"],
        "actions":[{"kind":"undo"},{"kind":"reopen"}],
    }));
    assert_eq!(undid["ok"],true,"{undid:?}");
    assert_eq!(undid["data"]["has_result"],false);
    assert_eq!(undid["data"]["redo_depth"],1);
    assert_eq!(undid["data"]["feature_graph"]["output_bindings"].as_array().unwrap().len(),115);
    let redid=run(json!({
        "operation":"g6_project","archive":undid["data"]["archive"],
        "actions":[{"kind":"redo"}],
    }));
    assert_eq!(redid["ok"],true,"{redid:?}");
    assert_eq!(redid["data"]["document_revision"],3);
    assert_eq!(redid["data"]["topology"]["shells"],2);
    assert_eq!(redid["data"]["semantic_names"],result["semantic_names"]);
}
#[test]
fn corrupt_schema1_archive_is_never_confused_with_g6_schema2(){
    let first=run(initial());
    let mut wrong=first["data"]["archive"].clone();
    wrong["schema_version"]=json!(1);
    let failure=run(json!({"operation":"g6_project","archive":wrong,"actions":[]}));
    assert_eq!(failure["ok"],false);
    assert!(failure["solid"].is_null());
    assert!(failure["data"].is_null());
    let mut tampered=first["data"]["archive"].clone();
    tampered["current"]["result"]["component_volumes"][0]=json!(999.);
    let failure=run(json!({"operation":"g6_project","archive":tampered,"actions":[]}));
    assert_eq!(failure["ok"],false);
    assert!(failure["data"].is_null());
}
#[test]
fn unsupported_g6_cut_refuses_commit_and_never_leaks_a_partial_shell(){
    let input=json!({
        "operation":"g6_project","feature_key":"four_not_two",
        "left":{"kind":"block","origin":[-1.5,1,-1.5],
            "width":3,"height":3,"depth":3},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":2,"height":5},
        "actions":[{"kind":"commit"}],
    });
    let result=run(input);
    assert_eq!(result["ok"],false);
    assert!(result["data"].is_null());
    assert!(result["solid"].is_null());
}
#[test]
fn unrelated_existing_g5_four_corner_api_retains_its_schema1_contract(){
    let legacy=run(json!({
        "operation":"g5_project","feature_key":"older_family",
        "left":{"kind":"block","origin":[-1.5,1,-1.5],
            "width":3,"height":3,"depth":3},
        "right":{"kind":"analytic_cylinder","origin":[0,0,0],
            "radius":2,"height":5},
        "actions":[{"kind":"commit"}],
    }));
    assert_eq!(legacy["ok"],true,"{legacy:?}");
    assert_eq!(legacy["data"]["archive_schema_version"],1);
    assert_eq!(legacy["data"]["topology"]["shells"],4);
    assert_eq!(legacy["data"]["semantic_names"].as_array().unwrap().len(),177);
}
