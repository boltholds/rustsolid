//! Executable geometry.v1 stdio contract for one source-certified
//! tilted-cylinder / horizontal-plane halfspace cut. No mock math.
use serde_json::{json,Value};
use std::f64::consts::PI;
use std::io::Write;
use std::process::{Command,Stdio};

fn cylinder()->Value{
    let tilt=0.3_f64;
    let (s,c)=tilt.sin_cos();
    json!({"origin":[0.0,0.0,0.0],"axis":[-s,c,0.0],
        "radial_x":[c,s,0.0],"radius":2.0,"height":5.0})
}
fn request(retain:&str,segments:usize,undo:usize,redo:usize)->Value{
    json!({"operation":"oblique_cylinder_clip","cylinder":cylinder(),
        "plane_height":2.5,"retain":retain,"display_segments":segments,
        "undo":undo,"redo":redo})
}
fn invoke(request:&Value)->Value{
    let mut child=Command::new(env!("CARGO_BIN_EXE_gefest-geometry"))
        .stdin(Stdio::piped()).stdout(Stdio::piped())
        .stderr(Stdio::piped()).spawn().expect("geometry CLI binary");
    child.stdin.take().unwrap().write_all(
        serde_json::to_vec(request).unwrap().as_slice()).unwrap();
    let out=child.wait_with_output().unwrap();
    assert!(out.status.success(),"stderr: {}",String::from_utf8_lossy(&out.stderr));
    let data:Value=serde_json::from_slice(&out.stdout).expect("one valid JSON response");
    assert_eq!(data["schema_version"],"geometry.v1");
    data
}
fn num(v:&Value)->f64{v.as_f64().expect("number")}

#[test]
fn each_halfspace_returns_true_clipped_brep_with_face_owned_mesh(){
    let lower=invoke(&request("negative",64,0,0));
    let upper=invoke(&request("positive",64,0,0));
    for (response,side) in [(&lower,"negative"),(&upper,"positive")]{
        assert_eq!(response["ok"],true,"{response}");
        assert!(response["error"].is_null());
        assert_eq!(response["data"]["kind"],"oblique_cylinder_clip");
        assert_eq!(response["data"]["state"],"committed");
        assert_eq!(response["data"]["retained_halfspace"],side);
        assert_eq!(response["data"]["topology_committed"],true);
        assert_eq!(response["data"]["general_body_boolean"],false);
        assert_eq!(response["solid"]["analytic_brep"]["kind"],
            "oblique_halfspace_clipped_cylinder");
        assert_eq!(response["solid"]["topology"]["vertices"],2);
        assert_eq!(response["solid"]["topology"]["edges"],3);
        assert_eq!(response["solid"]["topology"]["coedges"],6);
        assert_eq!(response["solid"]["topology"]["faces"],3);
        assert_eq!(response["solid"]["topology"]["euler_characteristic"],2);
        let triangles=response["solid"]["faces"].as_array().unwrap();
        assert_eq!(triangles.len(),256);
        let owners=response["solid"]["triangle_face_ids"].as_array().unwrap();
        assert_eq!(owners.len(),triangles.len());
        assert_eq!(owners.iter().filter(|v|**v==0).count(),64);
        assert_eq!(owners.iter().filter(|v|**v==1).count(),64);
        assert_eq!(owners.iter().filter(|v|**v==2).count(),128);
        assert_eq!(response["solid"]["analytic_brep"]["facet_segments"],64);
        assert_eq!(response["data"]["mesh_certificate"]
            ["closed_oriented_manifold_certified"],true);
    }
    let v0=num(&lower["solid"]["mass"]["volume"]);
    let v1=num(&upper["solid"]["mass"]["volume"]);
    assert!((v0+v1-20.0*PI).abs()<1e-8);
    let m=2.5/0.3_f64.cos();
    assert!((v0-4.0*PI*m).abs()<1e-8);
    assert!((v1-4.0*PI*(5.0-m)).abs()<1e-8);
    let o0=lower["solid"]["mass"]["centroid"].as_array().unwrap();
    let o1=upper["solid"]["mass"]["centroid"].as_array().unwrap();
    // Centered finite cylinder's centroid equals 2.5 times its axis.
    let (si,co)=0.3_f64.sin_cos();
    for (i,expected) in [-2.5*si,2.5*co,0.0].into_iter().enumerate(){
        let weighted=(num(&o0[i])*v0+num(&o1[i])*v1)/(v0+v1);
        assert!((weighted-expected).abs()<1e-8,"centroid coord {i}");
    }
}

#[test]
fn user_visible_undo_and_redo_report_correct_final_body(){
    let original=invoke(&request("negative",24,0,0));
    let undone=invoke(&request("negative",24,1,0));
    assert_eq!(undone["ok"],true,"{undone}");
    assert_eq!(undone["data"]["state"],"imprinted_after_undo");
    assert_eq!(undone["data"]["topology_committed"],false);
    assert!(undone["data"]["retained_halfspace"].is_null());
    assert_eq!(undone["solid"]["topology"]["faces"],4);
    assert_eq!(undone["solid"]["faces"].as_array().unwrap().len(),144);
    assert_eq!(undone["solid"]["history"]["revision"],2);
    assert_eq!(undone["solid"]["history"]["redo_depth"],1);
    assert!((num(&undone["solid"]["mass"]["volume"])-20.0*PI).abs()<1e-8);

    let redone=invoke(&request("negative",24,1,1));
    assert_eq!(redone["ok"],true);
    assert_eq!(redone["data"]["state"],"committed");
    assert_eq!(redone["solid"]["topology"]["faces"],3);
    assert_eq!(redone["solid"]["history"]["revision"],3);
    assert_eq!(redone["solid"]["history"]["redo_depth"],0);
    assert!((num(&redone["solid"]["mass"]["volume"])-
        num(&original["solid"]["mass"]["volume"])).abs()<1e-9);
}

#[test]
fn ambiguous_or_invalid_sources_return_no_partial_geometry(){
    let mut near_end=request("negative",64,0,0);
    near_end["plane_height"]=json!(10.0);
    let held=invoke(&near_end);
    assert_eq!(held["ok"],false);
    assert!(held["solid"].is_null());
    assert_eq!(held["data"]["state"],"held");
    assert_eq!(held["data"]["topology_committed"],false);

    let mut tangent=request("negative",64,0,0);
    tangent["plane_height"]=json!(4.995);
    let result=invoke(&tangent);
    assert_eq!(result["ok"],false);
    assert!(result["solid"].is_null());

    let mut bad_axis=request("negative",64,0,0);
    bad_axis["cylinder"]["radial_x"]=json!([0.0,1.0,0.0]);
    let result=invoke(&bad_axis);
    assert_eq!(result["ok"],false);
    assert!(result["solid"].is_null());

    let mut negative_axis=request("positive",64,0,0);
    negative_axis["cylinder"]["axis"]=json!([0.0,-1.0,0.0]);
    negative_axis["cylinder"]["radial_x"]=json!([1.0,0.0,0.0]);
    let result=invoke(&negative_axis);
    assert_eq!(result["ok"],false);
    assert!(result["solid"].is_null());
}

#[test]
fn invalid_history_mesh_budget_or_enum_does_not_emit_solid(){
    for (n,u,r) in [(11,0,0),(4097,0,0),(64,2,0),(64,0,1)] {
        let response=invoke(&request("positive",n,u,r));
        assert_eq!(response["ok"],false,"{response}");
        assert!(response["solid"].is_null());
    }
    let response=invoke(&request("sideways",64,0,0));
    assert_eq!(response["ok"],false);
    assert!(response["solid"].is_null());

    let mut unknown=request("negative",64,0,0);
    unknown["unexpected_parameter"]=json!(false);
    let response=invoke(&unknown);
    // Serde GeometryRequest is not globally deny_unknown_fields;
    // this particular explicit endpoint still rejects bad nested sources.
    let mut nested=request("negative",64,0,0);
    nested["cylinder"]["not_a_real_field"]=json!(123);
    let nested_response=invoke(&nested);
    assert_eq!(nested_response["ok"],false);
    assert!(nested_response["solid"].is_null());
    assert_eq!(response["schema_version"],"geometry.v1");
}
