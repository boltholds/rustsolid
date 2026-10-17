//! Gefest JSON adapter for the strictly bounded G5b project archive.
//! Stateless by design: each request carries an archive (optional), executes
//! a list of fully validated document actions and returns the new archive.
//! Callers may persist it with the G5ProjectArchive::save_atomic Rust API.
use super::{BrepInput,Response};
use rustsolid::{
    G5BooleanDocument,G5DocumentOperand,G5PrepareResult,G5ProjectArchive,
    GeometryError,GeometryTolerance,
};
use serde::Deserialize;
use serde_json::{json,Value};

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(s.into())
}
#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum G5ProjectAction{
    Commit,
    Undo,
    Redo,
    ReplaceSource {
        operand:G5SourceSelector,source:BrepInput,
    },
    Reopen,
}
#[derive(Debug,Deserialize)]
#[serde(rename_all="snake_case")]
pub(super) enum G5SourceSelector {First,Second}
impl G5SourceSelector{
    fn operand(&self)->G5DocumentOperand {
        match self {Self::First=>G5DocumentOperand::First,
            Self::Second=>G5DocumentOperand::Second}
    }
}
pub(super) fn evaluate(
    archive:Option<Value>,feature_key:Option<String>,
    left:Option<BrepInput>,right:Option<BrepInput>,
    tolerance:Option<GeometryTolerance>,actions:Vec<G5ProjectAction>,
)->Response{
    let answer=(||->Result<Value,GeometryError>{
        if actions.len()>32{
            return Err(invalid("G5b API has a maximum of 32 actions per call"));
        }
        let mut doc=if let Some(value)=archive {
            if left.is_some()||right.is_some()||feature_key.is_some()||tolerance.is_some(){
                return Err(invalid("G5b archive input may not be mixed with fresh source fields"));
            }
            let input=serde_json::to_string(&value)
                .map_err(|e|invalid(format!("G5b encoded archive invalid: {e}")))?;
            G5ProjectArchive::from_json(&input)?.reopen()?
        } else {
            let l=left.ok_or_else(||invalid("G5b new project requires left source"))?;
            let r=right.ok_or_else(||invalid("G5b new project requires right source"))?;
            let feature=feature_key
                .ok_or_else(||invalid("G5b new project requires feature_key"))?;
            let (a,_)=l.prepare()?;
            let (b,_)=r.prepare()?;
            G5BooleanDocument::new(&feature,a,b,tolerance.unwrap_or_default())?
        };
        let mut events=Vec::<Value>::new();
        for action in actions {
            match action {
                G5ProjectAction::Commit=>{
                    let proposal=match doc.prepare()?{
                        G5PrepareResult::Ready(x)=>x,
                        G5PrepareResult::Held(why)=>
                            return Err(invalid(format!(
                                "G5b bounded Boolean preparation Held: {why:?}"))),
                    };
                    let change=doc.commit(proposal)?;
                    events.push(json!({"kind":"commit",
                        "revision":change.revision_after,
                        "shells":change.shell_components,
                        "volume":change.committed_volume,
                        "semantic_entities":change.bound_entities}));
                },
                G5ProjectAction::Undo=>{
                    let changed=doc.undo()?.ok_or_else(||invalid(
                        "G5b Undo exceeds available history"))?;
                    events.push(json!({"kind":"undo",
                        "revision":changed.revision_after,
                        "has_result":changed.has_result}));
                },
                G5ProjectAction::Redo=>{
                    let changed=doc.redo()?.ok_or_else(||invalid(
                        "G5b Redo exceeds available history"))?;
                    events.push(json!({"kind":"redo",
                        "revision":changed.revision_after,
                        "has_result":changed.has_result}));
                },
                G5ProjectAction::ReplaceSource{operand,source}=>{
                    let (body,_)=source.prepare()?;
                    let revision=doc.replace_source(operand.operand(),body)?;
                    events.push(json!({"kind":"replace_source","revision":revision}));
                },
                G5ProjectAction::Reopen=>{
                    // Simulate a process boundary with no shared old handle
                    // registry: serialize, parse, reconstruct and revalidate.
                    let saved=G5ProjectArchive::from_document(&doc)?;
                    doc=G5ProjectArchive::from_json(&saved.to_json()?)?.reopen()?;
                    events.push(json!({"kind":"reopen",
                        "revision":doc.revision(),
                        "identity_count":doc.identity_count()}));
                }
            }
        }
        doc.validate()?;
        let output=G5ProjectArchive::from_document(&doc)?;
        let graph=output.graph.clone();
        let model=doc.model();
        let topology=model.map(|m|{
            let s=m.summary();
            json!({"vertices":s.vertices,"edges":s.edges,
                "coedges":s.coedges,"loops":s.loops,"faces":s.faces,
                "shells":s.shells,"euler_characteristic":s.euler_characteristic,
                "genus":s.genus})
        });
        let volume=model.map(|m|m.mass.volume);
        let names=doc.bindings().iter()
            .map(|b|b.name.as_str().to_string()).collect::<Vec<_>>();
        Ok(json!({
            "kind":"g5_persistent_boolean_document",
            "archive_schema_version":output.schema_version,
            "feature_key":doc.feature_key(),
            "document_revision":doc.revision(),
            "has_result":model.is_some(),
            "volume":volume,"topology":topology,
            "semantic_names":names,
            "undo_depth":doc.undo_depth(),
            "redo_depth":doc.redo_depth(),
            "events":events,
            "feature_graph":graph,
            "archive":output,
        }))
    })();
    match answer {
        Ok(data)=>Response{
            schema_version:"geometry.v1",ok:true,error:None,data:Some(data),solid:None,
        },
        Err(error)=>Response{
            schema_version:"geometry.v1",ok:false,
            error:Some(error.to_string()),data:None,solid:None,
        }
    }
}
