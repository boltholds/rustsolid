//! Stateless geometry.v1 bridge for audited G7 one-wall source-owned Boolean.
//! Does not change g5_project schema1 or g6_project schema2.
use super::{BrepInput,Response};
use rustsolid::{
    G7BooleanDocument,G7ProjectArchive,G7PrepareResult,
    G7Operand,GeometryError,GeometryTolerance,
};
use serde::Deserialize;
use serde_json::{json,Value};

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(s.into())
}
#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum G7Action {
    Commit,Undo,Redo,Reopen,
    ReplaceSource{operand:G7Source,source:BrepInput},
}
#[derive(Debug,Deserialize)]
#[serde(rename_all="snake_case")]
pub(super) enum G7Source{First,Second}
impl G7Source{
    fn operand(&self)->G7Operand{
        match self{Self::First=>G7Operand::First,Self::Second=>G7Operand::Second}
    }
}
pub(super) fn evaluate(
    archive:Option<Value>,
    feature_key:Option<String>,
    left:Option<BrepInput>,right:Option<BrepInput>,
    tolerance:Option<GeometryTolerance>,
    actions:Vec<G7Action>,
)->Response{
    let answer=(||->Result<Value,GeometryError>{
        if actions.len()>32{
            return Err(invalid("G7 maximum 32 actions per geometry.v1 request"));
        }
        let mut document=if let Some(value)=archive{
            if feature_key.is_some()||left.is_some()||right.is_some()
                ||tolerance.is_some(){
                return Err(invalid("G7 stored archive cannot be combined with new operands"));
            }
            let json=serde_json::to_string(&value)
                .map_err(|e|invalid(format!("G7 archive serialization failed: {e}")))?;
            G7ProjectArchive::from_json(&json)?.reopen()?
        }else{
            let key=feature_key.ok_or_else(||
                invalid("G7 new document requires feature_key"))?;
            let (first,_)=left.ok_or_else(||
                invalid("G7 new document requires a box"))?.prepare()?;
            let (second,_)=right.ok_or_else(||
                invalid("G7 new document requires a cylinder"))?.prepare()?;
            G7BooleanDocument::new(
                &key,first,second,tolerance.unwrap_or_default())?
        };
        let mut events=Vec::new();
        for action in actions{
            match action{
                G7Action::Commit=>{
                    let candidate=match document.prepare()?{
                        G7PrepareResult::Ready(x)=>x,
                        G7PrepareResult::Held(why)=>
                            return Err(invalid(format!(
                                "G7 original material/volume Held: {why:?}"
                            ))),
                    };
                    let report=document.commit(candidate)?;
                    events.push(json!({
                        "kind":"commit","revision":report.revision_after,
                        "volume":report.result_volume,
                        "shells":report.shell_components,
                        "semantic_entities":report.topology_names.len(),
                    }));
                },
                G7Action::Undo=>{
                    let change=document.undo()?.ok_or_else(||
                        invalid("G7 Undo exceeds available history"))?;
                    events.push(json!({
                        "kind":"undo","revision":change.revision_after,
                        "has_result":change.has_result,
                    }));
                },
                G7Action::Redo=>{
                    let change=document.redo()?.ok_or_else(||
                        invalid("G7 Redo exceeds available history"))?;
                    events.push(json!({
                        "kind":"redo","revision":change.revision_after,
                        "has_result":change.has_result,
                    }));
                },
                G7Action::ReplaceSource{operand,source}=>{
                    let (body,_)=source.prepare()?;
                    let revision=document.replace_source(operand.operand(),body)?;
                    events.push(json!({
                        "kind":"replace_source","revision":revision,
                    }));
                },
                G7Action::Reopen=>{
                    let saved=G7ProjectArchive::from_document(&document)?;
                    document=G7ProjectArchive::from_json(&saved.to_json()?)?
                        .reopen()?;
                    events.push(json!({
                        "kind":"reopen","revision":document.revision(),
                        "identity_count":document.bindings().len(),
                    }));
                },
            }
        }
        document.validate()?;
        let saved=G7ProjectArchive::from_document(&document)?;
        let topology=document.model().map(|model|{
            let s=model.summary();
            json!({
                "vertices":s.vertices,"edges":s.edges,"coedges":s.coedges,
                "loops":s.loops,"faces":s.faces,"shells":s.shells,
                "euler_characteristic":s.euler_characteristic,
                "genus":s.genus,
            })
        });
        let names=document.bindings().iter().map(|b|
            b.name.as_str().to_owned()).collect::<Vec<_>>();
        let graph=saved.graph.clone();
        Ok(json!({
            "kind":"g7_one_wall_boolean_document",
            "archive_schema_version":saved.schema_version,
            "feature_key":document.feature_key(),
            "document_revision":document.revision(),
            "has_result":document.model().is_some(),
            "volume":document.model().map(|x|x.mass.volume),
            "topology":topology,
            "semantic_names":names,
            "undo_depth":document.undo_depth(),
            "redo_depth":document.redo_depth(),
            "feature_graph":graph,
            "archive":saved,
            "events":events,
        }))
    })();
    match answer{
        Ok(data)=>Response{
            schema_version:"geometry.v1",ok:true,
            data:Some(data),solid:None,error:None,
        },
        Err(error)=>Response{
            schema_version:"geometry.v1",ok:false,
            data:None,solid:None,error:Some(error.to_string()),
        },
    }
}
