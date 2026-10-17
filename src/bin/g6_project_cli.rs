//! Stateless native geometry.v1 bridge for the new bounded G6c/G6d family.
//! Preserves the original G5b four-shell schema and G5Project CLI unmodified.
use super::{BrepInput,Response};
use rustsolid::{
    G5TwoWallDocument,G5TwoWallProjectArchive,G5TwoWallPrepareResult,
    G5TwoWallOperand,GeometryError,GeometryTolerance,
};
use serde::Deserialize;
use serde_json::{json,Value};

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(s.into())
}
#[derive(Debug,Deserialize)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub(super) enum G6Action {
    Commit,Undo,Redo,Reopen,
    ReplaceSource{operand:G6Source,source:BrepInput},
}
#[derive(Debug,Deserialize)]
#[serde(rename_all="snake_case")]
pub(super) enum G6Source{First,Second}
impl G6Source{
    fn operand(&self)->G5TwoWallOperand{
        match self {
            Self::First=>G5TwoWallOperand::First,
            Self::Second=>G5TwoWallOperand::Second,
        }
    }
}
pub(super) fn evaluate(
    archive:Option<Value>,feature_key:Option<String>,
    left:Option<BrepInput>,right:Option<BrepInput>,
    tolerance:Option<GeometryTolerance>,actions:Vec<G6Action>,
)->Response{
    let answer=(||->Result<Value,GeometryError>{
        if actions.len()>32{
            return Err(invalid("G6d maximum 32 actions per request"));
        }
        let mut doc=if let Some(archive)=archive {
            if feature_key.is_some()||left.is_some()||right.is_some()
                ||tolerance.is_some(){
                return Err(invalid("G6d cannot mix stored archive and fresh source fields"));
            }
            let encoded=serde_json::to_string(&archive)
                .map_err(|e|invalid(format!("G6d encoded archive: {e}")))?;
            G5TwoWallProjectArchive::from_json(&encoded)?.reopen()?
        }else{
            let key=feature_key
                .ok_or_else(||invalid("G6d new document requires feature_key"))?;
            let (first,_)=left.ok_or_else(||
                invalid("G6d new document requires left source"))?.prepare()?;
            let (second,_)=right.ok_or_else(||
                invalid("G6d new document requires right source"))?.prepare()?;
            G5TwoWallDocument::new(&key,first,second,tolerance.unwrap_or_default())?
        };
        let mut events=Vec::<Value>::new();
        for action in actions {
            match action {
                G6Action::Commit=>{
                    let proposal=match doc.prepare()?{
                        G5TwoWallPrepareResult::Ready(x)=>x,
                        G5TwoWallPrepareResult::Held(why)=>
                            return Err(invalid(format!(
                                "G6d certified material preparation Held: {why:?}"
                            ))),
                    };
                    let report=doc.commit(proposal)?;
                    events.push(json!({
                        "kind":"commit","revision":report.revision_after,
                        "volume":report.volume,"shells":report.shells,
                        "topology_names":report.identities.len(),
                    }));
                },
                G6Action::Undo=>{
                    let history=doc.undo()?.ok_or_else(||
                        invalid("G6d Undo exceeds history"))?;
                    events.push(json!({
                        "kind":"undo","revision":history.revision_after,
                        "has_result":history.has_result,
                    }));
                },
                G6Action::Redo=>{
                    let history=doc.redo()?.ok_or_else(||
                        invalid("G6d Redo exceeds history"))?;
                    events.push(json!({
                        "kind":"redo","revision":history.revision_after,
                        "has_result":history.has_result,
                    }));
                },
                G6Action::ReplaceSource{operand,source}=>{
                    let (new,_)=source.prepare()?;
                    let revision=doc.replace_source(operand.operand(),new)?;
                    events.push(json!({"kind":"replace_source","revision":revision}));
                },
                G6Action::Reopen=>{
                    let archive=G5TwoWallProjectArchive::from_document(&doc)?;
                    doc=G5TwoWallProjectArchive::from_json(&archive.to_json()?)?.reopen()?;
                    events.push(json!({
                        "kind":"reopen","revision":doc.revision(),
                        "topology_names":doc.bindings().len(),
                    }));
                },
            }
        }
        doc.validate()?;
        let archive=G5TwoWallProjectArchive::from_document(&doc)?;
        let model=doc.model();
        let topology=model.map(|m|{
            let s=m.summary();
            json!({
                "vertices":s.vertices,"edges":s.edges,"coedges":s.coedges,
                "loops":s.loops,"faces":s.faces,"shells":s.shells,
                "euler_characteristic":s.euler_characteristic,"genus":s.genus,
            })
        });
        let names=doc.bindings().iter()
            .map(|b|b.name.as_str().to_string()).collect::<Vec<_>>();
        let graph=archive.graph.clone();
        Ok(json!({
            "kind":"g6_two_wall_boolean_document",
            "archive_schema_version":archive.schema_version,
            "feature_key":doc.feature_key(),
            "document_revision":doc.revision(),
            "has_result":model.is_some(),
            "volume":model.map(|m|m.mass.volume),
            "topology":topology,
            "semantic_names":names,
            "undo_depth":doc.undo_depth(),
            "redo_depth":doc.redo_depth(),
            "feature_graph":graph,
            "archive":archive,
            "events":events,
        }))
    })();
    match answer{
        Ok(data)=>Response{
            schema_version:"geometry.v1",ok:true,error:None,
            data:Some(data),solid:None,
        },
        Err(e)=>Response{
            schema_version:"geometry.v1",ok:false,
            error:Some(e.to_string()),data:None,solid:None,
        },
    }
}
