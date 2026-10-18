//! G7d: strict schema-v3 persistent document for a one-wall circular notch.
//!
//! Only source primitive parameters, independent volume evidence, semantic
//! feature provenance and Undo/Redo states are serialized. Neither B-rep
//! geometry nor process-local handles may be accepted as material evidence.
//! Reopening recreates and re-certifies all G7a/G7b historical geometry.

use crate::{
    BrepBody,BooleanOperation,AnalyticOneWallPlan,AnalyticOneWallResult,
    AnalyticOneWallVolumeOracle,AnalyticOneWallVolumeResult,
    G5FeatureGraph,G5FeatureNode,G5GraphBinding,G5ProjectSource,
    G7BooleanDocument,GeometryError,GeometryTolerance,
};
use crate::g7_one_wall_document::{
    G7State,G7Accepted,G7IdentityRegistry,g7_valid_key,g7_restore_token,
};
use serde::{Serialize,Deserialize};
use std::fs::{self,File,OpenOptions};
use std::io::{Read,Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64,Ordering};

const SCHEMA_VERSION:u32=3;
const MAX_ARCHIVE_BYTES:usize=1024*1024;
const MAX_HISTORY_STATES:usize=64;
static NEXT_TEMP:AtomicU64=AtomicU64::new(1);

fn invalid(message:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(message.into())
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G7VolumeStamp {
    pub retained_volume:f64,
    pub source_lower:f64,
    pub source_upper:f64,
    pub flux_volume:f64,
    pub true_angular_span:f64,
    pub integration_panels:usize,
    pub source_numerical_guard:f64,
}
impl G7VolumeStamp{
    fn from_accepted(result:&G7Accepted)->Self{
        Self{
            retained_volume:result.candidate.brep().mass.volume,
            source_lower:result.volume.source.lower,
            source_upper:result.volume.source.upper,
            flux_volume:result.volume.result_flux.volume,
            true_angular_span:result.volume.result_flux.true_angle_span,
            integration_panels:result.volume.source.integration_panels,
            source_numerical_guard:result.volume.source.numerical_guard,
        }
    }
    fn valid(&self)->bool{
        self.retained_volume.is_finite()&&self.retained_volume>0.
            &&self.source_lower.is_finite()&&self.source_lower>0.
            &&self.source_upper.is_finite()&&self.source_upper>self.source_lower
            &&self.flux_volume.is_finite()&&self.flux_volume>0.
            &&self.true_angular_span.is_finite()
            &&self.true_angular_span>std::f64::consts::PI
            &&self.true_angular_span<std::f64::consts::TAU
            &&self.integration_panels>0&&self.integration_panels<=131072
            &&self.source_numerical_guard.is_finite()
            &&self.source_numerical_guard>0.
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G7ArchiveState {
    pub sources:[G5ProjectSource;2],
    pub result:Option<G7VolumeStamp>,
}
impl G7ArchiveState{
    fn from_state(state:&G7State)->Result<Self,GeometryError>{
        // The schema stores parametric primitive definitions, not intrinsic
        // Solid topology revision histories. Reject rather than silently
        // dropping any source-solid revision that cannot be reconstructed.
        if let BrepBody::Polyhedral(ref solid)=state.inputs[0]{
            if solid.revision!=0{
                return Err(invalid("G7d cannot persist nonzero source-solid revisions"));
            }
        }
        Ok(Self{
            sources:[
                G5ProjectSource::from_body(&state.inputs[0])?,
                G5ProjectSource::from_body(&state.inputs[1])?,
            ],
            result:state.accepted.as_ref().map(G7VolumeStamp::from_accepted),
        })
    }
    fn rebuild(&self,tolerance:GeometryTolerance)->Result<G7State,GeometryError>{
        let inputs=[
            self.sources[0].to_body()?,
            self.sources[1].to_body()?,
        ];
        let accepted=if let Some(saved)=&self.result{
            if !saved.valid(){
                return Err(invalid("G7d archive has an invalid numerical volume stamp"));
            }
            let candidate=match AnalyticOneWallPlan::build(
                &inputs[0],&inputs[1],BooleanOperation::Difference,tolerance)?
            {
                AnalyticOneWallResult::Ready(x)=>x,
                AnalyticOneWallResult::Held(_)=>
                    return Err(invalid("G7d archived geometry is no longer source-certified")),
            };
            let volume=match AnalyticOneWallVolumeOracle::build(
                &inputs[0],&inputs[1],&candidate,tolerance)?
            {
                AnalyticOneWallVolumeResult::Ready(x)=>x,
                AnalyticOneWallVolumeResult::Held(_)=>
                    return Err(invalid("G7d archived independent volume is no longer certified")),
            };
            let rebuilt=G7Accepted{candidate,volume};
            if &G7VolumeStamp::from_accepted(&rebuilt)!=saved {
                return Err(invalid("G7d archived volume witnesses were altered"));
            }
            Some(rebuilt)
        }else{None};
        let result=G7State{inputs,accepted};
        result.validate(tolerance)?;
        Ok(result)
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G7ProjectArchive{
    pub schema_version:u32,
    pub feature_key:String,
    pub document_revision:u64,
    pub tolerance:GeometryTolerance,
    pub graph:G5FeatureGraph,
    pub current:G7ArchiveState,
    pub undo:Vec<G7ArchiveState>,
    pub redo:Vec<G7ArchiveState>,
}
impl G7ProjectArchive{
    fn graph_from_document(doc:&G7BooleanDocument)
        ->Result<G5FeatureGraph,GeometryError>{
        let accepted=doc.state.accepted.iter()
            .chain(doc.undo.iter().filter_map(|x|x.accepted.as_ref()))
            .chain(doc.redo.iter().filter_map(|x|x.accepted.as_ref()))
            .next();
        let names=G7IdentityRegistry::names(
            accepted.map(|x|&x.candidate),&doc.feature_key)?;
        if !names.is_empty()&&names.len()!=84{
            return Err(invalid("G7d cannot persist a partial source-owned DAG"));
        }
        let mut bindings=names.into_values().map(|topology_name|G5GraphBinding{
            topology_name,producing_node:"difference".into(),
        }).collect::<Vec<_>>();
        bindings.sort_by(|a,b|a.topology_name.cmp(&b.topology_name));
        Ok(G5FeatureGraph{
            nodes:vec![
                G5FeatureNode::Source{id:"source_a".into(),operand_index:0},
                G5FeatureNode::Source{id:"source_b".into(),operand_index:1},
                G5FeatureNode::Difference{
                    id:"difference".into(),left:"source_a".into(),
                    right:"source_b".into(),
                    output_scope:format!("{}/difference",doc.feature_key),
                },
            ],
            output_bindings:bindings,
        })
    }
    pub fn from_document(doc:&G7BooleanDocument)->Result<Self,GeometryError>{
        doc.validate()?;
        if doc.undo.len()+doc.redo.len()>MAX_HISTORY_STATES{
            return Err(invalid("G7d archive history exceeds 64 states"));
        }
        let result=Self{
            schema_version:SCHEMA_VERSION,
            feature_key:doc.feature_key.clone(),
            document_revision:doc.revision,
            tolerance:doc.tolerance,
            graph:Self::graph_from_document(doc)?,
            current:G7ArchiveState::from_state(&doc.state)?,
            undo:doc.undo.iter().map(G7ArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
            redo:doc.redo.iter().map(G7ArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
        };
        let _=result.to_json()?;
        Ok(result)
    }
    pub fn reopen(&self)->Result<G7BooleanDocument,GeometryError>{
        if self.schema_version!=SCHEMA_VERSION
            ||!g7_valid_key(&self.feature_key)
            ||self.document_revision==u64::MAX
            ||self.document_revision<self.undo.len() as u64
            ||self.undo.len()+self.redo.len()>MAX_HISTORY_STATES{
            return Err(invalid("G7d unknown schema, revision or history policy"));
        }
        self.tolerance.validate()?;
        let current=self.current.rebuild(self.tolerance)?;
        let undo=self.undo.iter()
            .map(|x|x.rebuild(self.tolerance))
            .collect::<Result<Vec<_>,_>>()?;
        let redo=self.redo.iter()
            .map(|x|x.rebuild(self.tolerance))
            .collect::<Result<Vec<_>,_>>()?;
        // Under the narrow G7c lifecycle, two consecutive committed
        // source states require an intervening source replacement.
        let mut chronological=undo.iter().collect::<Vec<_>>();
        chronological.push(&current);
        chronological.extend(redo.iter().rev());
        for state in chronological.windows(2){
            if state[0].accepted.is_some()&&state[1].accepted.is_some(){
                return Err(invalid("G7d history skips an intervening source edit"));
            }
        }
        let token=g7_restore_token()?;
        let registry=G7IdentityRegistry::from_candidate(
            current.accepted.as_ref().map(|x|&x.candidate),
            &self.feature_key,token)?;
        let doc=G7BooleanDocument {
            token,feature_key:self.feature_key.clone(),
            revision:self.document_revision,tolerance:self.tolerance,
            state:current,registry,undo,redo,
        };
        doc.validate()?;
        if Self::graph_from_document(&doc)?!=self.graph{
            return Err(invalid("G7d restored semantic producer graph does not match archive"));
        }
        Ok(doc)
    }
    pub fn to_json(&self)->Result<String,GeometryError>{
        if self.schema_version!=SCHEMA_VERSION{
            return Err(invalid("G7d can serialize only archive schema 3"));
        }
        let result=serde_json::to_string_pretty(self)
            .map_err(|e|invalid(format!("G7d serialization failed: {e}")))?;
        if result.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G7d archive exceeds one MiB"));
        }
        Ok(result)
    }
    pub fn from_json(source:&str)->Result<Self,GeometryError>{
        if source.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G7d project exceeds bounded JSON size"));
        }
        let archive:Self=serde_json::from_str(source)
            .map_err(|e|invalid(format!("G7d invalid project JSON: {e}")))?;
        let _=archive.reopen()?;
        Ok(archive)
    }
    pub fn save_atomic(&self,target:impl AsRef<Path>)
        ->Result<(),GeometryError>{
        let encoded=self.to_json()?;
        let _=self.reopen()?;
        let target=target.as_ref();
        let parent=target.parent().filter(|p|!p.as_os_str().is_empty())
            .unwrap_or_else(||Path::new("."));
        if !parent.is_dir(){
            return Err(invalid("G7d target parent directory does not exist"));
        }
        let base=target.file_name().ok_or_else(||invalid("G7d target lacks filename"))?;
        let id=NEXT_TEMP.fetch_add(1,Ordering::Relaxed);
        if id==0||id==u64::MAX{
            return Err(invalid("G7d temporary file ID exhausted"));
        }
        let scratch=parent.join(format!(
            ".{}.g7d-{}-{id}.tmp",base.to_string_lossy(),std::process::id()));
        let status=(||->Result<(),GeometryError>{
            let mut output=OpenOptions::new().create_new(true)
                .write(true).open(&scratch)
                .map_err(|e|invalid(format!("G7d create temp failed: {e}")))?;
            output.write_all(encoded.as_bytes())
                .map_err(|e|invalid(format!("G7d write temp failed: {e}")))?;
            output.sync_all()
                .map_err(|e|invalid(format!("G7d sync temp failed: {e}")))?;
            drop(output);
            fs::rename(&scratch,target)
                .map_err(|e|invalid(format!("G7d atomic rename failed: {e}")))?;
            Ok(())
        })();
        if status.is_err(){
            let _=fs::remove_file(&scratch);
        }
        status?;
        // Directory fsync is best-effort/platform-dependent. Atomic rename
        // is not a distributed database or universal power-failure guarantee.
        let _=File::open(parent).and_then(|x|x.sync_all());
        Ok(())
    }
    pub fn load_file(source:impl AsRef<Path>)->Result<Self,GeometryError>{
        let mut input=File::open(source)
            .map_err(|e|invalid(format!("G7d open project failed: {e}")))?;
        let mut limited=(&mut input).take((MAX_ARCHIVE_BYTES+1) as u64);
        let mut bytes=Vec::new();
        limited.read_to_end(&mut bytes)
            .map_err(|e|invalid(format!("G7d read project failed: {e}")))?;
        if bytes.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G7d project exceeds 1 MiB"));
        }
        let string=std::str::from_utf8(&bytes)
            .map_err(|_|invalid("G7d project is not valid UTF-8 JSON"))?;
        Self::from_json(string)
    }
}
