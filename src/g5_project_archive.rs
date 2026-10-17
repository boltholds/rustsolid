//! G5b: versioned Gefest-bound four-shell project archive.
//!
//! Persist canonical source parameters, a typed three-node feature DAG,
//! compact certified-result stamps, and complete Undo/Redo state snapshots.
//! Never persist live handles or trust serialized B-rep meshes/vertices.
//! Reopen recomputes G4b2b3b/G4c against every historical source state.

use crate::{
    BrepBody,CylindricalBrep,G5BooleanDocument,GeometryError,
    GeometryTolerance,Point3,Solid,TopologyName,
};
use crate::atomic_boolean_document::{
    G5RecoveryHistory,G5RecoveryState,G5RecoveryStamp,
};
use serde::{Deserialize,Serialize};
use std::fs::{self,File,OpenOptions};
use std::io::{Read,Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64,Ordering};

const SCHEMA:u32=1;
const MAX_ARCHIVE_BYTES:usize=1024*1024;
const MAX_HISTORY:usize=64;
static NEXT_TEMP:AtomicU64=AtomicU64::new(1);

fn invalid(s:impl Into<String>)->GeometryError{
    GeometryError::InvalidEdit(s.into())
}
fn point(x:[f64;3])->Point3{
    Point3{x:x[0],y:x[1],z:x[2]}
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub enum G5ProjectSource {
    Block {
        origin:[f64;3],width:f64,height:f64,depth:f64,
        tolerance:GeometryTolerance,
    },
    UprightCylinder {
        origin:[f64;3],radius:f64,height:f64,
        tolerance:GeometryTolerance,
    },
}
impl G5ProjectSource {
    fn from_body(body:&BrepBody)->Result<Self,GeometryError>{
        body.validate()?;
        match body {
            BrepBody::Polyhedral(s)=>{
                crate::boolean::certify_box(s)
                    .map_err(|e|invalid(format!("G5b only canonical boxes: {e}")))?;
                let min=s.bbox.min;
                let delta=s.bbox.max.sub(min);
                let result=Self::Block{
                    origin:min.array(),width:delta.x,height:delta.y,
                    depth:delta.z,tolerance:s.tolerance,
                };
                result.to_body()?.validate()?;
                Ok(result)
            },
            BrepBody::Cylinder(c)=>{
                let rebuilt=CylindricalBrep::upright(
                    c.frame.origin,c.radius,c.height,c.tolerance)?;
                let budget=c.tolerance.length_at(
                    c.radius.max(c.height).max(1.0))*16.0;
                if c.frame.x.sub(rebuilt.frame.x).length_squared().sqrt()>budget
                    ||c.frame.y.sub(rebuilt.frame.y).length_squared().sqrt()>budget
                    ||c.frame.normal.sub(rebuilt.frame.normal)
                        .length_squared().sqrt()>budget{
                    return Err(invalid("G5b persistence requires upright cylinder chart"));
                }
                Ok(Self::UprightCylinder{
                    origin:c.frame.origin.array(),radius:c.radius,
                    height:c.height,tolerance:c.tolerance,
                })
            },
            _=>Err(invalid("G5b archive supports only original AABB/cylinder sources")),
        }
    }
    fn to_body(&self)->Result<BrepBody,GeometryError>{
        let result=match *self{
            Self::Block{origin,width,height,depth,tolerance}=>
                BrepBody::from(Solid::block_with_tolerance(
                    point(origin),width,height,depth,tolerance)?),
            Self::UprightCylinder{origin,radius,height,tolerance}=>
                BrepBody::from(CylindricalBrep::upright(
                    point(origin),radius,height,tolerance)?),
        };
        result.validate()?;
        Ok(result)
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5ResultStamp {
    pub mass_volume:f64,
    pub component_volumes:[f64;4],
    pub source_lower:f64,
    pub source_upper:f64,
    pub flux_volume:f64,
}
impl From<G5RecoveryStamp> for G5ResultStamp {
    fn from(x:G5RecoveryStamp)->Self{
        Self{mass_volume:x.mass_volume,component_volumes:x.component_volumes,
            source_lower:x.source_lower,source_upper:x.source_upper,
            flux_volume:x.flux_volume}
    }
}
impl From<G5ResultStamp> for G5RecoveryStamp {
    fn from(x:G5ResultStamp)->Self{
        Self{mass_volume:x.mass_volume,component_volumes:x.component_volumes,
            source_lower:x.source_lower,source_upper:x.source_upper,
            flux_volume:x.flux_volume}
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5ArchiveState {
    pub sources:[G5ProjectSource;2],
    pub result:Option<G5ResultStamp>,
}
impl G5ArchiveState {
    fn from_state(s:G5RecoveryState)->Result<Self,GeometryError>{
        Ok(Self{
            sources:[G5ProjectSource::from_body(&s.inputs[0])?,
                G5ProjectSource::from_body(&s.inputs[1])?],
            result:s.accepted.map(Into::into),
        })
    }
    fn to_state(&self)->Result<G5RecoveryState,GeometryError>{
        Ok(G5RecoveryState {
            inputs:[self.sources[0].to_body()?,self.sources[1].to_body()?],
            accepted:self.result.clone().map(Into::into),
        })
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq,Eq)]
#[serde(tag="kind",rename_all="snake_case",deny_unknown_fields)]
pub enum G5FeatureNode {
    Source { id:String,operand_index:u8 },
    Difference {
        id:String,left:String,right:String,
        output_scope:String,
    },
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq,Eq)]
#[serde(deny_unknown_fields)]
pub struct G5GraphBinding {
    /// Rebound to a fresh process-local handle after project open.
    pub topology_name:TopologyName,
    pub producing_node:String,
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq,Eq)]
#[serde(deny_unknown_fields)]
pub struct G5FeatureGraph {
    pub nodes:Vec<G5FeatureNode>,
    pub output_bindings:Vec<G5GraphBinding>,
}
impl G5FeatureGraph {
    fn from_document(doc:&G5BooleanDocument)->Result<Self,GeometryError>{
        let names=doc.archived_output_names()?;
        if names.len()!=0&&names.len()!=177 {
            return Err(invalid("G5b feature graph must own every result entity"));
        }
        let mut bindings=names.into_iter().map(|name|G5GraphBinding {
            topology_name:name,producing_node:"difference".into(),
        }).collect::<Vec<_>>();
        bindings.sort_by(|a,b|a.topology_name.cmp(&b.topology_name));
        Ok(Self{
            nodes:vec![
                G5FeatureNode::Source{id:"source_a".into(),operand_index:0},
                G5FeatureNode::Source{id:"source_b".into(),operand_index:1},
                G5FeatureNode::Difference{id:"difference".into(),
                    left:"source_a".into(),right:"source_b".into(),
                    output_scope:format!("{}/difference",doc.feature_key()),
                },
            ],
            output_bindings:bindings,
        })
    }
}
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
#[serde(deny_unknown_fields)]
pub struct G5ProjectArchive {
    pub schema_version:u32,
    pub feature_key:String,
    pub document_revision:u64,
    pub tolerance:GeometryTolerance,
    pub graph:G5FeatureGraph,
    pub current:G5ArchiveState,
    pub undo:Vec<G5ArchiveState>,
    pub redo:Vec<G5ArchiveState>,
}
impl G5ProjectArchive {
    pub fn from_document(doc:&G5BooleanDocument)
        ->Result<Self,GeometryError>{
        let saved=doc.export_recovery()?;
        if saved.undo.len()+saved.redo.len()>MAX_HISTORY{
            return Err(invalid("G5b project exceeded 64 history states"));
        }
        let result=Self{
            schema_version:SCHEMA,feature_key:saved.feature_key,
            document_revision:saved.revision,tolerance:saved.tolerance,
            graph:G5FeatureGraph::from_document(doc)?,
            current:G5ArchiveState::from_state(saved.current)?,
            undo:saved.undo.into_iter().map(G5ArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
            redo:saved.redo.into_iter().map(G5ArchiveState::from_state)
                .collect::<Result<Vec<_>,_>>()?,
        };
        let _=result.to_json()?; // enforce representability before export.
        Ok(result)
    }
    pub fn reopen(&self)->Result<G5BooleanDocument,GeometryError>{
        if self.schema_version!=SCHEMA{
            return Err(invalid(format!(
                "G5b unsupported project schema version {}",self.schema_version)));
        }
        if self.undo.len()+self.redo.len()>MAX_HISTORY
            ||self.document_revision==u64::MAX{
            return Err(invalid("G5b unsafe persisted history size or revision"));
        }
        self.tolerance.validate()?;
        let recovered=G5RecoveryHistory{
            feature_key:self.feature_key.clone(),
            tolerance:self.tolerance,
            revision:self.document_revision,
            current:self.current.to_state()?,
            undo:self.undo.iter().map(G5ArchiveState::to_state)
                .collect::<Result<Vec<_>,_>>()?,
            redo:self.redo.iter().map(G5ArchiveState::to_state)
                .collect::<Result<Vec<_>,_>>()?,
        };
        let document=G5BooleanDocument::restore_recovery(recovered)?;
        let expected=G5FeatureGraph::from_document(&document)?;
        if expected!=self.graph{
            return Err(invalid("G5b feature dependency graph/semantic names tampered or stale"));
        }
        document.validate()?;
        Ok(document)
    }
    /// Deterministic (single version) JSON: no B-rep arrays and NO handles.
    pub fn to_json(&self)->Result<String,GeometryError>{
        if self.schema_version!=SCHEMA{
            return Err(invalid("G5b cannot export unknown project schema version"));
        }
        let value=serde_json::to_string_pretty(self)
            .map_err(|e|invalid(format!("G5b JSON serialization failed: {e}")))?;
        if value.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G5b archive exceeds 1 MiB limit"));
        }
        Ok(value)
    }
    pub fn from_json(raw:&str)->Result<Self,GeometryError>{
        if raw.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G5b project exceeds maximum JSON size"));
        }
        let parsed:Self=serde_json::from_str(raw)
            .map_err(|e|invalid(format!("G5b invalid schema or JSON: {e}")))?;
        // Fail closed: every historical output and the graph must be fully
        // reconstructed independently before returning an archive.
        let _=parsed.reopen()?;
        Ok(parsed)
    }
    pub fn save_atomic(&self,target:impl AsRef<Path>)
        ->Result<(),GeometryError>{
        let encoded=self.to_json()?;
        // Disk effects are only attempted after successful in-memory G4c
        // reconstruction of EVERY archived state.
        self.reopen()?;
        let path=target.as_ref();
        let parent=path.parent().filter(|x|!x.as_os_str().is_empty())
            .unwrap_or_else(||Path::new("."));
        if !parent.is_dir(){
            return Err(invalid("G5b parent directory does not exist"));
        }
        let id=NEXT_TEMP.fetch_add(1,Ordering::Relaxed);
        if id==0||id==u64::MAX{
            return Err(invalid("G5b temporary file counter exhausted"));
        }
        let destination=path.file_name()
            .ok_or_else(||invalid("G5b invalid target filename"))?;
        let candidate=parent.join(format!(
            ".{}.g5b-{}-{id}.tmp",
            destination.to_string_lossy(),std::process::id()));
        let result=(||->Result<(),GeometryError>{
            let mut file=OpenOptions::new()
                .write(true).create_new(true).open(&candidate)
                .map_err(|e|invalid(format!("G5b create temp failed: {e}")))?;
            file.write_all(encoded.as_bytes())
                .map_err(|e|invalid(format!("G5b temp write failed: {e}")))?;
            file.sync_all()
                .map_err(|e|invalid(format!("G5b temp fsync failed: {e}")))?;
            drop(file);
            // Atomic replacement on platforms supporting same-directory
            // atomic rename. No application-wide database transaction.
            fs::rename(&candidate,path)
                .map_err(|e|invalid(format!("G5b atomic rename failed: {e}")))?;
            Ok(())
        })();
        if result.is_err(){let _=fs::remove_file(&candidate);}
        result?;
        // Directory durability is platform-dependent and best-effort; the
        // atomicity guarantee concerns readers of the final path.
        let _=File::open(parent).and_then(|dir|dir.sync_all());
        Ok(())
    }
    pub fn load_file(path:impl AsRef<Path>)->Result<Self,GeometryError>{
        let mut file=File::open(path)
            .map_err(|e|invalid(format!("G5b cannot open project: {e}")))?;
        let mut bounded=(&mut file).take((MAX_ARCHIVE_BYTES+1) as u64);
        let mut bytes=Vec::new();
        bounded.read_to_end(&mut bytes)
            .map_err(|e|invalid(format!("G5b project read failed: {e}")))?;
        if bytes.len()>MAX_ARCHIVE_BYTES{
            return Err(invalid("G5b file exceeds 1 MiB limit"));
        }
        let input=std::str::from_utf8(&bytes)
            .map_err(|_|invalid("G5b project is not UTF-8 JSON"))?;
        Self::from_json(input)
    }
}
