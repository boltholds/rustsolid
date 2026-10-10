//! Topological names and scoped generational handles.
//!
//! Names are deterministic for the same feature key, primitive indexing and edit
//! history. They are **not** geometric matching across arbitrary regeneration.
//! Handles are process-local: generations prevent ABA when an uncommitted
//! transaction is dropped and a numeric topology slot is allocated again.
use crate::{EditDelta, FaceRole, GeometryError, Solid, TopologyEntity};
use std::collections::BTreeMap;
use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_INCARNATION: AtomicU64 = AtomicU64::new(1);

fn fresh_incarnation() -> u64 {
    // A single monotonic allocator is shared by body tokens and entity generations.
    // Neither value is a persistent ID to save on disk or transfer between processes.
    let id = NEXT_INCARNATION.fetch_add(1, Ordering::Relaxed);
    assert!(id != 0 && id != u64::MAX, "RustSolid process-local identity counter exhausted");
    id
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TopologyName(pub String);

impl TopologyName {
    pub fn as_str(&self) -> &str { &self.0 }
}
impl fmt::Display for TopologyName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { self.0.fmt(f) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TopologyHandle {
    pub body_token: u64,
    pub entity: TopologyEntity,
    pub generation: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HandleError { ForeignBody, StaleHandle }
impl fmt::Display for HandleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignBody => write!(f, "topological handle belongs to another body"),
            Self::StaleHandle => write!(f, "topological handle refers to a retired incarnation"),
        }
    }
}
impl std::error::Error for HandleError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameChange {
    pub entity: TopologyEntity,
    pub name: TopologyName,
    pub parent: TopologyName,
    pub handle: TopologyHandle,
}

#[derive(Debug, Clone)]
struct NamedSlot {
    name: TopologyName,
    generation: u64,
    parent: Option<TopologyName>,
}

#[derive(Debug, Clone)]
pub struct TopologyIdentity {
    body_token: u64,
    feature_key: String,
    entries: BTreeMap<TopologyEntity, NamedSlot>,
    names: BTreeMap<TopologyName, TopologyEntity>,
    allocations: Vec<TopologyEntity>,
    next_operation: u64,
}

impl TopologyIdentity {
    pub(crate) fn empty() -> Self {
        Self { body_token: fresh_incarnation(), feature_key: "primitive".into(),
            entries: BTreeMap::new(), names: BTreeMap::new(),
            allocations: Vec::new(), next_operation: 0 }
    }

    fn insert(&mut self, entity: TopologyEntity, name: TopologyName, parent: Option<TopologyName>) -> Result<TopologyHandle, GeometryError> {
        if self.entries.contains_key(&entity) || self.names.contains_key(&name) {
            return Err(GeometryError::InvalidTopology(format!("duplicate topological identity {entity:?} / {name}")));
        }
        let generation = fresh_incarnation();
        let handle = TopologyHandle { body_token: self.body_token, entity, generation };
        self.entries.insert(entity, NamedSlot { name: name.clone(), generation, parent });
        self.names.insert(name, entity);
        self.allocations.push(entity);
        Ok(handle)
    }

    pub(crate) fn initialize(&mut self, body: &Solid) -> Result<(), GeometryError> {
        if !self.entries.is_empty() { return Err(GeometryError::InvalidTopology("identity store already initialized".into())); }
        let base = self.feature_key.clone();
        self.insert(TopologyEntity::Body, TopologyName(format!("{base}/body")), None)?;
        for v in &body.vertices {
            self.insert(TopologyEntity::Vertex(v.id), TopologyName(format!("{base}/vertex/{}", v.id.0)), None)?;
        }
        for face in &body.faces {
            let role = match face.role {
                FaceRole::BottomCap => "bottom".into(),
                FaceRole::TopCap => "top".into(),
                FaceRole::Side(i) => format!("side/{i}"),
            };
            self.insert(TopologyEntity::Face(face.id), TopologyName(format!("{base}/face/{role}")), None)?;
        }
        for e in &body.edges {
            self.insert(TopologyEntity::Edge(e.id), TopologyName(format!("{base}/edge/{}", e.id.0)), None)?;
        }
        for c in &body.coedges {
            self.insert(TopologyEntity::Coedge(c.id), TopologyName(format!("{base}/coedge/{}", c.id.0)), None)?;
        }
        for l in &body.loops {
            self.insert(TopologyEntity::Loop(l.id), TopologyName(format!("{base}/loop/{}", l.id.0)), None)?;
        }
        for s in &body.shells {
            self.insert(TopologyEntity::Shell(s.id), TopologyName(format!("{base}/shell/{}", s.id.0)), None)?;
        }
        Ok(())
    }

    pub fn feature_key(&self) -> &str { &self.feature_key }
    pub fn name_of(&self, entity: TopologyEntity) -> Option<&TopologyName> {
        self.entries.get(&entity).map(|s| &s.name)
    }
    pub fn parent_of(&self, entity: TopologyEntity) -> Option<&TopologyName> {
        self.entries.get(&entity).and_then(|s| s.parent.as_ref())
    }
    pub fn find(&self, name: &TopologyName) -> Option<TopologyHandle> {
        self.names.get(name).and_then(|&entity| self.handle(entity))
    }
    pub fn handle(&self, entity: TopologyEntity) -> Option<TopologyHandle> {
        self.entries.get(&entity).map(|slot| TopologyHandle {
            body_token: self.body_token, entity, generation: slot.generation,
        })
    }
    pub fn resolve(&self, handle: TopologyHandle) -> Result<TopologyEntity, HandleError> {
        if handle.body_token != self.body_token { return Err(HandleError::ForeignBody); }
        match self.entries.get(&handle.entity) {
            Some(slot) if slot.generation == handle.generation => Ok(handle.entity),
            _ => Err(HandleError::StaleHandle),
        }
    }

    /// May only be called on a pristine primitive, before edit-derived names
    /// or references have escaped into a feature history.
    pub(crate) fn rename_feature(&mut self, key: &str) -> Result<(), GeometryError> {
        if key.is_empty() || key.len() > 128 || !key.bytes().all(|b|
            b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b'/')) {
            return Err(GeometryError::InvalidEdit("feature key must be 1..128 ASCII letters/digits/._-/".into()));
        }
        if self.next_operation != 0 {
            return Err(GeometryError::InvalidEdit("feature key cannot change after an edit".into()));
        }
        let old = format!("{}/", self.feature_key);
        let prefix = format!("{key}/");
        self.names.clear();
        for (&entity, slot) in &mut self.entries {
            let suffix = slot.name.0.strip_prefix(&old)
                .ok_or_else(|| GeometryError::InvalidTopology("inconsistent feature identity".into()))?;
            slot.name = TopologyName(format!("{prefix}{suffix}"));
            self.names.insert(slot.name.clone(), entity);
        }
        self.feature_key = key.into();
        Ok(())
    }

    pub(crate) fn checkpoint(&self) -> (usize, u64) { (self.allocations.len(), self.next_operation) }
    pub(crate) fn rewind(&mut self, (len, ordinal): (usize, u64)) {
        for entity in self.allocations.drain(len..) {
            if let Some(slot) = self.entries.remove(&entity) { self.names.remove(&slot.name); }
        }
        self.next_operation = ordinal;
    }

    pub(crate) fn register_delta(&mut self, delta: EditDelta) -> Result<Vec<NameChange>, GeometryError> {
        let (parent_entity, additions): (TopologyEntity, Vec<(TopologyEntity, &str)>) = match delta {
            EditDelta::SplitEdge(s) => (TopologyEntity::Edge(s.original_edge), vec![
                (TopologyEntity::Vertex(s.inserted_vertex), "vertex"),
                (TopologyEntity::Edge(s.created_edge), "remainder-edge"),
                (TopologyEntity::Coedge(s.created_coedges[0]), "coedge-forward"),
                (TopologyEntity::Coedge(s.created_coedges[1]), "coedge-reverse"),
            ]),
            EditDelta::SplitFace(s) => (TopologyEntity::Face(s.original_face), vec![
                (TopologyEntity::Face(s.created_face), "child-face"),
                (TopologyEntity::Loop(s.created_loop), "child-loop"),
                (TopologyEntity::Edge(s.diagonal_edge), "diagonal-edge"),
                (TopologyEntity::Coedge(s.created_coedges[0]), "coedge-forward"),
                (TopologyEntity::Coedge(s.created_coedges[1]), "coedge-reverse"),
            ]),
        };
        let parent = self.name_of(parent_entity).cloned()
            .ok_or_else(|| GeometryError::InvalidEdit("missing edit parent topological name".into()))?;
        let next = self.next_operation.checked_add(1)
            .ok_or_else(|| GeometryError::InvalidEdit("edit operation counter overflow".into()))?;
        // These keys are deterministic if the feature key and ordered operation
        // sequence are the same. The parent's name remains attached to its
        // retained child, and new children receive explicit provenance.
        let prefix = format!("{}/edit-{next}", parent.as_str());
        for &(entity, suffix) in &additions {
            let name = TopologyName(format!("{prefix}/{suffix}"));
            if self.entries.contains_key(&entity) || self.names.contains_key(&name) {
                return Err(GeometryError::InvalidEdit("edit would collide with an existing topology name".into()));
            }
        }
        let mut result = Vec::with_capacity(additions.len());
        for (entity, suffix) in additions {
            let name = TopologyName(format!("{prefix}/{suffix}"));
            let handle = self.insert(entity, name.clone(), Some(parent.clone()))?;
            result.push(NameChange { entity, name, parent: parent.clone(), handle });
        }
        self.next_operation = next;
        Ok(result)
    }

    pub(crate) fn validate(&self, body: &Solid) -> Result<(), GeometryError> {
        let expected = 1 + body.vertices.len() + body.edges.len() + body.coedges.len()
            + body.loops.len() + body.faces.len() + body.shells.len();
        if expected != self.entries.len() || self.names.len() != self.entries.len()
            || self.allocations.len() != self.entries.len() {
            return Err(GeometryError::InvalidTopology("identity registry cardinality mismatch".into()));
        }
        let mut expected_entities = Vec::with_capacity(expected);
        expected_entities.push(TopologyEntity::Body);
        expected_entities.extend(body.vertices.iter().map(|e| TopologyEntity::Vertex(e.id)));
        expected_entities.extend(body.edges.iter().map(|e| TopologyEntity::Edge(e.id)));
        expected_entities.extend(body.coedges.iter().map(|e| TopologyEntity::Coedge(e.id)));
        expected_entities.extend(body.loops.iter().map(|e| TopologyEntity::Loop(e.id)));
        expected_entities.extend(body.faces.iter().map(|e| TopologyEntity::Face(e.id)));
        expected_entities.extend(body.shells.iter().map(|e| TopologyEntity::Shell(e.id)));
        for entity in expected_entities {
            let slot = self.entries.get(&entity).ok_or_else(|| GeometryError::InvalidTopology(format!("unnamed topology: {entity:?}")))?;
            if self.names.get(&slot.name) != Some(&entity) {
                return Err(GeometryError::InvalidTopology(format!("name index disagrees with {entity:?}")));
            }
        }
        Ok(())
    }
}

impl Solid {
    /// Change the feature key on a newly constructed body. Call before edits.
    /// Reconstruction with the same key and equivalent indexed profile yields
    /// the same names, even though process-local handles will be different.
    pub fn with_feature_key(mut self, key: &str) -> Result<Self, GeometryError> {
        if self.revision != 0 { return Err(GeometryError::InvalidEdit("feature key is only assignable before editing".into())); }
        self.identities.rename_feature(key)?;
        Ok(self)
    }
    pub fn topology_name(&self, entity: TopologyEntity) -> Option<&TopologyName> {
        self.identities.name_of(entity)
    }
    pub fn topology_parent_name(&self, entity: TopologyEntity) -> Option<&TopologyName> {
        self.identities.parent_of(entity)
    }
    pub fn topology_handle(&self, entity: TopologyEntity) -> Option<TopologyHandle> {
        self.identities.handle(entity)
    }
    pub fn resolve_topology_name(&self, name: &TopologyName) -> Option<TopologyHandle> {
        self.identities.find(name)
    }
    pub fn resolve_topology_handle(&self, handle: TopologyHandle) -> Result<TopologyEntity, HandleError> {
        self.identities.resolve(handle)
    }
}
