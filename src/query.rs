//! A deliberately bounded, read-only GraphQL-inspired language for B-rep queries.
//!
//! This is NOT the GraphQL specification: it has a typed selection grammar,
//! field arguments and nested traversal, but intentionally omits mutations,
//! directives, fragments, variables, subscriptions and custom resolvers.
//! Model validation and query complexity checks precede any traversal.
use crate::{BrepModel, BrepOrigin, CoedgeId, Curve2, Curve3, EdgeId, FaceId,
    GeometryError, LoopId, Point3, Surface3};
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct QueryLimits {
    pub max_bytes: usize,
    pub max_depth: usize,
    pub max_fields: usize,
    pub max_items: usize,
    pub max_visits: usize,
}
impl Default for QueryLimits {
    fn default() -> Self {
        Self { max_bytes: 4096, max_depth: 8, max_fields: 128,
            max_items: 128, max_visits: 2048 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryError {
    pub code: &'static str,
    pub offset: usize,
    pub detail: String,
}
impl fmt::Display for QueryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at byte {}: {}", self.code, self.offset, self.detail)
    }
}
impl std::error::Error for QueryError {}
fn qerr(code: &'static str, offset: usize, message: impl Into<String>) -> QueryError {
    QueryError { code, offset, detail: message.into() }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Literal { Text(String), Integer(u32), Boolean(bool) }
#[derive(Debug, Clone)]
struct Selection {
    name: String, at: usize,
    args: BTreeMap<String, Literal>,
    fields: Vec<Selection>,
}
#[derive(Debug, Clone)]
enum TokenKind { Word(String), Integer(u32), String(String), Open, Close, LParen, RParen, Colon, End }
#[derive(Debug, Clone)]
struct Token { kind: TokenKind, at: usize }

fn tokenize(src: &str, limits: QueryLimits) -> Result<Vec<Token>, QueryError> {
    if src.len() > limits.max_bytes { return Err(qerr("query.too_large", 0, "query exceeds byte limit")); }
    let bytes = src.as_bytes();
    let mut i = 0;
    let mut tokens = Vec::new();
    while i < bytes.len() {
        let start = i;
        match bytes[i] {
            b' ' | b'\t' | b'\r' | b'\n' | b',' => { i += 1; continue; }
            b'#' => { while i < bytes.len() && bytes[i] != b'\n' { i += 1; } continue; }
            b'{' | b'}' | b'(' | b')' | b':' => {
                let kind = match bytes[i] { b'{' => TokenKind::Open, b'}' => TokenKind::Close,
                    b'(' => TokenKind::LParen, b')' => TokenKind::RParen, _ => TokenKind::Colon };
                i += 1;
                tokens.push(Token { kind, at:start });
            }
            b'"' => {
                i += 1;
                let mut escaped = false;
                while i < bytes.len() {
                    let current = bytes[i];
                    i += 1;
                    if current == b'"' && !escaped { break; }
                    if current == b'\\' && !escaped { escaped = true; }
                    else { escaped = false; }
                }
                let raw = &src[start..i];
                let text: String = serde_json::from_str(raw)
                    .map_err(|_| qerr("query.bad_string", start, "invalid quoted string"))?;
                tokens.push(Token {kind:TokenKind::String(text),at:start});
            }
            b'0'..=b'9' => {
                while i < bytes.len() && bytes[i].is_ascii_digit() { i += 1; }
                let number = src[start..i].parse::<u32>()
                    .map_err(|_| qerr("query.bad_integer", start, "integer exceeds u32"))?;
                tokens.push(Token {kind:TokenKind::Integer(number),at:start});
            }
            b'A'..=b'Z' | b'a'..=b'z' | b'_' => {
                i += 1;
                while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i]==b'_') { i += 1; }
                tokens.push(Token {kind:TokenKind::Word(src[start..i].into()),at:start});
            }
            _ => return Err(qerr("query.invalid_character", start, "unexpected character")),
        }
        if tokens.len() > 512 { return Err(qerr("query.too_many_tokens", start, "token budget exhausted")); }
    }
    tokens.push(Token { kind:TokenKind::End, at:bytes.len() });
    Ok(tokens)
}

struct Parser { tokens: Vec<Token>, position: usize, fields: usize, limits: QueryLimits }
impl Parser {
    fn peek(&self) -> &Token { &self.tokens[self.position] }
    fn take(&mut self) -> Token { let t = self.tokens[self.position].clone(); self.position += 1; t }
    fn expect(&mut self, expected: fn(&TokenKind)->bool, detail: &'static str) -> Result<Token, QueryError> {
        if expected(&self.peek().kind) { Ok(self.take()) }
        else { Err(qerr("query.syntax", self.peek().at, detail)) }
    }
    fn word(&mut self) -> Result<(String, usize), QueryError> {
        let token = self.expect(|k| matches!(k,TokenKind::Word(_)), "expected field or argument name")?;
        let TokenKind::Word(value) = token.kind else { unreachable!() };
        Ok((value,token.at))
    }
    fn selections(&mut self, depth:usize) -> Result<Vec<Selection>, QueryError> {
        if depth > self.limits.max_depth {
            return Err(qerr("query.depth", self.peek().at, "maximum query nesting depth exceeded"));
        }
        self.expect(|k|matches!(k,TokenKind::Open), "expected '{'")?;
        let mut fields = Vec::new();
        while !matches!(self.peek().kind, TokenKind::Close | TokenKind::End) {
            let (name, at) = self.word()?;
            self.fields += 1;
            if self.fields > self.limits.max_fields {
                return Err(qerr("query.fields", at, "field selection budget exhausted"));
            }
            let mut args = BTreeMap::new();
            if matches!(self.peek().kind, TokenKind::LParen) {
                self.take();
                while !matches!(self.peek().kind,TokenKind::RParen | TokenKind::End) {
                    let (key,key_at) = self.word()?;
                    self.expect(|k|matches!(k,TokenKind::Colon), "expected ':' in argument")?;
                    let value = match self.take() {
                        Token { kind:TokenKind::Integer(n), .. } => Literal::Integer(n),
                        Token { kind:TokenKind::String(s), .. } => Literal::Text(s),
                        Token { kind:TokenKind::Word(w), at } if w=="true" || w=="false" => Literal::Boolean(w=="true"),
                        Token { kind:TokenKind::Word(w), .. } => Literal::Text(w),
                        Token { at, .. } => return Err(qerr("query.argument",at,"invalid argument literal")),
                    };
                    if args.insert(key.clone(),value).is_some() {
                        return Err(qerr("query.argument",key_at,format!("duplicate argument '{key}'")));
                    }
                }
                self.expect(|k|matches!(k,TokenKind::RParen),"expected ')' after arguments")?;
            }
            let children = if matches!(self.peek().kind,TokenKind::Open) {
                self.selections(depth+1)?
            } else { Vec::new() };
            if fields.iter().any(|f:&Selection|f.name==name) {
                return Err(qerr("query.field",at,format!("duplicate field '{name}' (aliases unsupported)")));
            }
            fields.push(Selection {name,at,args,fields:children});
        }
        self.expect(|k| matches!(k,TokenKind::Close), "expected '}'")?;
        if fields.is_empty() { return Err(qerr("query.empty",self.peek().at,"empty field selection")); }
        Ok(fields)
    }
}

#[derive(Debug, Clone, Copy)]
enum Kind { Root, Body, Topology, Mass, Bounds, Point, Face, Edge, Vertex, Shell, Loop, Coedge, Surface, Curve, Pcurve }
#[derive(Debug, Clone, Copy)]
enum Shape { Scalar, Object(Kind), List(Kind) }
fn shape(kind:Kind, field:&str) -> Option<Shape> {
    use Kind::*;
    use Shape::*;
    match kind {
        Root => match field {"body"=>Some(Object(Body)), _=>None},
        Body => match field {
            "kind"=>Some(Scalar),"topology"=>Some(Object(Topology)),"mass"=>Some(Object(Mass)),
            "bounds"=>Some(Object(Bounds)),"faces"=>Some(List(Face)),"edges"=>Some(List(Edge)),
            "vertices"=>Some(List(Vertex)),"shells"=>Some(List(Shell)),_=>None
        },
        Topology => match field {"vertexCount"|"edgeCount"|"coedgeCount"|"loopCount"|
            "faceCount"|"shellCount"|"eulerCharacteristic"|"genus"|"closedEdges"|
            "seamEdges"|"circularEdges"|"cylindricalFaces"=>Some(Scalar),_=>None},
        Mass => match field {"volume"|"surfaceArea"=>Some(Scalar),"centroid"=>Some(Object(Point)),_=>None},
        Bounds => match field {"min"|"max"=>Some(Object(Point)),_=>None},
        Point => match field {"x"|"y"|"z"=>Some(Scalar),_=>None},
        Face => match field {"id"=>Some(Scalar),"surface"=>Some(Object(Surface)),"loops"=>Some(List(Loop)),_=>None},
        Edge => match field {"id"|"startId"|"endId"|"closed"|"seam"=>Some(Scalar),
            "curve"=>Some(Object(Curve)),"coedges"=>Some(List(Coedge)),_=>None},
        Vertex => match field {"id"=>Some(Scalar),"position"=>Some(Object(Point)),_=>None},
        Shell => match field {"id"|"closed"=>Some(Scalar),"faces"=>Some(List(Face)),_=>None},
        Loop => match field {"id"|"role"=>Some(Scalar),"coedges"=>Some(List(Coedge)),_=>None},
        Coedge => match field {"id"|"edgeId"|"faceId"|"loopId"|"nextId"|"prevId"|"twinId"|
            "reversed"=>Some(Scalar),"edge"=>Some(Object(Edge)),"face"=>Some(Object(Face)),
            "pcurve"=>Some(Object(Pcurve)),_=>None},
        Surface => match field {"kind"|"radius"=>Some(Scalar),"origin"|"normal"=>Some(Object(Point)),_=>None},
        Curve => match field {"kind"|"radius"=>Some(Scalar),"origin"|"direction"=>Some(Object(Point)),_=>None},
        Pcurve => match field {"kind"=>Some(Scalar),_=>None},
    }
}
fn validate_fields(kind:Kind, fields:&[Selection]) -> Result<(),QueryError> {
    for field in fields {
        let Some(s) = shape(kind,&field.name) else {
            return Err(qerr("query.unknown_field",field.at,format!("'{0}' is not a field of {kind:?}",field.name)));
        };
        for (name,value) in &field.args {
            let valid = match (s,name.as_str(),value) {
                (Shape::List(_),"limit",Literal::Integer(_)) => true,
                (Shape::List(_),"id",Literal::Integer(_)) => true,
                (Shape::List(Kind::Face),"kind",Literal::Text(v)) if ["plane","cylinder","sphere"].contains(&v.as_str()) => true,
                (Shape::List(Kind::Edge),"curve",Literal::Text(v)) if ["line","circle"].contains(&v.as_str()) => true,
                (Shape::List(Kind::Edge),"seam"|"closed",Literal::Boolean(_)) => true,
                _=>false,
            };
            if !valid { return Err(qerr("query.invalid_argument",field.at,
                format!("unsupported or incorrectly typed argument '{name}' on '{}'",field.name))); }
        }
        match s {
            Shape::Scalar if !field.fields.is_empty() => return Err(qerr("query.scalar_selection",field.at,"scalar has no nested fields")),
            Shape::Scalar => {},
            Shape::Object(child)|Shape::List(child) if field.fields.is_empty() =>
                return Err(qerr("query.missing_selection",field.at,"object/list requires subfields")),
            Shape::Object(child)|Shape::List(child) => validate_fields(child,&field.fields)?,
        }
    }
    Ok(())
}
fn int_arg(s:&Selection,key:&str)->Option<usize> {
    match s.args.get(key) {Some(Literal::Integer(v))=>Some(*v as usize),_=>None}
}
fn str_arg<'a>(s:&'a Selection,key:&str)->Option<&'a str> {
    match s.args.get(key) {Some(Literal::Text(v))=>Some(v),_=>None}
}
fn bool_arg(s:&Selection,key:&str)->Option<bool> {
    match s.args.get(key) {Some(Literal::Boolean(v))=>Some(*v),_=>None}
}

#[derive(Clone, Copy)]
enum Node { Root, Body, Topology, Mass, Bounds, Point(Point3),
    Face(FaceId), Edge(EdgeId), Vertex(usize), Shell(usize), Loop(LoopId), Coedge(CoedgeId),
    Surface(Surface3), Curve(Curve3), Pcurve(Curve2) }
struct Execution<'a> { model:&'a BrepModel, limits:QueryLimits, visits:usize }
impl Execution<'_> {
    fn get(&mut self,node:Node, fields:&[Selection]) -> Result<Value,QueryError> {
        self.visits += 1;
        if self.visits > self.limits.max_visits {
            return Err(qerr("query.cost",0,"maximum materialized objects exceeded"));
        }
        let mut output=Map::new();
        for s in fields {
            let value=self.field(node,s)?;
            output.insert(s.name.clone(),value);
        }
        Ok(Value::Object(output))
    }
    fn collection(&mut self, values:Vec<Node>, field:&Selection) -> Result<Value,QueryError> {
        let limit=int_arg(field,"limit").unwrap_or(32);
        if limit==0 || limit>self.limits.max_items {
            return Err(qerr("query.limit",field.at,"list limit outside permitted range"));
        }
        let desired_id=int_arg(field,"id");
        let mut output=Vec::new();
        for (index,node) in values.into_iter().enumerate() {
            let id=match node {
                Node::Face(id)=>id.0 as usize, Node::Edge(id)=>id.0 as usize,
                Node::Vertex(id)|Node::Shell(id)=>id,
                Node::Loop(id)=>id.0 as usize,Node::Coedge(id)=>id.0 as usize,
                _=>index,
            };
            if desired_id.is_some_and(|want|want != id) { continue; }
            if let Some(kind)=str_arg(field,"kind") {
                if let Node::Face(id) = node {
                    if surface_kind(self.model.face_surface(id))!=kind {continue;}
                }
            }
            if let Some(curve)=str_arg(field,"curve") {
                if let Node::Edge(id) = node {
                    if curve_kind(self.model.edge_curve(id))!=curve {continue;}
                }
            }
            if let Some(seam)=bool_arg(field,"seam") {
                if let Node::Edge(id)=node {if edge_seam(self.model,id)!=seam{continue;}}
            }
            if let Some(closed)=bool_arg(field,"closed") {
                if let Node::Edge(id)=node {
                    let e=&self.model.edges[id.0 as usize];
                    if (e.start==e.end)!=closed {continue;}
                }
            }
            output.push(self.get(node,&field.fields)?);
            if output.len()==limit {break;}
        }
        Ok(Value::Array(output))
    }
    fn field(&mut self,node:Node,s:&Selection) -> Result<Value,QueryError> {
        let m=self.model;
        let bad=||qerr("query.internal",s.at,"field resolver and schema disagree");
        let result=match (node,s.name.as_str()) {
            (Node::Root,"body")=>self.get(Node::Body,&s.fields)?,
            (Node::Body,"kind")=>json!(match m.origin {BrepOrigin::Polyhedral=>"polyhedral",BrepOrigin::AnalyticCylinder=>"analytic_cylinder"}),
            (Node::Body,"topology")=>self.get(Node::Topology,&s.fields)?,
            (Node::Body,"mass")=>self.get(Node::Mass,&s.fields)?,
            (Node::Body,"bounds")=>self.get(Node::Bounds,&s.fields)?,
            (Node::Body,"faces")=>self.collection(m.faces.iter().map(|f|Node::Face(f.id)).collect(),s)?,
            (Node::Body,"edges")=>self.collection(m.edges.iter().map(|e|Node::Edge(e.id)).collect(),s)?,
            (Node::Body,"vertices")=>self.collection((0..m.vertices.len()).map(Node::Vertex).collect(),s)?,
            (Node::Body,"shells")=>self.collection((0..m.shells.len()).map(Node::Shell).collect(),s)?,
            (Node::Topology,k)=>{
                let a=m.summary();
                match k {
                    "vertexCount"=>json!(a.vertices),"edgeCount"=>json!(a.edges),"coedgeCount"=>json!(a.coedges),
                    "loopCount"=>json!(a.loops),"faceCount"=>json!(a.faces),"shellCount"=>json!(a.shells),
                    "eulerCharacteristic"=>json!(a.euler_characteristic),"genus"=>json!(a.genus),
                    "closedEdges"=>json!(a.closed_edges),"seamEdges"=>json!(a.seam_edges),
                    "circularEdges"=>json!(a.circular_edges),"cylindricalFaces"=>json!(a.cylindrical_faces),
                    _=>return Err(bad()),
                }
            }
            (Node::Mass,"volume")=>json!(m.mass.volume),
            (Node::Mass,"surfaceArea")=>json!(m.mass.surface_area),
            (Node::Mass,"centroid")=>self.get(Node::Point(m.mass.centroid),&s.fields)?,
            (Node::Bounds,"min")=>self.get(Node::Point(m.bbox.min),&s.fields)?,
            (Node::Bounds,"max")=>self.get(Node::Point(m.bbox.max),&s.fields)?,
            (Node::Point(p),"x")=>json!(p.x),(Node::Point(p),"y")=>json!(p.y),(Node::Point(p),"z")=>json!(p.z),
            (Node::Face(id),"id")=>json!(id.0),
            (Node::Face(id),"surface")=>self.get(Node::Surface(m.face_surface(id).ok_or_else(bad)?),&s.fields)?,
            (Node::Face(id),"loops")=>{
                let face=&m.faces[id.0 as usize];
                self.collection(face.loops.iter().copied().map(Node::Loop).collect(),s)?
            }
            (Node::Edge(id),"id")=>json!(id.0),
            (Node::Edge(id),"startId")=>json!(m.edges[id.0 as usize].start.0),
            (Node::Edge(id),"endId")=>json!(m.edges[id.0 as usize].end.0),
            (Node::Edge(id),"closed")=>json!(m.edges[id.0 as usize].start==m.edges[id.0 as usize].end),
            (Node::Edge(id),"seam")=>json!(edge_seam(m,id)),
            (Node::Edge(id),"curve")=>self.get(Node::Curve(m.edge_curve(id).ok_or_else(bad)?),&s.fields)?,
            (Node::Edge(id),"coedges")=>{
                let edge=&m.edges[id.0 as usize];
                self.collection(edge.coedges.iter().copied().map(Node::Coedge).collect(),s)?
            }
            (Node::Vertex(idx),"id")=>json!(m.vertices[idx].id.0),
            (Node::Vertex(idx),"position")=>self.get(Node::Point(m.vertices[idx].position),&s.fields)?,
            (Node::Shell(idx),"id")=>json!(m.shells[idx].id.0),
            (Node::Shell(idx),"closed")=>json!(m.shells[idx].closed),
            (Node::Shell(idx),"faces")=>self.collection(m.shells[idx].faces.iter().copied().map(Node::Face).collect(),s)?,
            (Node::Loop(id),"id")=>json!(id.0),
            (Node::Loop(id),"role")=>json!(match m.loops[id.0 as usize].role {crate::LoopRole::Outer=>"outer",crate::LoopRole::Inner=>"inner"}),
            (Node::Loop(id),"coedges")=>{
                let lp=&m.loops[id.0 as usize];
                let first=lp.first_coedge;
                let mut cur=first;
                let mut ids=Vec::new();
                for _ in 0..m.coedges.len() {
                    ids.push(Node::Coedge(cur));
                    cur=m.coedges[cur.0 as usize].next;
                    if cur==first {break;}
                }
                self.collection(ids,s)?
            }
            (Node::Coedge(id),"id")=>json!(id.0),
            (Node::Coedge(id),"edgeId")=>json!(m.coedges[id.0 as usize].edge.0),
            (Node::Coedge(id),"faceId")=>json!(m.coedges[id.0 as usize].face.0),
            (Node::Coedge(id),"loopId")=>json!(m.coedges[id.0 as usize].loop_id.0),
            (Node::Coedge(id),"nextId")=>json!(m.coedges[id.0 as usize].next.0),
            (Node::Coedge(id),"prevId")=>json!(m.coedges[id.0 as usize].prev.0),
            (Node::Coedge(id),"twinId")=>json!(m.coedges[id.0 as usize].twin.0),
            (Node::Coedge(id),"reversed")=>json!(m.coedges[id.0 as usize].reversed),
            (Node::Coedge(id),"edge")=>self.get(Node::Edge(m.coedges[id.0 as usize].edge),&s.fields)?,
            (Node::Coedge(id),"face")=>self.get(Node::Face(m.coedges[id.0 as usize].face),&s.fields)?,
            (Node::Coedge(id),"pcurve")=>self.get(Node::Pcurve(m.coedge_pcurve(id).ok_or_else(bad)?),&s.fields)?,
            (Node::Surface(surf),"kind")=>json!(surface_kind(Some(surf))),
            (Node::Surface(surf),"radius")=>match surf {Surface3::Cylinder(x)=>json!(x.radius),Surface3::Sphere(x)=>json!(x.radius),_=>Value::Null},
            (Node::Surface(surf),"origin")=>{
                let p=match surf {Surface3::Plane(x)=>x.frame.origin,Surface3::Cylinder(x)=>x.frame.origin,Surface3::Sphere(x)=>x.frame.origin};
                self.get(Node::Point(p),&s.fields)?
            }
            (Node::Surface(surf),"normal")=>{
                let p=match surf {Surface3::Plane(x)=>x.frame.normal,Surface3::Cylinder(x)=>x.frame.normal,Surface3::Sphere(x)=>x.frame.normal};
                self.get(Node::Point(p),&s.fields)?
            }
            (Node::Curve(curve),"kind")=>json!(curve_kind(Some(curve))),
            (Node::Curve(curve),"radius")=>match curve {Curve3::Circle(x)=>json!(x.radius),_=>Value::Null},
            (Node::Curve(curve),"origin")=>{
                let p=match curve {Curve3::Line(x)=>x.origin,Curve3::Circle(x)=>x.frame.origin};
                self.get(Node::Point(p),&s.fields)?
            }
            (Node::Curve(curve),"direction")=>{
                let p=match curve {Curve3::Line(x)=>x.direction,Curve3::Circle(x)=>x.frame.normal};
                self.get(Node::Point(p),&s.fields)?
            }
            (Node::Pcurve(trim),"kind")=>json!(match trim {Curve2::Line(_)=>"line",Curve2::Circle(_)=>"circle"}),
            _=>return Err(bad()),
        };
        Ok(result)
    }
}
fn surface_kind(surface:Option<Surface3>)->&'static str {
    match surface {Some(Surface3::Plane(_))=>"plane",Some(Surface3::Cylinder(_))=>"cylinder",
        Some(Surface3::Sphere(_))=>"sphere",None=>"missing"}
}
fn curve_kind(curve:Option<Curve3>)->&'static str {
    match curve {Some(Curve3::Line(_))=>"line",Some(Curve3::Circle(_))=>"circle",None=>"missing"}
}
fn edge_seam(m:&BrepModel,id:EdgeId)->bool {
    let Some(edge)=m.edges.get(id.0 as usize) else {return false};
    let [a,b]=edge.coedges;
    m.coedges[a.0 as usize].face==m.coedges[b.0 as usize].face
}

/// Evaluate a typed field selection against a *validated* shared B-rep view.
/// Dynamic complexity and output size remain bounded; no mutations or plugins.
pub fn geometry_query(model:&BrepModel, source:&str, limits:QueryLimits) -> Result<Value,QueryError> {
    model.validate().map_err(|e:GeometryError|qerr("query.invalid_model",0,e.to_string()))?;
    let tokens=tokenize(source,limits)?;
    let mut parser=Parser{tokens,position:0,fields:0,limits};
    if matches!(&parser.peek().kind,TokenKind::Word(s) if s=="query") {
        parser.take();
        // An optional operation name may follow `query`, like GraphQL.
        if matches!(parser.peek().kind,TokenKind::Word(_)) {parser.take();}
    }
    let selections=parser.selections(1)?;
    if !matches!(parser.peek().kind,TokenKind::End) {
        return Err(qerr("query.syntax",parser.peek().at,"unexpected tokens after selection"));
    }
    validate_fields(Kind::Root,&selections)?;
    let mut context=Execution{model,limits,visits:0};
    let value=context.get(Node::Root,&selections)?;
    // Defense-in-depth on unexpectedly large nested output.
    if serde_json::to_vec(&value).map_or(true, |v|v.len()>256*1024) {
        return Err(qerr("query.output_size",0,"response exceeds 256 KiB"));
    }
    Ok(value)
}
