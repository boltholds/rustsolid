//! Independent primitive parity experiment. No SolidCore geometry is installed into
//! the RustSolid canonical B-rep or committed to Gefest.
//
// NOTE: SolidCore is a port of NVIDIA USD BRep's SMLib, Apache-2.0.
// See https://github.com/tangency-labs/solidcore (NOTICE/ATTRIBUTION.md).

use rustsolid::{Point3 as RustPoint3, Solid};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{
    brep_properties, MassPropertyBehavior, PropertyAttributes,
};
use solidcore::math::Point3 as SmPoint3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::query::brep::{calculate_bounding_box, BoundingBoxOptions};
use solidcore::topology::brep::{p1 as brep_queries, p2::is_manifold_solid};

/// Coordinates and dimensions in common model units: X, Y, Z.
///
/// RustSolid::block(origin, width, height, depth) uses X, Y, Z respectively.
/// SolidCore::create_box(origin, length, width, height) also uses X, Y, Z.
#[derive(Debug, Clone, Copy)]
struct BoxCase {
    name: &'static str,
    origin: [f64; 3],
    size: [f64; 3],
}

const CASES: [BoxCase; 3] = [
    BoxCase {
        name: "origin-5x10x20",
        origin: [0.0, 0.0, 0.0],
        size: [5.0, 10.0, 20.0],
    },
    BoxCase {
        name: "translated-2x3x4",
        origin: [11.0, -17.0, 2.0],
        size: [2.0, 3.0, 4.0],
    },
    BoxCase {
        name: "offset-small",
        origin: [-250.0, 700.0, -30.0],
        size: [0.5, 0.75, 1.25],
    },
];

#[derive(Debug)]
struct ParityReport {
    name: &'static str,
    volume: (f64, f64),
    area: (f64, f64),
    topological_counts: (usize, usize, usize),
}

fn fail(op: &str, e: impl std::fmt::Debug) -> String {
    format!("{op}: {e:?}")
}

/// Check a value against an independent formula, using a strict tolerance in
/// the scale of the individual model rather than absolute world coordinates.
fn nearly_equal(label: &str, actual: f64, expected: f64, scale: f64) -> Result<(), String> {
    let threshold = 1.0e-8 * scale.max(1.0);
    if !actual.is_finite() || (actual - expected).abs() > threshold {
        return Err(format!(
            "{label}: actual={actual:?}, expected={expected:?}, tolerance={threshold:?}"
        ));
    }
    Ok(())
}

fn compare_box(case: &BoxCase) -> Result<ParityReport, String> {
    let [x, y, z] = case.size;
    let [ox, oy, oz] = case.origin;
    if !case.origin.iter().chain(case.size.iter()).all(|v| v.is_finite())
        || case.size.iter().any(|v| *v <= 0.0)
    {
        return Err(format!("{}: invalid source specification", case.name));
    }

    // RustSolid: its own source-built B-rep, validated before interrogation.
    let rust_solid = Solid::block(
        RustPoint3 { x: ox, y: oy, z: oz },
        x, y, z,
    )
    .map_err(|e| fail("RustSolid::block", e))?;
    rust_solid
        .validate()
        .map_err(|e| fail("RustSolid::validate", e))?;

    // SolidCore: independent NVIDIA SMLib port, never converted into RustSolid.
    let mut kernel = Kernel::default();
    let mut sm_solid = create_box(&mut kernel, SmPoint3::new(ox, oy, oz), x, y, z)
        .map_err(|e| fail("SolidCore::create_box", e))?;
    let manifold = is_manifold_solid(&sm_solid)
        .map_err(|e| fail("SolidCore::is_manifold_solid", e))?;
    if !manifold {
        return Err(format!("{}: SolidCore result is not manifold", case.name));
    }

    let sm_faces = brep_queries::faces(&sm_solid, None)
        .map_err(|e| fail("SolidCore::faces", e))?.len();
    let sm_edges = brep_queries::edges(&sm_solid, None)
        .map_err(|e| fail("SolidCore::edges", e))?.len();
    let sm_vertices = brep_queries::vertices(&sm_solid, None)
        .map_err(|e| fail("SolidCore::vertices", e))?.len();

    let rust_counts = (
        rust_solid.faces.len(),
        rust_solid.edges.len(),
        rust_solid.vertices.len(),
    );
    let sm_counts = (sm_faces, sm_edges, sm_vertices);
    if rust_counts != sm_counts || rust_counts != (6, 12, 8) {
        return Err(format!(
            "{}: differing B-rep cardinality: RustSolid={rust_counts:?}, SolidCore={sm_counts:?}",
            case.name
        ));
    }

    // Independently query the tight SolidCore bounding box, excluding tolerance
    // expansion. This may lazily populate internal SolidCore trim caches.
    let tight = calculate_bounding_box(
        &kernel,
        &mut sm_solid,
        BoundingBoxOptions {
            approx: false,
            tight: true,
            expand_by_zone_tol_3d: false,
        },
    )
    .map_err(|e| fail("SolidCore::bounding_box", e))?;

    for axis in 0..3 {
        let expected_min = case.origin[axis];
        let expected_max = case.origin[axis] + case.size[axis];
        nearly_equal("SolidCore bbox min", tight.min_at(axis), expected_min, case.size[axis])?;
        nearly_equal("SolidCore bbox max", tight.max_at(axis), expected_max, case.size[axis])?;
        let rust_min = rust_solid.bbox.min.array()[axis];
        let rust_max = rust_solid.bbox.max.array()[axis];
        nearly_equal("RustSolid bbox min", rust_min, expected_min, case.size[axis])?;
        nearly_equal("RustSolid bbox max", rust_max, expected_max, case.size[axis])?;
    }

    // SolidCore mass/area are integrated from its own B-rep. Attributes are
    // empty (uniform unit density); the queried extent provides the origin.
    let store = AttributeStore::new();
    let key_of = |_: solidcore::EntityRef| ObjectKey(0);
    let attributes = PropertyAttributes { store: &store, key_of: &key_of };
    let properties = brep_properties(
        &kernel,
        &mut sm_solid,
        1.0e-10,
        tight.mid(),
        MassPropertyBehavior::default(),
        &attributes,
    )
    .map_err(|e| fail("SolidCore::brep_properties", e))?;

    let analytic_volume = x * y * z;
    let analytic_area = 2.0 * (x * y + x * z + y * z);
    for (label, actual, expected) in [
        ("RustSolid volume", rust_solid.mass.volume, analytic_volume),
        ("SolidCore volume", properties.volume, analytic_volume),
        ("RustSolid area", rust_solid.mass.surface_area, analytic_area),
        ("SolidCore area", properties.area, analytic_area),
    ] {
        nearly_equal(label, actual, expected, expected.abs())?;
    }

    Ok(ParityReport {
        name: case.name,
        volume: (rust_solid.mass.volume, properties.volume),
        area: (rust_solid.mass.surface_area, properties.area),
        topological_counts: rust_counts,
    })
}

fn run() -> Result<(), String> {
    for case in &CASES {
        let report = compare_box(case)?;
        println!(
            "PASS {} | F/E/V {:?} | volume RustSolid={:.12} SolidCore={:.12} | area RustSolid={:.12} SolidCore={:.12}",
            report.name, report.topological_counts,
            report.volume.0, report.volume.1, report.area.0, report.area.1
        );
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("FAIL {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_primitive_parity() {
        run().unwrap();
    }

    #[test]
    fn rejects_invalid_source_before_mutating_either_backend() {
        let invalid = BoxCase { name: "invalid", origin: [0.0; 3], size: [1.0, 0.0, 1.0] };
        assert!(compare_box(&invalid).is_err());
    }
}
