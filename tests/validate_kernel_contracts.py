#!/usr/bin/env python3
"""Architecture-only schema and synthetic transition tests for RFC 0003.

No CAD geometry is constructed, no shape validity is asserted, and this script
does not replace G7/G8 independent mathematical and Rust integration tests.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
ARCH = ROOT / "docs" / "architecture"
EXPECTED_OWNERS = {
    "SourceSpec": "rustsolid-source",
    "BrepModel": "rustsolid-topology",
    "TopologyLineage": "rustsolid-document",
    "CertifiedCandidate": "rustsolid-boolean",
}
MANDATORY_SCENARIOS = {
    "g7_guarded_ready", "g8_guarded_ready", "formally_certified_ready",
    "same_operand_identity_difference_empty",
    "empty_with_disjointness_witness", "different_source_unit_no_transform",
    "explicit_unit_transform_allowed", "source_version_changed",
    "candidate_from_different_sources", "validation_for_different_sources",
    "source_changed_before_commit", "missing_actual_boundary_flux",
    "correct_volume_wrong_material", "proximity_only_vertex_merge",
    "periodic_seam_as_new_material_cut", "interval_claim_missing_rounding_proof",
    "ambiguous_used_reference", "ambiguous_unreferenced_new_feature",
    "foreign_document_candidate", "stale_expected_revision",
    "serialized_candidate_is_not_commit_authority", "failed_pre_commit_revalidation",
    "undo_checkpoint_unavailable", "already_committed", "empty_without_source_witness",
}

def read_json(filename: str) -> dict:
    return json.loads((ARCH / filename).read_text(encoding="utf-8"))

def checked_number(v: object) -> bool:
    return isinstance(v, (float, int)) and not isinstance(v, bool) and math.isfinite(v)

def legal_interval(v: object) -> bool:
    return (isinstance(v, list) and len(v) == 2
            and all(checked_number(n) for n in v) and v[0] <= v[1])

def patch_case(base: dict, edits: dict) -> dict:
    state = deepcopy(base)
    for path, value in edits.items():
        fields = path.split(".")
        target = state
        for field in fields[:-1]:
            if isinstance(target, list):
                assert field.isdigit(), f"non-index list path: {path}"
                target = target[int(field)]
            else:
                assert isinstance(target, dict) and field in target, f"unknown patch path {path}"
                target = target[field]
        field = fields[-1]
        if isinstance(target, list):
            assert field.isdigit() and int(field) < len(target), f"invalid list patch {path}"
            target[int(field)] = value
        else:
            assert isinstance(target, dict) and field in target, f"unknown patch key {path}"
            target[field] = value
    return state

@dataclass(frozen=True)
class Result:
    status: str
    reason: str

def evaluate_design(s: dict) -> Result:
    src, op, shape = s["source"], s["operation"], s["brep"]
    evidence, lineage, context = s["validation"], s["lineage"], s["commit"]
    def held(reason: str) -> Result:
        return Result("Held", reason)
    if not op["supported"] or shape["kind"] not in {"Empty", "Solids"} or not shape["supported_domain"]:
        return Result("Unsupported", "outside pinned supported modeling domain")
    if op["kind"] not in {"Union", "Intersection", "Difference"}:
        return Result("Unsupported", "unrecognized operation kind")
    specs = src["specs"]
    if len(specs) < 2 or len({x["id"] for x in specs}) != len(specs):
        return held("invalid source IDs")
    by_id = {x["id"]: x for x in specs}
    if len(src["operands"]) != 2:
        return held("invalid pair of operand references")
    for ref in src["operands"]:
        if ref["id"] not in by_id or ref["revision"] != by_id[ref["id"]]["revision"]:
            return held("source version or reference mismatch")
    if not all(x["finite"] and x["valid_frame"] and x["positive_dimensions"] for x in specs):
        return held("invalid source coordinates, frame or dimensions")
    if len({x["units"] for x in specs}) != 1 and not src["explicit_unit_transform"]:
        return held("implicit unit conversion")
    stamp = src["snapshot_stamp"]
    if not stamp or any(src[key] != stamp for key in ("candidate_stamp", "certificate_stamp", "current_stamp")):
        return held("snapshot/candidate/certificate/current source stamp mismatch")
    if shape["kind"] == "Empty":
        if not shape["empty_proof"]:
            return held("empty geometry without material proof")
    else:
        mandatory = [
            "closed_manifold", "oriented_shells", "valid_carrier_catalog",
            "pcurve_3d_agreement", "physical_vertex_binding",
            "seam_not_extra_cut", "inner_loops_valid", "shell_nesting_verified",
            "exact_analytic_carriers", "source_face_coverage_complete",
            "material_classification_complete",
        ]
        for flag in mandatory:
            if shape[flag] is not True:
                return held(f"missing canonical analytic B-rep proof: {flag}")
    if not evidence["certifier_owned"] or not evidence["independent_paths"]:
        return held("missing independent accepted validation authority")
    a, b = evidence["source_interval"], evidence["flux_interval"]
    if not legal_interval(a) or not legal_interval(b):
        return held("missing or invalid volume interval")
    guard = evidence["numerical_guard"]
    if not checked_number(guard) or guard < 0:
        return held("invalid numerical guard")
    if a[1] + guard < b[0] or b[1] + guard < a[0]:
        return held("source interval and actual flux incompatible")
    measured = evidence["measured_result_volume"]
    if (not checked_number(measured) or
            not (a[0] - guard <= measured <= a[1] + guard) or
            not (b[0] - guard <= measured <= b[1] + guard)):
        return held("result volume not independently enclosed")
    if shape["kind"] == "Empty" and (
        a != [0, 0] or b != [0, 0] or measured != 0
    ):
        return held("empty result must have zero certified volume")
    assurance = evidence["assurance"]
    if assurance == "GuardedF64":
        if guard <= 0 or evidence["scoped_family"] not in evidence["guarded_allowlist"]:
            return held("guarded f64 outside explicit domain acceptance policy")
    elif assurance == "CertifiedInterval":
        if not evidence["formal_interval_proof"] or guard != 0:
            return held("unsupported claim of formal interval certification")
    else:
        return held("unsupported numerical assurance tier")
    if not lineage["provenance_complete"] or not lineage["merge_source_proven"] or lineage["names_recycled"]:
        return held("invalid source or persistent topology lineage")
    if lineage["unresolved_names"] and lineage["dependent_reference_ambiguous"]:
        return held("ambiguous reference used by a dependent feature")
    if context["candidate_incarnation"] != context["current_incarnation"]:
        return held("foreign document proposal")
    if context["expected_revision"] != context["current_revision"]:
        return held("stale document revision")
    if context["candidate_deserialized"]:
        return held("serialized candidate is not a trusted runtime proof")
    if context["already_committed"]:
        return held("proposal already committed")
    if not context["revalidation_succeeded"]:
        return held("exclusive pre-commit revalidation failed")
    if not context["history_reserved"]:
        return held("cannot allocate atomic history checkpoint")
    return Result("Ready", "eligible in synthetic contract state machine")

def check_declarations(contract: dict, blueprint: dict, rfc: str) -> None:
    assert contract["schema"] == "rustsolid.kernel_contracts.v1"
    assert contract["status"] == "design_only_no_new_crates"
    assert set(contract["contracts"]) == set(EXPECTED_OWNERS)
    assert len(contract["flow"]) == 6
    assert contract["flow"][-1]["from"] == "CertifiedCandidate"
    assert contract["flow"][-1]["to"] == "DocumentRevision"
    assert rfc.count("## ") >= 7
    for name, owner in EXPECTED_OWNERS.items():
        decl = contract["contracts"][name]
        assert decl["owner"] == owner, f"{name}: wrong declared owner"
        assert name in rfc, f"{name}: missing RFC section"
        assert owner in blueprint["packages"], f"{name}: missing blueprint package"
        assert len(decl["required_fields"]) > 0
    deps = blueprint["packages"]
    for dep in contract["contracts"]["SourceSpec"]["dependencies"]:
        assert dep in deps["rustsolid-source"]["dependencies"]
    for source, forbidden in contract["invariants"]["dependency_forbidden"]:
        assert forbidden not in deps[source]["dependencies"], f"forbidden direct edge: {source}->{forbidden}"
    assert deps["rustsolid-validation"]["dependencies"].count("rustsolid-source") == 1
    assert "rustsolid-construct" not in deps["rustsolid-validation"]["dependencies"]
    support = contract["supporting_contracts"]
    assert support["LineageEvidence"]["owner"] == "rustsolid-assembly"
    assert support["ActualBoundaryFlux"]["owner"] == "rustsolid-validation"
    assert support["SourceVolumeEnclosure"]["owner"] == "rustsolid-validation"
    assert "BrepModel" in support["SourceVolumeEnclosure"]["must_not_read"]
    assert "CandidatePlan" in support["ActualBoundaryFlux"]["must_not_read"]
    assert contract["contracts"]["CertifiedCandidate"]["validation_owner"] == "rustsolid-validation"
    for f in ("Mesh","MassProperties","BrepOrigin"):
        assert f in contract["contracts"]["BrepModel"]["forbidden_values"]
    assert {"Ready","Held","Unsupported"} == set(contract["contracts"]["CertifiedCandidate"]["outcome_variants"])
    assert {"Preserved","SplitFrom","MergedFrom","GeneratedFrom","Deleted","Unresolved"} == set(contract["contracts"]["TopologyLineage"]["relation_variants"])

def check_scenarios(scenarios: dict) -> tuple[int, int]:
    assert scenarios["schema"] == "rustsolid.kernel_contract_scenarios.v1"
    assert scenarios["status"] == "synthetic_architecture_gate_not_geometric_execution"
    seen = set()
    accepted = 0
    rejected = 0
    base = scenarios["base"]
    for case in scenarios["cases"]:
        key = case["id"]
        assert key not in seen, f"duplicate scenario id {key}"
        seen.add(key)
        assert case["expected"] in {"Ready", "Held", "Unsupported"}
        modified = patch_case(base, case["patch"])
        untouched = deepcopy(modified)
        verdict = evaluate_design(modified)
        assert verdict.status == case["expected"], (
            f"{key}: expected {case['expected']}, got {verdict.status}: {verdict.reason}"
        )
        assert modified == untouched, f"{key}: design gate mutated staged source state"
        if verdict.status == "Ready":
            accepted += 1
        else:
            rejected += 1
    assert MANDATORY_SCENARIOS.issubset(seen), "essential positive/negative contract scenario missing"
    assert accepted >= 5 and rejected >= 30, "insufficient synthetic gate coverage"
    return accepted, rejected

def main() -> int:
    try:
        contract = read_json("kernel-contracts-v1.json")
        blueprint = read_json("kernel-crate-plan-v1.json")
        scenarios = read_json("kernel-contract-scenarios-v1.json")
        rfc = (ARCH / "kernel-contracts-v1.md").read_text(encoding="utf-8")
        check_declarations(contract, blueprint, rfc)
        accepted, rejected = check_scenarios(scenarios)
        print(
            f"PASS design only: {len(contract['contracts'])} cross-crate contracts, "
            f"{accepted} eligible and {rejected} rejected/unsupported synthetic scenarios."
        )
        print("No Rust CAD algorithm compiled/executed by this contract checker.")
        return 0
    except (AssertionError, KeyError, ValueError, TypeError, OSError, IndexError) as exc:
        print(f"FAIL: kernel contract gate: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
