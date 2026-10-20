#!/usr/bin/env python3
"""Planning-only design gate for RustSolid crate extraction.

Check architectural DAG and source-file primary owners without assuming
that planned crates have already been implemented or can be built.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "architecture" / "kernel-crate-plan-v1.json"

FORBIDDEN = {
    "rustsolid-geometry": {"rustsolid-topology", "rustsolid-boolean"},
    "rustsolid-topology": {"rustsolid-boolean", "rustsolid-modeling", "rustsolid-document", "rustsolid-meshing"},
    "rustsolid-validation": {
        "rustsolid-construct", "rustsolid-assembly", "rustsolid-boolean", "rustsolid-classification",
        "rustsolid-imprint", "rustsolid-modeling", "rustsolid-document",
        "rustsolid-meshing", "rustsolid-archive",
    },
    "rustsolid-boolean": {"rustsolid-modeling", "rustsolid-document", "rustsolid-meshing"},
    "rustsolid-meshing": {"rustsolid-boolean", "rustsolid-document"},
    "gefest-geometry": {"rustsolid-qualification"},
}

def check(plan: dict) -> None:
    assert plan.get("schema") == "rustsolid.crate_blueprint.v1", "wrong blueprint schema"
    packages = plan["packages"]
    owners = plan["source_file_primary_owner"]
    assert len(packages) >= 16, "missing target packages"
    assert len(owners) == plan["source_file_count"], "source ownership ledger truncated"
    assert len(owners) >= 80, "baseline unexpectedly shrank"
    assert "rustsolid" not in packages, "transitional monolith must not be final target"
    assert packages["rustsolid-validation"]["status"] == "planned", "status discrepancy"
    assert packages["rustsolid-source"]["dependencies"] == ["rustsolid-foundation", "rustsolid-geometry"], "immutable source must be a leaf below topology"
    assert "rustsolid-source" in packages["rustsolid-validation"]["dependencies"], "oracle source dependency missing"
    assert packages["gefest-geometry"]["crate_path"] == "apps/gefest-geometry"
    assert set(FORBIDDEN).issubset(packages), "missing core boundary"
    for name, entry in packages.items():
        direct = entry["dependencies"]
        assert len(direct) == len(set(direct)), f"duplicate dependency in {name}"
        assert set(direct).issubset(packages), f"unknown package dependency in {name}"
        assert name not in direct, f"self dependency in {name}"
        assert entry["status"] in {"partial", "planned"}, f"bad status {name}"
        if entry["status"] == "partial":
            assert (ROOT / entry["crate_path"] / "Cargo.toml").is_file(), f"{name} manifest missing"
    state = {}
    reach = {}
    def closure(package: str, path: tuple[str, ...] = ()) -> set[str]:
        if state.get(package) == 1:
            raise AssertionError("cycle: " + " -> ".join((*path, package)))
        if state.get(package) == 2:
            return reach[package]
        state[package] = 1
        descendants = set()
        for dependency in packages[package]["dependencies"]:
            descendants.add(dependency)
            descendants.update(closure(dependency, (*path, package)))
        reach[package] = descendants
        state[package] = 2
        return descendants
    for package in packages:
        closure(package)
    for package, forbidden in FORBIDDEN.items():
        violation = reach[package] & forbidden
        assert not violation, f"{package} has forbidden dependency path to {sorted(violation)}"
    assert set(owners.values()).issubset(set(packages) | {"retire"}), "unknown source owner"
    assert all(p.startswith("src/") and p.endswith(".rs") for p in owners), "non-source ledger entry"
    current_paths = {
        p.relative_to(ROOT).as_posix()
        for p in (ROOT / "src").rglob("*.rs")
    }
    unmapped = sorted(current_paths - set(owners))
    assert not unmapped, f"unmapped Rust files: {unmapped}"
    for path in current_paths:
        if path.startswith("src/bin/"):
            assert owners[path] == "gefest-geometry", f"CLI ownership error: {path}"
    splits = plan["mandatory_split_notes"]
    assert len(splits) >= 20, "incomplete mandatory mixed-layer split inventory"
    assert all(p in owners and explanation.strip() for p, explanation in splits.items())
    print(
        f"PASS: {len(packages)} target packages ({sum(p['status']=='partial' for p in packages.values())} partially extracted), "
        f"{len(owners)} baseline sources, {len(current_paths)} current Rust source files, "
        f"{len(splits)} mandatory source splits; DAG acyclic"
    )
    print("Design check only. Planned crates have NOT been compiled or implemented.")

def main() -> int:
    try:
        check(json.loads(SOURCE.read_text(encoding="utf-8")))
        return 0
    except (AssertionError, ValueError, KeyError, OSError, TypeError) as exc:
        print(f"FAIL: RustSolid crate blueprint: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
