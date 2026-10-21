#!/usr/bin/env python3
"""Fail CI when a kernel crate violates the one-way RustSolid dependency graph.

Uses Cargo's own metadata instead of parsing source imports or assuming
filesystem names are package names. Update ALLOWED only after an explicit
architecture review when adding a workspace member.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALLOWED: dict[str, set[str]] = {
    "rustsolid-foundation": set(),
    "rustsolid-predicates": {"rustsolid-foundation"},
    "rustsolid-geometry": {"rustsolid-foundation"},
    "rustsolid-source": {"rustsolid-foundation", "rustsolid-geometry"},
    "rustsolid-validation": {"rustsolid-foundation", "rustsolid-geometry", "rustsolid-source", "rustsolid-topology"},
    "rustsolid-boolean": {"rustsolid-foundation","rustsolid-geometry","rustsolid-source","rustsolid-topology","rustsolid-validation"},
    "rustsolid-document": {"rustsolid-boolean", "rustsolid-source", "rustsolid-topology", "rustsolid-foundation", "rustsolid-geometry"},
    "rustsolid-uv": {"rustsolid-geometry"},
    "rustsolid-topology": {"rustsolid-foundation", "rustsolid-geometry"},
    "rustsolid-intersect": {
        "rustsolid-foundation", "rustsolid-predicates", "rustsolid-geometry",
    },
    # Transitional monolith. Only this crate may depend on all lower levels.
    # Retire it as those subsystems become independent crates.
    "rustsolid": {
        "rustsolid-foundation", "rustsolid-predicates",
        "rustsolid-geometry", "rustsolid-intersect", "rustsolid-uv",
        "rustsolid-topology", "rustsolid-validation", "rustsolid-source",
    },
}


def main() -> int:
    try:
        raw = subprocess.check_output(
            ["cargo", "metadata", "--format-version=1", "--no-deps"],
            cwd=ROOT, text=True, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"ERROR: cargo metadata failed: {e}", file=sys.stderr)
        return 2
    metadata = json.loads(raw)
    member_ids = set(metadata["workspace_members"])
    pkgs = [p for p in metadata["packages"] if p["id"] in member_ids]
    by_name = {p["name"]: p for p in pkgs}
    present = set(by_name)
    if len(pkgs) != len(by_name):
        print("ERROR: duplicate workspace package names", file=sys.stderr)
        return 1
    if present != set(ALLOWED):
        print(f"ERROR: workspace drift: unregistered={sorted(present-set(ALLOWED))}, "
              f"absent={sorted(set(ALLOWED)-present)}", file=sys.stderr)
        return 1
    errors = []
    for name, pkg in sorted(by_name.items()):
        deps = set()
        for dep in pkg["dependencies"]:
            dep_name = dep["name"]
            if dep_name in present:
                deps.add(dep_name)
            elif dep_name.startswith("rustsolid-"):
                errors.append(f"{name}: unknown external RustSolid dependency {dep_name}")
        invalid = deps - ALLOWED[name]
        if invalid:
            errors.append(f"{name}: forbidden workspace dependencies {sorted(invalid)}")
        print(f"{name}: {', '.join(sorted(deps)) or '(leaf)'}")
    if errors:
        for error in errors:
            print("ERROR:", error, file=sys.stderr)
        return 1
    print(f"PASS: {len(pkgs)} workspace packages obey kernel layer boundaries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
