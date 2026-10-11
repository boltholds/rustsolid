"""Validate the versioned official Parasolid capability evidence ledger.

The test checks provenance/consistency, not third-party functionality or
private source code, and never treats Siemens output as a RustSolid oracle.
"""
import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "parasolid-official-capabilities.v1.json"

RELEASE_SOURCE = {
    "2018 public capability": "data_2018",
    "v32.0": "tracking_2019",
    "v33.0": "mixed_2020",
    "2020 public capability": "qa_2020",
    "v38.0": "v38_2025",
    "v39.0": "v39_2026",
}
STATUSES = {"partial", "missing", "implemented"}
PRIORITIES = {"P0", "P1", "P2"}


def validate_ledger(doc):
    if doc.get("schema") != "rustsolid.parasolid_official_capabilities.v1":
        raise ValueError("Unknown ledger schema")
    if doc.get("kernel_baseline") != "v0.14.0":
        raise ValueError("Re-evaluate comparison when the RustSolid baseline changes")
    reference = doc.get("previous_decomp_reference", {})
    if reference.get("available_for_direct_verification") is not False:
        raise ValueError("Decompiled source availability must be independently reverified")
    sources = doc.get("official_sources", {})
    rows = doc.get("capabilities", [])
    if not isinstance(rows, list) or len(rows) < 10:
        raise ValueError("Capability ledger is incomplete")
    unique_ids = set()
    for row in rows:
        identifier = row.get("id")
        if not isinstance(identifier, str) or not identifier or identifier in unique_ids:
            raise ValueError("Missing or duplicate capability ID")
        unique_ids.add(identifier)
        release = row.get("official_release")
        if release not in RELEASE_SOURCE:
            raise ValueError(f"{identifier}: unknown or unversioned release")
        url = row.get("official_url")
        if sources.get(RELEASE_SOURCE[release]) != url:
            raise ValueError(f"{identifier}: source does not support claimed release")
        if not isinstance(url, str) or not (
            url.startswith("https://blogs.sw.siemens.com/plm-components/")
            or url.startswith("https://developer.siemens.com/")
            or url.startswith("https://www.siemens.com/")
        ):
            raise ValueError(f"{identifier}: not an official Siemens URL")
        if not isinstance(row.get("official_behavior"), str) or len(row["official_behavior"]) < 30:
            raise ValueError(f"{identifier}: missing behavior statement")
        if row.get("priority") not in PRIORITIES:
            raise ValueError(f"{identifier}: unknown engineering priority")
        rust = row.get("rustsolid_v014", {})
        if rust.get("status") not in STATUSES or not rust.get("gap") or not rust.get("implemented"):
            raise ValueError(f"{identifier}: needs an explicit independent kernel gap")
        oracle = row.get("independent_acceptance_test", "")
        if not isinstance(oracle, str) or len(oracle) < 40:
            raise ValueError(f"{identifier}: missing falsifiable acceptance check")
        if "opencascade" in oracle.lower() or "parasolid" in oracle.lower():
            raise ValueError(f"{identifier}: another CAD engine must not be the correctness oracle")
        symbols = row.get("symbols_from_previous_decomp_notes")
        if not isinstance(symbols, list) or any(not isinstance(s, str) for s in symbols):
            raise ValueError(f"{identifier}: invalid symbol observations")
        allowed = "prior_notes_only_not_rechecked" if symbols else "not_established"
        if row.get("decomp_verification") != allowed:
            raise ValueError(f"{identifier}: impossible symbol verification claim")
        if release == "v39.0" and symbols:
            raise ValueError(f"{identifier}: cannot project v39 capabilities back to v38 export names")
    if "v38.0" not in {r["official_release"] for r in rows} or (
        "v39.0" not in {r["official_release"] for r in rows}
    ):
        raise ValueError("Need explicit v38/v39 version boundary")
    return len(rows)


class OfficialParasolidCapabilityAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(LEDGER.read_text(encoding="utf-8"))

    def test_official_source_mapping_and_independent_acceptance(self):
        self.assertGreaterEqual(validate_ledger(self.original), 10)

    def test_invalid_release_projection_is_rejected(self):
        doc = json.loads(json.dumps(self.original))
        candidate = next(r for r in doc["capabilities"] if r["official_release"] == "v39.0")
        candidate["official_release"] = "v38.0"
        with self.assertRaises(ValueError):
            validate_ledger(doc)

    def test_unverified_decompiler_names_cannot_be_labeled_verified(self):
        doc = json.loads(json.dumps(self.original))
        candidate = next(r for r in doc["capabilities"]
                         if r["symbols_from_previous_decomp_notes"])
        candidate["decomp_verification"] = "verified"
        with self.assertRaises(ValueError):
            validate_ledger(doc)

    def test_missing_independent_oracle_is_rejected(self):
        doc = json.loads(json.dumps(self.original))
        doc["capabilities"][0]["independent_acceptance_test"] = ""
        with self.assertRaises(ValueError):
            validate_ledger(doc)


if __name__ == "__main__":
    unittest.main()
