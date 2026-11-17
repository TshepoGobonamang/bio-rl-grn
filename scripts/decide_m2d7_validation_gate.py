"""M2D7: machine-readable decision gate for the real GSE216595 source audit.

This gate consumes observed audit artifacts, never makes up missing schema fields
or pretends that a GEO sample accession proves biological replication.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def decide(schema: dict, barcode: dict | None = None) -> dict:
    tfmap = schema.get("tfmap") or {}
    reasons = []
    if not tfmap.get("sha256") or not tfmap.get("columns") or tfmap.get("n_rows", 0) < 1:
        reasons.append("real_tfmap_schema_and_checksum_missing")
    if not schema.get("h5ad") or not (schema["h5ad"] or {}).get("sha256"):
        reasons.append("h5ad_schema_and_checksum_missing")
    if barcode is None or barcode.get("n_matched_cells", 0) < 1:
        reasons.append("no_verified_real_barcode_join")
    elif not barcode.get("all_expression_cells_mapped", False):
        reasons.append("incomplete_cell_barcode_coverage")
    if not schema.get("true_controls_identified", False):
        reasons.append("valid_negative_controls_unverified")
    if not schema.get("independent_replicates_identified", False):
        reasons.append("biologically_independent_replicates_unverified")
    if not schema.get("source_provenance_cleared", False):
        reasons.append("source_independence_not_certified")
    return {
        "milestone": "M2C-D7", "candidate_dataset": "GSE216595",
        "n_observed_tfmap_rows": tfmap.get("n_rows"),
        "n_observed_h5ad_cells": (schema.get("h5ad") or {}).get("n_cells"),
        "n_joined_cells": (barcode or {}).get("n_matched_cells"),
        "blocking_reasons": reasons,
        "schema_inspected": bool(tfmap.get("sha256") and tfmap.get("columns")),
        "ready_for_perturbation_inference": not reasons,
        "independent_regulatory_validation_certified": False,
        "note": "A passed intake gate is not an independent regulatory-accuracy result.",
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schema", type=Path, required=True)
    p.add_argument("--barcode-audit", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    schema = json.loads(a.schema.read_text(encoding="utf-8"))
    barcode = (json.loads(a.barcode_audit.read_text(encoding="utf-8"))
               if a.barcode_audit else None)
    result = decide(schema, barcode)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
