"""M2D7B evidence-backed experimental-design approval, independent of model training.

A JSON evidence dossier must link actual checksum-pinned source observations to
control, biological replication, barcode alignment, and provenance decisions.
Flags in earlier manifests alone never certify these conditions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def evaluate(schema: dict, join: dict | None, dossier: dict | None) -> dict:
    problems = []
    tfmap = schema.get("tfmap") or {}
    h5ad = schema.get("h5ad") or {}
    if not (tfmap.get("sha256") and tfmap.get("columns") and tfmap.get("n_rows", 0) > 0):
        problems.append("source_tfmap_audit_missing")
    if not (h5ad.get("sha256") and h5ad.get("n_cells", 0) > 0 and h5ad.get("obs_columns")):
        problems.append("observed_h5ad_audit_missing")
    if not join or join.get("n_matched_cells", 0) < 1:
        problems.append("barcode_join_missing")
    elif not join.get("all_expression_cells_mapped", False):
        problems.append("incomplete_barcode_join")
    if not dossier:
        problems.append("experiment_design_dossier_missing")
    else:
        required = {
            "geo_series", "assay_protocol_url", "control_definition",
            "control_source_reference", "biological_replicate_definition",
            "replicate_source_reference", "source_independence_review_url",
            "verified_by", "verified_on", "tfmap_sha256", "h5ad_sha256",
            "barcode_join_report_sha256",
        }
        missing = sorted(k for k in required if not isinstance(dossier.get(k), str)
                         or not dossier[k].strip())
        if missing:
            problems.append("dossier_fields_missing:" + ",".join(missing))
        if dossier.get("geo_series") != "GSE216595":
            problems.append("wrong_geo_series")
        if dossier.get("tfmap_sha256") != tfmap.get("sha256"):
            problems.append("tfmap_checksum_mismatch")
        if dossier.get("h5ad_sha256") != h5ad.get("sha256"):
            problems.append("h5ad_checksum_mismatch")
        if not join or dossier.get("barcode_join_report_sha256") != join.get("report_sha256"):
            problems.append("join_report_checksum_missing_or_mismatch")
        if dossier.get("controls_comparison_validated") is not True:
            problems.append("control_comparison_not_reviewed")
        if dossier.get("replicates_independent_validated") is not True:
            problems.append("biological_replicate_review_incomplete")
        if dossier.get("source_overlap_review_cleared") is not True:
            problems.append("source_provenance_review_incomplete")
    return {
        "milestone": "M2C-D7B",
        "schema_observed": bool(tfmap.get("sha256") and tfmap.get("columns")),
        "assay_inference_readiness": "eligible_for_review" if not problems else "blocked",
        "independent_biological_validation_certified": False,
        "blocking_reasons": problems,
        "note": "Even a fully populated dossier needs human methodological review; no automatic independent-ground-truth certification."
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schema", required=True, type=Path)
    p.add_argument("--join", type=Path)
    p.add_argument("--dossier", type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    load = lambda path: json.loads(path.read_text(encoding="utf-8")) if path else None
    result = evaluate(load(args.schema), load(args.join), load(args.dossier))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
