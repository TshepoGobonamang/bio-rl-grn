"""M2D3 provenance eligibility gate for an external regulatory validation assay.

Fails closed: missing accession-level provenance is NOT evidence of independence.
The audit is deliberately separate from reliability model training.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

PRIOR_REQUIRED = {"source", "target", "evidence_accessions"}
ASSAY_REQUIRED = {"source", "target", "label", "assay_accessions"}


def _ids(value: object) -> set[str]:
    """Semicolon-delimited experiment identifiers. Never guess ID equivalences."""
    if pd.isna(value) or not str(value).strip():
        return set()
    return {v.strip().upper() for v in str(value).split(";") if v.strip()}


def audit(prior: pd.DataFrame, assay: pd.DataFrame, manifest: dict) -> tuple[pd.DataFrame, dict]:
    """Return edge provenance status and fail-closed eligibility summary.

    Prior rows must supply evidence_accessions tied to the specific interaction;
    assay rows must supply independently documented assay_accessions.
    False negatives are accepted ONLY when the external assay explicitly labels 0.
    """
    missing_prior = PRIOR_REQUIRED - set(prior.columns)
    missing_assay = ASSAY_REQUIRED - set(assay.columns)
    if missing_prior or missing_assay:
        raise ValueError(f"Missing columns: prior={sorted(missing_prior)} assay={sorted(missing_assay)}")
    if not manifest.get("candidate_universe_expression_only", False):
        raise ValueError("Candidate universe must be selected without curated priors or reference labels")
    if not manifest.get("predictions_frozen_before_assay_access", False):
        raise ValueError("Predictions must be frozen before external assay access")
    if not manifest.get("distinct_assay_origin_documented", False):
        raise ValueError("Independent assay origin is undocumented")
    if not manifest.get("assay_species") or not manifest.get("assay_type"):
        raise ValueError("Assay organism and type are required")
    if prior[["source", "target"]].isna().any().any() or assay[["source", "target"]].isna().any().any():
        raise ValueError("Gene identifiers cannot be missing")
    if (prior.source == prior.target).any() or (assay.source == assay.target).any():
        raise ValueError("Self-edges forbidden in external validation")
    if not assay.label.isin([0, 1]).all():
        raise ValueError("Assay labels must be explicit binary 0/1; absent edges are NOT negatives")
    if assay.duplicated(["source", "target"]).any():
        raise ValueError("Duplicate assay edges need explicit conflict adjudication")
    if prior.duplicated(["source", "target"]).any():
        raise ValueError("Duplicate prior edges need accession aggregation before audit")

    prior_map = {(str(r.source), str(r.target)): _ids(r.evidence_accessions)
                 for r in prior.itertuples(index=False)}
    rows = []
    for r in assay.itertuples(index=False):
        edge = (str(r.source), str(r.target))
        source_ids = prior_map.get(edge, set())
        assay_ids = _ids(r.assay_accessions)
        is_prior_edge = edge in prior_map
        if not assay_ids:
            status = "unknown_assay_accession"
        elif is_prior_edge and not source_ids:
            status = "unknown_prior_accession"
        elif source_ids & assay_ids:
            status = "direct_accession_overlap"
        elif is_prior_edge:
            # Distinct accessions do not prove distinct underlying studies;
            # mandatory expert source review before treating as independent.
            status = "pending_source_level_review"
        else:
            status = "pending_global_prior_provenance_review"
        rows.append({
            "source": edge[0], "target": edge[1], "label": int(r.label),
            "prior_record_present": is_prior_edge,
            "accession_overlap_count": len(source_ids & assay_ids),
            "provenance_status": status,
            "eligible_for_independent_claim": False,
        })
    table = pd.DataFrame(rows, columns=[
        "source", "target", "label", "prior_record_present",
        "accession_overlap_count", "provenance_status",
        "eligible_for_independent_claim"
    ])
    counts = table.provenance_status.value_counts().to_dict() if len(table) else {}
    summary = {
        "stage": "M2D3_source_provenance_eligibility_audit",
        "assay_type": manifest["assay_type"],
        "assay_species": manifest["assay_species"],
        "n_assay_edges": len(table),
        "n_explicit_positive": int((table.label == 1).sum()) if len(table) else 0,
        "n_explicit_negative": int((table.label == 0).sum()) if len(table) else 0,
        "statuses": {str(k): int(v) for k, v in counts.items()},
        "independent_validation_certified": False,
        "reason": "Accession inequality alone cannot establish experiment/source independence; explicit external review required.",
    }
    return table, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", required=True, type=Path)
    parser.add_argument("--assay", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    table, summary = audit(pd.read_csv(args.prior, dtype={"evidence_accessions": "string"}),
                           pd.read_csv(args.assay, dtype={"assay_accessions": "string"}),
                           manifest)
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "m2d3_provenance_edges.csv", index=False)
    (args.output / "m2d3_provenance_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
