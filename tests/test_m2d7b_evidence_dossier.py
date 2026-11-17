"""M2D7B source-evidence gate regression tests."""
from scripts.decide_m2d7b_evidence_dossier import evaluate


def schema():
    return {"tfmap": {"sha256": "a", "columns": ["barcode"], "n_rows": 12},
            "h5ad": {"sha256": "b", "n_cells": 12, "obs_columns": ["sample"]}}


def join():
    return {"n_matched_cells": 12, "all_expression_cells_mapped": True,
            "report_sha256": "c"}


def dossier():
    return {
        "geo_series": "GSE216595",
        "assay_protocol_url": "https://example.org/protocol",
        "control_definition": "matched mock",
        "control_source_reference": "study methods",
        "biological_replicate_definition": "independent culture",
        "replicate_source_reference": "study methods table",
        "source_independence_review_url": "https://example.org/review",
        "verified_by": "reviewer", "verified_on": "2026-11-17",
        "tfmap_sha256": "a", "h5ad_sha256": "b",
        "barcode_join_report_sha256": "c",
        "controls_comparison_validated": True,
        "replicates_independent_validated": True,
        "source_overlap_review_cleared": True,
    }


def test_missing_real_evidence_fails_closed():
    result = evaluate(schema(), join(), None)
    assert result["assay_inference_readiness"] == "blocked"
    assert "experiment_design_dossier_missing" in result["blocking_reasons"]


def test_bad_hashes_and_unverified_controls_block():
    evidence = dossier()
    evidence["tfmap_sha256"] = "wrong"
    evidence["controls_comparison_validated"] = False
    result = evaluate(schema(), join(), evidence)
    assert "tfmap_checksum_mismatch" in result["blocking_reasons"]
    assert "control_comparison_not_reviewed" in result["blocking_reasons"]


def test_verified_dossier_still_never_certifies_biological_truth():
    result = evaluate(schema(), join(), dossier())
    assert result["assay_inference_readiness"] == "eligible_for_review"
    assert result["independent_biological_validation_certified"] is False


def test_missing_observed_h5ad_blocks_despite_dossier():
    result = evaluate({"tfmap": schema()["tfmap"]}, join(), dossier())
    assert result["assay_inference_readiness"] == "blocked"
    assert "observed_h5ad_audit_missing" in result["blocking_reasons"]
