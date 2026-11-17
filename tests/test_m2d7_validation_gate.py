import pytest
from scripts.decide_m2d7_validation_gate import decide


def test_no_artifacts_is_blocked():
    report = decide({})
    assert not report["ready_for_perturbation_inference"]
    assert "real_tfmap_schema_and_checksum_missing" in report["blocking_reasons"]
    assert "h5ad_schema_and_checksum_missing" in report["blocking_reasons"]


def test_real_tfmap_alone_cannot_pass():
    result = decide({"tfmap": {"sha256": "abc", "columns": ["barcode"], "n_rows": 100}})
    assert result["schema_inspected"]
    assert not result["ready_for_perturbation_inference"]
    assert "valid_negative_controls_unverified" in result["blocking_reasons"]


def test_complete_intake_can_pass_but_never_certifies_regulatory_truth():
    schema = {
        "tfmap": {"sha256": "abc", "columns": ["barcode"], "n_rows": 10},
        "h5ad": {"sha256": "def", "n_cells": 10},
        "true_controls_identified": True,
        "independent_replicates_identified": True,
        "source_provenance_cleared": True,
    }
    barcode = {"n_matched_cells": 10, "all_expression_cells_mapped": True}
    report = decide(schema, barcode)
    assert report["ready_for_perturbation_inference"]
    assert not report["independent_regulatory_validation_certified"]


def test_partial_join_blocks_inference():
    schema = {
        "tfmap": {"sha256": "abc", "columns": ["barcode"], "n_rows": 10},
        "h5ad": {"sha256": "def"},
        "true_controls_identified": True,
        "independent_replicates_identified": True,
        "source_provenance_cleared": True,
    }
    report = decide(schema, {"n_matched_cells": 3, "all_expression_cells_mapped": False})
    assert "incomplete_cell_barcode_coverage" in report["blocking_reasons"]
