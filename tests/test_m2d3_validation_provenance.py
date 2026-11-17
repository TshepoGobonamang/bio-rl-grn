"""M2D3 protects against false independence and implicit negatives."""
import pandas as pd
import pytest
from scripts.audit_m2d3_validation_provenance import audit


def _manifest(**changes):
    data = {
        "candidate_universe_expression_only": True,
        "predictions_frozen_before_assay_access": True,
        "distinct_assay_origin_documented": True,
        "assay_species": "Homo sapiens",
        "assay_type": "perturb-seq",
    }
    data.update(changes)
    return data


def _prior():
    return pd.DataFrame({"source": ["A", "C"], "target": ["B", "D"],
                         "evidence_accessions": ["GSE1; PMID:10", pd.NA]})


def _assay():
    return pd.DataFrame({
        "source": ["A", "C", "E", "F"], "target": ["B", "D", "G", "H"],
        "label": [1, 0, 1, 0],
        "assay_accessions": ["GSE1", "GSE2", "GSE3", pd.NA],
    })


def test_fail_closed_statuses_and_explicit_labels():
    table, report = audit(_prior(), _assay(), _manifest())
    assert list(table.provenance_status) == [
        "direct_accession_overlap", "unknown_prior_accession",
        "pending_global_prior_provenance_review", "unknown_assay_accession",
    ]
    assert not table.eligible_for_independent_claim.any()
    assert report["n_explicit_negative"] == 2
    assert not report["independent_validation_certified"]


@pytest.mark.parametrize("change", [
    {"candidate_universe_expression_only": False},
    {"predictions_frozen_before_assay_access": False},
    {"distinct_assay_origin_documented": False},
    {"assay_type": ""},
])
def test_disallow_unfrozen_or_prior_selected_evaluations(change):
    with pytest.raises(ValueError):
        audit(_prior(), _assay(), _manifest(**change))


def test_no_inferred_negatives_or_duplicate_conflicts():
    assay = _assay()
    assay.loc[0, "label"] = -1
    with pytest.raises(ValueError, match="explicit binary"):
        audit(_prior(), assay, _manifest())
    assay = pd.concat([_assay(), _assay().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate assay"):
        audit(_prior(), assay, _manifest())


def test_distinct_accessions_still_not_certified():
    prior = pd.DataFrame({"source": ["A"], "target": ["B"],
                          "evidence_accessions": ["GSE4"]})
    assay = pd.DataFrame({"source": ["A"], "target": ["B"],
                          "label": [1], "assay_accessions": ["GSE9"]})
    table, result = audit(prior, assay, _manifest())
    assert table.provenance_status.iloc[0] == "pending_source_level_review"
    assert result["independent_validation_certified"] is False
