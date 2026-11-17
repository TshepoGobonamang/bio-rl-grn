"""M2D6B exact barcode mapping guards."""
import pandas as pd
import pytest
from scripts.audit_m2d6b_barcode_join import audit_join


def check(mapping, cells):
    return audit_join(mapping, cells, map_barcode="barcode",
                      cell_barcode="cell", tf_col="construct",
                      map_sample="sample", cell_sample="library")


def test_same_barcode_different_samples_is_safe():
    mapping = pd.DataFrame({"sample":["s1","s2"],
                            "barcode":["AA-1","AA-1"],
                            "construct":["TF1","TF2"]})
    cells = pd.DataFrame({"library":["s2","s1","s3"],
                          "cell":["AA-1","AA-1","BB-1"]})
    joined, report = check(mapping, cells)
    assert joined.tf_assignment.tolist()[:2] == ["TF2", "TF1"]
    assert pd.isna(joined.tf_assignment.iloc[2])
    assert report["n_matched_cells"] == 2
    assert report["n_unmatched_expression_cells"] == 1
    assert not report["biological_replicates_verified"]


def test_conflicting_assignments_never_arbitrarily_resolved():
    mapping = pd.DataFrame({"sample":["s1","s1"],"barcode":["AA-1","AA-1"],
                            "construct":["TF1","TF2"]})
    cells = pd.DataFrame({"library":["s1"],"cell":["AA-1"]})
    with pytest.raises(ValueError, match="conflicting"):
        check(mapping, cells)


def test_missing_or_duplicate_cell_identifiers_rejected():
    mapping = pd.DataFrame({"sample":["s1"],"barcode":["AA-1"],
                            "construct":["TF1"]})
    cells = pd.DataFrame({"library":["s1","s1"],"cell":["AA-1","AA-1"]})
    with pytest.raises(ValueError, match="Duplicated"):
        check(mapping, cells)
    cells = pd.DataFrame({"library":["s1"],"cell":[None]})
    with pytest.raises(ValueError, match="missing"):
        check(mapping, cells)
