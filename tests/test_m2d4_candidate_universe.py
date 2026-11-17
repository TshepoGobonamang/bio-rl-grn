"""Tests for candidate freeze and non-use of perturbation/prior labels."""
import hashlib
import numpy as np
import pandas as pd
import pytest

from scripts.prepare_m2d4_candidate_universe import candidate_universe


def _expression():
    return pd.DataFrame({
        "TF_A": [0., 2., 1., 4., 3.],
        "TF_B": [2., 1., 5., 3., 0.],
        "GENE_A": [1., 2., 3., 4., 5.],
        "GENE_B": [5., 4., 3., 2., 1.],
        "CONSTANT": [2.] * 5,
    }, index=["c1", "c2", "c3", "c4", "c5"])


def test_only_expression_and_catalog_define_universe():
    a, am = candidate_universe(_expression(), ["TF_A", "TF_B"], 4, 2)
    b, bm = candidate_universe(_expression().sample(frac=1, random_state=91),
                               ["TF_B", "TF_A"], 4, 2)
    pd.testing.assert_frame_equal(a, b)
    assert am["candidate_sha256"] == bm["candidate_sha256"]
    assert set(a.source) == {"TF_A", "TF_B"}
    assert len(a) == 6
    assert (a.source != a.target).all()
    assert not am["external_labels_accessed"]


def test_nonfinite_and_insufficient_catalog_rejected():
    bad = _expression()
    bad.loc["c1", "TF_A"] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        candidate_universe(bad, ["TF_A", "TF_B"], 4, 2)
    with pytest.raises(ValueError, match="Insufficient"):
        candidate_universe(_expression(), ["TF_A", "MISSING"], 4, 2)


def test_duplicate_gene_and_repeated_catalog_rejected():
    bad = _expression()
    bad.columns = ["TF_A", "TF_A", "GENE_A", "GENE_B", "CONSTANT"]
    with pytest.raises(ValueError, match="unique"):
        candidate_universe(bad, ["TF_A", "TF_B"], 4, 2)
    with pytest.raises(ValueError, match="unique"):
        candidate_universe(_expression(), ["TF_A", "TF_A"], 4, 2)
