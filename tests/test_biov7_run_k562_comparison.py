"""Meaningful BIO-V7 boundary and action-invariance regressions."""
import json
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scripts import biov7_run_k562_comparison as driver

from scripts.biov7_run_k562_comparison import (
    Inputs, METHODS, TRANSFORMS, VARIANTS, invariant_action_evidence,
    load_protocol, paired_cell_interval, select_actions,
)


def test_heldout_test_rejected_before_matrix_read():
    # Deliberately unusable matrices prove the partition guard runs first.
    inputs = Inputs(None, None, [], None, [], [], [], [], {}, {})
    with pytest.raises(ValueError, match="prohibits test"):
        inputs.frame("untouched_test", "log1p_original")


def test_constant_source_does_not_choose_a_regulator_or_mutate_input():
    contexts = pd.DataFrame({
        "tf": ["Z", "A", "Z", "A"], "target": ["T1", "T1", "T2", "T2"],
        "graph_parents": ["P|Q"] * 4, "eligible_action": [True] * 4,
        "gain": [0.] * 4, "expr_rank": [.7, .7, .5, .6],
    })
    edges = list(zip(contexts.tf, contexts.target))
    source = np.array([.9, .9, 0., 0.])
    evidence = {"dnase": (source, np.ones(4, bool))}
    stable = invariant_action_evidence(evidence, edges)
    np.testing.assert_array_equal(stable["dnase"][0], np.zeros(4))
    np.testing.assert_array_equal(source, [.9, .9, 0., 0.])
    scores, selections = select_actions(contexts, evidence, pd.DataFrame(),
                                        ["expression_only", "fixed_sources"], .3)
    by_method = selections.pivot(index="target", columns="method", values="selected_tf")
    assert by_method.expression_only.tolist() == ["A", "A"]
    assert by_method.expression_only.equals(by_method.fixed_sources)
    # Raw cross-target scores retain their offset; decision scores remove it.
    fixed = scores.query("method == 'fixed_sources'")
    np.testing.assert_allclose(fixed.score - fixed.action_score, [.27, .27, 0., 0.])


def test_observed_zero_and_missing_are_not_implicitly_the_same_source():
    evidence = {"source": (np.array([.9, .9]), np.array([True, False]))}
    stable = invariant_action_evidence(evidence, [("A", "T"), ("B", "T")])
    np.testing.assert_array_equal(stable["source"][0], [.9, .9])
    np.testing.assert_array_equal(stable["source"][1], [True, False])


def test_paired_cell_bootstrap_preserves_shared_cell_pairing_and_sign():
    right = np.arange(131 * 3, dtype=float).reshape(131, 3)
    left = right - .25
    delta, low, high = paired_cell_interval(left, right, seed=12, n_resamples=100)
    assert delta == low == high == -.25
    with pytest.raises(ValueError, match="align"):
        paired_cell_interval(left[:-1], right[:-1], seed=12, n_resamples=100)


def test_protocol_forbids_test_driven_setting_change(tmp_path):
    protocol = {"protocol": "TEST", "seed": 20261127, "alpha": 1.,
                "max_parents": 2, "boost": .3, "n_folds": 3, "shrinkage": .05,
                "bootstrap_resamples": 10000, "transforms": TRANSFORMS,
                "variants": VARIANTS, "methods": METHODS}
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(protocol))
    assert load_protocol(path)["seed"] == 20261127
    protocol["boost"] = .8
    path.write_text(json.dumps(protocol))
    with pytest.raises(ValueError, match="boost"):
        load_protocol(path)


@pytest.mark.parametrize("tamper", ["decision", "implementation", "protocol"])
def test_evaluation_rejects_tampering_before_loading_expression(tmp_path, monkeypatch, tamper):
    output = tmp_path / "frozen"
    output.mkdir()
    decision = output / "decision.csv"
    decision.write_text("tf,target\nTF,T\n")
    protocol = tmp_path / "protocol.json"
    locked = Path(__file__).parents[1] / "configs" / "biov7_k562_comparison_2026_11_27.json"
    protocol.write_bytes(locked.read_bytes())
    code_hashes = {"synthetic_module.py": "original-code-hash"}
    manifest = {
        "stage": "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE",
        "protocol_sha256": driver.sha(protocol),
        "implementation_sha256": code_hashes,
        "test_used": False, "internal_score_cells_used": 0,
        "output_sha256": {"decision.csv": driver.sha(decision)},
    }
    (output / "freeze_manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(driver, "implementation_hashes", lambda: code_hashes)

    def prohibited_load(*args):
        raise AssertionError("Expression must not be loaded after a failed freeze integrity check")

    monkeypatch.setattr(driver, "load_inputs", prohibited_load)
    if tamper == "decision":
        decision.write_text("tf,target\nDIFFERENT,T\n")
    elif tamper == "implementation":
        monkeypatch.setattr(driver, "implementation_hashes",
                            lambda: {"synthetic_module.py": "changed-code-hash"})
    else:
        # Even a scientifically harmless serialization change invalidates the
        # exact protocol bytes associated with already-frozen decisions.
        protocol.write_bytes(protocol.read_bytes() + b"\n")
    args = Namespace(output=output, protocol=protocol,
                     cohort=tmp_path / "unused", original=tmp_path / "unused",
                     biov5=tmp_path / "unused", biov6=tmp_path / "unused")
    with pytest.raises(ValueError, match="[Ff]rozen|protocol"):
        driver.evaluate(args)


def test_normalization_uses_full_transcriptome_totals_only_for_requested_rows():
    class RestrictedRows:
        def __init__(self, data):
            self.data = sparse.csr_matrix(data)
            self.allowed = {0}
            self.accessed = []

        def __getitem__(self, key):
            rows, columns = key
            row_ids = np.asarray(rows).ravel().tolist()
            assert set(row_ids).issubset(self.allowed), "Unrequested expression rows accessed"
            self.accessed.append(tuple(row_ids))
            return self.data[rows, columns]

    # Poison the test row; only exact requested fit/score rows may be indexed.
    matrix = RestrictedRows([[2., 4.], [4., 2.], [np.nan, np.inf]])
    original = RestrictedRows([[2., 4., 3., 1.], [4., 2., 2., 12.],
                               [np.nan, np.inf, np.nan, np.nan]])
    partitions = pd.DataFrame({"cell_id": ["fit", "score", "test"],
                               "partition": ["fit", "internal_score", "untouched_test"]})
    data = Inputs(matrix, original, [0, 1], partitions, ["A", "B"],
                  ["A"], ["B"], [], {}, {})
    fit = data.frame("fit", "log1p_cp10k")
    np.testing.assert_allclose(fit.to_numpy(), np.log1p([[2000., 4000.]]))
    assert fit.index.tolist() == ["fit"]
    matrix.allowed = original.allowed = {1}
    score = data.frame("internal_score", "log1p_cp10k")
    np.testing.assert_allclose(score.to_numpy(), np.log1p([[2000., 1000.]]))
    assert score.index.tolist() == ["score"]
    assert set(matrix.accessed) == {(0,), (1,)}
    assert set(original.accessed) == {(0,), (1,)}


def test_wrong_original_source_identity_fails_before_any_expression_load(tmp_path):
    cohort = tmp_path / "cohort"
    cohort.mkdir()
    (cohort / "biov4b_fit_only_manifest.json").write_text(json.dumps({
        "original_K562_H5AD_sha256": "wrong-source-fingerprint"
    }))
    with pytest.raises(ValueError, match="source identity mismatch"):
        driver.load_inputs(cohort, tmp_path / "original", tmp_path / "biov5", tmp_path / "biov6")
