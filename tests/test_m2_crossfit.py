"""M2A tests: masks, disjoint folds, reproducibility and gold-free construction."""
import json
import numpy as np
import pandas as pd
import pytest
from biorlgrn.m11 import split_three
from biorlgrn.m2_crossfit import (
    CrossfitConfig, make_reliability_meta_dataset, _source_margin,
)
from scripts.build_m2_crossfit import build

def _dummy(n=240, seed=18):
    rng = np.random.default_rng(seed)
    x = pd.DataFrame(rng.normal(size=(n, 3)), columns=["A", "B", "C"],
                     index=[f"c{i}" for i in range(n)])
    x["B"] = x["A"].shift(1).fillna(0) + rng.normal(size=n) * .1
    pt = pd.DataFrame({"PseudoTime1": np.arange(n, dtype=float)}, index=x.index)
    return x, pt

def test_crossfit_manifest_and_schema():
    x, pt = _dummy()
    fit, score, test, _ = split_three(x)
    cfg = CrossfitConfig(folds=3, bootstrap_repeats=2)
    out, m = make_reliability_meta_dataset(fit, pt.loc[fit.index], list(fit), cfg, "toy")
    assert len(out) == 3 * 6
    assert m["uses_gold_reference"] is False
    assert m["n_temporal_labels_available"] > 0
    assert out.fold_id.nunique() == 3
    assert out.y_temporal_reproducibility_TARGET.between(0, 1).all()
    assert all(f["disjoint"] for f in m["folds"])

def test_outer_test_perturbation_has_no_effect():
    x, pt = _dummy()
    fit, score, test, _ = split_three(x)
    cfg = CrossfitConfig(bootstrap_repeats=0)
    a, _ = make_reliability_meta_dataset(fit, pt.loc[fit.index], list(fit), cfg)
    altered = x.copy()
    altered.loc[test.index, "A"] = 1e9
    altered.loc[score.index, "B"] = -1e9
    b, _ = make_reliability_meta_dataset(altered.loc[fit.index],
                                          pt.loc[fit.index], list(fit), cfg)
    pd.testing.assert_frame_equal(a, b)

def test_missing_branch_is_masked_not_negative():
    x, pt = _dummy(100)
    pt.iloc[:, :] = np.nan
    data, meta = _source_margin(x, pt, list(x), CrossfitConfig())
    assert (data.m_present == 0).all()
    assert data.z_signed.isna().all()
    assert meta["source_present"] is False

def test_no_gold_required(tmp_path):
    x, pt = _dummy()
    raw = tmp_path / "inputs"
    raw.mkdir()
    x.T.to_csv(raw / "ExpressionData.csv")
    pt.to_csv(raw / "PseudoTime.csv")
    assert not (raw / "GroundTruthNetwork.csv").exists()
    output = tmp_path / "out"
    meta = build(raw, output, folds=3, bootstrap_repeats=0)
    assert meta["uses_gold_reference"] is False
    assert meta["n_rows"] == 18
    assert not {"reference_label", "gold_label"}.intersection(
        pd.read_csv(output / "m2_crossfit_edges.csv").columns
    )
    assert json.loads((output / "m2_crossfit_manifest.json").read_text())[
        "outer_score_and_test_excluded_from_dataset"
    ]

def test_missing_validation_evidence_means_missing_target():
    x, pt = _dummy(240)
    pt["PseudoTime1"] = np.nan
    pt.iloc[:20, 0] = np.arange(20)
    f, _, _, _ = split_three(x)
    frame, meta = make_reliability_meta_dataset(
        f, pt.loc[f.index], list(f), CrossfitConfig(bootstrap_repeats=0)
    )
    assert frame.y_temporal_reproducibility_TARGET.isna().all()
    assert (frame.m_present == 0).all()
    assert meta["n_temporal_labels_available"] == 0

def test_bifurcating_branch_masks():
    x, pt = _dummy(500)
    pt["PseudoTime2"] = np.nan
    pt.iloc[:250, 0] = np.linspace(0, 1, 250)
    pt.iloc[250:, 0] = np.nan
    pt.iloc[250:, 1] = np.linspace(0, 1, 250)
    f, _, _, _ = split_three(x)
    frame, meta = make_reliability_meta_dataset(
        f, pt.loc[f.index], list(f), CrossfitConfig(bootstrap_repeats=1)
    )
    assert len(frame) == 18
    assert meta["n_temporal_labels_available"] == 18
    assert (frame.branch_n == 2).all()

def test_same_seed_same_dataset():
    x, pt = _dummy()
    f, _, _, _ = split_three(x)
    cfg = CrossfitConfig(bootstrap_repeats=2, seed=13)
    a, ma = make_reliability_meta_dataset(f, pt.loc[f.index], list(f), cfg)
    b, mb = make_reliability_meta_dataset(f, pt.loc[f.index], list(f), cfg)
    pd.testing.assert_frame_equal(a, b)
    assert ma == mb
