"""M2C-A gold-blind leakage and directional-gain regression tests."""
import json
import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import KFold
from biorlgrn.m2_crossfit import CrossfitConfig
from biorlgrn.m2_direction import (
    _directed_heldout_gain, heldout_directional_targets, make_directional_meta_dataset
)
from biorlgrn.m11 import split_three
from scripts.build_m2_direction import build


def _toy(n=240, seed=52):
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 1, n)
    data = pd.DataFrame({
        "A": np.sin(5 * np.pi * t) + rng.normal(scale=.10, size=n),
        "B": np.sin(5 * np.pi * (t - .07)) + rng.normal(scale=.13, size=n),
        "C": rng.normal(scale=.3, size=n) + t,
    }, index=[f"cell-{k}" for k in range(n)])
    pt = pd.DataFrame({"PseudoTime1": t}, index=data.index)
    return data, pt


def test_heldout_lagged_direction_differs_for_reverse_edges():
    rng = np.random.default_rng(91)
    def make_profile():
        a = rng.normal(size=20)
        b = np.r_[0, a[:-1]] + rng.normal(scale=.02, size=20)
        return np.column_stack([a, b])
    f, v = make_profile(), make_profile()
    forward, _, _ = _directed_heldout_gain(f, v, 0, 1, alpha=.1, lag=1)
    reverse, _, _ = _directed_heldout_gain(f, v, 1, 0, alpha=.1, lag=1)
    assert forward > reverse + .1


def test_directional_margin_antisymmetric_and_masks():
    x, pt = _toy()
    f, v = x.iloc[:150], x.iloc[150:]
    pf, pv = pt.loc[f.index], pt.loc[v.index]
    cfg = CrossfitConfig(bootstrap_repeats=0)
    target, audit = heldout_directional_targets(f, v, pf, pv, list(x), cfg)
    scores = dict(zip(zip(target.source, target.target),
                      target.y_pairwise_direction_margin_TARGET))
    assert len(target) == 6
    assert all(z["used"] for z in audit)
    for s, t in scores:
        assert np.isclose(scores[(s, t)], -scores[(t, s)])
    assert target.y_directional_lag_gain_TARGET.notna().all()
    absent = pv.copy()
    absent.loc[:, :] = np.nan
    missing, _ = heldout_directional_targets(f, v, pf, absent, list(x), cfg)
    assert missing.y_directional_lag_gain_TARGET.isna().all()
    assert (missing.m_direction_TARGET_ONLY == 0).all()


def test_feature_table_excludes_validation_targets_and_is_deterministic():
    x, pt = _toy()
    outer_fit, _, _, _ = split_three(x)
    cfg = CrossfitConfig(folds=3, bootstrap_repeats=0)
    f, t, m = make_directional_meta_dataset(
        outer_fit, pt.loc[outer_fit.index], list(x), cfg, "toy")
    f2, t2, m2 = make_directional_meta_dataset(
        outer_fit, pt.loc[outer_fit.index], list(x), cfg, "toy")
    assert len(f) == len(t) == 18
    assert set(f.columns).isdisjoint(set(m["target_columns"]))
    assert not any("_TARGET" in k for k in f.columns)
    assert m["uses_gold_reference"] is False
    assert all(z["disjoint"] for z in m["folds"])
    pd.testing.assert_frame_equal(f, f2)
    pd.testing.assert_frame_equal(t, t2)
    assert m == m2


def test_validation_fold_perturbation_does_not_change_its_fit_features():
    x, pt = _toy()
    outer_fit, _, _, _ = split_three(x)
    cfg = CrossfitConfig(folds=3, bootstrap_repeats=0)
    f1, t1, _ = make_directional_meta_dataset(
        outer_fit, pt.loc[outer_fit.index], list(x), cfg, "toy")
    split = KFold(n_splits=3, shuffle=True, random_state=cfg.seed)
    _, val_ix = next(split.split(outer_fit))
    changed = outer_fit.copy()
    changed.iloc[val_ix, 0] += 10 * np.linspace(0, 1, len(val_ix))
    f2, t2, _ = make_directional_meta_dataset(
        changed, pt.loc[changed.index], list(x), cfg, "toy")
    pd.testing.assert_frame_equal(f1[f1.fold_id.eq(0)].reset_index(drop=True),
                                  f2[f2.fold_id.eq(0)].reset_index(drop=True))
    assert not t1[t1.fold_id.eq(0)].equals(t2[t2.fold_id.eq(0)])


def test_raw_builder_requires_no_reference_network(tmp_path):
    x, pt = _toy()
    raw = tmp_path / "raw"
    raw.mkdir()
    x.T.to_csv(raw / "ExpressionData.csv")
    pt.to_csv(raw / "PseudoTime.csv")
    assert not (raw / "GroundTruthNetwork.csv").exists()
    m = build(raw, tmp_path / "out", dataset_id="toy", bootstrap_repeats=0)
    assert m["n_rows"] == 18
    assert m["no_gold_reference_file_opened"]
    assert m["n_directional_targets_available"] == 18
    features = pd.read_csv(tmp_path / "out" / "m2c_features.csv")
    targets = pd.read_csv(tmp_path / "out" / "m2c_targets.csv")
    assert len(features) == len(targets)
    assert not any("TARGET" in c for c in features)
    assert "y_directional_lag_gain_TARGET" in targets
    assert json.loads((tmp_path / "out" / "m2c_manifest.json").read_text())[
        "uses_outer_scoring_or_test_cells"] is False
