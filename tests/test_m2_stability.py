"""M2C-A2: reproducibility, negative control, and gold isolation tests."""
import json
import numpy as np
import pandas as pd
import pytest

from biorlgrn.m2_stability import (
    DEVELOPMENT, CONFIRMATION, DEV_BINS, DEV_SEEDS, LOCKED_BINS, LOCKED_SEED,
    aggregate_setting, run_network, run_setting, shuffle_pseudotime, summarize_setting,
)
from scripts.evaluate_m2ca2 import evaluate_one


def _toy(n=240):
    rng = np.random.default_rng(22)
    t = np.linspace(0, 1, n)
    x = pd.DataFrame({
        "g1": np.sin(10 * t) + rng.normal(0, .15, n),
        "g2": np.sin(10 * (t - .08)) + rng.normal(0, .12, n),
        "g3": rng.normal(0, .4, n) + .4 * t,
    }, index=[f"cell{i}" for i in range(n)])
    pt = pd.DataFrame({"PseudoTime1": t}, index=x.index)
    return x, pt


def _raw(tmp_path, name):
    x, p = _toy()
    d = tmp_path / name
    d.mkdir()
    x.T.to_csv(d / "ExpressionData.csv")
    p.to_csv(d / "PseudoTime.csv")
    return d


def test_shuffling_preserves_branch_missingness_and_ids():
    x, pt = _toy()
    pt.loc[pt.index[::2], "PseudoTime2"] = np.nan
    pt.loc[pt.index[1::2], "PseudoTime2"] = np.linspace(0, 1, len(pt.index[1::2]))
    a = shuffle_pseudotime(pt, 42)
    b = shuffle_pseudotime(pt, 42)
    pd.testing.assert_frame_equal(a, b)
    assert a.index.equals(pt.index)
    assert a.isna().equals(pt.isna())
    for col in pt:
        assert np.array_equal(
            np.sort(a[col].dropna().to_numpy()),
            np.sort(pt[col].dropna().to_numpy())
        )


def test_setting_is_deterministic_and_keeps_folds_disjoint():
    x, pt = _toy()
    a, fa = run_setting(x, pt, "toy", seed=42, bins=12)
    b, fb = run_setting(x, pt, "toy", seed=42, bins=12)
    pd.testing.assert_frame_equal(a, b)
    assert fa == fb
    assert len(a) == 18
    assert a.source.ne(a.target).all()
    assert not a.duplicated(["setting", "fold_id", "source", "target"]).any()
    assert all(f["disjoint"] for f in fa)
    assert a.target_available.eq(1).all()
    assert a.gain.between(-1, 1).all()
    stat = summarize_setting(a)
    assert -1 <= stat["gain_fold_spearman"] <= 1
    assert aggregate_setting(a).shape[0] == 6


def test_shuffled_time_changes_directional_scores():
    x, pt = _toy()
    clean, _ = run_setting(x, pt, "toy", seed=42, bins=12)
    bad, _ = run_setting(x, pt, "toy", seed=42, bins=12, shuffled=True)
    assert clean[["fold_id", "source", "target"]].equals(
        bad[["fold_id", "source", "target"]])
    assert not np.allclose(clean.gain, bad.gain)
    assert bad.shuffled.all() and not clean.shuffled.any()


def test_role_whitelist_prevents_silent_confirmation_reuse(tmp_path):
    dev = _raw(tmp_path, DEVELOPMENT[0])
    with pytest.raises(ValueError, match="confirmation"):
        run_network(dev, tmp_path / "out", role="confirmation")
    candidate = _raw(tmp_path, CONFIRMATION[0])
    with pytest.raises(ValueError, match="development"):
        run_network(candidate, tmp_path / "out2", role="development")


def test_reference_is_not_needed_for_locked_score_creation(tmp_path):
    dataset = CONFIRMATION[0]
    raw = _raw(tmp_path, dataset)
    results = tmp_path / "results"
    root = results / dataset
    assert not (raw / "GroundTruthNetwork.csv").exists()
    manifest = run_network(raw, root, role="confirmation")
    assert manifest["n_rows"] == 36
    assert manifest["n_settings"] == 2
    audit = json.loads((root / "m2ca2_manifest.json").read_text())
    assert audit["gold_standard_opened"] is False
    assert audit["n_outer_fit_cells"] == 144
    with pytest.raises(FileNotFoundError):
        evaluate_one(results, tmp_path, dataset)
    pd.DataFrame({"Gene1": ["g1", "g2"], "Gene2": ["g2", "g3"]}).to_csv(
        raw / "GroundTruthNetwork.csv", index=False)
    metrics = evaluate_one(results, tmp_path, dataset)
    assert set(metrics.method) == {
        "z_signed", "gain", "direction_margin", "fixed_blend",
        "shuffled_margin"
    }
    assert metrics.average_precision.between(0, 1).all()


def test_locked_settings_and_blend_are_predeclared():
    assert set(DEV_BINS) == {8, 12, 16}
    assert set(DEV_SEEDS) == {42, 73, 137}
    assert LOCKED_BINS == 12 and LOCKED_SEED == 42
    assert set(DEVELOPMENT).isdisjoint(CONFIRMATION)
