"""Tests for branch-aware pseudotime scoring, without reference-label access."""
import numpy as np
import pandas as pd
import pytest
from biorlgrn.pseudotime import load_pseudotime, temporal_margin, temporal_greedy
from biorlgrn.ridge import GraphRidge
from biorlgrn.demo import make_demo
from scripts import run_m12


def test_beeline_branch_mask_and_alignment(tmp_path):
    path = tmp_path / "PseudoTime.csv"
    pd.DataFrame({"PseudoTime1": [0., 1., None],
                  "PseudoTime2": [None, None, 4.]},
                 index=["a", "b", "c"]).to_csv(path)
    values = load_pseudotime(path, pd.Index(["c", "a", "b"]))
    assert values.index.tolist() == ["c", "a", "b"]
    assert np.isnan(values.loc["c", "PseudoTime1"])
    with pytest.raises(ValueError):
        load_pseudotime(path, pd.Index(["a", "missing"]))


def test_directional_scores_antisymmetric():
    rng = np.random.default_rng(3)
    n = 300
    t = np.linspace(0, 1, n)
    x = pd.DataFrame({
        "G1": np.sin(t * 9) + rng.normal(scale=.15, size=n),
        "G2": np.sin(np.maximum(t - .08, 0) * 9) + rng.normal(scale=.15, size=n),
    }, index=[f"cell-{i}" for i in range(n)])
    pt = pd.DataFrame({"PseudoTime1": t}, index=x.index)
    margins = temporal_margin(x, pt, ["G1", "G2"], bins=16)
    assert len(margins) == 2
    assert np.isclose(margins.score.iloc[0], -margins.score.iloc[1])
    eng = GraphRidge(x.iloc[:200], x.iloc[200:], alpha=1.)
    ranked, path = temporal_greedy(eng, ["G1", "G2"], margins, max_edges=1)
    assert len(ranked) == 2
    assert len(path) <= 1


def test_pseudotime_must_match_fit_rows():
    x = pd.DataFrame({"G1": range(50), "G2": range(50)})
    pt = pd.DataFrame({"PseudoTime1": range(50)}, index=range(50)[::-1])
    with pytest.raises(ValueError):
        temporal_margin(x, pt, ["G1", "G2"])


def test_no_gold_until_all_inference_outputs_exist(tmp_path, monkeypatch):
    toy = make_demo(tmp_path / "toy")
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.read_csv(toy / "expression.csv", index_col=0).T.to_csv(raw / "ExpressionData.csv")
    pd.read_csv(toy / "gold_edges.csv").rename(
        columns={"source": "Gene1", "target": "Gene2"}
    ).to_csv(raw / "GroundTruthNetwork.csv", index=False)
    x = pd.read_csv(toy / "expression.csv", index_col=0)
    pd.DataFrame({"PseudoTime1": np.arange(len(x))},
                 index=x.index).to_csv(raw / "PseudoTime.csv")
    out = tmp_path / "out"
    original = run_m12.audit_reference

    def assert_saved(*args, **kwargs):
        assert all((out / f"ranked_{name}.csv").exists() for name in [
            "greedy", "temporal_greedy", "shuffled_pseudotime", "temporal_only"
        ])
        return original(*args, **kwargs)
    monkeypatch.setattr(run_m12, "audit_reference", assert_saved)
    summary = run_m12.execute(raw, tmp_path / "prepared", out, bins=8, max_edges=4)
    assert summary["split"]["untouched_test_samples"] == 20
    assert len(pd.read_csv(out / "m12_metrics.csv")) == 4
