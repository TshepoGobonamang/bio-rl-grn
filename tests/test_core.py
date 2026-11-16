from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from biorlgrn.demo import make_demo
from biorlgrn.io import load_expression, prepare_dataset, read_prepared
from biorlgrn.ridge import GraphRidge
from biorlgrn.search import universe, greedy_rank, correlation_rank
from biorlgrn.evaluate import score_predictions
from biorlgrn.cli import run_experiment


def test_beeline_orientation_and_import(tmp_path):
    raw = pd.DataFrame({"cell1": [1, 2, 3], "cell2": [3, 4, 5]}, index=["TF1", "G2", "G3"])
    source = tmp_path / "ExpressionData.csv"
    raw.to_csv(source)
    X = load_expression(source, gene_axis="rows")
    assert X.shape == (2, 3)
    assert list(X.columns) == ["TF1", "G2", "G3"]
    gold = pd.DataFrame({"Gene1": ["TF1", "G2"], "Gene2": ["G2", "G3"]})
    gp = tmp_path / "GroundTruthNetwork.csv"
    gold.to_csv(gp, index=False)
    out = tmp_path / "prepared"
    manifest = prepare_dataset(source, out, gp)
    assert "all measured genes" in manifest["tf_universe_rule"]
    normalized = pd.read_csv(out / "expression.csv", index_col=0)
    assert normalized.shape == (2, 3)
    assert len(pd.read_csv(out / "gold_edges.csv")) == 2


def test_conditional_gain_changes_with_graph(tmp_path):
    make_demo(tmp_path / "demo")
    X, _ = read_prepared(tmp_path / "demo")
    engine = GraphRidge(X.iloc[:75], X.iloc[75:], alpha=0.1)
    empty = engine.gain(set(), ("G1", "G3"))
    changed = engine.gain({("G3", "G5")}, ("G1", "G5"))
    alone = engine.gain(set(), ("G1", "G5"))
    assert np.isfinite(empty)
    assert abs(changed - alone) > 1e-8


def test_no_gold_needed_for_training(tmp_path):
    p = make_demo(tmp_path / "demo")
    (p / "gold_edges.csv").unlink()
    out = tmp_path / "results"
    manifest = run_experiment(p, out, max_edges=4)
    assert manifest["ground_truth_accessed_during_fit"] is False
    assert (out / "ranked_greedy.csv").exists()
    with pytest.raises(FileNotFoundError):
        score_predictions(p, out / "ranked_greedy.csv")


def test_metrics_and_rank_coverage(tmp_path):
    p = make_demo(tmp_path / "demo")
    out = tmp_path / "results"
    run_experiment(p, out, max_edges=4)
    X, tf = read_prepared(p)
    n = len(universe(list(X.columns), tf))
    for method in ["pearson", "spearman", "greedy"]:
        scores = pd.read_csv(out / f"ranked_{method}.csv")
        assert len(scores) == n
        m = score_predictions(p, out / f"ranked_{method}.csv", k=3)
        assert 0 <= m["auprc"] <= 1
        assert 0 <= m["auroc"] <= 1
