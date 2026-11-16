import pandas as pd
import pytest
from biorlgrn.demo import make_demo
from biorlgrn.m11 import audit_reference, evaluate_ranked, orientation_metrics, split_three
from scripts import run_m11


def test_three_way_split():
    X = pd.DataFrame({"g": range(100)})
    fit, score, test, meta = split_three(X, 42)
    assert (len(fit), len(score), len(test)) == (60, 20, 20)
    assert set(fit.index).isdisjoint(score.index)
    assert set(fit.index).isdisjoint(test.index)
    assert set(score.index).isdisjoint(test.index)
    assert meta["split_indices_disjoint"]


def test_reference_audit(tmp_path):
    path = tmp_path / "gold.csv"
    pd.DataFrame({"Gene1": ["A", "A", "A", "Q", "B", "A", "C"],
                  "Gene2": ["B", "B", "A", "B", "Q", "C", "B"]}).to_csv(path, index=False)
    valid, raw, counts = audit_reference(path, ["A", "B", "C"], ["A"])
    assert list(zip(valid.source, valid.target)) == [("A", "B"), ("A", "C")]
    assert counts["raw_gold_rows"] == 7
    for label in ["duplicate", "self_loop_excluded", "source_not_measured",
                  "target_not_measured", "source_not_candidate_tf"]:
        assert counts["gold_exclusions"][label] == 1


def test_orientation_and_ties():
    ranked = pd.DataFrame([("A", "B", 1.), ("B", "A", 1.),
                           ("A", "C", 2.), ("C", "A", .1)],
                          columns=["source", "target", "score"])
    m = orientation_metrics(ranked, {("A", "B"), ("A", "C")}, 4)
    assert m["top_k_reciprocal_pairs"] == 2
    assert m["direction_correct"] == 1
    assert m["direction_ties"] == 1
    assert m["direction_accuracy_ties_half"] == .75


def test_reference_isolation(tmp_path, monkeypatch):
    toy = make_demo(tmp_path / "toy")
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.read_csv(toy / "expression.csv", index_col=0).T.to_csv(raw / "ExpressionData.csv")
    pd.read_csv(toy / "gold_edges.csv").rename(columns={
        "source": "Gene1", "target": "Gene2"
    }).to_csv(raw / "GroundTruthNetwork.csv", index=False)
    prepared, results = tmp_path / "prepared", tmp_path / "results"
    original = run_m11.audit_reference

    def checked(*args, **kwargs):
        assert (results / "ranked_greedy.csv").exists()
        assert (results / "ranked_reciprocal_ablation.csv").exists()
        return original(*args, **kwargs)

    monkeypatch.setattr(run_m11, "audit_reference", checked)
    summary = run_m11.execute(raw, prepared, results, max_edges=4)
    assert summary["split"]["untouched_test_samples"] == 20
    assert summary["reference_provenance"]["gold_eligible"] == 3
    assert (results / "m11_metrics.csv").exists()
    assert (results / "reference_row_audit.csv").exists()
    assert (results / "evaluation_annotated_greedy.csv").exists()


def test_complete_candidate_universe():
    ranked = pd.DataFrame([("A", "B", 1.), ("B", "A", .5)],
                          columns=["source", "target", "score"])
    result = evaluate_ranked(ranked, {("A", "B")}, ["A", "B"], ["A", "B"], 1)
    assert result["average_precision"] == 1.
    assert result["direction_accuracy_ties_half"] == 1.
    with pytest.raises(ValueError):
        evaluate_ranked(ranked.head(1), {("A", "B")}, ["A", "B"], ["A", "B"], 1)
