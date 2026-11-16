import pandas as pd
from biorlgrn.demo import make_demo
from scripts.run_e1 import execute


def test_e1_separates_inference_from_evaluation(tmp_path):
    toy = make_demo(tmp_path / "toy")
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.read_csv(toy / "expression.csv", index_col=0).T.to_csv(raw / "ExpressionData.csv")
    pd.read_csv(toy / "gold_edges.csv").rename(
        columns={"source": "Gene1", "target": "Gene2"}
    ).to_csv(raw / "GroundTruthNetwork.csv", index=False)
    summary = execute(raw, tmp_path / "prepared", tmp_path / "results", max_edges=4)
    assert summary["status"] == "exploratory_E1_not_confirmatory"
    assert summary["run_manifest"]["ground_truth_accessed_during_fit"] is False
    assert set(summary["ranked_methods"]) == {"pearson", "spearman", "greedy"}
    metrics = pd.read_csv(tmp_path / "results" / "metrics.csv")
    assert len(metrics) == 3
    assert metrics["n_candidates"].nunique() == 1
