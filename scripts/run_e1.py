"""Run exploratory E1 on an extracted BEELINE replicate.

Reference gold edges are accessed only after inference and saving rankings.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
from biorlgrn.cli import run_experiment
from biorlgrn.evaluate import score_predictions
from biorlgrn.io import prepare_dataset, sha256


def execute(raw_dir: Path, processed_dir: Path, results_dir: Path,
            tf_list: Path | None = None, seed: int = 42, max_edges: int = 19,
            alpha: float = 1.0, penalty: float = 0.0, holdout: float = 0.25) -> dict:
    raw_dir = Path(raw_dir)
    expression = raw_dir / "ExpressionData.csv"
    gold = raw_dir / "GroundTruthNetwork.csv"
    if not expression.is_file() or not gold.is_file():
        raise FileNotFoundError(f"Missing BEELINE input files in {raw_dir}")
    dataset = prepare_dataset(expression, processed_dir, gold, tf_list, gene_axis="rows")
    run = run_experiment(processed_dir, results_dir, seed=seed,
                         max_edges=max_edges, alpha=alpha,
                         sparsity_penalty=penalty, holdout=holdout)
    rows = []
    # The only stage that opens the evaluation reference network.
    for method in ("pearson", "spearman", "greedy"):
        rows.append({"method": method, **score_predictions(
            processed_dir, results_dir / f"ranked_{method}.csv", k=max_edges)})
    results_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(results_dir / "metrics.csv", index=False)
    summary = {
        "status": "exploratory_E1_not_confirmatory",
        "raw_dir": str(raw_dir.resolve()),
        "raw_expression_sha256": sha256(expression),
        "raw_gold_sha256": sha256(gold),
        "prepared_manifest": dataset,
        "run_manifest": run,
        "ranked_methods": [r["method"] for r in rows],
        "caveat": "Greedy selection repeatedly queries its internal scoring split; "
                  "retain untouched expression contexts for confirmatory evaluation.",
    }
    (results_dir / "e1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--processed-dir", default="data/processed/dyn_LL_500_1", type=Path)
    parser.add_argument("--results-dir", default="results/dyn_LL_500_1", type=Path)
    parser.add_argument("--tf-list", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-edges", type=int, default=19)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--penalty", type=float, default=0.0)
    parser.add_argument("--holdout", type=float, default=0.25)
    args = parser.parse_args()
    result = execute(args.raw_dir, args.processed_dir, args.results_dir,
                     tf_list=args.tf_list, seed=args.seed, max_edges=args.max_edges,
                     alpha=args.alpha, penalty=args.penalty, holdout=args.holdout)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
