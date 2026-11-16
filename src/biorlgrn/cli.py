"""Command-line interface: import, train baseline/greedy, evaluate."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from .demo import make_demo
from .evaluate import score_predictions
from .io import prepare_dataset, read_prepared
from .ridge import GraphRidge
from .search import correlation_rank, greedy_rank


def run_experiment(data_dir, output_dir, seed=42, max_edges=20, alpha=1.0,
                   sparsity_penalty=0., holdout=0.25):
    if not 0.05 <= holdout <= 0.5:
        raise ValueError("holdout must be in [0.05,0.5]")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    X, tfs = read_prepared(data_dir)
    indices = np.arange(len(X))
    fit_idx, score_idx = train_test_split(indices, test_size=holdout,
                                          random_state=seed, shuffle=True)
    fit = X.iloc[fit_idx]
    scoring = X.iloc[score_idx]
    # All fitting data transformations remain fit-only.
    engine = GraphRidge(fit, scoring, alpha=alpha)
    for method in ("pearson", "spearman"):
        correlation_rank(fit, tfs, method=method).to_csv(out / f"ranked_{method}.csv", index=False)
    edges, path = greedy_rank(engine, tfs, max_edges=max_edges,
                              sparsity_penalty=sparsity_penalty)
    edges.to_csv(out / "ranked_greedy.csv", index=False)
    path.to_csv(out / "greedy_trajectory.csv", index=False)
    m = {"seed": seed, "n_genes": len(X.columns), "n_samples": len(X),
         "n_tfs": len(tfs), "fit_samples": len(fit_idx), "score_samples": len(score_idx),
         "ridge_alpha": alpha, "sparsity_penalty": sparsity_penalty,
         "max_edges": max_edges, "selected_edges": len(path),
         "stage": "E1 prototype: static + greedy; no RL or adaptive reliability yet",
         "ground_truth_accessed_during_fit": False}
    (out / "run_manifest.json").write_text(json.dumps(m, indent=2))
    return m


def main(argv=None):
    parser = argparse.ArgumentParser(prog="biorlgrn", description="Bio-RL-GRN benchmark prototype")
    sub = parser.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo", help="Generate independent toy smoke-test data")
    demo.add_argument("--output", default="data/demo")
    demo.add_argument("--seed", type=int, default=42)
    prep = sub.add_parser("prepare", help="Convert BEELINE run to standardized format")
    prep.add_argument("--expression", required=True)
    prep.add_argument("--ground-truth", default=None)
    prep.add_argument("--tf-list", default=None)
    prep.add_argument("--gene-axis", choices=["rows", "columns"], default="rows")
    prep.add_argument("--output", required=True)
    run = sub.add_parser("run", help="Train baselines/greedy WITHOUT reading gold edges")
    run.add_argument("--data-dir", required=True)
    run.add_argument("--output-dir", required=True)
    run.add_argument("--seed", type=int, default=42)
    run.add_argument("--max-edges", type=int, default=20)
    run.add_argument("--ridge-alpha", type=float, default=1.0)
    run.add_argument("--sparsity-penalty", type=float, default=0.0)
    run.add_argument("--holdout", type=float, default=0.25)
    ev = sub.add_parser("evaluate", help="Evaluate completed rankings against gold edges")
    ev.add_argument("--data-dir", required=True)
    ev.add_argument("--results-dir", required=True)
    ev.add_argument("--k", type=int, default=20)
    args = parser.parse_args(argv)
    if args.command == "demo":
        print(f"Created toy smoke test: {make_demo(args.output, args.seed)}")
    elif args.command == "prepare":
        info = prepare_dataset(args.expression, args.output, args.ground_truth,
                               args.tf_list, args.gene_axis)
        print(json.dumps(info, indent=2))
    elif args.command == "run":
        print(json.dumps(run_experiment(args.data_dir, args.output_dir,
              args.seed, args.max_edges, args.ridge_alpha,
              args.sparsity_penalty, args.holdout), indent=2))
    elif args.command == "evaluate":
        out = Path(args.results_dir)
        rows = []
        for name in ("pearson", "spearman", "greedy"):
            path = out / f"ranked_{name}.csv"
            if path.exists():
                rows.append({"method": name, **score_predictions(args.data_dir, path, k=args.k)})
        if not rows:
            raise FileNotFoundError("No ranked_*.csv files to evaluate")
        pd.DataFrame(rows).to_csv(out / "metrics.csv", index=False)
        print(pd.DataFrame(rows).to_string(index=False))

if __name__ == "__main__":
    main()
