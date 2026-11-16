"""M1.1 controlled experiment: directional diagnostics and untouched expression test."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
from biorlgrn.io import prepare_dataset, read_prepared, sha256
from biorlgrn.m11 import (
    audit_reference, evaluate_ranked, greedy_reciprocity_ablation,
    split_three, independent_expression_loss,
)
from biorlgrn.ridge import GraphRidge
from biorlgrn.search import correlation_rank, greedy_rank


def execute(raw_dir: Path, processed_dir: Path, results_dir: Path,
            *, max_edges: int = 20, alpha: float = 1.0, seed: int = 42,
            reverse_penalty: float = .01, tf_list: Path | None = None) -> dict:
    raw_dir, processed_dir, results_dir = map(Path, (raw_dir, processed_dir, results_dir))
    results_dir.mkdir(parents=True, exist_ok=True)
    expression_path = raw_dir / "ExpressionData.csv"
    gold_path = raw_dir / "GroundTruthNetwork.csv"
    if not expression_path.is_file() or not gold_path.is_file():
        raise FileNotFoundError("Expected BEELINE ExpressionData.csv and GroundTruthNetwork.csv")
    (processed_dir / "gold_edges.csv").unlink(missing_ok=True)
    manifest = prepare_dataset(expression_path, processed_dir, None, tf_list, "rows")
    X, tfs = read_prepared(processed_dir)
    fit, internal, untouched, split = split_three(X, seed=seed)
    engine = GraphRidge(fit, internal, alpha=alpha)
    ranks = {
        method: correlation_rank(fit, tfs, method)
        for method in ("pearson", "spearman")
    }
    ranks["greedy"], path = greedy_rank(engine, tfs, max_edges=max_edges)
    ranks["reciprocal_ablation"], reverse_path = greedy_reciprocity_ablation(
        engine, tfs, max_edges=max_edges, reverse_penalty=reverse_penalty
    )
    for name, ranked in ranks.items():
        ranked.to_csv(results_dir / f"ranked_{name}.csv", index=False)
    path.to_csv(results_dir / "greedy_trajectory.csv", index=False)
    reverse_path.to_csv(results_dir / "reverse_trajectory.csv", index=False)
    chosen = set(zip(path.source, path.target))
    reverse_chosen = set(zip(reverse_path.source, reverse_path.target))
    # Test is assessed only after *both* graphs have been selected.
    # Compare all four fixed top-K graphs on exactly the same untouched cells.
    # No method is allowed to use this partition during graph selection.
    test_diagnostics = {}
    for name, ranked in ranks.items():
        if name == "greedy":
            selected = chosen
        elif name == "reciprocal_ablation":
            selected = reverse_chosen
        else:
            selected = set(zip(
                ranked.head(max_edges).source, ranked.head(max_edges).target
            ))
        test_diagnostics[name] = independent_expression_loss(
            fit, untouched, selected, alpha
        )
    # No reference gold is consulted until all ranks have been generated and saved.
    eligible, raw_audit, provenance = audit_reference(gold_path, list(X.columns), tfs)
    raw_audit.to_csv(results_dir / "reference_row_audit.csv", index=False)
    eligible.to_csv(processed_dir / "gold_edges.csv", index=False)
    gold = set(zip(eligible.source, eligible.target))
    metrics = []
    for name, ranked in ranks.items():
        m = evaluate_ranked(ranked, gold, list(X.columns), tfs, k=max_edges)
        metrics.append({
            "method": name,
            **{key: value for key, value in m.items() if key != "reciprocal_pairs"},
        })
        # Labels are attached only to evaluation artifacts (not inference rankings).
        annotated = ranked.copy()
        annotated["reference_label"] = [
            int((s, t) in gold) for s, t in zip(annotated.source, annotated.target)
        ]
        annotated["reverse_in_reference"] = [
            int((t, s) in gold) for s, t in zip(annotated.source, annotated.target)
        ]
        annotated.to_csv(
            results_dir / f"evaluation_annotated_{name}.csv", index=False
        )
    pd.DataFrame(metrics).to_csv(results_dir / "m11_metrics.csv", index=False)
    summary = {
        "stage": "M1.1 exploratory: no RL or learned biological priors",
        "dataset": str(raw_dir), "seed": seed, "n_genes": len(X.columns),
        "n_samples": len(X), "n_candidate_tfs": len(tfs),
        "candidate_rule": manifest["tf_universe_rule"],
        "max_edges": max_edges, "reciprocal_penalty": reverse_penalty,
        "split": split, "independent_expression_test": test_diagnostics,
        "reference_provenance": provenance,
        "source_expression_sha256": sha256(expression_path),
        "source_gold_sha256": sha256(gold_path),
        "caveat": (
            "Single-replicate exploratory analysis; do not tune based on test "
            "results or label-based metrics."
        ),
    }
    (results_dir / "m11_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed/m11"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/m11"))
    parser.add_argument("--tf-list", type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-edges", type=int, default=20)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--reverse-penalty", type=float, default=.01)
    args = parser.parse_args()
    print(json.dumps(execute(
        args.raw_dir, args.processed_dir, args.results_dir,
        max_edges=args.max_edges, alpha=args.alpha, seed=args.seed,
        reverse_penalty=args.reverse_penalty, tf_list=args.tf_list
    ), indent=2))


if __name__ == "__main__":
    main()
