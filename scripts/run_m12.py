"""Exploratory pseudotime-direction experiment with independent test and shuffled control."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from biorlgrn.io import prepare_dataset, read_prepared, sha256
from biorlgrn.ridge import GraphRidge
from biorlgrn.search import greedy_rank
from biorlgrn.m11 import (
    audit_reference, evaluate_ranked, split_three, independent_expression_loss
)
from biorlgrn.pseudotime import load_pseudotime, temporal_margin, temporal_greedy


def execute(raw: Path, prepared: Path, results: Path, *,
            seed: int = 42, bins: int = 16, temporal_weight: float = .015,
            max_edges: int = 20, alpha: float = 1.) -> dict:
    raw, prepared, results = Path(raw), Path(prepared), Path(results)
    results.mkdir(parents=True, exist_ok=True)
    (prepared / "gold_edges.csv").unlink(missing_ok=True)
    manifest = prepare_dataset(raw / "ExpressionData.csv", prepared, None, None, "rows")
    x, tfs = read_prepared(prepared)
    pseudo = load_pseudotime(raw / "PseudoTime.csv", x.index)
    fit, score, test, split = split_three(x, seed=seed)
    engine = GraphRidge(fit, score, alpha=alpha)
    evidence = temporal_margin(fit, pseudo.loc[fit.index], tfs, bins=bins)
    shuffle = pseudo.loc[fit.index].copy()
    rng = np.random.default_rng(seed + 830)
    for branch in shuffle:
        keep = shuffle[branch].notna()
        values = shuffle.loc[keep, branch].to_numpy(copy=True)
        rng.shuffle(values)
        shuffle.loc[keep, branch] = values
    negative = temporal_margin(fit, shuffle, tfs, bins=bins)
    baseline, base_path = greedy_rank(engine, tfs, max_edges=max_edges)
    combined, temporal_path = temporal_greedy(
        engine, tfs, evidence, max_edges=max_edges, temporal_weight=temporal_weight
    )
    shuffled_rank, shuffle_path = temporal_greedy(
        engine, tfs, negative, max_edges=max_edges, temporal_weight=temporal_weight
    )
    temporal_only = evidence.sort_values(
        ["score", "source", "target"], ascending=[False, True, True]
    ).reset_index(drop=True)
    ranks = {
        "greedy": baseline, "temporal_greedy": combined,
        "shuffled_pseudotime": shuffled_rank, "temporal_only": temporal_only
    }
    for name, ranked in ranks.items():
        ranked.to_csv(results / f"ranked_{name}.csv", index=False)
    evidence.to_csv(results / "temporal_margins.csv", index=False)
    base_path.to_csv(results / "greedy_trajectory.csv", index=False)
    temporal_path.to_csv(results / "temporal_trajectory.csv", index=False)
    shuffle_path.to_csv(results / "shuffled_trajectory.csv", index=False)
    # All graph selections are fixed before evaluating untouched expression cells.
    losses = {}
    for name, ranked in ranks.items():
        chosen = set(zip(ranked.head(max_edges).source, ranked.head(max_edges).target))
        losses[name] = independent_expression_loss(fit, test, chosen, alpha)
    # Only now read the gold-reference network (not available to search).
    gold_frame, audit, provenance = audit_reference(
        raw / "GroundTruthNetwork.csv", list(x.columns), tfs
    )
    audit.to_csv(results / "reference_audit.csv", index=False)
    gold = set(zip(gold_frame.source, gold_frame.target))
    metrics = []
    for name, ranked in ranks.items():
        scores = evaluate_ranked(ranked, gold, list(x.columns), tfs, max_edges)
        metrics.append({
            "method": name,
            **{k: v for k, v in scores.items() if k != "reciprocal_pairs"},
            **losses[name],
        })
        labeled = ranked.copy()
        labeled["reference_label"] = [
            int(e in gold) for e in zip(ranked.source, ranked.target)
        ]
        labeled.to_csv(results / f"evaluation_{name}.csv", index=False)
    pd.DataFrame(metrics).to_csv(results / "m12_metrics.csv", index=False)
    summary = {
        "status": "exploratory_pseudotime_not_causal",
        "seed": seed, "n_cells": len(x), "n_genes": len(x.columns),
        "n_pseudotime_branches": len(pseudo.columns),
        "candidate_rule": manifest["tf_universe_rule"],
        "split": split, "bins": bins, "temporal_weight": temporal_weight,
        "max_edges": max_edges, "ridge_alpha": alpha,
        "reference_audit": provenance,
        "expression_sha256": sha256(raw / "ExpressionData.csv"),
        "pseudotime_sha256": sha256(raw / "PseudoTime.csv"),
        "gold_sha256": sha256(raw / "GroundTruthNetwork.csv"),
        "warning": "Synthetic BoolODE pseudotime may encode actual simulation time; "
                   "this is an optimistic temporal-ordering experiment, not causality.",
    }
    (results / "m12_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", required=True, type=Path)
    p.add_argument("--processed-dir", required=True, type=Path)
    p.add_argument("--results-dir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--bins", type=int, default=16)
    p.add_argument("--max-edges", type=int, default=20)
    p.add_argument("--temporal-weight", type=float, default=.015)
    p.add_argument("--alpha", type=float, default=1.)
    args = p.parse_args()
    print(json.dumps(execute(
        args.raw_dir, args.processed_dir, args.results_dir,
        seed=args.seed, bins=args.bins, max_edges=args.max_edges,
        temporal_weight=args.temporal_weight, alpha=args.alpha
    ), indent=2))


if __name__ == "__main__":
    main()
