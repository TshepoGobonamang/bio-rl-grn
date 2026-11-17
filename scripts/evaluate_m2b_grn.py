"""Evaluate frozen reliability gates in matched BEELINE graph construction.

Trains NOTHING: the M2B gates in `m2b_inference_weights.csv` were fit
on other network topologies. Reference gold is opened only after rankings
and independent expression losses have already been saved.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from biorlgrn.io import load_expression, sha256
from biorlgrn.m11 import split_three, independent_expression_loss, audit_reference, evaluate_ranked
from biorlgrn.pseudotime import load_pseudotime, temporal_margin, temporal_greedy
from biorlgrn.search import greedy_rank, universe
from biorlgrn.ridge import GraphRidge

METHODS = ("expression_only", "fixed", "global", "stability", "adaptive",
           "shuffled_time")
WEIGHT = 0.015


def evaluate_network(raw_dir: Path, weights_csv: Path, output_dir: Path, *,
                     dataset_id: str | None = None, seed: int = 42,
                     max_edges: int = 20, bins: int = 16) -> dict:
    raw_dir, output_dir = Path(raw_dir), Path(output_dir)
    dataset_id = dataset_id or raw_dir.name
    output_dir.mkdir(parents=True, exist_ok=True)
    x = load_expression(raw_dir / "ExpressionData.csv", gene_axis="rows")
    pseudotime = load_pseudotime(raw_dir / "PseudoTime.csv", x.index)
    # Keep all genes as candidate regulators as in M1.2. No gold-derived TF list.
    tfs = list(x.columns)
    fit, graph_score, untouched, split = split_three(x, seed=seed)
    engine = GraphRidge(fit, graph_score, alpha=1.)
    temporal = temporal_margin(fit, pseudotime.loc[fit.index], tfs, bins=bins)
    universe_pairs = universe(list(x.columns), tfs)
    if list(zip(temporal.source, temporal.target)) != universe_pairs:
        raise ValueError("Temporal evidence and candidate universe mismatch")
    saved_weights = pd.read_csv(weights_csv)
    required = {"dataset_id", "source", "target", "method", "rho"}
    if not required.issubset(saved_weights):
        raise ValueError("Missing reliability-gate weight columns")
    ranks, trajectories = {}, {}
    ranks["expression_only"], trajectories["expression_only"] = greedy_rank(
        engine, tfs, max_edges=max_edges)
    for method in ("fixed", "global", "stability", "adaptive"):
        sub = saved_weights.loc[
            saved_weights.dataset_id.eq(dataset_id) & saved_weights.method.eq(method)
        ].copy()
        edges = list(zip(sub.source, sub.target))
        if len(edges) != len(universe_pairs) or len(set(edges)) != len(edges) or set(edges) != set(universe_pairs):
            raise ValueError(f"Missing or duplicate weights for {dataset_id}/{method}")
        if not sub.rho.between(0., 1.).all():
            raise ValueError("Reliability weights must be bounded")
        gate = dict(zip(edges, sub.rho.astype(float)))
        gated = temporal.copy()
        gated["score"] = [
            val * gate[(source, target)]
            for source, target, val in zip(temporal.source, temporal.target, temporal.score)
        ]
        ranks[method], trajectories[method] = temporal_greedy(
            engine, tfs, gated, max_edges=max_edges, temporal_weight=WEIGHT)
    # Negative control preserves branch missingness and is independent of gold.
    shuffled = pseudotime.loc[fit.index].copy()
    rng = np.random.default_rng(seed + 830)
    for branch in shuffled:
        mask = shuffled[branch].notna()
        values = shuffled.loc[mask, branch].to_numpy(copy=True)
        rng.shuffle(values)
        shuffled.loc[mask, branch] = values
    permuted_margin = temporal_margin(fit, shuffled, tfs, bins=bins)
    ranks["shuffled_time"], trajectories["shuffled_time"] = temporal_greedy(
        engine, tfs, permuted_margin, max_edges=max_edges, temporal_weight=WEIGHT)
    for name in METHODS:
        ranks[name].to_csv(output_dir / f"ranked_{name}.csv", index=False)
        trajectories[name].to_csv(output_dir / f"trajectory_{name}.csv", index=False)
    test_losses = {}
    # Only after all graphs have been selected do we inspect untouched cells.
    for name in METHODS:
        path = trajectories[name]
        selected = set(zip(path.source, path.target))
        test_losses[name] = independent_expression_loss(fit, untouched, selected, 1.)
    # Only after inference and untouched expression testing do we open gold.
    gold_path = raw_dir / "GroundTruthNetwork.csv"
    gold_df, row_audit, provenance = audit_reference(gold_path, list(x.columns), tfs)
    row_audit.to_csv(output_dir / "reference_audit.csv", index=False)
    gold = set(zip(gold_df.source, gold_df.target))
    metrics = []
    for name in METHODS:
        rank = ranks[name]
        m = evaluate_ranked(rank, gold, list(x.columns), tfs, max_edges)
        metrics.append({
            "dataset_id": dataset_id, "method": name,
            **{k: v for k, v in m.items() if k != "reciprocal_pairs"},
            **test_losses[name],
        })
        annotated = rank.copy()
        annotated["evaluation_reference_label"] = [
            int((s, t) in gold) for s, t in zip(rank.source, rank.target)
        ]
        annotated.to_csv(output_dir / f"evaluation_{name}.csv", index=False)
    pd.DataFrame(metrics).to_csv(output_dir / "m2b_grn_metrics.csv", index=False)
    manifest = {
        "dataset_id": dataset_id,
        "stage": "M2B exploratory held-out topology gate evaluation",
        "n_cells": len(x), "n_genes": len(x.columns),
        "candidate_universe": len(universe_pairs),
        "n_gold_eligible": len(gold),
        "fixed_temporal_weight": WEIGHT,
        "seed": seed, "bins": bins, "max_edges": max_edges,
        "split": split, "reference_provenance": provenance,
        "expression_sha256": sha256(raw_dir / "ExpressionData.csv"),
        "pseudotime_sha256": sha256(raw_dir / "PseudoTime.csv"),
        "gold_sha256": sha256(gold_path),
        "no_gold_in_reliability_or_graph_training": True,
        "warning": (
            "Temporal reproducibility is a surrogate, not biological truth; "
            "this uses previously inspected BEELINE topologies and is exploratory."
        ),
    }
    (output_dir / "m2b_grn_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--weights-csv", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--dataset-id")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-edges", type=int, default=20)
    p.add_argument("--bins", type=int, default=16)
    a = p.parse_args()
    print(json.dumps(evaluate_network(
        a.raw_dir, a.weights_csv, a.output_dir,
        dataset_id=a.dataset_id, seed=a.seed, max_edges=a.max_edges,
        bins=a.bins
    ), indent=2))


if __name__ == "__main__":
    main()
