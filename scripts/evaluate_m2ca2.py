"""Locked confirmation-only evaluation of M2C-A2 evidence scores.

Reads synthetic gold references ONLY after the gold-blind score files exist.
Never fit or choose blend coefficients using these evaluation results.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from biorlgrn.m11 import audit_reference, evaluate_ranked
from biorlgrn.m2_stability import (
    CONFIRMATION, LOCKED_SEED, LOCKED_BINS, BLEND_Z, BLEND_GAIN,
)

def _scale(v: np.ndarray) -> np.ndarray:
    if not np.isfinite(v).all():
        raise ValueError("Evidence score unavailable: do not silently impute it as negative")
    denominator = float(np.quantile(np.abs(v), .95))
    return np.clip(v / max(denominator, 1e-10), -1., 1.)


def evaluate_one(root: Path, raw_root: Path, name: str) -> pd.DataFrame:
    if name not in CONFIRMATION:
        raise ValueError("Gold reference evaluation restricted to locked confirmation networks")
    folder = Path(root) / name
    manifest = json.loads((folder / "m2ca2_manifest.json").read_text(encoding="utf-8"))
    if manifest["role"] != "confirmation" or manifest["gold_standard_opened"]:
        raise ValueError("Confirmation/provenance guard failed")
    if manifest["predeclared_blend_weights"] != {
            "z_signed": BLEND_Z, "margin": BLEND_GAIN}:
        raise ValueError("Blend weights changed after pre-registration")
    raw = Path(raw_root) / name
    plain = pd.read_csv(folder / "m2ca2_edge_scores.csv")
    original = plain.loc[
        (plain.split_seed == LOCKED_SEED) &
        (plain.bin_cap == LOCKED_BINS) & (~plain.shuffled.astype(bool))
    ].sort_values(["source", "target"]).reset_index(drop=True)
    shuffled = plain.loc[
        (plain.split_seed == LOCKED_SEED) &
        (plain.bin_cap == LOCKED_BINS) & (plain.shuffled.astype(bool))
    ].sort_values(["source", "target"]).reset_index(drop=True)
    if len(original) != len(shuffled) or not original[
            ["source", "target"]].equals(shuffled[["source", "target"]]):
        raise ValueError("Shuffled control must have identical candidate universe")
    if not (original.n_available == 3).all() or not (shuffled.n_available == 3).all():
        raise ValueError("Incomplete fold support; report missingness, do not silently fill")
    # Build all per-edge ranking files before touching the reference.
    scores = {
        "z_signed": original.z_signed.to_numpy(dtype=float),
        "gain": original.gain.to_numpy(dtype=float),
        "direction_margin": original.margin.to_numpy(dtype=float),
        "shuffled_margin": shuffled.margin.to_numpy(dtype=float),
    }
    scores["fixed_blend"] = (
        BLEND_Z * _scale(scores["z_signed"]) +
        BLEND_GAIN * _scale(scores["direction_margin"])
    )
    ranking = {}
    for method, vals in scores.items():
        rank = original[["source", "target"]].copy()
        rank["score"] = vals
        rank = rank.sort_values(
            ["score", "source", "target"], ascending=[False, True, True]
        ).reset_index(drop=True)
        rank.to_csv(folder / f"ranked_{method}.csv", index=False)
        ranking[method] = rank
    # Reference is opened only NOW, and not imported into any score-creation module.
    genes = sorted(set(original.source) | set(original.target))
    gold_df, audit, audit_stats = audit_reference(
        raw / "GroundTruthNetwork.csv", genes, genes
    )
    audit.to_csv(folder / "m2ca2_reference_audit.csv", index=False)
    gold = set(zip(gold_df.source, gold_df.target))
    metrics = []
    direction_correct = {}
    for name_method, ranked in ranking.items():
        m = evaluate_ranked(ranked, gold, genes, genes, k=min(20, len(ranked)))
        metrics.append({
            "dataset_id": name, "method": name_method,
            **{k: v for k, v in m.items() if k != "reciprocal_pairs"},
            "reference_only_evaluation": True,
        })
        if name_method in ("z_signed", "fixed_blend"):
            scores_map = dict(zip(zip(ranked.source, ranked.target), ranked.score))
            direction_correct[name_method] = {
                (s, t): float(scores_map[(s, t)]) > float(scores_map[(t, s)])
                for s, t in gold if (t, s) not in gold
            }
    old, new = direction_correct["z_signed"], direction_correct["fixed_blend"]
    changes = {
        "reference_eligible_unique_directions": len(old),
        "fixed_blend_corrects_old_errors": sum(not old[k] and new[k] for k in old),
        "fixed_blend_introduces_errors": sum(old[k] and not new[k] for k in old),
        "reference_audit": audit_stats,
        "note": "The 'fixed_blend' mixture was predeclared; no gold-based tuning.",
    }
    (folder / "m2ca2_orientation_changes.json").write_text(
        json.dumps(changes, indent=2), encoding="utf-8")
    return pd.DataFrame(metrics)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-root", type=Path, required=True)
    p.add_argument("--raw-root", type=Path, required=True)
    args = p.parse_args()
    missing = [name for name in CONFIRMATION
               if not (args.results_root / name / "m2ca2_edge_scores.csv").is_file()]
    if missing:
        raise FileNotFoundError(f"No gold-independent scores for {missing}")
    all_metrics = pd.concat([
        evaluate_one(args.results_root, args.raw_root, n) for n in CONFIRMATION
    ], ignore_index=True)
    all_metrics.to_csv(
        args.results_root / "m2ca2_locked_confirmation_metrics.csv", index=False
    )
    cols = ["dataset_id", "method", "average_precision", "auroc",
            "direction_accuracy_ties_half", "precision_at_k"]
    print(all_metrics[cols].to_string(index=False))


if __name__ == "__main__":
    main()
