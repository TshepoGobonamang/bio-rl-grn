"""Aggregate exploratory M1.1 replicate metrics; keep individual datasets identifiable."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


def aggregate(root: Path):
    runs = []
    for file in sorted(Path(root).glob("*/m11_metrics.csv")):
        frame = pd.read_csv(file)
        frame.insert(0, "replicate", file.parent.name)
        runs.append(frame)
    if not runs:
        raise FileNotFoundError(f"No replicate metrics found under {root}")
    combined = pd.concat(runs, ignore_index=True)
    combined.to_csv(root / "m11_replicate_metrics.csv", index=False)
    cols = ["auroc", "average_precision", "trapezoidal_pr_auc",
            "precision_at_k", "direction_accuracy_ties_half",
            "top_k_reciprocal_pairs"]
    summary = combined.groupby("method", as_index=False)[cols].agg(["mean", "std", "count"])
    summary.columns = ["_".join(x).rstrip("_") for x in summary.columns.to_flat_index()]
    summary.to_csv(root / "m11_descriptive_summary.csv", index=False)
    return combined, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    _, summary = aggregate(args.root)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
