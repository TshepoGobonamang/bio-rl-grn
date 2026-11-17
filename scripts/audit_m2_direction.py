"""Audit M2C direction-sensitive cross-fitted targets before any learner training."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

KEYS = ["dataset_id", "fold_id", "source", "target", "pair_id"]


def audit(root: Path):
    manifests = sorted(Path(root).glob("*/m2c_manifest.json"))
    if not manifests:
        raise FileNotFoundError("No M2C dataset manifests found")
    records = []
    for manifest_file in manifests:
        folder = manifest_file.parent
        meta = json.loads(manifest_file.read_text(encoding="utf-8"))
        features = pd.read_csv(folder / "m2c_features.csv")
        target = pd.read_csv(folder / "m2c_targets.csv")
        if features.columns.tolist() != KEYS + meta["feature_columns"]:
            raise ValueError("Unexpected fit-only features")
        if target.columns.tolist() != KEYS + meta["target_columns"]:
            raise ValueError("Unexpected held-out target columns")
        if len(features) != meta["n_rows"] or not features[KEYS].equals(target[KEYS]):
            raise ValueError("Feature/target row misalignment")
        if features.duplicated(KEYS).any() or target.duplicated(KEYS).any():
            raise ValueError("Duplicate edge-fold sample")
        if meta["uses_gold_reference"] or meta["uses_outer_scoring_or_test_cells"]:
            raise ValueError("Reference or outer holdout leakage indicated by manifest")
        if not all(f["disjoint"] for f in meta["folds"]):
            raise ValueError("Crossfit fold intersection is not empty")
        observed = target.m_direction_TARGET_ONLY.eq(1)
        if (observed & target.y_directional_lag_gain_TARGET.isna()).any():
            raise ValueError("Target missing despite available evidence")
        if ((~observed) & target.y_directional_lag_gain_TARGET.notna()).any():
            raise ValueError("Missing evidence interpreted as observed zero")
        if not target.loc[observed, "y_directional_lag_gain_TARGET"].between(-1., 1.).all():
            raise ValueError("Directional gains exceed bounded target scale")
        directed = dict(zip(
            zip(target.fold_id, target.source, target.target),
            target.y_pairwise_direction_margin_TARGET
        ))
        checked = 0
        nonzero = 0
        for row in target.loc[observed].itertuples(index=False):
            reverse = (row.fold_id, row.target, row.source)
            if reverse in directed and np.isfinite(directed[reverse]):
                checked += 1
                if not np.isclose(row.y_pairwise_direction_margin_TARGET,
                                  -directed[reverse], atol=1e-9):
                    raise ValueError("Pairwise orientation margin is not antisymmetric")
                nonzero += int(abs(row.y_pairwise_direction_margin_TARGET) > 1e-5)
        if checked == 0:
            raise ValueError("No direction-pair comparisons were available")
        usable = target.loc[observed, "y_directional_lag_gain_TARGET"]
        records.append({
            "dataset_id": meta["dataset_id"],
            "edge_fold_records": len(target),
            "available_directional_targets": int(observed.sum()),
            "unavailable_directional_targets": int((~observed).sum()),
            "mean_signed_heldout_gain": float(usable.mean()),
            "fraction_positive_gain": float((usable > 0).mean()),
            "fraction_negative_gain": float((usable < 0).mean()),
            "fraction_clipped_abs_gain_ge_1": float((usable.abs() >= 1).mean()),
            "mean_absolute_direction_margin": float(
                target.y_pairwise_direction_margin_TARGET.abs().mean()),
            "pairwise_checks": checked,
            "nonzero_margin_fraction": nonzero / checked,
            "has_label_leakage": False,
        })
    summary = pd.DataFrame(records)
    root = Path(root)
    summary.to_csv(root / "m2c_quality_summary.csv", index=False)
    (root / "m2c_quality_manifest.json").write_text(json.dumps({
        "status": "M2C-A target construction only, no reliability learner",
        "datasets": summary.dataset_id.tolist(),
        "edge_fold_rows_total": int(summary.edge_fold_records.sum()),
        "available_targets_total": int(summary.available_directional_targets.sum()),
        "reference_edges_seen": 0,
        "note": "Antisymmetry is a mathematical property of a signed forward-vs-reverse gain, not evidence of true causal direction.",
    }, indent=2), encoding="utf-8")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    args = p.parse_args()
    print(audit(args.root).to_string(index=False))


if __name__ == "__main__":
    main()
