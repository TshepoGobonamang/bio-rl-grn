"""Audit M2A predictor/target boundaries and summarize BEELINE pseudo-targets.

This is a quality-control utility, not reliability-model training or GRN evaluation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def audit(root: Path, require_targets: bool = True) -> pd.DataFrame:
    root = Path(root)
    manifests = sorted(root.glob("*/m2_crossfit_manifest.json"))
    if not manifests:
        raise FileNotFoundError(f"No M2A cross-fit manifests under {root}")
    reports = []
    for path in manifests:
        folder = path.parent
        m = json.loads(path.read_text(encoding="utf-8"))
        f = pd.read_csv(folder / "m2_features.csv")
        t = pd.read_csv(folder / "m2_targets.csv")
        keys = m["data_keys"]
        allowed = keys + m["feature_columns"]
        target_cols = [c for c in t if c not in keys]
        if list(f.columns) != allowed:
            raise ValueError(f"Unexpected predictor columns in {folder}")
        if not target_cols or not all(c.endswith(("_TARGET", "_TARGET_ONLY")) for c in target_cols):
            raise ValueError(f"Unsafe target schema in {folder}")
        if len(f) != m["n_rows"] or len(f) != len(t):
            raise ValueError(f"Incorrect row count in {folder}")
        if f.duplicated(keys).any() or t.duplicated(keys).any() or not f[keys].equals(t[keys]):
            raise ValueError(f"Unaligned / duplicate fold-edge keys in {folder}")
        if m["uses_gold_reference"] or not m["outer_score_and_test_excluded_from_dataset"]:
            raise ValueError(f"Gold or outer holdout exposure in {folder}")
        if not all(part["disjoint"] for part in m["folds"]):
            raise ValueError(f"Train/validation overlap in {folder}")
        if ((f.m_present == 0) & f.z_signed.notna()).any():
            raise ValueError(f"Missing evidence encoded as an observed score in {folder}")
        if ((t.validation_available_TARGET_ONLY == 0) &
                t.validation_z_signed_TARGET_ONLY.notna()).any():
            raise ValueError(f"Missing validation evidence encoded as a value in {folder}")
        valid = ((f.m_present == 1) & (t.validation_available_TARGET_ONLY == 1))
        if (valid & t.y_temporal_reproducibility_TARGET.isna()).any():
            raise ValueError(f"Missing reproducibility target for available evidence in {folder}")
        if ((~valid) & t.y_temporal_reproducibility_TARGET.notna()).any():
            raise ValueError(f"Target populated for unavailable evidence in {folder}")
        usable = t.loc[valid, "y_temporal_reproducibility_TARGET"]
        if require_targets and usable.empty:
            raise ValueError(f"No available reliability pseudo-targets in {folder}")
        if (usable.lt(0) | usable.gt(1)).any():
            raise ValueError(f"Pseudo-target outside 0–1 in {folder}")
        expr_gain = t.y_expression_conditional_gain_TARGET
        corr = usable.corr(expr_gain.loc[valid]) if len(usable) > 2 else np.nan
        reports.append({
            "dataset_id": m["dataset_id"],
            "n_features": len(f),
            "n_folds": int(f.fold_id.nunique()),
            "n_unique_edges": int(f[["source", "target"]].drop_duplicates().shape[0]),
            "train_only_evidence_fraction": float(f.m_present.mean()),
            "target_availability_fraction": float(valid.mean()),
            "mean_soft_target": float(usable.mean()) if len(usable) else np.nan,
            "fraction_target_extreme_002_098": float(
                ((usable < .02) | (usable > .98)).mean()
            ) if len(usable) else np.nan,
            "mean_conditional_expression_gain": float(expr_gain.mean()),
            "fraction_positive_expression_gain": float((expr_gain > 0).mean()),
            "temporal_expression_gain_correlation": float(corr),
            "gold_reference_accessed": False,
            "feature_target_separation_verified": True,
        })
    result = pd.DataFrame(reports)
    result.to_csv(root / "m2a_quality_summary.csv", index=False)
    report = {
        "datasets": result.dataset_id.tolist(),
        "n_datasets": len(result),
        "total_edge_fold_records": int(result.n_features.sum()),
        "n_reference_labels_used": 0,
        "status": "M2A audit only; not model fitting or biological validation",
        "warning": "Fold-level edge rows are dependent; do not use random-row split for meta-model validation.",
    }
    (root / "m2a_quality_manifest.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--allow-empty-targets", action="store_true")
    args = parser.parse_args()
    print(audit(args.root, require_targets=not args.allow_empty_targets).to_string(index=False))


if __name__ == "__main__":
    main()
