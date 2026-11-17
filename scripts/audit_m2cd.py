"""M2C-D: strict feature/target alignment and target reproducibility audit."""
from __future__ import annotations
import argparse
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

KEYS = ["dataset_id", "fold_id", "source", "target", "pair_id"]


def audit(root: Path) -> pd.DataFrame:
    entries = sorted(Path(root).glob("*/m2cd_manifest.json"))
    if not entries:
        raise FileNotFoundError("No M2C-D outputs")
    rows = []
    for filename in entries:
        folder = filename.parent
        m = json.loads(filename.read_text(encoding="utf-8"))
        f = pd.read_csv(folder / "m2cd_features.csv")
        t = pd.read_csv(folder / "m2cd_targets.csv")
        if list(f) != KEYS + m["feature_columns"]:
            raise ValueError("Unexpected feature columns or target leakage")
        if list(t) != KEYS + m["target_columns"]:
            raise ValueError("Unexpected target columns")
        if not f[KEYS].equals(t[KEYS]) or f.duplicated(KEYS).any() or t.duplicated(KEYS).any():
            raise ValueError("Feature/target key mismatch or duplication")
        if len(t) != m["n_rows"] or not all(x["disjoint"] for x in m["folds"]):
            raise ValueError("Incorrect records or overlapping inner folds")
        if not m["no_gold_reference_access"] or m["outer_score_test_cells_accessed"]:
            raise ValueError("Evaluation reference/outer test used")
        if ((f.m_independent_prior == 0) & f.independent_prior_score.notna()).any():
            raise ValueError("Missing biological source interpreted as value")
        if ((f.m_independent_prior == 1) & f.independent_prior_score.isna()).any():
            raise ValueError("Present source has missing value")
        if m["dataset_kind"] == "synthetic" and f.m_independent_prior.sum():
            raise ValueError("Synthetic network was given biological evidence")
        before = t.heldout_conditional_before_TARGET_ONLY.to_numpy(dtype=float)
        after = t.heldout_conditional_after_TARGET_ONLY.to_numpy(dtype=float)
        y = t.y_parent_conditional_gain_TARGET.to_numpy(dtype=float)
        recomputed = np.clip((before - after) / np.maximum(before, 1e-10), -1., 1.)
        if not np.isfinite(y).all() or not np.allclose(y, recomputed, rtol=1e-10, atol=1e-10):
            raise ValueError("Held-out conditional targets fail loss reconstruction")
        if not t.y_marginal_expression_gain_TARGET.between(-1, 1).all():
            raise ValueError("Out-of-bounds marginal gain")
        wide = t.pivot(index=["source", "target"], columns="fold_id",
                       values="y_parent_conditional_gain_TARGET")
        rhos = []
        for a, b in combinations(wide.columns, 2):
            if wide[a].nunique() > 1 and wide[b].nunique() > 1:
                rho = spearmanr(wide[a], wide[b]).statistic
                if np.isfinite(rho):
                    rhos.append(float(rho))
        rows.append({
            "dataset_id": m["dataset_id"], "records": len(t),
            "candidate_edges": len(wide), "folds": len(wide.columns),
            "independent_evidence_edges": int(f.loc[f.fold_id.eq(0), "m_independent_prior"].sum()),
            "mean_parent_conditional_gain": float(t.y_parent_conditional_gain_TARGET.mean()),
            "mean_marginal_gain": float(t.y_marginal_expression_gain_TARGET.mean()),
            "fraction_negative_conditional_gain": float((y < 0).mean()),
            "fraction_positive_conditional_gain": float((y > 0).mean()),
            "fraction_clipped_conditional_gain": float((np.abs(y) >= 1).mean()),
            "mean_abs_conditional_minus_marginal": float(np.abs(
                y - t.y_marginal_expression_gain_TARGET.to_numpy(dtype=float)).mean()),
            "conditional_fold_spearman": float(np.mean(rhos)) if rhos else np.nan,
            "gold_labels_used": False,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(Path(root) / "m2cd_quality_summary.csv", index=False)
    (Path(root) / "m2cd_quality_manifest.json").write_text(json.dumps({
        "stage": "M2C-D target quality only",
        "n_networks": len(summary), "total_rows": int(summary.records.sum()),
        "uses_gold_labels": False,
        "notes": "Cross-fold correlations are dependent descriptive statistics; no causality."
    }, indent=2), encoding="utf-8")
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    a = p.parse_args()
    print(audit(a.root).to_string(index=False))


if __name__ == "__main__":
    main()
