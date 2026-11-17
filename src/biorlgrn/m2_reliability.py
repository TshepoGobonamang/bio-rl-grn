"""M2B: gold-blind, leave-one-network-out reliability-gate comparisons.

Outputs measure stability pseudo-target prediction, NOT regulatory truth.
M2A meta-features are constructed from fitting cells only; validation labels are
read solely to fit models on non-held-out network structures and to evaluate
predictions *after* the held-out model has been frozen.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

KEYS = ["dataset_id", "fold_id", "source", "target", "pair_id"]
TARGET = "y_temporal_reproducibility_TARGET"
FEATURES = [
    "z_signed", "m_present", "coverage", "branch_n", "bin_n",
    "abs_z", "expression_baseline", "target_fit_variance",
    "n_genes", "scaffold_n", "bootstrap_variance",
    "bootstrap_sign_agreement",
]
ALPHA = 30.0


@dataclass
class M2AData:
    features: pd.DataFrame
    targets: pd.DataFrame
    sources: tuple[str, ...]


def load_m2a(root: Path) -> M2AData:
    """Read separate M2A files, never the mixed diagnostic CSV."""
    root = Path(root)
    inputs, labels, source_ids = [], [], []
    for manifest_path in sorted(root.glob("*/m2_crossfit_manifest.json")):
        folder = manifest_path.parent
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("uses_gold_reference") is not False:
            raise ValueError(f"Unverified gold-label isolation in {folder}")
        if manifest.get("outer_score_and_test_excluded_from_dataset") is not True:
            raise ValueError(f"Outer holdout was accessed in {folder}")
        if manifest["feature_columns"] != FEATURES:
            raise ValueError("Unexpected M2A feature allowlist")
        if not all(f.get("disjoint", False) for f in manifest["folds"]):
            raise ValueError("An M2A inner split overlaps")
        f = pd.read_csv(folder / "m2_features.csv")
        t = pd.read_csv(folder / "m2_targets.csv")
        if f.columns.tolist() != KEYS + FEATURES:
            raise ValueError("Feature file has unexpected or target-derived columns")
        if any(c.endswith(("_TARGET", "_TARGET_ONLY")) for c in f):
            raise ValueError("Target leakage into predictor features")
        if not set(KEYS + [TARGET]).issubset(t):
            raise ValueError("Target file does not expose the required soft label")
        if (t.columns[:len(KEYS)].tolist() != KEYS or
                not f[KEYS].equals(t[KEYS]) or
                f.duplicated(KEYS).any() or
                t.duplicated(KEYS).any()):
            raise ValueError("Target/feature edge-fold alignment failure")
        if len(f) != manifest["n_rows"] or f.dataset_id.nunique() != 1:
            raise ValueError("Invalid dataset cardinality")
        if ((f.m_present == 0) & f.z_signed.notna()).any():
            raise ValueError("Missing temporal evidence was encoded as present")
        if ((f.m_present == 0) & t[TARGET].notna()).any():
            raise ValueError("Missing temporal evidence has a soft target")
        if t[TARGET].dropna().lt(0).any() or t[TARGET].dropna().gt(1).any():
            raise ValueError("Soft labels must lie in [0,1]")
        inputs.append(f)
        labels.append(t.loc[:, KEYS + [TARGET]])
        source_ids.append(str(f.dataset_id.iloc[0]))
    if len(inputs) < 3 or len(set(source_ids)) != len(source_ids):
        raise ValueError("At least three distinct network contexts required for M2B")
    return M2AData(
        features=pd.concat(inputs, ignore_index=True),
        targets=pd.concat(labels, ignore_index=True),
        sources=tuple(source_ids),
    )


def _heuristic(x: pd.DataFrame) -> np.ndarray:
    """Train-only bootstrap sign agreement as a simple nonlearned gate."""
    r = x.bootstrap_sign_agreement.fillna(0.5).to_numpy(dtype=float)
    return np.where(x.m_present.to_numpy() > 0, np.clip(r, 0., 1.), 0.)


def _predictor() -> Pipeline:
    # All transforms are FIT on source networks only. Parameters are frozen.
    return Pipeline([
        ("imputer", SimpleImputer(strategy="constant", fill_value=0.0,
                                  keep_empty_features=True)),
        ("scale", StandardScaler()),
        ("ridge", Ridge(alpha=ALPHA)),
    ])


def leave_one_network_out(data: M2AData):
    """Fit three independent models: each held-out topology is never used for fitting.

    A model is fit on *two other network contexts*. All folds of the held-out
    context are evaluated together; edge-fold rows are not i.i.d. replications.
    """
    x, y = data.features, data.targets
    if not x[KEYS].equals(y[KEYS]):
        raise ValueError("Features/labels misaligned")
    records, summaries, weights = [], [], []
    for held_out in data.sources:
        is_test = x.dataset_id.eq(held_out).to_numpy()
        is_train = ~is_test
        valid_train = is_train & y[TARGET].notna().to_numpy() & x.m_present.eq(1).to_numpy()
        if valid_train.sum() < 10:
            raise ValueError(f"Not enough source-topology training labels for {held_out}")
        train_x = x.loc[valid_train, FEATURES].copy()
        train_y = y.loc[valid_train, TARGET].astype(float).to_numpy()
        # Each source topology contributes equal total sample weight despite
        # different candidate-universe sizes (306 vs 42).
        ids = x.loc[valid_train, "dataset_id"]
        counts = ids.value_counts()
        sw = ids.map(lambda v: 1.0 / counts[v]).to_numpy(dtype=float)
        sw *= len(sw) / sw.sum()
        model = _predictor()
        model.fit(train_x, train_y, ridge__sample_weight=sw)
        test_x = x.loc[is_test, FEATURES].copy()
        availability = test_x.m_present.to_numpy(dtype=int)
        learned = np.where(
            availability > 0,
            np.clip(model.predict(test_x), 0., 1.),
            0.,
        )
        global_rho = float(np.average(train_y, weights=sw))
        estimates = {
            "fixed": np.where(availability > 0, 1., 0.),
            "global": np.where(availability > 0, global_rho, 0.),
            "stability": _heuristic(test_x),
            "adaptive": learned,
        }
        true = y.loc[is_test, TARGET].to_numpy(dtype=float)
        for method, values in estimates.items():
            row = x.loc[is_test, KEYS].copy()
            row["method"] = method
            row["rho"] = values
            row["target_soft_for_evaluation_only"] = true
            records.append(row)
            observed = np.isfinite(true)
            summaries.append({
                "held_out_network": held_out, "method": method,
                "n_eval_rows": int(observed.sum()),
                "n_unique_gene_pairs": int(x.loc[is_test, "pair_id"].nunique()),
                "mse_reproducibility": float(np.mean((values[observed] - true[observed]) ** 2)),
                "mae_reproducibility": float(np.mean(np.abs(values[observed] - true[observed]))),
                "mean_predicted_reliability": float(values[observed].mean()),
                "mean_target_reproducibility": float(true[observed].mean()),
                "train_networks": ",".join(v for v in data.sources if v != held_out),
                "trained_on_held_out_network": False,
                "training_rows": int(valid_train.sum()),
            })
            # Average three inner-fold scores into one gate per directed edge.
            row_no_target = row.drop(columns=["target_soft_for_evaluation_only"])
            edge_scores = row_no_target.groupby(["dataset_id", "source", "target", "method"],
                                               as_index=False)["rho"].mean()
            weights.append(edge_scores)
    pred = pd.concat(records, ignore_index=True)
    scores = pd.DataFrame(summaries)
    edges = pd.concat(weights, ignore_index=True)
    if not pred.rho.between(0., 1.).all() or not edges.rho.between(0., 1.).all():
        raise RuntimeError("Reliability weights must be in [0,1]")
    if pred[KEYS + ["method"]].duplicated().any() or edges[
            ["dataset_id", "source", "target", "method"]].duplicated().any():
        raise RuntimeError("Duplicated reliability predictions")
    return pred, scores, edges


def export_m2b(root: Path, output: Path):
    """Export diagnostics separately from weights consumed by GRN graph search."""
    dataset = load_m2a(root)
    pred, score, weights = leave_one_network_out(dataset)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    pred.to_csv(out / "m2b_meta_evaluation_only.csv", index=False)
    score.to_csv(out / "m2b_lono_metrics.csv", index=False)
    weights.to_csv(out / "m2b_inference_weights.csv", index=False)
    meta = {
        "status": "exploratory_meta_reliability_no_RL",
        "held_out_networks": list(dataset.sources),
        "target": TARGET,
        "training_feature_columns": FEATURES,
        "alpha_fixed_in_advance": ALPHA,
        "n_meta_training_rows_total": len(dataset.features),
        "n_meta_validation_rows": len(pred),
        "no_gold_label_training": True,
        "all_networks_held_out_once": True,
        "experimental_note": (
            "Soft target is reproducibility of temporal association, not true "
            "TF regulation. Direction-reversed edges can share the same target."
        ),
        "warning": "No random edge-fold validation and no inference of regulatory truth.",
    }
    (out / "m2b_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return score, meta
