"""M2C-A: direction-sensitive out-of-fold temporal *predictive* usefulness.

A candidate s→t is assessed by fitting binned lagged regressions on F_b
and scoring them on disjoint V_b cells. This is NOT causal edge truth.
No gold regulatory network or outer scoring/test cells are used.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from .m2_crossfit import CrossfitConfig, make_reliability_meta_dataset
from .search import universe


def _hash_ids(ids: pd.Index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode()).hexdigest()


def _branch_bins(x: pd.DataFrame, p: pd.Series, bins: int) -> np.ndarray:
    """Return time-ordered *independent-cell* bin means, not trajectories."""
    time = p.to_numpy(dtype=float)
    ids = np.where(np.isfinite(time))[0]
    ids = ids[np.argsort(time[ids], kind="stable")]
    groups = np.array_split(ids, bins)
    return np.stack([x.iloc[g].to_numpy(dtype=float).mean(axis=0) for g in groups])


def _predict_ridge(fx: np.ndarray, fy: np.ndarray, vx: np.ndarray,
                   alpha: float) -> np.ndarray:
    """Fit only on F; apply F-scaling to V; do not penalize intercept."""
    means = fx.mean(axis=0)
    scales = np.maximum(fx.std(axis=0), 1e-8)
    train = (fx - means) / scales
    valid = (vx - means) / scales
    coef = np.linalg.solve(
        train.T @ train + alpha * np.eye(train.shape[1]),
        train.T @ (fy - fy.mean())
    )
    return fy.mean() + valid @ coef


def _directed_heldout_gain(train: np.ndarray, valid: np.ndarray,
                           s: int, t: int, alpha: float,
                           lag: int) -> tuple[float, float, float]:
    """Fit autoregression + source predictor on TRAIN bins; score on VAL bins.

    Both candidate and baseline predict the SAME target, making loss differences
    directed. Includes target's own lag and a shared bin-time trend.
    """
    n = len(train) - lag
    if n < 6 or len(valid) - lag < 6:
        return np.nan, np.nan, np.nan
    # Per-gene normalization learned from fitting bins only.
    center = train.mean(axis=0)
    scale = np.maximum(train.std(axis=0), 1e-8)
    f = (train - center) / scale
    v = (valid - center) / scale
    times_f = np.linspace(-1, 1, len(f) - lag)
    times_v = np.linspace(-1, 1, len(v) - lag)
    xb_f = np.column_stack([f[:-lag, t], times_f])
    xb_v = np.column_stack([v[:-lag, t], times_v])
    xf_f = np.column_stack([xb_f, f[:-lag, s]])
    xf_v = np.column_stack([xb_v, v[:-lag, s]])
    y_f = f[lag:, t]
    y_v = v[lag:, t]
    base = _predict_ridge(xb_f, y_f, xb_v, alpha)
    full = _predict_ridge(xf_f, y_f, xf_v, alpha)
    loss_base = float(np.mean((y_v - base) ** 2))
    loss_full = float(np.mean((y_v - full) ** 2))
    gain = (loss_base - loss_full) / max(loss_base, 1e-8)
    return float(np.clip(gain, -1., 1.)), loss_base, loss_full


def heldout_directional_targets(fit: pd.DataFrame, valid: pd.DataFrame,
                                 pt_fit: pd.DataFrame, pt_valid: pd.DataFrame,
                                 tfs: list[str], cfg: CrossfitConfig):
    """One target per admissible *directed* edge; branch-local masks preserved."""
    if not fit.index.equals(pt_fit.index) or not valid.index.equals(pt_valid.index):
        raise ValueError("Pseudotime must align to expression cell identifiers")
    if list(fit.columns) != list(valid.columns):
        raise ValueError("Gene order must match")
    if list(pt_fit.columns) != list(pt_valid.columns):
        raise ValueError("Pseudotime branch order must match")
    genes = list(fit.columns)
    weighted = np.zeros((len(genes), len(genes)))
    base_losses = np.zeros_like(weighted)
    full_losses = np.zeros_like(weighted)
    weights = np.zeros_like(weighted)
    branch_info = []
    for branch in pt_fit.columns:
        nf = int(pt_fit[branch].notna().sum())
        nv = int(pt_valid[branch].notna().sum())
        nbin = min(cfg.max_bins, nf // cfg.min_cells_per_bin,
                   nv // cfg.min_cells_per_bin)
        if nbin < cfg.min_bins or nbin - cfg.lag < 6:
            branch_info.append({"branch": str(branch), "used": False,
                                "fit_cells": nf, "validation_cells": nv})
            continue
        f = _branch_bins(fit, pt_fit[branch], nbin)
        v = _branch_bins(valid, pt_valid[branch], nbin)
        for si in range(len(genes)):
            for ti in range(len(genes)):
                if si == ti:
                    continue
                g, before, after = _directed_heldout_gain(
                    f, v, si, ti, cfg.temporal_ridge_alpha, cfg.lag)
                if np.isfinite(g):
                    weighted[si, ti] += nv * g
                    base_losses[si, ti] += nv * before
                    full_losses[si, ti] += nv * after
                    weights[si, ti] += nv
        branch_info.append({"branch": str(branch), "used": True,
                            "fit_cells": nf, "validation_cells": nv,
                            "bins": nbin})
    index = {g: i for i, g in enumerate(genes)}
    gain = np.divide(weighted, weights, out=np.full_like(weights, np.nan),
                     where=weights > 0)
    base = np.divide(base_losses, weights, out=np.full_like(weights, np.nan),
                     where=weights > 0)
    full = np.divide(full_losses, weights, out=np.full_like(weights, np.nan),
                     where=weights > 0)
    records = []
    for s, t in universe(genes, tfs):
        si, ti = index[s], index[t]
        available = bool(weights[si, ti] > 0 and weights[ti, si] > 0)
        records.append({
            "source": s, "target": t,
            "y_directional_lag_gain_TARGET": float(gain[si, ti]) if available else np.nan,
            "y_pairwise_direction_margin_TARGET": (
                float(gain[si, ti] - gain[ti, si]) if available else np.nan
            ),
            "validation_lag_base_loss_TARGET_ONLY": (
                float(base[si, ti]) if available else np.nan
            ),
            "validation_lag_full_loss_TARGET_ONLY": (
                float(full[si, ti]) if available else np.nan
            ),
            "m_direction_TARGET_ONLY": int(available),
        })
    return pd.DataFrame.from_records(records), branch_info


def make_directional_meta_dataset(fitting: pd.DataFrame, pt_fit: pd.DataFrame,
                                   tfs: list[str], cfg: CrossfitConfig,
                                   dataset_id: str):
    """Reuse M2A FIT-only feature schema; replace symmetric target with a held-out directional target."""
    original, m2a = make_reliability_meta_dataset(
        fitting, pt_fit, tfs, cfg, dataset_id=dataset_id
    )
    keys = ["dataset_id", "fold_id", "source", "target", "pair_id"]
    features = original.loc[:, keys + m2a["feature_columns"]].copy()
    if not features.columns.is_unique or any("_TARGET" in col for col in features.columns):
        raise RuntimeError("Target contamination in features")
    folds = KFold(n_splits=cfg.folds, shuffle=True, random_state=cfg.seed)
    all_targets = []
    fold_audit = []
    for fold_id, (fit_idx, val_idx) in enumerate(folds.split(fitting)):
        f, v = fitting.iloc[fit_idx], fitting.iloc[val_idx]
        a, b = pt_fit.loc[f.index], pt_fit.loc[v.index]
        directional, branches = heldout_directional_targets(f, v, a, b, tfs, cfg)
        directional.insert(0, "fold_id", fold_id)
        directional.insert(0, "dataset_id", dataset_id)
        directional["pair_id"] = [
            "::".join(sorted((s, t))) for s, t in zip(
                directional.source, directional.target
            )
        ]
        all_targets.append(directional)
        fold_audit.append({
            "fold_id": fold_id, "fit_n": len(f), "validation_n": len(v),
            "fit_ids_sha256": _hash_ids(f.index),
            "validation_ids_sha256": _hash_ids(v.index),
            "disjoint": f.index.intersection(v.index).empty,
            "branches": branches,
        })
    target = pd.concat(all_targets, ignore_index=True)
    target = target[keys + [c for c in target.columns if c not in keys]]
    if not features[keys].equals(target[keys]):
        raise RuntimeError("Crossfit feature/target identity mismatch")
    if features.duplicated(keys).any() or target.duplicated(keys).any():
        raise RuntimeError("Duplicate edge-fold records")
    eligible = target.m_direction_TARGET_ONLY.eq(1)
    if eligible.any():
        if target.loc[eligible, "y_directional_lag_gain_TARGET"].isna().any():
            raise RuntimeError("Missing directional usefulness target")
        if not target.loc[eligible, "y_directional_lag_gain_TARGET"].between(-1., 1.).all():
            raise RuntimeError("Directed gains out of supported range")
        diffs = dict(zip(zip(target.dataset_id, target.fold_id,
                              target.source, target.target),
                         target.y_pairwise_direction_margin_TARGET))
        for row in target.loc[eligible].itertuples(index=False):
            reverse = (row.dataset_id, row.fold_id, row.target, row.source)
            if reverse in diffs and np.isfinite(diffs[reverse]):
                if not np.isclose(row.y_pairwise_direction_margin_TARGET,
                                  -diffs[reverse], atol=1e-9):
                    raise RuntimeError("Direction margin is not antisymmetric")
    elif m2a["n_temporal_labels_available"]:
        raise RuntimeError("All held-out directional targets unexpectedly unavailable")
    manifest = {
        "dataset_id": dataset_id,
        "stage": "M2C_A_directional_crossfit_TARGET_CONSTRUCTION_ONLY",
        "configuration": asdict(cfg),
        "outer_fit_ids_sha256": _hash_ids(fitting.index),
        "n_outer_fit_cells": len(fitting),
        "n_rows": len(target),
        "n_candidates_per_fold": m2a["n_candidates_per_fold"],
        "n_directional_targets_available": int(eligible.sum()),
        "n_directional_targets_missing": int((~eligible).sum()),
        "feature_columns": m2a["feature_columns"],
        "target_columns": [c for c in target if c.endswith(("_TARGET", "_TARGET_ONLY"))],
        "data_keys": keys,
        "folds": fold_audit,
        "uses_gold_reference": False,
        "uses_outer_scoring_or_test_cells": False,
        "missing_is_not_zero": True,
        "method_note": (
            "Train-fold binned lagged models evaluated on independent validation "
            "cells, with target autocorrelation and time trend controlled. This "
            "is directional predictive association, NOT causal regulation."
        ),
    }
    return features, target, manifest
