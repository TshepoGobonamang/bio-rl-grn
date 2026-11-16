"""Gold-blind cross-fitted temporal-evidence reliability meta-dataset; NOT edge truth."""
from __future__ import annotations
import hashlib
from dataclasses import dataclass, asdict
from typing import Any
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from .pseudotime import temporal_margin
from .ridge import GraphRidge
from .search import universe

@dataclass(frozen=True)
class CrossfitConfig:
    folds: int = 3
    seed: int = 42
    max_bins: int = 16
    min_bins: int = 8
    min_cells_per_bin: int = 3
    lag: int = 1
    temporal_ridge_alpha: float = 1.0
    expression_ridge_alpha: float = 1.0
    stability_temperature: float = 0.05
    bootstrap_repeats: int = 3

    def __post_init__(self):
        if self.folds < 2 or self.min_bins < 8 or self.max_bins < self.min_bins:
            raise ValueError("Invalid fold or pseudotime-bin settings")
        if self.min_cells_per_bin < 2 or self.bootstrap_repeats < 0:
            raise ValueError("Invalid bin support or bootstrap settings")
        if self.stability_temperature <= 0 or self.expression_ridge_alpha <= 0:
            raise ValueError("Temperature and ridge alpha must be positive")

def _hash_ids(index: pd.Index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, index))).encode()).hexdigest()

def _source_margin(expression: pd.DataFrame, pseudotime: pd.DataFrame,
                   tfs: list[str], cfg: CrossfitConfig) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Branch-local support masks: unobserved evidence is NaN, never zero."""
    if not expression.index.equals(pseudotime.index):
        raise ValueError("Misaligned pseudotime and expression")
    edges = universe(list(expression.columns), tfs)
    pairs = pd.DataFrame(edges, columns=["source", "target"])
    weighted = np.zeros(len(edges), dtype=float)
    weight = 0.
    branches_used: list[str] = []
    coverage = np.zeros(len(expression), dtype=bool)
    bin_counts: list[int] = []
    for branch in pseudotime.columns:
        observed = np.isfinite(pseudotime[branch].to_numpy(dtype=float))
        n = int(observed.sum())
        bins = min(cfg.max_bins, n // cfg.min_cells_per_bin)
        if bins < cfg.min_bins:
            continue
        values = temporal_margin(
            expression, pseudotime.loc[:, [branch]], tfs, bins=bins,
            min_cells_per_bin=cfg.min_cells_per_bin,
            lag=cfg.lag, alpha=cfg.temporal_ridge_alpha
        )
        if list(zip(values.source, values.target)) != edges:
            raise ValueError("Temporal estimator changed candidate edge order")
        weighted += n * values.score.to_numpy(dtype=float)
        weight += n
        coverage |= observed
        branches_used.append(str(branch))
        bin_counts.append(int(bins))
    present = weight > 0
    pairs["z_signed"] = weighted / weight if present else np.nan
    pairs["m_present"] = int(present)
    pairs["coverage"] = float(coverage.mean())
    pairs["branch_n"] = len(branches_used)
    pairs["bin_n_min"] = min(bin_counts) if bin_counts else 0
    return pairs, {
        "branches_used": branches_used,
        "branch_n_available": int(len(pseudotime.columns)),
        "branch_n_used": len(branches_used),
        "min_bins_used": int(min(bin_counts)) if bin_counts else 0,
        "usable_cell_coverage": float(coverage.mean()),
        "source_present": bool(present),
    }

def _scaffold_parent(correlation: pd.DataFrame, source: str, target: str,
                     tfs: list[str]) -> tuple[str, ...]:
    choices = [p for p in tfs if p not in {source, target}]
    if not choices:
        return tuple()
    return (sorted(choices, key=lambda p: (-float(correlation.loc[p, target]), p))[0],)

def make_reliability_meta_dataset(
        fitting: pd.DataFrame, pt_fit: pd.DataFrame, tfs: list[str],
        cfg: CrossfitConfig, dataset_id: str = "unknown") -> tuple[pd.DataFrame, dict]:
    """F_b: evidence/scaffold/features. V_b: independent reliability surrogates."""
    if not fitting.index.equals(pt_fit.index):
        raise ValueError("Outer fitting cells and pseudotime must align exactly")
    if fitting.empty or len(fitting) < cfg.folds * 12:
        raise ValueError("Insufficient outer-fitting cells for cross-fitting")
    if fitting.index.has_duplicates or fitting.columns.has_duplicates:
        raise ValueError("Duplicate sample or gene identifiers")
    if set(tfs) - set(fitting.columns) or not tfs:
        raise ValueError("Candidate regulators must be independent measured genes")
    if not np.isfinite(fitting.to_numpy(dtype=float)).all():
        raise ValueError("Expression contains non-finite values")
    candidates = universe(list(fitting.columns), tfs)
    crossfit = KFold(n_splits=cfg.folds, shuffle=True, random_state=cfg.seed)
    results, folds = [], []
    for fold, (train_ix, val_ix) in enumerate(crossfit.split(fitting)):
        F, V = fitting.iloc[train_ix], fitting.iloc[val_ix]
        pF, pV = pt_fit.loc[F.index], pt_fit.loc[V.index]
        e_train, train_meta = _source_margin(F, pF, tfs, cfg)
        e_val, val_meta = _source_margin(V, pV, tfs, cfg)
        if list(zip(e_train.source, e_train.target)) != candidates:
            raise ValueError("Candidate universe mismatch")
        boot = []
        rng = np.random.default_rng(cfg.seed + 1009 * (fold + 1))
        for _ in range(cfg.bootstrap_repeats):
            sampled = rng.integers(0, len(F), len(F))
            boot_F, boot_P = F.iloc[sampled].copy(), pF.iloc[sampled].copy()
            synthetic_ids = [f"bootstrap-{i}" for i in range(len(sampled))]
            boot_F.index = synthetic_ids
            boot_P.index = synthetic_ids
            boost, _ = _source_margin(boot_F, boot_P, tfs, cfg)
            boot.append(boost.z_signed.to_numpy(dtype=float))
        if boot:
            boot_scores = np.stack(boot)
            boot_variance = np.array([
                np.var(col[np.isfinite(col)]) if np.isfinite(col).any() else np.nan
                for col in boot_scores.T
            ])
            signed = np.sign(e_train.z_signed.to_numpy(dtype=float))[None, :] * np.sign(boot_scores)
            bootstrap_agreement = np.array([
                np.mean(col[np.isfinite(col)] > 0) if np.isfinite(col).any() else np.nan
                for col in signed.T
            ])
        else:
            boot_variance = np.full(len(candidates), np.nan)
            bootstrap_agreement = np.full(len(candidates), np.nan)
        correlation = F.corr(method="pearson").abs().fillna(0.0)
        engine = GraphRidge(F, V, alpha=cfg.expression_ridge_alpha)
        fvar = F.var(ddof=0).clip(lower=1e-8)
        for idx, (s, t) in enumerate(candidates):
            base_parents = _scaffold_parent(correlation, s, t, tfs)
            before = engine.target_loss(t, base_parents)
            after = engine.target_loss(t, tuple(sorted((*base_parents, s))))
            rel_gain = (before - after) / max(before, 1e-8)
            z, v = float(e_train.z_signed.iloc[idx]), float(e_val.z_signed.iloc[idx])
            observed = bool(e_train.m_present.iloc[idx] and e_val.m_present.iloc[idx])
            y_stable = (
                float((1.0 + np.tanh(z * v / cfg.stability_temperature**2)) / 2)
                if observed else np.nan
            )
            results.append({
                "dataset_id": dataset_id, "source": s, "target": t,
                "source_id": "beeline_pseudotime", "fold_id": fold,
                "pair_id": "::".join(sorted((s, t))),
                "z_signed": z, "m_present": int(e_train.m_present.iloc[idx]),
                "coverage": float(e_train.coverage.iloc[idx]),
                "branch_n": int(e_train.branch_n.iloc[idx]),
                "bin_n": int(e_train.bin_n_min.iloc[idx]),
                "abs_z": abs(z) if np.isfinite(z) else np.nan,
                "expression_baseline": float(correlation.loc[s, t]),
                "target_fit_variance": float(fvar.loc[t]),
                "n_genes": len(fitting.columns),
                "scaffold_parent": base_parents[0] if base_parents else "",
                "scaffold_n": len(base_parents),
                "bootstrap_variance": float(boot_variance[idx]),
                "bootstrap_sign_agreement": (
                    float(bootstrap_agreement[idx]) if np.isfinite(z) else np.nan
                ),
                "validation_z_signed_TARGET_ONLY": v,
                "validation_available_TARGET_ONLY": int(e_val.m_present.iloc[idx]),
                "y_temporal_reproducibility_TARGET": y_stable,
                "y_expression_conditional_gain_TARGET": float(np.clip(rel_gain, -1., 1.)),
                "expression_loss_before_TARGET_ONLY": float(before),
                "expression_loss_after_TARGET_ONLY": float(after),
            })
        folds.append({
            "fold_id": fold, "fit_n": len(F), "validation_n": len(V),
            "fit_ids_sha256": _hash_ids(F.index),
            "validation_ids_sha256": _hash_ids(V.index),
            "train_evidence": train_meta, "validation_evidence": val_meta,
            "disjoint": set(F.index).isdisjoint(V.index),
        })
    frame = pd.DataFrame.from_records(results)
    features = ["z_signed", "m_present", "coverage", "branch_n", "bin_n",
                "abs_z", "expression_baseline", "target_fit_variance",
                "n_genes", "scaffold_n", "bootstrap_variance", "bootstrap_sign_agreement"]
    targets = ["y_temporal_reproducibility_TARGET", "y_expression_conditional_gain_TARGET"]
    manifest = {
        "dataset_id": dataset_id, "stage": "M2A_crossfit_dataset_only_NO_training",
        "configuration": asdict(cfg), "n_rows": len(frame),
        "n_candidates_per_fold": len(candidates), "n_outer_fit_cells": len(fitting),
        "outer_fit_ids_sha256": _hash_ids(fitting.index),
        "folds": folds, "feature_columns": features, "target_columns": targets,
        "n_temporal_labels_available": int(frame[targets[0]].notna().sum()),
        "n_temporal_evidence_missing": int((frame.m_present == 0).sum()),
        "missing_evidence_is_not_negative": True,
        "uses_gold_reference": False,
        "note": "Temporal target measures out-of-fold sign reproducibility; expression target measures held-out predictive gain. Neither is regulatory truth.",
    }
    if len(frame) != len(candidates) * cfg.folds or not all(f["disjoint"] for f in folds):
        raise RuntimeError("Crossfit integrity failure")
    return frame, manifest
