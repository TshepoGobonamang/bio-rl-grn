"""M2C-A2: gold-blind stability/sensitivity experiment for directional evidence.

Keep all evaluation references out of this module. Individual cells are
cross-sectional; lagged pseudotime-bin predictors do NOT establish causality.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import KFold

from .io import load_expression, sha256
from .m11 import split_three
from .m2_crossfit import CrossfitConfig, _source_margin
from .m2_direction import heldout_directional_targets
from .pseudotime import load_pseudotime
from .search import universe

KEYS = ["dataset_id", "setting", "split_seed", "bin_cap", "shuffled",
        "fold_id", "source", "target"]
DEVELOPMENT = ("dyn-LL-500-1", "dyn-LI-500-1", "dyn-BF-500-1")
CONFIRMATION = ("dyn-CY-500-1", "dyn-TF-500-1", "dyn-BFC-500-1")
DEV_SEEDS = (42, 73, 137)
DEV_BINS = (8, 12, 16)
LOCKED_SEED = 42
LOCKED_BINS = 12
BLEND_Z = 0.75
BLEND_GAIN = 0.25


def _ids_hash(ids: pd.Index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode()).hexdigest()


def shuffle_pseudotime(p: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Permute observed branch values, preserving missingness and cell IDs."""
    rng = np.random.default_rng(seed)
    out = p.copy()
    for c in out.columns:
        mask = out[c].notna()
        values = out.loc[mask, c].to_numpy(copy=True)
        rng.shuffle(values)
        out.loc[mask, c] = values
    return out


def run_setting(x: pd.DataFrame, pt: pd.DataFrame, dataset_id: str, *,
                seed: int, bins: int, shuffled: bool = False) -> tuple[pd.DataFrame, list]:
    """Evaluate all directed edges under one inner split, fold by fold."""
    if not x.index.equals(pt.index):
        raise ValueError("Pseudotime / expression cell IDs are not aligned")
    if x.index.has_duplicates or x.columns.has_duplicates:
        raise ValueError("Duplicate cell or gene identifiers")
    cfg = CrossfitConfig(folds=3, seed=seed, min_bins=8, max_bins=bins,
                         min_cells_per_bin=3, bootstrap_repeats=0)
    setting = f"s{seed}_b{bins}" + ("_shuffled" if shuffled else "")
    split = KFold(n_splits=3, shuffle=True, random_state=seed)
    candidates = universe(list(x.columns), list(x.columns))
    all_rows, fold_info = [], []
    for fold_id, (fidx, vidx) in enumerate(split.split(x)):
        f, v = x.iloc[fidx], x.iloc[vidx]
        pf, pv = pt.loc[f.index], pt.loc[v.index]
        if shuffled:
            pf = shuffle_pseudotime(pf, seed + 1709 * (fold_id + 1))
            pv = shuffle_pseudotime(pv, seed + 3793 * (fold_id + 1))
        evidence, evidence_info = _source_margin(f, pf, list(x.columns), cfg)
        target, branch_info = heldout_directional_targets(
            f, v, pf, pv, list(x.columns), cfg
        )
        if list(zip(evidence.source, evidence.target)) != candidates:
            raise RuntimeError("Candidate edge order changed during feature generation")
        if list(zip(target.source, target.target)) != candidates:
            raise RuntimeError("Candidate edge order changed during target construction")
        if not f.index.intersection(v.index).empty:
            raise RuntimeError("F and V share cells")
        row = evidence[["source", "target", "z_signed", "m_present"]].copy()
        row["gain"] = target["y_directional_lag_gain_TARGET"].to_numpy()
        row["margin"] = target["y_pairwise_direction_margin_TARGET"].to_numpy()
        row["target_available"] = target["m_direction_TARGET_ONLY"].to_numpy()
        for name, val in (("dataset_id", dataset_id), ("setting", setting),
                          ("split_seed", seed), ("bin_cap", bins),
                          ("shuffled", shuffled), ("fold_id", fold_id)):
            row.insert(0, name, val)
        row = row.loc[:, KEYS + ["z_signed", "m_present", "gain", "margin",
                                "target_available"]]
        all_rows.append(row)
        fold_info.append({
            "fold_id": fold_id, "fit_n": len(f), "validation_n": len(v),
            "fit_ids_sha256": _ids_hash(f.index),
            "validation_ids_sha256": _ids_hash(v.index),
            "disjoint": True, "source_info": evidence_info, "branches": branch_info,
        })
    output = pd.concat(all_rows, ignore_index=True)
    if output.duplicated(KEYS).any() or len(output) != 3 * len(candidates):
        raise RuntimeError("Duplicate or missing fold-edge records")
    return output, fold_info


def _pairwise_fold_rho(frame: pd.DataFrame, field: str) -> float:
    wide = frame.pivot(index=["source", "target"], columns="fold_id", values=field)
    correlations = []
    for f1, f2 in itertools.combinations(wide.columns, 2):
        a, b = wide[f1], wide[f2]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 4 or a[ok].nunique() < 2 or b[ok].nunique() < 2:
            continue
        r = spearmanr(a[ok], b[ok]).statistic
        if np.isfinite(r):
            correlations.append(float(r))
    return float(np.mean(correlations)) if correlations else np.nan


def _sign_agreement(frame: pd.DataFrame, field: str) -> float:
    wide = frame.pivot(index=["source", "target"], columns="fold_id", values=field)
    good = wide.notna().all(axis=1)
    a = np.sign(wide.loc[good].to_numpy(dtype=float))
    return float(((a > 0).all(axis=1) | (a < 0).all(axis=1)).mean()) if len(a) else np.nan


def aggregate_setting(frame: pd.DataFrame) -> pd.DataFrame:
    if frame[["setting", "dataset_id"]].drop_duplicates().shape[0] != 1:
        raise ValueError("Aggregate exactly one network and one setting")
    return frame.groupby(["dataset_id", "setting", "split_seed", "bin_cap",
                          "shuffled", "source", "target"], as_index=False).agg(
        z_signed=("z_signed", "mean"), gain=("gain", "mean"),
        margin=("margin", "mean"),
        n_available=("target_available", "sum"),
        n_observed=("m_present", "sum"),
    )


def summarize_setting(frame: pd.DataFrame) -> dict:
    out = {
        "dataset_id": str(frame.dataset_id.iloc[0]),
        "setting": str(frame.setting.iloc[0]),
        "bin_cap": int(frame.bin_cap.iloc[0]),
        "split_seed": int(frame.split_seed.iloc[0]),
        "shuffled": bool(frame.shuffled.iloc[0]),
        "edge_fold_rows": len(frame),
        "available_target_fraction": float(frame.target_available.mean()),
        "mean_gain": float(frame.gain.mean()),
        "negative_gain_fraction": float((frame.gain < 0).mean()),
        "gain_clipping_fraction": float((frame.gain.abs() >= 1).mean()),
    }
    for field in ("z_signed", "gain", "margin"):
        out[f"{field}_fold_spearman"] = _pairwise_fold_rho(frame, field)
        out[f"{field}_fold_sign_consistency"] = _sign_agreement(frame, field)
    return out


def run_network(raw_dir: Path, output_dir: Path, *, role: str) -> dict:
    """Generate development sensitivity or locked confirmation diagnostics."""
    raw = Path(raw_dir)
    name = raw.name
    if role == "development" and name not in DEVELOPMENT:
        raise ValueError(f"{name}: not a predefined development network")
    if role == "confirmation" and name not in CONFIRMATION:
        raise ValueError(f"{name}: not a predefined confirmation network")
    expr_file = raw / "ExpressionData.csv"
    pt_file = raw / "PseudoTime.csv"
    x = load_expression(expr_file, gene_axis="rows")
    pt = load_pseudotime(pt_file, x.index)
    outer, score, untouched, split = split_three(x, seed=42)
    settings = (
        tuple(itertools.product(DEV_SEEDS, DEV_BINS))
        if role == "development" else ((LOCKED_SEED, LOCKED_BINS),)
    )
    all_frames, metadata = [], []
    for seed, bins in settings:
        frame, folds = run_setting(
            outer, pt.loc[outer.index], name, seed=seed, bins=bins)
        all_frames.append(frame)
        metadata.append({"seed": seed, "bins": bins, "shuffled": False,
                         "folds": folds})
    # One fixed-per-network shuffled negative control. No gold-based tuning.
    control, folds = run_setting(
        outer, pt.loc[outer.index], name,
        seed=LOCKED_SEED, bins=LOCKED_BINS, shuffled=True)
    all_frames.append(control)
    metadata.append({"seed": LOCKED_SEED, "bins": LOCKED_BINS, "shuffled": True,
                     "folds": folds})
    folds_table = pd.concat(all_frames, ignore_index=True)
    summary = pd.DataFrame([
        summarize_setting(g) for _, g in folds_table.groupby("setting", sort=True)
    ])
    aggregate = pd.concat([
        aggregate_setting(g) for _, g in folds_table.groupby("setting", sort=True)
    ], ignore_index=True)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    folds_table.to_csv(out / "m2ca2_fold_scores.csv", index=False)
    summary.to_csv(out / "m2ca2_stability.csv", index=False)
    aggregate.to_csv(out / "m2ca2_edge_scores.csv", index=False)
    manifest = {
        "dataset_id": name, "role": role, "n_genes": len(x.columns),
        "n_outer_fit_cells": len(outer), "n_outer_score_cells_excluded": len(score),
        "n_outer_test_cells_excluded": len(untouched),
        "outer_split": split, "settings": metadata,
        "development_networks": DEVELOPMENT,
        "confirmation_networks": CONFIRMATION,
        "predeclared_blend_weights": {"z_signed": BLEND_Z, "margin": BLEND_GAIN},
        "expression_sha256": sha256(expr_file),
        "pseudotime_sha256": sha256(pt_file),
        "gold_standard_opened": False,
        "no_outer_score_test_access": True,
        "note": "Independent-cell pseudotime-bin association, not causal TF direction.",
    }
    import json
    (out / "m2ca2_manifest.json").write_text(json.dumps(
        manifest, indent=2, default=lambda o: list(o) if isinstance(o, tuple) else str(o)
    ), encoding="utf-8")
    return {"dataset_id": name, "role": role,
            "n_rows": len(folds_table),
            "n_settings": len(metadata)}
