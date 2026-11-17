"""M2C-D: cross-fitted source-blind *parent-conditional* expression usefulness.

No reference network edges, evaluation labels or external priors participate in
scaffold selection or target fitting. Biological evidence is a separate feature
with explicit availability; target is predictive association, not causality.
"""
from __future__ import annotations
import hashlib
from dataclasses import asdict, dataclass
from itertools import product
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from .m2_crossfit import CrossfitConfig, _source_margin
from .m2_external import align_evidence
from .ridge import GraphRidge
from .search import universe


@dataclass(frozen=True)
class ParentConfig:
    folds: int = 3
    seed: int = 42
    max_parents: int = 2
    ridge_alpha: float = 1.0
    max_bins: int = 12
    def __post_init__(self):
        if self.folds < 2 or self.max_parents < 0 or self.ridge_alpha <= 0:
            raise ValueError("Invalid parent-conditional configuration")


def _hash_ids(index: pd.Index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, index))).encode()).hexdigest()


def _scaffold(F: pd.DataFrame, source: str, target: str, n: int,
              correlation: pd.DataFrame | None = None) -> tuple[str, ...]:
    """Source/target excluded; correlation cache must originate in F only."""
    if n <= 0:
        return ()
    corr = F.corr().abs().fillna(0.) if correlation is None else correlation
    if not F.columns.equals(corr.index) or not F.columns.equals(corr.columns):
        raise ValueError("Scaffold correlation must use the inner-fit gene universe")
    choices = [g for g in F.columns if g not in {source, target}]
    ranked = sorted(choices, key=lambda g: (-float(corr.loc[g, target]), g))
    return tuple(sorted(ranked[:n]))


def parent_conditional_dataset(
    fit: pd.DataFrame, pseudo: pd.DataFrame | None,
    tfs: list[str], cfg: ParentConfig,
    dataset_id: str, *, biological_links: pd.DataFrame | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Generate *separate* fit-only feature and V-fold target files.

    Outer fitting cells only; KFold inner validation is not used for predictors
    or scaffold selection. Parent effect is direct ridge prediction of measured
    target at the same cell, not a lagged or causal temporal effect.
    """
    if fit.empty or fit.index.has_duplicates or fit.columns.has_duplicates:
        raise ValueError("Empty, duplicate-cell or duplicate-gene expression")
    if len(fit) < cfg.folds * 12 or not np.isfinite(fit.to_numpy(dtype=float)).all():
        raise ValueError("Insufficient or non-finite expression")
    if pseudo is not None and not fit.index.equals(pseudo.index):
        raise ValueError("Pseudotime must align to expression cells")
    if set(tfs) - set(fit.columns):
        raise ValueError("Unknown candidate TF")
    candidates = universe(list(fit.columns), tfs)
    aligned = align_evidence(
        pd.DataFrame(candidates, columns=["source", "target"]), biological_links
    )
    if list(zip(aligned.source, aligned.target)) != candidates:
        raise RuntimeError("Biological prior join changed candidate order")
    splits = KFold(n_splits=cfg.folds, shuffle=True, random_state=cfg.seed)
    feature_rows, target_rows, audits = [], [], []
    for fold, (i_f, i_v) in enumerate(splits.split(fit)):
        F, V = fit.iloc[i_f], fit.iloc[i_v]
        if not F.index.intersection(V.index).empty:
            raise RuntimeError("Overlapping inner fold")
        scorer = GraphRidge(F, V, alpha=cfg.ridge_alpha)
        if pseudo is not None:
            e, provenance = _source_margin(
                F, pseudo.loc[F.index], tfs,
                CrossfitConfig(folds=cfg.folds, seed=cfg.seed,
                               max_bins=cfg.max_bins, bootstrap_repeats=0)
            )
            if list(zip(e.source, e.target)) != candidates:
                raise RuntimeError("Candidate order not preserved")
        else:
            e = pd.DataFrame({"z_signed": np.full(len(candidates), np.nan),
                              "m_present": np.zeros(len(candidates), dtype=int)})
            provenance = {"source_present": False}
        correlations = F.corr().abs().fillna(0.)
        variance = F.var(ddof=0).clip(lower=1e-8)
        for k, (s, t) in enumerate(candidates):
            others = _scaffold(F, s, t, cfg.max_parents, correlations)
            before = scorer.target_loss(t, others)
            after = scorer.target_loss(t, tuple(sorted((*others, s))))
            # Parallel marginal comparator, not the learned target.
            marginal_before = scorer.target_loss(t, ())
            marginal_after = scorer.target_loss(t, (s,))
            gain = float(np.clip((before-after)/max(before, 1e-10), -1., 1.))
            marginal = float(np.clip(
                (marginal_before-marginal_after)/max(marginal_before, 1e-10), -1., 1.
            ))
            keys = {"dataset_id": dataset_id, "fold_id": fold,
                    "source": s, "target": t,
                    "pair_id": "::".join(sorted((s, t)))}
            feature_rows.append({
                **keys,
                "z_signed": float(e.z_signed.iloc[k]),
                "m_temporal": int(e.m_present.iloc[k]),
                "m_independent_prior": int(aligned.m_independent_prior.iloc[k]),
                "independent_prior_score": float(aligned.independent_prior_score.iloc[k]),
                "expression_abs_correlation": float(correlations.loc[s, t]),
                "target_fit_variance": float(variance.loc[t]),
                "scaffold_n": len(others),
                "scaffold_parents": "|".join(others),
            })
            target_rows.append({
                **keys,
                "y_parent_conditional_gain_TARGET": gain,
                "y_marginal_expression_gain_TARGET": marginal,
                "heldout_conditional_before_TARGET_ONLY": before,
                "heldout_conditional_after_TARGET_ONLY": after,
                "m_target_TARGET_ONLY": 1,
            })
        audits.append({
            "fold_id": fold, "fit_n": len(F), "validation_n": len(V),
            "fit_ids_sha256": _hash_ids(F.index),
            "validation_ids_sha256": _hash_ids(V.index),
            "disjoint": True, "temporal_evidence": provenance
        })
    keys = ["dataset_id", "fold_id", "source", "target", "pair_id"]
    f = pd.DataFrame.from_records(feature_rows)
    t = pd.DataFrame.from_records(target_rows)
    if not f[keys].equals(t[keys]) or f.duplicated(keys).any() or t.duplicated(keys).any():
        raise RuntimeError("Unsafe fold-edge alignment")
    manifest = {
        "stage": "M2C_D_parent_conditional_crossfit_NO_model_training",
        "dataset_id": dataset_id, "configuration": asdict(cfg),
        "n_outer_fit_cells": len(fit), "n_candidates_per_fold": len(candidates),
        "n_rows": len(f), "n_independent_links": int(aligned.m_independent_prior.sum()),
        "independent_biological_prior_present": biological_links is not None,
        "outer_fit_ids_sha256": _hash_ids(fit.index),
        "data_keys": keys,
        "feature_columns": [c for c in f if c not in keys],
        "target_columns": [c for c in t if c not in keys],
        "folds": audits, "no_gold_reference_access": True,
        "outer_score_test_cells_accessed": False,
        "target_is_causal_regulatory_truth": False,
        "warning": "Source-blind parent-conditioned predictive association, not biological causation."
    }
    return f, t, manifest
