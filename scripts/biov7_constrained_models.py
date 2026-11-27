"""BIO-V7: TF-constrained matched one-edit predictive experiments.

This is an explicit adaptation of ``grnbench.gcev``, not a redefinition of the
historical synthetic benchmark. Every target has one common TF-only scaffold;
every compared action inserts exactly one eligible TF. No routine repartitions
the study cohort or reads any expression file. Callers own the frozen fit,
internal-score and test boundaries.

Nested source-utility labels use genuinely disjoint training/search-score/
validation subsets of the supplied fitting cells. Final scaffolds are selected
on all supplied fitting cells, then coefficients/gains are cross-fitted while
holding those scaffolds fixed. Those final gains are selection statistics, not
unbiased cross-validation of the full scaffold-selection procedure. The outer
internal-score endpoint must remain unseen until decisions are frozen.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

from grnbench.core import _loss, _rank_unit
from grnbench.gcev import _corr_state, _winner

VARIANCE_FLOOR = 1e-9


def _validate(fit: pd.DataFrame, edges: list[tuple[str, str]]) -> None:
    if not fit.index.is_unique or not fit.columns.is_unique:
        raise ValueError("Expression identifiers must be unique")
    if not np.isfinite(fit.to_numpy(float)).all():
        raise ValueError("Expression must be finite")
    if len(edges) != len(set(edges)) or any(s == t for s, t in edges):
        raise ValueError("Candidate edges must be unique and exclude self-edges")
    if any(s not in fit or t not in fit for s, t in edges):
        raise ValueError("Candidate genes are missing from expression")
    if set(t for _, t in edges) != set(fit.columns):
        raise ValueError("Every expression target needs an eligible candidate set")


def _candidate_map(edges: list[tuple[str, str]]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for source, target in edges:
        result.setdefault(target, []).append(source)
    return {target: sorted(sources) for target, sources in result.items()}


def tf_scaffold(fit: pd.DataFrame, target: str,
                eligible_tfs: list[str], size: int) -> tuple[str, ...]:
    """Highest absolute fit-only correlations among this target's eligible TFs."""
    if size < 0:
        raise ValueError("Scaffold size cannot be negative")
    ranked = []
    y = fit[target].to_numpy(float)
    for tf in sorted(eligible_tfs):
        if tf == target:
            raise ValueError("Scaffold eligibility contains a self-edge")
        with np.errstate(divide="ignore", invalid="ignore"):
            value = np.corrcoef(fit[tf].to_numpy(float), y)[0, 1]
        ranked.append((abs(float(value)) if np.isfinite(value) else 0.0, tf))
    ranked.sort(key=lambda item: (-item[0], item[1]))
    if len(ranked) < size + 2:
        raise ValueError("Each target needs at least two residual candidate edits")
    return tuple(tf for _, tf in ranked[:size])


def _context_gains(search_fit: pd.DataFrame, search_score: pd.DataFrame,
                   target: str, state: tuple[str, ...], candidates: list[str],
                   alpha: float) -> np.ndarray:
    base = _loss(search_fit, search_score, target, state, alpha)
    return np.asarray([
        np.clip((base - _loss(search_fit, search_score, target,
                             state + (tf,), alpha)) / max(base, VARIANCE_FLOOR),
                -1.0, 1.0)
        for tf in candidates
    ], dtype=float)


def nested_valuation(fit_df: pd.DataFrame,
                     edges: list[tuple[str, str]],
                     evidence: Mapping[str, tuple[np.ndarray, np.ndarray]],
                     *, alpha: float = 1.0, max_parents: int = 2,
                     boost: float = .3, seed: int = 20261127,
                     n_folds: int = 3) -> pd.DataFrame:
    """Existing GCEV FEATURES-compatible labels under common TF-only scaffolds.

    An exactly constant effective source vector cannot alter the chosen edit.
    Its action and utility are explicitly kept equal, preventing floating-point
    rounding of a common offset from being mistaken for a biological effect.
    This is an algebraic invariant, not an outcome-dependent decision.
    """
    _validate(fit_df, edges)
    if len(fit_df) < 48 or n_folds < 3 or not 0 < boost <= 1:
        raise ValueError("Need >=48 fitting cells, >=3 folds and 0<boost<=1")
    if alpha <= 0 or not evidence:
        raise ValueError("Positive ridge alpha and at least one source required")
    for name, (values, mask) in evidence.items():
        if np.shape(values) != (len(edges),) or np.shape(mask) != (len(edges),):
            raise ValueError(f"Evidence shape mismatch: {name}")
        if not np.isfinite(values).all() or np.any(np.asarray(values) < 0) or np.any(np.asarray(values) > 1):
            raise ValueError(f"Evidence values must be finite in [0,1]: {name}")
    eligible = _candidate_map(edges)
    edge_ix = {edge: i for i, edge in enumerate(edges)}
    sources = sorted(evidence)
    rows = []
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    for fold, (train_index, valid_index) in enumerate(kf.split(fit_df)):
        train, validation = fit_df.iloc[train_index], fit_df.iloc[valid_index]
        perm = np.random.default_rng(seed + 51 * fold + 13).permutation(len(train))
        cut = max(10, int(len(train) * .75))
        if cut >= len(train) - 3:
            raise ValueError("Insufficient cells for nested source valuation")
        search_fit, search_score = train.iloc[perm[:cut]], train.iloc[perm[cut:]]
        contexts = {}
        for target in fit_df.columns:
            state = tf_scaffold(search_fit, target, eligible[target], max_parents)
            candidates = [tf for tf in eligible[target] if tf not in state]
            gains = _context_gains(search_fit, search_score, target, state,
                                   candidates, alpha)
            contexts[target] = (state, candidates, gains)
        ranks = _rank_unit(np.concatenate([v[2] for v in contexts.values()]))
        offset = 0
        for target, (state, candidates, gains) in contexts.items():
            base_rank = ranks[offset:offset + len(candidates)]
            offset += len(candidates)
            indices = np.asarray([edge_ix[(tf, target)] for tf in candidates])
            ev = {s: np.asarray(evidence[s][0], float)[indices] for s in sources}
            masks = {s: np.asarray(evidence[s][1], bool)[indices] for s in sources}
            for source in sources:
                vals, mask = ev[source], masks[source]
                other = [s for s in sources if s != source]
                others = sum((boost * masks[s] * ev[s] for s in other),
                             start=np.zeros(len(candidates)))
                minus_scores = base_rank + others
                effective = mask * vals
                constant = bool(np.ptp(effective) == 0.0)
                baseline_index = _winner(candidates, minus_scores)
                boosted_index = (baseline_index if constant else
                                 _winner(candidates, minus_scores + boost * effective))
                minus_loss = _loss(train, validation, target,
                                   state + (candidates[baseline_index],), alpha)
                plus_loss = (minus_loss if boosted_index == baseline_index else
                             _loss(train, validation, target,
                                   state + (candidates[boosted_index],), alpha))
                delta = float(np.clip((minus_loss - plus_loss) /
                                     max(minus_loss, VARIANCE_FLOOR), -1, 1))
                observed = vals[mask]
                other_mask = (np.logical_or.reduce([masks[s] for s in other])
                              if other else np.zeros(len(candidates), bool))
                disagreements = [abs(float(ev[s][i] - vals[i]))
                                 for s in other for i in range(len(candidates))
                                 if mask[i] and masks[s][i]]
                winning_prior = float(vals[boosted_index]) if mask[boosted_index] else 0.
                baseline_prior = float(vals[baseline_index]) if mask[baseline_index] else 0.
                rows.append({
                    "fold": fold, "target": str(target), "evidence_source": source,
                    "graph_parents": "|".join(state),
                    "minus_edge": candidates[baseline_index],
                    "plus_edge": candidates[boosted_index],
                    "different_edge": int(baseline_index != boosted_index),
                    "minus_validation_nmse": float(minus_loss),
                    "plus_validation_nmse": float(plus_loss),
                    "source_value": delta,
                    "equal_edit_budget": True,
                    "used_evaluation_reference": False,
                    "tf_only_scaffold": True,
                    "source_constant_within_target": constant,
                    "effective_prior_range": float(np.ptp(effective)),
                    "source_coverage": float(mask.mean()),
                    "source_mean": float(observed.mean()) if len(observed) else 0.,
                    "source_max": float(observed.max()) if len(observed) else 0.,
                    "source_std": float(observed.std()) if len(observed) else 0.,
                    "other_source_coverage": float(other_mask.mean()),
                    "source_disagreement": float(np.mean(disagreements)) if disagreements else 0.,
                    "candidate_count": float(len(candidates)),
                    "parent_count": float(len(state)),
                    "parent_target_correlation": _corr_state(search_fit, target, state),
                    "baseline_gain": float(gains[baseline_index]),
                    "baseline_margin": float(np.sort(minus_scores)[-1] - np.sort(minus_scores)[-2]),
                    "winning_source_prior": winning_prior,
                    "action_prior_contrast": winning_prior - baseline_prior,
                })
    result = pd.DataFrame(rows)
    if result.duplicated(["fold", "target", "evidence_source"]).any():
        raise AssertionError("Duplicate source valuation rows")
    return result


def final_contexts(fit_df: pd.DataFrame, edges: list[tuple[str, str]],
                   *, alpha: float = 1.0, max_parents: int = 2,
                   seed: int = 20261127, n_folds: int = 3) -> pd.DataFrame:
    """Common full-fit scaffold plus fit-only cross-fitted residual-edit gains.

    Return every input edge in its original order. Scaffold edges have NaN gain
    and rank because they are already present and are not candidate insertions.
    Callers MUST filter ``eligible_action`` before composing/choosing actions.
    """
    _validate(fit_df, edges)
    if len(fit_df) < 48 or n_folds < 3 or alpha <= 0:
        raise ValueError("Need >=48 fitting cells, >=3 folds, positive alpha")
    eligible = _candidate_map(edges)
    states = {t: tf_scaffold(fit_df, t, eligible[t], max_parents)
              for t in fit_df.columns}
    candidate_edges = [(s, t) for s, t in edges if s not in states[t]]
    fold_gains = np.zeros((n_folds, len(candidate_edges)))
    residual = {t: [s for s in eligible[t] if s not in states[t]]
                for t in fit_df.columns}
    by_edge: dict[tuple[str, str], int] = {edge: i for i, edge in enumerate(candidate_edges)}
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    for fold, (train_index, valid_index) in enumerate(kf.split(fit_df)):
        for target in fit_df.columns:
            values = _context_gains(fit_df.iloc[train_index], fit_df.iloc[valid_index],
                                    target, states[target], residual[target], alpha)
            for source, gain in zip(residual[target], values):
                fold_gains[fold, by_edge[(source, target)]] = gain
    mean_gains = fold_gains.mean(axis=0)
    ranks = _rank_unit(mean_gains)
    rows = []
    for tf, target in edges:
        scaffold_edge = tf in states[target]
        index = by_edge.get((tf, target))
        rows.append({"tf": tf, "target": target,
                     "graph_parents": "|".join(states[target]),
                     "scaffold_edge": scaffold_edge,
                     "eligible_action": not scaffold_edge,
                     "gain": np.nan if scaffold_edge else float(mean_gains[index]),
                     "expr_rank": np.nan if scaffold_edge else float(ranks[index]),
                     "scaffold_selected_on_full_fit": True,
                     "gain_coefficient_folds": n_folds})
    return pd.DataFrame(rows)


def evaluate_graphs(fit_df: pd.DataFrame, evaluation_df: pd.DataFrame,
                    selections: pd.DataFrame, *, alpha: float = 1.0
                    ) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Evaluate frozen matched graph choices, returning per-target and cell errors.

    ``selections`` requires one row per (method,target), with ``graph_parents``
    pipe-separated and ``selected_tf`` the single inserted TF. Every method must
    use identical target-specific scaffolds. Error matrices are evaluation
    cells by fit_df.columns order and normalized using fitting target variance.
    Caller may use paired cell resampling conditional on this one fitted cohort;
    that is not biological-replicate or cross-network uncertainty.
    """
    required = {"method", "target", "graph_parents", "selected_tf"}
    if not required.issubset(selections.columns):
        raise ValueError("Selections require method,target,graph_parents,selected_tf")
    if selections.duplicated(["method", "target"]).any():
        raise ValueError("Duplicate target selection for a method")
    if list(fit_df.columns) != list(evaluation_df.columns):
        raise ValueError("Evaluation genes/order differ from fitting genes")
    if not fit_df.index.is_unique or not evaluation_df.index.is_unique:
        raise ValueError("Cell identifiers must be unique")
    if len(fit_df.index.intersection(evaluation_df.index)):
        raise ValueError("Evaluation cells overlap fitting cells")
    if not np.isfinite(fit_df.to_numpy(float)).all() or not np.isfinite(evaluation_df.to_numpy(float)).all():
        raise ValueError("Nonfinite expression")
    if selections.groupby("target").graph_parents.nunique().max() != 1:
        raise ValueError("Compared methods use unequal scaffolds")
    target_order = list(fit_df.columns)
    rows, cell_errors = [], {}
    for method, frame in selections.groupby("method", sort=True):
        if set(frame.target) != set(target_order):
            raise ValueError("Each method must cover exactly the same targets")
        frame = frame.set_index("target")
        errors = np.zeros((len(evaluation_df), len(target_order)), dtype=float)
        for j, target in enumerate(target_order):
            row = frame.loc[target]
            state = tuple(s for s in str(row.graph_parents).split("|") if s)
            source = str(row.selected_tf)
            if source in state or target in state or source == target or len(set(state)) != len(state):
                raise ValueError("Graph action is a self-edge, duplicate, or already present")
            parents = state + (source,)
            if any(p not in fit_df for p in parents):
                raise ValueError("Graph parent missing from fitting data")
            y, yv = fit_df[target].to_numpy(float), evaluation_df[target].to_numpy(float)
            x, xv = fit_df[list(parents)].to_numpy(float), evaluation_df[list(parents)].to_numpy(float)
            mu, sd = x.mean(axis=0), np.maximum(x.std(axis=0), VARIANCE_FLOOR)
            model = Ridge(alpha=alpha).fit((x - mu) / sd, y)
            prediction = model.predict((xv - mu) / sd)
            sq_error = (yv - prediction) ** 2
            variance = float(np.var(y))
            denominator = max(variance, VARIANCE_FLOOR)
            errors[:, j] = sq_error / denominator
            rows.append({"method": method, "target": target,
                         "graph_parents": "|".join(state), "selected_tf": source,
                         "total_parent_count": len(parents),
                         "n_fit_cells": len(fit_df), "n_evaluation_cells": len(evaluation_df),
                         "fit_target_variance": variance,
                         "variance_floor_applied": variance < VARIANCE_FLOOR,
                         "nmse_denominator": denominator,
                         "mse": float(sq_error.mean()),
                         "nmse": float(errors[:, j].mean())})
        cell_errors[str(method)] = errors
    return pd.DataFrame(rows), cell_errors
