"""Graph-Conditional Evidence Valuation (GCEV), experimental v0.2.

Gold-independent, fit-only nested validation estimates *predictive usefulness*
of providing one evidence source to a fixed-budget graph edit. It does not
estimate a causal TF-target relationship, calibrated edge probability or a
biological source's reliability in unseen tissues.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold

from .core import _loss, _scaffold, _rank_unit


SOURCE_FEATURES = (
    'source_coverage', 'source_mean', 'source_max', 'source_std',
    'other_source_coverage', 'source_disagreement', 'candidate_count',
)
GRAPH_FEATURES = (
    'parent_count', 'parent_target_correlation', 'baseline_gain',
    'baseline_margin', 'winning_source_prior', 'action_prior_contrast',
)
FEATURES = SOURCE_FEATURES + GRAPH_FEATURES


@dataclass
class GCEVResult:
    per_target_source: pd.DataFrame
    paired_decisions: pd.DataFrame
    gate_fit_metadata: dict


def _winner(candidates: list[str], scores: np.ndarray) -> int:
    """Deterministic, source-lexicographic tie breaker; no gold labels."""
    return min(range(len(candidates)), key=lambda i: (-float(scores[i]), candidates[i]))


def _corr_state(train: pd.DataFrame, target: str, parents: tuple[str, ...]) -> float:
    if not parents:
        return 0.0
    ct = train[target].to_numpy(float)
    vals = []
    for name in parents:
        with np.errstate(invalid='ignore', divide='ignore'):
            v = np.corrcoef(train[name].to_numpy(float), ct)[0, 1]
        vals.append(abs(float(v)) if np.isfinite(v) else 0.0)
    return float(np.mean(vals))


def _one_context(train: pd.DataFrame, validation: pd.DataFrame,
                 edges: list[tuple[str, str]],
                 evidence: Mapping[str, tuple[np.ndarray, np.ndarray]],
                 *, alpha: float, max_parents: int, boost: float,
                 split_seed: int, fold: int) -> list[dict]:
    """Matched one-edit graph outcomes, source-ablated vs source-present.

    Choose edit using *training-only* inner score, measure source contribution
    on the distinct inner validation partition. Both competitors use the same
    parent scaffold and exactly one candidate insertion (same action budget).
    """
    rng = np.random.default_rng(split_seed)
    perm = rng.permutation(len(train))
    cut = max(10, int(len(train) * .75))
    if cut >= len(train) - 3:
        raise ValueError('Not enough cells for the nested source-valuation split')
    search_fit = train.iloc[perm[:cut]]
    search_score = train.iloc[perm[cut:]]
    edge_ix = {e: i for i, e in enumerate(edges)}
    sources = sorted(evidence)
    records = []
    # Normalize expression gains across the full fitted candidate universe,
    # matching the final ranking's global normalization. Normalizing within
    # a three-edge target would create rank jumps of 0.5 and make any source
    # boost <=0.3 mathematically incapable of changing the chosen edit.
    contexts = {}
    for target in train.columns:
        possible = sorted(s for s, t in edges if t == target)
        if len(possible) < 2:
            continue
        state = _scaffold(search_fit, '', target, max_parents)
        candidates = [s for s in possible if s not in state]
        if len(candidates) < 2:
            continue
        base_loss = _loss(search_fit, search_score, target, state, alpha)
        gains = np.asarray([
            float(np.clip((base_loss - _loss(search_fit, search_score,
                     target, state + (source,), alpha)) / max(base_loss, 1e-9), -1, 1))
            for source in candidates
        ])
        contexts[target] = (state, candidates, gains)
    if not contexts:
        return []
    ordered_targets = list(contexts)
    all_scores = np.concatenate([contexts[t][2] for t in ordered_targets])
    global_rank = _rank_unit(all_scores)
    rank_map = {}
    offset = 0
    for target in ordered_targets:
        n = len(contexts[target][2])
        rank_map[target] = global_rank[offset:offset+n]
        offset += n
    for target in ordered_targets:
        state, candidates, gains = contexts[target]
        base_rank = rank_map[target]
        evid_values, evid_mask = {}, {}
        for name in sources:
            vals, mask = evidence[name]
            idx = [edge_ix[(source, target)] for source in candidates]
            evid_values[name] = np.asarray(vals[idx], dtype=float)
            evid_mask[name] = np.asarray(mask[idx], dtype=bool)
        for name in sources:
            vals = evid_values[name]
            mask = evid_mask[name]
            others = np.zeros(len(candidates))
            for other in sources:
                if other != name:
                    others += boost * evid_mask[other] * evid_values[other]
            minus_scores = base_rank + others
            plus_scores = minus_scores + boost * mask * vals
            baseline_idx = _winner(candidates, minus_scores)
            boosted_idx = _winner(candidates, plus_scores)
            # Refit both graph actions on precisely the SAME training sample;
            # the source cannot see these validation losses while selecting.
            minus_loss = _loss(train, validation, target,
                               state + (candidates[baseline_idx],), alpha)
            plus_loss = _loss(train, validation, target,
                              state + (candidates[boosted_idx],), alpha)
            delta = float(np.clip((minus_loss - plus_loss) / max(minus_loss, 1e-9), -1, 1))
            observed_vals = vals[mask]
            other_mask = np.logical_or.reduce([evid_mask[n] for n in sources if n != name]) if len(sources)>1 else np.zeros(len(candidates),bool)
            # Missing source evidence stays masked; it is not negative evidence.
            disagreement = float(np.mean([
                abs(float(evid_values[n][i] - vals[i]))
                for n in sources if n != name
                for i in range(len(candidates)) if mask[i] and evid_mask[n][i]
            ])) if any(mask[i] and evid_mask[n][i] for n in sources if n != name for i in range(len(candidates))) else 0.0
            winning_prior = float(vals[boosted_idx]) if mask[boosted_idx] else 0.0
            baseline_prior = float(vals[baseline_idx]) if mask[baseline_idx] else 0.0
            features = {
                'source_coverage':float(mask.mean()),
                'source_mean':float(observed_vals.mean()) if len(observed_vals) else 0.,
                'source_max':float(observed_vals.max()) if len(observed_vals) else 0.,
                'source_std':float(observed_vals.std()) if len(observed_vals) else 0.,
                'other_source_coverage':float(other_mask.mean()),
                'source_disagreement':disagreement,
                'candidate_count':float(len(candidates)),
                'parent_count':float(len(state)),
                'parent_target_correlation':_corr_state(search_fit, target, state),
                'baseline_gain':float(gains[baseline_idx]),
                'baseline_margin':float(np.sort(minus_scores)[-1]-np.sort(minus_scores)[-2]),
                'winning_source_prior':winning_prior,
                'action_prior_contrast':float(winning_prior-baseline_prior),
            }
            records.append({
                'fold':fold,'target':str(target),'evidence_source':str(name),
                'graph_parents':'|'.join(state),'minus_edge':candidates[baseline_idx],
                'plus_edge':candidates[boosted_idx],
                'different_edge':int(baseline_idx!=boosted_idx),
                'minus_validation_nmse':float(minus_loss),
                'plus_validation_nmse':float(plus_loss),
                'source_value':delta,
                'equal_edit_budget':True,
                'used_evaluation_reference':False,
                **features,
            })
    return records


def paired_source_valuation(fit: pd.DataFrame,
                            edges: list[tuple[str, str]],
                            evidence: Mapping[str, tuple[np.ndarray, np.ndarray]],
                            *, alpha: float=1.0, max_parents: int=2,
                            boost: float=.3, seed: int=42,
                            n_folds: int=3) -> pd.DataFrame:
    """Source-ablation usefulness targets built *only* from outer fitting cells."""
    if not evidence:
        raise ValueError('GCEV requires at least one real or labeled synthetic source')
    if len(fit) < 48 or n_folds<3:
        raise ValueError('GCEV requires >=48 outer fitting cells and >=3 folds')
    if not (0 < boost <= 1):
        raise ValueError('Evidence weight must be in (0,1]')
    for name, (values, mask) in evidence.items():
        if len(values)!=len(edges) or len(mask)!=len(edges):
            raise ValueError(f'Length mismatch in evidence source {name}')
        if not np.isfinite(values).all() or (values<0).any() or (values>1).any():
            raise ValueError('Source values must be finite and within [0,1]')
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    rows = []
    for j,(a,b) in enumerate(kf.split(fit)):
        rows += _one_context(fit.iloc[a],fit.iloc[b],edges,evidence,
                             alpha=alpha,max_parents=max_parents,boost=boost,
                             split_seed=seed+51*j+13,fold=j)
    if not rows:
        raise ValueError('No eligible graph states to value: need >=2 candidate edits per target')
    result = pd.DataFrame(rows)
    if result.duplicated(['fold','target','evidence_source']).any():
        raise AssertionError('Duplicate graph evidence valuation target rows')
    return result


def learn_graph_conditional_weights(targets: pd.DataFrame,
                                    *, shrinkage: float=.05,
                                    seed: int=42) -> GCEVResult:
    """Target-held-out OOF learned weights; graph state ablation; fixed global gate.

    Each held-out target's labels are excluded from learning its weight; all
    source-value labels derive exclusively from an outer fit subset.  A global
    comparator is also target-out-of-fold. The graph state ablation receives
    identical fitting rows and excludes only GRAPH_FEATURES.
    """
    required = {'target','evidence_source','source_value', *FEATURES}
    if not required.issubset(targets.columns):
        raise ValueError('Missing GCEV training columns: '+str(sorted(required-set(targets.columns))))
    if shrinkage<=0:
        raise ValueError('Positive shrinkage required')
    groups = targets['target'].astype(str).to_numpy()
    unique = np.unique(groups)
    if len(unique)<3:
        raise ValueError('At least three distinct held-out target genes are required')
    sources = sorted(targets['evidence_source'].astype(str).unique())
    one_hot = np.column_stack([(targets['evidence_source']==s).to_numpy(float) for s in sources])
    y = targets['source_value'].to_numpy(float)
    estimated = {k:np.zeros(len(y),dtype=float) for k in ['gcev','no_graph_state','global_source']}
    cv = GroupKFold(n_splits=min(3,len(unique)))
    for a,b in cv.split(targets,y,groups):
        for method, cols in [('gcev',FEATURES),('no_graph_state',SOURCE_FEATURES)]:
            xx = np.column_stack([targets[list(cols)].to_numpy(float),one_hot])
            mu, sd = xx[a].mean(axis=0),xx[a].std(axis=0)
            sd = np.maximum(sd,1e-8)
            reg = Ridge(alpha=10.0)
            reg.fit((xx[a]-mu)/sd,y[a]);estimated[method][b] = reg.predict((xx[b]-mu)/sd)
        for source in sources:
            src_mask = (targets['evidence_source']==source).to_numpy()
            src_fit = a[src_mask[a]]
            vals = y[src_fit]
            estimated['global_source'][b[src_mask[b]]] = float(vals.mean()) if len(vals) else 0.0
    summary=targets[['target','evidence_source']].copy()
    for method,pred in estimated.items():
        summary[f'{method}_predicted_value'] = pred
        summary[f'{method}_weight'] = np.clip(pred/shrinkage,0.0,1.0)
    per_target_source=summary.groupby(['target','evidence_source'],as_index=False).agg({
        key:'mean' for key in summary if key.endswith('_weight') or key.endswith('_value')
    })
    # Include actual observed target utility only in the *diagnostic* table,
    # never pass observed labels from the held-out group to the learned gate.
    target_diagnostics=targets.copy()
    for key,pred in estimated.items():
        target_diagnostics[f'{key}_oof_predicted_value']=pred
    return GCEVResult(
        per_target_source=per_target_source,
        paired_decisions=target_diagnostics,
        gate_fit_metadata={
            'oof_by_target':True,
            'n_outer_fit_cells_are_the_only_training_input':True,
            'n_source_value_rows':len(targets),
            'n_target_groups':int(len(unique)),
            'n_sources':len(sources),
            'features_graph_conditional':list(FEATURES),
            'features_graph_ablated':list(SOURCE_FEATURES),
            'source_valuation':'paired_graph_edit_same_budget_inner_validation',
            'reliability_not_calibrated_edge_probability':True,
            'no_gold_used':True,
            'negative_effect_gate':'clip(value/shrinkage,0,1)',
            'shrinkage':shrinkage,
            'seed':seed,
        }
    )


def compose_graph_conditional_rank(expr_rank: np.ndarray,
                                   edges: list[tuple[str, str]],
                                   evidence: Mapping[str, tuple[np.ndarray,np.ndarray]],
                                   learned: pd.DataFrame,
                                   method: str,boost: float=.3) -> np.ndarray:
    """Fixed source utility, group-independent predicted gate, or GCEV gate."""
    if method not in ('gcev','no_graph_state','global_source','fixed_sources'):
        raise ValueError('Unsupported source weighting method')
    score=np.asarray(expr_rank,dtype=float).copy()
    weights = {
        (str(row.target),str(row.evidence_source)):float(getattr(row,f'{method}_weight'))
        for row in learned.itertuples(index=False)
    } if method!='fixed_sources' else {}
    for source,(values,mask) in sorted(evidence.items()):
        for i,(_,target) in enumerate(edges):
            w = 1.0 if method=='fixed_sources' else weights.get((str(target),source))
            if w is None:
                raise ValueError(f'Missing out-of-fold source weighting for {target}/{source}')
            if mask[i]:
                score[i] += boost*w*values[i]
    if not np.isfinite(score).all():
        raise ValueError('GCEV generated nonfinite ranking scores')
    return score
