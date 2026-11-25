"""Leakage-aware expression-only and evidence-assisted GRN benchmark.

No reference-network input enters the fitting, graph selection, or adaptive
reliability methods in this module. Predictive gains are not causal effects.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold


@dataclass(frozen=True)
class Partition:
    fit: np.ndarray
    score: np.ndarray
    test: np.ndarray


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for piece in iter(lambda: f.read(1024 * 1024), b''):
            h.update(piece)
    return h.hexdigest()


def partition(n: int, seed: int = 42) -> Partition:
    if n < 30:
        raise ValueError('At least 30 cells are required for this pilot')
    indices = np.random.default_rng(seed).permutation(n)
    a, b = int(0.6 * n), int(0.8 * n)
    return Partition(indices[:a], indices[a:b], indices[b:])


def read_expression(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, index_col=0)
    if data.empty or not data.index.is_unique or not data.columns.is_unique:
        raise ValueError('Expression matrix must have unique cells and genes')
    if any(str(x).startswith('Unnamed:') for x in data.columns):
        raise ValueError('Expression matrix appears to contain an extra index column')
    values = data.apply(pd.to_numeric, errors='raise')
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError('Expression data contain nonfinite values')
    if values.shape[1] < 3:
        raise ValueError('At least three genes are required')
    return values


def candidate_pairs(genes: list[str], regulators: list[str]) -> list[tuple[str, str]]:
    return [(s, t) for t in genes for s in regulators if s != t]


def read_prior(path: Path | None, edges: list[tuple[str, str]]) -> tuple[np.ndarray, np.ndarray]:
    score = np.zeros(len(edges), dtype=float)
    mask = np.zeros(len(edges), dtype=bool)
    if path is None:
        return score, mask
    df = pd.read_csv(path)
    if not {'source', 'target', 'score'}.issubset(df):
        raise ValueError('Prior requires source,target,score')
    if df[['source', 'target']].isna().any().any() or df.duplicated(['source', 'target']).any():
        raise ValueError('Prior contains empty or duplicate directed edge identifiers')
    vals = pd.to_numeric(df['score'], errors='raise')
    if vals.isna().any() or not np.isfinite(vals).all() or not vals.between(0, 1).all():
        raise ValueError('Prior scores must be finite and within [0,1]')
    lookup = {(str(s), str(t)): float(v) for s, t, v in zip(df.source, df.target, vals)}
    for i, edge in enumerate(edges):
        if edge in lookup:
            score[i] = lookup[edge]
            mask[i] = True
    return score, mask


def evidence_sensitivity(scores: np.ndarray, mask: np.ndarray, missing: float,
                         corrupt: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if not 0 <= missing <= 1 or not 0 <= corrupt <= 1:
        raise ValueError('Dropout/corruption fractions must be in [0,1]')
    new_score, new_mask = scores.copy(), mask.copy()
    rng = np.random.default_rng(seed)
    available = np.flatnonzero(mask)
    if len(available):
        drop = rng.choice(available, size=int(round(missing * len(available))), replace=False)
        new_mask[drop] = False
        kept = np.flatnonzero(new_mask)
        affected = rng.choice(kept, size=int(round(corrupt * len(kept))), replace=False)
        # Clearly declared wrong-strength perturbation, not a biological negative.
        new_score[affected] = 1.0 - new_score[affected]
    return new_score, new_mask


def _loss(x_fit: pd.DataFrame, x_score: pd.DataFrame, target: str,
          parents: tuple[str, ...], alpha: float) -> float:
    y = x_fit[target].to_numpy(dtype=float)
    yv = x_score[target].to_numpy(dtype=float)
    if parents:
        tr, te = x_fit[list(parents)].to_numpy(float), x_score[list(parents)].to_numpy(float)
        mu, sd = tr.mean(axis=0), tr.std(axis=0)
        sd = np.maximum(sd, 1e-9)
        model = Ridge(alpha=alpha)
        model.fit((tr - mu) / sd, y)
        pred = model.predict((te - mu) / sd)
    else:
        pred = np.full_like(yv, y.mean())
    return float(np.mean((yv - pred) ** 2) / max(float(np.var(y)), 1e-9))


def _scaffold(x: pd.DataFrame, source: str, target: str, size: int) -> tuple[str, ...]:
    cols = [c for c in x.columns if c not in {source, target}]
    if not cols or size <= 0:
        return ()
    target_values = x[target].to_numpy(float)
    ranked = []
    for candidate in cols:
        v = x[candidate].to_numpy(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            c = np.corrcoef(v, target_values)[0, 1]
        ranked.append((0.0 if not np.isfinite(c) else abs(float(c)), candidate))
    ranked.sort(key=lambda z: (-z[0], z[1]))
    return tuple(name for _, name in ranked[:size])


def conditional_gains(fit: pd.DataFrame, valid: pd.DataFrame,
                      edges: list[tuple[str, str]], alpha: float,
                      max_parents: int) -> np.ndarray:
    gains = []
    for source, target in edges:
        scaffold = _scaffold(fit, source, target, max_parents)
        before = _loss(fit, valid, target, scaffold, alpha)
        after = _loss(fit, valid, target, scaffold + (source,), alpha)
        gains.append(float(np.clip((before - after) / max(before, 1e-9), -1, 1)))
    return np.asarray(gains, float)


def correlations(fit: pd.DataFrame, edges: list[tuple[str, str]],
                 method: str) -> np.ndarray:
    c = fit.corr(method=method).fillna(0.0)
    return np.asarray([abs(float(c.loc[s, t])) for s, t in edges])


def crossfit_targets(fit: pd.DataFrame, edges: list[tuple[str, str]],
                     alpha: float, max_parents: int, seed: int) -> np.ndarray:
    kf = KFold(n_splits=3, shuffle=True, random_state=seed)
    values = np.zeros((3, len(edges)), float)
    for j, (a, b) in enumerate(kf.split(fit)):
        values[j] = conditional_gains(fit.iloc[a], fit.iloc[b], edges, alpha, max_parents)
    return values.mean(axis=0)


def _rank_unit(values: np.ndarray) -> np.ndarray:
    if len(values) == 0:
        return values
    return (rankdata(values, method='average') - 1) / max(len(values) - 1, 1)


def adaptive_oof(fit: pd.DataFrame, edges: list[tuple[str, str]],
                 corr: np.ndarray, prior: np.ndarray, mask: np.ndarray,
                 alpha: float, max_parents: int, seed: int) -> np.ndarray:
    y = crossfit_targets(fit, edges, alpha, max_parents, seed)
    # Biological masks and source values are never interpreted as truth labels.
    features = np.column_stack([corr, prior * mask, mask.astype(float),
                                _rank_unit(corr) * mask.astype(float)])
    groups = np.array([target for _, target in edges])
    if len(np.unique(groups)) < 3:
        raise ValueError('Adaptive target-group crossfit needs three distinct targets')
    prediction = np.zeros(len(edges), dtype=float)
    cv = GroupKFold(n_splits=min(3, len(np.unique(groups))))
    for tr, te in cv.split(features, y, groups):
        mu, sd = features[tr].mean(axis=0), features[tr].std(axis=0)
        sd = np.maximum(sd, 1e-9)
        learned = Ridge(alpha=10.0)
        learned.fit((features[tr] - mu) / sd, y[tr])
        prediction[te] = learned.predict((features[te] - mu) / sd)
    return prediction


def greedy_ranking(fit: pd.DataFrame, score: pd.DataFrame,
                   edges: list[tuple[str, str]], alpha: float,
                   k: int, budget: int, fallback: np.ndarray):
    chosen: list[tuple[str, str]] = []
    parents = {target: tuple() for target in fit.columns}
    cached_loss = {target: _loss(fit, score, target, (), alpha) for target in fit.columns}
    trajectory = []
    queries = 0
    remaining = set(edges)
    for step in range(min(k, len(edges))):
        best, best_gain, best_new_loss = None, -float('inf'), None
        if queries + len(remaining) > budget:
            break
        for source, target in sorted(remaining):
            after = _loss(fit, score, target, parents[target] + (source,), alpha)
            queries += 1
            gain = cached_loss[target] - after
            if gain > best_gain:
                best, best_gain, best_new_loss = (source, target), gain, after
        if best is None or best_gain <= 1e-12:
            break
        remaining.remove(best)
        chosen.append(best)
        s, t = best
        parents[t] += (s,)
        cached_loss[t] = best_new_loss
        trajectory.append({'step':step+1,'source':s,'target':t,
                           'internal_gain':float(best_gain),'utility_queries':queries,
                           'score_loss_total':float(sum(cached_loss.values()))})
    leftovers = sorted(remaining, key=lambda e: (-float(fallback[edges.index(e)]), e))
    return chosen + leftovers, trajectory, queries, len(chosen)


def heldout_graph_loss(fit_and_score: pd.DataFrame, test: pd.DataFrame,
                       selected: list[tuple[str, str]], alpha: float) -> float:
    by_target = {t:[] for t in fit_and_score.columns}
    for source, target in selected:
        by_target[target].append(source)
    return float(np.mean([_loss(fit_and_score, test, target,
                                tuple(by_target[target]), alpha)
                          for target in fit_and_score.columns]))


def run_fit(expression: Path, out: Path, *, tf_list: Path | None = None,
            prior: Path | None = None, seed: int = 42, max_genes: int = 12,
            max_parents: int = 2, k: int = 10, alpha: float = 1.0,
            graph_query_budget: int = 5000, missing: float = 0.0,
            corrupt: float = 0.0, allow_all_genes: bool = False,
            gcev: bool = False, prior_secondary: Path | None = None):
    if k < 1 or max_genes < 3 or max_parents < 0 or graph_query_budget < 1:
        raise ValueError('Invalid positive model configuration')
    expression = Path(expression)
    x = read_expression(expression)
    split = partition(len(x), seed)
    if len(x.columns) > max_genes:
        # Gene universe selected using fitting cells ONLY, never references/priors.
        variance = x.iloc[split.fit].var().sort_values(ascending=False, kind='stable')
        selected_genes = variance.index[:max_genes].tolist()
        x = x[selected_genes]
    if tf_list is None:
        if not allow_all_genes:
            raise ValueError('Provide a TF whitelist; --allow-all-genes is permitted for toy smoke runs only')
        regulators = list(x.columns)
    else:
        tf = pd.read_csv(tf_list)
        if 'tf' not in tf.columns:
            raise ValueError('TF list requires tf header')
        if tf['tf'].isna().any() or tf['tf'].duplicated().any():
            raise ValueError('TF list has invalid or duplicate gene names')
        regulators = [g for g in x.columns if g in set(tf['tf'].astype(str))]
        if not regulators:
            raise ValueError('No TF whitelist genes matched expression genes')
    fit, score, test = x.iloc[split.fit], x.iloc[split.score], x.iloc[split.test]
    edges = candidate_pairs(list(x.columns), regulators)
    if not edges:
        raise ValueError('No eligible directed candidate edges')
    evidence, mask = read_prior(prior, edges)
    evidence, mask = evidence_sensitivity(evidence, mask, missing, corrupt, seed + 500)
    if (missing or corrupt) and prior is None:
        raise ValueError('Evidence corruption/dropout requires a genuine explicitly provided prior file')
    pearson = correlations(fit, edges, 'pearson')
    spearman = correlations(fit, edges, 'spearman')
    conditional = conditional_gains(fit, score, edges, alpha, max_parents)
    expr_rank = _rank_unit(conditional)
    rankings = {'pearson':pearson,'spearman':spearman,'conditional':conditional}
    if mask.any():
        fixed = expr_rank + 0.3 * mask * evidence
        heuristic = expr_rank + 0.3 * mask * evidence * np.maximum(conditional, 0)
        adaptive = expr_rank + 0.3 * _rank_unit(adaptive_oof(fit, edges, pearson,
                                                            evidence, mask, alpha,
                                                            max_parents, seed))
        rankings.update(fixed=fixed, heuristic=heuristic, adaptive=adaptive)
    # GCEV cross-fitting and counterfactual source-ablation targets see ONLY
    # the outer fit subset; graph-scoring and outer test cells are not accessed.
    gcev_result = None
    gcev_elapsed_seconds = None
    gcev_fit_loss_evaluations = None
    source_files = {}
    if gcev:
        if prior is None:
            raise ValueError('GCEV requires a provided (real or explicitly synthetic) evidence source')
        from .gcev import (paired_source_valuation, learn_graph_conditional_weights,
                           compose_graph_conditional_rank)
        source_files['primary'] = str(prior)
        source_matrices = {'primary':(evidence,mask)}
        if prior_secondary is not None:
            secondary, secondary_mask = read_prior(prior_secondary, edges)
            secondary, secondary_mask = evidence_sensitivity(secondary, secondary_mask,
                                                              missing, corrupt, seed+501)
            source_matrices['secondary'] = (secondary,secondary_mask)
            source_files['secondary'] = str(prior_secondary)
        gcev_start = time.perf_counter()
        labels = paired_source_valuation(fit,edges,source_matrices,
                                         alpha=alpha,max_parents=max_parents,
                                         boost=.3,seed=seed)
        gcev_result = learn_graph_conditional_weights(labels,shrinkage=.05,seed=seed)
        gcev_elapsed_seconds = round(time.perf_counter()-gcev_start,6)
        # Every source uses two held-out graph losses. Each target/fold graph
        # scaffold needs one baseline loss and one loss per eligible action.
        unique_contexts = labels.drop_duplicates(['fold','target'])
        gcev_fit_loss_evaluations = int((unique_contexts.candidate_count+1).sum()+2*len(labels))
        for method in ('fixed_sources','global_source','no_graph_state','no_parent_state','no_action_state','gcev'):
            rankings[method] = compose_graph_conditional_rank(expr_rank,edges,
                         source_matrices,gcev_result.per_target_source,method,boost=.3)
    elif prior_secondary is not None:
        raise ValueError('A secondary source requires --gcev')
    order, trajectory, queries, accepted = greedy_ranking(fit, score, edges, alpha,
                                                          k, graph_query_budget, conditional)
    # Greedy order preserves actual edit selections. Remaining edges are fallback
    # rankings only and must not be described as accepted graph edits.
    out.mkdir(parents=True, exist_ok=True)
    (out/'rankings').mkdir(exist_ok=True)
    if gcev_result is not None:
        # Diagnostic labels are kept apart from predictor scores. Evaluation
        # cannot ingest them; ranking seal includes GCEV methods below.
        gcev_dir=out/'gcev'
        gcev_dir.mkdir(exist_ok=True)
        gcev_result.paired_decisions.to_csv(gcev_dir/'paired_source_value_targets.csv',index=False)
        gcev_result.per_target_source.to_csv(gcev_dir/'oof_source_weights.csv',index=False)
        (gcev_dir/'training_contract.json').write_text(json.dumps(
            gcev_result.gate_fit_metadata,indent=2)+'\n')
    rows = []
    for method, values in rankings.items():
        ids = sorted(range(len(edges)), key=lambda i: (-float(values[i]), edges[i]))
        ordered = [edges[i] for i in ids]
        pd.DataFrame({'source':[e[0] for e in ordered],
                      'target':[e[1] for e in ordered],
                      'score':[float(values[i]) for i in ids]}).to_csv(out/'rankings'/f'{method}.csv', index=False)
        heldout = heldout_graph_loss(x.iloc[np.concatenate((split.fit, split.score))],
                                     test, ordered[:k], alpha)
        rows.append({'method':method,'test_nmse':heldout,'top_k':min(k,len(edges)),
                     'selected_edges':min(k,len(edges)),'utility_queries':0})
    # Explicit greedy rank position; scores are monotonically descending rank,
    # not comparable to regression utility gains.
    pd.DataFrame({'source':[e[0] for e in order], 'target':[e[1] for e in order],
                  'score':np.arange(len(order),0,-1)}).to_csv(out/'rankings'/'greedy.csv', index=False)
    rows.append({'method':'greedy',
                 'test_nmse':heldout_graph_loss(x.iloc[np.concatenate((split.fit, split.score))],
                                               test, order[:k], alpha),
                 'top_k':min(k,len(edges)),'selected_edges':accepted,'utility_queries':queries})
    pd.DataFrame(rows).to_csv(out/'expression_test_results.csv', index=False)
    pd.DataFrame(trajectory).to_csv(out/'greedy_trajectory.csv', index=False)
    manifest = {'status':'EVALUATION_PENDING','source_expression':str(expression),
                'expression_sha256':sha256(expression),
                'prior_sha256':sha256(prior) if prior is not None else None,
                'secondary_prior_sha256':sha256(prior_secondary) if prior_secondary else None,
                'gcev_enabled':bool(gcev),
                'gcev_training_elapsed_seconds':gcev_elapsed_seconds,
                'gcev_nested_loss_evaluations':gcev_fit_loss_evaluations,
                'gcev_models_have_unequal_training_compute_vs_fixed':bool(gcev),
                'gcev_source_files':source_files,
                'gcev_methods':['fixed_sources','global_source','no_graph_state','no_parent_state','no_action_state','gcev'] if gcev else [],
                'gcev_training_contract_sha256':sha256(out/'gcev'/'training_contract.json') if gcev else None,
                'gcev_paired_targets_sha256':sha256(out/'gcev'/'paired_source_value_targets.csv') if gcev else None,
                'gcev_weights_sha256':sha256(out/'gcev'/'oof_source_weights.csv') if gcev else None,
                'tf_list_sha256':sha256(tf_list) if tf_list else None,
                'n_cells':len(x),'n_genes':len(x.columns),'genes':list(x.columns),
                'n_regulators':len(regulators),'regulators':regulators,
                'n_candidate_edges':len(edges),'n_prior_present':int(mask.sum()),
                'seed':seed,'alpha':alpha,'k':k,'max_parents':max_parents,
                'graph_query_budget':graph_query_budget,'prior_missing_fraction':missing,
                'prior_corruption_fraction':corrupt,
                'fit_cells':len(split.fit),'internal_scoring_cells':len(split.score),
                'untouched_test_cells':len(split.test),
                'partition_sha256':hashlib.sha256(np.concatenate((split.fit,split.score,split.test)).tobytes()).hexdigest(),
                'gold_reference_accessed':False,
                'ranking_sha256':{p.stem:sha256(p) for p in sorted((out/'rankings').glob('*.csv'))},
                'limitation':'Predictive association not causation; selected model comparisons are exploratory until independent networks and complete evaluation reference are verified.'}
    (out/'fit_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def evaluate(fit_dir: Path, reference: Path, out: Path, *, reference_complete: bool = False,
             k: int = 10):
    if not reference_complete:
        raise ValueError('AP/AUROC require an explicitly asserted complete benchmark reference; absence from an incomplete biological reference is not a negative')
    from sklearn.metrics import average_precision_score, roc_auc_score
    fit_dir, reference, out = Path(fit_dir), Path(reference), Path(out)
    manifest = json.loads((fit_dir/'fit_manifest.json').read_text())
    if manifest.get('gcev_enabled'):
        checks = [
            ('paired_source_value_targets.csv','gcev_paired_targets_sha256'),
            ('oof_source_weights.csv','gcev_weights_sha256'),
            ('training_contract.json','gcev_training_contract_sha256'),
        ]
        for name,key in checks:
            p = fit_dir/'gcev'/name
            if not p.exists() or sha256(p) != manifest.get(key):
                raise ValueError(f'GCEV artifact mutated after fitting: {name}')
    reference_df = pd.read_csv(reference)
    if not {'source','target'}.issubset(reference_df):
        raise ValueError('Reference requires source,target')
    if reference_df[['source','target']].isna().any().any():
        raise ValueError('Reference contains missing edge identifiers')
    gold = set(zip(reference_df.source.astype(str),reference_df.target.astype(str)))
    candidate = set(candidate_pairs(manifest['genes'],manifest['regulators']))
    gold_eligible = gold & candidate
    if not gold_eligible or len(gold_eligible) == len(candidate):
        raise ValueError('AUROC requires positive and negative eligible references')
    rows = []
    for path in sorted((fit_dir/'rankings').glob('*.csv')):
        if sha256(path) != manifest['ranking_sha256'].get(path.stem):
            raise ValueError(f'Ranking mutated after fitting: {path}')
        ranked = pd.read_csv(path)
        edges = list(zip(ranked.source.astype(str),ranked.target.astype(str)))
        if len(edges)!=len(candidate) or len(set(edges))!=len(candidate) or set(edges)!=candidate:
            raise ValueError(f'Candidate universe mismatch for {path.stem}')
        y = np.asarray([int(e in gold_eligible) for e in edges])
        v = ranked.score.to_numpy(dtype=float)
        if not np.isfinite(v).all():
            raise ValueError('Non-finite ranking scores')
        # Per-pair orientation: only for one-sided gold where both arrows eligible.
        scores = dict(zip(edges,v))
        unique = [e for e in gold_eligible if (e[1],e[0]) in candidate and
                  (e[1],e[0]) not in gold_eligible]
        oriented = [int(scores[e]>scores[(e[1],e[0])]) for e in unique
                    if scores[e]!=scores[(e[1],e[0])]]
        rows.append({'method':path.stem,'ap':float(average_precision_score(y,v)),
                     'auroc':float(roc_auc_score(y,v)),
                     'precision_at_k':float(y[:k].mean()),
                     'true_positives_at_k':int(y[:k].sum()),
                     'eligible_gold_edges':len(gold_eligible),'candidate_edges':len(candidate),
                     'reference_prevalence':len(gold_eligible)/len(candidate),
                     'orientation_correct':sum(oriented),
                     'orientation_assessable_untied':len(oriented),
                     'orientation_tied':len(unique)-len(oriented)})
    if not rows:
        raise ValueError('Missing sealed rankings')
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows).merge(pd.read_csv(fit_dir/'expression_test_results.csv'),on='method',validate='one_to_one')
    df.to_csv(out/'metrics.csv', index=False)
    (out/'evaluation_manifest.json').write_text(json.dumps({
        'fit_manifest_sha256':sha256(fit_dir/'fit_manifest.json'),
        'reference_sha256':sha256(reference),'reference_complete_explicit':True,
        'candidate_edges':len(candidate),'eligible_reference_edges':len(gold_eligible),
        'evaluation_uses_gold_only_after_fit':True,
        'metric':'sklearn noninterpolated average_precision_score (AP), not trapezoidal AUPRC'
    },indent=2)+'\n')
    return df


def generate_smoke(out: Path, seed:int=42, n_cells:int=180, n_genes:int=8):
    if n_cells<30 or n_genes<3:
        raise ValueError('Smoke requires >=30 cells and >=3 genes')
    rng=np.random.default_rng(seed)
    w=np.zeros((n_genes,n_genes))
    for target in range(1,n_genes):
        for source in range(target):
            if rng.random()<0.28:
                w[source,target]=rng.choice([-1.,1.])*rng.uniform(0.3,0.8)
    x=np.zeros((n_cells,n_genes))
    for target in range(n_genes):
        x[:,target]=x@w[:,target]+rng.standard_normal(n_cells)
    out.mkdir(parents=True,exist_ok=True)
    genes=[f'g{i:02d}' for i in range(n_genes)]
    expression=out/'expression.csv'
    pd.DataFrame(x,index=[f'cell{i:04d}' for i in range(n_cells)],columns=genes).to_csv(expression,index_label='cell_id')
    pd.DataFrame({'tf':genes}).to_csv(out/'tfs.csv',index=False)
    gold=[{'source':genes[i],'target':genes[j]} for i,j in zip(*np.nonzero(w))]
    if not gold:
        gold=[{'source':genes[0],'target':genes[1]}]
    pd.DataFrame(gold).to_csv(out/'reference_complete.csv',index=False)
    # Artificial modality, independent random values; never biological evidence.
    toy_prior=[]
    for target in genes:
        for source in genes:
            if source!=target and rng.random()<0.35:
                toy_prior.append({'source':source,'target':target,'score':rng.random()})
    pd.DataFrame(toy_prior).to_csv(out/'TOY_ONLY_random_prior.csv',index=False)
    (out/'DATA_WARNING.txt').write_text('SMOKE FIXTURE ONLY: synthetic expression, fully known generated graph and random fake prior; NOT experimental biology or publishable independent validation.\n')
    return expression
