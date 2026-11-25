"""Pre-registered 48-network x four-perturbation GCEV benchmark (synthetic only).

Protocol was committed to GitHub BEFORE any confirmation seed was evaluated.
This command does not alter the preregistration and never selects favorable networks.
"""
from __future__ import annotations

import argparse
import json
import shutil
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from grnbench.core import candidate_pairs, evaluate, run_fit, sha256


def generate_network(path: Path, spec: dict) -> dict:
    """Generate fully known DAG, including correctly labeled nonlinear edges."""
    seed = int(spec['seed'])
    family = str(spec['family'])
    n_genes, n_cells = int(spec['n_genes']), int(spec['n_cells'])
    if family not in ('linear', 'nonlinear') or n_genes not in (7, 9):
        raise ValueError('Network family or size violates preregistration')
    rng = np.random.default_rng(seed)
    w = np.zeros((n_genes, n_genes))
    for attempt in range(1, 101):
        w.fill(0)
        for target in range(1, n_genes):
            for source in range(target):
                if rng.random() < .30:
                    w[source, target] = rng.choice([-1., 1.]) * rng.uniform(.35, .80)
        if np.count_nonzero(w) >= 2:
            break
    else:
        raise RuntimeError('Could not obtain >=2 edges in 100 predeclared redraws')
    expr = np.zeros((n_cells, n_genes), dtype=float)
    for target in range(n_genes):
        parents = np.flatnonzero(w[:, target])
        effect = np.zeros(n_cells)
        if family == 'linear':
            for source in parents:
                effect += w[source, target] * expr[:, source]
        else:
            for source in parents:
                effect += 1.1 * w[source, target] * np.tanh(expr[:, source])
            if len(parents) >= 2:
                effect += .30 * np.tanh(expr[:, parents[0]] * expr[:, parents[1]])
        expr[:, target] = effect + rng.normal(size=n_cells)
    path.mkdir(parents=True, exist_ok=True)
    genes = [f'g{i:02d}' for i in range(n_genes)]
    pd.DataFrame(expr, index=[f'cell{i:04d}' for i in range(n_cells)], columns=genes).to_csv(path / 'expression.csv', index_label='cell_id')
    pd.DataFrame({'tf': genes}).to_csv(path / 'tfs.csv', index=False)
    reference = [{'source': genes[i], 'target': genes[j]} for i,j in zip(*np.nonzero(w))]
    pd.DataFrame(reference).to_csv(path / 'reference_complete.csv', index=False)
    metadata = {'network_id': spec['network_id'], 'seed': seed, 'family': family,
                'n_genes': n_genes, 'n_cells': n_cells,
                'n_truth_edges': len(reference), 'graph_redraw_attempts': attempt,
                'reference_complete': True, 'source_truth_used_only_to_manufacture_synthetic_inputs': True}
    (path / 'generation_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return metadata


def generate_synthetic_evidence(path: Path, spec: dict) -> dict:
    """Simulator-only oracle-derived source CSVs: NOT external biological evidence."""
    genes = list(pd.read_csv(path / 'tfs.csv')['tf'].astype(str))
    gold = set(map(tuple, pd.read_csv(path / 'reference_complete.csv')[["source", "target"]].itertuples(index=False, name=None)))
    rng = np.random.default_rng(int(spec['seed']) + 888)
    regime = spec['evidence_regime']
    sources = {}
    for index, name in enumerate(('SIM_A', 'SIM_B')):
        rows = []
        for source, target in candidate_pairs(genes, genes):
            if rng.random() > .77:
                continue
            original_truth = (source, target) in gold
            if regime == 'target_switch':
                reliable = (genes.index(target) % 2 == index)
                effective = original_truth if reliable else not original_truth
            elif regime == 'global_stable':
                flip_probability = .12 if index == 0 else .25
                effective = original_truth if rng.random() > flip_probability else not original_truth
            else:
                raise ValueError('Unknown locked evidence regime')
            raw = (.90 if effective else .10) + rng.normal(scale=.07)
            rows.append({'source': source, 'target': target, 'score': float(np.clip(raw,0,1))})
        f = path / f'SYNTHETIC_ONLY_{name}_truth_derived.csv'
        pd.DataFrame(rows, columns=['source','target','score']).to_csv(f, index=False)
        sources[name] = {'path': str(f), 'sha256': sha256(f), 'n_records': len(rows)}
    return sources


def one_network(spec: dict, settings: list[dict], root_str: str, protocol_sha: str) -> dict:
    root = Path(root_str)
    net = root / 'networks' / spec['network_id']
    raw = net / 'data'
    stage = 'generate'
    try:
        generation = generate_network(raw, spec)
        prior = generate_synthetic_evidence(raw, spec)
        rows = []
        for setting in settings:
            stage = 'fit_' + setting['id']
            dest = net / 'conditions' / setting['id']
            fit = dest / 'fit'
            run_fit(raw / 'expression.csv', fit,
                    tf_list=raw / 'tfs.csv',
                    prior=Path(prior['SIM_A']['path']),
                    prior_secondary=Path(prior['SIM_B']['path']),
                    seed=int(spec['seed']), max_genes=int(spec['n_genes']),
                    max_parents=2, k=7, alpha=1.0,
                    graph_query_budget=1600, gcev=True,
                    missing=float(setting['missing']), corrupt=float(setting['corrupt']))
            stage = 'evaluate_' + setting['id']
            ev = evaluate(fit, raw / 'reference_complete.csv', dest / 'evaluation',
                          reference_complete=True, k=7)
            fit_meta = json.loads((fit / 'fit_manifest.json').read_text())
            if fit_meta['gold_reference_accessed'] is not False:
                raise AssertionError('Gold was used during fitting')
            if len(ev) < 9 or not {'gcev','no_graph_state','fixed_sources','global_source'}.issubset(set(ev.method)):
                raise AssertionError('Missing mandatory matched comparison method')
            for _, metric in ev.iterrows():
                d = metric.to_dict()
                for key in ['network_id','seed','family','n_genes','n_cells','evidence_regime']:
                    d[key] = spec[key]
                d['scenario'] = setting['id']
                d['missing'] = setting['missing']
                d['corrupt'] = setting['corrupt']
                d['prereg_sha256'] = protocol_sha
                d['gcev_training_seconds'] = fit_meta.get('gcev_training_elapsed_seconds')
                d['gcev_nested_loss_evaluations'] = fit_meta.get('gcev_nested_loss_evaluations')
                d['fit_manifest_sha256'] = sha256(fit / 'fit_manifest.json')
                rows.append(d)
        pd.DataFrame(rows).to_csv(net/'network_metrics.csv', index=False)
        return {'network_id':spec['network_id'], 'success':True,'metrics_path':str(net/'network_metrics.csv'),
                'generation':generation, 'source_hashes':prior}
    except Exception as e:
        return {'network_id':spec['network_id'], 'success':False, 'stage':stage,
                'error':str(e), 'traceback':traceback.format_exc()[-3200:]}


def stratified_bootstrap(paired: pd.DataFrame, metric: str, draws: int, seed: int):
    rng = np.random.default_rng(seed)
    strata = [g[metric].to_numpy(dtype=float) for _,g in paired.groupby(['family','n_genes','evidence_regime'], sort=True)]
    if len(strata) != 8 or any(len(group)!=6 for group in strata):
        raise ValueError('Expected eight preregistered strata, 6 independent DAGs each')
    sampled = np.column_stack([g[rng.integers(len(g),size=(draws,len(g)))].mean(axis=1) for g in strata]).mean(axis=1)
    values = paired[metric].to_numpy(dtype=float)
    return {'mean':float(values.mean()), 'ci95_low':float(np.quantile(sampled,.025)),
            'ci95_high':float(np.quantile(sampled,.975)), 'n_networks':len(values),
            'bootstrap_draws':draws}


def summarize(metrics: pd.DataFrame, protocol: dict, output: Path):
    selected = ['gcev','no_graph_state','fixed_sources','global_source','conditional','greedy']
    p = metrics.loc[metrics.method.isin(selected)].copy()
    expected = len(protocol['network_design']) * len(protocol['scenarios']) * len(selected)
    if len(p) != expected or p.duplicated(['network_id','scenario','method']).any():
        raise ValueError(f'Incomplete method grid: {len(p)} vs {expected}')
    wide = p.pivot(index=['network_id','family','n_genes','evidence_regime','scenario'],columns='method',values=['ap','auroc','test_nmse','precision_at_k'])
    out = wide.reset_index()
    out.columns=['_'.join(map(str,c)).strip('_') for c in out.columns.to_flat_index()]
    for baseline in ['no_graph_state','fixed_sources','global_source','conditional','greedy']:
        out[f'delta_ap_vs_{baseline}'] = out['ap_gcev']-out[f'ap_{baseline}']
    out['delta_relative_test_nmse_vs_no_graph_state'] = (out['test_nmse_gcev']-out['test_nmse_no_graph_state'])/out['test_nmse_no_graph_state'].clip(lower=1e-12)
    out.to_csv(output/'paired_network_effects.csv',index=False)
    bootstrap=[]
    for scenario in [x['id'] for x in protocol['scenarios']]:
        condition=out.loc[out.scenario==scenario].copy()
        for col in [k for k in condition.columns if k.startswith('delta_')]:
            b=stratified_bootstrap(condition,col,draws=4000,seed=20261125)
            bootstrap.append({'scenario':scenario,'comparison':col,**b})
    pd.DataFrame(bootstrap).to_csv(output/'network_bootstrap_95.csv',index=False)
    primary = out.loc[out.scenario=='clean']
    ap=stratified_bootstrap(primary,'delta_ap_vs_no_graph_state',draws=4000,seed=20261125)
    nmse=stratified_bootstrap(primary,'delta_relative_test_nmse_vs_no_graph_state',draws=4000,seed=20261125)
    decision = bool(ap['ci95_low'] > 0 and nmse['ci95_high'] <= .02)
    method_means=(metrics.groupby(['scenario','method'],as_index=False)
                          .agg(mean_ap=('ap','mean'), mean_auroc=('auroc','mean'),
                               mean_test_nmse=('test_nmse','mean'), n_networks=('network_id','nunique')))
    method_means.to_csv(output/'method_means.csv',index=False)
    result={'protocol_id':protocol['protocol_id'], 'status':'COMPLETED_SYNTHETIC_ONLY_NOT_BIOLOGICAL',
            'n_independent_network_realizations':len(protocol['network_design']),
            'n_repeated_condition_runs':len(protocol['network_design'])*len(protocol['scenarios']),
            'primary_ap':ap, 'primary_relative_nmse':nmse,
            'preregistered_primary_success':decision,
            'decision':'PRIMARY_CRITERION_MET' if decision else 'PRIMARY_CRITERION_NOT_MET',
            'interpretation':'Even success only establishes these two simulated DAG families with oracle-derived evidence, not independent biology or RL.',
            'no_network_selected_posthoc':True}
    (output/'primary_decision.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def validate_protocol(protocol: dict):
    if protocol['protocol_id'] != 'GCEV-LCK-2026-11-25-001':
        raise ValueError('Wrong locked preregistration')
    seeds=[int(n['seed']) for n in protocol['network_design']]
    if seeds != list(range(2001,2049)):
        raise ValueError('Seed schedule differs from the frozen protocol')
    scenario=[(x['id'],x['missing'],x['corrupt']) for x in protocol['scenarios']]
    if scenario != [('clean',0.0,0.0),('dropout50',.5,0.0),('corrupt25',0.0,.25),('combined',.5,.25)]:
        raise ValueError('Perturbation design differs from frozen protocol')
    return True


def run(protocol_file: Path, output: Path, workers: int=4):
    protocol=json.loads(protocol_file.read_text())
    validate_protocol(protocol)
    if sha256(protocol_file) != '556c6745948cc89b10c9b7ee92a64d4b5637bcefbe180c0c062082b7d2f0272c':
        raise ValueError('Preregistered protocol file bytes differ from original frozen SHA256')
    output.mkdir(parents=True,exist_ok=True)
    frozen=output/'frozen_protocol.json'
    if frozen.exists() and sha256(frozen) != sha256(protocol_file):
        raise ValueError('Cannot resume with modified protocol')
    shutil.copyfile(protocol_file,frozen)
    protocol_sha=sha256(frozen)
    start=time.perf_counter()
    jobs=[]
    if workers == 1:
        for spec in protocol['network_design']:
            item=one_network(spec, protocol['scenarios'], str(output), protocol_sha)
            jobs.append(item)
            print(f"{item['network_id']}: {'OK' if item['success'] else 'FAILED'}",flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            tasks=[pool.submit(one_network,spec,protocol['scenarios'],str(output),protocol_sha) for spec in protocol['network_design']]
            for future in as_completed(tasks):
                item=future.result()
                jobs.append(item)
                print(f"{len(jobs)}/{len(tasks)} {item['network_id']}: {'OK' if item['success'] else 'FAILED'}",flush=True)
    jobs.sort(key=lambda x:x['network_id'])
    (output/'network_execution_ledger.json').write_text(json.dumps(jobs,indent=2)+'\n')
    bad=[item for item in jobs if not item['success']]
    if bad:
        (output/'failure_ledger.json').write_text(json.dumps(bad,indent=2)+'\n')
        raise RuntimeError(f'{len(bad)} network realizations failed; see failure_ledger; no cherry-picking permitted')
    all_metrics=pd.concat([pd.read_csv(x['metrics_path']) for x in jobs],ignore_index=True)
    all_metrics.to_csv(output/'per_network_method_metrics.csv',index=False)
    result=summarize(all_metrics,protocol,output)
    result['protocol_sha256']=protocol_sha
    result['elapsed_wall_seconds']=round(time.perf_counter()-start,3)
    result['n_failed_networks']=0
    result['source_ref_accessed_only_post_fit']=True
    (output/'primary_decision.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CONFIRMATORY_SYNTHETIC_RESULTS',json.dumps(result,sort_keys=True),flush=True)
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    opts=parser.parse_args()
    run(opts.protocol,opts.output,opts.workers)
