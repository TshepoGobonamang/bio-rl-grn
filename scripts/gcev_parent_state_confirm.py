"""Preregistered GCEV mechanism confirmation: 32 unseen DAGs x four mechanisms x two stresses: direct parent-feature ablation.

The original synthetic evidence is deliberately oracle-derived: this is a computation
and mechanistic simulation, NOT independent biological evidence.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from grnbench.core import run_fit, evaluate, sha256
from scripts.confirmatory_gcev import generate_synthetic_evidence

ID = 'GCEV-PARENT-LOCK-2026-11-25-001'
LEVELS = [('linear',0,0),('saturation',1,0),('interaction',0,1),('combined_mechanism',1,1)]
SCENARIOS = [('clean',0.,0.),('source_stress',.5,.25)]
SEEDS = tuple(range(7201,7233))

def design(seed):
    if seed not in SEEDS: raise ValueError('Seed outside frozen confirmation set')
    i=seed-7201
    return {'network_id':f'dag_{seed}', 'seed':seed,
            'n_genes':7 if i%2==0 else 9,
            'n_cells':140 if i%2==0 else 180,
            'evidence_regime':'target_switch' if (i//2)%2==0 else 'global_stable'}

def draw_dag(seed,n_genes,n_cells):
    rng=np.random.default_rng(seed)
    weights=np.zeros((n_genes,n_genes))
    for attempt in range(1,101):
        weights[:]=0
        for target in range(1,n_genes):
            for src in range(target):
                if rng.random()<.30:
                    weights[src,target]=rng.choice([-1.,1.])*rng.uniform(.35,.80)
        if np.count_nonzero(weights)>=2: break
    else: raise RuntimeError('100 attempts with <2 edges')
    return weights,rng.standard_normal((n_cells,n_genes)),attempt

def expression_for(weights,noise,saturation,interaction):
    if saturation not in (0,1) or interaction not in (0,1): raise ValueError('Locked binary mechanisms only')
    x=np.zeros_like(noise)
    for t in range(noise.shape[1]):
        p=np.flatnonzero(weights[:,t])
        v=noise[:,t].copy()
        for s in p:
            v+=weights[s,t]*((1-saturation)*x[:,s]+saturation*1.1*np.tanh(x[:,s]))
        if interaction and len(p)>=2:
            v+=.30*np.tanh(x[:,p[0]]*x[:,p[1]])
        x[:,t]=v
    return x

def generate_pair(root,spec):
    root.mkdir(parents=True,exist_ok=True)
    n=spec['n_genes']; w,noise,attempt=draw_dag(spec['seed'],n,spec['n_cells'])
    genes=[f'g{i:02d}' for i in range(n)]
    common=root/'common';common.mkdir(exist_ok=True)
    pd.DataFrame({'tf':genes}).to_csv(common/'tfs.csv',index=False)
    refs=[{'source':genes[s],'target':genes[t]} for s,t in zip(*np.nonzero(w))]
    pd.DataFrame(refs,columns=['source','target']).to_csv(common/'reference_complete.csv',index=False)
    priors=generate_synthetic_evidence(common,spec)
    outputs={}
    for name,sat,inter in LEVELS:
        path=root/'mechanisms'/name;path.mkdir(parents=True,exist_ok=True)
        x=expression_for(w,noise,sat,inter)
        pd.DataFrame(x,columns=genes,index=[f'cell{i:04d}' for i in range(len(x))]).to_csv(path/'expression.csv',index_label='cell_id')
        outputs[name]={'saturation':sat,'interaction':inter,'expression_sha256':sha256(path/'expression.csv')}
    # Only simulation metadata: never pass generating oracle weights into grnbench.run_fit.
    multireg=int(sum(np.count_nonzero(w[:,t])>=2 for t in range(n)))
    attrs={'seed':spec['seed'],'n_true_edges':len(refs),'graph_redraw_attempts':attempt,
           'multi_parent_targets':multireg,'fraction_multi_parent_targets':multireg/n,
           'reference_sha256':sha256(common/'reference_complete.csv'),
           'source_hashes':{k:v['sha256'] for k,v in priors.items()},
           'mechanisms':outputs,'shared_exogenous_noise':True,'shared_dag_and_evidence':True}
    (root/'generation_manifest.json').write_text(json.dumps(attrs,indent=2)+'\n')
    return attrs

def one_dag(spec, output, protocol_sha):
    root=Path(output)/'networks'/spec['network_id']
    stage='generate'; rows=[]; diagnostics=[]; seals={}
    try:
        attrs=generate_pair(root,spec)
        common=root/'common'
        for name,sat,inter in LEVELS:
            for scenario,missing,corrupt in SCENARIOS:
                stage=f'{name}_{scenario}'
                dest=root/'conditions'/name/scenario
                inp=root/'mechanisms'/name/'expression.csv'
                manifest=run_fit(inp,dest/'fit',tf_list=common/'tfs.csv',
                    prior=common/'SYNTHETIC_ONLY_SIM_A_truth_derived.csv',
                    prior_secondary=common/'SYNTHETIC_ONLY_SIM_B_truth_derived.csv',
                    seed=spec['seed'],max_genes=spec['n_genes'],max_parents=2,k=7,
                    alpha=1.0,graph_query_budget=1600,gcev=True,
                    missing=missing,corrupt=corrupt)
                if manifest['gold_reference_accessed'] is not False:raise AssertionError('Fit saw gold')
                if len(manifest['ranking_sha256'])<9:raise AssertionError('Methods missing')
                if manifest['partition_sha256']!=seals.setdefault('partition',manifest['partition_sha256']):
                    raise AssertionError('Paired mechanisms used different outer splits')
                for field in ('prior_sha256','secondary_prior_sha256'):
                    if manifest[field]!=seals.setdefault(field,manifest[field]):
                        raise AssertionError('Source evidence hashes differ across paired mechanism')
                ev=evaluate(dest/'fit',common/'reference_complete.csv',dest/'eval',reference_complete=True,k=7)
                for _,m in ev.iterrows():
                    d=m.to_dict();d.update(spec);d.update({'mechanism':name,'saturation':sat,
                        'interaction':inter,'scenario':scenario,'missing':missing,'corrupt':corrupt,
                        'protocol_sha256':protocol_sha,
                        'gcev_training_seconds':manifest['gcev_training_elapsed_seconds'],
                        'gcev_nested_loss_evaluations':manifest['gcev_nested_loss_evaluations']})
                    rows.append(d)
                vals=pd.read_csv(dest/'fit'/'gcev'/'paired_source_value_targets.csv')
                diagnostics.append({'network_id':spec['network_id'],'mechanism':name,
                    'scenario':scenario,'n_paired_source_values':len(vals),
                    'parent_count_mean':float(vals.parent_count.mean()),
                    'parent_count_std':float(vals.parent_count.std(ddof=0)),
                    'parent_target_correlation_mean':float(vals.parent_target_correlation.mean()),
                    'parent_target_correlation_std':float(vals.parent_target_correlation.std(ddof=0)),
                    'parent_count_unique':int(vals.parent_count.nunique()),
                    'oof_nmse_gcev':float(((vals.source_value-vals.gcev_oof_predicted_value)**2).mean()),
                    'oof_nmse_no_parent':float(((vals.source_value-vals.no_parent_state_oof_predicted_value)**2).mean()),
                    'oof_nmse_no_action':float(((vals.source_value-vals.no_action_state_oof_predicted_value)**2).mean()),
                    'changed_edit_fraction':float(vals.different_edge.mean()),
                    'mean_abs_source_value':float(vals.source_value.abs().mean()),
                    'mean_source_disagreement':float(vals.source_disagreement.mean()),
                    'mean_parent_target_corr':float(vals.parent_target_correlation.mean()),
                    'training_seconds':float(manifest['gcev_training_elapsed_seconds']),
                    'nested_loss_evaluations':int(manifest['gcev_nested_loss_evaluations'])})
        pd.DataFrame(rows).to_csv(root/'metrics.csv',index=False)
        pd.DataFrame(diagnostics).to_csv(root/'mechanism_diagnostics.csv',index=False)
        return {'network_id':spec['network_id'],'success':True,
                'metrics':str(root/'metrics.csv'),'diagnostics':str(root/'mechanism_diagnostics.csv'),
                'generation':attrs,'partition_hash':seals['partition']}
    except Exception as e:
        return {'network_id':spec['network_id'],'success':False,'stage':stage,
                'error':str(e),'traceback':traceback.format_exc()[-3200:]}

def bootstrap(values,strata,draws=4000,seed=20261125):
    values=np.asarray(values,dtype=float)
    strata=np.asarray(strata)
    if len(values)!=32 or not np.isfinite(values).all():raise ValueError('Expected 32 finite paired DAG observations')
    rng=np.random.default_rng(seed)
    parts=[]
    for stratum in sorted(set(strata)):
        group=values[strata==stratum]
        if len(group)!=8:raise ValueError('Expected 8 DAGs per stratum')
        parts.append(group[rng.integers(0,8,size=(draws,8))].mean(axis=1))
    means=np.mean(np.column_stack(parts),axis=1)
    return {'mean':float(values.mean()),'low95':float(np.quantile(means,.025)),
            'high95':float(np.quantile(means,.975)),'unit':'independent DAG','n_dag':32}

def summarize(metrics, diag, protocol, out):
    # A single strict, balanced grid; do not silently drop a failed method or network.
    wanted = ['gcev','no_parent_state','no_action_state','no_graph_state',
              'fixed_sources','global_source','conditional','greedy']
    m = metrics.loc[metrics.method.isin(wanted)].copy()
    if len(m) != 32*4*2*len(wanted) or m.duplicated(['network_id','mechanism','scenario','method']).any():
        raise ValueError('Incomplete 32 x 4 x 2 x 8 ablation method grid')
    wide=m.pivot(index=['network_id','seed','n_genes','evidence_regime','mechanism','saturation','interaction','scenario'],
                 columns='method',values=['ap','test_nmse','auroc'])
    wide.columns=['_'.join(col) for col in wide.columns]
    wide=wide.reset_index()
    wide['dAP_parent']=wide.ap_gcev-wide.ap_no_parent_state
    wide['dAP_action_without_parent']=wide.ap_no_parent_state-wide.ap_no_graph_state
    wide['dAP_action_with_parent']=wide.ap_gcev-wide.ap_no_action_state
    wide['dAP_fixed']=wide.ap_gcev-wide.ap_fixed_sources
    wide['dNMSE_parent']=(wide.test_nmse_gcev-wide.test_nmse_no_parent_state)/wide.test_nmse_no_parent_state.clip(lower=1e-10)
    wide.to_csv(out/'per_dag_mechanism_comparisons.csv',index=False)
    clean=wide[wide.scenario=='clean'].copy()
    pairing=clean.groupby(['network_id','n_genes','evidence_regime'],as_index=False).agg(
         primary_parent_delta_ap=('dAP_parent','mean'),
         primary_parent_delta_nmse=('dNMSE_parent','mean'),
         action_delta_ap_without_parent=('dAP_action_without_parent','mean'),
         action_delta_ap_with_parent=('dAP_action_with_parent','mean'),
         gcev_vs_fixed_ap=('dAP_fixed','mean'))
    pairing['stratum']=pairing.n_genes.astype(str)+'_'+pairing.evidence_regime
    if len(pairing)!=32 or pairing.primary_parent_delta_ap.isna().any():
        raise ValueError('Primary DAG aggregation is incomplete')
    pairing.to_csv(out/'per_dag_parent_effects.csv',index=False)
    st=pairing.stratum.to_numpy()
    primary=bootstrap(pairing.primary_parent_delta_ap,st)
    nmse=bootstrap(pairing.primary_parent_delta_nmse,st)
    other={k:bootstrap(pairing[k],st) for k in
           ['action_delta_ap_without_parent','action_delta_ap_with_parent','gcev_vs_fixed_ap']}
    stress=wide.groupby('scenario').agg(mean_parent_delta_ap=('dAP_parent','mean'),
                                          mean_delta_test_nmse=('dNMSE_parent','mean')).to_dict('index')
    diag=diag.copy()
    # Diagnostics are already one set per DAG×mechanism×stress; report as descriptive (not IID).
    diag.to_csv(out/'per_condition_source_decisions.csv',index=False)
    parents={'n_distinct_parent_counts':int(diag.parent_count_unique.max()),
             'fraction_contexts_with_single_parent_count':float((diag.parent_count_unique==1).mean()),
             'mean_parent_target_correlation_std':float(diag.parent_target_correlation_std.mean()),
             'mean_value_prediction_mse_gcev':float(diag.oof_nmse_gcev.mean()),
             'mean_value_prediction_mse_no_parent':float(diag.oof_nmse_no_parent.mean()),
             'mean_value_prediction_mse_no_action':float(diag.oof_nmse_no_action.mean()),
             'mean_changed_edit_rate':float(diag.changed_edit_fraction.mean()),
             'mean_training_seconds':float(diag.training_seconds.mean())}
    verdict=bool(primary['low95']>0 and primary['mean']>=.01 and nmse['high95']<=.02)
    result={'protocol_id':ID,'independent_dags':32,'conditions':256,
            'n_method_evaluations':int(len(metrics)),
            'primary_parent_state_delta_ap':primary,
            'primary_parent_state_relative_test_nmse':nmse,
            'primary_parent_state_supported':verdict,
            'decision':'PARENT_FEATURE_ADVANTAGE_SUPPORTED_IN_SIMULATION' if verdict else 'PRIMARY_PARENT_FEATURE_ADVANTAGE_NOT_ESTABLISHED',
            'secondary_descriptive':other, 'stress_descriptive':stress,
            'parent_feature_diagnostics':parents,
            'independent_biological_validation':False,
            'scope':'Direct parent feature block only; action features still depend on fitted parent scaffold.'}
    (out/'parent_state_decision.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def run(protocol_file,out,workers):
    proto=json.loads(Path(protocol_file).read_text())
    if proto['protocol_id']!=ID or [v['seed'] for v in proto['networks']]!=list(SEEDS):
        raise ValueError('Protocol identity/seed schedule was modified')
    if proto.get('ablation_feature_contract',{}).get('parent_state') != ['parent_count','parent_target_correlation']:
        raise ValueError('Frozen parent-state feature contract has changed')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    frozen=out/'protocol.json'
    if frozen.exists() and sha256(frozen)!=sha256(protocol_file):raise ValueError('Cannot overwrite frozen protocol')
    shutil.copyfile(protocol_file,frozen);hashval=sha256(frozen)
    start=time.perf_counter();ledger=[]
    if workers==1:
        for spec in proto['networks']:
            v=one_dag(spec,str(out),hashval);ledger.append(v)
            print(len(ledger),'/',32,v['network_id'],v['success'],flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(one_dag,s,str(out),hashval) for s in proto['networks']]
            for future in as_completed(futures):
                v=future.result();ledger.append(v)
                print(len(ledger),'/',32,v['network_id'],v['success'],flush=True)
    ledger.sort(key=lambda x:x['network_id']);(out/'execution_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
    failures=[x for x in ledger if not x['success']]
    if failures:
        (out/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
        raise RuntimeError(f'{len(failures)} failed DAGs; no exclusions permitted')
    metrics=pd.concat([pd.read_csv(x['metrics']) for x in ledger],ignore_index=True)
    diag=pd.concat([pd.read_csv(x['diagnostics']) for x in ledger],ignore_index=True)
    metrics.to_csv(out/'all_method_metrics.csv',index=False)
    result=summarize(metrics,diag,proto,out);result['protocol_sha256']=hashval
    result['elapsed_seconds']=round(time.perf_counter()-start,2)
    (out/'parent_state_decision.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();run(a.protocol,a.output,a.workers)
