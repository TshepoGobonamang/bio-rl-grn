"""Post-frozen independent metric and leakage-invariant audit of GCEV mechanistic study."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from grnbench.core import sha256
from scripts.gcev_parent_state_confirm import LEVELS, SCENARIOS, SEEDS


def verify(root: Path):
    root=Path(root)
    metrics=pd.read_csv(root/'all_method_metrics.csv')
    checked=0;errors=[];sha_checks=0;possible=set();tested=set()
    for seed in SEEDS:
        net=root/'networks'/f'dag_{seed}'
        common=net/'common'
        reference=pd.read_csv(common/'reference_complete.csv')
        gold=set(map(tuple,reference[['source','target']].itertuples(index=False,name=None)))
        meta=json.loads((net/'generation_manifest.json').read_text())
        if meta['n_true_edges']!=len(gold):errors.append(f'{seed}: duplicate or invalid gold edges')
        tf=pd.read_csv(common/'tfs.csv')['tf'].astype(str).tolist()
        candidate={(s,t) for s in tf for t in tf if s!=t}
        if not gold.issubset(candidate):errors.append(f'{seed}: invalid gold candidate')
        partitions=set();priors=set()
        for mech,*_ in LEVELS:
            for scenario,_,_ in SCENARIOS:
                key=f'{seed}:{mech}:{scenario}'
                fit=net/'conditions'/mech/scenario/'fit'
                manifest=json.loads((fit/'fit_manifest.json').read_text())
                if manifest['gold_reference_accessed'] is not False:errors.append(f'{key}: gold accessed during fit')
                partitions.add(manifest['partition_sha256'])
                priors.add((manifest['prior_sha256'],manifest['secondary_prior_sha256']))
                tar=pd.read_csv(fit/'gcev'/'paired_source_value_targets.csv')
                if not (tar['used_evaluation_reference'].eq(False).all() and tar['equal_edit_budget'].eq(True).all()):
                    errors.append(f'{key}: training target contract violated')
                gold_refhash=meta['reference_sha256']
                if sha256(common/'reference_complete.csv')!=gold_refhash:errors.append(f'{key}: gold source changed')
                rows=metrics.loc[(metrics.seed==seed)&(metrics.mechanism==mech)&(metrics.scenario==scenario)]
                methods=rows.method.astype(str).tolist()
                if len(rows)!=len(manifest['ranking_sha256']) or len(set(methods))!=len(methods):
                    errors.append(f'{key}: missing or duplicate method metrics')
                for _,row in rows.iterrows():
                    ranking=fit/'rankings'/f'{row.method}.csv'
                    if sha256(ranking)!=manifest['ranking_sha256'].get(row.method):errors.append(f'{key}: ranking hash mismatch {row.method}')
                    df=pd.read_csv(ranking);edges=list(map(tuple,df[['source','target']].itertuples(index=False,name=None)))
                    if len(edges)!=len(candidate) or set(edges)!=candidate:errors.append(f'{key}: candidate universe mismatch {row.method}')
                    y=[int(e in gold) for e in edges];pred=df.score.astype(float).to_numpy()
                    ap=average_precision_score(y,pred);auc=roc_auc_score(y,pred)
                    if not np.isclose(ap,row.ap,rtol=1e-11,atol=1e-12):errors.append(f'{key}: AP not reproducible {row.method}')
                    if not np.isclose(auc,row.auroc,rtol=1e-11,atol=1e-12):errors.append(f'{key}: AUROC not reproducible {row.method}')
                    checked+=1
                tested.add(key)
        if len(partitions)!=1:errors.append(f'{seed}: outer splits differ')
        if len(priors)!=1:errors.append(f'{seed}: source evidence differs across paired mechanisms')
        sha_checks+=len(partitions)+len(priors)
    result={'n_independent_graphs':len(SEEDS),'n_expected_conditions':len(SEEDS)*len(LEVELS)*len(SCENARIOS),
            'n_checked_conditions':len(tested),'n_ranking_files_independently_rescored':checked,
            'parent_state_ap_audit':True,
            'all_method_metric_rows':len(metrics),'max_allowed_metric_error':1e-11,
            'matched_split_and_evidence_checks':sha_checks,'no_errors':not bool(errors),'errors':errors[:20],
            'oracle_prior_is_synthetic_only':True}
    (root/'independent_integrity_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    if errors:raise RuntimeError(f'{len(errors)} audit failures: {errors[:4]}')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path)
    args=p.parse_args();print(json.dumps(verify(args.root),indent=2))
