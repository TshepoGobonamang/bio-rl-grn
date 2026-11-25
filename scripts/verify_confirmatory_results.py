"""Independent post-fit audit of every sealed ranking and synthetic reference metric.

Recomputes sklearn AP/AUROC and checks that no source/reference files mutate.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from grnbench.core import candidate_pairs, sha256


def audit(root: Path, tolerance: float=1e-11) -> dict:
    p = root/'per_network_method_metrics.csv'
    data = pd.read_csv(p)
    checks=0
    differences=[]
    failures=[]
    manifests=0
    for seed in range(2001,2049):
        network=f'net_{seed}'
        raw=root/'networks'/network/'data'
        true=pd.read_csv(raw/'reference_complete.csv')
        gold=set(map(tuple,true[['source','target']].itertuples(index=False,name=None)))
        for scenario in ['clean','dropout50','corrupt25','combined']:
            fit=root/'networks'/network/'conditions'/scenario/'fit'
            manifest=json.loads((fit/'fit_manifest.json').read_text())
            assert manifest['gold_reference_accessed'] is False
            assert (root/'networks'/network/'conditions'/scenario/'evaluation'/'evaluation_manifest.json').exists()
            if manifest['expression_sha256'] != sha256(raw/'expression.csv'):
                failures.append(f'{network}/{scenario}: expression mutated')
            if manifest['prior_sha256'] != sha256(raw/'SYNTHETIC_ONLY_SIM_A_truth_derived.csv'):
                failures.append(f'{network}/{scenario}: primary prior mutated')
            if manifest['secondary_prior_sha256'] != sha256(raw/'SYNTHETIC_ONLY_SIM_B_truth_derived.csv'):
                failures.append(f'{network}/{scenario}: secondary prior mutated')
            universe=set(candidate_pairs(manifest['genes'],manifest['regulators']))
            eligible=universe & gold
            if not (0<len(eligible)<len(universe)):
                failures.append(f'{network}/{scenario}: unscorable reference')
            subset=data[(data.network_id==network)&(data.scenario==scenario)]
            if len(subset)!=11: failures.append(f'{network}/{scenario}: incomplete method grid')
            for _,row in subset.iterrows():
                fname=fit/'rankings'/f'{row.method}.csv'
                if sha256(fname) != manifest['ranking_sha256'].get(row.method):
                    failures.append(f'{network}/{scenario}/{row.method}: ranking hash mismatch')
                    continue
                r=pd.read_csv(fname)
                ids=list(zip(r.source,r.target))
                if len(ids)!=len(universe) or set(ids)!=universe:
                    failures.append(f'{network}/{scenario}/{row.method}: candidate mismatch')
                    continue
                y=np.asarray([int(e in eligible) for e in ids])
                ap=float(average_precision_score(y,r.score));auc=float(roc_auc_score(y,r.score))
                differences.append(max(abs(ap-row.ap),abs(auc-row.auroc)))
                if abs(ap-row.ap)>tolerance or abs(auc-row.auroc)>tolerance:
                    failures.append(f'{network}/{scenario}/{row.method}: metric mismatch')
                checks+=1
            manifests+=1
    summary={'status':'PASS' if not failures and checks==2112 else 'FAIL',
             'prereg_protocol_hash':json.loads((root/'primary_decision.json').read_text())['protocol_sha256'],
             'source_runs_checked':manifests,
             'method_rankings_checked':checks,
             'max_abs_metric_recomputation_difference':max(differences) if differences else None,
             'failure_count':len(failures),'failures':failures[:60],
             'exact_reference_is_synthetic_oracle':True,
             'no_biological_accuracy_claim':True}
    (root/'independent_ranking_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    z=argparse.ArgumentParser()
    z.add_argument('--root',required=True,type=Path)
    args=z.parse_args()
    r=audit(args.root)
    print(json.dumps(r,indent=2))
    if r['status']!='PASS':raise SystemExit(1)
