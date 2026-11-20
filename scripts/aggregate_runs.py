"""Aggregate independent evaluation runs without treating edges/cells as IID replicates."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def aggregate(evaluations: list[Path], out: Path, seed=2026, draws=2000):
    if len(evaluations) < 2:
        raise ValueError('At least two independent network realization runs needed for descriptive aggregation')
    frames=[]
    names=set()
    for path in evaluations:
        if path.name in names:
            raise ValueError('Duplicate run identifier; supply independent realization IDs')
        names.add(path.name)
        if not (path/'evaluation_manifest.json').exists():
            raise ValueError(f'Missing evaluation provenance manifest: {path}')
        data=pd.read_csv(path/'metrics.csv')
        required={'method','ap','auroc','test_nmse','reference_prevalence'}
        if not required.issubset(data.columns):
            raise ValueError(f'Missing essential benchmark columns in {path}')
        if data.method.duplicated().any():
            raise ValueError(f'Duplicate method values in {path}')
        data.insert(0,'realization_id',path.name)
        frames.append(data)
    per_network=pd.concat(frames,ignore_index=True)
    methods=[set(v.method) for v in frames]
    if any(m != methods[0] for m in methods[1:]):
        raise ValueError('Cannot pool incomparable method sets across networks')
    n=len(evaluations)
    summary=[]
    rng=np.random.default_rng(seed)
    for method in sorted(methods[0]):
        x=per_network.loc[per_network.method==method].set_index('realization_id')
        for metric in ('ap','auroc','test_nmse'):
            vals=x[metric].to_numpy(float)
            sampled=vals[rng.integers(n,size=(draws,n))].mean(axis=1)
            summary.append({'method':method,'metric':metric,'n_independent_networks_reported':n,
                            'mean':float(vals.mean()),'low_95':float(np.quantile(sampled,.025)),
                            'high_95':float(np.quantile(sampled,.975)),
                            'unit':'network_realization_provided_by_user',
                            'note':'Bootstrap resamples network realizations only; topology clustering and true independence require source audit'})
    out.mkdir(parents=True,exist_ok=True)
    per_network.to_csv(out/'per_network.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'uncertainty.csv',index=False)
    (out/'aggregate_manifest.json').write_text(json.dumps({'realization_ids':sorted(names),
         'seed':seed,'bootstrap_draws':draws,'provenance_warning':'Analyst must verify independent graph realizations, non-overlapping test contexts and reference completeness; results not automatically confirmatory.'},indent=2)+'\n')
    return pd.DataFrame(summary)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('evaluations',nargs='+',type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--seed',type=int,default=2026)
    a=p.parse_args()
    print(aggregate(a.evaluations,a.output,a.seed).to_string(index=False))

if __name__=='__main__':main()
