"""Controlled engineering pilot only: synthetic graphs and fabricated evidence.

NEVER use these metrics as an empirical human-gene-regulatory validation.
This fixture explicitly uses simulation truth to manufacture sources and then
runs the estimator with the reference inaccessible until evaluation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from grnbench.core import candidate_pairs, evaluate, generate_smoke, run_fit


def simulate_evidence(data_dir: Path, seed: int) -> tuple[Path,Path]:
    rng=np.random.default_rng(seed + 888)
    gold_df=pd.read_csv(data_dir/'reference_complete.csv')
    gold=set(zip(gold_df.source,gold_df.target))
    genes=list(pd.read_csv(data_dir/'tfs.csv').tf)
    src_paths=[]
    for source_idx,name in enumerate(('SIM_A','SIM_B')):
        records=[]
        for s,t in candidate_pairs(genes,genes):
            if rng.random()>.77:
                continue
            target_id=genes.index(t)
            reliable = (target_id%2==source_idx)
            edge_truth=(s,t) in gold
            if not reliable: edge_truth=not edge_truth
            value=np.clip((.9 if edge_truth else .1)+rng.normal(scale=.07),0,1)
            records.append({'source':s,'target':t,'score':float(value)})
        path=data_dir/f'TOY_ONLY_{name}_gold_derived.csv'
        pd.DataFrame(records).to_csv(path,index=False)
        src_paths.append(path)
    return src_paths[0],src_paths[1]


def run_pilot(root:Path,seeds:list[int],n_cells:int=130,n_genes:int=7):
    outputs=[]
    for seed in seeds:
        base=root/f'network_{seed}'
        data=base/'data'
        generate_smoke(data,seed=seed,n_cells=n_cells,n_genes=n_genes)
        p1,p2=simulate_evidence(data,seed)
        fit=base/'fit'
        run_fit(data/'expression.csv',fit,tf_list=data/'tfs.csv',
                prior=p1,prior_secondary=p2,gcev=True,seed=seed,
                max_genes=n_genes,max_parents=2,k=min(7,n_genes+1),
                graph_query_budget=1600)
        ev=evaluate(fit,data/'reference_complete.csv',base/'evaluation',
                    reference_complete=True,k=min(7,n_genes+1))
        ev.insert(0,'seed',seed)
        outputs.append(ev)
    root.mkdir(parents=True,exist_ok=True)
    table=pd.concat(outputs,ignore_index=True)
    table.to_csv(root/'SYNTHETIC_ONLY_pilot_metrics.csv',index=False)
    contract={'status':'SYNTHETIC ENGINEERING PILOT — NOT INDEPENDENT BIOLOGICAL EVIDENCE',
              'seeds':seeds,'n_cells_per_network':n_cells,'n_genes':n_genes,
              'priors_derived_from_generation_truth':True,
              'causal_inference_claim':False,'real_biological_source':False,
              'gold_read_during_estimator_fit':False,
              'completion_claim':False}
    (root/'PILOT_CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n')
    print(table[table.method.isin(['conditional','fixed_sources','global_source','no_graph_state','gcev'])][['seed','method','ap','test_nmse']].to_string(index=False))
    print('\nSYNTHETIC ENGINEERING PILOT ONLY — NOT A PUBLICATION VALIDATION')
    return table


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seeds',type=int,nargs='+',default=[11,42,77])
    p.add_argument('--cells',type=int,default=130)
    p.add_argument('--genes',type=int,default=7)
    a=p.parse_args()
    run_pilot(a.output,a.seeds,a.cells,a.genes)
