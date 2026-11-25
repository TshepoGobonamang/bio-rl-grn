"""Produce descriptive bulk CRISPRi responses from independently identified replicate TSVs.

This script does NOT compute direct regulatory edges, significance, or biological AP.
Use ONLY after the original ENCODE file/source and replicate audits have passed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda: f.read(1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def compute_response(config: dict) -> tuple[pd.DataFrame, dict]:
    if config.get('source_metadata_audit_passed') is not True:
        raise ValueError('Original ENCODE source metadata and per-replicate quantification audit required')
    if not config.get('annotation') or not config.get('tf'):
        raise ValueError('Explicit common annotation and TF identity required')
    rows=config.get('samples') or []
    if len(set(str(row.get('biosample_id','')) for row in rows)) != len(rows):
        raise ValueError('Duplicate biosamples: no technical pseudoreplication allowed')
    arms={arm:[r for r in rows if r.get('arm')==arm] for arm in ('perturbed','control')}
    if len(arms['perturbed']) < 2 or len(arms['control']) < 2:
        raise ValueError('At least two distinct biological biosamples per arm required')
    if len(arms['perturbed'])+len(arms['control'])!=len(rows):
        raise ValueError('Unsupported arm name')
    names=[str(r.get('replicate_number')) for r in rows]
    if any(not r.get('biosample_id') or not r.get('experiment_accession') or not r.get('file_accession') for r in rows):
        raise ValueError('Missing biosample/experiment/file provenance')
    if len({(r['arm'],r['replicate_number']) for r in rows}) != len(rows):
        raise ValueError('Repeated arm/biological replicate number')
    if any(r.get('annotation')!=config['annotation'] for r in rows):
        raise ValueError('Genome annotation mismatch')
    expected={ 'perturbed': config.get('treatment_accession'),
               'control': config.get('control_accession') }
    if any(r['experiment_accession'] != expected[r['arm']] for r in rows):
        raise ValueError('Sample accession does not match frozen treatment/control')
    col_gene=config.get('gene_column','gene_id')
    col_expr=config.get('abundance_column','TPM')
    if col_expr not in ('TPM','expected_count'):
        raise ValueError('Unsupported abundance scale; use predeclared TPM or expected_count')
    values={}
    hashes={}
    for row in rows:
        path=Path(row['path'])
        frame=pd.read_csv(path,sep='\t',dtype={col_gene:'string'})
        if col_gene not in frame or col_expr not in frame:
            raise ValueError(f'{path}: wrong quantification columns')
        if frame[col_gene].isna().any() or frame[col_gene].duplicated().any():
            raise ValueError(f'{path}: missing/duplicate gene IDs')
        x=pd.to_numeric(frame[col_expr],errors='raise')
        if not np.isfinite(x).all() or (x<0).any():
            raise ValueError(f'{path}: abundance must be finite and nonnegative')
        values[(row['arm'],row['biosample_id'])]=pd.Series(x.to_numpy(),index=frame[col_gene].astype(str).to_numpy())
        hashes[row['file_accession']]=sha256(path)
    gene_sets=[set(x.index) for x in values.values()]
    overlap=set.intersection(*gene_sets)
    if len(overlap)<2:
        raise ValueError('Insufficient exact gene overlap')
    genes=sorted(overlap)
    matrices={arm:np.vstack([values[(arm,r['biosample_id'])].reindex(genes).to_numpy()
                                for r in arms[arm]]) for arm in arms}
    x=matrices['perturbed']; y=matrices['control']; pseudo=1.0
    # Descriptive log2 expression shift, not DE testing or a complete regulatory reference.
    effect=np.log2((x.mean(axis=0)+pseudo)/(y.mean(axis=0)+pseudo))
    ctrl_mean=y.mean(axis=0)
    consistent=(np.all(x>ctrl_mean,axis=0)|np.all(x<ctrl_mean,axis=0))
    result=pd.DataFrame({'gene_id':genes, 'perturbed_mean':x.mean(axis=0),
                         'control_mean':ctrl_mean, 'log2_mean_abundance_ratio':effect,
                         'perturbed_replicate_direction_consistent':consistent,
                         'n_perturbed_biosamples':len(x),'n_control_biosamples':len(y)})
    summary={'tf':config['tf'],'treatment_accession':config['treatment_accession'],
             'control_accession':config['control_accession'], 'annotation':config['annotation'],
             'abundance_column':col_expr,'n_common_genes':len(genes),
             'n_biological_biosamples_per_arm':{k:len(v) for k,v in arms.items()},
             'file_sha256':hashes,'direct_binding_verified':False,
             'regulatory_edge_truth_labels_generated':False,
             'fdr_significance_computed':False,
             'model_predictions_unsealed':False,
             'note':'Descriptive RNA response only; no per-gene false-negative or direct TF edge claims.'}
    return result, summary


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--config',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    df,summary=compute_response(json.loads(a.config.read_text()))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(a.output,index=False)
    a.output.with_suffix('.metadata.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
