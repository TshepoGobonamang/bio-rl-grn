"""Prepare deterministic outer fit/internal score/untouched test splits from original K562 matrix.

Select features ONLY on outer fitting cells. This tool does not train a model
or access perturbation outcome data, regulatory prior sources, or gold edges.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()

def prepare(x, cells, genes, *, seed=20261126, n_features=512, min_fraction=0.05):
    if not sparse.issparse(x):
        raise ValueError("Original K562 sparse expression matrix required")
    x=x.tocsr()
    if x.shape != (len(cells),len(genes)):
        raise ValueError("Original matrix and identifier lengths mismatch")
    if len(set(cells))!=len(cells) or len(set(genes))!=len(genes):
        raise ValueError("Duplicate original cell or gene identifiers")
    if len(cells)<100 or n_features<=0:
        raise ValueError("Insufficient source cells or target features")
    if not np.isfinite(x.data).all() or (x.data<0).any():
        raise ValueError("Original matrix has nonfinite/negative source values")
    rng=np.random.default_rng(seed)
    indices=rng.permutation(len(cells))
    n_fit=int(0.60*len(cells))
    n_score=int(0.20*len(cells))
    folds={
        "fit":indices[:n_fit],
        "internal_score":indices[n_fit:n_fit+n_score],
        "untouched_test":indices[n_fit+n_score:],
    }
    fit=x[folds["fit"],:].copy()
    prevalence=np.asarray(fit.getnnz(axis=0)).ravel()/len(folds["fit"])
    fit.data=np.log1p(fit.data)
    mean=np.asarray(fit.mean(axis=0)).ravel()
    squared=fit.copy();squared.data**=2
    variance=np.maximum(0,np.asarray(squared.mean(axis=0)).ravel()-mean**2)
    eligible=np.flatnonzero((prevalence>=min_fraction)&(variance>0))
    if len(eligible)<n_features:
        raise ValueError("Insufficient genes meeting fitting-cell-only detection/variance criteria")
    symbols=np.asarray(genes,dtype=str)
    sort=np.lexsort((symbols[eligible],-variance[eligible]))
    selected=eligible[sort[:n_features]]
    partition=np.full(len(cells),"",dtype=object)
    for label, ix in folds.items():partition[ix]=label
    if not partition.all():
        raise ValueError("Incomplete cell split assignment")
    manifest={
        "stage":"FIT_ONLY_K562_GENE_SELECTION_NO_BIOLOGICAL_PREDICTIONS",
        "seed":seed,"n_cells":len(cells),"n_source_genes":len(genes),
        "n_selected_genes":n_features,
        "fit_cells":len(folds["fit"]),
        "internal_score_cells":len(folds["internal_score"]),
        "untouched_test_cells":len(folds["untouched_test"]),
        "minimum_fit_detected_fraction":min_fraction,
        "feature_selection":"fit-only log1p original-scale variance; stable symbol tie-break",
        "source_expression_scale":"original count-like X; true source raw-count processing not separately certified",
        "no_perturbation_response_or_gold_accessed":True,
        "source_prior_verified":False,"gcev_model_fitted":False,
    }
    return x[:,selected].tocsr(),pd.DataFrame({"cell_id":cells,"partition":partition}),pd.DataFrame({"gene_symbol":symbols[selected]}),manifest

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--seed",type=int,default=20261126)
    p.add_argument("--n-features",type=int,default=512)
    a=p.parse_args()
    source=json.loads((a.input/"k562_source_audit.json").read_text())
    if source.get("source_name")!="GSM9224207_Fresh_K562.h5ad" or not source.get("all_nonnegative_and_finite"):
        raise ValueError("Original K562-only source audit does not qualify")
    for name,sha in source["export_hashes"].items():
        if digest(a.input/name)!=sha:raise ValueError("Source export checksum mismatch: "+name)
    x=sparse.load_npz(a.input/"k562_expression_original_scale.npz")
    cells=pd.read_csv(a.input/"k562_cell_ids.csv",dtype="string")["cell_id"].tolist()
    genes=pd.read_csv(a.input/"k562_gene_ids.csv",dtype="string")["gene_id"].tolist()
    reduced,splits,selected,manifest=prepare(x,cells,genes,seed=a.seed,n_features=a.n_features)
    a.output.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(a.output/"k562_selected_original_scale.npz",reduced)
    splits.to_csv(a.output/"cell_partitions.csv",index=False)
    selected.to_csv(a.output/"selected_gene_symbols.csv",index=False)
    manifest["original_H5AD_sha256"]=source["source_sha256"]
    manifest["source_audit_sha256"]=digest(a.input/"k562_source_audit.json")
    manifest["export_sha256"]={
        p.name:digest(p) for p in
        [a.output/"k562_selected_original_scale.npz",
         a.output/"cell_partitions.csv",a.output/"selected_gene_symbols.csv"]
    }
    (a.output/"split_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(manifest,indent=2))

if __name__=="__main__":main()
