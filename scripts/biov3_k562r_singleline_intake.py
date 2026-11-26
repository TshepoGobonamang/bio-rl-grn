"""Qualify single-cell-line 10x K562-r filtered matrix without ENCODE outcomes.

K562-r (ATCC CRL-3344; fixed RNA) is NOT asserted identical to the
K562 culture used for ENCODE CRISPRi responses. This is training-candidate
qualification, not independent biological confirmation.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

SOURCE_PAGE = 'https://www.10xgenomics.com/datasets/10k-human-k562-r-cells-singleplex-sample-1-standard'


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def qualify_matrix(raw: csc_matrix, gene_ids: list[str], gene_names: list[str], cell_ids: list[str],
                   n_cells=800, n_genes=256, seed=20261126):
    if not isinstance(raw, csc_matrix): raw=raw.tocsc()
    if raw.shape!=(len(gene_ids),len(cell_ids)):
        raise ValueError('Incompatible 10x matrix dimensions')
    if len(set(cell_ids))!=len(cell_ids) or not all(cell_ids):
        raise ValueError('Duplicate or missing original cell IDs')
    if len(set(gene_ids))!=len(gene_ids) or not all(gene_ids):
        raise ValueError('Duplicate or missing original Ensembl gene IDs')
    if not all(x.startswith('ENSG') for x in gene_ids):
        raise ValueError('10x source lacks human Ensembl gene identifiers')
    if len(cell_ids)<n_cells or len(gene_ids)<n_genes:
        raise ValueError('Insufficient source cells/genes')
    if raw.data.size and (np.any(raw.data<0) or not np.isfinite(raw.data).all() or np.any(raw.data!=np.floor(raw.data))):
        raise ValueError('Not valid nonnegative integer UMI counts')
    rng=np.random.default_rng(seed)
    chosen=np.sort(rng.choice(raw.shape[1],size=n_cells,replace=False))
    rng2=np.random.default_rng(seed+1)
    perm=rng2.permutation(n_cells)
    nfit=int(0.6*n_cells);nscore=int(0.2*n_cells)
    fit=perm[:nfit]; score=perm[nfit:nfit+nscore]; test=perm[nfit+nscore:]
    counts=raw[:,chosen].toarray().astype(float).T
    train=counts[fit]
    var=train.var(axis=0)
    expressed=(train>0).sum(axis=0)>=max(1,min(10,int(nfit*0.05)))
    var[~expressed]=-np.inf
    if int(expressed.sum())<n_genes:
        raise ValueError('Insufficient fit-expressed genes')
    # Fit-only feature selection; final test and internal-score cells not consulted.
    selected=np.sort(np.argsort(-var,kind='stable')[:n_genes])
    expression=pd.DataFrame(counts[:,selected].astype(np.int64),columns=[gene_ids[i] for i in selected])
    expression.insert(0,'cell_id',[cell_ids[i] for i in chosen])
    assignment=pd.DataFrame({'cell_id':expression['cell_id'],
                             'partition':['fit' if i in fit else 'score' if i in score else 'test' for i in range(n_cells)]})
    features=pd.DataFrame({'gene_id':[gene_ids[i] for i in selected],
                           'gene_symbol':[gene_names[i] for i in selected],
                           'fit_only_variance':[float(var[i]) for i in selected]})
    if len(set(expression.columns))!=len(expression.columns):
        raise ValueError('Duplicate selected identifiers')
    return expression,assignment,features


def execute(h5: Path, output: Path, *, seed=20261126, cells=800, genes=256):
    import h5py
    with h5py.File(h5,'r') as f:
        m=f['matrix']
        shape=tuple(int(x) for x in m['shape'][:])
        data=m['data'][:];indices=m['indices'][:];indptr=m['indptr'][:]
        def asstr(x):return [v.decode() if isinstance(v,bytes) else str(v) for v in x]
        ids=asstr(m['features']['id'][:]); names=asstr(m['features']['name'][:]); cells_id=asstr(m['barcodes'][:])
        types=asstr(m['features']['feature_type'][:]) if 'feature_type' in m['features'] else ['Gene Expression']*len(ids)
        if len(types)!=len(ids):raise ValueError('Missing gene-expression feature types')
        raw=csc_matrix((data,indices,indptr),shape=shape)
        take=np.where(np.array(types)=='Gene Expression')[0]
        if len(take)==0:raise ValueError('No gene expression features')
        matrix=raw[take,:].tocsc()
        gene_ids=[ids[i] for i in take]; gene_names=[names[i] for i in take]
    X,part,features=qualify_matrix(matrix,gene_ids,gene_names,cells_id, n_cells=cells,n_genes=genes,seed=seed)
    output.mkdir(parents=True,exist_ok=True)
    xfile=output/'K562r_fit_only_selected_expression.csv.gz'
    pfile=output/'cell_partitions.csv'
    ffile=output/'original_feature_mapping.csv'
    X.to_csv(xfile,index=False,compression='gzip')
    part.to_csv(pfile,index=False)
    features.to_csv(ffile,index=False)
    result={'dataset':'10x 10k Human K562-r Cells Singleplex Sample',
            'provider_source_page':SOURCE_PAGE,
            'cell_line':'K562-r','atcc':'CRL-3344','preservation':'fixed RNA / Flex',
            'not_same_line_provenance_as_ENCODE_K562_CRISPRi':True,
            'raw_source_h5_sha256':sha256(h5),
            'source_shape_genes_by_cells':list(shape),
            'n_cells_exported':len(X),
            'n_ensg_features':len(features),
            'partition_n':part['partition'].value_counts().to_dict(),
            'feature_selection':'fit-only variance; no test cell gene selection',
            'selected_feature_id_type':'Ensembl stable IDs from original 10x HDF5',
            'original_expression_matrix_md5_verified':False,
            'expression_sha256':sha256(xfile),
            'partition_sha256':sha256(pfile),
            'feature_mapping_sha256':sha256(ffile),
            'any_ENCODE_perturbation_outcome_loaded':False,
            'regulatory_prior_qualified':False,
            'GCEV_rankings_fitted':False,
            'source_identity_status':'SINGLE_LINE_K562_R_CANDIDATE_NOT_SAME_AS_K562',
            'decision':'EXPLORATORY_TRAINING_CANDIDATE_NEEDS_SOURCE_COHORT_COMPARABILITY_REVIEW'}
    (output/'source_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--h5',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cells',type=int,default=800)
    p.add_argument('--genes',type=int,default=256)
    args=p.parse_args()
    execute(args.h5,args.out,cells=args.cells,genes=args.genes)
if __name__=='__main__':main()
