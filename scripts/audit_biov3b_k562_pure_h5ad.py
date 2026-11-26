"""BIO-V3B authentic K562-only H5AD source inspection; never read ENCODE outcomes.

Audit the entire expression source with sparse preservation and export only when
the original scale is demonstrably usable. No TF selection or final split yet.
"""
from __future__ import annotations
import argparse, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()


def inspect_array(x, cells, genes, source_hash: str) -> dict:
    if x.shape != (len(cells),len(genes)):
        raise ValueError("Gene-by-cell matrix shape inconsistent with the H5AD dimensions")
    if len(set(cells))!=len(cells) or len(set(genes))!=len(genes):
        raise ValueError("Duplicate cell or gene identifiers")
    if not cells or not genes or not all(cells) or not all(genes):
        raise ValueError("Empty cells or genes")
    nonzeros=x.data if sparse.issparse(x) else np.asarray(x).ravel()
    if not np.isfinite(nonzeros).all() or (nonzeros<0).any():
        raise ValueError("Source expression has nonfinite or negative values")
    if nonzeros.size == 0:
        raise ValueError("Empty source expression matrix")
    integerish=float(np.mean(np.isclose(nonzeros,np.rint(nonzeros),atol=1e-6)))
    symbol_like=float(np.mean([bool(re.fullmatch(r"[A-Za-z0-9_.-]+",g)) for g in genes]))
    ensg=int(sum(g.startswith("ENSG") for g in genes))
    return {
        "source_sha256":source_hash,
        "n_cells":len(cells),"n_genes":len(genes),
        "matrix_is_sparse":bool(sparse.issparse(x)),
        "nonzero_entries":int(np.count_nonzero(nonzeros)),
        "min_nonzero":float(np.min(nonzeros)), "max_nonzero":float(np.max(nonzeros)),
        "fraction_nonzero_entries_near_integer":integerish,
        "all_nonnegative_and_finite":True,
        "ensembl_like_gene_ids":ensg,
        "gene_identifier_format_simple_fraction":symbol_like,
        "cell_index_sha256":hashlib.sha256("\n".join(cells).encode()).hexdigest(),
        "gene_index_sha256":hashlib.sha256("\n".join(genes).encode()).hexdigest(),
        "assumed_unperturbed_from_source_metadata_only":True,
        "encode_perturbation_data_accessed":False,
        "gene_annotation_version_independently_verified":False,
        "gcev_training_performed":False,
        "regulatory_prior_qualified":False,
        "cell_type_from_numeric_mixture_inferred":False,
        "matrix_scale_classification":"COUNT_LIKE_NEEDS_SOURCE_CONFIRMATION" if integerish>=0.999 else "NONINTEGER_EXPRESSION_SCALE_NEEDS_SOURCE_CONFIRMATION",
        "decision":"SOURCE_ARRAY_VERIFIED_GENE_AND_SCALE_REVIEW_REQUIRED",
    }


def audit(source: Path, output: Path, export: bool=False) -> dict:
    import anndata
    ad=anndata.read_h5ad(source)
    cells=list(map(str,ad.obs_names))
    genes=list(map(str,ad.var_names))
    findings=inspect_array(ad.X,cells,genes,sha256(source))
    findings["source_name"]=source.name
    findings["obs_metadata_fields"]=list(map(str,ad.obs.columns))
    findings["var_metadata_fields"]=list(map(str,ad.var.columns))
    findings["layers"]=list(map(str,ad.layers.keys()))
    findings["raw_layer_present"]=ad.raw is not None
    output.mkdir(parents=True,exist_ok=True)
    if export:
        x=ad.X if sparse.issparse(ad.X) else sparse.csr_matrix(ad.X)
        sparse.save_npz(output/"k562_expression_original_scale.npz",x.tocsr())
        pd.DataFrame({"cell_id":cells}).to_csv(output/"k562_cell_ids.csv",index=False)
        pd.DataFrame({"gene_id":genes}).to_csv(output/"k562_gene_ids.csv",index=False)
        findings["exported_original_scale_sparse_matrix"]=True
        findings["export_hashes"]={
            p.name:sha256(p) for p in (output/"k562_expression_original_scale.npz",
                                      output/"k562_cell_ids.csv",output/"k562_gene_ids.csv")}
        findings["interpretation"]="Original source expression scale and IDs preserved; do not assume raw counts or normalize before verifying source pipeline."
    else:
        findings["exported_original_scale_sparse_matrix"]=False
    (output/"k562_source_audit.json").write_text(json.dumps(findings,indent=2)+"\n")
    return findings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--export",action="store_true")
    args=p.parse_args()
    results=audit(args.h5ad,args.output,export=args.export)
    print(json.dumps({k:results[k] for k in (
        "source_sha256","n_cells","n_genes","matrix_scale_classification",
        "ensembl_like_gene_ids","decision")},indent=2))


if __name__=="__main__":main()
