"""BIO-V4B: source-anchored TF-aware fit-only selection for a 64-gene K562 cohort.

This is a *development redesign* after the original 512-HVG set contained only
3 JASPAR motif-covered TFs. Same split, never use ENCODE response or GRN truth.
"""
from __future__ import annotations
import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from scripts.biov4_qualify_gene_universe import (
    read_gencode, read_jaspar_meme, sha
)

H5AD_SOURCE_SHA = "14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6"
GENCODE_SHA = "e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9"
JASPAR_SHA = "ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15"


def choose(genes: list[str], detected: np.ndarray, variance: np.ndarray,
           gencode: dict, jaspar: dict, n_tf: int=16, n_target: int=48):
    if len(genes)!=len(detected) or len(genes)!=len(variance):
        raise ValueError("Full gene list and fitting-only statistics must be aligned")
    if len(genes)!=len(set(genes)):
        raise ValueError("Duplicate source gene symbols")
    if not np.isfinite(detected).all() or not np.isfinite(variance).all():
        raise ValueError("Unusable fitting-only variance/detection")
    tf_candidates, target_candidates = [],[]
    for i,symbol in enumerate(genes):
        if (detected[i]<.05 or variance[i]<=0 or
            re.match(r"^(MT-|RPL[0-9]|RPS[0-9])",symbol)):
            continue
        matches=gencode.get(symbol.casefold(),[])
        unique={v["gencode_gene_id"]:v for v in matches}
        if len(unique)!=1:
            continue
        g=next(iter(unique.values()))
        if g["gene_type"]!="protein_coding":
            continue
        motifs=sorted({x["motif_id"] for x in jaspar.get(symbol.casefold(),[])})
        row={"gene_symbol":symbol,"gencode_gene_id":g["gencode_gene_id"],
             "fit_detection":float(detected[i]),
             "fit_log1p_variance":float(variance[i]),
             "motif_ids":"|".join(motifs),
             "motif_count":len(motifs),
             "gene_index":i}
        if motifs:
            tf_candidates.append(row)
        else:
            target_candidates.append(row)
    order=lambda x:(-x["fit_log1p_variance"],x["gene_symbol"])
    tf_candidates.sort(key=order)
    target_candidates.sort(key=order)
    if len(tf_candidates)<n_tf or len(target_candidates)<n_target:
        raise ValueError(f"Insufficient fit-only motif-eligible TFs / nonTF targets: "
                         f"{len(tf_candidates)}/{len(target_candidates)} for {n_tf}/{n_target}")
    chosen_tfs=tf_candidates[:n_tf]
    chosen_targets=target_candidates[:n_target]
    all_selected=chosen_tfs+chosen_targets
    return chosen_tfs,chosen_targets,all_selected,{
        "full_genes":len(genes),
        "fit_detected_and_annotated_motif_TF_candidates":len(tf_candidates),
        "fit_detected_and_annotated_nonTF_target_candidates":len(target_candidates),
        "n_regulators":len(chosen_tfs),"n_additional_targets":len(chosen_targets),
        "n_selected_genes":len(all_selected),
        "directed_candidate_edges_no_self":len(chosen_tfs)*(len(all_selected)-1),
    }


def prepare(source: Path, gtf: Path, meme:Path, dest:Path):
    audited=json.loads((source/"biov3b_pure_k562/k562_source_audit.json").read_text())
    if audited.get("source_sha256")!=H5AD_SOURCE_SHA:
        raise ValueError("Original K562 H5AD source checksum mismatch")
    for name,digest in audited["export_hashes"].items():
        if sha(source/"biov3b_pure_k562"/name)!=digest:
            raise ValueError("Original source matrix / identifier checksum mismatch: "+name)
    original=source/"biov3b_K562_fit_only_512"
    split_manifest=json.loads((original/"split_manifest.json").read_text())
    if (split_manifest["fit_cells"],split_manifest["internal_score_cells"],
        split_manifest["untouched_test_cells"],split_manifest["seed"])!=(394,131,133,20261126):
        raise ValueError("Frozen BIO-V3B partitions unexpectedly modified")
    for name,digest in split_manifest["export_sha256"].items():
        if sha(original/name)!=digest:
            raise ValueError("Original BIO-V3B split checksum mismatch")
    if sha(gtf)!=GENCODE_SHA or sha(meme)!=JASPAR_SHA:
        raise ValueError("Original GENCODE / JASPAR source hashes differ from frozen BIO-V4")
    x=sparse.load_npz(source/"biov3b_pure_k562/k562_expression_original_scale.npz").tocsr()
    cells=pd.read_csv(source/"biov3b_pure_k562/k562_cell_ids.csv",dtype="string")["cell_id"].tolist()
    genes=pd.read_csv(source/"biov3b_pure_k562/k562_gene_ids.csv",dtype="string")["gene_id"].tolist()
    split=pd.read_csv(original/"cell_partitions.csv",dtype="string")
    if x.shape!=(658,60671) or len(cells)!=658 or len(genes)!=60671:
        raise ValueError("Full original K562 matrix shape differs")
    if split["cell_id"].tolist()!=cells or split["partition"].value_counts().to_dict()!={"fit":394,"internal_score":131,"untouched_test":133}:
        raise ValueError("Exact ordered partition identities do not match original K562 cells")
    fit=np.flatnonzero(split["partition"].to_numpy()=="fit")
    xf=x[fit,:].copy()
    xf.eliminate_zeros()
    detected=np.asarray(xf.getnnz(axis=0)).ravel()/len(fit)
    xf.data=np.log1p(xf.data)
    means=np.asarray(xf.mean(axis=0)).ravel()
    sq=xf.copy();sq.data**=2
    variance=np.maximum(0,np.asarray(sq.mean(axis=0)).ravel()-means**2)
    gc=read_gencode(gtf)
    jm=read_jaspar_meme(meme)
    tfs,targets,all_genes,stats=choose(genes,detected,variance,gc,jm)
    selected_idx=[r["gene_index"] for r in all_genes]
    subset=x[:,selected_idx].tocsr()
    dest.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(dest/"k562_tfaware64_original_scale.npz",subset)
    tf_cols=["gene_symbol","gencode_gene_id","fit_detection","fit_log1p_variance","motif_ids","motif_count"]
    all_cols=tf_cols+["role"]
    pd.DataFrame(tfs)[tf_cols].to_csv(dest/"selected_16_TFs.csv",index=False)
    pd.DataFrame(targets)[tf_cols].to_csv(dest/"selected_48_targets.csv",index=False)
    pd.DataFrame([{**{c:r[c] for c in tf_cols},
                   "role":"regulator" if j<16 else "target_only"}
                  for j,r in enumerate(all_genes)])[all_cols].to_csv(dest/"selected_gene_universe.csv",index=False)
    shutil.copyfile(original/"cell_partitions.csv",dest/"cell_partitions.csv")
    provenance={
        "stage":"DEVELOPMENT_TF_AWARE_K562_GENE_COHORT_NOT_CONFIRMATORY",
        "protocol":"GCEV-BIOV4B-K562-TF16-TARGET48-2026-11-26-001",
        "original_K562_H5AD_sha256":H5AD_SOURCE_SHA,
        "original_biov3b_split_sha256":sha(original/"split_manifest.json"),
        "GENCODE_v29_sha256":GENCODE_SHA, "JASPAR_2026_sha256":JASPAR_SHA,
        "n_cells":len(cells),"n_fit_cells":394,"n_internal_score_cells":131,
        "n_untouched_test_cells":133,
        "original_512_HVG_set_retained_not_modified":True,
        "reselected_due_to_observed_original_3_TF_eligibility_failure":True,
        "feature_selection":"fit-only log1p variance; protein-coding unique GENCODE v29; >=5% detection; 16 JASPAR-matched TF + 48 nonmotif target genes",
        "source_expression_scale":"original count-like X not independently verified as raw UMI",
        "no_ENCODE_response_data_opened":True,
        "no_reference_edge_data_opened":True,
        "no_motif_to_target_edges_generated":True,
        "no_GCEV_models_fitted":True,
        "not_prospectively_uninspected_gene_universe":True,
        "candidate_universe_biased_towards_motif_available_TFs":True,
        "counts":stats,
    }
    files=[dest/"k562_tfaware64_original_scale.npz",
           dest/"selected_16_TFs.csv",dest/"selected_48_targets.csv",
           dest/"selected_gene_universe.csv",dest/"cell_partitions.csv"]
    provenance["sha256"]={p.name:sha(p) for p in files}
    (dest/"biov4b_fit_only_manifest.json").write_text(json.dumps(provenance,indent=2)+"\n")
    return provenance


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True,type=Path)
    ap.add_argument("--gencode",required=True,type=Path)
    ap.add_argument("--jaspar",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    a=ap.parse_args()
    report=prepare(a.input,a.gencode,a.jaspar,a.output)
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
