"""BIO-V4: qualify real K562 fit-only gene symbols and JASPAR TF motifs.

No reference/CRISPRi response access, promoter scan, TF-target prediction or
model training. Only GENCODE gene annotation, MEME motif identities and
expression detection on the already frozen outer-fitting cells.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
GENE_ID = re.compile(r"^ENSG[0-9]+(?:\.[0-9]+)?$")
MEME_ID = re.compile(r"^MA[0-9]+\.[0-9]+$")
SINGLE_TF = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda:fh.read(1048576),b""):
            h.update(chunk)
    return h.hexdigest()


def attr(text: str) -> dict:
    return dict(re.findall(r'(\w+)\s+"([^"]+)"', text))


def read_gencode(path: Path) -> dict[str, list[dict]]:
    """Only genuine primary-chromosome gene rows, keeping ambiguity explicit."""
    records: dict[str,list[dict]]=defaultdict(list)
    with gzip.open(path,"rt",encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)!=9 or c[2]!="gene" or c[0] not in PRIMARY:
                continue
            d=attr(c[8])
            name=d.get("gene_name","")
            gid=d.get("gene_id","")
            kind=d.get("gene_type",d.get("gene_biotype",""))
            if not name or not GENE_ID.fullmatch(gid) or not kind:
                continue
            records[name.casefold()].append({
                "gene_symbol":name,"gencode_gene_id":gid,
                "gene_type":kind, "chromosome":c[0],
                "start_1based":int(c[3]),"end_1based":int(c[4]),
                "strand":c[6],
            })
    if len(records)<1000:
        raise ValueError("GENCODE gene annotation unexpectedly sparse or incompatible")
    return records


def read_jaspar_meme(path: Path) -> dict[str,list[dict]]:
    """Exact single-factor names only; multi-factor motifs stay excluded."""
    motifs: dict[str,list[dict]]=defaultdict(list)
    with path.open("r",encoding="utf-8") as f:
        header=f.readline()
        if not header.startswith("MEME version"):
            raise ValueError("Not a source MEME motif file")
        for line in f:
            if not line.startswith("MOTIF "):continue
            tokens=line.rstrip("\n").split(maxsplit=2)
            if len(tokens)!=3 or not MEME_ID.fullmatch(tokens[1]):
                continue
            name=tokens[2].strip()
            if not SINGLE_TF.fullmatch(name) or "::" in name:
                continue
            motifs[name.casefold()].append({
                "motif_id":tokens[1],"jaspar_factor_name":name})
    if len(motifs)<100:
        raise ValueError("Unexpectedly few JASPAR 2026 single-factor motifs")
    return motifs


def inspect_universe(selected: list[str], fit_detected: np.ndarray,
                     gencode: dict, motifs: dict, min_tfs: int=4):
    if len(selected)!=len(fit_detected) or len(set(selected))!=len(selected):
        raise ValueError("Feature order/uniqueness mismatch")
    if len(selected)<5:
        raise ValueError("Too few expression-selected genes")
    annotation=[]
    eligible=[]
    matched_motifs=[]
    for symbol, fraction in zip(selected, fit_detected):
        record=gencode.get(symbol.casefold(),[])
        distinct={r["gencode_gene_id"]:r for r in record}
        status="not_found" if not distinct else ("ambiguous_gene_id" if len(distinct)>1 else "unique_gene_id")
        g=next(iter(distinct.values())) if status=="unique_gene_id" else None
        tf_profiles=motifs.get(symbol.casefold(),[])
        motif_ids=sorted({m["motif_id"] for m in tf_profiles})
        annotation.append({
            "selected_gene_symbol":symbol,
            "gene_annotation_status":status,
            "gencode_gene_id":g["gencode_gene_id"] if g else None,
            "gencode_official_symbol":g["gene_symbol"] if g else None,
            "gene_type":g["gene_type"] if g else None,
            "chromosome":g["chromosome"] if g else None,
            "fit_detected_fraction":float(fraction),
            "single_factor_motif_count":len(motif_ids),
        })
        if g and motif_ids:
            matched_motifs += [
                {"tf_symbol":symbol,"gencode_gene_id":g["gencode_gene_id"],
                 "motif_id":m["motif_id"],"jaspar_factor_name":m["jaspar_factor_name"]}
                for m in tf_profiles
            ]
        if (g and g["gene_type"]=="protein_coding"
            and float(fraction)>=0.05 and motif_ids):
            eligible.append({
                "tf_symbol":symbol, "gencode_gene_id":g["gencode_gene_id"],
                "fit_detected_fraction":float(fraction),
                "n_single_factor_motifs":len(motif_ids),
                "jaspar_motif_ids":"|".join(motif_ids),
            })
    ann=pd.DataFrame(annotation)
    tf=pd.DataFrame(eligible,columns=[
        "tf_symbol","gencode_gene_id","fit_detected_fraction",
        "n_single_factor_motifs","jaspar_motif_ids"])
    motif=pd.DataFrame(matched_motifs,columns=[
        "tf_symbol","gencode_gene_id","motif_id","jaspar_factor_name"])
    counts={
        "selected_genes":len(selected),
        "unique_gencode_gene_ids":int((ann.gene_annotation_status=="unique_gene_id").sum()),
        "ambiguous_gencode_gene_ids":int((ann.gene_annotation_status=="ambiguous_gene_id").sum()),
        "unmatched_genes":int((ann.gene_annotation_status=="not_found").sum()),
        "unique_protein_coding_genes":int((ann.gene_type=="protein_coding").sum()),
        "mitochondrial_prefix_genes":sum(str(s).startswith("MT-") for s in selected),
        "ribosomal_prefix_genes":sum(bool(re.match(r"^RP[LS][0-9]",str(s))) for s in selected),
        "genes_with_single_factor_motif_profiles":int((ann.single_factor_motif_count>0).sum()),
        "eligible_regulators":len(eligible),
        "eligible_matrix_ids":int(len(motif)),
    }
    decision="MOTIF_TF_ELIGIBILITY_OBSERVED_NEEDS_PROMOTER_SCANNING" if len(eligible)>=min_tfs else "NO_GO_INSUFFICIENT_TF_COVERAGE_REVISE_FIT_ONLY_SELECTION"
    return ann, motif, tf, counts, decision


def run(input_root:Path, gtf:Path, meme:Path, out:Path, min_tfs:int=4):
    source=input_root/"biov3b_K562_fit_only_512"
    manifest=json.loads((source/"split_manifest.json").read_text())
    if manifest.get("original_H5AD_sha256")!="14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6":
        raise ValueError("Unexpected original K562-only source")
    if manifest.get("fit_cells")!=394 or manifest.get("internal_score_cells")!=131 or manifest.get("untouched_test_cells")!=133:
        raise ValueError("Frozen K562 fit/score/test split not preserved")
    if manifest.get("n_selected_genes")!=512 or manifest.get("seed")!=20261126:
        raise ValueError("Unexpected original fit-only gene universe")
    for filename,expected in manifest["export_sha256"].items():
        if sha(source/filename)!=expected:
            raise ValueError("BIOV3B output checksum mismatch: "+filename)
    selected=pd.read_csv(source/"selected_gene_symbols.csv",dtype="string")["gene_symbol"].tolist()
    split=pd.read_csv(source/"cell_partitions.csv",dtype="string")
    x=sparse.load_npz(source/"k562_selected_original_scale.npz").tocsr()
    if x.shape!=(658,512) or len(selected)!=512 or len(split)!=658 or split["cell_id"].duplicated().any():
        raise ValueError("Frozen expression matrix shape/cell identifiers invalid")
    fit=np.flatnonzero(split["partition"].to_numpy()=="fit")
    if len(fit)!=394:
        raise ValueError("Unverified fitting-only split")
    xf=x[fit,:].copy()
    xf.eliminate_zeros()
    detected=np.asarray(xf.getnnz(axis=0)).ravel()/len(fit)
    gc=read_gencode(gtf)
    jm=read_jaspar_meme(meme)
    ann,motifs,tfs,counts,decision=inspect_universe(selected,detected,gc,jm,min_tfs)
    out.mkdir(parents=True,exist_ok=True)
    ann.to_csv(out/"biov4_qualifying_gene_annotations.csv",index=False)
    motifs.to_csv(out/"biov4_single_tf_motifs.csv",index=False)
    tfs.to_csv(out/"biov4_eligible_regulators.csv",index=False)
    report={
        "protocol":"GCEV-BIOV4-GENE-TF-JASPAR-2026-11-26-001",
        "original_gse306708_h5ad_sha256":manifest["original_H5AD_sha256"],
        "original_split_manifest_sha256":sha(source/"split_manifest.json"),
        "gencode_release":"29", "gencode_source_sha256":sha(gtf),
        "jaspar_release":"2026_CORE_vertebrates_nonredundant",
        "jaspar_meme_sha256":sha(meme),
        "total_gencode_unambiguous_names":sum(len({v["gencode_gene_id"] for v in vv})==1 for vv in gc.values()),
        "jaspar_single_factor_name_count":len(jm),
        "fit_cells_only_used_for_detection":394,
        "unmodified_original_512_genes":True,
        "cell_partitions_unchanged":True,
        "no_encode_response_access":True,
        "tf_binding_prior_edges_computed":False,
        "no_pwm_promoter_scanning_performed":True,
        "no_gcev_models_fitted":True,
        "no_biological_accuracy_evaluated":True,
        "counts":counts,
        "decision":decision,
        "next_step":"Pin GRCh38 promoter coordinates/sequences and source-specific motif scans; if regulator coverage inadequate, freeze a new expression-fit-only TF-anchored gene selection before modelling.",
        "results_sha256": {
            p.name:sha(p) for p in [out/"biov4_qualifying_gene_annotations.csv",
                                   out/"biov4_single_tf_motifs.csv",
                                   out/"biov4_eligible_regulators.csv"]
        }
    }
    (out/"biov4_source_eligibility.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--gencode",type=Path,required=True)
    p.add_argument("--jaspar",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=run(a.input,a.gencode,a.jaspar,a.output)
    print(json.dumps({"counts":r["counts"],"decision":r["decision"]},indent=2))


if __name__=="__main__":main()
