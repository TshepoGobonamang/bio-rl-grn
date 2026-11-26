"""BIO-V2: independent integrity checks of real DESeq2 result artifacts.

Recomputes published counts of tested/significant genes from saved gene-wise
tables. It does not rerun the statistical model, validate direct TF binding,
or compare GCEV with biological outcomes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
import numpy as np

TF_SET = ("TFDP1", "KLF2", "NR2F2", "RELA", "SP1", "NR2C2")
REQUIRED = ("gene_id", "baseMean", "log2FoldChange", "lfcSE",
            "stat", "pvalue", "padj", "bh_fdr_005")


def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda:fh.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def _as_bool(col: pd.Series) -> pd.Series:
    x=col.astype("string").str.lower().str.strip()
    if not x.isin(["true","false","1","0"]).all():
        raise ValueError("Invalid BH-FDR boolean values")
    return x.isin(["true","1"])


def verify_one(tf: str, frame: pd.DataFrame, summary: pd.Series, expected_input: int) -> dict:
    if not set(REQUIRED).issubset(frame):
        raise ValueError(f"{tf}: missing gene-wise DESeq2 output columns")
    if len(frame)==0 or frame["gene_id"].isna().any() or frame["gene_id"].duplicated().any():
        raise ValueError(f"{tf}: empty/duplicate gene identifiers")
    if int(summary["genes_input"]) != int(expected_input):
        raise ValueError(f"{tf}: input gene count differs from original QC")
    if int(summary["genes_retained"]) != len(frame) or len(frame)>expected_input:
        raise ValueError(f"{tf}: gene counts inconsistent with DESeq2 summary")
    for col in ("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"):
        # Missing DESeq2 p/adjusted p is permitted because of filtering/outlier flags.
        x=pd.to_numeric(frame[col],errors="coerce")
        if frame[col].notna().sum()!=x.notna().sum():
            raise ValueError(f"{tf}: {col} has nonnumeric values")
        if (x.dropna().abs()==float("inf")).any():
            raise ValueError(f"{tf}: {col} has infinite values")
        if col=="baseMean" and ((x.dropna()<0).any() or x.isna().any()):
            raise ValueError(f"{tf}: invalid base mean")
        if col=="lfcSE" and (x.dropna()<0).any():
            raise ValueError(f"{tf}: negative standard error")
        if col in ("padj","pvalue") and ((x.dropna()<0).any() or (x.dropna()>1).any()):
            raise ValueError(f"{tf}: p values outside [0,1]")
    padj=pd.to_numeric(frame["padj"],errors="coerce")
    pval=pd.to_numeric(frame["pvalue"],errors="coerce")
    lfc=pd.to_numeric(frame["log2FoldChange"],errors="coerce")
    sig=padj.notna() & (padj<0.05)
    calls=_as_bool(frame["bh_fdr_005"])
    if not calls.equals(sig.astype(bool)):
        raise ValueError(f"{tf}: gene-wise FDR calls inconsistent with saved p-adj values")
    found={
        "tested_pvalues":int(pval.notna().sum()),
        "tested_padj":int(padj.notna().sum()),
        "fdr_005":int(sig.sum()),
        "fdr_005_up":int((sig & (lfc>0)).sum()),
        "fdr_005_down":int((sig & (lfc<0)).sum()),
    }
    for key,value in found.items():
        if int(summary[key])!=value:
            raise ValueError(f"{tf}: summary {key} mismatch ({summary[key]} != {value})")
    if int(summary["n_controls"])<2 or int(summary["n_perturbed"])<2:
        raise ValueError(f"{tf}: insufficient original biosample-level replication")
    if str(summary["is_direct_TF_binding_truth"]).lower() not in ("false","0"):
        raise ValueError(f"{tf}: false direct binding claim")
    if str(summary["is_independent_validation_of_GCEV"]).lower() not in ("false","0"):
        raise ValueError(f"{tf}: false GCEV external validation claim")
    return {"tf":tf, "n_input":int(expected_input), "n_retained":len(frame),
            **found,"direct_regulatory_truth_verified":False,
            "gcev_biological_accuracy_test_performed":False}


def verify(qc:Path, deseq2:Path, expected_tfs=TF_SET) -> dict:
    source=json.loads((qc/"qc_manifest.json").read_text())
    if source.get("stage")!="BIOV2_EXPLORATORY_REPLICATE_QC_ONLY":
        raise ValueError("Not a source-authenticated BIO-V2 quality manifest")
    n_ensg=int(source["n_common_ensembl_gene_ids"])
    if n_ensg<1000 or int(source["n_individual_biosamples"])!=17:
        raise ValueError("Unexpected original ENCODE source gene/biosample coverage")
    if int(source["n_control_biosamples"])!=2:
        raise ValueError("Unexpected original shared control design")
    summary_file=deseq2/"deseq2_summary_all_six.csv"
    summary=pd.read_csv(summary_file)
    if summary["tf"].duplicated().any() or set(summary["tf"])!=set(expected_tfs):
        raise ValueError("DESeq2 summary is missing frozen TFs or includes replacements")
    rows,files=[],{}
    for tf in expected_tfs:
        individual=deseq2/tf/"summary.csv"
        row=pd.read_csv(individual)
        if len(row)!=1 or str(row["tf"].iloc[0])!=tf:
            raise ValueError(f"{tf}: invalid per-TF DESeq2 summary")
        pooled=summary.set_index("tf").loc[tf]
        for col in row.columns:
            if col!="tf" and str(row[col].iloc[0])!=str(pooled[col]):
                raise ValueError(f"{tf}: aggregate summary differs from per-TF result {col}")
        gene_file=deseq2/tf/"deseq2_genewise.csv.gz"
        table=pd.read_csv(gene_file,compression="gzip",dtype={"gene_id":"string"})
        rows.append(verify_one(tf,table,pooled,n_ensg))
        files[tf]={"gene_wise_sha256":digest(gene_file),
                   "summary_sha256":digest(individual)}
    version=(deseq2/"deseq2_environment.txt").read_text()
    if "DESeq2_version=" not in version or "tximport_version=" not in version:
        raise ValueError("DESeq2/tximport version metadata missing")
    return {
        "stage":"BIOV2_REAL_DESEQ2_ARTIFACT_CONSISTENCY_VERIFIED",
        "original_qc_sha256":digest(qc/"qc_manifest.json"),
        "deseq2_summary_sha256":digest(summary_file),
        "method":"DESeq2 original RSEM estimated counts with tximport effective lengths",
        "software_environment":version,
        "six_tf_response_stats":rows,
        "source_files":files,
        "multiple_testing_note":"BH-FDR 0.05 within each TF contrast, not six-TF global FDR",
        "control_reuse_note":"Two biological control biosamples shared across all six contrasts; TF results not independent.",
        "knockdown_efficacy_review_pending":True,
        "full_encode_qc_review_pending":True,
        "source_disjointness_verified":False,
        "true_direct_regulatory_edges":False,
        "gcev_biological_superiority_evaluated":False,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--qc",required=True,type=Path)
    ap.add_argument("--deseq2",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    result=verify(args.qc,args.deseq2)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"stage":result["stage"],
                      "tf_counts":result["six_tf_response_stats"]},indent=2))


if __name__=="__main__":
    main()
