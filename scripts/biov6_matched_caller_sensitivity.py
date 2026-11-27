"""BIO-V6S source-version sensitivity: two matched ENCODE K562 DNase peak callers.

Exploratory process robustness, *not* post-hoc replacement of the original
locked lexicographic cross-pipeline result. No TF-target response or graph gold.
"""
from __future__ import annotations
import argparse
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.audit_encode_crispri_panel import fetch_original
from scripts.biov6_k562_dnase_accessibility import (
    EXPERIMENT, inspect_dnase_original, load_frozen_motif, download_checked,
    parse_and_merge_narrowpeak, covered_bases, compute_joint, digest
)

PIPELINES={
    "v1":"dnase-peak-calling-step-v-1-0",
    "v2":"dnase-call-hotspots-se-step-v-2-1",
}


def select_matched(meta:dict, code:str) -> dict[int,dict]:
    if code not in PIPELINES:raise ValueError("Unregistered pipeline variant")
    ending="/analysis-step-versions/"+PIPELINES[code]+"/"
    chosen={}
    for rep in (1,2):
        candidates=[f for f in meta["eligible_files_by_replicate"][str(rep)]
                    if f.get("pipeline")==ending]
        if not candidates:
            raise ValueError(f"NO-GO: matching source {code} missing replicate {rep}")
        chosen[rep]=sorted(candidates,key=lambda f:f["accession"])[0]
    if chosen[1]["biosample_id"]==chosen[2]["biosample_id"]:
        raise ValueError("Same biosample in both original DNase replicates")
    return chosen


def promoter_stats(coordinates:pd.DataFrame, merged:dict[int,dict]) -> pd.DataFrame:
    rows=[]
    for g in coordinates.itertuples(index=False):
        a,b=int(g.start_0based),int(g.end_0based_exclusive)
        if b-a!=1201:raise ValueError("Source promoter length changed")
        r1=covered_bases(merged[1],str(g.chromosome),a,b)/1201
        r2=covered_bases(merged[2],str(g.chromosome),a,b)/1201
        rows.append({
            "gene_symbol":g.gene_symbol,
            "gencode_gene_id":g.gencode_gene_id,
            "promoter_start0":a,"promoter_end0_exclusive":b,
            "dnase_replicate1_coverage":r1,
            "dnase_replicate2_coverage":r2,
            "dnase_accessibility_mean":(r1+r2)/2,
            "dnase_accessibility_min":min(r1,r2),
            "dnase_both_replicates_peak_present":bool(r1>0 and r2>0),
            "dnase_replicate_abs_disagreement":abs(r1-r2)
        })
    return pd.DataFrame(rows)


def run(biov5:Path,original_primary:Path,out:Path) -> dict:
    edges,coordinates,motif_manifest=load_frozen_motif(biov5)
    obj,source_sha=fetch_original(EXPERIMENT)
    meta=inspect_dnase_original(obj,source_sha)
    out.mkdir(parents=True,exist_ok=True)
    report={
        "protocol":"BIOV6_SECONDARY_MATCHED_PIPELINE_SENSITIVITY_2026_11_27",
        "stage":"DEVELOPMENT_SECONDARY_CALLER_SENSITIVITY_NOT_INDEPENDENT_VALIDATION",
        "original_source_experiment":"ENCSR000EOT",
        "original_experiment_json_sha256":source_sha,
        "parent_biov5_prior_sha256":motif_manifest["file_sha256"]["k562_tf_to_promoter_motif_prior_1008.csv"],
        "original_locked_first_analysis_unmodified":True,
        "ENCODE_CRISPRi_response_outcomes_accessed":False,
        "GCEV_method_fitted":False,
        "two_matched_peak_callers":{},
        "file_sha256":{}
    }
    originals={}
    if original_primary.exists():
        p=pd.read_csv(original_primary)
        if len(p)!=64:raise ValueError("Original primary 64-promoter output missing")
        originals=p.set_index("gene_symbol")
    summaries={}
    with tempfile.TemporaryDirectory(prefix="biov6_matched_caller_") as td:
        for code in ("v1","v2"):
            selected=select_matched(meta,code)
            sources=[]
            merged={}
            for rep in (1,2):
                f=selected[rep]
                target=Path(td)/(f"{code}_{f['accession']}.bed.gz")
                source=download_checked(f,target)
                merged[rep],n=parse_and_merge_narrowpeak(target)
                sources.append({**source,**n})
            coverage=promoter_stats(coordinates,merged)
            joint=compute_joint(edges,coverage)
            coverage_path=out/f"k562_dnase_{code}_matched_promoter_coverage_64.csv"
            joint_path=out/f"k562_dnase_{code}_matched_joint_context_1008.csv"
            coverage.to_csv(coverage_path,index=False)
            joint.to_csv(joint_path,index=False)
            report["file_sha256"][coverage_path.name]=digest(coverage_path)
            report["file_sha256"][joint_path.name]=digest(joint_path)
            metrics={
                "selected_original_files":sources,
                "promoters_peak_rep1":int(coverage.dnase_replicate1_coverage.gt(0).sum()),
                "promoters_peak_rep2":int(coverage.dnase_replicate2_coverage.gt(0).sum()),
                "promoters_both_peak":int(coverage.dnase_both_replicates_peak_present.sum()),
                "promoters_no_peak_both":int(coverage.dnase_accessibility_mean.eq(0).sum()),
                "mean_promoter_coverage":float(coverage.dnase_accessibility_mean.mean()),
                "mean_absolute_replicate_disagreement":float(coverage.dnase_replicate_abs_disagreement.mean()),
                "replicate_coverage_spearman":float(coverage.dnase_replicate1_coverage.corr(coverage.dnase_replicate2_coverage,method="spearman")),
                "positive_joint_motif_dnase_product":int(joint.joint_motif_dnase_product.gt(0).sum()),
                "mean_joint_motif_dnase_product":float(joint.joint_motif_dnase_product.mean())
            }
            if originals:
                c=coverage.set_index("gene_symbol")
                metrics["mean_abs_difference_vs_original_heterogeneous_coverage"]=float(
                    (c.loc[originals.index,"dnase_accessibility_mean"]-originals["dnase_accessibility_mean"]).abs().mean())
            report["two_matched_peak_callers"][code]=metrics
            summaries[code]=coverage.set_index("gene_symbol")
            print("MATCHED_PIPELINE",code,
                  "BOTH",metrics["promoters_both_peak"],
                  "MEAN",metrics["mean_promoter_coverage"],flush=True)
    a=summaries["v1"];b=summaries["v2"]
    report["v1_v2_mean_abs_promoter_coverage_difference"]=float(
        (a["dnase_accessibility_mean"]-b["dnase_accessibility_mean"]).abs().mean())
    report["v1_v2_promoter_coverage_spearman"]=float(
        a["dnase_accessibility_mean"].corr(b["dnase_accessibility_mean"],method="spearman"))
    report["two_distinct_independent_assay_types"]=False
    report["accessibility_is_tf_specific"]=False
    report["joint_motif_dnase_product_is_third_independent_source"]=False
    report["not_binding_truth"]=True
    (out/"biov6_secondary_pipeline_sensitivity_manifest.json").write_text(
        json.dumps(report,indent=2)+"\n")
    return report


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--biov5",required=True,type=Path)
    ap.add_argument("--original-primary",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    x=ap.parse_args()
    report=run(x.biov5,x.original_primary,x.output)
    print(json.dumps({k:report[k] for k in (
        "stage","two_matched_peak_callers",
        "v1_v2_mean_abs_promoter_coverage_difference",
        "v1_v2_promoter_coverage_spearman")},indent=2))


if __name__=="__main__":main()
