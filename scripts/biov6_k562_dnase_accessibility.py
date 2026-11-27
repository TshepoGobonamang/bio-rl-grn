"""BIO-V6: authentic K562 DNase-seq replicate accessibility at frozen GENCODE promoters.

Independent *assay-type* evidence, but target-level (not TF binding). No
CRISPRi outcomes, no gold edges and no GCEV model fitting.
"""
from __future__ import annotations
import argparse
from bisect import bisect_right
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re
import tempfile
import urllib.request

import numpy as np
import pandas as pd

from scripts.audit_encode_crispri_panel import fetch_original, _biosample, record_id

EXPERIMENT = "ENCSR000EOT"
HOST = "https://www.encodeproject.org"
EXPECTED_MOTIF_SHA = "ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0"
EXPECTED_STAGE = "BIOV5_IN_SILICO_SEQUENCE_BINDING_POTENTIAL_NO_BIOLOGICAL_GOLD"
PROTOCOL = "GCEV-BIOV6-K562-DNASE-2026-11-27-001"
CHROMS = {f"chr{i}" for i in range(1,23)} | {"chrX","chrY"}


def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20), b""):
            h.update(block)
    return h.hexdigest()


def inspect_dnase_original(obj: dict, original_sha: str) -> dict:
    """Choose exactly one peak file per genuine independent biosample."""
    if (obj.get("accession")!=EXPERIMENT or obj.get("status")!="released"
        or "dnase" not in str(obj.get("assay_title",obj.get("assay_term_name",""))).lower()
        or "k562" not in str(obj.get("biosample_summary","")).lower()
        or str(obj.get("replication_type","")).lower()!="isogenic"):
        raise ValueError("Original ENCODE DNase assay identity, status or sample invalid")
    reps=defaultdict(set)
    for item in obj.get("replicates",[]):
        if not isinstance(item,dict):
            raise ValueError("Unresolved original replicate metadata")
        n=item.get("biological_replicate_number")
        b=_biosample(item)
        if n is None or not b:
            raise ValueError("Cannot establish independent biological biosample")
        reps[int(n)].add(b)
    if not all(len(reps.get(k,set()))==1 for k in (1,2)):
        raise ValueError("Required K562 biological replicates 1 and 2 not established")
    if next(iter(reps[1]))==next(iter(reps[2])):
        raise ValueError("Biological replicate biosamples are duplicated")
    found={1:[],2:[]}
    for item in obj.get("files",[]):
        if not isinstance(item,dict): continue
        if item.get("status")!="released" or item.get("assembly")!="GRCh38":
            continue
        if str(item.get("file_format","")).lower()!="bed" or (
            "narrowpeak" not in str(item.get("file_format_type","")).lower()):
            continue
        if str(item.get("output_type","")).lower()!="peaks":
            continue
        numbers=item.get("biological_replicates")
        if not isinstance(numbers,list) or len(numbers)!=1:continue
        try: n=int(numbers[0])
        except (ValueError,TypeError):continue
        if n not in found:continue
        a=record_id(item)
        href=item.get("href") or ""
        md5=str(item.get("md5sum") or "").lower()
        size=item.get("file_size")
        if not (re.fullmatch(r"ENCFF[A-Z0-9]{6}",a) and
            href.startswith(f"/files/{a}/@@download/") and href.endswith(".bed.gz")
            and re.fullmatch(r"[0-9a-f]{32}",md5)
            and isinstance(size,int) and 0<size<=50_000_000):
            continue
        found[n].append({"accession":a,"bio_replicate":n,
            "biosample_id":next(iter(reps[n])),"href":href,
            "md5sum":md5,"file_size":size,
            "source_file_format_type":item.get("file_format_type"),
            "source_output_type":item.get("output_type"),
            "source_assembly":item.get("assembly"),
            "pipeline":record_id(item.get("analysis_step_version") or item.get("analysis") or "")})
    if not all(found.values()):
        raise ValueError("NO-GO: missing released GRCh38 narrowPeak per biological replicate; do not substitute")
    selected={n:sorted(files,key=lambda r:r["accession"])[0] for n,files in found.items()}
    return {"protocol":PROTOCOL,"original_experiment_accession":EXPERIMENT,
        "experiment_url":f"{HOST}/experiments/{EXPERIMENT}/",
        "original_experiment_json_sha256":original_sha,
        "biosample_ids_by_replicate":{str(n):list(reps[n])[0] for n in (1,2)},
        "original_replicate_and_library_sha256":obj.get("_gcev_original_replicate_sha256",{}),
        "eligible_files_by_replicate":{str(n):found[n] for n in (1,2)},
        "selected_files_by_replicate":{str(n):selected[n] for n in (1,2)},
        "unreviewed_source_audits":obj.get("audit",obj.get("audits",None)),
        "biological_replicates_are_independent_cultures":False,
        "replication_type":"isogenic",
        "source_file_discovery_pass":True}


def download_checked(info: dict, dest: Path) -> dict:
    url=HOST+info["href"]
    req=urllib.request.Request(url,headers={"User-Agent":"GCEV-BIOV6-ENCODE-DNase/0.1",
        "Accept":"application/octet-stream"})
    h=hashlib.md5()
    n=0
    with urllib.request.urlopen(req,timeout=120) as response, dest.open("wb") as out:
        for chunk in iter(lambda:response.read(1<<20),b""):
            n+=len(chunk)
            if n>50_000_000: raise ValueError("Original DNase BED too large")
            h.update(chunk);out.write(chunk)
    if n!=info["file_size"] or h.hexdigest()!=info["md5sum"]:
        raise ValueError("Original ENCODE DNase file MD5/size mismatch")
    with gzip.open(dest,"rb") as f:
        if not f.read(1):raise ValueError("Original peak BED is empty")
    return {"download_url":url,"sha256":digest(dest),"verified_size_bytes":n,
            "md5":h.hexdigest(),**info}


def parse_and_merge_narrowpeak(path: Path) -> tuple[dict[str,list[tuple[int,int]]],dict]:
    source=defaultdict(list)
    n=0
    with gzip.open(path,"rt") as f:
        for line in f:
            if not line.strip() or line.startswith(("#","track","browser")):continue
            c=line.rstrip("\n").split("\t")
            if len(c)<10:raise ValueError("Original narrowPeak row lacks required 10 BED fields")
            chrom=c[0]
            if chrom not in CHROMS:
                # Explicitly discard reference contigs not in the frozen 64 primary promoters.
                continue
            try:
                start,end=int(c[1]),int(c[2])
                # Require all numerical narrowPeak fields, without interpreting their scores.
                for idx in (4,6,7,8,9):
                    float(c[idx])
            except (ValueError,IndexError) as exc:
                raise ValueError("Invalid original ENCODE narrowPeak row") from exc
            if start<0 or end<=start:
                raise ValueError("Invalid original genomic BED half-open coordinates")
            source[chrom].append((start,end));n+=1
    if n<1000:raise ValueError("Unexpectedly sparse primary-chromosome DNase peak file")
    result={}
    for chrom,intervals in source.items():
        intervals.sort()
        merged=[]
        for a,b in intervals:
            if merged and a<=merged[-1][1]:
                merged[-1]=(merged[-1][0],max(merged[-1][1],b))
            else:
                merged.append((a,b))
        result[chrom]=merged
    return result,{"original_primary_peak_intervals":n,
        "merged_nonoverlapping_intervals":sum(len(v) for v in result.values()),
        "original_peaks_by_primary_chrom":{ch:len(v) for ch,v in sorted(source.items())}}


def covered_bases(merged: dict, chrom: str, start: int, end: int) -> int:
    if chrom not in merged:return 0
    intervals=merged[chrom]
    # Not assumed nonoverlapping in the original BED: merged input is required.
    hi=[b for a,b in intervals]
    idx=bisect_right(hi,start)
    covered=0
    while idx<len(intervals) and intervals[idx][0]<end:
        a,b=intervals[idx]
        covered+=max(0,min(b,end)-max(a,start))
        idx+=1
    if covered>end-start:raise AssertionError("Union coverage exceeds promoter length")
    return covered


def load_frozen_motif(prior_root:Path) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    m=json.loads((prior_root/"biov5_motif_prior_manifest.json").read_text())
    if (m.get("stage")!=EXPECTED_STAGE or m.get("tested_directed_no_self_edges")!=1008 or
        m.get("ENCODE_perturbation_response_labels_accessed") is not False or
        m.get("gcev_model_fitted") is not False):
        raise ValueError("Not the frozen, outcome-blind BIO-V5 source")
    e=prior_root/"k562_tf_to_promoter_motif_prior_1008.csv"
    p=prior_root/"gencode29_promoter_coordinates_and_sources.csv"
    if digest(e)!=EXPECTED_MOTIF_SHA:
        raise ValueError("BIO-V5 regulatory motif scores changed from the locked cohort")
    for path in (e,p):
        if m["file_sha256"].get(path.name)!=digest(path):
            raise ValueError("BIO-V5 original output file checksum mismatch: "+path.name)
    edge=pd.read_csv(e)
    coords=pd.read_csv(p)
    if (len(edge)!=1008 or edge.duplicated(["tf","target"]).any() or
        len(coords)!=64 or coords["gene_symbol"].duplicated().any() or
        not (coords["length"]==1201).all()):
        raise ValueError("Prior edge/promoter dimensions incorrect")
    return edge,coords,m


def compute_joint(edge: pd.DataFrame, promoters: pd.DataFrame) -> pd.DataFrame:
    result=edge[["tf","target","tf_gencode_gene_id","target_gencode_gene_id",
         "prior_score","empirical_p","bh_q"]].rename(
             columns={"prior_score":"motif_prior_score","empirical_p":"motif_empirical_p",
                      "bh_q":"motif_bh_q"}).copy()
    p=promoters.set_index("gene_symbol")
    if result.target.isin(p.index).sum()!=len(result):
        raise ValueError("Missing promoter in fixed sequence universe")
    for column in ("dnase_replicate1_coverage","dnase_replicate2_coverage",
                   "dnase_accessibility_mean","dnase_accessibility_min",
                   "dnase_both_replicates_peak_present","dnase_replicate_abs_disagreement"):
        result[column]=result["target"].map(p[column])
    result["joint_motif_dnase_product"]=(
        result["motif_prior_score"]*result["dnase_accessibility_mean"])
    result["accessibility_source_observed"]=True
    result["direct_tf_binding_verified"]=False
    result["independent_regulatory_edge_truth"]=False
    if (result.isna().any().any() or not result.dnase_accessibility_mean.between(0,1).all()
        or not result.joint_motif_dnase_product.between(0,1).all()
        or result.duplicated(["tf","target"]).any() or (result.tf==result.target).any()
        or len(result)!=1008):
        raise ValueError("Incomplete or invalid 1008-edge contextual prior")
    if result.groupby("target").dnase_accessibility_mean.nunique().max()!=1:
        raise AssertionError("Accessibility must remain TF-independent within target")
    return result


def run(prior_root:Path, out:Path) -> dict:
    edges,coords,motif_meta=load_frozen_motif(prior_root)
    obj,original_sha=fetch_original(EXPERIMENT)
    metadata=inspect_dnase_original(obj,original_sha)
    out.mkdir(parents=True,exist_ok=True)
    source_files=[]
    with tempfile.TemporaryDirectory(prefix="biov6_dnase_") as temp:
        merged={}
        for n in (1,2):
            info=metadata["selected_files_by_replicate"][str(n)]
            path=Path(temp)/(info["accession"]+".bed.gz")
            file_info=download_checked(info,path)
            merged[n],peak_counts=parse_and_merge_narrowpeak(path)
            source_files.append({**file_info,**peak_counts})
            print("DNASE_SOURCE_VERIFIED",n,info["accession"],file_info["sha256"],flush=True)
    rows=[]
    for g in coords.itertuples(index=False):
        chrom=str(g.chromosome)
        a,b=int(g.start_0based),int(g.end_0based_exclusive)
        if b-a!=1201:raise ValueError("Original BIO-V5 promoter window changed")
        c1=covered_bases(merged[1],chrom,a,b)/1201
        c2=covered_bases(merged[2],chrom,a,b)/1201
        rows.append({"gene_symbol":g.gene_symbol,
            "gencode_gene_id":g.gencode_gene_id,"chromosome":chrom,
            "promoter_start0":a,"promoter_end0_exclusive":b,
            "promoter_gene_strand":g.gencode_strand,
            "dnase_replicate1_coverage":c1,
            "dnase_replicate2_coverage":c2,
            "dnase_accessibility_mean":(c1+c2)/2,
            "dnase_accessibility_min":min(c1,c2),
            "dnase_both_replicates_peak_present":bool(c1>0 and c2>0),
            "dnase_replicate_abs_disagreement":abs(c1-c2)})
    promoters=pd.DataFrame(rows)
    joint=compute_joint(edges,promoters)
    promoters.to_csv(out/"k562_dnase_promoter_accessibility_64.csv",index=False)
    joint.to_csv(out/"k562_tf_target_accessibility_and_motif_context_1008.csv",index=False)
    report={
      **metadata,"stage":"BIOV6_REPLICATE_DNASE_ACCESSIBILITY_CONTEXT_NO_TF_BINDING_GOLD",
      "source_bed_file_summaries":source_files,
      "original_biov5_motif_prior_sha256":EXPECTED_MOTIF_SHA,
      "original_biov5_source_manifest_sha256":digest(prior_root/"biov5_motif_prior_manifest.json"),
      "n_promoters":len(promoters),"n_directed_edges":len(joint),
      "promoters_with_peak_replicate1":int(promoters.dnase_replicate1_coverage.gt(0).sum()),
      "promoters_with_peak_replicate2":int(promoters.dnase_replicate2_coverage.gt(0).sum()),
      "promoters_with_peak_both_replicates":int(promoters.dnase_both_replicates_peak_present.sum()),
      "mean_coverage_over_promoters":float(promoters.dnase_accessibility_mean.mean()),
      "mean_absolute_replicate_disagreement":float(promoters.dnase_replicate_abs_disagreement.mean()),
      "nonzero_joint_prior_scores":int(joint.joint_motif_dnase_product.gt(0).sum()),
      "joint_prior_is_independent_third_source":False,
      "accessibility_is_TF_specific_evidence":False,
      "zero_accessibility_is_biologically_proven_closed":False,
      "original_expression_cells_split_unchanged":True,
      "source_disjoint_from_BIOV2_perturbation_response_data_by_assay":True,
      "independent_TF_target_gold_verified":False,
      "ENCODE_perturbation_expression_outcomes_accessed":False,
      "GCEV_predictions_fitted":False,"GCEV_biological_accuracy_compared":False,
      "file_sha256":{p.name:digest(p) for p in (
        out/"k562_dnase_promoter_accessibility_64.csv",
        out/"k562_tf_target_accessibility_and_motif_context_1008.csv")}
    }
    (out/"biov6_source_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in (
        "stage","n_promoters","n_directed_edges",
        "promoters_with_peak_both_replicates",
        "mean_coverage_over_promoters","nonzero_joint_prior_scores")},indent=2))
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--biov5",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    run(args.biov5,args.output)


if __name__=="__main__":main()
