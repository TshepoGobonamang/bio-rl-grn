"""BIO-V5: source-verified GENCODE gene-start promoters with actual JASPAR motif scans.

Sequence-only developmental binding-potential evidence. No ENCODE perturbation
labels, ChIP labels, fitted GCEV model or TF-target causal truth are used.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

H5AD_SHA = "14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6"
GENCODE_SHA = "e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9"
JASPAR_SHA = "ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15"
PROTOCOL = "GCEV-BIOV5-MOTIF-GRCH38-2026-11-27-001"
NULL_SHUFFLES = 128
PROMOTER_UPSTREAM = 1000
PROMOTER_DOWNSTREAM = 200
PROMOTER_LENGTH = 1201
DNA = set("ACGT")
COMP = str.maketrans("ACGT", "TGCA")
ASCII = np.full(256, 4, dtype=np.uint8)
for symbol, idx in (("A",0),("C",1),("G",2),("T",3)):
    ASCII[ord(symbol)] = idx


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def parse_gtf_attributes(text: str) -> dict:
    return dict(re.findall(r'(\w+)\s+"([^"]+)"', text))


def locate_gencode_gene_features(gtf: Path, selected: pd.DataFrame) -> dict[str,dict]:
    """Fail closed if any original gene has no unique exact ID/symbol gene row."""
    keep = {str(r.gencode_gene_id):str(r.gene_symbol) for r in selected.itertuples()}
    if len(keep)!=len(selected):
        raise ValueError("Duplicated original selected GENCODE IDs")
    found=defaultdict(list)
    with gzip.open(gtf, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)!=9 or c[2]!="gene": continue
            a=parse_gtf_attributes(c[8])
            gid=a.get("gene_id","")
            if gid not in keep: continue
            if c[0] not in {f"chr{i}" for i in range(1,23)}|{"chrX","chrY"}:
                raise ValueError("Selected gene assigned to non-primary chromosome: "+gid)
            if a.get("gene_name")!=keep[gid]:
                raise ValueError("GENCODE gene ID/symbol mismatch: "+gid)
            if c[6] not in ("+","-"):
                raise ValueError("Invalid original gene strand: "+gid)
            found[gid].append({
                "gene_symbol":keep[gid], "gencode_gene_id":gid,
                "chromosome":c[0], "gencode_strand":c[6],
                "gene_start_1based":int(c[3]),"gene_end_1based":int(c[4]),
                "gene_type":a.get("gene_type", ""),
            })
    if set(found)!=set(keep) or any(len(v)!=1 for v in found.values()):
        raise ValueError("Each selected gene requires exactly one primary GENCODE gene-level TSS")
    result={}
    for v in found.values():
        item=v[0]
        if item["gene_type"]!="protein_coding":
            raise ValueError("Non-protein-coding selected gene: "+item["gene_symbol"])
        if item["gene_start_1based"]<1 or item["gene_end_1based"]<item["gene_start_1based"]:
            raise ValueError("Invalid selected GENCODE gene interval")
        result[item["gene_symbol"]]=item
    return result


def promoter_coordinates(g: dict) -> dict:
    """0-based, half-open 1201-bp window around GENCODE gene-level TSS."""
    tss0=(g["gene_start_1based"]-1 if g["gencode_strand"]=="+" else
          g["gene_end_1based"]-1)
    if g["gencode_strand"]=="+":
        begin=tss0-PROMOTER_UPSTREAM
        end=tss0+PROMOTER_DOWNSTREAM+1
    else:
        begin=tss0-PROMOTER_DOWNSTREAM
        end=tss0+PROMOTER_UPSTREAM+1
    if begin<0 or end-begin!=PROMOTER_LENGTH:
        raise ValueError("GENCODE TSS promoter falls outside expected chromosome range")
    url=(f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
         f"chrom={g['chromosome']};start={begin};end={end}")
    return {**g, "tss_0based":tss0,"start_0based":begin,
            "end_0based_exclusive":end,"length":PROMOTER_LENGTH,
            "ucsc_reference_genome":"hg38","url":url}


def fetch_ucsc_promoter(region: dict, timeout: int=45, attempts: int=5) -> tuple[str,dict]:
    """Accept only UCSC hg38 response matching exact requested region and length."""
    payload=None
    for attempt in range(attempts):
        try:
            req=urllib.request.Request(region["url"],headers={
                "Accept":"application/json",
                "User-Agent":"GCEV-BIOV5-research-sequence-provenance/0.1"
            })
            with urllib.request.urlopen(req,timeout=timeout) as response:
                payload=response.read(100_001)
            if len(payload)>100_000: raise ValueError("UCSC response unexpectedly large")
            doc=json.loads(payload)
            if doc.get("genome")!="hg38" or doc.get("chrom")!=region["chromosome"]:
                raise ValueError("UCSC sequence source genome/chromosome mismatch")
            if int(doc.get("start",-1))!=region["start_0based"] or int(doc.get("end",-1))!=region["end_0based_exclusive"]:
                raise ValueError("UCSC sequence source coordinate mismatch")
            raw=str(doc.get("dna","")).upper()
            if len(raw)!=PROMOTER_LENGTH or set(raw)-DNA:
                raise ValueError("UCSC promoter contains unknown bases or invalid length")
            oriented=raw if region["gencode_strand"]=="+" else raw.translate(COMP)[::-1]
            return oriented, {
                **region,
                "ucsc_raw_response_sha256":digest_bytes(payload),
                "genomic_forward_dna_sha256":digest_bytes(raw.encode("ascii")),
                "gene_oriented_dna_sha256":digest_bytes(oriented.encode("ascii")),
                "n_GC":raw.count("G")+raw.count("C"),
                "GC_fraction":(raw.count("G")+raw.count("C"))/len(raw),
                "sequence_pass":True,
            }
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            if attempt==attempts-1: raise RuntimeError(
                f"UCSC original sequence request failed: {region['gene_symbol']}: {e}") from e
            time.sleep(min(2**attempt,12))
    raise RuntimeError("Original UCSC sequence unavailable")


def parse_meme_selected(path: Path, selected: pd.DataFrame) -> dict[str,list[dict]]:
    """JASPAR 2026 MEME letter-probability matrices, frozen selected IDs only."""
    require={}
    for r in selected.itertuples(index=False):
        if r.role!="regulator": continue
        require[r.gene_symbol]=str(r.motif_ids).split("|")
    if len(require)!=16 or any(not ids or any(not re.fullmatch(r"MA\d+\.\d+",i) for i in ids)
                               for ids in require.values()):
        raise ValueError("Expected 16 source-identified TFs and complete JASPAR motif IDs")
    ids_to_tf={}
    for tf, ids in require.items():
        for motif_id in ids:
            if motif_id in ids_to_tf and ids_to_tf[motif_id]!=tf:
                raise ValueError("Selected motif ID unexpectedly shared across TFs")
            ids_to_tf[motif_id]=tf
    selected_rows={}
    lines=path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("MEME version"):
        raise ValueError("JASPAR source is not in expected MEME format")
    k=0
    while k<len(lines):
        line=lines[k].strip()
        if not line.startswith("MOTIF "):
            k+=1;continue
        parts=line.split(maxsplit=2)
        if len(parts)!=3:
            raise ValueError("Unrecognized JASPAR MEME MOTIF header")
        motif_id, motif_name=parts[1:]
        if motif_id not in ids_to_tf:
            k+=1;continue
        if motif_id in selected_rows:
            raise ValueError("Selected original JASPAR motif ID appears twice")
        if motif_name.casefold()!=ids_to_tf[motif_id].casefold():
            raise ValueError(f"Motif ID/name not original TF: {motif_id} {motif_name}")
        j=k+1
        while j<len(lines) and not lines[j].strip().startswith(("letter-probability matrix:","MOTIF ")):
            j+=1
        if j>=len(lines) or not lines[j].strip().startswith("letter-probability matrix:"):
            raise ValueError("Original selected motif lacks probability matrix: "+motif_id)
        meta=lines[j]
        match=re.search(r"\bw\s*=\s*(\d+)",meta)
        if not match:
            raise ValueError("Motif width not specified: "+motif_id)
        width=int(match.group(1))
        if width<4 or width>40:
            raise ValueError("Invalid JASPAR matrix width")
        matrix=[]
        for line0 in lines[j+1:j+1+width]:
            vals=line0.split()
            if len(vals)!=4: raise ValueError("JASPAR probability matrix row shape mismatch")
            matrix.append([float(x) for x in vals])
        x=np.asarray(matrix,dtype=float)
        if x.shape!=(width,4) or not np.isfinite(x).all() or (x<0).any() or (
            np.max(np.abs(x.sum(axis=1)-1))>0.01):
            raise ValueError("Invalid JASPAR PWM frequencies")
        x=x / x.sum(axis=1,keepdims=True)
        pwm=np.log2(np.maximum(x,1e-6)/0.25)
        min_bits=float(pwm.min(axis=1).sum())
        max_bits=float(pwm.max(axis=1).sum())
        if max_bits-min_bits<=1e-8:
            raise ValueError("Uninformative constant JASPAR motif")
        selected_rows[motif_id]={
            "motif_id":motif_id,"tf_symbol":ids_to_tf[motif_id],
            "jaspar_factor_name":motif_name,"width":width,
            "probabilities":x,"pwm_bits":pwm,
            "min_bits":min_bits, "max_bits":max_bits,
            "motif_probability_sha256":digest_bytes(
                x.astype("<f8").tobytes()),
        }
        k=j+1+width
    if set(selected_rows)!=set(ids_to_tf):
        raise ValueError(f"Missing required JASPAR motifs: {sorted(set(ids_to_tf)-set(selected_rows))}")
    out=defaultdict(list)
    for motif_id, data in sorted(selected_rows.items()):
        out[data["tf_symbol"]].append(data)
    return dict(out)


def encode(seq: str) -> np.ndarray:
    b=np.frombuffer(seq.encode("ascii"),dtype=np.uint8)
    result=ASCII[b]
    if len(result)!=PROMOTER_LENGTH or np.any(result>3):
        raise ValueError("Non-ACGT or incorrectly sized promoter sequence")
    return result


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def scan_one_motif(fwd: np.ndarray, rev: np.ndarray, motif: dict) -> dict:
    """Exhaustive PWM scan in both strands; motif relative score scaled 0–1."""
    weight=motif["pwm_bits"]
    n=len(fwd)
    width=weight.shape[0]
    if len(rev)!=n or n<width:
        raise ValueError("Original sequence too short for JASPAR motif")
    best_score=-np.inf
    best_start=-1
    best_strand=None
    for strand,codes in (("+",fwd),("-",rev)):
        scores=np.zeros(n-width+1,dtype=float)
        for i in range(width):
            scores+=weight[i,codes[i:i+len(scores)]]
        best=int(np.argmax(scores))
        val=float(scores[best])
        if val>best_score:
            best_score=val
            best_start=best if strand=="+" else n-(best+width)
            best_strand=strand
    rel=float(np.clip((best_score-motif["min_bits"])/
                     (motif["max_bits"]-motif["min_bits"]),0,1))
    return {"relative":rel,"raw_bits":best_score,
            "position_0based":best_start,"strand":best_strand,
            "motif_id":motif["motif_id"]}


def dinucleotide_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Random Eulerian trail through original directed nucleotide transitions."""
    if len(seq)<2 or set(seq)-DNA: raise ValueError("Cannot shuffle invalid DNA")
    outgoing={b:[] for b in DNA}
    for left,right in zip(seq,seq[1:]):
        outgoing[left].append(right)
    for adj in outgoing.values():
        rng.shuffle(adj)
    stack=[seq[0]]
    trail=[]
    while stack:
        node=stack[-1]
        if outgoing[node]:
            stack.append(outgoing[node].pop())
        else:
            trail.append(stack.pop())
    result="".join(reversed(trail))
    if (len(result)!=len(seq) or result[0]!=seq[0] or result[-1]!=seq[-1]
        or Counter(zip(result,result[1:]))!=Counter(zip(seq,seq[1:]))):
        raise AssertionError("Null shuffle failed exact dinucleotide-preserving invariant")
    return result


def one_target_edges(target: str, seq: str, selected: pd.DataFrame,
                     motifs: dict[str,list[dict]]) -> list[dict]:
    seed=int.from_bytes(hashlib.sha256(
        (PROTOCOL+"|"+target).encode("utf-8")).digest()[:8],"big")
    rng=np.random.default_rng(seed)
    sequences=[seq]+[dinucleotide_shuffle(seq,rng) for _ in range(NULL_SHUFFLES)]
    encoded=[(encode(s),encode(reverse_complement(s))) for s in sequences]
    # Compute each unique motif once per original/null sequence; re-use across TFs.
    records={}
    for tf, group in motifs.items():
        for motif in group:
            records[motif["motif_id"]]=[scan_one_motif(a,b,motif)
                                       for a,b in encoded]
    rows=[]
    for r in selected.itertuples(index=False):
        if r.role!="regulator" or r.gene_symbol==target:
            continue
        tf=r.gene_symbol
        group=motifs[tf]
        best=max((records[m["motif_id"]][0] for m in group),
                 key=lambda x:x["relative"])
        null=np.asarray([max(records[m["motif_id"]][i]["relative"] for m in group)
                         for i in range(1,NULL_SHUFFLES+1)])
        observed=float(best["relative"])
        median=float(np.median(null))
        pval=float((1+int(np.count_nonzero(null>=observed-1e-12)))/
                   (NULL_SHUFFLES+1))
        percentile95=float(np.quantile(null,.95))
        strength=(max(0.0,(observed-median)/(1-median))
                  if median < 1-1e-10 else 0.0)
        rows.append({
            "tf":tf,"target":target,
            "tf_gencode_gene_id":r.gencode_gene_id,
            "target_gencode_gene_id":str(
                selected.loc[selected.gene_symbol==target,"gencode_gene_id"].iloc[0]),
            "prior_score":float(np.clip(strength,0,1)),
            "observed_best_pwm_relative":observed,
            "observed_best_pwm_raw_bits":best["raw_bits"],
            "null_median_pwm_relative":median,
            "null_95th_pwm_relative":percentile95,
            "observed_minus_null95":observed-percentile95,
            "higher_than_null95":bool(observed>percentile95),
            "empirical_p":pval,
            "best_motif_id":best["motif_id"],
            "best_position_0based":best["position_0based"],
            "best_strand_gene_orientation":best["strand"],
            "motif_profiles_tested":len(group),
            "null_shuffles":NULL_SHUFFLES,
            "source_observed":True,
            "direct_tf_binding_verified":False,
            "regulatory_edge_truth":False,
        })
    return rows


def bh_q(p):
    p=np.asarray(p,dtype=float)
    if not len(p) or not np.isfinite(p).all() or np.any((p<0)|(p>1)):
        raise ValueError("Invalid empirical p-values for BH reporting")
    n=len(p)
    order=np.argsort(p,kind="stable")
    ranked=p[order]*n/np.arange(1,n+1)
    q=np.minimum.accumulate(ranked[::-1])[::-1]
    out=np.empty(n,dtype=float)
    out[order]=np.minimum(q,1.0)
    return out


def source_inputs(cohort_root:Path, gtf:Path, meme:Path) -> tuple[pd.DataFrame,dict]:
    if sha(gtf)!=GENCODE_SHA or sha(meme)!=JASPAR_SHA:
        raise ValueError("Source GENCODE or JASPAR hash fails frozen protocol")
    file=cohort_root/"biov4b_fit_only_manifest.json"
    manifest=json.loads(file.read_text())
    if (manifest.get("original_K562_H5AD_sha256")!=H5AD_SHA
        or manifest.get("n_fit_cells")!=394
        or manifest.get("n_internal_score_cells")!=131
        or manifest.get("n_untouched_test_cells")!=133
        or manifest.get("counts",{}).get("directed_candidate_edges_no_self")!=1008
        or manifest.get("GENCODE_v29_sha256")!=GENCODE_SHA
        or manifest.get("JASPAR_2026_sha256")!=JASPAR_SHA):
        raise ValueError("Original BIO-V4B source selection is not the locked 16x64 cohort")
    for name,expected in manifest["sha256"].items():
        if sha(cohort_root/name)!=expected:
            raise ValueError("Original BIO-V4B file altered: "+name)
    selected=pd.read_csv(cohort_root/"selected_gene_universe.csv",dtype="string")
    if (len(selected)!=64 or selected["gene_symbol"].duplicated().any()
        or selected["gencode_gene_id"].duplicated().any()
        or int((selected["role"]=="regulator").sum())!=16
        or int((selected["role"]=="target_only").sum())!=48):
        raise ValueError("Malformed original BIO-V4B measured candidate gene universe")
    return selected, manifest


def run(cohort_root:Path, gtf:Path, meme:Path, out:Path) -> dict:
    selected,manifest=source_inputs(cohort_root,gtf,meme)
    coordinates=locate_gencode_gene_features(gtf,selected)
    motifs=parse_meme_selected(meme,selected)
    if len(motifs)!=16 or sum(map(len,motifs.values()))!=17:
        raise ValueError("Unexpected selected single TF MEME profile count")
    out.mkdir(parents=True,exist_ok=True)
    coords=[]
    sequences={}
    for row in selected.itertuples(index=False):
        gene=row.gene_symbol
        region=promoter_coordinates(coordinates[gene])
        seq,metadata=fetch_ucsc_promoter(region)
        sequences[gene]=seq
        coords.append(metadata)
        print(f"ORIGINAL_GRCH38_PROMOTER {gene} {region['chromosome']}:{region['start_0based']}-{region['end_0based_exclusive']}",flush=True)
    with (out/"gene_oriented_hg38_promoters.fa").open("w") as f:
        for row in selected.itertuples(index=False):
            g=row.gene_symbol
            meta=next(m for m in coords if m["gene_symbol"]==g)
            f.write(f">{g}|{row.gencode_gene_id}|{meta['chromosome']}:{meta['start_0based']}-{meta['end_0based_exclusive']}|{meta['gencode_strand']}\n{sequences[g]}\n")
    pd.DataFrame(coords).to_csv(out/"genecode29_promoter_coordinates_and_sources.csv",index=False)
    source_motifs=[]
    for tf,models in motifs.items():
        for m in models:
            source_motifs.append({"tf":tf,"motif_id":m["motif_id"],
                 "motif_width":m["width"],"motif_name":m["jaspar_factor_name"],
                 "matrix_probability_sha256":m["motif_probability_sha256"],
                 "pwm_min_bits":m["min_bits"],"pwm_max_bits":m["max_bits"]})
    pd.DataFrame(source_motifs).to_csv(out/"jaspar_selected_motifs_sha256.csv",index=False)
    edge_records=[]
    for row in selected.itertuples(index=False):
        edge_records+=one_target_edges(row.gene_symbol,sequences[row.gene_symbol],
                                      selected,motifs)
        print(f"SCANNED_{len(edge_records)} {row.gene_symbol}",flush=True)
    edge=pd.DataFrame(edge_records)
    if (len(edge)!=1008 or edge.duplicated(["tf","target"]).any() or
        (edge.tf==edge.target).any() or not set(edge.tf)==set(motifs) or
        not set(edge.target)==set(selected.gene_symbol)):
        raise AssertionError("Incomplete or invalid 1,008-edge promoter motif scan")
    edge["bh_q"]=bh_q(edge.empirical_p.to_numpy(float))
    if not np.isfinite(edge.select_dtypes(include="number").to_numpy(float)).all():
        raise AssertionError("Nonnumeric/nan in motif prior")
    if not edge.prior_score.between(0,1).all():
        raise AssertionError("Motif prior outside 0..1")
    edge.to_csv(out/"k562_tf_to_promoter_motif_prior_1008.csv",index=False)
    results={
        "protocol":PROTOCOL,"stage":"BIOV5_IN_SILICO_SEQUENCE_BINDING_POTENTIAL_NO_BIOLOGICAL_GOLD",
        "original_biov4b_manifest_sha256":sha(cohort_root/"biov4b_fit_only_manifest.json"),
        "source_H5AD_sha256":H5AD_SHA,"gencode_v29_gtf_sha256":GENCODE_SHA,
        "jaspar2026_CORE_MEME_sha256":JASPAR_SHA,"reference_assembly":"UCSC hg38 (GRCh38)",
        "promoter_TSS":"GENCODE v29 gene feature start/end, NOT K562-active transcript TSS",
        "promoter_upstream_bp":PROMOTER_UPSTREAM,"promoter_downstream_bp":PROMOTER_DOWNSTREAM,
        "promoter_length":PROMOTER_LENGTH,"TFs":16,"target_genes":64,
        "tested_directed_no_self_edges":len(edge),
        "motif_profiles":sum(map(len,motifs.values())),
        "null_type":"exact dinucleotide-preserving randomized Eulerian trail",
        "null_sequences_per_gene":NULL_SHUFFLES,
        "p_value_resolution":1/(NULL_SHUFFLES+1),
        "prior_method":"max relative PWM log2 odds over registered JASPAR motifs and two gene-oriented strands; positive excess above dinucleotide-null median scaled to [0,1]",
        "prior_positive_scores":int((edge.prior_score>0).sum()),
        "empirical_p_below_005":int((edge.empirical_p<.05).sum()),
        "BH_q_below_005":int((edge.bh_q<.05).sum()),
        "mean_prior_strength":float(edge.prior_score.mean()),
        "median_prior_strength":float(edge.prior_score.median()),
        "n_source_promoters_fully_retrieved":len(sequences),
        "n_MOTIF_sourced_TF_regulators":len(motifs),
        "all_candidate_scores_observed":bool(edge.source_observed.all()),
        "original_cell_partitions_retained_not_recomputed":True,
        "K562_expression_used_in_prior_score":False,
        "ENCODE_perturbation_response_labels_accessed":False,
        "direct_tf_binding_verified":False,"gcev_model_fitted":False,
        "regulatory_accuracy_or_superiority_evaluated":False,
        "reference_sources_checked_for_future_test_independence":False,
        "scientific_status":"MOTIF_BINDING_POTENTIAL_NOT_DIRECT_REGULATION",
        "file_sha256":{
            p.name:sha(p) for p in (out/"gene_oriented_hg38_promoters.fa",
                   out/"genecode29_promoter_coordinates_and_sources.csv",
                   out/"jaspar_selected_motifs_sha256.csv",
                   out/"k562_tf_to_promoter_motif_prior_1008.csv")
        }
    }
    (out/"biov5_motif_prior_manifest.json").write_text(
        json.dumps(results,indent=2)+"\n")
    print(json.dumps({key:results[key] for key in (
        "stage","tested_directed_no_self_edges","prior_positive_scores",
        "mean_prior_strength","empirical_p_below_005",
        "BH_q_below_005")},indent=2))
    return results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cohort",type=Path,required=True)
    p.add_argument("--gencode",type=Path,required=True)
    p.add_argument("--jaspar",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.cohort,args.gencode,args.jasper,args.output)


if __name__=="__main__":main()
