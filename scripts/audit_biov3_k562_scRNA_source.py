"""BIO-V3: inspect public mixed K562/BT142 GEO single-cell sources without training.

Never infer a cell type from a filename. Only report schema, source hashes,
explicit metadata label frequencies and exact identifier overlap. Do not open
ENCODE perturbation responses or model reference edges.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def delim(sample: str) -> str:
    counts={c:sample.count(c) for c in ("\t", ",", ";")}
    best=max(counts,key=counts.get)
    if counts[best] < 1:
        raise ValueError("GEO text file has no recognized field separators")
    return best


def read_header(path: Path) -> tuple[list[str], str]:
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        line=f.readline()
    if not line:
        raise ValueError("Empty GEO supplementary data")
    sep=delim(line)
    header=next(csv.reader([line],delimiter=sep))
    if not header or len(set(header))!=len(header):
        raise ValueError("Empty or duplicated source header")
    return header,sep


def audit(expression: Path, metadata: Path) -> dict:
    expr_head, expr_sep=read_header(expression)
    meta_head, meta_sep=read_header(metadata)
    if len(expr_head)<2 or len(meta_head)<2:
        raise ValueError("Both matrix and metadata need multiple fields")
    low_card={k:Counter() for k in meta_head}
    values={k:set() for k in meta_head}
    n_meta=0
    with gzip.open(metadata,"rt",encoding="utf-8-sig",newline="") as f:
        rows=csv.DictReader(f,delimiter=meta_sep)
        for row in rows:
            n_meta+=1
            for k in meta_head:
                v=(row.get(k) or "").strip()
                if v:
                    low_card[k][v]+=1
                    if len(values[k])<=100000:
                        values[k].add(v)
    if n_meta==0: raise ValueError("Empty GEO metadata rows")
    n_expr=0
    with gzip.open(expression,"rt",encoding="utf-8-sig",newline="") as f:
        next(f)
        for line in f:
            if line.strip():n_expr+=1
    expr_header=set(expr_head[1:])
    overlap={k:len(v & expr_header) for k,v in values.items()}
    candidate_cell_fields=[k for k in meta_head if
        any(word in k.lower().replace("_","").replace(" ","") for word in
            ("cellline","celltype","identity","sampletype","origin"))]
    label_fields={}
    for k in candidate_cell_fields:
        if len(low_card[k])<=30:
            label_fields[k]={
                "nonmissing":sum(low_card[k].values()),
                "values":[{"value":x,"rows":int(n)} for x,n in low_card[k].most_common(30)]
            }
    k562_labels={k:int(sum(n for v,n in low_card[k].items() if v.strip().upper()=="K562"))
                 for k in candidate_cell_fields}
    matched_keys=max(overlap.values()) if overlap else 0
    return {
      "stage":"ORIGINAL_GEO_SOURCE_SCHEMA_INSPECTION_ONLY",
      "original_sample":"GSM5534702",
      "expression_file":expression.name,"expression_sha256":sha256(expression),
      "metadata_file":metadata.name,"metadata_sha256":sha256(metadata),
      "matrix_shape_as_text_rows_after_first":int(n_expr),
      "matrix_header_column_count":len(expr_head),
      "matrix_header_first_field":expr_head[0],
      "metadata_rows":n_meta,"metadata_columns":meta_head,
      "metadata_categorical_cell_line_candidates":label_fields,
      "metadata_exact_k562_label_count_by_candidate_field":k562_labels,
      "metadata_fields_matching_expression_header":overlap,
      "largest_exact_cell_key_header_overlap":matched_keys,
      "source_expression_cell_orientation_verified":False,
      "independent_training_expression_qualified":False,
      "independent_regulatory_prior_qualified":False,
      "encode_perturbation_response_accessed":False,
      "decision":"SCHEMA_OBSERVED_MANUAL_CELL_LINE_AND_JOIN_REVIEW_REQUIRED",
      "note":"No input is interpreted as a K562 training subset without exact per-cell labels and ID joins; GEO source may mix K562 and BT142."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expression",type=Path,required=True)
    p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    r=audit(x.expression,x.metadata)
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps(r,indent=2))

if __name__=="__main__":
    main()
