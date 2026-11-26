"""BIO-V3: source-audited matrix orientation and marker diagnostics for K562–BT142.

Records marker statistics by *original numeric CellType*, never autoassigns
numeric labels to K562. Requires original authors' code and direct marker
consistency review before accepting a K562 training subset.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from itertools import chain
from pathlib import Path


MARKERS=("HBG2","PTPRZ1")


def audit(expression: Path, metadata: Path) -> dict:
    with gzip.open(metadata, "rt", encoding="utf-8-sig",newline="") as f:
        meta_reader=csv.DictReader(f,delimiter="\t")
        if not {"orig.ident","CellType"}.issubset(set(meta_reader.fieldnames or [])):
            raise ValueError("Original metadata must expose exact 'orig.ident' and 'CellType'")
        label_for={}
        for row in meta_reader:
            cell=(row["orig.ident"] or "").strip()
            label=(row["CellType"] or "").strip()
            if not cell or not label or cell in label_for:
                raise ValueError("Missing or duplicated cell ID / numeric cell identity")
            label_for[cell]=label
    with gzip.open(expression,"rt", encoding="utf-8-sig", newline="") as f:
        rows=csv.reader(f,delimiter="\t")
        header=next(rows)
        first=next(rows)
        if len(header)!=len(set(header)) or len(first)!=len(header)+1:
            raise ValueError("Gene-by-cell matrix orientation not confirmed: expected first data row gene ID + cell values")
        matched=[(i, cell) for i,cell in enumerate(header) if cell in label_for]
        if len(matched)<0.95*len(header):
            raise ValueError("Too many cells lack an exact match to original metadata")
        cells_by_label={}
        for i,cell in matched:
            cells_by_label.setdefault(label_for[cell],[]).append(i)
        if any(len(ix)<2 for ix in cells_by_label.values()):
            raise ValueError("Insufficient cells in a source metadata category")
        found={}
        for record in chain((first,),rows):
            if not record: continue
            gene=record[0]
            if gene not in MARKERS:continue
            if len(record)!=len(header)+1:
                raise ValueError(f"Expression row shape differs for marker {gene}")
            x=[]
            for value in record[1:]:
                try: v=float(value)
                except ValueError as e:raise ValueError(f"Non-numeric {gene} source abundance") from e
                if not math.isfinite(v) or v<0:
                    raise ValueError(f"Invalid {gene} original abundance")
                x.append(v)
            found[gene]={label:{
                "n_cells":len(indices),
                "mean_abundance":sum(x[i] for i in indices)/len(indices),
                "fraction_positive":sum(x[i]>0 for i in indices)/len(indices)
            } for label,indices in cells_by_label.items()}
            if len(found)==len(MARKERS):break
    return {
      "source":"GSM5534702 original GEO expression + metadata",
      "expression_sha256":hashlib.sha256(expression.read_bytes()).hexdigest(),
      "metadata_sha256":hashlib.sha256(metadata.read_bytes()).hexdigest(),
      "n_exactly_joined_cells":len(matched),
      "n_unmatched_expression_cells":len(header)-len(matched),
      "n_unmatched_metadata_cells":len(set(label_for)-set(header)),
      "n_metadata_cells":len(label_for),
      "celltype_numeric_label_counts":{k:len(v) for k,v in cells_by_label.items()},
      "gene_by_cell_orientation_observed":True,
      "expression_first_header_is_cell_id":True,
      "row_extra_gene_index_present":True,
      "marker_expression_by_original_numeric_label":found,
      "author_original_code":"https://github.com/petervangalen/MAESTER-2021/blob/fe8001bfd242d4440dcb72be42a28f346adbc047/2_Cell_line_mixes_decontX/2.2_TenX_CellLineMix_decontX.R",
      "author_cluster_labeling":"RenameIdents: cluster 0->BT142, cluster 1->K562; factor levels c(BT142,K562). Does not automatically certify numeric GEO export labels.",
      "label_to_K562_assigned":False,
      "fit_model_performed":False,
      "encode_crispri_outcomes_accessed":False,
      "decision":"SOURCE_MARKER_DIAGNOSTICS_ONLY_AUTHOR_EXPORT_MAPPING_REVIEW_REQUIRED"
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expression",required=True,type=Path)
    p.add_argument("--metadata",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    x=p.parse_args()
    result=audit(x.expression,x.metadata)
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
