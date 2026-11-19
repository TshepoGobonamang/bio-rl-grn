"""M2D15: source-aware, explicit-key TFmap/H5AD reconciliation.

No suffix stripping, implicit sample mapping, inferred control labels or
biological replicates. Exact cell + sample keys only when explicitly configured.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
import pandas as pd

def normalized(x):
    return x.astype("string").str.strip()

def reconcile(obs: pd.DataFrame, mapping: pd.DataFrame,
              obs_cell: str, map_cell: str, obs_tf: str, map_tf: str,
              obs_sample: str | None = None, map_sample: str | None = None) -> dict:
    if bool(obs_sample) != bool(map_sample):
        raise ValueError("Sample keys must be provided on both sides")
    for df, cols, label in [(obs,[obs_cell,obs_tf]+([obs_sample] if obs_sample else []),"H5AD"),
                             (mapping,[map_cell,map_tf]+([map_sample] if map_sample else []),"TFmap")]:
        if len(set(cols)) != len(cols):
            raise ValueError(f"{label}: column roles must be distinct")
        absent=[s for s in cols if s not in df.columns]
        if absent: raise ValueError(f"{label} columns absent: {absent}")
    lk=["cell"]+(["sample"] if obs_sample else [])
    a=pd.DataFrame({"cell":normalized(obs[obs_cell]),"h5ad_tf":normalized(obs[obs_tf])})
    b=pd.DataFrame({"cell":normalized(mapping[map_cell]),"map_tf":normalized(mapping[map_tf])})
    if obs_sample:
        a["sample"]=normalized(obs[obs_sample])
        b["sample"]=normalized(mapping[map_sample])
    for label,frame in [("H5AD",a),("TFmap",b)]:
        if frame[lk].isna().any().any() or (frame[lk]=="").any().any():
            raise ValueError(f"{label} contains blank join identifiers")
    # Distinct source records and conflicting labels must not be collapsed silently.
    if a.duplicated(lk).any():
        raise ValueError("Duplicate H5AD join keys")
    duplicates=b.duplicated(lk,keep=False)
    conflicts=(b.loc[duplicates].groupby(lk,dropna=False)["map_tf"].nunique(dropna=False)>1)
    if conflicts.any():
        raise ValueError(f"Conflicting TFmap TF assignments for {int(conflicts.sum())} join keys")
    identical_duplicate_rows=int(b.duplicated(lk).sum())
    b=b.drop_duplicates(lk,keep="first")
    merged=a.merge(b,on=lk,how="left",validate="one_to_one",indicator=True)
    matched=merged["_merge"].eq("both")
    comparable=matched & merged["map_tf"].notna() & merged["h5ad_tf"].notna()
    concordant=comparable & merged["map_tf"].eq(merged["h5ad_tf"]).fillna(False)
    mismatch=comparable & ~merged["map_tf"].eq(merged["h5ad_tf"]).fillna(False)
    missing_map_labels=matched & merged["map_tf"].isna()
    unused=b.merge(a[lk],on=lk,how="left",indicator=True)["_merge"].eq("left_only").sum()
    return {
        "join_type":"EXACT_SAMPLE_CELL" if obs_sample else "EXACT_CELL_ONLY",
        "sample_aware":bool(obs_sample),
        "n_h5ad_cells":int(len(a)),"n_tfmap_unique_keys":int(len(b)),
        "tfmap_duplicate_same_key_rows_dropped":identical_duplicate_rows,
        "matched_keys":int(matched.sum()),"unmatched_h5ad":int((~matched).sum()),
        "unmatched_tfmap":int(unused),
        "missing_tfmap_label_at_match":int(missing_map_labels.sum()),
        "comparable_labels":int(comparable.sum()),
        "concordant_labels":int(concordant.sum()),
        "discordant_labels":int(mismatch.sum()),
        "fraction_h5ad_joined":float(matched.mean()) if len(a) else 0,
        "fraction_concordant_among_comparable":float(concordant.sum()/comparable.sum()) if comparable.sum() else None,
        "matched_control_verified":False,"biological_replicates_verified":False,
        "direct_regulatory_edges_verified":False,
        "warning": "Exact matching does not prove independent cell assignment or experimental control design."
    }

def inspect_h5ad(source:Path):
    import anndata
    with tempfile.TemporaryDirectory() as t:
        actual=source
        if source.suffix==".gz":
            actual=Path(t)/"source.h5ad"
            with gzip.open(source,"rb") as inp, actual.open("wb") as out:
                shutil.copyfileobj(inp,out)
        ad=anndata.read_h5ad(actual,backed="r")
        try:
            o=ad.obs.copy()
            o.insert(0,"_h5ad_index_exact",o.index.astype(str))
            return o
        finally:
            ad.file.close()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--h5ad",type=Path,required=True)
    p.add_argument("--tfmap",type=Path,required=True)
    p.add_argument("--h5ad-cell-column",default="_h5ad_index_exact")
    p.add_argument("--tfmap-cell-column",required=True)
    p.add_argument("--h5ad-tf-column",default="TF")
    p.add_argument("--tfmap-tf-column",required=True)
    p.add_argument("--h5ad-sample-column")
    p.add_argument("--tfmap-sample-column")
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    observed=inspect_h5ad(x.h5ad)
    mapping=pd.read_csv(x.tfmap,compression="infer",dtype="string")
    result=reconcile(observed,mapping,x.h5ad_cell_column,x.tfmap_cell_column,x.h5ad_tf_column,x.tfmap_tf_column,x.h5ad_sample_column,x.tfmap_sample_column)
    for key,path in [("h5ad_sha256",x.h5ad),("tfmap_sha256",x.tfmap)]:
        h=hashlib.sha256()
        with path.open("rb") as src:
            for block in iter(lambda:src.read(1024*1024),b""):h.update(block)
        result[key]=h.hexdigest()
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
