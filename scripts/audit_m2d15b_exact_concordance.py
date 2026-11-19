"""M2D15B: real headerless TF-map vs H5AD exact barcode evidence audit.

No implicit barcode transformations, culture assumptions or control inference.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

def audit_records(mapping: pd.DataFrame, obs: pd.DataFrame):
    if mapping.shape[1] != 2:
        raise ValueError("Expected exactly 2 TFmap columns; manual inspection required")
    m = mapping.copy()
    m.columns = ["barcode", "map_tf"]
    m["barcode"] = m["barcode"].astype("string")
    m["map_tf"] = m["map_tf"].astype("string")
    if m["barcode"].isna().any() or m["barcode"].eq("").any():
        raise ValueError("TFmap empty barcodes")
    if m["barcode"].duplicated().any():
        raise ValueError("Duplicate TFmap barcodes; cannot silently resolve")
    if not obs.index.is_unique or "TF" not in obs:
        raise ValueError("H5AD must have unique cell index and TF field")
    h = pd.DataFrame({"barcode":pd.Index(obs.index).astype("string"),"h5ad_tf":obs["TF"].astype("string").to_numpy()})
    j = h.merge(m,on="barcode",how="left",validate="one_to_one",indicator=True)
    matched = j["_merge"].eq("both")
    comparable = matched & j["map_tf"].notna() & j["h5ad_tf"].notna()
    agree = (j["map_tf"].eq(j["h5ad_tf"]) & comparable).fillna(False)
    map_barcodes = set(m["barcode"])
    h_barcodes = set(h["barcode"])
    return {
        "n_h5ad_cells":len(h),
        "n_tfmap_records_headerless":len(m),
        "tfmap_unique_barcodes":m["barcode"].nunique(),
        "tfmap_missing_assignment":int(m["map_tf"].isna().sum()),
        "tfmap_distinct_nonmissing_labels":int(m["map_tf"].nunique(dropna=True)),
        "exact_matched_cells":int(matched.sum()),
        "exact_unmatched_h5ad":int((~matched).sum()),
        "exact_unmatched_tfmap":len(map_barcodes-h_barcodes),
        "exact_comparable_tf_labels":int(comparable.sum()),
        "exact_concordant_tf_labels":int(agree.sum()),
        "exact_discordant_tf_labels":int(comparable.sum()-agree.sum()),
        "exact_concordance_fraction":float(agree.sum()/comparable.sum()) if comparable.sum() else None,
        "sample_aware_join":False,
        "barcode_transform_applied":False,
        "controls_verified":False,
        "culture_replicates_verified":False,
        "independent_regulatory_labels_verified":False,
        "decision":"EXACT_DIAGNOSTIC_ONLY",
        "caution":"Even exact matches may not be unique across experimental samples; full sample provenance required."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tfmap",required=True,type=Path)
    p.add_argument("--h5ad",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    from scripts.reconcile_m2d15_tfmap_h5ad import inspect_h5ad
    m=pd.read_csv(a.tfmap,header=None,compression="infer",dtype="string")
    o=inspect_h5ad(a.h5ad)
    d=audit_records(m,o)
    for name,path in (("tfmap_sha256",a.tfmap),("h5ad_sha256",a.h5ad)):
        h=hashlib.sha256()
        with path.open("rb") as f:
            for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
        d[name]=h.hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2)+"\n")
    print(json.dumps(d,indent=2))
if __name__=="__main__":main()
