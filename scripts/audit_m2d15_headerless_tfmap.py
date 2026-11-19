"""Audit whether source TFmap CSV is headerless, without guessing join keys."""
import argparse, hashlib, json, re
from pathlib import Path
import pandas as pd

DNA16=re.compile(r"^[ACGT]{16}$")
ORF=re.compile(r"^[A-Za-z0-9_-]+-NM_[0-9]+$")
def audit(path):
    raw=pd.read_csv(path,header=None,dtype="string",compression="infer",keep_default_na=True)
    if raw.shape[1]!=2: raise ValueError("Expected exactly two source fields; audit manually")
    first=[None if pd.isna(x) else str(x) for x in raw.iloc[0].tolist()]
    probable=bool(first[0] and DNA16.fullmatch(first[0]) and first[1] and ORF.fullmatch(first[1]))
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    return {
      "geo_accession":"GSM6681047","source_sha256":h,
      "source_records_if_headerless":int(len(raw)),
      "source_records_if_header":int(len(raw)-1),
      "first_record":dict(zip(["barcode","TF"],first)),
      "first_record_looks_like_data":probable,
      "column_0_barcode_16nt_fraction":float(raw.iloc[:,0].dropna().str.fullmatch(DNA16).mean()),
      "column_0_unique":int(raw.iloc[:,0].nunique(dropna=True)),
      "column_1_unique_nonmissing":int(raw.iloc[:,1].nunique(dropna=True)),
      "column_1_missing":int(raw.iloc[:,1].isna().sum()),
      "column_1_empty":int(raw.iloc[:,1].fillna("").str.strip().eq("").sum()),
      "join_key_validated":False,"h5ad_join_executed":False,
      "note":"A data-looking first record strongly suggests headerless CSV but does not certify sample-aware cell identity."
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",required=True,type=Path);p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();r=audit(a.source);a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n");print(json.dumps(r,indent=2))
if __name__=="__main__":main()
