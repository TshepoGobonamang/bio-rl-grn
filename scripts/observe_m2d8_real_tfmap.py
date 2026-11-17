"""M2D8: inspect actual GSE216595 TF-map records and emit privacy-safe summaries.

A schema audit is not experimental control or replicate certification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def audit_source(path: Path) -> dict:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    # Preserve barcode strings, leading zeroes and original column names.
    frame = pd.read_csv(path, dtype="string", compression="infer")
    if frame.empty or frame.columns.has_duplicates:
        raise ValueError("Empty source or duplicate column names")
    cols = list(map(str, frame.columns))
    result = {
        "geo_accession": "GSM6681047",
        "file_name": path.name,
        "sha256": h.hexdigest(),
        "n_rows": int(len(frame)),
        "n_columns": len(cols),
        "columns": cols,
        "field_profiles": {
            col: {
                "non_missing": int(frame[col].notna().sum()),
                "unique_non_missing": int(frame[col].nunique(dropna=True)),
                "empty_string_count": int(frame[col].fillna("").str.strip().eq("").sum()),
            }
            for col in cols
        },
        "real_file_parsed": True,
        "tf_column_verified": False,
        "control_identity_verified": False,
        "replicate_independence_verified": False,
        "barcode_to_h5ad_join_verified": False,
        "source_provenance_cleared": False,
        "biological_labels_generated": False,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit_source(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
