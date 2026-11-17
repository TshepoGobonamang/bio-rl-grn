"""M2D6B: exact, sample-aware TF barcode-to-expression join quality gate.

Requires operator-selected columns from inspected *actual* GEO schema.
Never trims suffixes, guesses identifiers, chooses conflicting perturbations,
or promotes library/sample IDs into biological replicates.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def audit_join(map_df: pd.DataFrame, cell_df: pd.DataFrame,
               *, map_barcode: str, cell_barcode: str, tf_col: str,
               map_sample: str, cell_sample: str) -> tuple[pd.DataFrame, dict]:
    for df, required, tag in (
        (map_df, [map_barcode, tf_col, map_sample], "mapping"),
        (cell_df, [cell_barcode, cell_sample], "expression")):
        missing = set(required) - set(df.columns)
        if missing:
            raise ValueError(f"{tag} missing declared columns {sorted(missing)}")
        if df[required].isna().any().any() or any(
            df[c].astype(str).str.strip().eq("").any() for c in required):
            raise ValueError(f"{tag} contains missing or blank IDs")
    mapping = map_df[[map_sample, map_barcode, tf_col]].copy()
    mapping.columns = ["sample_id", "barcode", "tf_assignment"]
    cells = cell_df[[cell_sample, cell_barcode]].copy()
    cells.columns = ["sample_id", "barcode"]
    for df in [mapping, cells]:
        df["sample_id"] = df.sample_id.astype(str)
        df["barcode"] = df.barcode.astype(str)
    mapping["tf_assignment"] = mapping.tf_assignment.astype(str)
    if cells.duplicated(["sample_id", "barcode"]).any():
        raise ValueError("Duplicated expression sample/barcode pairs")
    unique_mapping = mapping.drop_duplicates()
    conflict = unique_mapping.groupby(["sample_id", "barcode"]).tf_assignment.nunique()
    n_conflicts = int((conflict > 1).sum())
    # Unresolved conflicts are always blocked, never arbitrarily first-assigned.
    if n_conflicts:
        raise ValueError(f"{n_conflicts} barcode mappings have conflicting TF assignments")
    unique_mapping = unique_mapping.drop_duplicates(["sample_id", "barcode"])
    merged = cells.merge(unique_mapping, how="left", on=["sample_id", "barcode"],
                         indicator=True, validate="one_to_one")
    unmapped = int(merged.tf_assignment.isna().sum())
    unreferenced = int(len(unique_mapping) - int((merged["_merge"] == "both").sum()))
    report = {
        "stage": "M2D6B_exact_sample_barcode_join_audit",
        "n_expression_cells": len(cells),
        "n_distinct_mapping_keys": len(unique_mapping),
        "n_matched_cells": int((merged["_merge"] == "both").sum()),
        "n_unmatched_expression_cells": unmapped,
        "n_unreferenced_mapping_keys": unreferenced,
        "n_ambiguous_assignments": n_conflicts,
        "all_expression_cells_mapped": unmapped == 0,
        "control_identity_verified": False,
        "biological_replicates_verified": False,
        "ready_for_perturbation_inference": False,
    }
    merged = merged.drop(columns=["_merge"])
    return merged, report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tfmap", type=Path, required=True)
    p.add_argument("--cells", type=Path, required=True,
                   help="Actual H5AD obs cell-index export with sample ID, not expression values")
    p.add_argument("--map-barcode", required=True)
    p.add_argument("--cell-barcode", required=True)
    p.add_argument("--tf-column", required=True)
    p.add_argument("--map-sample", required=True)
    p.add_argument("--cell-sample", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    mapping = pd.read_csv(a.tfmap, dtype=str)
    cells = pd.read_csv(a.cells, dtype=str)
    joined, report = audit_join(mapping, cells,
        map_barcode=a.map_barcode, cell_barcode=a.cell_barcode,
        tf_col=a.tf_column, map_sample=a.map_sample, cell_sample=a.cell_sample)
    for key, file in [("mapping_sha256", a.tfmap), ("cells_sha256", a.cells)]:
        report[key] = hashlib.sha256(file.read_bytes()).hexdigest()
    a.output.mkdir(parents=True, exist_ok=True)
    joined.to_csv(a.output / "barcode_join_audit.csv", index=False)
    (a.output / "barcode_join_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
