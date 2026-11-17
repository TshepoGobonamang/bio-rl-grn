"""M2D5: conservative perturbation intake and paired replicate response summaries.

Inputs are NORMALIZED, AUDITED CSV exports (not raw GEO h5ad); the converter
for actual GSE216595 remains blocked until its exact obs/var schema, controls
and genuine biological replication are verified.

Never promote cells or 10x sequencing lanes into biological replicates.
No binary response/nonresponse labels are generated here.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

METADATA = {"cell_id", "tf", "condition", "replicate_id", "batch"}
CONTROL, TREATMENT = "control", "overexpression"


def inspect_metadata(meta: pd.DataFrame) -> dict:
    if not METADATA.issubset(meta.columns):
        raise ValueError(f"Missing metadata: {sorted(METADATA - set(meta.columns))}")
    if meta[list(METADATA)].isna().any().any():
        raise ValueError("Missing cell, condition, batch or biological replicate assignment")
    if meta.cell_id.duplicated().any():
        raise ValueError("Duplicate cell identifiers")
    if not meta.condition.isin([CONTROL, TREATMENT]).all():
        raise ValueError("Only explicit control/overexpression assignments permitted")
    if not all(meta[c].astype(str).str.strip().ne("").all() for c in METADATA):
        raise ValueError("Blank metadata identifiers")
    # No assumption that a sequencing library (GSM) is a biological replicate.
    return {"n_cells": len(meta), "n_tfs": int(meta.tf.nunique()),
            "n_batches": int(meta.batch.nunique()),
            "n_declared_replicates": int(meta.replicate_id.nunique()),
            "replicate_independence_verified": False,
            "control_identity_verified": False}


def paired_responses(counts: pd.DataFrame, meta: pd.DataFrame,
                     *, biological_replicates_verified: bool = False,
                     min_cells: int = 10, min_pairs: int = 2) -> tuple[pd.DataFrame, dict]:
    """Pseudobulk by TF/condition/batch/biological replicate, pair controls.

    Each TF's controls must be explicitly annotated with the same TF and same
    batch/replicate ID. IDs must denote verified *biological* replicates, not
    split wells, cell barcodes, technical lanes or samples from a single culture.
    No p-values, significance calls or negative-label assignments.
    """
    audit = inspect_metadata(meta)
    if not biological_replicates_verified:
        raise ValueError("Biological replicate identity is not independently verified")
    if not counts.index.is_unique or counts.columns.has_duplicates:
        raise ValueError("Count matrix must have unique cells and genes")
    if set(counts.index) != set(meta.cell_id):
        raise ValueError("Count cells must exactly match metadata cell IDs")
    x = counts.loc[meta.cell_id].to_numpy(dtype=float)
    if not np.isfinite(x).all() or (x < 0).any() or not np.allclose(x, np.rint(x)):
        raise ValueError("Expected raw finite nonnegative integer UMI counts")
    if (x.sum(axis=1) <= 0).any():
        raise ValueError("Zero-library cell")
    ordered = counts.loc[meta.cell_id].copy()
    ordered.index = pd.RangeIndex(len(ordered))
    annotated = meta.reset_index(drop=True)
    combined = pd.concat([annotated, ordered], axis=1)
    gene_cols = list(counts.columns)
    if METADATA.intersection(gene_cols):
        raise ValueError("Gene columns conflict with metadata columns")
    keys = ["tf", "batch", "replicate_id", "condition"]
    cell_sizes = combined.groupby(keys, sort=True).size()
    totals = combined.groupby(keys, sort=True)[gene_cols].sum()
    valid = cell_sizes >= min_cells
    totals = totals.loc[valid]
    cell_sizes = cell_sizes.loc[valid]
    lib = totals.sum(axis=1)
    totals = totals.loc[lib > 0]
    log_cpm = np.log2(1 + totals.div(totals.sum(axis=1), axis=0) * 1e6)
    records = []
    for tf in sorted(annotated.tf.unique()):
        try:
            tf_rows = log_cpm.xs(tf, level="tf")
        except KeyError:
            continue
        for batch in sorted(annotated.loc[annotated.tf.eq(tf), "batch"].unique()):
            try:
                bt = tf_rows.xs(batch, level="batch")
            except KeyError:
                continue
            pairs = []
            for rep in sorted(annotated.loc[annotated.tf.eq(tf) & annotated.batch.eq(batch),
                                            "replicate_id"].unique()):
                if (rep, CONTROL) not in bt.index or (rep, TREATMENT) not in bt.index:
                    continue
                pairs.append((rep, bt.loc[(rep, TREATMENT)] - bt.loc[(rep, CONTROL)]))
            if len(pairs) < min_pairs:
                continue
            diffs = pd.DataFrame([d for _, d in pairs], index=[str(rep) for rep, _ in pairs])
            for gene in gene_cols:
                vals = diffs[gene].to_numpy(dtype=float)
                records.append({"source": tf, "target": gene, "batch": batch,
                                "n_biological_pairs": len(pairs),
                                "mean_log2_cpm_difference": float(np.mean(vals)),
                                "std_log2_cpm_difference": float(np.std(vals, ddof=1)),
                                "sign_agreement_fraction": float(max(np.mean(vals > 0),
                                                                      np.mean(vals < 0))),
                                "tested": True, "binary_response_label": pd.NA})
    result = pd.DataFrame.from_records(records, columns=[
        "source", "target", "batch", "n_biological_pairs",
        "mean_log2_cpm_difference", "std_log2_cpm_difference",
        "sign_agreement_fraction", "tested", "binary_response_label"])
    audit.update({"n_paired_response_rows": len(result),
                  "biological_replicates_verified": True,
                  "control_identity_verified": False,
                  "labels_generated": False,
                  "status": "effect_summaries_only_pending_assay_controls_and_provenance_review"})
    return result, audit


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--biological-replicates-verified", action="store_true",
                   help="Explicit assertion after checking actual independent culture replicates")
    p.add_argument("--min-cells", type=int, default=10)
    p.add_argument("--min-pairs", type=int, default=2)
    a = p.parse_args()
    if a.min_cells < 1 or a.min_pairs < 2:
        p.error("min-cells must be >=1 and min-pairs >=2")
    meta = pd.read_csv(a.metadata, dtype=str)
    counts = pd.read_csv(a.counts, index_col=0)
    counts.index = counts.index.astype(str)
    effects, audit = paired_responses(
        counts, meta, biological_replicates_verified=a.biological_replicates_verified,
        min_cells=a.min_cells, min_pairs=a.min_pairs)
    audit["counts_sha256"] = hashlib.sha256(a.counts.read_bytes()).hexdigest()
    audit["metadata_sha256"] = hashlib.sha256(a.metadata.read_bytes()).hexdigest()
    a.output.mkdir(parents=True, exist_ok=True)
    effects.to_csv(a.output / "m2d5_effect_summaries.csv", index=False)
    (a.output / "m2d5_intake_manifest.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
