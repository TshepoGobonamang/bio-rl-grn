"""M2C-D12 descriptive metadata audit: no inferred controls or replicates."""
import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path


def summarize(obs):
    import pandas as pd
    required = {"TF", "batch", "temp"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise ValueError(f"Missing metadata fields: {missing}")
    if not obs.index.is_unique:
        raise ValueError("Cell IDs are not unique")
    fields = ["TF", "batch", "temp"]
    frame = obs[fields].copy()
    for key in fields:
        frame[key] = frame[key].astype("string").fillna("<MISSING>")
    counts = {}
    for key in fields:
        c = frame.groupby(key, dropna=False, observed=True).size().sort_values(ascending=False)
        counts[key] = [{"value": str(v), "cells": int(n)} for v, n in c.items()]
    joint = frame.groupby(fields, dropna=False, observed=True).size().sort_values(ascending=False)
    tf_batch = frame.groupby(["TF", "batch"], dropna=False, observed=True).size()
    tf_diversity = tf_batch.groupby(level=0).size().sort_values(ascending=False)
    return {
        "n_cells": len(frame), "observed_fields": fields,
        "counts_by_field": counts,
        "tf_batch_temp_cells": [{"TF": str(t), "batch": str(b), "temp": str(x), "cells": int(n)}
                                for (t, b, x), n in joint.items()],
        "batch_count_per_tf": [{"TF": str(k), "n_batches": int(v)} for k, v in tf_diversity.items()],
        "possible_control_label_candidates_not_verified": [
            row["value"] for row in counts["TF"]
            if any(term in row["value"].casefold() for term in
                   ("gfp", "control", "empty", "vector", "ntc", "non-target", "mock"))
        ],
        "tf_assignments_validated_against_tfmap": False,
        "experimental_controls_verified": False,
        "biological_replicates_verified": False,
        "differential_expression_computed": False,
        "regulatory_edge_truth_inferred": False,
    }


def audit_h5ad(source: Path):
    import anndata
    h = hashlib.sha256()
    with source.open("rb") as stream:
        for chunk in iter(lambda: stream.read(2 * 1024 * 1024), b""):
            h.update(chunk)
    with tempfile.TemporaryDirectory() as td:
        uncompressed = Path(td) / "metadata.h5ad"
        if source.suffix == ".gz":
            with gzip.open(source, "rb") as inp, uncompressed.open("wb") as out:
                shutil.copyfileobj(inp, out)
        else:
            uncompressed = source
        adata = anndata.read_h5ad(uncompressed, backed="r")
        try:
            result = summarize(adata.obs)
            result["source_sha256"] = h.hexdigest()
            result["source_name"] = source.name
            return result
        finally:
            adata.file.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    result = audit_h5ad(args.h5ad)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    print(f"Observed {result['n_cells']} cells; report saved to {args.output}")


if __name__ == "__main__":
    main()
