"""M2D6: inspect actual GEO TF barcode table and optional compressed h5ad.

Audit only: no TF response labels, invented controls, or replicate assertions.
Handles source data without committing human-cell expression to Git.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inspect_tfmap(path: Path) -> dict:
    # All fields are strings: cell barcodes, constructs and sample IDs are identifiers.
    df = pd.read_csv(path, dtype=str, compression="infer", low_memory=False)
    if df.empty:
        raise ValueError("TF mapping file empty")
    names = [str(c) for c in df.columns]
    if len(names) != len(set(names)):
        raise ValueError("Duplicated TF-map column names")
    return {
        "filename": path.name, "sha256": sha256(path),
        "n_rows": len(df), "n_columns": len(names),
        "columns": names,
        "missing_per_column": {str(k): int(v) for k, v in df.isna().sum().items()},
        "distinct_per_column": {str(k): int(v) for k, v in df.nunique().items()},
        "barcode_key_verified": False,
        "tf_assignment_key_verified": False,
        "control_identity_verified": False,
        "biological_replicates_verified": False,
    }


def inspect_h5ad(path: Path) -> dict:
    try:
        import anndata as ad
    except ImportError as ex:
        raise RuntimeError("Install anndata to inspect h5ad") from ex
    import gzip
    import tempfile
    import shutil
    # Decompress disk-backed; no dense expression array or count matrix loaded.
    if path.name.endswith(".h5ad.gz"):
        with tempfile.TemporaryDirectory() as td:
            extracted = Path(td) / "source.h5ad"
            with gzip.open(path, "rb") as src, extracted.open("wb") as out:
                shutil.copyfileobj(src, out, length=1024 * 1024)
            return _inspect_h5ad_backed(extracted, source_hash=sha256(path), name=path.name)
    return _inspect_h5ad_backed(path, source_hash=sha256(path), name=path.name)


def _inspect_h5ad_backed(path: Path, source_hash: str, name: str) -> dict:
    import anndata as ad
    data = ad.read_h5ad(path, backed="r")
    try:
        obs = data.obs
        var = data.var
        return {
            "filename": name, "sha256": source_hash,
            "n_cells": int(data.n_obs), "n_genes": int(data.n_vars),
            "obs_columns": list(map(str, obs.columns)),
            "var_columns": list(map(str, var.columns)),
            "layers": list(map(str, data.layers.keys())),
            "obsm": list(map(str, data.obsm.keys())),
            "obs_index_unique": bool(obs.index.is_unique),
            "var_index_unique": bool(var.index.is_unique),
            "cell_ids_hash": hashlib.sha256("\n".join(map(str, obs.index)).encode()).hexdigest(),
            "gene_ids_hash": hashlib.sha256("\n".join(map(str, var.index)).encode()).hexdigest(),
            "control_identity_verified": False,
            "biological_replicates_verified": False,
        }
    finally:
        data.file.close()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tfmap", type=Path, required=True)
    p.add_argument("--h5ad", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    audit = {"stage": "M2D6_real_schema_only",
             "tfmap": inspect_tfmap(a.tfmap),
             "h5ad": inspect_h5ad(a.h5ad) if a.h5ad else None,
             "true_controls_identified": False,
             "independent_replicates_identified": False,
             "perturbation_response_labels_created": False,
             "ready_for_confirmatory_validation": False}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
