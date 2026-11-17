"""M2D10: inspect actual processed GSE216595 h5ad without interpreting biological design.

Outputs metadata schema/counts only, not identifiable cell-level records.
Compressed input is streamed to scratch disk before backed reading.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inspect(path: Path) -> dict:
    try:
        import anndata
    except ImportError as exc:
        raise RuntimeError("anndata is required; install with pip install anndata") from exc
    original_hash = digest(path)
    with tempfile.TemporaryDirectory() as scratch:
        local = path
        if path.name.endswith(".gz"):
            local = Path(scratch) / "source.h5ad"
            with gzip.open(path, "rb") as source, local.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
        adata = anndata.read_h5ad(local, backed="r")
        try:
            obs = adata.obs
            var = adata.var
            if not obs.index.is_unique or not var.index.is_unique:
                raise ValueError("Nonunique cell or gene identifiers")
            summary = {
                "source_name": path.name,
                "source_sha256": original_hash,
                "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
                "obs_columns": list(map(str, obs.columns)),
                "var_columns": list(map(str, var.columns)),
                "layers": list(map(str, adata.layers.keys())),
                "obs_index_unique": True, "var_index_unique": True,
                "obs_id_digest": hashlib.sha256("\n".join(map(str, obs.index)).encode()).hexdigest(),
                "var_id_digest": hashlib.sha256("\n".join(map(str, var.index)).encode()).hexdigest(),
                "obs_field_profiles": {
                    str(name): {
                        "missing": int(obs[name].isna().sum()),
                        "unique_nonmissing": int(obs[name].nunique(dropna=True)),
                    } for name in obs.columns
                },
                "tf_barcode_join_verified": False,
                "controls_verified": False,
                "independent_culture_replicates_verified": False,
                "regulatory_labels_computed": False,
            }
            return summary
        finally:
            adata.file.close()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    result = inspect(a.h5ad)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
