"""M2D4: label-blind candidate universe, frozen before perturbation data access.

TF catalog MUST be externally predefined and NOT derived from CollecTRI or gold.
Input expression: cells x measured genes, first CSV column unique cell identifier.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def candidate_universe(expression: pd.DataFrame, tf_catalog: list[str],
                       n_genes: int = 128, n_tfs: int = 24) -> tuple[pd.DataFrame, dict]:
    if expression.empty or expression.columns.has_duplicates or expression.index.has_duplicates:
        raise ValueError("Expression needs unique cell and gene IDs")
    if not all(isinstance(x, str) and x.strip() for x in expression.columns):
        raise ValueError("Gene symbols must be nonempty strings")
    if not np.isfinite(expression.to_numpy(dtype=float)).all():
        raise ValueError("Expression contains nonfinite values")
    if len(tf_catalog) != len(set(tf_catalog)) or not tf_catalog:
        raise ValueError("TF catalog must have unique nonempty genes")
    if n_tfs < 1 or n_genes <= n_tfs:
        raise ValueError("Need a positive TF count smaller than selected genes")
    variance = expression.var(axis=0, ddof=1)
    measured = set(expression.columns)
    eligible = [t for t in tf_catalog if t in measured and variance[t] > 0]
    ranked_tfs = sorted(eligible, key=lambda g: (-float(variance[g]), g))[:n_tfs]
    if len(ranked_tfs) != n_tfs:
        raise ValueError("Insufficient variable TFs in predefined catalog")
    ranked = sorted((g for g in expression.columns if variance[g] > 0 and g not in ranked_tfs),
                    key=lambda g: (-float(variance[g]), g))
    target_genes = ranked[:n_genes - n_tfs]
    if len(target_genes) != n_genes - n_tfs:
        raise ValueError("Insufficient genes")
    genes = sorted(set(ranked_tfs + target_genes))
    edges = pd.DataFrame([(s, t) for s in sorted(ranked_tfs) for t in genes if s != t],
                         columns=["source", "target"])
    digest = hashlib.sha256(edges.to_csv(index=False).encode()).hexdigest()
    return edges, {"n_cells": len(expression), "n_measured_genes": len(measured),
                   "n_selected_genes": len(genes), "n_tfs": len(ranked_tfs),
                   "selected_tfs": sorted(ranked_tfs),
                   "n_edges": len(edges), "candidate_sha256": digest,
                   "selection": "expression variance only; predeclared TF catalog",
                   "external_labels_accessed": False,
                   "independent_validation_certified": False}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--expression", type=Path, required=True)
    p.add_argument("--tf-catalog", type=Path, required=True)
    p.add_argument("--tf-catalog-origin", required=True,
                   help="Documented non-prior/non-gold TF catalog provenance")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--n-genes", type=int, default=128)
    p.add_argument("--n-tfs", type=int, default=24)
    args = p.parse_args()
    if not args.tf_catalog_origin.strip():
        raise ValueError("Catalog origin must be documented")
    expr = pd.read_csv(args.expression, index_col=0)
    if not expr.index.is_unique:
        raise ValueError("Duplicate cell IDs")
    catalog = pd.read_csv(args.tf_catalog)
    if "gene" not in catalog or catalog.gene.isna().any():
        raise ValueError("Catalog requires nonmissing gene column")
    edges, manifest = candidate_universe(expr, catalog.gene.astype(str).tolist(),
                                         args.n_genes, args.n_tfs)
    manifest["tf_catalog_origin"] = args.tf_catalog_origin
    manifest["expression_sha256"] = hashlib.sha256(args.expression.read_bytes()).hexdigest()
    manifest["tf_catalog_sha256"] = hashlib.sha256(args.tf_catalog.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    edges.to_csv(args.output / "m2d4_frozen_candidates.csv", index=False)
    (args.output / "m2d4_candidate_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
