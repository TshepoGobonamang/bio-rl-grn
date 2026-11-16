"""Build gold-blind M2A cross-fitted reliability data from official BEELINE inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from biorlgrn.io import load_expression, load_tfs, sha256
from biorlgrn.m11 import split_three
from biorlgrn.pseudotime import load_pseudotime
from biorlgrn.m2_crossfit import CrossfitConfig, make_reliability_meta_dataset

def _hash_ids(index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, index))).encode()).hexdigest()

def build(raw_dir: Path, output_dir: Path, *, dataset_id: str | None = None,
          tf_list: Path | None = None, seed: int = 42, folds: int = 3,
          bins: int = 16, bootstrap_repeats: int = 3) -> dict:
    raw_dir, output_dir = Path(raw_dir), Path(output_dir)
    expression_path, pseudotime_path = (
        raw_dir / "ExpressionData.csv", raw_dir / "PseudoTime.csv"
    )
    if not expression_path.is_file() or not pseudotime_path.is_file():
        raise FileNotFoundError("Both ExpressionData.csv and PseudoTime.csv are required")
    X = load_expression(expression_path, gene_axis="rows")
    pt = load_pseudotime(pseudotime_path, X.index)
    tfs = load_tfs(tf_list, list(X.columns)) if tf_list else list(X.columns)
    outer_fit, outer_score, outer_test, split = split_three(X, seed=seed)
    cfg = CrossfitConfig(folds=folds, seed=seed, max_bins=bins,
                         bootstrap_repeats=bootstrap_repeats)
    rows, manifest = make_reliability_meta_dataset(
        outer_fit, pt.loc[outer_fit.index], tfs, cfg, dataset_id or raw_dir.name
    )
    manifest.update({
        "outer_split": split,
        "outer_score_ids_sha256": _hash_ids(outer_score.index),
        "outer_test_ids_sha256": _hash_ids(outer_test.index),
        "expression_sha256": sha256(expression_path),
        "pseudotime_sha256": sha256(pseudotime_path),
        "tf_annotation_source": str(tf_list) if tf_list else "all_measured_genes_not_gold",
        "tf_annotation_sha256": sha256(tf_list) if tf_list else None,
        "outer_score_and_test_excluded_from_dataset": True,
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(output_dir / "m2_crossfit_edges.csv", index=False)
    (output_dir / "m2_crossfit_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--dataset-id")
    p.add_argument("--tf-list", type=Path)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--folds", type=int, default=3)
    p.add_argument("--bins", type=int, default=16)
    p.add_argument("--bootstrap-repeats", type=int, default=3)
    a = p.parse_args()
    manifest = build(a.raw_dir, a.output_dir, dataset_id=a.dataset_id,
                     tf_list=a.tf_list, seed=a.seed, folds=a.folds,
                     bins=a.bins, bootstrap_repeats=a.bootstrap_repeats)
    print(json.dumps({key: manifest[key] for key in [
        "dataset_id", "n_rows", "n_candidates_per_fold",
        "n_temporal_labels_available", "n_temporal_evidence_missing"
    ]}, indent=2))

if __name__ == "__main__":
    main()
