"""Export gold-blind M2C-A direction-sensitive cross-fitted evidence targets."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from biorlgrn.io import load_expression, load_tfs, sha256
from biorlgrn.m11 import split_three
from biorlgrn.pseudotime import load_pseudotime
from biorlgrn.m2_crossfit import CrossfitConfig
from biorlgrn.m2_direction import make_directional_meta_dataset


def build(raw_dir: Path, output_dir: Path, *, dataset_id: str | None = None,
          tf_list: Path | None = None, seed: int = 42, folds: int = 3,
          bins: int = 16, bootstrap_repeats: int = 3):
    raw_dir, out = Path(raw_dir), Path(output_dir)
    exp_file = raw_dir / "ExpressionData.csv"
    pt_file = raw_dir / "PseudoTime.csv"
    if not exp_file.is_file() or not pt_file.is_file():
        raise FileNotFoundError("M2C requires ExpressionData.csv and PseudoTime.csv")
    expression = load_expression(exp_file, gene_axis="rows")
    pt = load_pseudotime(pt_file, expression.index)
    tfs = load_tfs(tf_list, list(expression.columns)) if tf_list else list(expression.columns)
    outer_fit, outer_score, outer_test, split = split_three(expression, seed=seed)
    cfg = CrossfitConfig(folds=folds, seed=seed, max_bins=bins,
                         bootstrap_repeats=bootstrap_repeats)
    features, targets, manifest = make_directional_meta_dataset(
        outer_fit, pt.loc[outer_fit.index], tfs, cfg, dataset_id or raw_dir.name
    )
    manifest.update({
        "outer_split": split,
        "outer_score_cells_excluded": len(outer_score),
        "outer_test_cells_excluded": len(outer_test),
        "expression_sha256": sha256(exp_file),
        "pseudotime_sha256": sha256(pt_file),
        "tf_annotation_source": str(tf_list) if tf_list else "all_measured_genes_not_gold",
        "tf_annotation_sha256": sha256(tf_list) if tf_list else None,
        "feature_file": "m2c_features.csv",
        "target_file": "m2c_targets.csv",
        "training_target": "y_directional_lag_gain_TARGET",
        "orientation_target": "y_pairwise_direction_margin_TARGET",
        "unavailable_target_is_NaN": True,
        "no_gold_reference_file_opened": True,
    })
    out.mkdir(parents=True, exist_ok=True)
    features.to_csv(out / "m2c_features.csv", index=False)
    targets.to_csv(out / "m2c_targets.csv", index=False)
    (out / "m2c_manifest.json").write_text(
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
    m = build(a.raw_dir, a.output_dir, dataset_id=a.dataset_id,
              tf_list=a.tf_list, seed=a.seed, folds=a.folds,
              bins=a.bins, bootstrap_repeats=a.bootstrap_repeats)
    print(json.dumps({k: m[k] for k in [
        "dataset_id", "n_rows", "n_directional_targets_available",
        "n_directional_targets_missing", "uses_gold_reference",
    ]}, indent=2))


if __name__ == "__main__":
    main()
