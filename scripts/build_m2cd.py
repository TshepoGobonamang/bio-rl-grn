"""Export gold-blind parent-conditional targets with optional independent biology."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from biorlgrn.io import load_expression, sha256
from biorlgrn.m11 import split_three
from biorlgrn.pseudotime import load_pseudotime
from biorlgrn.m2_external import read_external_evidence
from biorlgrn.m2_parent import ParentConfig, parent_conditional_dataset


def build(raw_dir: Path, output: Path, *, dataset_id: str | None = None,
          dataset_kind: str = "synthetic", evidence: Path | None = None,
          evidence_manifest: Path | None = None, max_parents: int = 2,
          seed: int = 42) -> dict:
    raw, out = Path(raw_dir), Path(output)
    x = load_expression(raw / "ExpressionData.csv")
    p = raw / "PseudoTime.csv"
    pseudotime = load_pseudotime(p, x.index) if p.exists() else None
    links, prior_audit = None, None
    if (evidence is None) != (evidence_manifest is None):
        raise ValueError("Both biological edge scores and provenance are required")
    if evidence is not None:
        links, prior_audit = read_external_evidence(
            evidence, evidence_manifest, list(x.columns), dataset_kind=dataset_kind
        )
    outer_fit, outer_score, outer_test, split = split_three(x, seed=seed)
    features, targets, manifest = parent_conditional_dataset(
        outer_fit, pseudotime.loc[outer_fit.index] if pseudotime is not None else None,
        list(x.columns), ParentConfig(seed=seed, max_parents=max_parents),
        dataset_id or raw.name, biological_links=links
    )
    manifest.update({
        "dataset_kind": dataset_kind,
        "expression_sha256": sha256(raw / "ExpressionData.csv"),
        "pseudotime_sha256": sha256(p) if p.exists() else None,
        "outer_split": split,
        "outer_score_cells_excluded": len(outer_score),
        "outer_test_cells_excluded": len(outer_test),
        "external_prior_provenance": prior_audit,
        "no_gold_file_opened": True,
        "no_independent_biology_claim_for_synthetic_network": links is None,
    })
    out.mkdir(parents=True, exist_ok=True)
    features.to_csv(out / "m2cd_features.csv", index=False)
    targets.to_csv(out / "m2cd_targets.csv", index=False)
    (out / "m2cd_manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--dataset-id")
    p.add_argument("--dataset-kind", choices=["synthetic", "biological"], default="synthetic")
    p.add_argument("--evidence", type=Path)
    p.add_argument("--evidence-manifest", type=Path)
    p.add_argument("--max-parents", type=int, default=2)
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()
    m = build(a.raw_dir, a.output_dir, dataset_id=a.dataset_id,
              dataset_kind=a.dataset_kind, evidence=a.evidence,
              evidence_manifest=a.evidence_manifest,
              max_parents=a.max_parents, seed=a.seed)
    print(json.dumps({k: m[k] for k in (
        "dataset_id", "n_rows", "n_independent_links",
        "independent_biological_prior_present", "no_gold_reference_access"
    )}, indent=2))


if __name__ == "__main__":
    main()
