"""M2D2 real hESC + CollecTRI curated-interaction integration SMOKE TEST.

Cohort is prior-enriched by construction and cannot be used to claim unbiased
GRN recovery, external regulatory validation, or independence from ChIP-seq gold.
No gold-standard network is read, and no model selection uses gold reference.
"""
from __future__ import annotations
import argparse
import io
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from biorlgrn.io import sha256
from biorlgrn.m11 import split_three
from biorlgrn.m2_parent import ParentConfig, parent_conditional_dataset
from biorlgrn.m2_external import align_evidence
from biorlgrn.pseudotime import load_pseudotime
from scripts.match_hesc_collecTRI import hesc_gene_inventory


def choose_prior_enriched_subset(links: pd.DataFrame,
                                  tf_count: int = 8,
                                  target_count: int = 24) -> list[str]:
    """Deterministic 8 regulator + 24 target genes; prior-driven intentionally."""
    tf_rank = links.groupby("source").target.nunique().sort_values(ascending=False)
    ranked_tfs = sorted(tf_rank.index, key=lambda s: (-tf_rank[s], s))
    tfs = ranked_tfs[:tf_count]
    if len(tfs) < tf_count:
        raise ValueError("Too few observed human TFs for pilot")
    partial = links.loc[links.source.isin(tfs)]
    target_rank = partial.groupby("target").source.nunique()
    targets = [x for x in sorted(target_rank.index,
                                 key=lambda t: (-target_rank[t], t))
               if x not in tfs][:target_count]
    if len(targets) < target_count:
        raise ValueError("Too few available target genes for matched pilot")
    genes = tfs + targets
    if len(genes) != tf_count + target_count or len(set(genes)) != len(genes):
        raise ValueError("Duplicate exploratory genes")
    return genes


def _extract_subset(archive: Path, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    members = {
        "expression": "BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv",
        "pseudotime": "BEELINE-data/inputs/scRNA-Seq/hESC/PseudoTime.csv",
    }
    frames = []
    wanted = set(genes)
    with ZipFile(archive) as zf:
        for chunk in pd.read_csv(zf.open(members["expression"]),
                                 index_col=0, chunksize=500):
            pick = chunk.loc[chunk.index.astype(str).isin(wanted)].copy()
            if not pick.empty:
                frames.append(pick)
        if not frames:
            raise ValueError("None of the matched human genes exists in expression")
        exp = pd.concat(frames)
        exp.index = exp.index.astype(str)
        if not set(genes).issubset(exp.index) or exp.index.has_duplicates:
            raise ValueError("Cannot uniquely align selected genes to expression")
        exp = exp.loc[genes].astype(float)
        if not np.isfinite(exp.to_numpy()).all():
            raise ValueError("Non-finite expression")
        pt = pd.read_csv(zf.open(members["pseudotime"]), index_col=0)
    x = exp.T
    x.index = x.index.map(str)
    pt.index = pt.index.map(str)
    if x.index.has_duplicates or not x.index.isin(pt.index).all():
        raise ValueError("hESC cells have missing/duplicated pseudotime IDs")
    return x, pt.reindex(x.index)


def run(archive: Path, matched_prior_csv: Path, out: Path) -> dict:
    genes, inventory = hesc_gene_inventory(archive)
    prior = pd.read_csv(matched_prior_csv)
    if not {"source", "target", "observed_curated_association"}.issubset(prior):
        raise ValueError("Expected vetted CollecTRI gene-overlap source")
    if (prior.observed_curated_association != 1).any():
        raise ValueError("Curated observation is presence, not estimated probability")
    cohort = choose_prior_enriched_subset(prior)
    x, pt = _extract_subset(archive, cohort)
    fit, score, test, split = split_three(x, seed=42)
    aligned = prior.loc[prior.source.isin(cohort) &
                        prior.target.isin(cohort) &
                        prior.source.ne(prior.target), ["source", "target"]].copy()
    aligned["evidence_score"] = 1.0
    # A 1 is a binary *recorded interaction*, not calibrated trust.
    cfg = ParentConfig(seed=42, max_parents=2)
    f, t, manifest = parent_conditional_dataset(
        fit, pt.loc[fit.index], list(x.columns), cfg, "hESC_COLLECTRI_EXPLORATORY",
        biological_links=aligned
    )
    coverage = {
        "stage": "real_biological_evidence_feature_integration_smoke_test_only",
        "sample": "BEELINE human embryonic stem cell scRNA-seq (hESC)",
        "prior": "CollecTRI v2.0 curated human TF-target interactions (2023)",
        "n_all_measured_genes": len(genes),
        "n_pilot_genes": len(cohort),
        "gene_selection": "explicitly prior-enriched 8 curated regulators + 24 targets; not unbiased",
        "n_cells": len(x),
        "n_matched_prior_directed_edges": len(aligned),
        "n_candidate_edges": int(manifest["n_candidates_per_fold"]),
        "fraction_of_candidate_edges_with_prior": len(aligned) / manifest["n_candidates_per_fold"],
        "n_meta_records": len(f),
        "n_outer_fit": len(fit), "n_internal_score_excluded": len(score),
        "n_outer_test_excluded": len(test),
        "prior_has_unverified_overlap_with_chip_reference": True,
        "no_benchmark_gold_or_chip_reference_opened": True,
        "no_grn_auroc_or_ap_claim": True,
        "source_info": inventory,
        "dataset_manifest": manifest,
    }
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    f.to_csv(out / "hesc_real_features.csv", index=False)
    t.to_csv(out / "hesc_real_targets.csv", index=False)
    # Gene list is intentionally transparent about prior-driven enrichment.
    pd.DataFrame({"gene": cohort, "cohort_source": (
        ["TF_prior_driven"] * 8 + ["target_prior_driven"] * 24)
    }).to_csv(out / "hesc_cohort_gene_list.csv", index=False)
    (out / "hesc_exploratory_manifest.json").write_text(
        json.dumps(coverage, indent=2), encoding="utf-8"
    )
    return {k: coverage[k] for k in (
        "n_all_measured_genes", "n_cells", "n_pilot_genes",
        "n_matched_prior_directed_edges",
        "n_meta_records", "n_outer_fit",
        "fraction_of_candidate_edges_with_prior"
    )}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--matched-prior", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    print(json.dumps(run(a.archive, a.matched_prior, a.output), indent=2))


if __name__ == "__main__":
    main()
