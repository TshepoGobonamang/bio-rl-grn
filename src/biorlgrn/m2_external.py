"""Strict import of independently derived, organism-matched edge evidence.

No TF-target priors can be inferred from synthetic BEELINE gene names.
Binding motifs are not regulatory labels; use scores only as weak priors.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .io import sha256

COLUMNS = ["source", "target", "score", "evidence_id"]
FORBIDDEN_FILENAMES = {
    "groundtruthnetwork.csv", "gold_edges.csv", "reference_audit.csv",
    "evaluation_labels.csv", "m2b_grn_metrics.csv"
}
EVIDENCE_TYPES = {"motif_promoter_link", "chip_binding", "atac_motif_link",
                  "perturbation_effect", "curated_binding"}
SPECIES = {"human", "mouse", "rat", "zebrafish", "fruit_fly", "yeast"}
SYNTHETIC_IDS = __import__("re").compile(r"^(?:g|G)\d+$")


def read_external_evidence(path: str | Path, metadata_path: str | Path,
                           measured_genes: list[str], *,
                           dataset_kind: str) -> tuple[pd.DataFrame, dict]:
    """Load independent source provenance WITHOUT touching evaluation gold.

    Missing row denotes unknown evidence, NOT a known non-interaction.
    A motif matrix on its own is NOT a TF-target edge: an explicit genome/promoter
    linking procedure and provenance are required for motif-based edge files.
    """
    p, mpath = Path(path), Path(metadata_path)
    if p.name.lower() in FORBIDDEN_FILENAMES or mpath.name.lower() in FORBIDDEN_FILENAMES:
        raise ValueError("Reference/evaluation files cannot be used as evidence")
    if dataset_kind != "biological":
        raise ValueError("Independent biological priors cannot be assigned to synthetic BEELINE genes")
    if not p.is_file() or not mpath.is_file():
        raise FileNotFoundError("Evidence data and provenance manifest required")
    meta = json.loads(mpath.read_text(encoding="utf-8"))
    required = {"database", "release", "organism", "genome_assembly", "source_url",
                "assay", "derivation", "role", "independent_of_evaluation_reference"}
    if not required.issubset(meta) or any(not str(meta[k]).strip() for k in required):
        raise ValueError("Incomplete independent-evidence provenance")
    if meta["role"] != "predictor_prior_only" or meta["independent_of_evaluation_reference"] is not True:
        raise ValueError("Source cannot be certified as independent of evaluation references")
    if meta["organism"] not in SPECIES or meta["assay"] not in EVIDENCE_TYPES:
        raise ValueError("Unsupported organism or binding evidence type")
    if not str(meta["source_url"]).startswith(("https://", "http://")):
        raise ValueError("Evidence source URL required")
    if "GroundTruthNetwork" in json.dumps(meta) or "BEELINE_GOLD" in json.dumps(meta):
        raise ValueError("Provenance indicates evaluation reference reuse")
    genes = set(measured_genes)
    if genes and all(SYNTHETIC_IDS.fullmatch(str(g)) for g in genes):
        raise ValueError("Gene identities are synthetic, cannot join to biological evidence")
    rows = pd.read_csv(p)
    if not set(COLUMNS).issubset(rows) or rows.empty:
        raise ValueError("Require source, target, score, evidence_id")
    rows = rows[COLUMNS].copy()
    if rows[["source", "target", "evidence_id"]].isna().any().any():
        raise ValueError("Unresolved or missing gene/evidence IDs")
    if rows.duplicated(["source", "target", "evidence_id"]).any():
        raise ValueError("Duplicate source-target-evidence records")
    rows["score"] = pd.to_numeric(rows["score"], errors="raise")
    if not np.isfinite(rows.score).all() or not rows.score.between(0, 1).all():
        raise ValueError("Biological evidence confidence must be finite in [0,1]")
    rows["source"] = rows.source.astype(str)
    rows["target"] = rows.target.astype(str)
    # For retained records: do not silently treat unmatched genes as zero.
    retained = rows.loc[rows.source.isin(genes) & rows.target.isin(genes)
                        & rows.source.ne(rows.target)].copy()
    if retained.empty:
        raise ValueError("No evidence links map to independently measured gene identifiers")
    grouped = retained.groupby(["source", "target"], as_index=False).agg(
        evidence_score=("score", "max"), n_evidence_records=("evidence_id", "nunique")
    )
    audit = {
        "dataset_kind": dataset_kind, "source_file_sha256": sha256(p),
        "metadata_sha256": sha256(mpath), "source_name": meta["database"],
        "source_release": meta["release"], "organism": meta["organism"],
        "genome_assembly": meta["genome_assembly"], "assay": meta["assay"],
        "source_url": meta["source_url"],
        "n_input_records": len(rows), "n_measured_links": len(grouped),
        "n_unmatched_or_self_loop_records": len(rows) - len(retained),
        "independence_asserted_by_provider": True,
        "reference_overlap_checked": False,
        "reference_labels_opened_during_loading": False,
        "warning": "Provider's independence assertion does not prove no empirical overlap."
    }
    return grouped, audit


def align_evidence(candidates: pd.DataFrame, links: pd.DataFrame | None) -> pd.DataFrame:
    """Explicit missingness; no source record is unknown, not a negative edge."""
    frame = candidates[["source", "target"]].copy()
    if links is None:
        frame["m_independent_prior"] = 0
        frame["independent_prior_score"] = np.nan
        return frame
    if links.duplicated(["source", "target"]).any():
        raise ValueError("Duplicated mapped biological evidence")
    frame = frame.merge(links[["source", "target", "evidence_score"]],
                        on=["source", "target"], how="left",
                        validate="one_to_one", sort=False)
    frame["m_independent_prior"] = frame.evidence_score.notna().astype(int)
    frame = frame.rename(columns={"evidence_score": "independent_prior_score"})
    return frame


def audit_evaluation_overlap(links: pd.DataFrame, gold: pd.DataFrame):
    """Post-inference diagnostic only; cannot feed back to construct priors."""
    for f in (links, gold):
        if not {"source", "target"}.issubset(f):
            raise ValueError("Missing directed gene identifiers")
    A = set(zip(links.source, links.target))
    B = set(zip(gold.source, gold.target))
    return {"prior_directed_edges": len(A), "reference_directed_edges": len(B),
            "overlapping_directed_edges": len(A & B),
            "overlap_is_evidence_leakage": "undetermined: provenance review required"}
