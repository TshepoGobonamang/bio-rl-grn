"""Evaluation-only overlap between a FROZEN curated prior and hESC ChIP reference.

Not a true independent validation. CollecTRI source literature could share
experimental provenance with the ChIP reference. No source tuning or training
is allowed after these counts are inspected.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

NETWORKS_MD5 = "a7dca84d302bb127bad962ae52815662"


def digest(path: Path) -> str:
    h = hashlib.md5()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_hesc_chip(archive: Path) -> tuple[pd.DataFrame, dict]:
    if digest(archive) != NETWORKS_MD5:
        raise ValueError("BEELINE-Networks checksum mismatch; refusing reference")
    with ZipFile(archive) as zf:
        names = [n for n in zf.namelist()
                 if n.lower().endswith("/human/hesc-chip-seq-network.csv")]
        if len(names) != 1:
            raise ValueError(f"Expected one human hESC ChIP reference, found {names}")
        content = zf.read(names[0])
    # BEELINE input schemas can be named or headerless; detect explicitly.
    first = pd.read_csv(io.BytesIO(content), nrows=3)
    standard = [
        ("source", "target"), ("Gene1", "Gene2"), ("gene1", "gene2"),
        ("TF", "Target"), ("regulator", "target")
    ]
    pair = next((p for p in standard if set(p).issubset(first.columns)), None)
    if pair is not None:
        frame = pd.read_csv(io.BytesIO(content))
        edges = frame[list(pair)].copy()
        edges.columns = ["source", "target"]
    elif first.shape[1] == 2:
        # Two unnamed fields? Read all rows with no header, then discard
        # literal source/target header tokens if present.
        edges = pd.read_csv(io.BytesIO(content), header=None)
        edges.columns = ["source", "target"]
        if str(edges.source.iloc[0]).lower() in {
                "source", "gene1", "tf", "regulator"}:
            edges = edges.iloc[1:]
    else:
        raise ValueError(f"Unsupported reference header; columns {list(first.columns)}")
    edges["source"] = edges.source.astype(str).str.strip()
    edges["target"] = edges.target.astype(str).str.strip()
    if edges[["source", "target"]].isna().any().any() or len(edges) == 0:
        raise ValueError("Reference contains invalid rows")
    unique = edges.drop_duplicates().reset_index(drop=True)
    unique = unique.loc[unique.source.ne(unique.target)].copy()
    return unique, {"archive_md5": NETWORKS_MD5, "archive_member": names[0],
                    "n_unique_nonself_reference_edges": len(unique),
                    "reference_version": "BEELINE-Networks v4 (published 2026)"}


def audit(prior_csv: Path, networks: Path, output: Path):
    prior = pd.read_csv(prior_csv)
    if not {"source", "target", "observed_curated_association"}.issubset(prior):
        raise ValueError("Missing frozen prior edge table")
    if (prior.observed_curated_association != 1).any():
        raise ValueError("Prior must be a curated association indicator")
    if prior.duplicated(["source", "target"]).any():
        raise ValueError("Prior contains duplicate candidate edges")
    gold, provenance = load_hesc_chip(networks)
    a = set(zip(prior.source, prior.target))
    b = set(zip(gold.source, gold.target))
    measured_genes = set(prior.source) | set(prior.target)
    eligible = {(s,t) for s,t in b
                if s in measured_genes and t in measured_genes}
    overlap = a & b
    report = {
        "status": "post_prior_freeze_reference_overlap_diagnostic",
        "prior_edges": len(a),
        "reference_total_edges": len(b),
        "reference_edges_with_endpoints_in_matched_prior_gene_set": len(eligible),
        "overlapping_directed_edges": len(overlap),
        "fraction_prior_edges_in_reference": len(overlap) / max(len(a),1),
        "fraction_reference_eligible_edges_in_prior": len(overlap) / max(len(eligible),1),
        "prior_retained_unchanged": True,
        "gold_used_for_feature_selection": False,
        "source_level_provenance_disjointness_proved": False,
        "interpretation": (
            "Edge-set overlap is descriptive agreement, NOT proof of independent "
            "validation. Underlying ChIP experiments and curated literature may overlap."
        ),
        "reference_provenance": provenance,
    }
    out = Path(output);out.mkdir(parents=True,exist_ok=True)
    (out / "m2d2_hesc_chip_overlap_diagnostic.json").write_text(
        json.dumps(report,indent=2),encoding="utf-8"
    )
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior",type=Path,required=True)
    p.add_argument("--networks-archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(audit(a.prior,a.networks_archive,a.output),indent=2))


if __name__=="__main__":
    main()
