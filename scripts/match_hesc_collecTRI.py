"""Map a pinned CollecTRI v2.0 human TF-target regulon to real BEELINE hESC genes.

This is an independent-of-EXPRESSION prior coverage audit. It does NOT verify
independence from BEELINE's ChIP-seq reference, or constitute biological
validation of predicted interactions. No evaluation reference is opened.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

BEELINE_MD5 = "9ea93065004f2f12022adeeaf3728edd"
COLLECTRI_MD5 = "cee4a3943c059e6dd8796ced6dec44f6"
SOURCE_MD5 = "c8622a46a49e21e955182f000cf1c067"
REGULON_URL = "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"


def digest(path: Path, algorithm="sha256") -> str:
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def hesc_gene_inventory(archive: Path) -> tuple[list[str], dict]:
    if digest(archive, "md5") != BEELINE_MD5:
        raise ValueError("BEELINE archive checksum does not match pinned release")
    with ZipFile(archive) as zf:
        names = [n for n in zf.namelist() if
                 "/scRNA-Seq/hESC/" in n and n.endswith("/ExpressionData.csv")]
        if len(names) != 1:
            raise ValueError(f"Expected one human hESC expression file, got {names}")
        genes = []
        with zf.open(names[0]) as binary:
            text = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
            reader = csv.reader(text)
            header = next(reader)
            if len(header) < 3:
                raise ValueError("Expression matrix lacks cells")
            for row in reader:
                if row:
                    genes.append(str(row[0]).strip())
        if not genes or len(set(genes)) != len(genes) or any(not g for g in genes):
            raise ValueError("Ambiguous or duplicated hESC gene IDs; no silent deduplication")
        return genes, {
            "zip_member": names[0], "n_expression_genes": len(genes),
            "n_expression_cells": len(header)-1,
            "gene_id_example": genes[:12],
            "zip_uncompressed_expression_bytes": zf.getinfo(names[0]).file_size,
        }


def map_regulons(genes: list[str], regulons: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    # Fail closed if source schema unexpectedly changes.
    required = {"source", "target"}
    if not required.issubset(regulons.columns):
        raise ValueError(f"Unrecognized CollecTRI schema {list(regulons.columns)}")
    rows = regulons.copy()
    if rows[["source", "target"]].isna().any().any():
        raise ValueError("CollecTRI contains missing regulatory gene identifiers")
    rows["source"] = rows["source"].astype(str).str.strip()
    rows["target"] = rows["target"].astype(str).str.strip()
    rows = rows.loc[rows.source.ne("") & rows.target.ne("") & rows.source.ne(rows.target)]
    measured = set(genes)
    matched = rows.loc[rows.source.isin(measured) & rows.target.isin(measured)].copy()
    grouped = matched.groupby(["source", "target"], as_index=False).size().rename(
        columns={"size": "n_source_rows"})
    # One is an *indicator of a curated association*, not probability or binding strength.
    grouped["observed_curated_association"] = 1
    grouped["prior_source"] = "CollecTRI_2023_v2"
    if "weight" in matched.columns:
        w = pd.to_numeric(matched["weight"], errors="raise")
        if not w.dropna().isin([-1, 1]).all():
            raise ValueError("Unexpected CollecTRI regulation signs")
        tmp = matched[["source", "target"]].copy()
        tmp["sign"] = np.sign(w).astype(float)
        signs = tmp.groupby(["source", "target"]).sign.agg(
            lambda x: int(x.iloc[0]) if x.nunique() == 1 else 0
        ).reset_index(name="signed_regulatory_effect")
        grouped = grouped.merge(signs, on=["source", "target"], validate="one_to_one")
    else:
        grouped["signed_regulatory_effect"] = np.nan
    coverage = {
        "n_collecTRI_raw_rows": len(regulons),
        "n_collecTRI_gene_pairs": len(
            rows[["source", "target"]].drop_duplicates()
        ),
        "n_hesc_measured_genes": len(genes),
        "n_hesc_measured_regulators_in_collectri": int(
            len(set(matched.source))
        ),
        "n_hesc_measured_targets_in_collectri": int(
            len(set(matched.target))
        ),
        "n_matched_directed_pairs": len(grouped),
        "n_matched_source_rows": len(matched),
        "gene_match_method": "exact, case-sensitive human gene identifiers; no alias or orthology mapping",
        "confidence_semantics": "binary association presence, NOT calibrated binding probability",
        "reference_overlap_checked": False,
        "independent_from_expression_measurement": True,
        "independent_from_BEEELINE_evaluation_reference": "unverified",
        "evaluation_reference_accessed": False,
    }
    return grouped, coverage


def execute(archive: Path, regulons_file: Path, source_file: Path, output: Path):
    for path, md5 in [(regulons_file, COLLECTRI_MD5),
                      (source_file, SOURCE_MD5)]:
        if digest(path, "md5") != md5:
            raise ValueError(f"Unexpected source release checksum: {path}")
    genes, expression = hesc_gene_inventory(archive)
    regulons = pd.read_csv(regulons_file)
    matched, summary = map_regulons(genes, regulons)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    matched.to_csv(out / "hesc_collecTRI_matched_prior_edges.csv", index=False)
    report = {
        "status": "real-human-expression-and-curated-prior-gene-overlap_verified_only",
        "BEELINE": expression,
        "COLLECTRI": summary,
        "archive_beeline_sha256": digest(archive),
        "collectri_v2_regulons_sha256": digest(regulons_file),
        "collectri_v2_source_provenance_sha256": digest(source_file),
        "regulons_source_url": REGULON_URL,
        "regulons_release": "CollecTRI 2.0, Zenodo 8192729, July 2023",
        "prior_assay": "multi-source curated human TF-target interactions",
        "warning": (
            "Matching human gene names supports a biological evidence channel, "
            "NOT nonoverlap with ChIP-seq references. No model trained; "
            "no independent regulatory recovery claim. Absent rows mean unknown."
        ),
    }
    (out / "m2d2_hesc_prior_match_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--regulons", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(execute(a.archive, a.regulons, a.source, a.output),
                     indent=2))


if __name__ == "__main__":
    main()
