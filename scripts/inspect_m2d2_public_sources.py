"""Inspect pinned public human single-cell expression and CollecTRI source *schemas*.

No benchmark evaluation references are opened. No TF-target edges inferred or
model fit in this stage. The output is safe, aggregate source metadata.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def inspect(archive: Path, regulons: Path, source: Path, output: Path) -> dict:
    if md5(archive) != "9ea93065004f2f12022adeeaf3728edd":
        raise ValueError("Unexpected BEELINE v4 archive version")
    if md5(regulons) != "cee4a3943c059e6dd8796ced6dec44f6":
        raise ValueError("Unexpected CollecTRI v2 regulon version")
    if md5(source) != "c8622a46a49e21e955182f000cf1c067":
        raise ValueError("Unexpected CollecTRI v2 source version")
    with ZipFile(archive) as zf:
        names = zf.namelist()
        selected = {}
        for part in ("ExpressionData.csv", "PseudoTime.csv", "GeneOrdering.csv"):
            matches = [n for n in names if "/scRNA-Seq/hESC/" in n
                       and n.endswith("/" + part) and not n.startswith("__MACOSX")]
            if len(matches) != 1:
                raise ValueError(f"Expected one hESC {part}; found {matches}")
            name = matches[0]
            info = zf.getinfo(name)
            with zf.open(name) as handle:
                frame = pd.read_csv(handle, nrows=5, index_col=0)
            selected[part] = {
                "zip_member": name, "uncompressed_bytes": info.file_size,
                "first_5_index": list(map(str, frame.index)),
                "first_12_columns": list(map(str, frame.columns[:12])),
                "n_header_columns": int(frame.shape[1]+1),
            }
        other_ref = sum(("BEELINE-Networks" in n) or ("GroundTruthNetwork" in n)
                        for n in names)
    reg_head = pd.read_csv(regulons, nrows=5)
    src_head = pd.read_csv(source, sep="\t", nrows=5, low_memory=False)
    result = {
        "status": "source_schema_inspection_only",
        "hesc": selected,
        "collectri_regulons": {
            "sha256": sha256(regulons),
            "columns": list(map(str, reg_head.columns)),
            "first_5_rows": reg_head.where(pd.notna(reg_head), None).to_dict(orient="records"),
            "n_header_columns": len(reg_head.columns),
        },
        "collectri_source": {
            "sha256": sha256(source),
            "columns": list(map(str, src_head.columns))[:45],
            "n_header_columns": len(src_head.columns),
        },
        "archive_reference_file_count_visible_but_not_opened": other_ref,
        "expression_sha256": sha256(archive),
        "zero_gold_references_opened": True,
        "note": "CollecTRI curated regulatory associations may overlap ChIP-seq evaluation sources. No independence established yet.",
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "m2d2_source_inventory.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--regulons", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    meta = inspect(a.archive, a.regulons, a.source, a.output)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
