"""Data import with explicit orientation and test-label separation."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_expression(path: str | Path, gene_axis: str = "rows") -> pd.DataFrame:
    """BEELINE uses genes in rows, cells in columns; internal format is samples x genes."""
    raw = pd.read_csv(path, index_col=0)
    if gene_axis not in {"rows", "columns"}:
        raise ValueError("gene_axis must be 'rows' or 'columns'")
    frame = raw.T if gene_axis == "rows" else raw
    frame.index = frame.index.map(str)
    frame.columns = frame.columns.map(str)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError("Duplicate sample or gene identifiers")
    if frame.empty:
        raise ValueError("Empty expression matrix")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("Expression matrix contains missing or non-finite values")
    frame.index.name = "sample"
    return frame.astype(float)


def load_gold(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    options = [("source", "target"), ("Gene1", "Gene2"), ("gene1", "gene2")]
    for a, b in options:
        if a in df.columns and b in df.columns:
            result = df[[a, b]].rename(columns={a: "source", b: "target"}).copy()
            result["source"] = result["source"].astype(str)
            result["target"] = result["target"].astype(str)
            return result.drop_duplicates().reset_index(drop=True)
    # For headerless two-column reference networks.
    test = pd.read_csv(path, header=None)
    if test.shape[1] in (2, 3):
        test = test.iloc[:, :2]
        test.columns = ["source", "target"]
        return test.astype(str).drop_duplicates().reset_index(drop=True)
    raise ValueError("Gold standard must have Gene1/Gene2 or source/target columns")


def load_tfs(path: str | Path, genes: list[str]) -> list[str]:
    df = pd.read_csv(path)
    col = "tf" if "tf" in df.columns else df.columns[0]
    tfs = list(dict.fromkeys(df[col].dropna().astype(str)))
    unknown = set(tfs) - set(genes)
    if unknown:
        raise ValueError(f"TF annotations contain genes not in expression: {sorted(unknown)[:5]}")
    if not tfs:
        raise ValueError("Empty TF list")
    return tfs


def prepare_dataset(expression_path: str | Path, output_dir: str | Path,
                    ground_truth_path: str | Path | None = None,
                    tf_path: str | Path | None = None,
                    gene_axis: str = "rows") -> dict:
    """Never infers allowed TFs from the gold standard; absent list -> all genes."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    expression_path = Path(expression_path)
    X = load_expression(expression_path, gene_axis)
    X.to_csv(out / "expression.csv")
    tfs = load_tfs(tf_path, list(X.columns)) if tf_path else list(X.columns)
    pd.DataFrame({"tf": tfs}).to_csv(out / "tf_list.csv", index=False)
    manifest = {
        "source_expression": str(expression_path.resolve()),
        "source_expression_sha256": sha256(expression_path),
        "gene_axis_original": gene_axis,
        "shape_samples_by_genes": list(X.shape),
        "tf_source": str(Path(tf_path).resolve()) if tf_path else None,
        "tf_universe_rule": "independent file" if tf_path else "all measured genes (not gold-derived)",
        "evaluation_label_access": "gold_edges.csv is for evaluate command only",
    }
    if tf_path:
        manifest["tf_source_sha256"] = sha256(Path(tf_path))
    if ground_truth_path:
        path = Path(ground_truth_path)
        gold = load_gold(path)
        keep = set(X.columns)
        gold = gold[gold.source.isin(keep) & gold.target.isin(keep)]
        gold = gold[(gold.source != gold.target) & gold.source.isin(set(tfs))]
        gold.to_csv(out / "gold_edges.csv", index=False)
        manifest["gold_source"] = str(path.resolve())
        manifest["gold_source_sha256"] = sha256(path)
        manifest["gold_edges_in_candidate_universe"] = len(gold)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def read_prepared(data_dir: str | Path) -> tuple[pd.DataFrame, list[str]]:
    p = Path(data_dir)
    X = load_expression(p / "expression.csv", gene_axis="columns")
    tfs = load_tfs(p / "tf_list.csv", list(X.columns))
    return X, tfs
