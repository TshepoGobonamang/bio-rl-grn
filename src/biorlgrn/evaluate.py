"""Evaluation performed ONLY after predictions have been saved."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from .io import load_gold, read_prepared
from .search import universe


def score_predictions(data_dir: str | Path, ranked_path: str | Path, k: int = 20) -> dict:
    X, tfs = read_prepared(data_dir)
    pairs = universe(list(X.columns), tfs)
    pred = pd.read_csv(ranked_path)
    if not {"source", "target", "score"}.issubset(pred.columns):
        raise ValueError("Prediction schema must contain source,target,score")
    got = list(zip(pred.source.astype(str), pred.target.astype(str)))
    if len(got) != len(set(got)) or set(got) != set(pairs):
        raise ValueError("Prediction must cover each candidate edge exactly once")
    if not np.isfinite(pred.score.astype(float).to_numpy()).all():
        raise ValueError("Non-finite predictions")
    gold_path = Path(data_dir) / "gold_edges.csv"
    if not gold_path.exists():
        raise FileNotFoundError("Gold standard not available; run has still produced edge rankings")
    gold = load_gold(gold_path)
    true = set(zip(gold.source.astype(str), gold.target.astype(str)))
    y = np.array([int(e in true) for e in pairs])
    if not 0 < y.sum() < len(y):
        raise ValueError("Need at least one gold edge and one non-gold candidate")
    m = dict(zip(got, pred.score.astype(float)))
    scores = np.array([m[e] for e in pairs], dtype=float)
    top = pred.sort_values("score", ascending=False).head(k)
    n_correct = sum((s, t) in true for s, t in zip(top.source, top.target))
    return {
        "auroc": float(roc_auc_score(y, scores)),
        "auprc": float(average_precision_score(y, scores)),
        "random_auprc_baseline": float(y.mean()),
        "precision_at_k": n_correct / max(len(top), 1),
        "recall_at_k": n_correct / int(y.sum()),
        "k": k, "n_candidates": len(pairs), "n_gold": int(y.sum()),
        "evaluation_assumption": "closed-world simulated gold standard",
    }
