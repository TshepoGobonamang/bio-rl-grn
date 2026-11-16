"""M1.1: split safeguards, orientation diagnostics and evaluation-only gold audit."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, auc
from .search import universe
from .ridge import GraphRidge


def split_three(X: pd.DataFrame, seed: int = 42, test_fraction: float = 0.2,
                dev_score_fraction: float = 0.25):
    """Fit / iterative graph scoring / untouched expression test, in that order."""
    if not 0 < test_fraction < 0.5 or not 0 < dev_score_fraction < 0.5:
        raise ValueError("Fractions must be in (0, 0.5)")
    ids = np.arange(len(X))
    dev, test = train_test_split(ids, test_size=test_fraction, random_state=seed)
    fit, score = train_test_split(dev, test_size=dev_score_fraction, random_state=seed + 1)
    if min(len(fit), len(score), len(test)) < 3:
        raise ValueError("Every split must contain at least three samples")
    return X.iloc[fit], X.iloc[score], X.iloc[test], {
        "fit_samples": len(fit), "internal_score_samples": len(score),
        "untouched_test_samples": len(test), "seed": seed,
        "test_fraction": test_fraction, "dev_score_fraction": dev_score_fraction,
        "split_indices_disjoint": len(set(fit) | set(score) | set(test)) == len(X),
        "test_evaluated_after_graph_selection": True,
    }


def audit_reference(gold_path: str | Path, genes: list[str], tfs: list[str]):
    """Audit gold rows only after inference; never create candidates from gold."""
    raw = pd.read_csv(gold_path)
    col_pairs = [("Gene1", "Gene2"), ("source", "target"), ("gene1", "gene2")]
    found = next(((a, b) for a, b in col_pairs if a in raw and b in raw), None)
    if found is None:
        raw = pd.read_csv(gold_path, header=None)
        if raw.shape[1] < 2:
            raise ValueError("Reference requires at least two columns")
        edges = raw.iloc[:, :2].copy()
    else:
        edges = raw.loc[:, list(found)].copy()
    edges.columns = ["source", "target"]
    genes, tfs = set(genes), set(tfs)
    status, seen = [], set()
    for s, t in edges.itertuples(index=False, name=None):
        if pd.isna(s) or pd.isna(t):
            label = "missing_endpoint"
        else:
            s, t = str(s), str(t)
            if s not in genes:
                label = "source_not_measured"
            elif t not in genes:
                label = "target_not_measured"
            elif s not in tfs:
                label = "source_not_candidate_tf"
            elif s == t:
                label = "self_loop_excluded"
            elif (s, t) in seen:
                label = "duplicate"
            else:
                label = "eligible"
                seen.add((s, t))
        status.append(label)
    edges = edges.assign(status=status)
    filtered = edges.loc[edges.status == "eligible", ["source", "target"]].copy()
    filtered["source"] = filtered["source"].astype(str)
    filtered["target"] = filtered["target"].astype(str)
    summary = {
        "raw_gold_rows": len(edges), "gold_eligible": len(filtered),
        "gold_exclusions": edges.status.value_counts().to_dict(),
    }
    return filtered.reset_index(drop=True), edges, summary


def orientation_metrics(ranked: pd.DataFrame, true_edges: set[tuple[str, str]], k: int):
    """Direction correctness is conditional on single-direction gold pairs."""
    ordered = ranked.sort_values(["score", "source", "target"], ascending=[False, True, True])
    selected = set(zip(ordered.head(k).source, ordered.head(k).target))
    reciprocal = sorted({
        tuple(sorted((s, t))) for s, t in selected if (t, s) in selected
    })
    scores = dict(zip(zip(ranked.source, ranked.target), ranked.score.astype(float)))
    correct = wrong = ties = unavailable = 0
    for s, t in sorted(true_edges):
        if (t, s) in true_edges:
            continue
        if (s, t) not in scores or (t, s) not in scores:
            unavailable += 1
            continue
        delta = scores[(s, t)] - scores[(t, s)]
        if abs(delta) < 1e-12:
            ties += 1
        elif delta > 0:
            correct += 1
        else:
            wrong += 1
    n = correct + wrong + ties
    return {
        "top_k_reciprocal_pairs": len(reciprocal),
        "top_k_edges_in_reciprocal_pairs": 2 * len(reciprocal),
        "reciprocal_pairs": [f"{a}<->{b}" for a, b in reciprocal],
        "direction_reference_unique_pairs": n,
        "direction_correct": correct, "direction_wrong": wrong,
        "direction_ties": ties, "direction_unavailable": unavailable,
        "direction_accuracy_ties_half": (correct + .5 * ties) / n if n else None,
    }


def evaluate_ranked(ranked: pd.DataFrame, true_edges: set[tuple[str, str]],
                    genes: list[str], tfs: list[str], k: int):
    candidates = universe(genes, tfs)
    scores = dict(zip(zip(ranked.source, ranked.target), ranked.score.astype(float)))
    if len(scores) != len(ranked) or set(scores) != set(candidates):
        raise ValueError("Ranking must contain every admissible edge exactly once")
    y = np.array([int(edge in true_edges) for edge in candidates], dtype=int)
    if not 0 < y.sum() < len(y):
        raise ValueError("Reference requires positive and negative candidates")
    pred = np.array([scores[e] for e in candidates], dtype=float)
    if not np.isfinite(pred).all():
        raise ValueError("Non-finite ranking score")
    precision, recall, _ = precision_recall_curve(y, pred)
    top = ranked.sort_values(
        ["score", "source", "target"], ascending=[False, True, True]
    ).head(k)
    top_tp = sum((s, t) in true_edges for s, t in zip(top.source, top.target))
    return {
        "auroc": float(roc_auc_score(y, pred)),
        "average_precision": float(average_precision_score(y, pred)),
        "trapezoidal_pr_auc": float(auc(recall, precision)),
        "positive_prevalence": float(y.mean()),
        "precision_at_k": top_tp / len(top),
        "recall_at_k": top_tp / int(y.sum()),
        "tp_at_k": top_tp, "n_candidates": len(candidates),
        "n_gold": int(y.sum()), **orientation_metrics(ranked, true_edges, k),
    }


def greedy_reciprocity_ablation(engine: GraphRidge, tfs: list[str], max_edges: int,
                                 reverse_penalty: float = .01):
    """Experimental anti-reciprocity heuristic; does NOT establish causal direction."""
    if reverse_penalty < 0:
        raise ValueError("reverse_penalty cannot be negative")
    candidates = universe(engine.genes, tfs)
    chosen, remaining = set(), set(candidates)
    path = []
    for step in range(min(max_edges, len(candidates))):
        gains = [
            (engine.gain(chosen, e) - reverse_penalty * int((e[1], e[0]) in chosen), e)
            for e in remaining
        ]
        value, edge = max(gains, key=lambda x: (x[0], -candidates.index(x[1])))
        if value <= 0:
            break
        before = engine.graph_loss(chosen)
        chosen.add(edge)
        remaining.remove(edge)
        path.append({
            "step": step + 1, "source": edge[0], "target": edge[1],
            "utility_gain_with_penalty": value, "loss_before": before,
            "loss_after": engine.graph_loss(chosen),
        })
    leftover = sorted(remaining, key=lambda e: (
        -(engine.gain(chosen, e) - reverse_penalty * int((e[1], e[0]) in chosen)), e
    ))
    order = [(r["source"], r["target"]) for r in path] + leftover
    ranked = pd.DataFrame(order, columns=["source", "target"])
    ranked["score"] = np.arange(len(order), 0, -1, dtype=float)
    return ranked, pd.DataFrame(path)


def independent_expression_loss(fit: pd.DataFrame, test: pd.DataFrame,
                                chosen: set[tuple[str, str]], alpha: float):
    """Instantiate test engine only AFTER graph selection."""
    engine = GraphRidge(fit, test, alpha=alpha)
    before = engine.graph_loss(set())
    after = engine.graph_loss(chosen)
    return {
        "untouched_null_graph_loss": before,
        "untouched_selected_graph_loss": after,
        "untouched_test_loss_change": before - after,
    }
