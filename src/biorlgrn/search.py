"""Non-RL reference methods using only expression data and an independent TF list."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .ridge import GraphRidge


def universe(genes: list[str], tfs: list[str]) -> list[tuple[str, str]]:
    return [(s, t) for s in tfs for t in genes if s != t]


def correlation_rank(fit: pd.DataFrame, tfs: list[str], method: str = "pearson") -> pd.DataFrame:
    if method not in ("pearson", "spearman"):
        raise ValueError("Unknown correlation method")
    cor = fit.corr(method=method).abs().fillna(0.)
    pairs = universe(list(fit.columns), tfs)
    out = pd.DataFrame(pairs, columns=["source", "target"])
    out["score"] = [float(cor.loc[s, t]) for s, t in pairs]
    return out.sort_values(["score", "source", "target"],
                           ascending=[False, True, True]).reset_index(drop=True)


def greedy_rank(engine: GraphRidge, tfs: list[str], max_edges: int,
                sparsity_penalty: float = 0.0):
    """Rank edges by sequential marginal graph-utility improvement.

    Utility: -graph_loss - sparsity_penalty * (edge_count / candidate_count).
    This prototype is GREEDY, not an RL policy, and does not use biological priors.
    """
    if max_edges < 0 or sparsity_penalty < 0:
        raise ValueError("max_edges and sparsity_penalty must be nonnegative")
    candidates = universe(engine.genes, tfs)
    chosen: set[tuple[str, str]] = set()
    path = []
    remaining = set(candidates)
    per_edit_cost = sparsity_penalty / max(len(candidates), 1)
    for step in range(min(max_edges, len(candidates))):
        scored = [(engine.gain(chosen, e) - per_edit_cost, e) for e in remaining]
        best_gain, best_edge = max(scored, key=lambda item: (item[0], -candidates.index(item[1])))
        if best_gain <= 0:
            break
        before = engine.graph_loss(chosen)
        chosen.add(best_edge)
        remaining.remove(best_edge)
        after = engine.graph_loss(chosen)
        path.append({"step": step + 1, "source": best_edge[0],
                     "target": best_edge[1], "utility_gain": best_gain,
                     "loss_before": before, "loss_after": after})
    # Never-selected edges are ordered by marginal value in the final graph.
    # This is a heuristic ranking, not a calibrated regulatory probability.
    leftover = sorted(remaining, key=lambda e: (-engine.gain(chosen, e), e))
    ranking = [(r["source"], r["target"]) for r in path] + leftover
    n = len(ranking)
    ranked = pd.DataFrame(ranking, columns=["source", "target"])
    ranked["score"] = np.arange(n, 0, -1, dtype=float)
    return ranked, pd.DataFrame(path, columns=[
        "step", "source", "target", "utility_gain", "loss_before", "loss_after"
    ])
