"""Pseudotime directional evidence for exploratory GRN inference; NOT causal evidence."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .search import universe
from .ridge import GraphRidge


def load_pseudotime(path: str | Path, samples: pd.Index) -> pd.DataFrame:
    """BEELINE: cell IDs as row index, one or more branch pseudotime columns."""
    raw = pd.read_csv(path, index_col=0, na_values=["NA", "NaN", "nan", ""])
    raw.index = raw.index.map(str)
    ids = pd.Index(samples).map(str)
    if raw.index.has_duplicates or ids.has_duplicates:
        raise ValueError("Duplicate cell identifiers")
    missing = ids.difference(raw.index)
    if len(missing):
        raise ValueError(f"Pseudotime missing {len(missing)} expression cells")
    if raw.empty or raw.shape[1] < 1:
        raise ValueError("Empty pseudotime file")
    pt = raw.reindex(ids).apply(pd.to_numeric, errors="coerce")
    if np.isinf(pt.to_numpy(dtype=float)).any():
        raise ValueError("Infinite pseudotime")
    if pt.isna().all(axis=1).any():
        raise ValueError("Expression cells with no valid pseudotime branch")
    return pt


def _lagged_gain(source: np.ndarray, target: np.ndarray, lag: int = 1,
                 alpha: float = 1.0) -> float:
    """Extra in-sample ridge fit for target[b+lag] given source[b],
    target[b] and the pseudotime-bin index. Heuristic, not causal.
    """
    n = len(source) - lag
    if n < 6:
        return 0.
    y = target[lag:].astype(float)
    time = np.linspace(-1., 1., n)
    base = np.column_stack((np.ones(n), target[:-lag], time))
    extended = np.column_stack((base, source[:-lag]))
    def mse(model):
        z = model[:, 1:]
        z0 = z - z.mean(axis=0)
        z0 /= np.maximum(z0.std(axis=0), 1e-10)
        yc = y - y.mean()
        coef = np.linalg.solve(z0.T @ z0 + alpha * np.eye(z0.shape[1]), z0.T @ yc)
        pred = y.mean() + z0 @ coef
        return float(np.mean((y - pred) ** 2))
    before, after = mse(base), mse(extended)
    return float(np.clip((before - after) / max(before, 1e-10), 0., 1.))


def temporal_margin(fit: pd.DataFrame, pt_fit: pd.DataFrame, tfs: list[str],
                    bins: int = 16, min_cells_per_bin: int = 3,
                    lag: int = 1, alpha: float = 1.) -> pd.DataFrame:
    """Signed antisymmetric forward-vs-reverse gain over each valid branch.
    Temporal features are built exclusively on fitting cells.
    """
    if not 8 <= bins <= 100 or lag < 1 or alpha <= 0:
        raise ValueError("Invalid temporal hyperparameters")
    if not fit.index.equals(pt_fit.index):
        raise ValueError("Pseudotime index must align exactly to fitting samples")
    genes = list(fit.columns)
    candidates = universe(genes, tfs)
    weights = np.zeros((len(genes), len(genes)), dtype=float)
    values = np.zeros_like(weights)
    x = fit.to_numpy(dtype=float)
    for branch in pt_fit.columns:
        times = pt_fit[branch].to_numpy(dtype=float)
        indices = np.where(np.isfinite(times))[0]
        if len(indices) < bins * min_cells_per_bin:
            continue
        indices = indices[np.argsort(times[indices], kind="stable")]
        groups = np.array_split(indices, bins)
        if any(len(group) < min_cells_per_bin for group in groups):
            continue
        means = np.stack([x[group].mean(axis=0) for group in groups])
        means = (means - means.mean(axis=0)) / np.maximum(means.std(axis=0), 1e-10)
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                fwd = _lagged_gain(means[:, i], means[:, j], lag, alpha)
                rev = _lagged_gain(means[:, j], means[:, i], lag, alpha)
                values[i, j] += len(indices) * (fwd - rev)
                values[j, i] -= len(indices) * (fwd - rev)
                weights[i, j] += len(indices)
                weights[j, i] += len(indices)
    if not np.any(weights):
        raise ValueError("No pseudotime branch has enough fitting cells")
    index = {g: i for i, g in enumerate(genes)}
    return pd.DataFrame([
        (s, t, float(values[index[s], index[t]] / weights[index[s], index[t]]))
        for s, t in candidates
    ], columns=["source", "target", "score"])


def temporal_greedy(engine: GraphRidge, tfs: list[str], margins: pd.DataFrame,
                    max_edges: int = 20, temporal_weight: float = .015):
    """Graph-conditioned greedy baseline plus an independently computed time margin."""
    if temporal_weight < 0 or max_edges < 0:
        raise ValueError("Invalid temporal weight or edit budget")
    candidates = universe(engine.genes, tfs)
    m = dict(zip(zip(margins.source, margins.target), margins.score))
    if len(m) != len(margins) or set(m) != set(candidates):
        raise ValueError("Temporal scores must cover every candidate exactly once")
    chosen, remaining, history = set(), set(candidates), []
    for step in range(min(max_edges, len(candidates))):
        gain, edge = max(
            ((engine.gain(chosen, e) + temporal_weight * m[e], e) for e in remaining),
            key=lambda v: (v[0], -candidates.index(v[1]))
        )
        if gain <= 0:
            break
        before = engine.graph_loss(chosen)
        chosen.add(edge)
        remaining.remove(edge)
        after = engine.graph_loss(chosen)
        history.append({"step": step + 1, "source": edge[0], "target": edge[1],
                        "combined_gain": gain, "expression_gain": before - after,
                        "temporal_margin": m[edge], "loss_before": before,
                        "loss_after": after})
    leftover = sorted(
        remaining,
        key=lambda e: (-(engine.gain(chosen, e) + temporal_weight * m[e]), e)
    )
    ordered = [(v["source"], v["target"]) for v in history] + leftover
    ranked = pd.DataFrame(ordered, columns=["source", "target"])
    ranked["score"] = np.arange(len(ordered), 0, -1, dtype=float)
    return ranked, pd.DataFrame(history)
