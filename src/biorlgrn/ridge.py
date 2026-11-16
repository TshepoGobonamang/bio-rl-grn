"""Graph-conditioned ridge loss evaluated on a partition not used for fitting."""
from __future__ import annotations
from functools import lru_cache
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


class GraphRidge:
    """Loss for a TF-target graph G; changes when parent sets change.

    For each target j: ridge fitted on fit samples; squared error measured on
    score samples, scaled by fit variance. Empty parent set predicts fit mean.
    """
    def __init__(self, fit: pd.DataFrame, score: pd.DataFrame, alpha: float = 1.0):
        if alpha <= 0:
            raise ValueError("alpha must be positive")
        if list(fit.columns) != list(score.columns):
            raise ValueError("fit/score genes must match and have same order")
        if len(fit) < 3 or len(score) < 1:
            raise ValueError("Insufficient fit or score samples")
        self.genes = list(fit.columns)
        self.index = {g: i for i, g in enumerate(self.genes)}
        self.F = fit.to_numpy(dtype=float)
        self.S = score.to_numpy(dtype=float)
        self.n_fit = len(fit)
        self.alpha = alpha
        self.fmean = self.F.mean(axis=0)
        self.fvar = np.maximum(self.F.var(axis=0), 1e-8)

    @lru_cache(maxsize=100000)
    def target_loss(self, target: str, parents: tuple[str, ...]) -> float:
        if target in parents:
            raise ValueError("Self-regulation not allowed in prototype")
        j = self.index[target]
        y_fit = self.F[:, j]
        if not parents:
            pred = np.full(len(self.S), self.fmean[j], dtype=float)
        else:
            ix = [self.index[p] for p in parents]
            # sklearn Ridge objective ||y-Xb||^2 + alpha*n_fit ||b||^2
            model = Ridge(alpha=self.alpha * self.n_fit, fit_intercept=True)
            model.fit(self.F[:, ix], y_fit)
            pred = model.predict(self.S[:, ix])
        return float(np.mean((self.S[:, j] - pred) ** 2) / self.fvar[j])

    def graph_loss(self, edges: set[tuple[str, str]]) -> float:
        parents = {g: [] for g in self.genes}
        for src, tgt in edges:
            parents[tgt].append(src)
        return float(sum(self.target_loss(g, tuple(sorted(parents[g])))
                         for g in self.genes) / len(self.genes))

    def gain(self, edges: set[tuple[str, str]], edge: tuple[str, str],
             action: str = "add") -> float:
        """Positive gain means the requested change reduces held-out loss."""
        src, tgt = edge
        selected = {s for s, t in edges if t == tgt}
        if action == "add" and edge in edges:
            raise ValueError("Edge already present")
        if action == "remove" and edge not in edges:
            raise ValueError("Edge not present")
        before = self.target_loss(tgt, tuple(sorted(selected)))
        if action == "add":
            selected.add(src)
        elif action == "remove":
            selected.remove(src)
        else:
            raise ValueError("Unknown action")
        after = self.target_loss(tgt, tuple(sorted(selected)))
        return (before - after) / len(self.genes)
