"""Tiny synthetic SMOKE TEST. NOT a BEELINE dataset and NOT a publication result."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def make_demo(output: str | Path, seed: int = 42):
    p = Path(output)
    p.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    n = 100
    g1, g2 = rng.normal(size=(2, n))
    g3 = 0.9 * g1 + rng.normal(scale=0.2, size=n)
    g4 = -0.8 * g2 + rng.normal(scale=0.2, size=n)
    g5 = 0.5 * g3 + rng.normal(scale=0.3, size=n)
    X = pd.DataFrame({"G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": g5})
    X.index = [f"cell-{i}" for i in range(n)]
    X.index.name = "sample"
    X.to_csv(p / "expression.csv")
    pd.DataFrame({"tf": ["G1", "G2", "G3"]}).to_csv(p / "tf_list.csv", index=False)
    pd.DataFrame({"source": ["G1", "G2", "G3"],
                  "target": ["G3", "G4", "G5"]}).to_csv(p / "gold_edges.csv", index=False)
    (p / "manifest.json").write_text(json.dumps({
        "kind": "toy_smoke_test_not_beeline", "seed": seed,
        "warning": "Use only to check that the pipeline executes"}, indent=2))
    return p
