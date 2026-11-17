"""A real human/prior integration pilot must disclose its selection bias."""
import numpy as np
import pandas as pd
import pytest
from scripts.run_m2d2_hesc_pilot import choose_prior_enriched_subset


def test_prior_enriched_gene_cohort_is_deterministic_and_transparent():
    tfs = [f"TF{i}" for i in range(10)]
    targets = [f"TG{i}" for i in range(30)]
    edges = pd.DataFrame([
        (s,t) for s in tfs for t in targets[:25]
    ], columns=["source","target"])
    a = choose_prior_enriched_subset(edges)
    b = choose_prior_enriched_subset(edges)
    assert a == b and len(a) == 32 and len(set(a)) == 32
    assert set(a[:8]).issubset(tfs)
    assert set(a[8:]).issubset(targets)


def test_not_enough_source_coverage_fails_not_silently_filled():
    p = pd.DataFrame({"source": ["SOX2", "POU5F1"],
                      "target": ["NANOG", "GATA6"]})
    with pytest.raises(ValueError, match="Too few"):
        choose_prior_enriched_subset(p)
