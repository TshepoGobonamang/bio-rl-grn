"""M2D2 actual source-validation contracts, with small fabricated fixtures only."""
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import pytest
from scripts.match_hesc_collecTRI import map_regulons


def test_mapping_returns_real_gene_overlap_and_no_fake_confidence():
    data = pd.DataFrame([
        {"source": "SOX2", "target": "POU5F1", "weight": 1},
        {"source": "SOX2", "target": "POU5F1", "weight": 1},
        {"source": "NANOG", "target": "POU5F1", "weight": -1},
        {"source": "MYC", "target": "DNMT3B", "weight": 1},
    ])
    matched, audit = map_regulons(
        ["SOX2", "POU5F1", "NANOG", "OCT4", "GAPDH"], data
    )
    assert len(matched) == 2
    assert audit["n_matched_directed_pairs"] == 2
    assert audit["reference_overlap_checked"] is False
    assert audit["confidence_semantics"].startswith("binary")
    assert matched.observed_curated_association.eq(1).all()
    assert matched.loc[
        (matched.source == "SOX2") & (matched.target == "POU5F1"),
        "n_source_rows"
    ].iloc[0] == 2


def test_gene_mapping_not_alias_normalized():
    data = pd.DataFrame({"source": ["OCT4"], "target": ["SOX2"],
                         "weight": [1]})
    out, audit = map_regulons(["POU5F1", "SOX2"], data)
    assert out.empty
    assert audit["n_matched_directed_pairs"] == 0


def test_missing_gene_identifiers_rejected():
    data = pd.DataFrame({"source": ["SOX2", None],
                         "target": ["NANOG", "SOX2"],
                         "weight": [1, 1]})
    with pytest.raises(ValueError, match="missing"):
        map_regulons(["SOX2", "NANOG"], data)


def test_binary_curated_presence_not_motif_probability():
    data = pd.DataFrame({"source": ["SOX2"], "target": ["NANOG"],
                         "weight": [-1]})
    out, meta = map_regulons(["SOX2", "NANOG"], data)
    assert out.signed_regulatory_effect.iloc[0] == -1
    assert out.observed_curated_association.iloc[0] == 1
    assert "NOT calibrated" in meta["confidence_semantics"]
