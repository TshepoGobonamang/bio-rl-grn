"""M2C-D tests: source-blind parent selection, no gold labels, prior provenance."""
import json
import numpy as np
import pandas as pd
import pytest
from biorlgrn.m2_parent import ParentConfig, _scaffold, parent_conditional_dataset
from biorlgrn.m2_external import (
    read_external_evidence, align_evidence, audit_evaluation_overlap
)
from scripts.build_m2cd import build
from scripts.audit_m2cd import audit


def _toy(n=180):
    rng = np.random.default_rng(12)
    latent = rng.normal(size=n)
    x = pd.DataFrame({
        "A": latent + rng.normal(scale=.2, size=n),
        "B": latent + rng.normal(scale=.2, size=n),
        "C": .8 * latent + rng.normal(scale=.2, size=n),
        "D": rng.normal(size=n),
    }, index=[f"cell_{i}" for i in range(n)])
    p = pd.DataFrame({"PseudoTime1": np.linspace(0, 1, n)}, index=x.index)
    return x, p


def test_source_excluded_from_scaffold_and_fit_only_selection():
    x, _ = _toy()
    assert "A" not in _scaffold(x, "A", "C", 2)
    assert "C" not in _scaffold(x, "A", "C", 2)
    assert _scaffold(x, "A", "C", 0) == ()
    assert _scaffold(x, "A", "C", 2) == ("B", "D")


def test_conditional_gain_differs_from_marginal_on_confounding_toy():
    x, p = _toy()
    f, t, m = parent_conditional_dataset(
        x, p, list(x), ParentConfig(max_parents=1), "test")
    selected = f.loc[(f.source == "A") & (f.target == "C")]
    assert (selected.scaffold_n == 1).all()
    j = t.loc[(t.source == "A") & (t.target == "C")]
    assert (j.y_parent_conditional_gain_TARGET -
            j.y_marginal_expression_gain_TARGET).abs().max() > 1e-3
    assert len(f) == len(t) == 36
    assert m["no_gold_reference_access"]
    assert not f.independent_prior_score.notna().any()
    assert not f.m_independent_prior.any()


def test_validation_perturbation_does_not_change_fold_features():
    x, p = _toy()
    f1, t1, _ = parent_conditional_dataset(
        x, p, list(x), ParentConfig(), "test")
    from sklearn.model_selection import KFold
    _, val = next(KFold(3, shuffle=True, random_state=42).split(x))
    changed = x.copy()
    changed.iloc[val, 0] += 12
    f2, t2, _ = parent_conditional_dataset(
        changed, p, list(x), ParentConfig(), "test")
    pd.testing.assert_frame_equal(
        f1[f1.fold_id.eq(0)].reset_index(drop=True),
        f2[f2.fold_id.eq(0)].reset_index(drop=True))
    assert not t1[t1.fold_id.eq(0)].equals(t2[t2.fold_id.eq(0)])


def test_real_evidence_requires_external_provenance(tmp_path):
    evidence = tmp_path / "motifs.csv"
    pd.DataFrame({
        "source": ["STAT3", "STAT3", "JUN"],
        "target": ["SOCS3", "SOCS3", "FOS"],
        "score": [.75, .8, .4],
        "evidence_id": ["site_1", "site_2", "site_3"],
    }).to_csv(evidence, index=False)
    metadata = tmp_path / "source.json"
    meta = {
        "database": "test_fixture_not_real_binding_evidence",
        "release": "test-v1", "organism": "human",
        "genome_assembly": "GRCh38",
        "source_url": "https://jaspar.elixir.no/downloads/",
        "assay": "motif_promoter_link",
        "derivation": "TEST ONLY: hypothetical matched promoter scan",
        "role": "predictor_prior_only", "independent_of_evaluation_reference": True,
    }
    metadata.write_text(json.dumps(meta))
    links, prov = read_external_evidence(
        evidence, metadata, ["STAT3", "SOCS3", "JUN", "FOS"],
        dataset_kind="biological")
    assert len(links) == 2
    assert links.loc[links.target.eq("SOCS3"), "evidence_score"].iloc[0] == .8
    assert prov["reference_overlap_checked"] is False
    candidates = pd.DataFrame({"source": ["STAT3", "STAT3"],
                               "target": ["SOCS3", "FOS"]})
    mapped = align_evidence(candidates, links)
    assert mapped.m_independent_prior.tolist() == [1, 0]
    assert np.isnan(mapped.independent_prior_score.iloc[1])
    with pytest.raises(ValueError, match="synthetic BEELINE"):
        read_external_evidence(evidence, metadata, ["g1", "g2"],
                               dataset_kind="synthetic")
    with pytest.raises(ValueError, match="synthetic"):
        read_external_evidence(evidence, metadata, ["g1", "g2"],
                               dataset_kind="biological")
    with pytest.raises(ValueError, match="independent"):
        meta["independent_of_evaluation_reference"] = False
        metadata.write_text(json.dumps(meta))
        read_external_evidence(evidence, metadata, ["STAT3", "SOCS3"],
                               dataset_kind="biological")


def test_gold_filename_cannot_be_used_as_prior(tmp_path):
    path = tmp_path / "GroundTruthNetwork.csv"
    path.write_text("source,target,score,evidence_id\nX,Y,1,xx\n")
    prov = tmp_path / "source.json"
    prov.write_text("{}")
    with pytest.raises(ValueError, match="Reference"):
        read_external_evidence(path, prov, ["X", "Y"], dataset_kind="biological")


def test_end_to_end_without_any_gold_or_independent_priors(tmp_path):
    x, pt = _toy(150)
    raw = tmp_path / "synthetic"
    raw.mkdir()
    x.T.to_csv(raw / "ExpressionData.csv")
    pt.to_csv(raw / "PseudoTime.csv")
    assert not (raw / "GroundTruthNetwork.csv").exists()
    root = tmp_path / "results"
    m = build(raw, root / "synthetic", max_parents=1)
    assert m["no_gold_file_opened"]
    assert not m["independent_biological_prior_present"]
    assert m["outer_score_cells_excluded"] == 30
    assert m["outer_test_cells_excluded"] == 30
    assert len(pd.read_csv(root / "synthetic" / "m2cd_targets.csv")) == 36
    metrics = audit(root)
    assert len(metrics) == 1 and metrics.gold_labels_used.iloc[0] == False
    assert metrics.independent_evidence_edges.iloc[0] == 0
    assert np.isfinite(metrics.conditional_fold_spearman.iloc[0])
