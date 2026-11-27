"""Boundary, score-integrity, and algebraic regression tests using synthetic data."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts import biov8_freeze_and_audit as audit


@pytest.fixture(scope="module")
def synthetic_frozen_inputs(tmp_path_factory):
    """Build a complete fake 16-TF freeze; no biological input is read."""
    root = tmp_path_factory.mktemp("biov8_synthetic")
    results = root / "results"
    results.mkdir()
    tfs = [f"TF{i:02d}" for i in range(16)]
    genes = tfs + [f"TG{i:02d}" for i in range(48)]
    ids = {g: f"SYNTHETIC{i:04d}.1" for i, g in enumerate(genes)}
    rows, contexts = [], []
    for j, target in enumerate(genes):
        candidates = [tf for tf in tfs if tf != target]
        parents = candidates[:2]
        a = 0. if j % 4 == 0 else (j % 5 + 1) / 10
        for i, tf in enumerate(candidates):
            m = ((i * 7 + j) % 19) / 18
            rows.append({"tf": tf, "target": target, "tf_gencode_gene_id": ids[tf],
                         "target_gencode_gene_id": ids[target], "motif_prior_score": repr(m),
                         "dnase_accessibility_mean": repr(a), "joint_motif_dnase_product": repr(m * a),
                         "accessibility_source_observed": "True", "direct_tf_binding_verified": "False",
                         "independent_regulatory_edge_truth": "False"})
            e = (i + 1) / 18
            scaffold = tf in parents
            contexts.append({"tf": tf, "target": target, "graph_parents": "|".join(parents),
                             "scaffold_edge": str(scaffold), "eligible_action": str(not scaffold),
                             "gain": "" if scaffold else format(e / 10, ".17g"),
                             "expr_rank": "" if scaffold else format(e, ".17g"),
                             "scaffold_selected_on_full_fit": "True", "gain_coefficient_folds": "3"})
    prior = pd.DataFrame(rows)
    contexts = pd.DataFrame(contexts)
    priors = {}
    input_hashes = {}
    provenance = {}
    for variant, (config_key, input_key) in audit.PRIOR_KEYS.items():
        path = root / f"{variant}.csv"
        prior.to_csv(path, index=False)
        priors[variant] = path
        input_hashes[input_key] = audit.sha(path)
        provenance[config_key] = audit.sha(path)
    scores, selections = [], []
    residual = contexts.loc[contexts.eligible_action == "True"]
    for method_index, method in enumerate(audit.METHODS):
        for target, group in residual.groupby("target", sort=True):
            group_rows = []
            for row in group.to_dict("records"):
                # Target-level constants preserve decisions but alter legacy
                # across-target scores; this deliberately exercises both fields.
                offset = method_index * (.001 if target.endswith("0") else .002)
                record = {k: row[k] for k in ["tf", "target", "graph_parents", "gain", "expr_rank"]}
                record.update({"method": method, "score": format(float(row["expr_rank"]) + offset, ".17g"),
                               "action_score": row["expr_rank"],
                               "constant_target_offsets_removed_for_argmax": "True",
                               "global_score_is_not_source_utility_endpoint": "True"})
                group_rows.append(record)
                scores.append(record)
            best = min(group_rows, key=lambda x: (-float(x["action_score"]), x["tf"]))
            selections.append({"method": method, "target": target, "graph_parents": best["graph_parents"],
                "selected_tf": best["tf"], "selected_score": best["score"],
                "selected_action_score": best["action_score"], "parent_count_before": "2",
                "inserted_edge_count": "1", "total_parent_count": "3"})
    scores, selections = pd.DataFrame(scores), pd.DataFrame(selections)
    output_hashes = {}
    for transform in audit.TRANSFORMS:
        path = results / transform / "candidate_contexts.csv"
        path.parent.mkdir(parents=True)
        contexts.to_csv(path, index=False)
        output_hashes[str(path.relative_to(results))] = audit.sha(path)
        for variant in audit.VARIANTS:
            case = results / transform / variant
            case.mkdir()
            for name, frame in [("candidate_scores.csv", scores), ("frozen_selections.csv", selections)]:
                path = case / name
                frame.to_csv(path, index=False)
                output_hashes[str(path.relative_to(results))] = audit.sha(path)
    v7 = {"transforms": audit.TRANSFORMS, "variants": audit.VARIANTS, "methods": audit.METHODS,
          "boost": .3, "frozen_genes": 64, "frozen_regulators": 16, "frozen_edges": 1008,
          "eligible_residual_actions": 880, "provenance": provenance}
    v7_path = root / "biov7_protocol.json"
    v7_path.write_text(json.dumps(v7))
    freeze = {
        "stage": "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE",
        "protocol_sha256": audit.sha(v7_path), "candidate_prior_edges": 1008,
        "eligible_single_edits": 880, "n_genes": 64, "n_TFs": 16, "scaffold_edges": 128,
        "inserted_edges_per_method": 64, "total_graph_edges_per_method": 192,
        "fit_cells_used": 394, "internal_score_cells_used": 0, "test_cells_used": 0,
        "test_used": False, "test_expression_analytics_performed": False,
        "perturbation_outcomes_accessed": False, "methods": audit.METHODS,
        "transforms": audit.TRANSFORMS, "variants": audit.VARIANTS,
        "input_sha256": input_hashes, "output_sha256": output_hashes,
    }
    freeze_path = results / "freeze_manifest.json"
    freeze_path.write_text(json.dumps(freeze))
    protocol = {"protocol_id": "SYNTHETIC-BOUNDARY-TEST", "boost": .3,
                "diagnostic_permutation_seed": 20261128, "minimum_measured_targets_per_tf": 20,
                "minimum_qualified_tfs": 6, "biov7_protocol_sha256": audit.sha(v7_path),
                "biov7_freeze_manifest_sha256": audit.sha(freeze_path), "input_sha256": provenance,
                "regulators": tfs, "regulator_ensembl_ids": {tf: ids[tf].split(".")[0] for tf in tfs},
                "constraints": {name: True for name in ["no_expression_access", "no_internal_score_access",
                    "untouched_test_remains_sealed", "no_outcome_access", "no_new_prediction_fitting",
                    "no_new_performance_evaluation"]}}
    protocol_path = root / "protocol.json"
    protocol_path.write_text(json.dumps(protocol))
    return {"root": root, "results": results, "priors": priors, "v7": v7_path,
            "protocol": protocol_path, "contexts": contexts, "scores": scores,
            "selections": selections, "prior": prior}


def test_constant_offsets_cannot_create_actions_and_product_is_attenuation():
    contexts = pd.DataFrame({"tf": ["A", "B", "A", "B"], "target": ["T0", "T0", "T1", "T1"],
                             "eligible_action": [True] * 4, "expr_rank": [.5, .49, .5, .49],
                             "graph_parents": ["P|Q"] * 4})
    prior = pd.DataFrame({"tf": contexts.tf, "target": contexts.target,
                          "motif_prior_score": [0., 1., 0., 1.],
                          "dnase_accessibility_mean": [0., 0., .5, .5]})
    diag, permutation = audit.diagnose_contexts(contexts, prior, seed=20261128, boost=.3)
    by = diag.set_index(["target", "representation"])
    for target in ["T0", "T1"]:
        assert by.loc[(target, "accessibility_only"), "selected_tf"] == "A"
        assert by.loc[(target, "legacy_additive_motif_accessibility"), "selected_tf"] == "B"
        assert by.loc[(target, "motif"), "selected_tf"] == "B"
    assert by.loc[("T0", "matched_accessibility_motif"), "selected_tf"] == "A"
    assert by.loc[("T0", "global_mean_accessibility_motif"), "selected_tf"] == "B"
    assert by.loc[("T1", "matched_accessibility_motif"), "selected_tf"] == "B"
    assert sorted(permutation.original_accessibility) == sorted(permutation.permuted_accessibility)
    assert not diag.predictive_usefulness_evaluated.any()
    assert not diag.independent_third_source.any()


def test_coefficient_boundary_and_lexical_ties():
    boundary, challenger, wins_tie = audit.beta_boundary(["B", "A"], np.array([.8, .7]), np.array([0., 1.]))
    assert boundary == pytest.approx(.1)
    assert challenger == "A" and wins_tie
    assert audit.beta_boundary(["B", "A"], np.array([.8, .7]), np.array([.2, .2])) == (None, "", False)
    assert audit.winner(["Z", "A"], np.array([.5, .5])) == 1


def test_observed_zero_cannot_be_reclassified_as_missing_or_gold(synthetic_frozen_inputs):
    table = synthetic_frozen_inputs["prior"].copy()
    table.loc[0, "accessibility_source_observed"] = "False"
    with pytest.raises(ValueError, match="remain observed"):
        audit.validate_prior(table)
    table.loc[0, "accessibility_source_observed"] = "True"
    table.loc[0, "independent_regulatory_edge_truth"] = "True"
    with pytest.raises(ValueError, match="biological truth"):
        audit.validate_prior(table)


def test_no_fabricated_scaffold_score_or_inconsistent_graph_allowed(synthetic_frozen_inputs):
    fixture = synthetic_frozen_inputs
    prior = audit.validate_prior(fixture["prior"])
    bad = fixture["contexts"].copy()
    scaffold_row = bad.index[bad.scaffold_edge == "True"][0]
    bad.loc[scaffold_row, "gain"] = "0"
    with pytest.raises(ValueError, match="remain absent"):
        audit.validate_contexts(bad, prior)
    contexts = audit.validate_contexts(fixture["contexts"], prior)
    selections = fixture["selections"].copy()
    selections.loc[0, "selected_tf"] = selections.loc[0, "graph_parents"].split("|")[0]
    with pytest.raises(ValueError, match="argmax"):
        audit.validate_scores_and_selections(fixture["scores"], selections, contexts)


def test_complete_synthetic_export_preserves_exact_scores_and_never_reads_nearby_outcomes(
        synthetic_frozen_inputs, tmp_path, monkeypatch):
    fixture = synthetic_frozen_inputs
    forbidden = {"evaluation_manifest.json", "per_target_results.csv", "per_cell_errors.npz", "source_utility_labels.csv"}
    for name in forbidden:
        (fixture["results"] / name).write_text("DO NOT READ\n")
    actual_sha = audit.sha
    accessed = []

    def checked_sha(path):
        assert Path(path).name not in forbidden, "Outcome/utility file was accessed"
        accessed.append(str(path))
        return actual_sha(path)

    monkeypatch.setattr(audit, "sha", checked_sha)
    output = tmp_path / "frozen"
    result = audit.run(fixture["protocol"], fixture["results"], fixture["priors"], output,
                       biov7_protocol_path=fixture["v7"])
    assert result["frozen_score_rows"] == 49280
    assert result["frozen_graph_decisions"] == 3584
    assert result["source_action_diagnostic_rows"] == 2688
    assert result["source_action_summary_rows"] == 42
    assert result["expression_matrices_read"] == result["new_models_fitted"] == 0
    assert result["untouched_test_cells_accessed"] == 0
    assert not result["biological_validation_performed"]
    exported = audit.read_csv(output / "frozen_biov7_candidate_scores.csv")
    first = exported.loc[(exported["transform"] == "log1p_original") & (exported.variant == "primary_mixed")]
    assert first.score.tolist() == fixture["scores"].score.tolist()
    assert first.action_score.tolist() == fixture["scores"].action_score.tolist()
    assert (first.score != first.action_score).any(), "Legacy target offsets were lost"
    coverage = pd.read_csv(output / "tf_candidate_coverage.csv")
    assert len(coverage) == 16 and coverage.external_measured_target_count.isna().all()
    summary = pd.read_csv(output / "source_action_summary.csv")
    assert summary.loc[summary.representation == "accessibility_only", "changed_vs_expression"].eq(0).all()
    assert summary.loc[summary.representation == "legacy_additive_motif_accessibility", "changed_vs_motif"].eq(0).all()
    assert len(result["input_sha256"]) == 24


def test_bad_frozen_hash_fails_before_tables_or_output(synthetic_frozen_inputs, tmp_path, monkeypatch):
    fixture = synthetic_frozen_inputs
    protocol = json.loads(fixture["protocol"].read_text())
    protocol["biov7_freeze_manifest_sha256"] = "0" * 64
    path = tmp_path / "bad_protocol.json"
    path.write_text(json.dumps(protocol))
    output = tmp_path / "must_not_exist"

    def prohibited(*args, **kwargs):
        raise AssertionError("Tables read after upstream integrity failure")

    monkeypatch.setattr(audit, "read_csv", prohibited)
    with pytest.raises(ValueError, match="checksum mismatch"):
        audit.run(path, fixture["results"], fixture["priors"], output,
                  biov7_protocol_path=fixture["v7"])
    assert not output.exists()


def test_relaxing_outcome_boundary_is_rejected(synthetic_frozen_inputs, tmp_path):
    protocol = json.loads(synthetic_frozen_inputs["protocol"].read_text())
    protocol["constraints"]["no_outcome_access"] = False
    path = tmp_path / "bad_protocol.json"
    path.write_text(json.dumps(protocol))
    with pytest.raises(ValueError, match="no_outcome_access"):
        audit.load_protocol(path)
