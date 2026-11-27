"""Scientific contracts for the BIO-V7 TF-constrained one-edit adapter."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grnbench.core import _loss
from grnbench.gcev import FEATURES


_path = Path(__file__).parents[1] / "scripts" / "biov7_constrained_models.py"
_spec = importlib.util.spec_from_file_location("biov7_constrained_models", _path)
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)


def fixture():
    rng = np.random.default_rng(923)
    genes = ["TF0", "TF1", "TF2", "TF3", "TF4", "T0", "T1", "T2"]
    values = rng.normal(size=(96, len(genes)))
    values[:, 5] = .8 * values[:, 0] + rng.normal(size=96)
    # An ineligible target is more correlated with T0 than any TF.
    values[:, 6] = values[:, 5] + rng.normal(scale=.01, size=96)
    data = pd.DataFrame(values, columns=genes,
                        index=[f"cell_{i}" for i in range(96)])
    edges = [(tf, target) for target in genes for tf in genes[:5] if tf != target]
    motif = np.linspace(0, 1, len(edges))
    target_values = {target: i / len(genes) for i, target in enumerate(genes)}
    dnase = np.asarray([target_values[t] for _, t in edges])
    evidence = {"motif": (motif, np.ones(len(edges), bool)),
                "dnase": (dnase, np.ones(len(edges), bool))}
    return data.iloc[:72], data.iloc[72:], edges, evidence


def test_constant_target_source_has_no_action_or_utility_and_scaffolds_are_tfs():
    fit, _, edges, evidence = fixture()
    labels = model.nested_valuation(fit, edges, evidence)
    assert set(FEATURES).issubset(labels.columns)
    assert len(labels) == 3 * 8 * 2
    dnase = labels[labels.evidence_source == "dnase"]
    assert dnase.source_constant_within_target.all()
    assert dnase.different_edge.eq(0).all()
    assert dnase.source_value.eq(0).all()
    assert dnase.minus_edge.equals(dnase.plus_edge)
    assert dnase.minus_validation_nmse.equals(dnase.plus_validation_nmse)
    for row in labels.itertuples():
        state = row.graph_parents.split("|")
        assert len(state) == 2
        assert all(tf.startswith("TF") and tf != row.target for tf in state)
        assert row.minus_edge not in state and row.plus_edge not in state


def test_final_contexts_preserve_universe_and_equal_action_budget():
    fit, _, edges, _ = fixture()
    contexts = model.final_contexts(fit, edges)
    assert list(zip(contexts.tf, contexts.target)) == edges
    assert contexts.scaffold_edge.sum() == 2 * fit.shape[1]
    assert contexts.eligible_action.sum() == len(edges) - 2 * fit.shape[1]
    assert contexts.loc[contexts.scaffold_edge, ["gain", "expr_rank"]].isna().all().all()
    assert np.isfinite(contexts.loc[contexts.eligible_action, ["gain", "expr_rank"]]).all().all()
    for row in contexts.itertuples():
        assert all(tf.startswith("TF") for tf in row.graph_parents.split("|"))


def test_evaluation_matches_training_variance_loss_and_rejects_cell_overlap():
    fit, evaluation, edges, _ = fixture()
    contexts = model.final_contexts(fit, edges)
    eligible = contexts.loc[contexts.eligible_action]
    selections = (eligible.sort_values(["target", "expr_rank", "tf"],
                                       ascending=[True, False, True])
                  .drop_duplicates("target")
                  .rename(columns={"tf": "selected_tf"}))
    selections["method"] = "expression"
    per_target, errors = model.evaluate_graphs(fit, evaluation, selections)
    assert per_target.total_parent_count.eq(3).all()
    assert errors["expression"].shape == (24, 8)
    for row in per_target.itertuples():
        parents = tuple(row.graph_parents.split("|")) + (row.selected_tf,)
        assert row.nmse == pytest.approx(_loss(fit, evaluation, row.target, parents, 1.0))
    assert errors["expression"].mean() == pytest.approx(per_target.nmse.mean())
    with pytest.raises(ValueError, match="overlap"):
        model.evaluate_graphs(fit, fit.iloc[:10], selections)


def test_evaluation_rejects_unequal_scaffold_between_methods():
    fit, evaluation, edges, _ = fixture()
    contexts = model.final_contexts(fit, edges)
    selected = contexts.loc[contexts.eligible_action].drop_duplicates("target").copy()
    selected = selected.rename(columns={"tf": "selected_tf"})
    selected["method"] = "first"
    other = selected.copy()
    other["method"] = "second"
    other.loc[other.index[0], "graph_parents"] = "TF3|TF4"
    assert selected.loc[selected.index[0], "graph_parents"] != "TF3|TF4"
    with pytest.raises(ValueError, match="unequal scaffolds"):
        model.evaluate_graphs(fit, evaluation, pd.concat([selected, other]))
