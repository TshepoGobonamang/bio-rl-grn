"""M2B leakage tests and genuine network-held-out reliability-gate checks."""
import json
import numpy as np
import pandas as pd
import pytest
from biorlgrn.m2_reliability import (
    ALPHA, FEATURES, KEYS, TARGET, load_m2a, leave_one_network_out,
)
from scripts.evaluate_m2b_grn import evaluate_network
from biorlgrn.demo import make_demo


def _fake_m2a(root):
    for network, g in (("net-a", 3), ("net-b", 4), ("net-c", 5)):
        d = root / network
        d.mkdir(parents=True, exist_ok=True)
        rows, ys = [], []
        for fold in range(3):
            for i in range(g):
                for j in range(g):
                    if i == j:
                        continue
                    rng = np.random.default_rng(100 + fold * 8 + i * 7 + j * 3 + g * 13)
                    z = float(rng.uniform(-0.7, .7))
                    agree = float(rng.uniform(.05, .95))
                    keys = {"dataset_id": network, "fold_id": fold,
                            "source": f"G{i}", "target": f"G{j}",
                            "pair_id": "::".join(sorted([f"G{i}", f"G{j}"]))}
                    vals = dict(zip(FEATURES, [
                        z, 1, 1., 1, 16, abs(z), .3, 1., g, 1, .005, agree
                    ]))
                    rows.append({**keys, **vals})
                    ys.append({**keys, TARGET: np.clip(.15 + .65 * agree + .10 * z, 0., 1.)})
        f = pd.DataFrame(rows, columns=KEYS + FEATURES)
        t = pd.DataFrame(ys, columns=KEYS + [TARGET])
        f.to_csv(d / "m2_features.csv", index=False)
        t.to_csv(d / "m2_targets.csv", index=False)
        (d / "m2_crossfit_manifest.json").write_text(json.dumps({
            "uses_gold_reference": False,
            "outer_score_and_test_excluded_from_dataset": True,
            "feature_columns": FEATURES, "n_rows": len(f),
            "folds": [{"disjoint": True} for _ in range(3)],
        }))
    return root


def test_leave_one_topology_out_outputs_bounded_gates(tmp_path):
    dataset = load_m2a(_fake_m2a(tmp_path / "m2a"))
    pred, score, weights = leave_one_network_out(dataset)
    assert len(score) == 3 * 4
    assert not score.trained_on_held_out_network.any()
    assert pred.rho.between(0, 1).all()
    assert weights.rho.between(0, 1).all()
    assert len(weights) == (3 * 2 + 4 * 3 + 5 * 4) * 4
    assert set(score.method) == {"fixed", "global", "stability", "adaptive"}
    assert score.loc[score.method == "adaptive", "n_eval_rows"].min() > 0


def test_held_out_target_perturbation_cannot_change_its_weights(tmp_path):
    source = _fake_m2a(tmp_path / "m2a")
    original = load_m2a(source)
    _, _, before = leave_one_network_out(original)
    path = source / "net-b" / "m2_targets.csv"
    t = pd.read_csv(path)
    t[TARGET] = np.clip(1 - t[TARGET], 0, 1)
    t.to_csv(path, index=False)
    modified = load_m2a(source)
    _, _, after = leave_one_network_out(modified)
    keys = ["dataset_id", "source", "target", "method"]
    lhs = before.loc[before.dataset_id.eq("net-b")].sort_values(keys).reset_index(drop=True)
    rhs = after.loc[after.dataset_id.eq("net-b")].sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(lhs, rhs)


def test_reject_accidental_target_column_in_features(tmp_path):
    source = _fake_m2a(tmp_path / "m2a")
    path = source / "net-a" / "m2_features.csv"
    f = pd.read_csv(path)
    f[TARGET] = .99
    f.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Feature file"):
        load_m2a(source)


def test_missing_temporal_evidence_has_zero_gate(tmp_path):
    source = _fake_m2a(tmp_path / "m2a")
    fp, tp = source / "net-c" / "m2_features.csv", source / "net-c" / "m2_targets.csv"
    f, t = pd.read_csv(fp), pd.read_csv(tp)
    f.loc[0, "z_signed"] = np.nan
    f.loc[0, "m_present"] = 0
    t.loc[0, TARGET] = np.nan
    f.to_csv(fp, index=False)
    t.to_csv(tp, index=False)
    data = load_m2a(source)
    pred, score, weights = leave_one_network_out(data)
    missing = pred.loc[pred.dataset_id.eq("net-c") & pred.fold_id.eq(0)
                       & pred.source.eq(f.loc[0, "source"]) &
                       pred.target.eq(f.loc[0, "target"])]
    assert (missing.rho == 0).all()


def test_grn_evaluation_requires_weights_for_all_edges(tmp_path):
    toy = make_demo(tmp_path / "toy")
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.read_csv(toy / "expression.csv", index_col=0).T.to_csv(
        raw / "ExpressionData.csv")
    pd.read_csv(toy / "gold_edges.csv").rename(
        columns={"source": "Gene1", "target": "Gene2"}
    ).to_csv(raw / "GroundTruthNetwork.csv", index=False)
    x = pd.read_csv(toy / "expression.csv", index_col=0)
    pd.DataFrame({"PseudoTime1": np.arange(len(x))}, index=x.index).to_csv(
        raw / "PseudoTime.csv")
    bad = tmp_path / "weights.csv"
    pd.DataFrame([{"dataset_id": "raw", "source": "G1", "target": "G2",
                   "method": "adaptive", "rho": .5}]).to_csv(bad, index=False)
    with pytest.raises(ValueError, match="Missing or duplicate weights"):
        evaluate_network(raw, bad, tmp_path / "outputs", bins=8, max_edges=4)
