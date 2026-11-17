import pandas as pd
import pytest
from scripts.prepare_m2d5_perturbation_responses import inspect_metadata, paired_responses


def fixture_data():
    rows, values = [], []
    for tf in ["TF_A"]:
        for rep in ["bio1", "bio2"]:
            for condition in ["control", "overexpression"]:
                for cell in range(2):
                    cid = f"{tf}_{rep}_{condition}_{cell}"
                    rows.append({"cell_id": cid, "tf": tf, "condition": condition,
                                 "replicate_id": rep, "batch": "batch_1"})
                    values.append({"cell_id": cid,
                                   "GENE_A": 40 if condition == "overexpression" else 10,
                                   "GENE_B": 10 if condition == "overexpression" else 40})
    meta = pd.DataFrame(rows)
    counts = pd.DataFrame(values).set_index("cell_id")
    return counts, meta


def test_requires_explicit_verification():
    counts, meta = fixture_data()
    with pytest.raises(ValueError, match="not independently verified"):
        paired_responses(counts, meta, min_cells=2)


def test_paired_effect_and_no_binary_labels():
    counts, meta = fixture_data()
    effects, manifest = paired_responses(counts, meta,
        biological_replicates_verified=True, min_cells=2)
    assert len(effects) == 2
    assert (effects.n_biological_pairs == 2).all()
    assert effects.loc[effects.target == "GENE_A",
                       "mean_log2_cpm_difference"].iloc[0] > 0
    assert effects.binary_response_label.isna().all()
    assert manifest["labels_generated"] is False


def test_missing_pair_cannot_count_cells_as_replicates():
    counts, meta = fixture_data()
    meta = meta.loc[~((meta.replicate_id == "bio2") &
                      (meta.condition == "control"))]
    counts = counts.loc[meta.cell_id]
    effects, _ = paired_responses(counts, meta,
        biological_replicates_verified=True, min_cells=2)
    assert effects.empty


def test_no_implicit_control_or_negative_counts():
    counts, meta = fixture_data()
    meta.loc[0, "condition"] = "unassigned"
    with pytest.raises(ValueError, match="explicit control"):
        inspect_metadata(meta)
    counts, meta = fixture_data()
    counts.iloc[0, 0] = -1
    with pytest.raises(ValueError, match="raw finite"):
        paired_responses(counts, meta, biological_replicates_verified=True, min_cells=2)
