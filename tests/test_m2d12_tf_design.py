import pandas as pd
import pytest
from scripts.audit_m2d12_tf_design import summarize


def test_frequency_and_joint_counts_are_observed_not_inferred():
    obs = pd.DataFrame(
        {"TF": ["SOX2", "SOX2", "GFP", "NANOG"],
         "batch": ["A", "B", "A", "A"],
         "temp": [1, 1, 2, 2]},
        index=["c1", "c2", "c3", "c4"])
    result = summarize(obs)
    assert result["n_cells"] == 4
    assert len(result["tf_batch_temp_cells"]) == 4
    assert result["batch_count_per_tf"][0] == {"TF": "SOX2", "n_batches": 2}
    assert "GFP" in result["possible_control_label_candidates_not_verified"]
    assert result["experimental_controls_verified"] is False
    assert result["biological_replicates_verified"] is False


def test_reject_missing_columns_and_duplicate_cells():
    data = pd.DataFrame({"TF": ["A"], "batch": ["B"]}, index=["x"])
    with pytest.raises(ValueError, match="Missing"):
        summarize(data)
    dup = pd.DataFrame({"TF": ["A", "B"], "batch": ["X", "X"],
                        "temp": [1, 2]}, index=["x", "x"])
    with pytest.raises(ValueError, match="not unique"):
        summarize(dup)
