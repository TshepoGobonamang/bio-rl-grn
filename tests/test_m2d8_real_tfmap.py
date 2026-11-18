import gzip
import pandas as pd
import pytest
from scripts.observe_m2d8_real_tfmap import audit_source


def test_gzip_source_observation_preserves_schema_and_no_biology_claim(tmp_path):
    source = tmp_path / "source.csv.gz"
    with gzip.open(source, "wt") as f:
        f.write("cell,TF_construct,sample\n001,SOX2,A\n002,OCT4,A\n003,,B\n")
    result = audit_source(source)
    assert result["n_rows"] == 3
    assert result["columns"] == ["cell", "TF_construct", "sample"]
    assert result["field_profiles"]["TF_construct"]["non_missing"] == 2
    assert len(result["sha256"]) == 64
    assert result["real_file_parsed"]
    assert not result["replicate_independence_verified"]
    assert not result["biological_labels_generated"]


def test_empty_source_rejected(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("cell,tf\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Empty"):
        audit_source(p)
