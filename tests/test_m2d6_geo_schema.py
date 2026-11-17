import gzip
import json
import pandas as pd
import pytest

from scripts.inspect_m2d6_geo_schema import inspect_tfmap


def test_actual_header_preserved_without_guessing(tmp_path):
    path = tmp_path / "map.csv.gz"
    with gzip.open(path, "wt") as out:
        out.write("cell_barcode,TF_construct,putative_batch\n")
        out.write("AA-1,SOX2,libraryA\nBB-1,OCT4,libraryB\n")
    result = inspect_tfmap(path)
    assert result["n_rows"] == 2
    assert result["columns"] == ["cell_barcode", "TF_construct", "putative_batch"]
    assert len(result["sha256"]) == 64
    assert not result["control_identity_verified"]
    assert not result["biological_replicates_verified"]
    assert not result["barcode_key_verified"]


def test_empty_input_fails_closed(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("barcode,tf\n", encoding="utf-8")
    with pytest.raises(ValueError, match="empty"):
        inspect_tfmap(path)


def test_preserve_identifier_strings(tmp_path):
    path = tmp_path / "map.csv"
    path.write_text("barcode,tf\n001,TF01\n002,TF02\n", encoding="utf-8")
    result = inspect_tfmap(path)
    assert result["distinct_per_column"]["barcode"] == 2
