import gzip
import shutil
import numpy as np
import pandas as pd
import pytest

from scripts.inspect_m2d10_h5ad import inspect


def test_real_h5ad_structure_not_fake_assay_truth(tmp_path):
    ad = pytest.importorskip("anndata")
    obj = ad.AnnData(X=np.array([[1, 2], [3, 4]], dtype=np.float32),
                     obs=pd.DataFrame({"library": ["s1", "s2"]}, index=["c1", "c2"]),
                     var=pd.DataFrame(index=["A", "B"]))
    raw = tmp_path / "test.h5ad"
    obj.write_h5ad(raw)
    gz = tmp_path / "test.h5ad.gz"
    with raw.open("rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    report = inspect(gz)
    assert report["n_cells"] == 2
    assert report["n_genes"] == 2
    assert report["obs_columns"] == ["library"]
    assert len(report["source_sha256"]) == 64
    assert not report["controls_verified"]
    assert not report["tf_barcode_join_verified"]
    assert not report["regulatory_labels_computed"]


def test_duplicate_cell_ids_rejected(tmp_path):
    ad = pytest.importorskip("anndata")
    obj = ad.AnnData(X=np.ones((2, 1)),
                     obs=pd.DataFrame(index=["duplicate", "duplicate"]),
                     var=pd.DataFrame(index=["G"]))
    raw = tmp_path / "bad.h5ad"
    obj.write_h5ad(raw)
    with pytest.raises(ValueError, match="Nonunique"):
        inspect(raw)
