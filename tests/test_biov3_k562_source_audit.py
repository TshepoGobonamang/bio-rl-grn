import gzip
from pathlib import Path
import pytest
from scripts.audit_biov3_k562_scRNA_source import audit


def write_gz(path: Path, text: str):
    with gzip.open(path,"wt") as f:
        f.write(text)


def test_mixed_cell_line_metadata_reports_explicit_k562_and_no_training_claim(tmp_path):
    expression=tmp_path/"expr.txt.gz"
    metadata=tmp_path/"meta.txt.gz"
    write_gz(expression,"gene\tcell1\tcell2\tcell3\nG1\t1\t0\t2\nG2\t0\t4\t0\n")
    write_gz(metadata,"cell_id\tcell_line\tbatch\ncell1\tK562\tA\ncell2\tBT142\tA\ncell3\tK562\tB\n")
    r=audit(expression,metadata)
    assert r["matrix_shape_as_text_rows_after_first"]==2
    assert r["metadata_rows"]==3
    assert r["metadata_exact_k562_label_count_by_candidate_field"]["cell_line"]==2
    assert r["metadata_fields_matching_expression_header"]["cell_id"]==3
    assert not r["independent_training_expression_qualified"]
    assert not r["encode_perturbation_response_accessed"]


def test_uninformative_sample_does_not_infer_k562_from_filename(tmp_path):
    expression=tmp_path/"K562_expression.txt.gz"
    metadata=tmp_path/"meta.txt.gz"
    write_gz(expression,"gene\ta\tb\nG\t1\t2\n")
    write_gz(metadata,"barcode\tcluster\na\t1\nb\t2\n")
    r=audit(expression,metadata)
    assert not r["metadata_categorical_cell_line_candidates"]
    assert not r["independent_training_expression_qualified"]


def test_reject_duplicate_source_headers(tmp_path):
    expression=tmp_path/"expr.txt.gz"
    metadata=tmp_path/"meta.txt.gz"
    write_gz(expression,"gene\ta\ta\nG\t1\t2\n")
    write_gz(metadata,"barcode\tcell_line\na\tK562\n")
    with pytest.raises(ValueError,match="duplicated"):
        audit(expression,metadata)
