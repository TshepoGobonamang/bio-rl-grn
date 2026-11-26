import gzip
from pathlib import Path
import pytest

from scripts.diagnose_biov3_k562_numeric_labels import audit


def gz(path,content):
    with gzip.open(path,"wt") as f:f.write(content)


def test_original_numeric_labels_not_assumed_to_be_cell_lines(tmp_path):
    expr=tmp_path/"expr.gz"
    meta=tmp_path/"meta.gz"
    gz(expr,"c1\tc2\tc3\tc4\nHBG2\t0\t10\t0\t12\nPTPRZ1\t8\t0\t7\t0\n")
    gz(meta,"orig.ident\tCellType\nc1\t1\nc2\t2\nc3\t1\nc4\t2\n")
    r=audit(expr,meta)
    assert r["n_exactly_joined_cells"]==4
    assert r["marker_expression_by_original_numeric_label"]["HBG2"]["2"]["mean_abundance"]==11
    assert r["marker_expression_by_original_numeric_label"]["PTPRZ1"]["1"]["mean_abundance"]==7.5
    assert not r["label_to_K562_assigned"]
    assert not r["encode_crispri_outcomes_accessed"]

def test_one_unmatched_barcode_is_not_silently_corrected(tmp_path):
    expr=tmp_path/"expr.gz"
    meta=tmp_path/"meta.gz"
    header="\t".join(f"c{i}" for i in range(100))
    hbg="\t".join(["HBG2"]+["1"]*100)
    pt="\t".join(["PTPRZ1"]+["2"]*100)
    gz(expr,header+"\n"+hbg+"\n"+pt+"\n")
    text="orig.ident\tCellType\n"+"".join(f"c{i}\t{1 if i<50 else 2}\n" for i in range(1,100))
    text+="missing1\t1\n"
    gz(meta,text)
    r=audit(expr,meta)
    assert r["n_exactly_joined_cells"]==99
    assert r["n_unmatched_expression_cells"]==1
    assert r["n_unmatched_metadata_cells"]==1

def test_wrong_matrix_orientation_fails_closed(tmp_path):
    expr=tmp_path/"expr.gz";meta=tmp_path/"meta.gz"
    gz(expr,"gene\tc1\tc2\nHBG2\t1\t2\n")
    gz(meta,"orig.ident\tCellType\nc1\t1\nc2\t2\n")
    with pytest.raises(ValueError,match="orientation"):
        audit(expr,meta)
