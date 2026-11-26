import gzip
from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.biov3_reconcile_k562_labels import audit_identity,export_matrix

def test_exact_author_labels_select_k562(tmp_path):
    m={'a':'1','b':'2','c':'1','d':'2'}
    a={'a':'BT142','b':'K562','c':'BT142','d':'K562'}
    d=audit_identity(m,a)
    assert d['code_to_line']=={'1':'BT142','2':'K562'}
    expr=tmp_path/'source.gz'
    with gzip.open(expr,'wt') as f:
        f.write('a\tb\tc\td\nGENE1\t1\t2\t3\t4\nGENE2\t0\t1\t0\t2\n')
    r=export_matrix(expr,m,a,tmp_path/'out')
    assert r['k562_cell_count']==2
    assert r['n_genes_unmapped_symbols']==2
    assert r['eligible_for_gcev_model_fit'] is False
    with gzip.open(tmp_path/'out'/'K562_gene_by_cell_original_counts.tsv.gz','rt') as f:
        assert f.readline().strip()=='gene\tb\td'

def test_unknown_numeric_assignment_not_inferred():
    with pytest.raises(ValueError,match='barcode mismatch'):
        audit_identity({'a':'1','b':'2'},{'a':'K562','c':'BT142'})

def test_mixed_numeric_code_records_disagreement_without_relabeling():
    r=audit_identity({'a':'1','b':'1','c':'2','d':'2'},
                     {'a':'K562','b':'BT142','c':'K562','d':'BT142'})
    assert r['numeric_code_maps_directly_to_cell_line'] is False
    assert r['code_to_line'] is None
    assert r['author_counts']=={'K562':2,'BT142':2}

def test_wrong_gene_row_length_rejected(tmp_path):
    expr=tmp_path/'source.gz'
    with gzip.open(expr,'wt') as f: f.write('a\tb\nc1\t2\nc2\t3\t4\n')
    with pytest.raises(ValueError,match='row length mismatch'):
        export_matrix(expr,{'a':'1','b':'2'},{'a':'BT142','b':'K562'},tmp_path/'out')
