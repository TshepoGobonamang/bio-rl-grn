import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.sparse import csc_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.biov3_k562r_singleline_intake import qualify_matrix

def test_fit_only_selection_and_disjoint_partitions():
    rng=np.random.default_rng(1)
    m=csc_matrix(rng.poisson(2,size=(10,20)))
    x,part,genes=qualify_matrix(m,[f'ENSG{i:09d}' for i in range(10)],
                                 [f'g{i}' for i in range(10)],
                                 [f'cell-{i}' for i in range(20)],n_cells=15,n_genes=4)
    assert len(x)==15 and len(genes)==4
    assert part['partition'].value_counts().to_dict()=={'fit':9,'score':3,'test':3}
    assert x['cell_id'].is_unique

def test_bad_source_counts_rejected():
    a=csc_matrix(np.array([[2.,-1.,2.],[0.,2.,1.],[1.,1.,1.]]))
    with pytest.raises(ValueError,match='integer UMI'):
        qualify_matrix(a,['ENSG1','ENSG2','ENSG3'],['a','b','c'],['a','b','c'],n_cells=3,n_genes=1)

def test_non_ensembl_labels_rejected():
    with pytest.raises(ValueError,match='Ensembl'):
        qualify_matrix(csc_matrix(np.ones((4,6))),['a','b','c','d'],['a']*4,[f'cell{i}' for i in range(6)],n_cells=6,n_genes=2)
