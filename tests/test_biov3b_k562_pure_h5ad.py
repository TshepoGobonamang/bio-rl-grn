import numpy as np
import pytest
from scipy import sparse
from scripts.audit_biov3b_k562_pure_h5ad import inspect_array

def test_raw_sparse_counts_preserve_source_unchanged():
    x=sparse.csr_matrix(np.array([[0,1,3],[2,0,0]],dtype=float))
    r=inspect_array(x,["cell1","cell2"],["ENSG1","ENSG2","ENSG3"],"0"*64)
    assert r["n_cells"]==2 and r["n_genes"]==3
    assert r["matrix_is_sparse"]
    assert r["fraction_nonzero_entries_near_integer"]==1
    assert not r["gcev_training_performed"]
    assert not r["cell_type_from_numeric_mixture_inferred"]

def test_normalized_float_scale_is_not_falsely_called_counts():
    x=sparse.csr_matrix(np.array([[0,.23,2.65],[1.44,0,0]]))
    r=inspect_array(x,["c1","c2"],["G1","G2","G3"],"1"*64)
    assert r["matrix_scale_classification"]=="NONINTEGER_EXPRESSION_SCALE_NEEDS_SOURCE_CONFIRMATION"

@pytest.mark.parametrize("x",[
    np.array([[np.nan,1],[2,3]]),np.array([[np.inf,1],[2,3]]),
    np.array([[-1.,1.],[2.,3.]])])
def test_invalid_original_expression_is_blocked(x):
    with pytest.raises(ValueError,match="nonfinite or negative"):
        inspect_array(x,["c1","c2"],["g1","g2"],"2"*64)

def test_duplicate_cell_identity_rejected():
    with pytest.raises(ValueError,match="Duplicate"):
        inspect_array(np.array([[1,2],[3,4]]),["x","x"],["a","b"],"3"*64)

def test_inconsistent_shape_rejected():
    with pytest.raises(ValueError,match="shape"):
        inspect_array(np.array([[1,2]]),["c1","c2"],["g1","g2"],"3"*64)
