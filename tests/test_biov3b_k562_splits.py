import numpy as np
from scipy import sparse
from scripts.prepare_biov3b_k562_splits import prepare

def example(seed=123):
    x=sparse.csr_matrix(np.random.default_rng(seed).poisson(1,size=(120,32)))
    return x,[f"c{i}" for i in range(120)],[f"G{i}" for i in range(32)]

def test_partition_and_fit_only_source_features():
    x,c,g=example()
    matrix,part,genes,manifest=prepare(x,c,g,n_features=8)
    assert matrix.shape==(120,8)
    assert manifest["fit_cells"]==72 and manifest["internal_score_cells"]==24
    assert manifest["untouched_test_cells"]==24
    assert set(part["partition"])=={"fit","internal_score","untouched_test"}
    assert len(set(genes["gene_symbol"]))==8
    assert manifest["no_perturbation_response_or_gold_accessed"]

def test_untouched_test_cells_cannot_change_feature_selection():
    x,c,g=example()
    _,part,g1,_=prepare(x,c,g,n_features=8)
    holdout=np.flatnonzero(part["partition"]=="untouched_test")
    xx=x.tolil()
    for i in holdout:
        xx[i,:]=np.arange(32)*1000
    _,part2,g2,_=prepare(xx.tocsr(),c,g,n_features=8)
    assert part["partition"].tolist()==part2["partition"].tolist()
    assert g1["gene_symbol"].tolist()==g2["gene_symbol"].tolist()
