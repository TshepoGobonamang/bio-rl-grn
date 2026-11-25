import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scripts.gcev_mechanism_factorial import (
    SEEDS, LEVELS, SCENARIOS, bootstrap, design, draw_dag, expression_for,
    generate_pair, summarize,
)

def test_frozen_design_strata():
    specs=[design(s) for s in SEEDS]
    assert len(specs)==32 and set(SEEDS).isdisjoint({11,42,77,*range(2001,2049)})
    counts=pd.DataFrame(specs).groupby(['n_genes','evidence_regime']).size()
    assert len(counts)==4 and all(counts==8)
    assert len(LEVELS)==4 and len(SCENARIOS)==2
    with pytest.raises(ValueError): design(2001)

def test_identical_graph_noise_and_prior_across_factorial(tmp_path):
    spec={'network_id':'development_9001','seed':9001,'n_genes':7,'n_cells':140,'evidence_regime':'target_switch'}
    meta=generate_pair(tmp_path,spec)
    assert meta['shared_exogenous_noise'] is True
    assert meta['n_true_edges']>=2
    x=[pd.read_csv(tmp_path/'mechanisms'/m/'expression.csv').iloc[:,1:].to_numpy() for m,_,_ in LEVELS]
    assert all(np.array_equal(a[:,0],x[0][:,0]) for a in x)
    assert not np.array_equal(x[0],x[1])
    assert not np.array_equal(x[0],x[2]) if meta['multi_parent_targets']>0 else True
    assert set(meta['source_hashes'])=={'SIM_A','SIM_B'}
    assert len({meta['reference_sha256']})==1

def test_mechanism_terms_require_two_parents():
    w=np.zeros((3,3));w[0,1]=.4;w[1,2]=.5
    noise=np.arange(18).reshape(6,3).astype(float)/9
    a=expression_for(w,noise,0,0);b=expression_for(w,noise,0,1)
    assert np.array_equal(a,b)  # no target has two true parents
    c=expression_for(w,noise,1,0)
    assert not np.array_equal(a,c)

def test_bootstrap_uses_whole_dag_and_fails_on_wrong_design():
    strata=np.repeat(['a','b','c','d'],8)
    vals=np.tile(np.linspace(-.1,.15,8),4)
    r=bootstrap(vals,strata,draws=50,seed=10)
    assert r['n_dag']==32 and r['unit']=='independent DAG'
    with pytest.raises(ValueError): bootstrap(vals[:-1],strata[:-1])

def test_incomplete_method_grid_refused(tmp_path):
    with pytest.raises(ValueError,match='Incomplete'):
        summarize(pd.DataFrame({'method':['gcev'],'network_id':['a'],'mechanism':['linear'],'scenario':['clean']}),
                  pd.DataFrame(),{},tmp_path)
