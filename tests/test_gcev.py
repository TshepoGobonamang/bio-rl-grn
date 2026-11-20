"""GCEV acceptance tests: paired decisions, genuine target hold-out and no gold."""
from __future__ import annotations
import json

import numpy as np
import pandas as pd
import pytest

from grnbench.core import candidate_pairs, generate_smoke, run_fit, evaluate
from grnbench.gcev import (paired_source_valuation, learn_graph_conditional_weights,
                           compose_graph_conditional_rank, GRAPH_FEATURES)


def _fixture():
    rng=np.random.default_rng(33)
    x=pd.DataFrame(rng.normal(size=(115,6)),columns=list('ABCDEF'))
    x['B']=1.5*x['A']+rng.normal(scale=.6,size=len(x))
    x['C']=x['A']-x['D']+rng.normal(scale=.7,size=len(x))
    edges=candidate_pairs(list(x),list(x))
    src_a=rng.random(len(edges));src_b=rng.random(len(edges))
    full=np.ones(len(edges),bool)
    return x,edges,{'A':(src_a,full),'B':(src_b,full)}


def test_paired_counterfactual_val_is_not_a_rank_stability_proxy():
    x,edges,src=_fixture()
    labels=paired_source_valuation(x,edges,src,boost=.8,max_parents=1,seed=42)
    assert len(labels)==3*len(x.columns)*len(src)
    assert labels['different_edge'].sum()>0, 'Controlled fixture must change graph actions'
    assert (labels.loc[labels.different_edge.eq(0),'source_value']==0).all()
    assert labels['source_value'].between(-1,1).all()
    assert labels['used_evaluation_reference'].eq(False).all()
    assert labels['equal_edit_budget'].eq(True).all()
    assert (labels['graph_parents'].str.len()>0).all()


def test_target_oof_predictions_vary_with_graph_features():
    x,edges,src=_fixture()
    labels=paired_source_valuation(x,edges,src,boost=.8,max_parents=1,seed=42)
    fit=learn_graph_conditional_weights(labels,shrinkage=.05)
    assert fit.gate_fit_metadata['oof_by_target']
    assert set(['gcev_weight','no_graph_state_weight','global_source_weight']).issubset(fit.per_target_source)
    assert fit.per_target_source.gcev_weight.between(0,1).all()
    assert set(GRAPH_FEATURES).issubset(fit.gate_fit_metadata['features_graph_conditional'])
    # OOF estimates are predictively different from both simple-source and
    # graph-ablated controls on this deliberately heterogeneous fixture.
    assert not np.allclose(fit.per_target_source.gcev_weight,
                           fit.per_target_source.no_graph_state_weight)
    score=compose_graph_conditional_rank(np.ones(len(edges))*.5,edges,src,
                                        fit.per_target_source,'gcev',boost=.3)
    assert np.isfinite(score).all()


def test_unobserved_source_is_missing_not_zero_as_negative():
    x,edges,src=_fixture()
    no=np.zeros(len(edges),bool)
    empty={'A':(np.ones(len(edges)),no)}
    labels=paired_source_valuation(x,edges,empty,seed=42)
    assert labels.different_edge.sum()==0
    assert labels.source_value.eq(0).all()
    fit=learn_graph_conditional_weights(labels)
    score=compose_graph_conditional_rank(np.ones(len(edges))*.5,edges,empty,
                                        fit.per_target_source,'gcev')
    assert np.allclose(score,.5)


def test_no_reference_entrypoint_and_separate_evaluation(tmp_path):
    d=tmp_path/'raw';generate_smoke(d,seed=29,n_cells=110,n_genes=6)
    a=tmp_path/'run'
    manifest=run_fit(d/'expression.csv',a,tf_list=d/'tfs.csv',
                     prior=d/'TOY_ONLY_random_prior.csv',gcev=True,max_genes=6,
                     k=4,seed=29,graph_query_budget=1000)
    assert manifest['gold_reference_accessed'] is False
    assert manifest['gcev_enabled']
    assert (a/'gcev'/'training_contract.json').exists()
    assert (a/'gcev'/'paired_source_value_targets.csv').exists()
    weights=pd.read_csv(a/'gcev'/'oof_source_weights.csv')
    assert weights.gcev_weight.between(0,1).all()
    assert json.loads((a/'gcev'/'training_contract.json').read_text())['no_gold_used']
    ev=evaluate(a,d/'reference_complete.csv',tmp_path/'evaluation',reference_complete=True)
    assert {'fixed_sources','global_source','no_graph_state','gcev'}.issubset(ev.method)
    assert 'gcev' in manifest['ranking_sha256']


def test_gcev_cannot_proceed_without_evidence_or_sufficient_cells(tmp_path):
    d=tmp_path/'raw';generate_smoke(d,seed=31,n_cells=65,n_genes=4)
    with pytest.raises(ValueError,match='provided'):
        run_fit(d/'expression.csv',tmp_path/'fail',tf_list=d/'tfs.csv',
                gcev=True,max_genes=4)
    x,edges,src=_fixture()
    with pytest.raises(ValueError,match='>=48'):
        paired_source_valuation(x.iloc[:40],edges,src)


def test_gcev_crossfit_is_independent_of_internal_and_final_test_cells(tmp_path):
    from grnbench.core import partition, read_expression
    data=tmp_path/'raw';generate_smoke(data,seed=88,n_cells=110,n_genes=6)
    expression=data/'expression.csv'
    orig=read_expression(expression)
    altered=orig.copy()
    split=partition(len(orig),seed=88)
    altered.iloc[np.concatenate((split.score,split.test)),:] += 40.0
    changed=data/'expression_altered.csv'
    altered.to_csv(changed,index_label='cell_id')
    orig_dir,alter_dir=tmp_path/'orig',tmp_path/'changed'
    opts=dict(tf_list=data/'tfs.csv',prior=data/'TOY_ONLY_random_prior.csv',
              gcev=True,max_genes=6,k=3,graph_query_budget=750,seed=88)
    run_fit(expression,orig_dir,**opts)
    run_fit(changed,alter_dir,**opts)
    for name in ('paired_source_value_targets.csv','oof_source_weights.csv'):
        one=pd.read_csv(orig_dir/'gcev'/name)
        two=pd.read_csv(alter_dir/'gcev'/name)
        pd.testing.assert_frame_equal(one,two,check_exact=False,atol=1e-10,rtol=1e-10)


def test_tampered_training_value_blocks_reference_evaluation(tmp_path):
    d=tmp_path/'raw';generate_smoke(d,seed=81,n_cells=90,n_genes=5)
    fit=tmp_path/'fit'
    run_fit(d/'expression.csv',fit,tf_list=d/'tfs.csv',
            prior=d/'TOY_ONLY_random_prior.csv',gcev=True,max_genes=5,k=4,
            seed=81,graph_query_budget=750)
    target_file=fit/'gcev'/'paired_source_value_targets.csv'
    target_file.write_text(target_file.read_text()+'\n')
    with pytest.raises(ValueError,match='GCEV artifact mutated'):
        evaluate(fit,d/'reference_complete.csv',tmp_path/'ev',reference_complete=True)
