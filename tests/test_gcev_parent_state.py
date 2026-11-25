"""Isolated explicit-parent-feature block: matched training rows and action features."""
import numpy as np
import pandas as pd
from grnbench.gcev import (SOURCE_FEATURES, PARENT_STATE_FEATURES, ACTION_FEATURES,
  GRAPH_FEATURES, paired_source_valuation, learn_graph_conditional_weights)
from grnbench.core import candidate_pairs
from scripts.gcev_parent_state_confirm import SEEDS, LEVELS, SCENARIOS, bootstrap, design, summarize
import pytest


def _targets():
    rng=np.random.default_rng(91)
    x=pd.DataFrame(rng.normal(size=(125,6)),columns=list('ABCDEF'))
    x['C']=x['A']+0.8*x['B']+rng.normal(size=len(x))*.5
    edges=candidate_pairs(list(x),list(x))
    evidence={'A':(rng.random(len(edges)),np.ones(len(edges),bool)),
              'B':(rng.random(len(edges)),np.ones(len(edges),bool))}
    return paired_source_valuation(x,edges,evidence,alpha=1.,max_parents=2,boost=.6,seed=91)

def test_feature_partition_does_not_remove_action_context():
    assert PARENT_STATE_FEATURES==('parent_count','parent_target_correlation')
    assert set(PARENT_STATE_FEATURES).isdisjoint(ACTION_FEATURES)
    assert set(PARENT_STATE_FEATURES).isdisjoint(SOURCE_FEATURES)
    assert GRAPH_FEATURES==PARENT_STATE_FEATURES+ACTION_FEATURES

def test_only_gcev_depends_on_explicit_parent_columns():
    t=_targets(); original=learn_graph_conditional_weights(t)
    changed=t.copy()
    changed['parent_target_correlation'] = changed.parent_target_correlation * 17. + np.arange(len(changed))*.07
    changed['parent_count'] = changed.parent_count + (np.arange(len(changed))%3)
    fit=learn_graph_conditional_weights(changed)
    for name in ['no_parent_state','no_graph_state','global_source']:
        np.testing.assert_allclose(original.paired_decisions[name+'_oof_predicted_value'],
                                   fit.paired_decisions[name+'_oof_predicted_value'],atol=1e-12)
    assert not np.allclose(original.paired_decisions.gcev_oof_predicted_value,
                           fit.paired_decisions.gcev_oof_predicted_value)
    assert original.gate_fit_metadata['features_no_parent_state']==list(SOURCE_FEATURES+ACTION_FEATURES)

def test_all_methods_are_group_oof_and_bounded():
    fit=learn_graph_conditional_weights(_targets())
    for name in ['gcev','no_parent_state','no_action_state','no_graph_state','global_source']:
        assert fit.per_target_source[name+'_weight'].between(0,1).all()
    assert fit.gate_fit_metadata['oof_by_target'] is True

def test_paired_seed_grid_and_statistics():
    specs=[design(s) for s in SEEDS]
    assert len(specs)==32 and set(SEEDS).isdisjoint({11,42,77,*range(2001,2049),*range(6101,6133)})
    assert len(LEVELS)==4 and len(SCENARIOS)==2
    assert len(pd.DataFrame(specs).groupby(['n_genes','evidence_regime']).size())==4
    assert all(pd.DataFrame(specs).groupby(['n_genes','evidence_regime']).size()==8)
    with pytest.raises(ValueError): design(6101)
    with pytest.raises(ValueError, match='Expected 32'):bootstrap(np.zeros(31),np.array(['x']*31))

def test_incomplete_ablation_method_grid_rejected(tmp_path):
    with pytest.raises(ValueError,match='Incomplete'):
        summarize(pd.DataFrame({'method':['gcev'],'network_id':['x'],'scenario':['clean']}),pd.DataFrame(),{},tmp_path)
