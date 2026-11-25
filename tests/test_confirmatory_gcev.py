"""Tests use development seeds; never preregistered confirmation seeds."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def load_runner():
    module_file = Path(__file__).parents[1] / 'scripts' / 'confirmatory_gcev.py'
    spec = importlib.util.spec_from_file_location('confirmatory_gcev', module_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spec_for(family='linear', regime='target_switch'):
    return {'network_id': 'DEV_11','seed': 11,'family': family,
            'evidence_regime': regime, 'n_genes': 7, 'n_cells': 140}


@pytest.mark.parametrize('family',['linear','nonlinear'])
def test_generation_has_valid_oracle_dag_and_data_only_for_development(tmp_path,family):
    mod = load_runner()
    root = tmp_path/'data'
    meta = mod.generate_network(root,spec_for(family))
    ref = pd.read_csv(root/'reference_complete.csv')
    x = pd.read_csv(root/'expression.csv',index_col=0)
    assert (x.shape,meta['n_truth_edges']) == ((140,7),len(ref))
    assert len(ref) >= 2
    assert all(int(r.source[1:]) < int(r.target[1:]) for r in ref.itertuples())
    assert np.isfinite(x.to_numpy(float)).all()


@pytest.mark.parametrize('regime',['target_switch','global_stable'])
def test_qualifying_synthetic_source_is_explicitly_oracle_derived(tmp_path,regime):
    mod=load_runner()
    s=spec_for(regime=regime)
    mod.generate_network(tmp_path,s)
    summaries=mod.generate_synthetic_evidence(tmp_path,s)
    for k in ('SIM_A','SIM_B'):
        p=Path(summaries[k]['path'])
        a=pd.read_csv(p)
        assert a.score.between(0,1).all()
        assert not a[['source','target']].duplicated().any()
        assert a.shape[0] >= 10
        assert summaries[k]['sha256'] == mod.sha256(p)


def test_primary_bootstrap_resamples_networks_not_edges():
    mod=load_runner()
    frames=[]
    for family in ('linear','nonlinear'):
        for size in (7,9):
            for regime in ('target_switch','global_stable'):
                for rep in range(6):
                    frames.append({'family':family,'n_genes':size,'evidence_regime':regime,
                                   'network_id':f'{family}_{size}_{regime}_{rep}',
                                   'diff':0.1 + rep*0.01})
    data=pd.DataFrame(frames)
    summary=mod.stratified_bootstrap(data,'diff',1000,123)
    assert summary['n_networks']==48
    assert summary['ci95_low']>0
    assert summary['ci95_high']<0.2
    with pytest.raises(ValueError,match='eight'):
        mod.stratified_bootstrap(data.iloc[:46],'diff',100,123)


def test_locked_protocol_seed_and_scenario_schedule(tmp_path):
    mod=load_runner()
    p=Path(__file__).parents[1]/'configs'/'GCEV_BENCHMARK_PREREG_2026_11_25.json'
    proto=json.loads(p.read_text())
    assert mod.validate_protocol(proto)
    assert len(proto['network_design'])==48
    from copy import deepcopy
    wrong=deepcopy(proto);wrong['network_design'][0]['seed']=42
    with pytest.raises(ValueError,match='Seed schedule'):
        mod.validate_protocol(wrong)


def test_real_pipeline_evaluation_development_seed_only(tmp_path):
    mod=load_runner()
    example=mod.one_network(spec_for(),[{'id':'clean','missing':0.,'corrupt':0.}],str(tmp_path),'DEV_NOT_PREREG')
    assert example['success'],example.get('error')
    t=pd.read_csv(example['metrics_path'])
    assert {'gcev','no_graph_state','fixed_sources','global_source','conditional','greedy'}.issubset(t.method)
    assert not t.duplicated(['method','network_id','scenario']).any()
    fitjson=json.loads((tmp_path/'networks'/'DEV_11'/'conditions'/'clean'/'fit'/'fit_manifest.json').read_text())
    assert fitjson['gold_reference_accessed'] is False


def test_deterministic_source_generation(tmp_path):
    mod=load_runner();s=spec_for(family='nonlinear',regime='global_stable')
    x=tmp_path/'x';y=tmp_path/'y'
    for p in (x,y):
        mod.generate_network(p,s)
        mod.generate_synthetic_evidence(p,s)
    for name in ('expression.csv','reference_complete.csv','SYNTHETIC_ONLY_SIM_A_truth_derived.csv','SYNTHETIC_ONLY_SIM_B_truth_derived.csv'):
        assert mod.sha256(x/name)==mod.sha256(y/name)
