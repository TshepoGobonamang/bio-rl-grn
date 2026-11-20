import json
import pandas as pd
import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from aggregate_runs import aggregate


def test_network_level_bootstrap(tmp_path):
    evals=[]
    for i in range(2):
        p=tmp_path/f'network_{i}'
        p.mkdir()
        (p/'evaluation_manifest.json').write_text(json.dumps({'source':'unit-test'}))
        pd.DataFrame([{'method':'pearson','ap':.2+i*.1,'auroc':.6,
                       'test_nmse':1+i*.1,'reference_prevalence':.1}]).to_csv(p/'metrics.csv',index=False)
        evals.append(p)
    report=aggregate(evals,tmp_path/'out',seed=42,draws=100)
    assert report.loc[report.metric=='ap','mean'].iloc[0]==pytest.approx(.25)
    assert (tmp_path/'out'/'per_network.csv').exists()
    assert set(report['n_independent_networks_reported'])=={2}


def test_unbalanced_methods_rejected(tmp_path):
    evals=[]
    for i,method in enumerate(('pearson','spearman')):
        p=tmp_path/f'net_{i}';p.mkdir()
        (p/'evaluation_manifest.json').write_text('{}')
        pd.DataFrame([{'method':method,'ap':.4,'auroc':.6,
                       'test_nmse':1,'reference_prevalence':.1}]).to_csv(p/'metrics.csv',index=False)
        evals.append(p)
    with pytest.raises(ValueError,match='incomparable'):
        aggregate(evals,tmp_path/'out')
