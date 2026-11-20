import json
import numpy as np
import pandas as pd
import pytest
from grnbench.core import (
    generate_smoke, partition, candidate_pairs, read_prior, evidence_sensitivity,
    run_fit, evaluate, sha256,
)


def test_partition_is_disjoint():
    p=partition(100,42)
    assert len(p.fit)==60 and len(p.score)==20 and len(p.test)==20
    assert len(set(p.fit)&set(p.score))==0
    assert len(set(p.test)&set(p.score))==0
    assert len(set(np.concatenate([p.fit,p.score,p.test])))==100


def test_prior_absence_is_masked(tmp_path):
    edges=candidate_pairs(['A','B','C'],['A','B'])
    prior=tmp_path/'prior.csv'
    pd.DataFrame([{'source':'A','target':'C','score':0.7}]).to_csv(prior,index=False)
    scores,mask=read_prior(prior,edges)
    assert int(mask.sum())==1
    assert scores[edges.index(('A','C'))]==0.7
    scores2,mask2=evidence_sensitivity(scores,mask,1,0,42)
    assert mask2.sum()==0
    assert np.isfinite(scores2).all()


def test_fit_is_gold_blind_and_evaluation_is_separate(tmp_path):
    data=tmp_path/'raw'
    generate_smoke(data,seed=7,n_cells=90,n_genes=6)
    run=tmp_path/'run'
    manifest=run_fit(data/'expression.csv',run,tf_list=data/'tfs.csv',
                     prior=data/'TOY_ONLY_random_prior.csv',max_genes=6,k=4,
                     graph_query_budget=1000,seed=7)
    assert manifest['gold_reference_accessed'] is False
    assert manifest['fit_cells']==54
    assert len(manifest['ranking_sha256'])==7
    assert (run/'greedy_trajectory.csv').exists()
    with pytest.raises(ValueError,match='complete benchmark'):
        evaluate(run,data/'reference_complete.csv',tmp_path/'eval')
    m=evaluate(run,data/'reference_complete.csv',tmp_path/'eval',
               reference_complete=True,k=4)
    assert set(m.method)=={'pearson','spearman','conditional','fixed','heuristic','adaptive','greedy'}
    assert all((m.ap>=0)&(m.ap<=1))
    assert all(np.isfinite(m.test_nmse))


def test_any_rank_mutation_blocks_gold_evaluation(tmp_path):
    data=tmp_path/'raw'
    generate_smoke(data,seed=2,n_cells=70,n_genes=5)
    fit=tmp_path/'fit'
    run_fit(data/'expression.csv',fit,tf_list=data/'tfs.csv',
            max_genes=5,k=3,graph_query_budget=300)
    f=fit/'rankings'/'pearson.csv'
    f.write_text(f.read_text()+'\n')
    with pytest.raises(ValueError,match='mutated'):
        evaluate(fit,data/'reference_complete.csv',tmp_path/'eval',reference_complete=True)


def test_all_gene_mode_must_be_explicit(tmp_path):
    data=tmp_path/'data'
    generate_smoke(data,n_cells=40,n_genes=4)
    with pytest.raises(ValueError,match='TF whitelist'):
        run_fit(data/'expression.csv',tmp_path/'fit',max_genes=4,k=2)


def test_no_prior_means_no_adaptive_claim(tmp_path):
    data=tmp_path/'data'
    generate_smoke(data,n_cells=60,n_genes=4)
    out=tmp_path/'fit'
    run_fit(data/'expression.csv',out,tf_list=data/'tfs.csv',max_genes=4,k=2)
    assert set(json.loads((out/'fit_manifest.json').read_text())['ranking_sha256'])=={'pearson','spearman','conditional','greedy'}
