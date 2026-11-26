import pandas as pd
import pytest

from scripts.verify_biov2_deseq2 import verify_one

def fixture():
    frame=pd.DataFrame({
        'gene_id':['ENSG00001.1','ENSG00002.1','ENSG00003.1'],
        'baseMean':[12.,5.,19.],
        'log2FoldChange':[1.5,-0.4,0.2],
        'lfcSE':[.3,.8,.4],
        'stat':[5.,-.5,.5],
        'pvalue':[0.0001,0.7,0.6],
        'padj':[0.0003,None,0.9],
        'bh_fdr_005':['TRUE','FALSE','FALSE']})
    row=pd.Series({'genes_input':3,'genes_retained':3,
        'tested_pvalues':3,'tested_padj':2,
        'fdr_005':1,'fdr_005_up':1,'fdr_005_down':0,
        'n_controls':2,'n_perturbed':2,
        'is_direct_TF_binding_truth':False,
        'is_independent_validation_of_GCEV':False})
    return frame,row

def test_original_deseq2_saved_count_verification():
    frame,row=fixture()
    result=verify_one('TFDP1',frame,row,expected_input=3)
    assert result['fdr_005']==1
    assert not result['direct_regulatory_truth_verified']
    assert not result['gcev_biological_accuracy_test_performed']

def test_reject_summary_inconsistent_with_gene_values():
    frame,row=fixture()
    row['fdr_005']=2
    with pytest.raises(ValueError,match='summary fdr_005 mismatch'):
        verify_one('TFDP1',frame,row,3)

def test_reject_false_significance_labels():
    frame,row=fixture()
    frame.loc[1,'bh_fdr_005']='TRUE'
    with pytest.raises(ValueError,match='gene-wise FDR calls'):
        verify_one('TFDP1',frame,row,3)

def test_reject_duplicate_genes():
    frame,row=fixture()
    frame.loc[1,'gene_id']=frame.loc[0,'gene_id']
    with pytest.raises(ValueError,match='duplicate'):
        verify_one('TFDP1',frame,row,3)

def test_reject_deseq2_called_direct_binding_or_independent_biovalidation():
    frame,row=fixture()
    row['is_direct_TF_binding_truth']=True
    with pytest.raises(ValueError,match='false direct binding'):
        verify_one('TFDP1',frame,row,3)
    frame,row=fixture()
    row['is_independent_validation_of_GCEV']=True
    with pytest.raises(ValueError,match='false GCEV'):
        verify_one('TFDP1',frame,row,3)

def test_reject_using_technical_not_biological_replication():
    frame,row=fixture()
    row['n_controls']=1
    with pytest.raises(ValueError,match='insufficient original'):
        verify_one('TFDP1',frame,row,3)
