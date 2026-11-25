import importlib.util
from pathlib import Path
import pandas as pd
import pytest

script=Path(__file__).resolve().parents[1] / 'scripts/prepare_encode_replicate_responses.py'
spec=importlib.util.spec_from_file_location('responses',script)
responses=importlib.util.module_from_spec(spec)
spec.loader.exec_module(responses)


def fixture(tmp_path):
    config={'source_metadata_audit_passed':True,'tf':'TFDP1',
            'annotation':'V29','abundance_column':'TPM',
            'treatment_accession':'ENCSR552HTC','control_accession':'ENCSR095PIC','samples':[]}
    for arm, acc, data in [('perturbed','ENCSR552HTC',[1,3]),('control','ENCSR095PIC',[6,6])]:
        for i in range(2):
            path=tmp_path / f'{arm}{i}.tsv'
            pd.DataFrame({'gene_id':['ENSG1','ENSG2'],'TPM':[data[i],2]}).to_csv(path,sep='\t',index=False)
            config['samples'].append({'arm':arm,'experiment_accession':acc,
                         'biosample_id':f'{arm}-{i}','replicate_number':i+1,
                         'file_accession':f'file-{arm}{i}','annotation':'V29','path':str(path)})
    return config


def test_quantification_only_not_ground_truth(tmp_path):
    result,meta=responses.compute_response(fixture(tmp_path))
    assert len(result)==2
    assert result.loc[result.gene_id=='ENSG1','log2_mean_abundance_ratio'].iloc[0]<0
    assert not meta['regulatory_edge_truth_labels_generated']
    assert meta['n_biological_biosamples_per_arm']=={'perturbed':2,'control':2}


def test_reject_pseudoreplication(tmp_path):
    cfg=fixture(tmp_path)
    cfg['samples'][1]['biosample_id']=cfg['samples'][0]['biosample_id']
    with pytest.raises(ValueError,match='Duplicate biosamples'):
        responses.compute_response(cfg)


def test_reject_bad_control(tmp_path):
    cfg=fixture(tmp_path)
    cfg['samples'][-1]['experiment_accession']='ENCSRBAD'
    with pytest.raises(ValueError,match='does not match'):
        responses.compute_response(cfg)


def test_reject_unverified_source(tmp_path):
    cfg=fixture(tmp_path)
    cfg['source_metadata_audit_passed']=False
    with pytest.raises(ValueError,match='audit required'):
        responses.compute_response(cfg)


def test_reject_incorrect_annotation(tmp_path):
    cfg=fixture(tmp_path)
    cfg['samples'][0]['annotation']='V30'
    with pytest.raises(ValueError,match='annotation mismatch'):
        responses.compute_response(cfg)
