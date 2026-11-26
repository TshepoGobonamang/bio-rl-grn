import copy
import pytest
from scripts.execute_encode_real_responses import selected_files


def example():
    accessions=['ENCSR095PIC','ENCSR552HTC','ENCSR620DOP',
                'ENCSR688UTN','ENCSR052ZPX','ENCSR650QAK','ENCSR808HQN']
    tfs=[None,'TFDP1','KLF2','NR2F2','RELA','SP1','NR2C2']
    panel={'protocol_id':'test','genome_assembly':'GRCh38','genome_annotation':'V29',
           'control':{'accession':accessions[0]},
           'perturbations':[{'accession':a,'tf':t} for a,t in zip(accessions[1:],tfs[1:])]}
    experiments=[]
    for j,(a,t) in enumerate(zip(accessions,tfs)):
        reps=[{'biological_replicate_number':k,'biosample_id':f'ENCBS{j}TEST{k}'}
              for k in (1,2)]
        files=[]
        for k in (1,2):
            ident=f'ENCFF{j*2+k:06d}'
            files.append({'accession':ident,'assembly':'GRCh38','genome_annotation':'V29',
                'biological_replicates':[k],
                'md5sum':'a'*32,'file_size':100,
                'href':f'/files/{ident}/@@download/{ident}.tsv'})
        experiments.append({'accession':a,'expected_tf':t,
              'source_quantification_pass':True,'source_json_sha256':'b'*64,
              'replicate_to_biosample':reps,
              'per_replicate_grch38_quantifications':files})
    audit={'data_origin':'ORIGINAL_ENCODE_API',
           'panel_per_replicate_quantification_pass':True,'errors':{},
           'retrieved_experiments':7,'experiments':experiments}
    return audit,panel


def test_frozen_v29_exact_biosample_file_selection():
    a,p=example()
    r=selected_files(a,p)
    assert len(r)==7
    assert sum(len(v) for v in r.values())==14
    assert [v['biosample_id'] for v in r['ENCSR095PIC']]==['ENCBS0TEST1','ENCBS0TEST2']


def test_unqualified_source_or_test_fixture_rejected():
    a,p=example()
    a['data_origin']='TEST_FIXTURES_NOT_REAL'
    with pytest.raises(ValueError,match='original ENCODE API'):
        selected_files(a,p)
    a,p=example()
    a['panel_per_replicate_quantification_pass']=False
    with pytest.raises(ValueError,match='not qualified'):
        selected_files(a,p)


def test_fails_if_multiple_versions_or_duplicate_replicate_files():
    a,p=example()
    clone=copy.deepcopy(a['experiments'][0]['per_replicate_grch38_quantifications'][0])
    clone['accession']='ENCFF999999'
    clone['href']='/files/ENCFF999999/@@download/ENCFF999999.tsv'
    a['experiments'][0]['per_replicate_grch38_quantifications'].append(clone)
    with pytest.raises(ValueError,match='multiple V29'):
        selected_files(a,p)


def test_unverified_url_and_md5_block_download():
    a,p=example()
    a['experiments'][0]['per_replicate_grch38_quantifications'][0]['href']='https://other.example/file'
    with pytest.raises(ValueError,match='unsafe'):
        selected_files(a,p)
    a,p=example()
    a['experiments'][0]['per_replicate_grch38_quantifications'][0]['md5sum']=None
    with pytest.raises(ValueError,match='MD5'):
        selected_files(a,p)
