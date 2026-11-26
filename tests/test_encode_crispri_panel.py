from pathlib import Path
import importlib.util

path = Path(__file__).resolve().parents[1] / 'scripts/audit_encode_crispri_panel.py'
spec = importlib.util.spec_from_file_location('encode_audit', path)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def exp(accession, tf, control='ENCSR095PIC', reps=(1,2), bios=('B1','B2'), quant=True):
    obj = {'accession':accession,'status':'released','assay_title':'CRISPRi RNA-seq',
           'biosample_summary':'Homo sapiens K562','replication_type':'isogenic',
           'replicates':[{'biological_replicate_number':r,'biosample':b}
                         for r,b in zip(reps,bios)],'files':[], 'audit':{}}
    if tf is None:
        obj['control_type'] = 'non-specific target'
    else:
        obj.update({'target':f'/targets/{tf}-human/',
                    'possible_controls':[f'/experiments/{control}/']})
    if quant:
        obj['files']=[{'accession': f'ENCFF{i:06d}', 'status':'released',
                       'file_format':'tsv','output_type':'gene quantifications',
                       'assembly':'GRCh38','genome_annotation':'V29',
                       'biological_replicates':[r]} for i,r in enumerate(reps)]
    return obj


def test_qualified_pair_and_control():
    m={'protocol_id':'test', 'control':{'accession':'ENCSR095PIC'},
       'perturbations':[{'tf':'TFDP1','accession':'ENCSR552HTC'}]}
    docs={'ENCSR095PIC':exp('ENCSR095PIC',None),
          'ENCSR552HTC':exp('ENCSR552HTC','TFDP1')}
    r=audit.evaluate_panel(m,lambda s:(docs[s], 'f'*64))
    assert r['panel_design_metadata_pass']
    assert r['panel_per_replicate_quantification_pass']
    assert r['decision']=='SOURCE_METADATA_QUALIFIED_FOR_NEXT_REVIEW'
    assert not r['independent_biological_validation_completed']


def test_shared_biosamples_rejected():
    x=audit.inspect_experiment(exp('ENCSR552HTC','TFDP1',bios=('B1','B1')),
                                expected_tf='TFDP1',expected_control='ENCSR095PIC', raw_sha256='a')
    assert not x['biosample_replicates_verified_from_metadata']


def test_merged_quantification_not_replication():
    x=exp('ENCSR552HTC','TFDP1')
    x['files']=[{'status':'released','file_format':'tsv','output_type':'gene quantifications',
                 'assembly':'GRCh38','genome_annotation':'V29','biological_replicates':[1,2]}]
    r=audit.inspect_experiment(x,expected_tf='TFDP1', expected_control='ENCSR095PIC',raw_sha256='a')
    assert r['biosample_replicates_verified_from_metadata']
    assert not r['per_replicate_quantification_available']
    assert not r['source_quantification_pass']


def test_unlinked_control_blocks():
    x=exp('ENCSR552HTC','TFDP1',control='ENCSRWRONG')
    r=audit.inspect_experiment(x,expected_tf='TFDP1',expected_control='ENCSR095PIC',raw_sha256='a')
    assert not r['design_metadata_pass']


def test_missing_replicate_info_fails_closed():
    x=exp('ENCSR552HTC','TFDP1')
    x['replicates']=['/replicates/id-1/','/replicates/id-2/']
    r=audit.inspect_experiment(x,expected_tf='TFDP1',expected_control='ENCSR095PIC',raw_sha256='a')
    assert not r['design_metadata_pass']
    assert r['unresolved_replicate_references']==2


def test_missing_accession_retained_as_error():
    m={'protocol_id':'test', 'control':{'accession':'ENCSR095PIC'},
       'perturbations':[{'tf':'TFDP1','accession':'ENCSR552HTC'}]}
    r=audit.evaluate_panel(m,lambda s: (exp(s,None),'f') if s=='ENCSR095PIC'
                           else (_ for _ in ()).throw(FileNotFoundError(s)))
    assert r['requested_experiments']==2
    assert r['retrieved_experiments']==1
    assert r['decision']=='NO_GO_FOR_BIOLOGICAL_SCORING'
    assert 'ENCSR552HTC' in r['errors']


def test_qc_flags_remain_visible():
    x=exp('ENCSR552HTC','TFDP1')
    x['audit']={'NOT_COMPLIANT':['example QC concern']}
    r=audit.inspect_experiment(x,expected_tf='TFDP1', expected_control='ENCSR095PIC',raw_sha256='a')
    assert r['manual_audit_review_required']


def test_original_linked_replicates_resolved_and_has_source_hashes():
    p1='/replicates/aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa/'
    p2='/replicates/bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb/'
    x=exp('ENCSR552HTC','TFDP1')
    x['replicates']=[p1,p2]
    records={
        p1:{'@id':p1,'experiment':'/experiments/ENCSR552HTC/',
            'biological_replicate_number':1,'biosample':'/biosamples/ENCBSAAA123/'},
        p2:{'@id':p2,'experiment':'/experiments/ENCSR552HTC/',
            'biological_replicate_number':2,'biosample':'/biosamples/ENCBSBBB234/'},
    }
    def fetch(url):
        path=url.replace(audit.HOST,'').split('?')[0]
        return records[path], 'e'*64
    resolved=audit.resolve_replicates(x,fetch)
    assert list(resolved['_gcev_original_replicate_sha256'])==[p1,p2]
    r=audit.inspect_experiment(resolved,expected_tf='TFDP1',
                               expected_control='ENCSR095PIC',raw_sha256='a')
    assert r['design_metadata_pass']
    assert r['source_quantification_pass']
    assert len(r['original_replicate_record_sha256'])==2


def test_original_replicate_wrong_experiment_rejected():
    p='/replicates/aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa/'
    x=exp('ENCSR552HTC','TFDP1')
    x['replicates']=[p]
    import pytest
    with pytest.raises(ValueError,match='different experiment'):
        audit.resolve_replicates(x,lambda url:(
            {'@id':p,'experiment':'/experiments/ENCSR620DOP/',
             'biological_replicate_number':1,'biosample':'/biosamples/ENCBSAAA123/'},
            'f'*64))


def test_missing_original_qc_annotations_requires_manual_review():
    x=exp('ENCSR552HTC','TFDP1')
    del x['audit']
    r=audit.inspect_experiment(x,expected_tf='TFDP1',
                               expected_control='ENCSR095PIC',raw_sha256='a')
    assert r['manual_audit_review_required']
    assert r['original_qc_metadata_missing']
