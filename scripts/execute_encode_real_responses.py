"""Run original ENCODE K562 per-biosample V29 expression response summaries.

Consumes only the authenticated live source audit, never fixtures. Downloads the
17 frozen original-file quantifications with MD5 validation and does not publish
raw original expression TSVs or generate direct TF-binding labels.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

from scripts.prepare_encode_replicate_responses import compute_response, sha256

HOST = 'https://www.encodeproject.org'
FILE_RE = re.compile(r'^/files/(ENCFF[0-9A-Z]{6})/@@download/\1\.tsv$')


def selected_files(audit: dict, panel: dict) -> dict:
    if audit.get('data_origin') != 'ORIGINAL_ENCODE_API':
        raise ValueError('Actual original ENCODE API report required; fixtures are prohibited')
    if not audit.get('panel_per_replicate_quantification_pass'):
        raise ValueError('ENCODE panel source and per-replicate metadata not qualified')
    if audit.get('errors') or audit.get('retrieved_experiments') != 7:
        raise ValueError('Missing/failed original source record')
    if panel.get('genome_assembly') != 'GRCh38' or panel.get('genome_annotation') != 'V29':
        raise ValueError('Frozen GRCh38 V29 protocol mismatch')
    expected = {panel['control']['accession']: None}
    expected.update({item['accession']: item['tf'] for item in panel['perturbations']})
    original = {x['accession']: x for x in audit['experiments']}
    if len(expected) != 7 or set(original) != set(expected):
        raise ValueError('Unexpected ENCODE accession set')
    selection = {}
    for accession, tf in expected.items():
        obj = original[accession]
        if obj.get('expected_tf') != tf or not obj.get('source_quantification_pass'):
            raise ValueError(f'{accession}: original source qualification mismatch')
        if not obj.get('source_json_sha256') or len(obj['source_json_sha256']) != 64:
            raise ValueError(f'{accession}: source hash missing')
        samples = obj.get('replicate_to_biosample') or []
        if len(samples) < 2 or len({x['biosample_id'] for x in samples}) != len(samples):
            raise ValueError(f'{accession}: independent biosample records missing')
        ids = {int(x['biological_replicate_number']): str(x['biosample_id']) for x in samples}
        if len(ids) != len(samples):
            raise ValueError(f'{accession}: repeated replicate number')
        q = obj.get('per_replicate_grch38_quantifications') or []
        byrep = {}
        for file in q:
            if file.get('genome_annotation') != 'V29' or file.get('assembly') != 'GRCh38':
                raise ValueError('Mixed annotation/assembly')
            nums = file.get('biological_replicates') or []
            if len(nums) != 1 or nums[0] not in ids:
                raise ValueError('File-to-biosample assignment not verifiable')
            n = int(nums[0])
            if n in byrep:
                raise ValueError(f'{accession}: multiple V29 quantification files for replicate {n}')
            acc = file['accession']
            href = file.get('href') or ''
            m = FILE_RE.fullmatch(href)
            if not m or m.group(1) != acc:
                raise ValueError(f'{accession}: unsafe/unverified ENCODE file location')
            md5 = file.get('md5sum')
            if not isinstance(md5, str) or not re.fullmatch(r'[0-9a-fA-F]{32}',md5):
                raise ValueError(f'{acc}: missing authentic MD5 checksum')
            size = file.get('file_size')
            if not isinstance(size, int) or not 0 < size <= 60_000_000:
                raise ValueError(f'{acc}: unsafe or unknown file size')
            byrep[n] = {'experiment_accession':accession,'replicate_number':n,
                        'biosample_id':ids[n],'file_accession':acc,
                        'annotation':'V29','href':href,
                        'md5sum':md5.lower(),'file_size':size}
        if set(byrep) != set(ids):
            raise ValueError(f'{accession}: incomplete V29 replicate quantification coverage')
        selection[accession]=[byrep[n] for n in sorted(byrep)]
    all_biosamples=[s['biosample_id'] for files in selection.values() for s in files]
    if len(set(all_biosamples)) != len(all_biosamples):
        raise ValueError('A biosample is reused across treatment and control experiments')
    return selection


def download_checked(file: dict, dest: Path, max_retries: int = 3) -> dict:
    url = HOST + file['href']
    err = None
    for attempt in range(max_retries):
        try:
            request = urllib.request.Request(
                url, headers={'User-Agent':'GCEV-ENCODE-source-reproducibility/0.1',
                              'Accept':'text/tab-separated-values,application/octet-stream'})
            with urllib.request.urlopen(request,timeout=120) as response, dest.open('wb') as output:
                h = hashlib.md5()
                n = 0
                for chunk in iter(lambda:response.read(1024*1024),b''):
                    n += len(chunk)
                    if n > 60_000_000:
                        raise ValueError('Downloaded gene file exceeds safe size limit')
                    h.update(chunk)
                    output.write(chunk)
            if n != file['file_size'] or h.hexdigest() != file['md5sum']:
                raise ValueError(f"Source content size/MD5 differs from ENCODE file {file['file_accession']}")
            return {'accession':file['file_accession'],'url':url,
                    'md5sum':h.hexdigest(),'sha256':sha256(dest),
                    'size_bytes':n,'experiment':file['experiment_accession'],
                    'biosample_id':file['biosample_id'],
                    'biological_replicate_number':file['replicate_number']}
        except (OSError, ValueError, TimeoutError) as exc:
            err = exc
            dest.unlink(missing_ok=True)
            if isinstance(exc, ValueError):
                raise
    raise RuntimeError(f"Could not obtain original file {file['file_accession']}: {err}")


def execute(audit: dict, panel: dict, output: Path) -> dict:
    selection = selected_files(audit,panel)
    output.mkdir(parents=True,exist_ok=True)
    control=panel['control']['accession']
    provenance=[]
    with tempfile.TemporaryDirectory(prefix='gcev_encode_v29_') as td:
        root=Path(td)
        originals={}
        for accession, items in selection.items():
            for file in items:
                dest=root/(file['file_accession']+'.tsv')
                provenance.append(download_checked(file,dest))
                originals[file['file_accession']]=dest
                print('VERIFIED', file['file_accession'], 'biosample',file['biosample_id'], flush=True)
        per_tf=[]
        for t in panel['perturbations']:
            accession=t['accession']
            rows=[]
            for arm,acc in [('control',control),('perturbed',accession)]:
                for f in selection[acc]:
                    rows.append({k:f[k] for k in ('experiment_accession','replicate_number',
                                                 'biosample_id','file_accession','annotation') } |
                                {'arm':arm,'path':str(originals[f['file_accession']])})
            c={'source_metadata_audit_passed':True,'annotation':'V29',
               'tf':t['tf'],'treatment_accession':accession,
               'control_accession':control,'gene_column':'gene_id',
               'abundance_column':'TPM','samples':rows}
            response, metadata=compute_response(c)
            f=output/(t['tf']+'_response_descriptive.csv.gz')
            response.to_csv(f,index=False,compression='gzip')
            metadata['exploratory_only']=True
            metadata['original_encode_qc_manually_reviewed']=False
            metadata['tf_knockdown_efficacy_verified']=False
            metadata['model_rankings_compared']=False
            metadata['source_disjointness_verified']=False
            metadata['response_sha256']=sha256(f)
            (output/(t['tf']+'_response_metadata.json')).write_text(json.dumps(metadata,indent=2)+'\n')
            per_tf.append({'tf':t['tf'],'treatment':accession,
                           'n_common_genes':len(response),'response_file':f.name,
                           'response_sha256':metadata['response_sha256'],
                           'mean_abs_log2_response':float(np.abs(response['log2_mean_abundance_ratio']).mean()),
                           'n_controls':metadata['n_biological_biosamples_per_arm']['control'],
                           'n_perturbed':metadata['n_biological_biosamples_per_arm']['perturbed']})
            print('DESCRIPTIVE',t['tf'],'genes',len(response), flush=True)
    report={'protocol_id':panel['protocol_id'],
            'data_origin':'ORIGINAL_ENCODE_REPLICATE_TSV',
            'annotation':'GRCh38/GENCODE V29',
            'source_qc_review_pending':True,
            'all_six_predeclared_tfs_present':len(per_tf)==6,
            'per_tf':per_tf,'original_file_checksums':provenance,
            'shared_control_accession':control,
            'perturbation_response_labels_are_direct_binding_truth':False,
            'model_prediction_rankings_evaluated':False,
            'knockdown_efficacy_verified':False,
            'stage':'EXPLORATORY_REAL_EXPRESSION_RESPONSES_NOT_BIOLOGICAL_VALIDATION'}
    (output/'original_source_response_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--audit',required=True,type=Path)
    p.add_argument('--panel',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    x=p.parse_args()
    report=execute(json.loads(x.audit.read_text()),json.loads(x.panel.read_text()),x.output)
    print(json.dumps({'stage':report['stage'],'per_tf':report['per_tf']},indent=2))


if __name__ == '__main__':
    main()
