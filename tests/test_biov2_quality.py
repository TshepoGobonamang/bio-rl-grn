from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import pytest
from scripts.biov2_quality import parse_v29_gtf, parse_quantification, build_qc, pairwise_correlations


def test_genecode_gtf_exact_gene_symbol(tmp_path):
    path=tmp_path/'v29.gtf.gz'
    with gzip.open(path,'wt') as out:
        out.write('# version 29\n')
        out.write('chr1\tGENCODE\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG00000000003.14"; gene_name "TFDP1"; gene_type "protein_coding";\n')
        out.write('chr1\tGENCODE\ttranscript\t100\t200\t.\t+\t.\tgene_id "ENSG00000000003.14"; gene_name "TFDP1";\n')
    assert parse_v29_gtf(path)=={'ENSG00000000003.14':'TFDP1'}


def test_identical_gene_replicate_synthetic_qc(tmp_path):
    names=['TFDP1','KLF2','NR2F2','RELA','SP1','NR2C2']
    ids=[f'ENSG{i:011d}.1' for i in range(1,1202)]
    symbol={id: name for id,name in zip(ids[:6],names)}
    accessions=['ENCSR095PIC','ENCSR552HTC','ENCSR620DOP','ENCSR688UTN',
                'ENCSR052ZPX','ENCSR650QAK','ENCSR808HQN']
    selection={};downloaded={}
    rng=np.random.default_rng(13)
    for k, acc in enumerate(accessions):
        s=[]
        for n in (1,2):
            fname=f'ENCFF{k:03d}{n:03d}'
            p=tmp_path/(fname+'.tsv')
            values=np.maximum(0, rng.normal(loc=20,scale=6,size=len(ids)))
            v=pd.DataFrame({'gene_id':ids,'expected_count':values*10,
                            'TPM':values,'effective_length':np.full(len(ids),800)})
            v.to_csv(p,sep='\t',index=False)
            downloaded[fname]=p
            s.append({'file_accession':fname,'replicate_number':n,
                      'biosample_id':f'B{k}_{n}'})
        selection[acc]=s
    report=build_qc(selection,downloaded,symbol,tmp_path/'output',annotation_sha256='x'*64)
    assert report['n_individual_biosamples']==14
    assert report['n_common_ensembl_gene_ids']==len(ids)
    assert len(pd.read_csv(tmp_path/'output'/'knockdown_mrna.csv'))==6
    assert (tmp_path/'output'/'SP1'/'rsem_expected_counts.csv.gz').exists()
    assert not report['gcev_predictions_evaluated']
    assert report['no_significance_computed_in_python_qc']


def test_quant_duplicate_gene_id_fails(tmp_path):
    path=tmp_path/'bad.tsv'
    pd.DataFrame({'gene_id':['ENSG00000000003.14']*2,'TPM':[1,2],
                  'expected_count':[5,6],'effective_length':[10,20]}).to_csv(path,sep='\t',index=False)
    with pytest.raises(ValueError,match='duplicate gene IDs'):
        parse_quantification(path)
