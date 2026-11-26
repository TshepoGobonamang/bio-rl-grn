"""BIO-V2: original ENCODE K562 CRISPRi replicate QC and DESeq2 inputs.

Only expression metadata, direct mRNA target knockdown checks, sample QC and
raw estimated-count/TPM/length matrices. No causal edges, DE significance,
or GCEV model performance are inferred by this module.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


GENE = re.compile(r'ENSG[0-9]{11}(?:\.[0-9]+)?(?:_PAR_Y)?$')
ATTR = re.compile(r'([A-Za-z_]+) "([^"]*)";')
TARGET_TFS = ('TFDP1', 'KLF2', 'NR2F2', 'RELA', 'SP1', 'NR2C2')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as inp:
        for part in iter(lambda: inp.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def parse_v29_gtf(path: Path) -> dict[str, str]:
    """Use source release gene_id and gene_name exactly; no alias guessing."""
    result = {}
    with gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else path.open('rt', encoding='utf-8') as inp:
        for line in inp:
            if not line or line[0] == '#':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != 'gene':
                continue
            attrs = dict(ATTR.findall(cols[8]))
            gene_id, name = attrs.get('gene_id'), attrs.get('gene_name')
            if not gene_id or not name or not GENE.fullmatch(gene_id):
                continue
            if gene_id in result and result[gene_id] != name:
                raise ValueError('GENCODE V29 conflicting gene symbol for gene id: ' + gene_id)
            result[gene_id] = name
    if not result:
        raise ValueError('No GENCODE gene records parsed')
    return result


def parse_quantification(path: Path) -> pd.DataFrame:
    columns = ('gene_id', 'expected_count', 'TPM', 'effective_length')
    table = pd.read_csv(path, sep='\t', dtype={'gene_id': 'string'}, usecols=list(columns))
    if any(x not in table for x in columns) or table['gene_id'].isna().any():
        raise ValueError(f'{path.name}: missing gene identifiers or expected fields')
    if table['gene_id'].duplicated().any():
        raise ValueError(f'{path.name}: duplicate gene IDs')
    for col in columns[1:]:
        table[col] = pd.to_numeric(table[col], errors='raise')
        if not np.isfinite(table[col]).all() or (table[col] < 0).any():
            raise ValueError(f'{path.name}: invalid {col} values')
    return table.set_index('gene_id')


def pairwise_correlations(tpm: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    # Rows = genes; columns = independent biosamples. Spearman across all ENSG only.
    corr = np.log2(1 + tpm).corr(method='spearman')
    rows = []
    for i in range(len(meta)):
        for j in range(i + 1, len(meta)):
            a, b = meta.iloc[i], meta.iloc[j]
            rows.append({'sample_a': a['sample_id'], 'sample_b': b['sample_id'],
                         'experiment_a': a['experiment'], 'experiment_b': b['experiment'],
                         'same_experiment': bool(a['experiment'] == b['experiment']),
                         'r_spearman_log2_tpm': float(corr.loc[a['sample_id'],b['sample_id']])})
    return pd.DataFrame(rows)


def build_qc(selected: dict, downloaded: dict, gene_symbol: dict[str,str], output: Path,
             *, annotation_sha256: str) -> dict:
    """Compute QC without selecting/excluding individual biosamples post hoc."""
    output.mkdir(parents=True, exist_ok=True)
    rows, frames = [], {}
    for accession, items in selected.items():
        for source in items:
            file_accession = source['file_accession']
            path = Path(downloaded[file_accession])
            frame = parse_quantification(path)
            sample_id = f'{accession}_R{source["replicate_number"]}'
            if sample_id in frames:
                raise ValueError('Duplicate sample IDs')
            frames[sample_id] = frame
            gene_idx = frame.index.astype(str)
            is_ensg = np.asarray(gene_idx.str.match(GENE.pattern), dtype=bool)
            counts = frame['expected_count'].to_numpy()
            tpm = frame['TPM'].to_numpy()
            rows.append({'sample_id':sample_id, 'experiment':accession,
                         'biological_replicate_number': int(source['replicate_number']),
                         'biosample_id':source['biosample_id'],
                         'source_file_accession':file_accession,
                         'source_sha256':sha256(path),
                         'total_rsem_expected_count':float(counts.sum()),
                         'n_biological_ensg_ids':int(is_ensg.sum()),
                         'nonzero_ensg_tpm':int(np.sum((tpm>0)&is_ensg)),
                         'nonzero_ensg_count':int(np.sum((counts>0)&is_ensg)),
                         'n_total_quantified_rows':int(len(frame))})
    meta = pd.DataFrame(rows)
    if meta['biosample_id'].duplicated().any():
        raise ValueError('Biosample repeated: cannot treat technical repeats as independent')
    gene_sets = [set(frame.index.astype(str)) for frame in frames.values()]
    common = sorted(set.intersection(*gene_sets))
    ens = [x for x in common if GENE.fullmatch(x)]
    if len(ens) < 1000:
        raise ValueError(f'Only {len(ens)} common biological Ensembl genes')
    tpm = pd.DataFrame({k:v.loc[ens,'TPM'].to_numpy(dtype=float) for k,v in frames.items()}, index=ens)
    pairs = pairwise_correlations(tpm, meta)
    within = pairs[pairs['same_experiment']]
    if len(within) < 7:
        raise ValueError('Missing within-experiment biological comparisons')
    mean_expr = np.log2(tpm+1)
    var = mean_expr.var(axis=1)
    selected_features = var.nlargest(min(2000, len(var))).index
    x = mean_expr.loc[selected_features].to_numpy().T
    x -= x.mean(axis=0, keepdims=True)
    # PCA via SVD, samples are biosamples; eigenvectors only exploratory
    u, singular, vt = np.linalg.svd(x, full_matrices=False)
    coords = u[:, :2] * singular[:2]
    if coords.shape[1] == 1:
        coords = np.column_stack((coords, np.zeros(len(coords))))
    pca = meta[['sample_id','experiment','biosample_id']].copy()
    pca['PC1'] = coords[:,0]; pca['PC2']=coords[:,1]
    variance_fraction = (singular[:2]**2 / max(float((singular**2).sum()), 1e-12)).tolist()
    # Direct mRNA target reduction is assay QC, not a guarantee of protein knockdown
    id_by_symbol = {}
    for gid,symbol in gene_symbol.items():
        if symbol in TARGET_TFS and gid in tpm.index:
            id_by_symbol.setdefault(symbol,[]).append(gid)
    control = 'ENCSR095PIC'
    ctrl_cols = meta.loc[meta['experiment']==control,'sample_id'].tolist()
    if len(ctrl_cols)<2:
        raise ValueError('Fewer than two distinct controls')
    knockdown=[]
    per_tf=[]
    for symbol,accession in [('TFDP1','ENCSR552HTC'),('KLF2','ENCSR620DOP'),
                             ('NR2F2','ENCSR688UTN'),('RELA','ENCSR052ZPX'),
                             ('SP1','ENCSR650QAK'),('NR2C2','ENCSR808HQN')]:
        treatment_cols = meta.loc[meta['experiment']==accession,'sample_id'].tolist()
        if len(treatment_cols) < 2:
            raise ValueError(f'{symbol}: fewer than two replicate biosamples')
        gid_list = id_by_symbol.get(symbol,[])
        if len(gid_list)!=1:
            knockdown.append({'tf':symbol,'experiment':accession,'gene_id':None,
                              'status':'UNRESOLVED_GENE_SYMBOL_MAPPING',
                              'candidate_matches':len(gid_list)})
        else:
            gid = gid_list[0]
            treated = tpm.loc[gid,treatment_cols]
            control_meas = tpm.loc[gid,ctrl_cols]
            ctl=float(control_meas.mean()); trt=float(treated.mean())
            kd_ratio=(trt+1)/(ctl+1)
            knockdown.append({'tf':symbol,'experiment':accession,'gene_id':gid,
                              'status':'DESCRIPTIVE_MRNA_CHECK_REQUIRES_QC',
                              'control_tpm_mean':ctl,'perturbed_tpm_mean':trt,
                              'control_TPM_by_biosample_json':json.dumps({str(meta.set_index('sample_id').loc[col,'biosample_id']):float(tpm.loc[gid,col]) for col in ctrl_cols},sort_keys=True),
                              'perturbed_TPM_by_biosample_json':json.dumps({str(meta.set_index('sample_id').loc[col,'biosample_id']):float(tpm.loc[gid,col]) for col in treatment_cols},sort_keys=True),
                              'log2_tpm_ratio_pseudo1':float(np.log2(kd_ratio)),
                              'percent_mrna_reduction_pseudocount':float(100*(1-kd_ratio)),
                              'treatment_replicates_below_control_mean':int((treated<ctl).sum()),
                              'treatment_biosamples':len(treatment_cols),
                              'control_biosamples':len(ctrl_cols)})
        sample_cols=ctrl_cols+treatment_cols
        sid_to_label={col:'C_R'+str(i+1) for i,col in enumerate(ctrl_cols)} | {col:'T_R'+str(i+1) for i,col in enumerate(treatment_cols)}
        fnames = {}
        for measure,filename in [('expected_count','rsem_expected_counts.csv.gz'),
                                 ('TPM','rsem_TPM.csv.gz'),
                                 ('effective_length','rsem_effective_length.csv.gz')]:
            # Export full common ENSG set, then R DESeq2 explicitly filters
            mat = pd.DataFrame({sid_to_label[col]:frames[col].loc[ens,measure].to_numpy(dtype=float) for col in sample_cols}, index=ens)
            mat.index.name='gene_id'
            dest=output/symbol/filename
            dest.parent.mkdir(parents=True, exist_ok=True)
            mat.to_csv(dest, compression='gzip')
            fnames[measure]={'path':str(dest.relative_to(output)), 'sha256':sha256(dest)}
        sample_map = pd.DataFrame([{'sample_id':sid_to_label[col],
                 'original_biosample_id':meta.set_index('sample_id').loc[col,'biosample_id'],
                 'original_replicate_number':int(meta.set_index('sample_id').loc[col,'biological_replicate_number']),
                 'condition': 'control' if col in ctrl_cols else 'perturbed'} for col in sample_cols])
        csvpath=output/symbol/'coldata.csv'
        sample_map.to_csv(csvpath,index=False)
        per_tf.append({'tf':symbol,'experiment':accession,'control':control,
                       'n_controls':len(ctrl_cols),'n_treated':len(treatment_cols),
                       'n_common_ensg':len(ens), 'estimated_counts_are_rsem':True,
                       'count_length_bias_correction_required':True,
                       'files':fnames, 'coldata_sha256':sha256(csvpath)})
    for name,data in [('sample_quality.csv',meta),('pairwise_correlations.csv',pairs),
                      ('sample_pca.csv',pca),('knockdown_mrna.csv',pd.DataFrame(knockdown))]:
        data.to_csv(output/name,index=False)
    summary={'stage':'BIOV2_EXPLORATORY_REPLICATE_QC_ONLY',
             'annotation':'GENCODE V29 GRCh38','annotation_gtf_sha256':annotation_sha256,
             'n_individual_biosamples':len(meta), 'n_control_biosamples':len(ctrl_cols),
             'n_common_genes_all_types':len(common),'n_common_ensembl_gene_ids':len(ens),
             'n_gencode_genes_parsed':len(gene_symbol),'n_ensg_mapped_to_gencode':sum(x in gene_symbol for x in ens),
             'n_within_experiment_pairs':len(within),
             'within_experiment_spearman_min':float(within['r_spearman_log2_tpm'].min()),
             'within_experiment_spearman_median':float(within['r_spearman_log2_tpm'].median()),
             'pca_fraction_variance_PC1_PC2':variance_fraction,
             'knockdown_efficacy_qc_pending':True,'encode_original_qc_review_pending':True,
             'source_disjointness_verified':False,'gcev_predictions_evaluated':False,
             'regulatory_direct_binding_verified':False,
             'no_significance_computed_in_python_qc':True,'tf_contrasts':per_tf}
    (output/'qc_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def run(panel_file:Path, output:Path, gtf_file:Path) -> dict:
    from scripts.audit_encode_crispri_panel import evaluate_panel, fetch_original
    from scripts.execute_encode_real_responses import selected_files, download_checked
    panel=json.loads(panel_file.read_text())
    if panel.get('genome_annotation')!='V29' or panel.get('genome_assembly')!='GRCh38':
        raise ValueError('Frozen annotation mismatch')
    report=evaluate_panel(panel,lambda a:fetch_original(a))
    # evaluate_panel is a pure-function validator; the existing ENCODE intake
    # requires the authenticated input-origin tag normally set by its CLI.
    report['data_origin']='ORIGINAL_ENCODE_API'
    output.mkdir(parents=True,exist_ok=True)
    (output/'original_source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    selected=selected_files(report,panel)
    mapping=parse_v29_gtf(gtf_file)
    with tempfile.TemporaryDirectory(prefix='encode_biov2_') as temp:
        files={}
        for acc,entries in selected.items():
            for q in entries:
                destination=Path(temp)/(q['file_accession']+'.tsv')
                provenance=download_checked(q,destination)
                files[q['file_accession']]=destination
                print('VERIFIED',provenance['accession'],provenance['biosample_id'],flush=True)
        return build_qc(selected,files,mapping,output,annotation_sha256=sha256(gtf_file))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--panel',type=Path,required=True)
    parser.add_argument('--gtf',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=run(args.panel,args.output,args.gtf)
    print(json.dumps({k:report[k] for k in ('stage','n_individual_biosamples','n_common_ensembl_gene_ids','within_experiment_spearman_min')},indent=2))

if __name__=='__main__':main()
