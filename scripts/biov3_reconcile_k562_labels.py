"""Verify publisher-derived K562 barcodes and extract a K562-only matrix.

This script will NEVER infer label 1/2 from markers, numeric factor codes, or
expression-based clustering. Trusted author cell annotations are essential.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _read_metadata(path: Path) -> dict[str, str]:
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not {'orig.ident', 'CellType'}.issubset(reader.fieldnames or []):
            raise ValueError('GEO metadata missing original orig.ident / CellType')
        obs = {}
        for r in reader:
            cell = (r['orig.ident'] or '').strip()
            code = (r['CellType'] or '').strip()
            if not cell or code not in {'1', '2'} or cell in obs:
                raise ValueError('Blank, duplicated or unexpected source cell/CellType code')
            obs[cell] = code
    if not obs:
        raise ValueError('Empty source metadata')
    return obs


def _read_author_labels(path: Path) -> dict[str, str]:
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if set(reader.fieldnames or []) != {'cell_id', 'cell_line'}:
            raise ValueError('Original author labels require exactly cell_id,cell_line columns')
        labels = {}
        for r in reader:
            cell = (r['cell_id'] or '').strip()
            label = (r['cell_line'] or '').strip()
            if not cell or label not in {'K562', 'BT142'} or cell in labels:
                raise ValueError('Author cell labels missing/duplicated or unsupported')
            labels[cell] = label
    if not labels:
        raise ValueError('Empty author label file')
    return labels


def audit_identity(meta: dict[str, str], labels: dict[str, str]) -> dict:
    if set(labels) != set(meta):
        raise ValueError(f'Author-to-GEO barcode mismatch: author_only={len(set(labels)-set(meta))}, GEO_only={len(set(meta)-set(labels))}')
    counts = Counter(labels.values())
    if min(counts.get('K562', 0), counts.get('BT142', 0)) < 2:
        raise ValueError('Both K562 and BT142 lines must be present in original author barcode reference')
    x = defaultdict(Counter)
    for cell, code in meta.items():
        x[code][labels[cell]] += 1
    mapping = {}
    for code, c in sorted(x.items()):
        if len(c) != 1:
            raise ValueError('Original GEO CellType code is not uniquely associated with author cell identities')
        mapping[code] = next(iter(c))
    if set(mapping.values()) != {'K562', 'BT142'}:
        raise ValueError('Cannot certify two distinct K562/BT142 numeric mappings')
    return {'author_counts': dict(counts),
            'code_to_line': mapping,
            'n_original_cells': len(meta),
            'n_author_matched_cells': len(labels),
            'original_labels_used_not_marker_predictions': True,
            'original_author_annotation_required': True}


def export_matrix(expression: Path, meta: dict[str, str], labels: dict[str, str], out: Path) -> dict:
    selected = {cell for cell, line in labels.items() if line == 'K562'}
    out.mkdir(parents=True, exist_ok=True)
    matrix = out / 'K562_gene_by_cell_original_counts.tsv.gz'
    genes = set()
    nnz = 0
    noninteger = 0
    with gzip.open(expression, 'rt', encoding='utf-8-sig', newline='') as source, \
         gzip.open(matrix, 'wt', encoding='utf-8', newline='') as dst:
        rows = csv.reader(source, delimiter='\t')
        header = next(rows)
        if len(header) != len(set(header)) or set(header) != set(meta):
            raise ValueError('Original expression cell header fails exact GEO metadata join')
        keep_idx = [i for i, cell in enumerate(header) if cell in selected]
        if len(keep_idx) != len(selected):
            raise ValueError('K562 cell IDs do not exactly match source header')
        writer = csv.writer(dst, delimiter='\t')
        writer.writerow(['gene'] + [header[i] for i in keep_idx])
        n_genes = 0
        for row in rows:
            if len(row) != len(header) + 1:
                raise ValueError(f'Original row length mismatch at gene row {n_genes + 1}')
            gene = row[0].strip()
            if not gene or gene in genes:
                raise ValueError('Original expression duplicate/empty gene symbol')
            genes.add(gene)
            values = [row[i+1] for i in keep_idx]
            for v in values:
                try:
                    z = float(v)
                except ValueError as exc:
                    raise ValueError('Expression matrix contains nonnumeric counts') from exc
                if not (z >= 0 and z < float('inf')):
                    raise ValueError('Expression counts must be finite and nonnegative')
                nnz += int(z > 0)
                noninteger += int(z != int(z))
            writer.writerow([gene] + values)
            n_genes += 1
    if n_genes == 0:
        raise ValueError('Original expression contains no gene rows')
    return {'k562_cell_count': len(selected),'n_genes_unmapped_symbols': n_genes,
            'nonzero_expression_entries': nnz,
            'noninteger_entries': noninteger,
            'matrix_file': matrix.name,
            'matrix_sha256': sha256(matrix),
            'matrix_orientation': 'gene_by_cell',
            'gene_id_provenance_status': 'ORIGINAL_GEO_SYMBOLS_NOT_YET_MAPPED_TO_GENCODE_V29',
            'eligible_for_gcev_model_fit': False}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--expression', type=Path, required=True)
    ap.add_argument('--metadata', type=Path, required=True)
    ap.add_argument('--author-labels', type=Path, required=True,
                    help='Exact original author RDS-derived cell_id,cell_line CSV (not inferred)')
    ap.add_argument('--author-rds-sha256', required=True,
                    help='SHA256 of actual original author Seurat RDS inspected in a separate source step')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--expected-expression-sha256', default='8bf3f175fe8635fec316263901bab18a4201d0199783482122142dd4387ffa75')
    ap.add_argument('--expected-metadata-sha256', default='ef87b7605684e66f52a9f53c8feb0863b9bf20a3aed37aa1081f6b18489558d3')
    args = ap.parse_args()
    if len(args.author_rds_sha256) != 64 or any(x not in '0123456789abcdef' for x in args.author_rds_sha256.lower()):
        raise ValueError('Original author RDS must have a recorded full SHA256')
    actual = {'expression': sha256(args.expression), 'metadata': sha256(args.metadata)}
    expected = {'expression':args.expected_expression_sha256, 'metadata': args.expected_metadata_sha256}
    if actual != expected:
        raise ValueError(f'GEO original source SHA256 differs from frozen audit: {actual}')
    meta = _read_metadata(args.metadata)
    author = _read_author_labels(args.author_labels)
    decision = audit_identity(meta,author)
    subset = export_matrix(args.expression,meta,author,args.output)
    result = {**decision, **subset,
              'original_source_sha256': actual,
              'author_rds_sha256': args.author_rds_sha256,
              'author_label_csv_sha256':sha256(args.author_labels),
              'encode_perturbation_responses_used':False,
              'regulatory_priors_mapped': False,
              'decision':'K562_MEMBERSHIP_VERIFIED_FROM_AUTHOR_BARCODE_LABELS_GENE_MAPPING_PENDING'}
    (args.output/'identity_and_export_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result, indent=2))

if __name__=='__main__': main()
