"""Inspect original ENCODE experiment JSON and fail closed on missing controls or replicates.

This audit NEVER reads expression matrices and NEVER constructs regulatory truth labels.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable

HOST = 'https://www.encodeproject.org'
ACCESSION = re.compile(r'ENCSR[A-Z0-9]{6}$')


def record_id(value: Any) -> str:
    """Return a stable ENCODE accession/UUID/path string, not a guessed label."""
    if isinstance(value, dict):
        return str(value.get('accession') or value.get('@id') or value.get('uuid') or '')
    return str(value or '')


def accession_of(value: Any) -> str:
    s = record_id(value)
    if ACCESSION.fullmatch(s):
        return s
    match = re.search(r'(ENCSR[A-Z0-9]{6})', s)
    return match.group(1) if match else ''


def fetch_json(url: str, *, timeout: int = 45) -> tuple[dict, str]:
    """Read one original public ENCODE JSON record, preserving its raw SHA256."""
    if not url.startswith(HOST + '/') or not (('/experiments/' in url) or ('/replicates/' in url) or ('/libraries/' in url)):
        raise ValueError('Unexpected original ENCODE record URL')
    req = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': 'GCEV-ENCODE-provenance-audit/0.2 academic metadata checker',
    })
    with urllib.request.urlopen(req, timeout=timeout) as response:
        data = response.read(8_000_001)
    if len(data) > 8_000_000:
        raise ValueError('ENCODE record exceeds safe metadata size')
    obj = json.loads(data)
    if not isinstance(obj, dict):
        raise ValueError('ENCODE response is not a JSON record')
    return obj, hashlib.sha256(data).hexdigest()


def resolve_replicates(obj: dict, fetch: Callable[[str], tuple[dict, str]]) -> dict:
    """Resolve original ENCODE replicate links; never infer biosample IDs from file counts.

    ENCODE experiment frame=embedded can still return string replicate paths;
    those MUST be fetched individually before auditing biological units.
    """
    resolved = []
    record_hashes = {}
    for item in (obj.get('replicates') or []):
        if isinstance(item, dict) and item.get('biological_replicate_number') is not None and _biosample(item):
            resolved.append(item)
            continue
        path = item.get('@id') if isinstance(item, dict) else str(item)
        if not isinstance(path, str) or not re.fullmatch(r'/replicates/[a-fA-F0-9-]{36}/', path):
            raise ValueError(f'Unverifiable ENCODE replicate reference: {path!r}')
        record, sha = fetch(HOST + path + '?format=json&frame=embedded')
        if record.get('@id') and record.get('@id') != path:
            raise ValueError(f'Replicate identity mismatch: {path}')
        # ENCODE's canonical model is replicate -> library -> biosample.
        # Do not invent an independent unit from a replicate number or file.
        if not _biosample(record) and isinstance(record.get('library'), str):
            library_path = record['library']
            if not re.fullmatch(r'/libraries/ENCLB[0-9A-Z]{6}/', library_path):
                raise ValueError(f'Unverifiable library reference: {library_path!r}')
            library, library_sha = fetch(HOST + library_path + '?format=json&frame=embedded')
            record['library'] = library
            record_hashes[library_path] = library_sha
        if not _biosample(record) or record.get('biological_replicate_number') is None:
            raise ValueError(f'Original replicate missing verified library biosample or number: {path}')
        if record.get('experiment'):
            exp = accession_of(record['experiment'])
            if exp and exp != accession_of(obj):
                raise ValueError(f'Replicate linked to a different experiment: {path}')
        record_hashes[path] = sha
        resolved.append(record)
    result = dict(obj)
    result['replicates'] = resolved
    result['_gcev_original_replicate_sha256'] = record_hashes
    return result


def fetch_original(accession: str, *, timeout: int = 45) -> tuple[dict, str]:
    url = f'{HOST}/experiments/{accession}/?format=json&frame=embedded'
    obj, source_sha = fetch_json(url, timeout=timeout)
    if accession_of(obj) != accession:
        raise ValueError(f'Unexpected ENCODE experiment identity for {accession}')
    return resolve_replicates(obj, lambda u: fetch_json(u, timeout=timeout)), source_sha


def jsonfile(directory: Path, accession: str) -> tuple[dict, str]:
    path = directory / (accession + '.json')
    content = path.read_bytes()
    obj = json.loads(content)
    if accession_of(obj) != accession:
        raise ValueError(f'Fixture identity mismatch: {accession}')
    return obj, hashlib.sha256(content).hexdigest()


def _biosample(r: dict) -> str:
    """Stable source biosample ID; support ENCODE replicate.library.biosample."""
    direct = r.get('biosample')
    if direct:
        return record_id(direct)
    library = r.get('library')
    if isinstance(library, dict) and library.get('biosample'):
        return record_id(library['biosample'])
    return ''


def _replicate(r: dict) -> int | None:
    v = r.get('biological_replicate_number')
    try:
        return int(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def inspect_experiment(obj: dict, *, expected_tf: str | None,
                       expected_control: str, raw_sha256: str) -> dict:
    accession = accession_of(obj)
    assays = ' '.join(str(obj.get(k) or '') for k in ('assay_title', 'assay_term_name')).lower()
    assay_ok = 'crispri' in assays and 'rna-seq' in assays
    context_ok = 'k562' in str(obj.get('biosample_summary') or obj.get('simple_biosample_summary') or '').lower()
    released = obj.get('status') == 'released'
    control_type = str(obj.get('control_type') or '')
    if expected_tf is None:
        source_control_ok = control_type.lower().strip() == 'non-specific target'
        target_ok = not obj.get('target')
        linked_control = None
    else:
        target_obj = obj.get('target')
        target_text = ((' '.join(str(target_obj.get(k) or '') for k in ('@id','name','label','accession')))
                       if isinstance(target_obj, dict) else record_id(target_obj))
        target_ok = expected_tf in re.split(r'[^A-Za-z0-9]+', target_text)
        control_paths = obj.get('possible_controls', []) or obj.get('controls', []) or []
        control_acc = {accession_of(x) for x in control_paths}
        source_control_ok = expected_control in control_acc
        linked_control = sorted(control_acc)

    reps = obj.get('replicates') or []
    ids, biosamples = [], []
    replicate_map = []
    unresolved = 0
    for item in reps:
        if not isinstance(item, dict):
            unresolved += 1  # do not count unresolved link-only replicate rows
            continue
        num, bio = _replicate(item), _biosample(item)
        if num is None or not bio:
            unresolved += 1
            continue
        ids.append(num)
        biosamples.append(bio)
        replicate_map.append({'biological_replicate_number': num, 'biosample_id': bio})
    biological = (
        str(obj.get('replication_type') or '').lower() == 'isogenic'
        and len(set(ids)) >= 2
        and len(set(biosamples)) >= 2
        and len(ids) == len(set(ids)) == len(set(biosamples))
    )
    files = obj.get('files') or []
    quant_files = []
    unresolved_files = 0
    for item in files:
        if not isinstance(item, dict):
            unresolved_files += 1
            continue
        if item.get('status') != 'released':
            continue
        if str(item.get('file_format') or '').lower() != 'tsv':
            continue
        if 'gene quantification' not in str(item.get('output_type') or '').lower():
            continue
        replicate_numbers = item.get('biological_replicates') or []
        if not isinstance(replicate_numbers, list):
            replicate_numbers = []
        try:
            replicate_numbers = [int(n) for n in replicate_numbers]
        except (TypeError, ValueError):
            replicate_numbers = []
        q = {
            'accession': record_id(item),
            'file_format': item.get('file_format'),
            'output_type': item.get('output_type'),
            'assembly': item.get('assembly'),
            'genome_annotation': item.get('genome_annotation'),
            'biological_replicates': replicate_numbers,
            'file_size': item.get('file_size'),
            'md5sum': item.get('md5sum'),
            'href': item.get('href'),
            'one_biological_replicate': len(replicate_numbers) == 1,
        }
        quant_files.append(q)
    qualified = [q for q in quant_files
                 if q['one_biological_replicate'] and q['assembly'] == 'GRCh38'
                 and q['biological_replicates'][0] in set(ids)]
    file_reps = {q['biological_replicates'][0] for q in qualified}
    annotations = {q['genome_annotation'] for q in qualified}
    quant_ok = len(file_reps) >= 2 and len(annotations) == 1 and None not in annotations

    audits = obj.get('audit') or obj.get('audits') or {}
    if not isinstance(audits, dict):
        audits = {'unparsed': str(audits)}
    audit_counts = {k: len(v) if isinstance(v, list) else 1 for k, v in audits.items()}
    findings = {
        'source': f'{HOST}/experiments/{accession}/',
        'source_json_sha256': raw_sha256,
        'accession': accession,
        'expected_tf': expected_tf,
        'assay_match': assay_ok,
        'cell_line_match': context_ok,
        'released': released,
        'target_match': target_ok,
        'linked_expected_control': source_control_ok,
        'linked_control_accessions': linked_control,
        'control_type': control_type,
        'replication_type': obj.get('replication_type'),
        'biological_replicate_numbers': sorted(set(ids)),
        'biosample_identifiers': sorted(set(biosamples)),
        'replicate_to_biosample': sorted(replicate_map, key=lambda x: x['biological_replicate_number']),
        'unresolved_replicate_references': unresolved,
        'original_replicate_record_sha256': obj.get('_gcev_original_replicate_sha256', {}),
        'biosample_replicates_verified_from_metadata': biological,
        'quantification_candidates': quant_files,
        'per_replicate_grch38_quantifications': qualified,
        'per_replicate_quantification_available': quant_ok,
        'unresolved_file_references': unresolved_files,
        'audit_counts': audit_counts,
        'manual_audit_review_required': bool(audit_counts) or not ('audit' in obj or 'audits' in obj),
        'original_qc_metadata_missing': not ('audit' in obj or 'audits' in obj),
    }
    findings['design_metadata_pass'] = bool(released and assay_ok and context_ok and target_ok
                                           and source_control_ok and biological)
    findings['source_quantification_pass'] = bool(findings['design_metadata_pass'] and quant_ok)
    return findings


def evaluate_panel(manifest: dict, fetch: Callable[[str], tuple[dict, str]]) -> dict:
    required = manifest['control']['accession']
    plan = [(None, required)] + [(p['tf'], p['accession']) for p in manifest['perturbations']]
    findings, errors = [], {}
    for tf, accession in plan:
        try:
            obj, checksum = fetch(accession)
            findings.append(inspect_experiment(obj, expected_tf=tf,
                                              expected_control=required,
                                              raw_sha256=checksum))
        except (ValueError, KeyError, TypeError, FileNotFoundError, urllib.error.URLError,
                TimeoutError, json.JSONDecodeError) as exc:
            errors[accession] = f'{type(exc).__name__}: {exc}'
    control = next((x for x in findings if x['accession'] == required), None)
    quant_annotations = {q['genome_annotation'] for x in findings
                         for q in x['per_replicate_grch38_quantifications']}
    panel_design = bool(len(findings) == len(plan) and not errors
                        and control and control['design_metadata_pass']
                        and all(x['design_metadata_pass'] for x in findings))
    panel_quant = bool(panel_design and all(x['source_quantification_pass'] for x in findings)
                       and len(quant_annotations) == 1 and None not in quant_annotations)
    return {
        'protocol_id': manifest['protocol_id'],
        'stage': 'original_source_metadata_only',
        'requested_experiments': len(plan),
        'retrieved_experiments': len(findings),
        'errors': errors,
        'experiments': findings,
        'shared_gene_annotation': next(iter(quant_annotations)) if len(quant_annotations) == 1 else None,
        'panel_design_metadata_pass': panel_design,
        'panel_per_replicate_quantification_pass': panel_quant,
        'quality_control_review_pending': any(x['manual_audit_review_required'] for x in findings),
        'perturbation_effects_computed': False,
        'source_disjointness_confirmed': False,
        'model_rankings_compared': False,
        'independent_biological_validation_completed': False,
        'next_step': 'inspect released per-replicate quantification, replicate quality, knockdown and source provenance',
        'decision': 'SOURCE_METADATA_QUALIFIED_FOR_NEXT_REVIEW' if panel_quant else 'NO_GO_FOR_BIOLOGICAL_SCORING',
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--offline-root', type=Path,
                        help='development fixture directory, never an actual-data source')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    if args.offline_root:
        fetch = lambda acc: jsonfile(args.offline_root, acc)
    else:
        fetch = fetch_original
    report = evaluate_panel(manifest, fetch)
    report['data_origin'] = 'TEST_FIXTURES_NOT_REAL' if args.offline_root else 'ORIGINAL_ENCODE_API'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('data_origin','requested_experiments',
          'retrieved_experiments','decision','errors')}, indent=2))
    # Nonzero stops a workflow if metadata missing/unqualified; JSON is retained for diagnosis.
    return 0 if report['panel_per_replicate_quantification_pass'] else 2


if __name__ == '__main__':
    sys.exit(main())
