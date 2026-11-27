"""BIO-V8 ENCODE metadata discovery; never fetch an assay data payload.

This module qualifies metadata candidates, not regulatory truth.  A successful
screen does not authorize outcome retrieval, certify knockdown efficacy, or
establish independence from motif-derivation experiments.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HOST = "https://www.encodeproject.org"
ASSAYS = ("CRISPRi RNA-seq", "shRNA RNA-seq", "knockdown RNA-seq", "TF ChIP-seq")
COLLECTIONS = {
    "experiments", "replicates", "libraries", "biosamples", "targets", "genes",
    "files", "treatments", "genetic-modifications",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def validate_metadata_url(url: str) -> str:
    """A transport boundary: even redirects must stay on JSON metadata URLs."""
    p = urllib.parse.urlsplit(url)
    if p.scheme != "https" or p.hostname not in {"www.encodeproject.org", "encodeproject.org"}:
        raise ValueError("Only original ENCODE HTTPS metadata is permitted")
    if p.username or p.password or p.port not in (None, 443) or p.fragment:
        raise ValueError("Credentials, alternate ports, and fragments are forbidden")
    decoded = p.path
    for _ in range(3):
        decoded = urllib.parse.unquote(decoded)
    if "@@" in decoded or ".." in decoded or "\\" in decoded or "." in decoded:
        raise ValueError("Assay/download payloads and non-record paths are sealed")
    pieces = decoded.strip("/").split("/")
    is_search = pieces == ["search"]
    if not is_search and not (
        len(pieces) == 2 and pieces[0] in COLLECTIONS
        and re.fullmatch(r"[A-Za-z0-9_:+-]+", pieces[1])
    ):
        raise ValueError("URL is not an allowed metadata collection record")
    query = urllib.parse.parse_qs(p.query, keep_blank_values=True)
    if query.get("format") != ["json"]:
        raise ValueError("Explicit format=json is required")
    if not is_search and set(query) - {"format", "frame"}:
        raise ValueError("Unexpected metadata record query parameter")
    if any("@@" in value or "download" in value.lower()
           for values in query.values() for value in values):
        raise ValueError("Download query is forbidden")
    return url


class MetadataRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        validate_metadata_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def record_url(value: Any, collection: str | None = None) -> str:
    if isinstance(value, dict):
        value = value.get("@id") or value.get("accession")
    if not isinstance(value, str):
        raise ValueError("Metadata reference lacks an original identifier")
    if value.startswith("https://"):
        raw = value.split("?")[0]
    elif value.startswith("/"):
        raw = HOST + value
    elif collection:
        raw = HOST + "/" + collection + "/" + value + "/"
    else:
        raise ValueError("Cannot infer metadata collection")
    return validate_metadata_url(raw.rstrip("/") + "/?format=json&frame=embedded")


def identity(value: Any) -> str:
    if isinstance(value, dict):
        value = value.get("accession") or value.get("@id") or value.get("uuid") or ""
    text = str(value or "")
    return text.strip("/").split("/")[-1]


class MetadataClient:
    def __init__(self, output: Path, cache: Path | None = None, timeout: int = 30):
        self.raw_dir = output / "raw_json"
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.cache = cache
        self.timeout = timeout
        self.ledger: list[dict] = []
        self.records: dict[str, dict] = {}
        self.opener = urllib.request.build_opener(MetadataRedirect())

    def get(self, url: str) -> dict:
        validate_metadata_url(url)
        if url in self.records:
            return self.records[url]
        key = hashlib.sha256(url.encode()).hexdigest()
        entry = {"url": url, "requested_at_utc": utc_now(), "request_kind": "metadata_json",
                 "payload_opened": False, "cache_key": key}
        try:
            source = self.cache / (key + ".json") if self.cache else None
            sidecar = self.cache / (key + ".request.json") if self.cache else None
            if source and source.exists() and sidecar and sidecar.exists():
                provenance = json.loads(sidecar.read_text())
                if provenance.get("url") != url:
                    raise ValueError("Cache URL does not match requested original source")
                raw = source.read_bytes()
                if hashlib.sha256(raw).hexdigest() != provenance.get("sha256"):
                    raise ValueError("Cached original metadata checksum mismatch")
                entry["retrieved_at_utc"] = provenance.get("retrieved_at_utc")
                if not entry["retrieved_at_utc"]:
                    raise ValueError("Cached metadata lacks original retrieval time")
                entry["transport"] = "verified_metadata_cache"
            else:
                req = urllib.request.Request(url, headers={
                    "Accept": "application/json", "User-Agent": "GCEV-BIOV8-metadata-only/1.0",
                })
                try:
                    response = self.opener.open(req, timeout=self.timeout)
                except urllib.error.HTTPError as exc:
                    # ENCODE may return HTTP 404 for a valid search with total=0.
                    # A failed query without that original JSON evidence is NOT empty.
                    if exc.code != 404 or urllib.parse.urlsplit(url).path != "/search/":
                        raise
                    if "json" not in exc.headers.get("Content-Type", "").lower():
                        raise
                    raw = exc.read(32_000_001)
                    empty = json.loads(raw)
                    if empty.get("total") != 0 or empty.get("@graph") != []:
                        raise ValueError("HTTP 404 search did not establish a valid zero-result query")
                    entry["http_status"] = 404
                else:
                    with response:
                        validate_metadata_url(response.geturl())
                        content_type = response.headers.get("Content-Type", "").lower()
                        if "json" not in content_type:
                            raise ValueError("Response is not JSON metadata; refusing payload")
                        raw = response.read(32_000_001)
                        entry["http_status"] = response.status
                entry["retrieved_at_utc"] = utc_now()
                entry["transport"] = "original_encode_https"
            if len(raw) > 32_000_000:
                raise ValueError("Metadata record exceeded 32 MB limit")
            obj = json.loads(raw)
            if not isinstance(obj, dict):
                raise ValueError("Metadata response must be a JSON object")
            if obj.get("status") == "error" or obj.get("@type") == ["Error"]:
                raise ValueError("ENCODE returned an error record")
            entry.update({"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw), "status": "ok"})
            (self.raw_dir / (key + ".json")).write_bytes(raw)
            (self.raw_dir / (key + ".request.json")).write_text(json.dumps(entry, indent=2) + "\n")
            self.records[url] = obj
            self.ledger.append(entry)
            return obj
        except Exception as exc:
            entry.update({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
            self.ledger.append(entry)
            raise


def list_value(value: Any) -> list:
    return value if isinstance(value, list) else ([] if value is None else [value])


def hydrate_experiment(obj: dict, client: MetadataClient, *, include_files: bool = True,
                       identity_only: bool = False) -> dict:
    """Resolve source metadata links, with bounded recursion and no data URLs."""
    reference_keys = {
        "target": "targets", "targets": "targets", "genes": "genes", "gene": "genes",
        "modified_site_by_target_id": "targets", "modified_gene": "genes",
        "modified_genes": "genes", "replicates": "replicates", "library": "libraries",
        "biosample": "biosamples", "treatments": "treatments",
        "genetic_modifications": "genetic-modifications",
        "applied_modifications": "genetic-modifications",
    }
    if identity_only:
        reference_keys = {k: v for k, v in reference_keys.items()
                          if k not in {"replicates", "library", "biosample"}}
    errors: list[str] = []

    def expand(value: Any, collection: str, depth: int) -> Any:
        if depth > 7:
            errors.append("Metadata reference recursion limit reached")
            return value
        if isinstance(value, list):
            return [expand(v, collection, depth) for v in value]
        if isinstance(value, str):
            try:
                value = client.get(record_url(value, collection))
            except Exception as exc:
                errors.append(str(exc))
                return {"@id": value, "_unresolved": True}
        if not isinstance(value, dict):
            return value
        required = {"genes": "symbol", "targets": "genes", "replicates": "biological_replicate_number",
                    "libraries": "biosample"}.get(collection)
        partial = (required is not None and required not in value) or ("@id" in value and len(value) <= 2)
        if partial and (value.get("@id") or value.get("accession")):
            try:
                value = client.get(record_url(value, collection))
            except Exception as exc:
                errors.append(str(exc))
                return {**value, "_unresolved": True}
        result = dict(value)
        for key, child_collection in reference_keys.items():
            if key in value:
                result[key] = expand(value[key], child_collection, depth + 1)
        return result

    result = expand(obj, "experiments", 0)
    if not include_files:
        result["_metadata_resolution_errors"] = errors
        return result
    resolved_files = []
    for value in obj.get("files", []):
        try:
            if isinstance(value, str):
                value = client.get(record_url(value, "files"))
            elif (isinstance(value, dict) and value.get("status") == "released"
                  and value.get("file_format") in {"tsv", "bed"}):
                # Full file metadata carries audit/replicate provenance absent in search frames.
                value = client.get(record_url(value, "files"))
            resolved_files.append(value)
        except Exception as exc:
            errors.append(str(exc))
            resolved_files.append({"@id": identity(value), "_unresolved": True})
    result["files"] = resolved_files
    result["_metadata_resolution_errors"] = errors
    return result


def target_identity(obj: dict) -> dict:
    canonical: set[str] = set()
    gene_ids: dict[str, set[str]] = defaultdict(set)
    labels: set[str] = set()
    paths: list[str] = []
    target_keys = {"target", "targets", "modified_site_by_target_id", "modified_gene", "modified_genes"}
    container_keys = {"replicates", "library", "biosample", "treatments", "genetic_modifications", "applied_modifications"}

    def target(value: Any, path: str) -> None:
        for item in list_value(value):
            if not isinstance(item, dict):
                if isinstance(item, str):
                    labels.add(identity(item).removesuffix("-human"))
                continue
            paths.append(path)
            for key in ("label", "name"):
                if item.get(key):
                    labels.add(str(item[key]).removesuffix("-human"))
            for gene in ([item] if item.get("symbol") else []) + list_value(item.get("genes")) + list_value(item.get("gene")):
                if isinstance(gene, dict) and gene.get("symbol"):
                    symbol = str(gene["symbol"])
                    canonical.add(symbol)
                    for key in ("geneid", "ensembl_id", "@id", "name"):
                        gene_ids[symbol].update(re.findall(r"ENSG\d{11}", str(gene.get(key, ""))))

    def walk(value: Any, path: str) -> None:
        if isinstance(value, list):
            for index, item in enumerate(value):
                walk(item, path + f"[{index}]")
        elif isinstance(value, dict):
            for key in target_keys:
                if key in value:
                    target(value[key], path + "." + key)
            for key in container_keys:
                if key in value:
                    walk(value[key], path + "." + key)
    walk(obj, "experiment")
    alias = "TCF4" in labels and ("TCF7L2" in canonical or not canonical)
    return {"canonical_symbols": sorted(canonical), "target_labels": sorted(labels),
            "canonical_gene_ids": {k: sorted(v) for k, v in sorted(gene_ids.items())},
            "target_metadata_paths": sorted(set(paths)), "tcf4_alias_ambiguity": alias,
            "canonical_identity_verified": bool(canonical) and len(canonical) == 1}


def biosample_id(replicate: dict) -> str:
    bio = replicate.get("biosample")
    if not bio and isinstance(replicate.get("library"), dict):
        bio = replicate["library"].get("biosample")
    if isinstance(bio, dict) and bio.get("_unresolved"):
        return ""
    return identity(bio)


def replicate_metadata(obj: dict) -> dict:
    mapping: dict[int, set[str]] = defaultdict(set)
    unresolved = 0
    wrong_experiment = 0
    for r in obj.get("replicates", []):
        if not isinstance(r, dict) or r.get("_unresolved"):
            unresolved += 1
            continue
        if r.get("experiment") and identity(r["experiment"]) != identity(obj):
            wrong_experiment += 1
            continue
        try:
            number = int(r["biological_replicate_number"])
        except (KeyError, ValueError, TypeError):
            unresolved += 1
            continue
        bio = biosample_id(r)
        if number < 1 or not bio:
            unresolved += 1
            continue
        mapping[number].add(bio)
    unique_bios = set().union(*mapping.values()) if mapping else set()
    one_to_one = all(len(v) == 1 for v in mapping.values()) and len(unique_bios) == len(mapping)
    return {"biological_replicate_numbers": sorted(mapping), "biosample_identifiers": sorted(unique_bios),
            "replicate_biosample_map": {str(k): sorted(v) for k, v in sorted(mapping.items())},
            "unresolved_replicates": unresolved,
            "wrong_experiment_replicates": wrong_experiment,
            "independent_biosample_metadata_pass": len(mapping) >= 2 and len(unique_bios) >= 2 and one_to_one and not unresolved and not wrong_experiment}


def assay_family(obj: dict) -> str:
    title = str(obj.get("assay_title", ""))
    term = str(obj.get("assay_term_name", ""))
    text = (title + " " + term).lower()
    if "chip-seq" in text and "histone" not in title.lower():
        return "TF ChIP-seq"
    if "rna-seq" in text:
        if "crispri" in text:
            return "CRISPRi RNA-seq"
        if "shrna" in text:
            return "shRNA RNA-seq"
        if "knockdown" in text:
            return "knockdown RNA-seq"
    return "unresolved"


def is_k562(obj: dict) -> bool:
    ontology = obj.get("biosample_ontology") or {}
    if isinstance(ontology, dict) and ontology.get("term_name"):
        return str(ontology["term_name"]).upper() == "K562"
    return bool(re.search(r"(?<![A-Za-z0-9])K562(?![A-Za-z0-9])", str(obj.get("biosample_summary", "")), re.I))


def qc_metadata(obj: dict) -> dict:
    severity: dict[str, list] = defaultdict(list)
    missing = "audit" not in obj and "audits" not in obj
    for source, value in [("experiment", obj)] + [("file:" + identity(f), f) for f in obj.get("files", []) if isinstance(f, dict)]:
        audit = value.get("audit", value.get("audits", {}))
        if not isinstance(audit, dict):
            severity["UNPARSED"].append({"source": source, "detail": audit})
            continue
        for level, entries in audit.items():
            for entry in list_value(entries):
                severity[str(level).upper().replace(" ", "_")].append({"source": source, "detail": entry})
    blocked = missing or any(severity.get(k) for k in ("ERROR", "NOT_COMPLIANT", "UNPARSED"))
    return {"qc_audits": dict(severity), "qc_audit_counts": {k: len(v) for k, v in severity.items()},
            "qc_metadata_missing": missing, "blocking_qc_flag": blocked,
            "manual_qc_review_required": missing or any(severity.values())}


def file_metadata(obj: dict) -> list[dict]:
    fields = ("accession", "@id", "status", "file_format", "file_format_type", "output_type", "assembly",
              "genome_annotation", "biological_replicates", "technical_replicates", "md5sum", "file_size",
              "dataset", "derived_from", "quality_metrics", "controlled_by", "replicate", "audit")
    rows = []
    for value in obj.get("files", []):
        if isinstance(value, dict):
            rows.append({key: value.get(key) for key in fields})
    return rows


def quantification_annotations(obj: dict, assembly: str, annotation: str | None) -> dict[str, set[int]]:
    valid_reps = set(replicate_metadata(obj)["biological_replicate_numbers"])
    by_annotation: dict[str, set[int]] = defaultdict(set)
    for f in file_metadata(obj):
        reps = f.get("biological_replicates") or []
        if not isinstance(reps, list) or len(reps) != 1:
            continue
        try:
            rep = int(reps[0])
        except (TypeError, ValueError):
            continue
        ann = f.get("genome_annotation")
        if (f.get("status") == "released" and f.get("assembly") == assembly
                and f.get("file_format") == "tsv" and str(f.get("output_type", "")).lower() == "gene quantifications"
                and ann and rep in valid_reps and (annotation is None or ann == annotation)
                and (not f.get("dataset") or identity(f["dataset"]) == identity(obj))):
            by_annotation[str(ann)].add(rep)
    return dict(by_annotation)


def inspect_experiment(obj: dict, controls: list[dict], regulators: list[str], protocol: dict) -> dict:
    target = target_identity(obj)
    exact = sorted(set(target["canonical_symbols"]) & set(regulators))
    label_only = sorted(set(target["target_labels"]) & set(regulators))
    family = assay_family(obj)
    reps = replicate_metadata(obj)
    qc = qc_metadata(obj)
    excluded = set(protocol.get("excluded_outcome_experiments", []))
    linked = {identity(c) for c in list_value(obj.get("possible_controls")) + list_value(obj.get("controls"))}
    assembly = protocol.get("assembly", "GRCh38")
    annotation = protocol.get("expected_genome_annotation")
    treatment_ann = quantification_annotations(obj, assembly, annotation)
    control_records = []
    control_passes = []
    for control in controls:
        cid = identity(control)
        cr = replicate_metadata(control)
        cq = qc_metadata(control)
        control_type = str(control.get("control_type", "")).lower().strip()
        ct = target_identity(control)
        linked_ok = cid in linked
        same_context = is_k562(control) and control.get("status") == "released"
        control_family = assay_family(control)
        fresh = cid not in excluded
        distinct = not (set(cr["biosample_identifiers"]) & set(reps["biosample_identifiers"]))
        no_target = not ct["canonical_symbols"] and not ct["target_labels"]
        if family == "TF ChIP-seq":
            appropriate = ("input" in control_type or "control" in str(control.get("assay_title", "")).lower()) and no_target
            matches = []
        else:
            appropriate = control_type in {"non-specific target", "non-targeting control", "non-targeting", "non-specific-target"} and no_target and control_family == family
            control_ann = quantification_annotations(control, assembly, annotation)
            matches = sorted(a for a, rr in treatment_ann.items() if len(rr) >= 2 and len(control_ann.get(a, set())) >= 2)
        ok = (linked_ok and same_context and fresh and distinct and appropriate
              and cr["independent_biosample_metadata_pass"] and not cq["blocking_qc_flag"]
              and not control.get("_metadata_resolution_errors"))
        if family != "TF ChIP-seq":
            ok = ok and bool(matches)
        control_passes.append(bool(ok))
        control_records.append({"accession": cid, "linked": linked_ok, "control_type": control_type,
                                "assay_family": control_family, "appropriate_control": appropriate,
                                "same_released_k562_context": same_context, "fresh_control_metadata": fresh,
                                "distinct_treatment_control_biosamples": distinct, "matching_annotations": matches,
                                "metadata_pass": bool(ok), "replicates": cr, "qc": cq,
                                "files": file_metadata(control), "resolution_errors": control.get("_metadata_resolution_errors", [])})
    expected_ids = protocol.get("regulator_ensembl_ids", {})
    gene_id_ok = (len(exact) == 1 and (not expected_ids or
                  set(target["canonical_gene_ids"].get(exact[0], [])) == {expected_ids.get(exact[0])}))
    identity_ok = target["canonical_identity_verified"] and len(exact) == 1 and not target["tcf4_alias_ambiguity"] and gene_id_ok
    allowed_assay = family in protocol.get("allowed_assay_families", ASSAYS)
    base = bool(obj.get("status") == "released" and is_k562(obj) and identity_ok and allowed_assay
                and reps["independent_biosample_metadata_pass"] and not qc["blocking_qc_flag"]
                and identity(obj) not in excluded and not obj.get("_metadata_resolution_errors"))
    matched_control = any(control_passes)
    pooled_idr = []
    if family == "TF ChIP-seq":
        for f in file_metadata(obj):
            fr = f.get("biological_replicates") or []
            if (f.get("status") == "released" and f.get("assembly") == assembly
                    and f.get("file_format") == "bed" and "idr" in str(f.get("output_type", "")).lower()
                    and len(set(fr)) >= 2):
                pooled_idr.append(f.get("accession"))
    response_pass = base and matched_control and family != "TF ChIP-seq"
    binding_candidate = base and matched_control and family == "TF ChIP-seq" and bool(pooled_idr)
    status = ("RESPONSE_METADATA_SCREEN_PASS_OUTCOMES_SEALED" if response_pass else
              "BINDING_METADATA_CANDIDATE_ONLY" if binding_candidate else "NO_GO_PENDING_METADATA_OR_QC")
    return {"accession": identity(obj), "assay_family": family, "assay_title": obj.get("assay_title"),
            "status": obj.get("status"), "cell_line_match": is_k562(obj), **target,
            "matched_regulators": exact, "label_only_candidates": label_only,
            "target_identity_pass": identity_ok, "canonical_gene_id_match": gene_id_ok, "replicates": reps, "qc": qc,
            "files": file_metadata(obj), "linked_control_accessions": sorted(linked),
            "control_metadata": control_records, "control_metadata_pass": matched_control,
            "treatment_quantification_annotations": {a: sorted(v) for a, v in treatment_ann.items()},
            "pooled_idr_file_metadata_candidates": pooled_idr,
            "metadata_resolution_errors": obj.get("_metadata_resolution_errors", []),
            "response_metadata_pass": bool(response_pass), "binding_metadata_candidate": bool(binding_candidate),
            "screen_decision": status, "previously_inspected_outcomes_excluded": identity(obj) in excluded,
            "source_disjointness_verified": False, "knockdown_efficacy_verified": False,
            "biological_validation_completed": False, "outcome_payloads_opened": 0}


def search_queries(protocol: dict) -> list[dict]:
    cfg = protocol.get("search", {})
    shared = {"type": "Experiment", "status": "released", "biosample_ontology.term_name": "K562"}
    queries = [{**shared, "assay_title": name} for name in protocol.get("allowed_assay_families", ASSAYS)]
    # A broad assay-term fallback catches records with older assay_title vocabulary.
    if cfg.get("include_assay_term_fallbacks", True):
        queries.extend([{**shared, "assay_term_name": term} for term in ("RNA-seq", "ChIP-seq")])
    return queries


def discover(client: MetadataClient, protocol: dict) -> tuple[dict[str, dict], list[dict]]:
    cfg = protocol.get("search", {})
    page_size = int(cfg.get("page_size", 100))
    max_pages = int(cfg.get("max_pages", 25))
    experiments: dict[str, dict] = {}
    searches = []
    for query in search_queries(protocol):
        seen: set[str] = set()
        trace = {"query": query, "urls": [], "complete": False, "errors": [], "reported_total": None}
        for page in range(max_pages):
            params = {**query, "limit": page_size, "from": page * page_size, "format": "json", "frame": "embedded"}
            url = HOST + "/search/?" + urllib.parse.urlencode(params)
            trace["urls"].append(url)
            try:
                obj = client.get(url)
                graph = obj.get("@graph")
                total = obj.get("total")
                if not isinstance(graph, list) or not isinstance(total, int):
                    raise ValueError("Search response lacks @graph or integer total")
                if trace["reported_total"] is None:
                    trace["reported_total"] = total
                elif trace["reported_total"] != total:
                    raise ValueError("Search total changed during pagination")
                before = len(seen)
                for record in graph:
                    accession = identity(record)
                    if not re.fullmatch(r"ENCSR[A-Z0-9]{6}", accession):
                        raise ValueError("Search result lacks experiment accession")
                    seen.add(accession)
                    experiments[accession] = record
                if len(seen) >= total:
                    trace["complete"] = True
                    break
                if not graph or len(seen) == before:
                    raise ValueError("Search pagination failed to advance before declared total")
            except Exception as exc:
                trace["errors"].append(f"{type(exc).__name__}: {exc}")
                break
        if not trace["complete"] and not trace["errors"]:
            trace["errors"].append("Search page limit reached; coverage is incomplete")
        trace["unique_returned"] = len(seen)
        searches.append(trace)
    return experiments, searches


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items() if k in fields})


def run(protocol: dict, output: Path, cache: Path | None = None) -> dict:
    regulators = protocol["regulators"]
    if len(regulators) != 16 or len(set(regulators)) != 16 or not all(isinstance(x, str) for x in regulators):
        raise ValueError("BIO-V8 requires the exact frozen list of 16 unique TF symbols")
    if set(protocol.get("allowed_assay_families", ASSAYS)) - set(ASSAYS):
        raise ValueError("Protocol contains an unapproved assay family")
    output.mkdir(parents=True, exist_ok=True)
    client = MetadataClient(output, cache, int(protocol.get("fetch_timeout_seconds", 30)))
    discovered, searches = discover(client, protocol)
    results, errors, rejected_search_records = [], [], []
    expanded_cache: dict[str, dict] = {}

    def original(accession: str) -> dict:
        if accession not in expanded_cache:
            obj = client.get(record_url(accession, "experiments"))
            if identity(obj) != accession:
                raise ValueError("Original experiment identity mismatch")
            expanded_cache[accession] = hydrate_experiment(obj, client)
        return expanded_cache[accession]

    for accession, stub in sorted(discovered.items()):
        # Resolve canonical target identities before expensive replicate/file metadata.
        # The cache reuses a target/gene record across all experiments for that TF.
        stub = hydrate_experiment(stub, client, include_files=False, identity_only=True)
        if stub.get("_metadata_resolution_errors"):
            errors.append({"accession": accession, "stage": "search_target_identity_resolution",
                           "error": stub["_metadata_resolution_errors"]})
        ti = target_identity(stub)
        possible = set(ti["canonical_symbols"] + ti["target_labels"]) & set(regulators)
        # Unresolved targets need original JSON; known different canonical genes do not.
        if ti["canonical_symbols"] and not possible:
            continue
        if not possible and assay_family(stub) == "TF ChIP-seq" and not ti["target_labels"]:
            # Input/control ChIP experiments without a target are not TF assays.
            continue
        try:
            obj = original(accession)
            ti = target_identity(obj)
            possible = set(ti["canonical_symbols"] + ti["target_labels"]) & set(regulators)
            if not possible:
                rejected_search_records.append({"accession": accession, "canonical_symbols": ti["canonical_symbols"],
                                                "reason": "no_exact_frozen_regulator_match", "resolution_errors": obj.get("_metadata_resolution_errors", [])})
                continue
            controls = []
            for value in list_value(obj.get("possible_controls")) + list_value(obj.get("controls")):
                cid = identity(value)
                if not re.fullmatch(r"ENCSR[A-Z0-9]{6}", cid):
                    errors.append({"accession": accession, "error": "Unresolvable control reference " + cid})
                    continue
                try:
                    controls.append(original(cid))
                except Exception as exc:
                    errors.append({"accession": accession, "control": cid, "error": str(exc)})
            results.append(inspect_experiment(obj, controls, regulators, protocol))
        except Exception as exc:
            errors.append({"accession": accession, "error": f"{type(exc).__name__}: {exc}"})
    complete = bool(searches) and all(s["complete"] for s in searches) and not errors
    # Any unresolved target metadata invalidates a global claim of complete coverage.
    unresolved = any(x.get("metadata_resolution_errors") for x in results) or any(x.get("resolution_errors") for x in rejected_search_records)
    complete = complete and not unresolved
    coverage = []
    for tf in regulators:
        exact_rows = [r for r in results if tf in r["matched_regulators"]]
        ambiguous = [r for r in results if tf in r["label_only_candidates"] and not r["target_identity_pass"]]
        response = [r["accession"] for r in exact_rows if r["response_metadata_pass"]]
        binding = [r["accession"] for r in exact_rows if r["binding_metadata_candidate"]]
        coverage.append({"tf": tf, "exact_canonical_experiments": len(exact_rows),
                         "experiment_accessions": [r["accession"] for r in exact_rows],
                         "ambiguous_label_candidates": [r["accession"] for r in ambiguous],
                         "response_metadata_pass_accessions": response, "binding_metadata_candidate_accessions": binding,
                         "coverage_state": "metadata_candidates_found" if exact_rows else ("no_candidate_in_completed_queries" if complete else "unknown_search_or_resolution_incomplete"),
                         "measured_candidate_target_count": None, "minimum_target_coverage_evaluated": False})
    response_tfs = [r["tf"] for r in coverage if r["response_metadata_pass_accessions"]]
    binding_tfs = [r["tf"] for r in coverage if r["binding_metadata_candidate_accessions"]]
    minimum = int(protocol.get("min_panel_tfs", 6))
    panel_pass = complete and len(response_tfs) >= minimum
    summary = {"protocol_id": protocol["protocol_id"], "completed_at_utc": utc_now(), "phase": "BIO-V8 metadata-only qualification",
               "search_complete": complete, "searches": searches, "unique_search_experiments": len(discovered),
               "exact_regulators": regulators, "min_panel_tfs": minimum, "min_targets_per_tf": int(protocol.get("min_targets_per_tf", 20)),
               "response_metadata_pass_tfs": response_tfs, "binding_metadata_candidate_tfs": binding_tfs,
               "response_panel_metadata_threshold_pass": panel_pass,
               "decision": "METADATA_PANEL_CANDIDATE_PENDING_PROVENANCE_QC_AND_TARGET_COVERAGE" if panel_pass else "NO_GO_FOR_BIOLOGICAL_SCORING",
               "errors": errors, "unresolved_metadata_present": unresolved, "request_count": len(client.ledger),
               "failed_requests": sum(x["status"] == "error" for x in client.ledger),
               "new_assay_outcome_payloads_retrieved": 0, "test_cells_accessed": 0,
               "independent_biological_validation_completed": False, "source_reference_independence_verified": False,
               "target_coverage_and_knockdown_efficacy_verified": False,
               "interpretation": "Metadata eligibility only. Search absence is scoped to recorded queries. Uninspected outcomes remain sealed; binding metadata are not complete regulatory truth."}
    for name, value in (("summary.json", summary), ("experiments.json", results), ("request_ledger.json", client.ledger),
                        ("search_records_excluded.json", rejected_search_records), ("tf_coverage.json", coverage)):
        (output / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    write_csv(output / "experiments.csv", results, ["accession", "assay_family", "status", "matched_regulators", "label_only_candidates",
              "target_identity_pass", "tcf4_alias_ambiguity", "response_metadata_pass", "binding_metadata_candidate", "screen_decision",
              "replicates", "qc", "linked_control_accessions", "control_metadata_pass", "metadata_resolution_errors"])
    write_csv(output / "tf_coverage.csv", coverage, list(coverage[0]))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    result = run(protocol, args.output, args.cache)
    print(json.dumps({k: result[k] for k in ("decision", "search_complete", "response_metadata_pass_tfs", "binding_metadata_candidate_tfs", "failed_requests")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
