"""Meaningful fail-closed metadata and blinded-transport tests; no live ENCODE."""
import copy
import importlib.util
import io
import json
import urllib.error
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/biov8_encode_metadata_audit.py"
SPEC = importlib.util.spec_from_file_location("biov8_metadata", PATH)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def experiment(accession="ENCSR123ABC", symbol="RUNX1", bios=("ENCBS111AAA", "ENCBS222BBB")):
    x = {"accession": accession, "status": "released", "assay_title": "CRISPRi RNA-seq",
         "assay_term_name": "RNA-seq", "biosample_ontology": {"term_name": "K562"},
         "replication_type": "isogenic", "audit": {}, "files": [],
         "replicates": [{"biological_replicate_number": i + 1, "technical_replicate_number": 1,
                          "library": {"biosample": {"accession": b}}} for i, b in enumerate(bios)]}
    if symbol:
        x["target"] = {"label": symbol, "genes": [{"symbol": symbol, "geneid": "ENSG00000159216"}]}
        x["possible_controls"] = ["/experiments/ENCSR456DEF/"]
    else:
        x["control_type"] = "non-specific target"
    for i in (1, 2):
        x["files"].append({"accession": f"ENCFF000AA{i}", "status": "released", "file_format": "tsv",
                           "output_type": "gene quantifications", "assembly": "GRCh38", "genome_annotation": "V29",
                           "biological_replicates": [i], "audit": {}})
    return x


def control():
    return experiment("ENCSR456DEF", None, ("ENCBS333CCC", "ENCBS444DDD"))


def protocol():
    return {"assembly": "GRCh38", "expected_genome_annotation": "V29",
            "regulator_ensembl_ids": {"RUNX1": "ENSG00000159216"}}


def inspect(x=None, c=None, p=None):
    return audit.inspect_experiment(x or experiment(), [c or control()], ["RUNX1"], p or protocol())


def test_pass_is_metadata_only_not_biological_validation():
    r = inspect()
    assert r["response_metadata_pass"]
    assert r["screen_decision"] == "RESPONSE_METADATA_SCREEN_PASS_OUTCOMES_SEALED"
    assert not r["biological_validation_completed"]
    assert not r["source_disjointness_verified"]
    assert r["outcome_payloads_opened"] == 0


def test_tcf4_label_must_not_be_misidentified_tcf7l2():
    x = experiment(symbol="TCF4")
    x["target"]["genes"] = [{"symbol": "TCF7L2", "geneid": "ENSG00000148737"}]
    p = {**protocol(), "regulator_ensembl_ids": {"TCF4": "ENSG00000196628"}}
    r = audit.inspect_experiment(x, [control()], ["TCF4"], p)
    assert r["tcf4_alias_ambiguity"]
    assert r["canonical_symbols"] == ["TCF7L2"]
    assert not r["target_identity_pass"]
    assert not r["response_metadata_pass"]


def test_unresolved_tcf4_label_cannot_pass_identity():
    x = experiment(symbol="TCF4")
    x["target"].pop("genes")
    r = audit.target_identity(x)
    assert r["tcf4_alias_ambiguity"]
    assert not r["canonical_identity_verified"]


def test_wrong_gene_id_blocks_even_correct_label_and_symbol():
    x = experiment()
    x["target"]["genes"][0]["geneid"] = "ENSG00000148737"
    r = inspect(x)
    assert not r["canonical_gene_id_match"]
    assert not r["target_identity_pass"]


def test_duplicate_biosample_is_not_two_biological_replicates():
    r = inspect(experiment(bios=("ENCBS111AAA", "ENCBS111AAA")))
    assert not r["replicates"]["independent_biosample_metadata_pass"]
    assert not r["response_metadata_pass"]


def test_technical_duplicate_does_not_invalidate_two_real_biological_units():
    x = experiment()
    extra = copy.deepcopy(x["replicates"][0])
    extra["technical_replicate_number"] = 2
    x["replicates"].append(extra)
    assert inspect(x)["response_metadata_pass"]


def test_one_rep_number_cannot_be_reassigned_to_multiple_biosamples():
    x = experiment()
    x["replicates"][1]["biological_replicate_number"] = 1
    assert not inspect(x)["response_metadata_pass"]


def test_replicate_or_file_from_wrong_experiment_is_not_borrowed():
    x = experiment()
    x["replicates"][0]["experiment"] = "/experiments/ENCSR999ZZZ/"
    assert not inspect(x)["response_metadata_pass"]
    x = experiment()
    x["files"][0]["dataset"] = "/experiments/ENCSR999ZZZ/"
    assert not inspect(x)["response_metadata_pass"]


def test_unlinked_or_targeted_control_blocks():
    x = experiment()
    x["possible_controls"] = ["/experiments/ENCSR999ZZZ/"]
    assert not inspect(x)["response_metadata_pass"]
    c = control()
    c["target"] = {"label": "RUNX1", "genes": [{"symbol": "RUNX1"}]}
    assert not inspect(c=c)["response_metadata_pass"]


def test_preinspected_control_cannot_be_fresh():
    p = {**protocol(), "excluded_outcome_experiments": ["ENCSR456DEF"]}
    r = inspect(p=p)
    assert not r["control_metadata"][0]["fresh_control_metadata"]
    assert not r["response_metadata_pass"]


def test_annotation_mismatch_and_pooled_quantification_block():
    c = control()
    for f in c["files"]:
        f["genome_annotation"] = "V24"
    assert not inspect(c=c)["response_metadata_pass"]
    x = experiment()
    for f in x["files"]:
        f["biological_replicates"] = [1, 2]
    assert not inspect(x)["response_metadata_pass"]


def test_error_or_not_compliant_audit_fails_closed():
    for level in ("ERROR", "NOT_COMPLIANT"):
        x = experiment()
        x["audit"] = {level: [{"category": "unresolved concern"}]}
        assert not inspect(x)["response_metadata_pass"]
        c = control()
        c["files"][0]["audit"] = {level: [{"category": "file concern"}]}
        assert not inspect(c=c)["response_metadata_pass"]


def test_missing_audit_or_metadata_resolution_is_not_silent_pass():
    x = experiment()
    x.pop("audit")
    assert not inspect(x)["response_metadata_pass"]
    x = experiment()
    x["_metadata_resolution_errors"] = ["HTTP 403 resolving a source"]
    assert not inspect(x)["response_metadata_pass"]


def test_target_in_biosample_genetic_modification_is_found():
    x = experiment()
    target = x.pop("target")
    for rep in x["replicates"]:
        rep["library"]["biosample"]["applied_modifications"] = [{"modified_site_by_target_id": target}]
    r = inspect(x)
    assert r["target_identity_pass"]
    assert any("applied_modifications" in p for p in r["target_metadata_paths"])


@pytest.mark.parametrize("url", [
    "https://www.encodeproject.org/files/ENCFF123ABC/@@download/ENCFF123ABC.tsv?format=json",
    "https://www.encodeproject.org/files/ENCFF123ABC/%40%40download/?format=json",
    "https://www.encodeproject.org/files/ENCFF123ABC/%2540%2540download/?format=json",
    "https://www.encodeproject.org/files/ENCFF123ABC/",
    "https://www.encodeproject.org/files/ENCFF123ABC/?format=tsv",
    "https://www.encodeproject.org/files/ENCFF123ABC/?format=json&download=true",
    "https://www.encodeproject.org.attacker.example/files/ENCFF123ABC/?format=json",
    "https://user:password@www.encodeproject.org/files/ENCFF123ABC/?format=json",
    "http://www.encodeproject.org/experiments/ENCSR123ABC/?format=json",
    "https://www.encodeproject.org/search/../files/ENCFF123ABC/?format=json",
])
def test_metadata_allowlist_blocks_sealed_and_foreign_urls(url):
    with pytest.raises(ValueError):
        audit.validate_metadata_url(url)


def test_metadata_allowlist_accepts_file_record_but_not_payload():
    url = "https://www.encodeproject.org/files/ENCFF123ABC/?format=json&frame=embedded"
    assert audit.validate_metadata_url(url) == url


def test_search_failure_and_truncation_are_never_reported_as_absence():
    class FailingClient:
        def get(self, url):
            raise ConnectionError("Network unavailable")
    _, traces = audit.discover(FailingClient(), {"allowed_assay_families": ["TF ChIP-seq"],
                                               "search": {"include_assay_term_fallbacks": False}})
    assert len(traces) == 1
    assert not traces[0]["complete"]
    assert traces[0]["reported_total"] is None
    class TruncatedClient:
        def get(self, url):
            return {"total": 2, "@graph": [{"accession": "ENCSR123ABC"}]}
    _, traces = audit.discover(TruncatedClient(), {"allowed_assay_families": ["TF ChIP-seq"],
                        "search": {"include_assay_term_fallbacks": False, "max_pages": 1}})
    assert not traces[0]["complete"]
    assert traces[0]["unique_returned"] == 1


def test_empty_original_search_json_can_establish_scoped_absence():
    class EmptyClient:
        def get(self, url):
            return {"total": 0, "@graph": []}
    records, traces = audit.discover(EmptyClient(), {"allowed_assay_families": ["TF ChIP-seq"],
                                      "search": {"include_assay_term_fallbacks": False}})
    assert records == {}
    assert traces[0]["complete"]
    assert traces[0]["reported_total"] == 0


def test_original_404_search_requires_explicit_empty_json(tmp_path):
    url = "https://www.encodeproject.org/search/?type=Experiment&format=json"
    class Empty404:
        def open(self, req, timeout):
            raise urllib.error.HTTPError(url, 404, "Not found", {"Content-Type": "application/json"},
                                         io.BytesIO(json.dumps({"total": 0, "@graph": []}).encode()))
    client = audit.MetadataClient(tmp_path)
    client.opener = Empty404()
    assert client.get(url)["total"] == 0
    assert client.ledger[0]["http_status"] == 404
    assert client.ledger[0]["status"] == "ok"
    class Error404:
        def open(self, req, timeout):
            raise urllib.error.HTTPError(url, 404, "Not found", {"Content-Type": "application/json"},
                                         io.BytesIO(b'{"status":"error","detail":"query failed"}'))
    client2 = audit.MetadataClient(tmp_path / "failed")
    client2.opener = Error404()
    with pytest.raises(ValueError, match="did not establish"):
        client2.get(url)
    assert client2.ledger[0]["status"] == "error"


def test_binding_metadata_never_becomes_regulatory_truth():
    x = experiment()
    x["assay_title"] = "TF ChIP-seq"
    x["assay_term_name"] = "ChIP-seq"
    x["files"] = [{"accession": "ENCFF999AAA", "status": "released", "assembly": "GRCh38",
                    "file_format": "bed", "output_type": "optimal IDR thresholded peaks",
                    "biological_replicates": [1, 2], "audit": {}}]
    c = control()
    c.update({"assay_title": "Control ChIP-seq", "assay_term_name": "ChIP-seq", "control_type": "input library"})
    r = inspect(x, c)
    assert r["binding_metadata_candidate"]
    assert not r["response_metadata_pass"]
    assert not r["biological_validation_completed"]
    assert r["screen_decision"] == "BINDING_METADATA_CANDIDATE_ONLY"
