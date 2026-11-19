import copy
import pytest
from scripts.decide_m2d13_design_readiness import decide

OBS = {
 "n_cells":4, "source_sha256":"abc",
 "counts_by_field":{"TF":[{"value":"SOX2","cells":2},{"value":"GFP","cells":2}],
 "batch":[{"value":"A","cells":2},{"value":"B","cells":2}],
 "temp":[{"value":"1","cells":4}]},
 "tf_batch_temp_cells":[
 {"TF":"SOX2","batch":"A","temp":"1","cells":1},
 {"TF":"SOX2","batch":"B","temp":"1","cells":1},
 {"TF":"GFP","batch":"A","temp":"1","cells":1},
 {"TF":"GFP","batch":"B","temp":"1","cells":1}],
 "batch_count_per_tf":[{"TF":"SOX2","n_batches":2},{"TF":"GFP","n_batches":2}],
}
def test_unknown_experiment_is_blocked():
    r = decide(OBS, {})
    assert r["decision"] == "NO_GO"
    assert r["n_tfs_present_in_multiple_batches_not_replicates"] == 2
    assert not r["independent_culture_replication_pass"]

def test_control_name_alone_cannot_qualify():
    r = decide(OBS, {"validated_control_tf_labels":["GFP"]})
    assert r["decision"] == "NO_GO"
    assert not r["control_evidence_pass"]

def test_unsupported_control_fails():
    with pytest.raises(ValueError, match="actual D12 TF"):
        decide(OBS, {"validated_control_tf_labels":["FAKE"]})

def test_inconsistent_source_counts_fail():
    v=copy.deepcopy(OBS)
    v["tf_batch_temp_cells"][0]["cells"]=3
    with pytest.raises(ValueError, match="reconcile"):
        decide(v, {})

def test_fully_evidenced_status_is_only_design_review():
    e={"validated_control_tf_labels":["GFP"],"control_identity_verified":True,
       "control_identity_source":"publication:method-section",
       "biological_culture_replicate_provenance":{"verified":True,"source":"sample sheet","independent_culture_id_field":"culture_id","sequencing_batch_not_used_as_replicate":True},
       "tf_barcode_mapping_verified":True,"tf_barcode_mapping_source":"checksum-confirmed join",
       "matched_cell_context_verified":True,"matched_cell_context_source":"documented cell types"}
    r=decide(OBS,e)
    assert r["decision"]=="ELIGIBLE_FOR_FURTHER_DESIGN_REVIEW"
    assert not r["decision"]=="BIOLOGICALLY_VALIDATED"
