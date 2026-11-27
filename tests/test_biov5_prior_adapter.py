import hashlib
import json
import numpy as np
import pandas as pd
import pytest

from scripts.biov5_prior_adapter import SOURCE_NAME, load_biov5_sequence_prior


def fixture(tmp_path):
    tfs=[f"TF{i}" for i in range(16)]
    genes=tfs+[f"G{i}" for i in range(48)]
    edges=[(tf,g) for tf in tfs for g in genes if tf!=g]
    assert len(edges)==1008
    path=tmp_path/"k562_tf_to_promoter_motif_prior_1008.csv"
    table=pd.DataFrame({"tf":[a for a,b in edges],
                        "target":[b for a,b in edges],
                        "prior_score":[0.0 if i%2 else .5 for i in range(len(edges))],
                        "source_observed":[True]*len(edges),
                        "empirical_p":[.25]*len(edges),
                        "bh_q":[.8]*len(edges)})
    table.to_csv(path,index=False)
    manifest={
       "stage":"BIOV5_IN_SILICO_SEQUENCE_BINDING_POTENTIAL_NO_BIOLOGICAL_GOLD",
       "tested_directed_no_self_edges":1008,
       "ENCODE_perturbation_response_labels_accessed":False,
       "gcev_model_fitted":False,
       "file_sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest()}
    }
    (tmp_path/"biov5_motif_prior_manifest.json").write_text(json.dumps(manifest))
    return edges,path


def test_observed_zero_not_missing_and_order_reindexed(tmp_path):
    edges,_=fixture(tmp_path)
    inverse=list(reversed(edges))
    result=load_biov5_sequence_prior(tmp_path,inverse)
    values,mask=result[SOURCE_NAME]
    assert len(values)==1008 and mask.all()
    assert values[0]==0.0 and values[1]==.5
    assert mask.dtype==bool


def test_missing_extra_or_repeated_directed_pairs_fail(tmp_path):
    edges,_=fixture(tmp_path)
    with pytest.raises(ValueError,match="candidate identities"):
        load_biov5_sequence_prior(tmp_path,edges[:-1]+[("UNMATCHED","G1")])


def test_mutated_source_prior_checksum_fails(tmp_path):
    edges,path=fixture(tmp_path)
    x=pd.read_csv(path)
    x.loc[0,"prior_score"]=.9
    x.to_csv(path,index=False)
    with pytest.raises(ValueError,match="checksum"):
        load_biov5_sequence_prior(tmp_path,edges)


def test_unobserved_source_and_missing_dna_must_not_be_false_zero(tmp_path):
    edges,path=fixture(tmp_path)
    x=pd.read_csv(path)
    x.loc[1,"source_observed"]=False
    x.to_csv(path,index=False)
    manifest_path=tmp_path/"biov5_motif_prior_manifest.json"
    manifest=json.loads(manifest_path.read_text())
    manifest["file_sha256"][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match="Unobserved"):
        load_biov5_sequence_prior(tmp_path,edges)
