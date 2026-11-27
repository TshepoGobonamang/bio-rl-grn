import hashlib
import json
import numpy as np
import pandas as pd
import pytest
from scripts.biov6_two_evidence_adapter import (
    load_two_biov6_evidence_sources,MOTIF_NAME,DNASE_NAME,EXPECTED_MOTIF_SHA
)


def original(tmp_path):
    tfs=[f"TF{i}" for i in range(16)]
    genes=tfs+[f"G{i}" for i in range(48)]
    pairs=[(a,b) for a in tfs for b in genes if a!=b]
    df=pd.DataFrame({"tf":[a for a,b in pairs],"target":[b for a,b in pairs],
        "motif_prior_score":[.75 if i%2==0 else 0.0 for i in range(1008)],
        "dnase_accessibility_mean":[.5]*1008,
        "joint_motif_dnase_product":[.375 if i%2==0 else 0.0 for i in range(1008)],
        "accessibility_source_observed":[True]*1008})
    p=tmp_path/"k562_tf_target_accessibility_and_motif_context_1008.csv"
    df.to_csv(p,index=False)
    m={"stage":"BIOV6_REPLICATE_DNASE_ACCESSIBILITY_CONTEXT_NO_TF_BINDING_GOLD",
        "n_promoters":64,"n_directed_edges":1008,
        "original_biov5_motif_prior_sha256":EXPECTED_MOTIF_SHA,
        "ENCODE_perturbation_expression_outcomes_accessed":False,
        "GCEV_predictions_fitted":False,
        "joint_prior_is_independent_third_source":False,
        "source_bed_file_summaries":[{},{}],
        "file_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()}}
    (tmp_path/"biov6_source_audit.json").write_text(json.dumps(m))
    return pairs,df,p,m


def test_two_sources_not_three_and_known_zero_not_missing(tmp_path):
    pairs,df,p,m=original(tmp_path)
    output=load_two_biov6_evidence_sources(tmp_path,list(reversed(pairs)))
    assert set(output)=={MOTIF_NAME,DNASE_NAME}
    a,mask=output[MOTIF_NAME]
    assert a[0]==0 and mask.all()
    assert output[DNASE_NAME][1].all()
    assert np.allclose(output[DNASE_NAME][0],.5)


def test_mismatched_candidate_universe_rejected(tmp_path):
    pairs,*_=original(tmp_path)
    with pytest.raises(ValueError,match="Mismatch"):
        load_two_biov6_evidence_sources(tmp_path,pairs[:-1]+[("UNLISTED","G1")])


def test_forged_third_source_or_outcome_leakage_blocks(tmp_path):
    pairs,df,p,m=original(tmp_path)
    m["joint_prior_is_independent_third_source"]=True
    (tmp_path/"biov6_source_audit.json").write_text(json.dumps(m))
    with pytest.raises(ValueError,match="not qualified"):
        load_two_biov6_evidence_sources(tmp_path,pairs)


def test_original_file_hash_stays_mandatory(tmp_path):
    pairs,df,p,m=original(tmp_path)
    df.loc[0,"dnase_accessibility_mean"]=.7
    df.to_csv(p,index=False)
    with pytest.raises(ValueError,match="checksum"):
        load_two_biov6_evidence_sources(tmp_path,pairs)
