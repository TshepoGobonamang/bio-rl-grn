import pandas as pd
import pytest
from scripts.reconcile_m2d15_tfmap_h5ad import reconcile

def test_exact_cell_and_sample_join():
    obs=pd.DataFrame({"cell":["AA","AA","BB"],"batch":["S1","S2","S1"],"TF":["SOX2","NANOG","PAX6"]})
    mapping=pd.DataFrame({"barcode":["AA","AA","BB","CC"],"sample":["S1","S2","S1","S1"],"ORF":["SOX2","NANOG","OCT4","EMX2"]})
    r=reconcile(obs,mapping,"cell","barcode","TF","ORF","batch","sample")
    assert (r["matched_keys"],r["unmatched_h5ad"],r["concordant_labels"],r["discordant_labels"],r["unmatched_tfmap"])==(3,0,2,1,1)
    assert not r["matched_control_verified"]

def test_reject_ambiguous_tf_map():
    obs=pd.DataFrame({"c":["A"],"tf":["TF1"]})
    mapping=pd.DataFrame({"c":["A","A"],"tf":["TF1","TF2"]})
    with pytest.raises(ValueError,match="Conflicting"):
        reconcile(obs,mapping,"c","c","tf","tf")

def test_identical_duplicates_are_counted():
    obs=pd.DataFrame({"c":["A"],"tf":["TF1"]})
    mapping=pd.DataFrame({"c":["A","A"],"tf":["TF1","TF1"]})
    r=reconcile(obs,mapping,"c","c","tf","tf")
    assert r["tfmap_duplicate_same_key_rows_dropped"]==1
    assert r["concordant_labels"]==1

def test_missing_sample_pair_rejected():
    obs=pd.DataFrame({"c":["A"],"s":["S1"],"tf":["TF1"]})
    mapping=pd.DataFrame({"c":["A"],"tf":["TF1"]})
    with pytest.raises(ValueError,match="Sample keys"):
        reconcile(obs,mapping,"c","c","tf","tf","s",None)

def test_no_implicit_barcode_normalization():
    obs=pd.DataFrame({"c":["AA-1"],"tf":["TF1"]})
    mapping=pd.DataFrame({"c":["AA"],"tf":["TF1"]})
    r=reconcile(obs,mapping,"c","c","tf","tf")
    assert r["matched_keys"]==0
    assert r["unmatched_h5ad"]==1
