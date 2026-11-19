import pandas as pd
from scripts.audit_m2d15b_exact_concordance import audit_records

def test_exact_join_with_headerless_first_record_preserved():
    m=pd.DataFrame([["AAAA","TF1"],["BBBB","TF2"],["CCCC",None]])
    obs=pd.DataFrame({"TF":["TF1","OTHER","TF0"]},index=["AAAA","BBBB","DDDD"])
    r=audit_records(m,obs)
    assert r["n_tfmap_records_headerless"]==3
    assert r["exact_matched_cells"]==2
    assert r["exact_concordant_tf_labels"]==1
    assert r["exact_discordant_tf_labels"]==1
    assert r["exact_unmatched_h5ad"]==1
    assert r["exact_unmatched_tfmap"]==1
    assert r["sample_aware_join"] is False
    assert r["barcode_transform_applied"] is False

def test_no_suffix_rewriting():
    m=pd.DataFrame([["AAAA","TF1"]])
    obs=pd.DataFrame({"TF":["TF1"]},index=["AAAA-1"])
    r=audit_records(m,obs)
    assert r["exact_matched_cells"]==0
    assert r["exact_concordance_fraction"] is None
