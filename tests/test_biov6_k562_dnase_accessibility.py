import hashlib
import json
import numpy as np
import pandas as pd
import pytest
from scripts.biov6_k562_dnase_accessibility import (
    inspect_dnase_original, covered_bases, compute_joint, EXPERIMENT
)


def source_experiment():
    result={"accession":EXPERIMENT,"status":"released",
            "assay_title":"DNase-seq","biosample_summary":"Homo sapiens K562",
            "replication_type":"isogenic",
            "replicates":[
                {"biological_replicate_number":1,"library":{"biosample":{"accession":"ENCBS001ABC"}}},
                {"biological_replicate_number":2,"library":{"biosample":{"accession":"ENCBS002ABC"}}}
            ],"files":[]}
    for i in (1,2):
        acc=f"ENCFF{i:06d}"
        result["files"].append({"accession":acc,"status":"released",
            "assembly":"GRCh38","file_format":"bed","file_format_type":"narrowPeak",
            "output_type":"peaks","biological_replicates":[i],
            "md5sum":"a"*32,"file_size":2222,
            "href":f"/files/{acc}/@@download/{acc}.bed.gz"})
    return result


def test_two_source_replicates_not_assumed_from_two_files():
    x=inspect_dnase_original(source_experiment(),"f"*64)
    assert x["selected_files_by_replicate"]["1"]["biosample_id"]=="ENCBS001ABC"
    assert x["selected_files_by_replicate"]["2"]["biosample_id"]=="ENCBS002ABC"
    assert x["source_file_discovery_pass"]


def test_one_biosample_across_two_replicates_fails():
    obj=source_experiment()
    obj["replicates"][1]["library"]["biosample"]["accession"]="ENCBS001ABC"
    with pytest.raises(ValueError,match="duplicated"):
        inspect_dnase_original(obj,"f"*64)


def test_nonmatching_assembly_or_assay_blocks():
    obj=source_experiment()
    obj["files"][1]["assembly"]="hg19"
    with pytest.raises(ValueError,match="NO-GO"):
        inspect_dnase_original(obj,"f"*64)
    obj=source_experiment()
    obj["assay_title"]="ChIP-seq"
    with pytest.raises(ValueError,match="assay identity"):
        inspect_dnase_original(obj,"f"*64)


def test_unions_never_double_count_overlapping_peak_intervals():
    # These intervals are assumed merged by upstream BED parser.
    merged={"chr1":[(100,180),(200,240),(800,850)]}
    assert covered_bases(merged,"chr1",120,230)==(180-120)+(230-200)
    assert covered_bases(merged,"chr1",100,900)==80+40+50
    assert covered_bases(merged,"chr2",100,900)==0
    assert covered_bases(merged,"chr1",240,800)==0


def test_joint_prior_is_context_not_independent_third_source():
    tfs=[f"TF{i}" for i in range(16)]
    genes=tfs+[f"G{i}" for i in range(48)]
    edges=[(tf,g) for tf in tfs for g in genes if tf!=g]
    df=pd.DataFrame({"tf":[p[0] for p in edges],"target":[p[1] for p in edges],
       "tf_gencode_gene_id":["ENSG1"]*1008,
       "target_gencode_gene_id":["ENSG2"]*1008,
       "prior_score":[.8]*1008,"empirical_p":[.2]*1008,"bh_q":[.6]*1008})
    promoter=pd.DataFrame({"gene_symbol":genes,
       "dnase_replicate1_coverage":[.3]*64,
       "dnase_replicate2_coverage":[.5]*64,
       "dnase_accessibility_mean":[.4]*64,
       "dnase_accessibility_min":[.3]*64,
       "dnase_both_replicates_peak_present":[True]*64,
       "dnase_replicate_abs_disagreement":[.2]*64})
    results=compute_joint(df,promoter)
    assert len(results)==1008
    assert np.allclose(results.joint_motif_dnase_product,.32)
    assert results.dnase_accessibility_mean.nunique()==1
    assert not results.direct_tf_binding_verified.any()
    assert not results.independent_regulatory_edge_truth.any()
    assert results.accessibility_source_observed.all()
    assert (results.groupby("target").dnase_accessibility_mean.nunique()==1).all()


def test_unknown_promoter_fails_not_imputed_as_zero():
    tfs=[f"TF{i}" for i in range(16)]
    genes=tfs+[f"G{i}" for i in range(48)]
    edges=[(tf,g) for tf in tfs for g in genes if tf!=g]
    df=pd.DataFrame({"tf":[a for a,b in edges],"target":[b for a,b in edges],
        "tf_gencode_gene_id":["ENSG1"]*1008,"target_gencode_gene_id":["ENSG2"]*1008,
        "prior_score":[0.0]*1008,"empirical_p":[.2]*1008,"bh_q":[.8]*1008})
    p=pd.DataFrame({"gene_symbol":["TF0"]})
    with pytest.raises(ValueError,match="Missing promoter"):
        compute_joint(df,p)
