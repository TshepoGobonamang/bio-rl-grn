import pytest
from scripts.biov6_matched_caller_sensitivity import PIPELINES,select_matched


def example():
    return {
      "eligible_files_by_replicate":{
        "1":[{"accession":"ENCFF821KDJ","pipeline":"/analysis-step-versions/dnase-call-hotspots-se-step-v-2-1/","biosample_id":"B1"},
             {"accession":"ENCFF185XRG","pipeline":"/analysis-step-versions/dnase-peak-calling-step-v-1-0/","biosample_id":"B1"}],
        "2":[{"accession":"ENCFF433TIR","pipeline":"/analysis-step-versions/dnase-call-hotspots-se-step-v-2-1/","biosample_id":"B2"},
             {"accession":"ENCFF807ICZ","pipeline":"/analysis-step-versions/dnase-peak-calling-step-v-1-0/","biosample_id":"B2"},
             {"accession":"ENCFF784DIM","pipeline":"/analysis-step-versions/dnase-peak-calling-step-v-1-0/","biosample_id":"B2"}]
      }
    }


def test_two_pipeline_versions_fixed_and_no_posthoc_choice():
    m=example()
    assert select_matched(m,"v2")[1]["accession"]=="ENCFF821KDJ"
    assert select_matched(m,"v2")[2]["accession"]=="ENCFF433TIR"
    assert select_matched(m,"v1")[1]["accession"]=="ENCFF185XRG"
    assert select_matched(m,"v1")[2]["accession"]=="ENCFF784DIM"


def test_missing_rep_does_not_get_silently_substituted():
    m=example()
    m["eligible_files_by_replicate"]["2"]=[f for f in m["eligible_files_by_replicate"]["2"] if PIPELINES["v2"] not in f["pipeline"]]
    with pytest.raises(ValueError,match="NO-GO"):
        select_matched(m,"v2")


def test_reused_biosample_rejected():
    m=example()
    for f in m["eligible_files_by_replicate"]["2"]:
        f["biosample_id"]="B1"
    with pytest.raises(ValueError,match="Same biosample"):
        select_matched(m,"v1")
