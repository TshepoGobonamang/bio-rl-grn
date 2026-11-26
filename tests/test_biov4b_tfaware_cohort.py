import numpy as np
import pytest
from scripts.biov4b_select_tfaware_cohort import choose


def candidate_data():
    genes=[f"TF{i}" for i in range(20)]+[f"GENE{i}" for i in range(60)]+["MT-CO1","RPL4","DUP"]
    detected=np.full(len(genes),.30)
    variance=np.array([1.0+0.01*i for i in range(len(genes))])
    annotations={}
    for n,i in zip(genes,range(len(genes))):
        annotations[n.casefold()]=[{
            "gencode_gene_id":f"ENSG{i:011d}.1",
            "gene_symbol":n,"gene_type":"protein_coding","chromosome":"chr1"}]
    annotations["dup"] += [{"gencode_gene_id":"ENSG99999999999.1",
                             "gene_symbol":"DUP","gene_type":"protein_coding",
                             "chromosome":"chr2"}]
    motif={f"tf{i}".casefold():[{"motif_id":f"MA{i:04d}.1","jaspar_factor_name":f"TF{i}"}]
           for i in range(20)}
    return genes,detected,variance,annotations,motif


def test_fit_only_tf_enrichment_is_deterministic_and_no_gold():
    args=candidate_data()
    t,tar,selected,stats=choose(*args)
    assert len(t)==16 and len(tar)==48
    assert stats["directed_candidate_edges_no_self"]==16*63
    assert t[0]["gene_symbol"]=="TF19"
    assert tar[0]["gene_symbol"]=="GENE59"
    assert len({x["gene_symbol"] for x in selected})==64
    assert all(x["gene_symbol"].startswith("TF") for x in t)
    assert all(x["gene_symbol"].startswith("GENE") for x in tar)


def test_preserves_original_list_and_no_replacement_when_missing_regulators():
    args=list(candidate_data())
    gene_copy=args[0][:]
    args[1][0:16]=0.0
    with pytest.raises(ValueError,match="Insufficient fit-only"):
        choose(*args)
    assert args[0]==gene_copy


def test_excludes_mito_ribosomal_and_ambiguous():
    t,tar,selected,_=choose(*candidate_data())
    assert not set(["MT-CO1","RPL4","DUP"]) & {s["gene_symbol"] for s in selected}


def test_eligible_tf_requires_motif_not_symbol_name():
    args=list(candidate_data())
    args[-1].clear()
    with pytest.raises(ValueError,match="Insufficient fit-only"):
        choose(*args)


def test_source_data_malformed_fails_closed():
    g,d,v,gc,j=candidate_data()
    with pytest.raises(ValueError,match="Duplicate"):
        choose(["A"]*len(g),d,v,gc,j)
    with pytest.raises(ValueError,match="aligned"):
        choose(g,d[:-1],v,gc,j)
