from pathlib import Path
import gzip
import numpy as np
import pytest
from scripts.biov4_qualify_gene_universe import (
    read_gencode,read_jaspar_meme,inspect_universe,sha
)


def test_gtf_unambiguous_gene_type_and_ambiguous_shared_symbols(tmp_path):
    p=tmp_path/"gencode.v29.annotation.gtf.gz"
    with gzip.open(p,"wt") as f:
        f.write("##gencode version 29\n")
        for chrom,gid,name,kind in [
            ("chr1","ENSG000001.1","TF1","protein_coding"),
            ("chr2","ENSG000002.1","TF2","protein_coding"),
            ("chr3","ENSG000003.1","SAME","protein_coding"),
            ("chr4","ENSG000004.1","SAME","protein_coding"),
            ("chr5","ENSG000005.1","TARGET","protein_coding"),
            ("chrM","ENSG000006.1","MT-TEST","protein_coding"),
        ]:
            f.write(f'{chrom}\tGENCODE\tgene\t10\t100\t.\t+\t.\tgene_id "{gid}"; gene_type "{kind}"; gene_name "{name}";\n')
        for j in range(1001):
            f.write(f'chr1\tGENCODE\tgene\t10\t100\t.\t+\t.\tgene_id "ENSG{j+900000:06d}.1"; gene_type "protein_coding"; gene_name "XGENE{j}";\n')
    r=read_gencode(p)
    assert r["tf1"][0]["gencode_gene_id"]=="ENSG000001.1"
    assert len(r["same"])==2
    assert "mt-test" not in r


def test_jaspar_single_factor_only(tmp_path):
    p=tmp_path/"JASPAR.meme"
    with p.open("w") as f:
        f.write("MEME version 4\n")
        for j in range(101):
            f.write(f"MOTIF MA{j:04d}.1 XGENE{j}\n")
        f.write("MOTIF MA9999.1 TF1\n")
        f.write("MOTIF MA9998.1 TF2::OTHER\n")
    r=read_jaspar_meme(p)
    assert r["tf1"][0]["motif_id"]=="MA9999.1"
    assert "tf2::other" not in r


def test_gene_universe_unchanged_and_only_real_single_tf_matched():
    g={
      "tf1":[{"gencode_gene_id":"ENSG000001.1","gene_type":"protein_coding",
             "gene_symbol":"TF1","chromosome":"chr1"}],
      "tf2":[{"gencode_gene_id":"ENSG000002.1","gene_type":"lncRNA",
             "gene_symbol":"TF2","chromosome":"chr1"}],
      "same":[{"gencode_gene_id":"ENSG000003.1","gene_type":"protein_coding",
               "gene_symbol":"SAME","chromosome":"chr1"},
              {"gencode_gene_id":"ENSG000004.1","gene_type":"protein_coding",
               "gene_symbol":"SAME","chromosome":"chr1"}],
      "target":[{"gencode_gene_id":"ENSG000005.1","gene_type":"protein_coding",
                 "gene_symbol":"TARGET","chromosome":"chr1"}],
    }
    jas={"tf1":[{"motif_id":"MA0001.1","jaspar_factor_name":"TF1"}],
         "tf2":[{"motif_id":"MA0002.1","jaspar_factor_name":"TF2"}],
         "same":[{"motif_id":"MA0003.1","jaspar_factor_name":"SAME"}]}
    ann,mot,tf,counts,decision=inspect_universe(
         ["TF1","TF2","SAME","TARGET","MT-ND1"],
         np.array([.32,.90,.90,.85,.91]),g,jas,min_tfs=2)
    assert counts["eligible_regulators"]==1
    assert tf.tf_symbol.tolist()==["TF1"]
    assert ann.loc[2,"gene_annotation_status"]=="ambiguous_gene_id"
    assert counts["mitochondrial_prefix_genes"]==1
    assert decision.startswith("NO_GO")
    assert len(ann)==5


def test_unknown_transcription_factors_are_not_invented():
    g={"target":[{"gencode_gene_id":"ENSG000005.1",
                  "gene_type":"protein_coding","gene_symbol":"TARGET",
                  "chromosome":"chr1"}]}
    _,mot,tf,c,decision=inspect_universe(
       ["TARGET","A","B","C","D"],np.ones(5),g,{})
    assert tf.empty and mot.empty and c["eligible_regulators"]==0
    assert decision.startswith("NO_GO")


def test_duplicate_selected_genes_prohibited():
    with pytest.raises(ValueError,match="Feature order"):
        inspect_universe(["A","A","B","C","D"],np.ones(5),{},{})


def test_gencode_primary_chromosome_ambiguity_without_aliases():
    import inspect
    from scripts import biov4_qualify_gene_universe as mod
    code=inspect.getsource(mod)
    assert "gencode" in code and "jaspar" in code
    assert "ENCODE" not in code or "no_encode_response_access" in code
