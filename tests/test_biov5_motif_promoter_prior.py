"""Tests for actual sequence-backed promoter prior; never rely on assay outcomes."""
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from scripts.biov5_scan_grch38_promoters import (
    PROMOTER_LENGTH, bh_q, dinucleotide_shuffle, encode,
    one_target_edges, promoter_coordinates, reverse_complement,
    scan_one_motif, parse_meme_selected
)


def gene(strand, start, end):
    return {"gene_symbol":"TGENE","gencode_gene_id":"ENSG0000001.1",
            "gene_start_1based":start,"gene_end_1based":end,
            "chromosome":"chr2","gencode_strand":strand}


def test_exact_gene_oriented_promoters_positive_and_negative():
    f=promoter_coordinates(gene("+",1201,9000))
    assert (f["tss_0based"],f["start_0based"],
            f["end_0based_exclusive"])==(1200,200,1401)
    r=promoter_coordinates(gene("-",1201,2001))
    assert (r["tss_0based"],r["start_0based"],
            r["end_0based_exclusive"])==(2000,1800,3001)
    assert f["length"]==r["length"]==PROMOTER_LENGTH
    assert reverse_complement("ACGTTGCA")=="TGCAACGT"


def test_promoter_boundary_blocks_implicit_shortening():
    with pytest.raises(ValueError,match="falls outside"):
        promoter_coordinates(gene("+",100,999))


def test_exact_dinucleotide_and_terminal_base_preservation():
    seq=("ACGTTGCGAATCCCAGCTTG"*60)+"A"
    assert len(seq)==PROMOTER_LENGTH
    rng=np.random.default_rng(9)
    nulls=[dinucleotide_shuffle(seq,rng) for _ in range(8)]
    for n in nulls:
        assert len(n)==len(seq)
        assert n[0]==seq[0] and n[-1]==seq[-1]
        assert Counter(zip(n,n[1:]))==Counter(zip(seq,seq[1:]))
    other=np.random.default_rng(9)
    assert nulls==[dinucleotide_shuffle(seq,other) for _ in range(8)]


def motif(model="MA0001.1"):
    x=np.tile(np.array([[.97,.01,.01,.01]]),(5,1))
    pwm=np.log2(x/.25)
    return {"motif_id":model,"width":5,"pwm_bits":pwm,
            "min_bits":float(pwm.min(axis=1).sum()),
            "max_bits":float(pwm.max(axis=1).sum())}


def test_both_strands_pwm_hits_have_gene_oriented_positions():
    m=motif()
    seq="C"*500+"TTTTT"+"C"*696
    assert len(seq)==PROMOTER_LENGTH
    hit=scan_one_motif(encode(seq),encode(reverse_complement(seq)),m)
    assert hit["relative"]>0.999
    assert hit["strand"]=="-"
    assert hit["position_0based"]==500
    fw="C"*500+"AAAAA"+"C"*696
    hit2=scan_one_motif(encode(fw),encode(reverse_complement(fw)),m)
    assert hit2["strand"]=="+" and hit2["position_0based"]==500


def test_null_statistics_and_no_self_loops_or_gold_reads():
    rows=pd.DataFrame([
        {"gene_symbol":"TF1","gencode_gene_id":"ENSG0000001.1","role":"regulator"},
        {"gene_symbol":"TF2","gencode_gene_id":"ENSG0000002.1","role":"regulator"},
        {"gene_symbol":"TARGET","gencode_gene_id":"ENSG0000003.1","role":"target_only"}])
    motifs={"TF1":[motif("MA0001.1")],
            "TF2":[motif("MA0002.1")]}
    seq=("ACGTCGATTCGACGGTTAAC"*60)+"A"
    result=one_target_edges("TF1",seq,rows,motifs)
    assert len(result)==1
    assert result[0]["tf"]=="TF2" and result[0]["target"]=="TF1"
    assert 0<=result[0]["prior_score"]<=1
    assert result[0]["null_shuffles"]==128
    assert 0<result[0]["empirical_p"]<=1
    assert not result[0]["regulatory_edge_truth"]
    assert result==one_target_edges("TF1",seq,rows,motifs)


def test_invalid_seq_and_wrong_length_are_not_silently_scanned():
    with pytest.raises(ValueError,match="Non-ACGT"):
        encode("N"*PROMOTER_LENGTH)
    with pytest.raises(ValueError,match="Non-ACGT"):
        encode("A"*12)
    with pytest.raises(ValueError,match="invalid DNA"):
        dinucleotide_shuffle("ACNGT",np.random.default_rng(7))


def test_bh_q_is_not_interpreted_as_direct_regulatory_binding():
    q=bh_q([.01,.04,.04,.5])
    assert np.allclose(q,[.04,.053333333333,.053333333333,.5])
    with pytest.raises(ValueError,match="Invalid empirical"):
        bh_q([.1,np.nan])


def test_selected_exact_motif_names_and_probabilities(tmp_path):
    rows=[]
    for i in range(16):
        rows.append({"gene_symbol":f"TF{i}","gencode_gene_id":f"ENSG{i:011d}.1",
                     "role":"regulator","motif_ids":f"MA{i:04d}.1"})
    rows[8]["motif_ids"]="MA0008.1|MA0016.1"
    selected=pd.DataFrame(rows)
    motif_file=tmp_path/"motifs.meme"
    with motif_file.open("w") as f:
        f.write("MEME version 4\n\nALPHABET= ACGT\n\n")
        for i in range(17):
            name="TF8" if i==16 else f"TF{i}"
            f.write(f"MOTIF MA{i:04d}.1 {name}\n")
            f.write("letter-probability matrix: alength= 4 w= 5 nsites= 22 E= 0\n")
            for _ in range(5):f.write(".7 .1 .1 .1\n")
            f.write("\n")
    res=parse_meme_selected(motif_file,selected)
    assert len(res)==16 and len(res["TF8"])==2
    assert sum(map(len,res.values()))==17
    assert all(all(m["width"]==5 for m in group) for group in res.values())
    text=motif_file.read_text().replace("MOTIF MA0000.1 TF0","MOTIF MA0000.1 TF9")
    motif_file.write_text(text)
    with pytest.raises(ValueError,match="not original TF"):
        parse_meme_selected(motif_file,selected)



def test_dinucleotide_null_is_deterministic_across_python_hash_seeds():
    """Guard against randomized set-iteration order changing RNG consumption."""
    import os
    import subprocess
    import sys
    code=(
        "import numpy as np,hashlib;"
        "from scripts.biov5_scan_grch38_promoters import dinucleotide_shuffle;"
        "s=('ACGTTGCGAATCCCAGCTTG'*60)+'A';"
        "r=np.random.default_rng(20261127);"
        "out='|'.join(dinucleotide_shuffle(s,r) for _ in range(4));"
        "print(hashlib.sha256(out.encode()).hexdigest())"
    )
    results=[]
    for seed in ("0","1","77","random"):
        env={**os.environ,"PYTHONHASHSEED":seed}
        value=subprocess.check_output(
            [sys.executable,"-c",code],env=env,text=True).strip()
        results.append(value)
    assert len(set(results))==1
