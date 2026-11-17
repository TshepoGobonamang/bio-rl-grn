"""BEELINE hESC ChIP reference is audit-only, never a prior/training input."""
from zipfile import ZipFile
import pandas as pd
from scripts import audit_m2d2_chip_overlap as mod


def test_posthoc_overlap_counts_and_independence_guard(tmp_path, monkeypatch):
    data = tmp_path / "BEELINE-Networks.zip"
    with ZipFile(data, "w") as zf:
        zf.writestr(
            "BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv",
            "Gene1,Gene2\nSOX2,NANOG\nPOU5F1,NANOG\nGATA6,NANOG\n"
        )
    # This fixture doesn't have the publisher's archive hash.
    monkeypatch.setattr(mod, "digest", lambda p: mod.NETWORKS_MD5)
    prior = tmp_path / "prior.csv"
    pd.DataFrame({
        "source": ["SOX2", "SOX2", "POU5F1"],
        "target": ["NANOG", "GATA6", "NANOG"],
        "observed_curated_association": [1, 1, 1],
    }).to_csv(prior,index=False)
    report = mod.audit(prior,data,tmp_path / "result")
    assert report["overlapping_directed_edges"] == 2
    assert report["prior_edges"] == 3
    assert report["reference_total_edges"] == 3
    assert report["source_level_provenance_disjointness_proved"] is False
    assert report["gold_used_for_feature_selection"] is False
