"""BIO-V6: load motif and promoter-accessibility sources into the GCEV interface.

DNase accessibility is *target-specific*, broadcast across regulators. The
motif*accessibility product is derived and intentionally NOT exported as a third
independent source. No model fitting or biological scoring occurs here.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

MOTIF_NAME="jaspar2026_GRCh38_promoter_motif_potential"
DNASE_NAME="encode_ENCSR000EOT_K562_dnase_promoter_coverage"
EXPECTED_MOTIF_SHA="ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0"


def load_two_biov6_evidence_sources(
    root:Path, edges:list[tuple[str,str]]
) -> dict[str,tuple[np.ndarray,np.ndarray]]:
    root=Path(root)
    manifest=json.loads((root/"biov6_source_audit.json").read_text())
    if (manifest.get("stage")!="BIOV6_REPLICATE_DNASE_ACCESSIBILITY_CONTEXT_NO_TF_BINDING_GOLD"
        or manifest.get("n_promoters")!=64 or manifest.get("n_directed_edges")!=1008
        or manifest.get("original_biov5_motif_prior_sha256")!=EXPECTED_MOTIF_SHA
        or manifest.get("ENCODE_perturbation_expression_outcomes_accessed") is not False
        or manifest.get("GCEV_predictions_fitted") is not False
        or manifest.get("joint_prior_is_independent_third_source") is not False
        or len(manifest.get("source_bed_file_summaries",[]))!=2):
        raise ValueError("BIO-V6 source provenance or study stage not qualified")
    path=root/"k562_tf_target_accessibility_and_motif_context_1008.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest["file_sha256"][path.name]:
        raise ValueError("Original K562 DNase evidence source checksum mismatch")
    data=pd.read_csv(path,dtype={"tf":"string","target":"string"})
    cols={"tf","target","motif_prior_score","dnase_accessibility_mean",
          "accessibility_source_observed","joint_motif_dnase_product"}
    if len(data)!=1008 or not cols.issubset(data) or data.duplicated(["tf","target"]).any():
        raise ValueError("BIO-V6 exact 1008-edge candidate source table invalid")
    original=list(zip(data.tf.astype(str),data.target.astype(str)))
    if (len(edges)!=1008 or len(set(edges))!=1008 or set(original)!=set(edges)
        or any(a==b for a,b in original)):
        raise ValueError("Mismatch with the original directed source candidate universe")
    if not data.accessibility_source_observed.astype(str).str.lower().eq("true").all():
        raise ValueError("An observed zero must not be recoded as missing accessibility")
    motif=pd.to_numeric(data.motif_prior_score,errors="coerce").to_numpy(dtype=float)
    dnase=pd.to_numeric(data.dnase_accessibility_mean,errors="coerce").to_numpy(dtype=float)
    joint=pd.to_numeric(data.joint_motif_dnase_product,errors="coerce").to_numpy(dtype=float)
    if not all(np.isfinite(v).all() and ((v>=0)&(v<=1)).all() for v in (motif,dnase,joint)):
        raise ValueError("Nonfinite or out-of-range original DNase/motif scores")
    if not np.allclose(joint,motif*dnase,atol=1e-12):
        raise ValueError("Source joint context inconsistent with its two parent sources")
    if data.groupby("target").dnase_accessibility_mean.nunique().max()!=1:
        raise ValueError("DNase promoter accessibility unexpectedly differs by TF")
    lookup={e:i for i,e in enumerate(original)}
    mask=np.ones(1008,dtype=bool)
    index=np.asarray([lookup[e] for e in edges],dtype=int)
    return {
        MOTIF_NAME:(motif[index].copy(),mask.copy()),
        DNASE_NAME:(dnase[index].copy(),mask.copy())
    }
