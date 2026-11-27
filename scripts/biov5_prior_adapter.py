"""Load genuine BIO-V5 promoter binding-potential scores in GCEV evidence format.

This adapter strictly aligns source-observed motif scores to an existing locked
directed candidate universe. It does NOT fit GCEV or access ENCODE outcomes.
A score of zero is still observed evidence (mask=True), not missingness.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_NAME = "jaspar2026_grch38_gene_promoter_motif_dinucleotide_null"
EXPECTED_EDGES = 1008


def load_biov5_sequence_prior(
    root: Path,
    edges: list[tuple[str,str]],
) -> dict[str,tuple[np.ndarray,np.ndarray]]:
    root=Path(root)
    manifest_path=root/"biov5_motif_prior_manifest.json"
    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("tested_directed_no_self_edges")!=EXPECTED_EDGES or
        manifest.get("stage")!="BIOV5_IN_SILICO_SEQUENCE_BINDING_POTENTIAL_NO_BIOLOGICAL_GOLD" or
        manifest.get("ENCODE_perturbation_response_labels_accessed") is not False or
        manifest.get("gcev_model_fitted") is not False):
        raise ValueError("BIO-V5 source evidence manifest does not pass qualification")
    path=root/"k562_tf_to_promoter_motif_prior_1008.csv"
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if manifest.get("file_sha256",{}).get(path.name)!=h:
        raise ValueError("BIO-V5 source score file checksum mismatch")
    table=pd.read_csv(path,dtype={"tf":"string","target":"string"})
    if (len(table)!=EXPECTED_EDGES or
        not {"tf","target","prior_score","source_observed","empirical_p","bh_q"}.issubset(table) or
        table.duplicated(["tf","target"]).any() or
        (table.tf==table.target).any()):
        raise ValueError("Expected 1,008 unique non-self BIO-V5 directed sources")
    scores=pd.to_numeric(table.prior_score,errors="coerce").to_numpy(float)
    p=pd.to_numeric(table.empirical_p,errors="coerce").to_numpy(float)
    q=pd.to_numeric(table.bh_q,errors="coerce").to_numpy(float)
    observed=table.source_observed.astype(str).str.lower()
    if (not np.isfinite(scores).all() or not np.all((scores>=0)&(scores<=1)) or
        not np.isfinite(p).all() or not np.isfinite(q).all() or
        not np.all((p>=0)&(p<=1)) or not np.all((q>=0)&(q<=1)) or
        not observed.eq("true").all()):
        raise ValueError("Unobserved or invalid original promoter motif evidence")
    if len(edges)!=EXPECTED_EDGES or len(set(edges))!=len(edges):
        raise ValueError("Model candidate edge universe differs from pinned 1,008 edges")
    observed_edges=list(zip(table.tf.astype(str),table.target.astype(str)))
    if set(observed_edges)!=set(edges):
        raise ValueError("Model TF/target candidate identities do not match the frozen source")
    index={pair:i for i,pair in enumerate(observed_edges)}
    aligned=np.asarray([scores[index[e]] for e in edges],dtype=float)
    mask=np.ones(len(edges),dtype=bool)
    return {SOURCE_NAME:(aligned,mask)}
