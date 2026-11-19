"""M2D13: fail-closed experimental-design readiness from observed D12 summary.

The evidence manifest is a human-audited record of source-backed design facts,
NOT an inferred metadata mapping. Unknown or absent evidence blocks inference.
"""
import argparse
import json
from pathlib import Path


def decide(d12, evidence):
    required = {"n_cells", "counts_by_field", "tf_batch_temp_cells", "batch_count_per_tf", "source_sha256"}
    if not required.issubset(d12):
        raise ValueError("D12 actual-source observation missing fields: " + str(sorted(required - set(d12))))
    if not d12["tf_batch_temp_cells"] or int(d12["n_cells"]) <= 0:
        raise ValueError("D12 observations empty")
    if sum(int(v["cells"]) for v in d12["tf_batch_temp_cells"]) != int(d12["n_cells"]):
        raise ValueError("D12 joint counts do not reconcile to cell total")
    for field in ("TF", "batch", "temp"):
        if sum(int(v["cells"]) for v in d12["counts_by_field"][field]) != int(d12["n_cells"]):
            raise ValueError(f"D12 {field} counts do not reconcile")
    observed_tf = {str(v["value"]) for v in d12["counts_by_field"]["TF"]}
    observed_batches = {str(v["value"]) for v in d12["counts_by_field"]["batch"]}
    controls = evidence.get("validated_control_tf_labels", [])
    if not isinstance(controls, list) or any(x not in observed_tf for x in controls):
        raise ValueError("Any confirmed control label must appear in actual D12 TF values")
    design = evidence.get("biological_culture_replicate_provenance", {})
    controlled = bool(controls) and bool(evidence.get("control_identity_source")) and evidence.get("control_identity_verified") is True
    replicate_ok = (design.get("verified") is True and
                    bool(design.get("source")) and
                    bool(design.get("independent_culture_id_field")) and
                    design.get("sequencing_batch_not_used_as_replicate") is True)
    mapping_ok = evidence.get("tf_barcode_mapping_verified") is True and bool(evidence.get("tf_barcode_mapping_source"))
    same_context = evidence.get("matched_cell_context_verified") is True and bool(evidence.get("matched_cell_context_source"))
    n_tf_multibatch = sum(int(v["n_batches"]) >= 2 for v in d12["batch_count_per_tf"])
    ready = bool(controlled and replicate_ok and mapping_ok and same_context)
    return {
        "source_sha256": d12["source_sha256"],
        "n_cells_reconciled": int(d12["n_cells"]),
        "n_observed_tfs": len(observed_tf),
        "n_observed_batches": len(observed_batches),
        "n_tfs_present_in_multiple_batches_not_replicates": n_tf_multibatch,
        "verified_control_labels": controls if controlled else [],
        "control_evidence_pass": controlled,
        "independent_culture_replication_pass": replicate_ok,
        "tf_map_independence_pass": mapping_ok,
        "matched_context_pass": same_context,
        "confirmatory_comparison_ready": ready,
        "decision": "ELIGIBLE_FOR_FURTHER_DESIGN_REVIEW" if ready else "NO_GO",
        "note": "Even PASS does not establish direct TF-target edges or independent reference-source disjointness."
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--d12", required=True, type=Path)
    p.add_argument("--evidence", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    result = decide(json.loads(a.d12.read_text()), json.loads(a.evidence.read_text()))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
