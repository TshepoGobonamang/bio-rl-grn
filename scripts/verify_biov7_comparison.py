"""Independent saved-output checks for the BIO-V7 K562 comparison.

This verifier does not import the fitting driver or fit models. It recomputes
metrics and paired-cell uncertainty from saved error matrices and checks the
frozen graph decisions. Optional reference comparison establishes reproduction
in another environment using the SAME source dataset, not independent biology.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


TRANSFORMS = ("log1p_original", "log1p_cp10k")
VARIANTS = ("primary_mixed", "matched_v1", "matched_v2", "motif_only")
METHODS = ("expression_only", "fixed_sources", "global_source", "no_graph_state",
           "no_parent_state", "no_action_state", "gcev")
COMPARATORS = ("no_graph_state", "fixed_sources", "expression_only")
DNASE = "encode_ENCSR000EOT_K562_dnase_promoter_coverage"
TOLERANCE = 1e-10


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(actual, expected, label: str) -> float:
    actual, expected = np.asarray(actual, float), np.asarray(expected, float)
    require(actual.shape == expected.shape, label + ": numeric shape mismatch")
    require(bool(np.isfinite(actual).all() and np.isfinite(expected).all()),
            label + ": nonfinite value")
    difference = float(np.max(np.abs(actual - expected))) if actual.size else 0.
    require(difference <= TOLERANCE, f"{label}: maximum error {difference} exceeds {TOLERANCE}")
    return difference


def hash_snapshot(root: Path, files: dict[str, str], label: str) -> int:
    require(bool(files), label + ": empty hash snapshot")
    for name, digest in files.items():
        rel = Path(name)
        require(not rel.is_absolute() and ".." not in rel.parts,
                label + ": unsafe relative path")
        require((root / rel).is_file(), label + ": missing " + name)
        require(sha(root / rel) == digest, label + ": checksum mismatch " + name)
    return len(files)


def keys(frame: pd.DataFrame, columns: list[str]) -> set[tuple]:
    require(not frame.duplicated(columns).any(), "Duplicate key: " + ",".join(columns))
    return set(frame[columns].itertuples(index=False, name=None))


def compare_table(actual: pd.DataFrame, reference: pd.DataFrame,
                  index: list[str], label: str) -> float:
    require(set(actual.columns) == set(reference.columns), label + ": column mismatch")
    require(keys(actual, index) == keys(reference, index), label + ": key mismatch")
    actual = actual.sort_values(index).reset_index(drop=True)
    reference = reference.sort_values(index).reset_index(drop=True)
    numeric = list(actual.select_dtypes(include=np.number).columns)
    error = close(actual[numeric], reference[numeric], label)
    for column in actual.columns.difference(numeric):
        require(actual[column].astype(str).equals(reference[column].astype(str)),
                label + ": categorical mismatch in " + column)
    return error


def verify(root: Path, reference: Path | None = None) -> dict:
    root = Path(root)
    evaluation = json.loads((root / "evaluation_manifest.json").read_text())
    frozen = json.loads((root / "freeze_manifest.json").read_text())
    audit = json.loads((root / "audit_summary.json").read_text())
    axes = json.loads((root / "per_cell_errors_axes.json").read_text())
    require(evaluation.get("stage") == "BIOV7_INTERNAL_SCORE_EVALUATED_TEST_SEALED",
            "Unexpected evaluation stage")
    require(frozen.get("stage") == "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE",
            "Unexpected freeze stage")
    require(evaluation.get("freeze_manifest_sha256") == sha(root / "freeze_manifest.json"),
            "Evaluation does not reference the current frozen snapshot")
    require(audit.get("freeze_manifest_sha256") == evaluation["freeze_manifest_sha256"],
            "Audit and evaluation freeze fingerprints disagree")
    for field in ("input_sha256", "protocol_sha256", "implementation_sha256"):
        require(evaluation.get(field) == frozen.get(field), "Freeze/evaluation mismatch: " + field)
    for record in (evaluation, frozen, audit):
        require(record.get("test_used") is False, "Test expression used")
        require(record.get("perturbation_outcomes_accessed") is False,
                "Perturbation response accessed")
    require(frozen.get("internal_score_cells_used") == 0 and frozen.get("fit_cells_used") == 394,
            "Freeze stage did not use only 394 fitting cells")
    require(audit.get("internal_score_cells_evaluated") == 131 and audit.get("test_cells_used") == 0,
            "Evaluation cell boundary mismatch")
    require(frozen.get("partition_counts") == {"fit": 394, "internal_score": 131, "untouched_test": 133},
            "Frozen partitions changed")
    checked_hashes = hash_snapshot(root, evaluation["output_sha256"], "evaluation")
    checked_hashes += hash_snapshot(root, frozen["output_sha256"], "freeze")

    summary = pd.read_csv(root / "summary.csv")
    contrasts = pd.read_csv(root / "contrasts.csv")
    per_target = pd.read_csv(root / "per_target_results.csv")
    source_audit = pd.read_csv(root / "source_action_audit.csv")
    expected_methods = {(t, v, m) for t in TRANSFORMS for v in VARIANTS for m in METHODS}
    require(len(summary) == 56 and keys(summary, ["transform", "variant", "method"]) == expected_methods,
            "Expected exactly 56 unique method summaries")
    expected_contrasts = {(t, v, "gcev_minus_" + c) for t in TRANSFORMS for v in VARIANTS for c in COMPARATORS}
    require(len(contrasts) == 24 and keys(contrasts, ["transform", "variant", "contrast"]) == expected_contrasts,
            "Expected exactly 24 unique paired contrasts")
    genes, cells = axes["columns"], axes["rows"]
    require(len(genes) == len(set(genes)) == 64 and len(cells) == len(set(cells)) == 131,
            "Cell-error axes must be 131 unique cells by 64 unique targets")
    expected_targets = {(t, v, m, gene) for t, v, m in expected_methods for gene in genes}
    require(len(per_target) == 3584 and keys(per_target, ["transform", "variant", "method", "target"]) == expected_targets,
            "Expected exactly 3,584 unique target metrics")
    require(per_target.total_parent_count.eq(3).all() and summary.total_graph_edges.eq(192).all(),
            "Graph budgets differ")
    require(summary.test_used.eq(False).all(), "Test use claimed in summary")
    summary_index = summary.set_index(["transform", "variant", "method"]).sort_index()
    contrast_index = contrasts.set_index(["transform", "variant", "contrast"]).sort_index()
    target_index = per_target.set_index(["transform", "variant", "method", "target"]).sort_index()
    expected_arrays = {"__".join(key) for key in expected_methods}
    maximum_error = 0.
    selection_tables = {}
    dnase_rows = 0
    with np.load(root / "per_cell_errors.npz", allow_pickle=False) as archive:
        require(set(archive.files) == expected_arrays, "Cell-error archive has missing/extra methods")
        draws = np.random.default_rng(20261127).integers(0, 131, size=(10000, 131))
        for transform in TRANSFORMS:
            for variant in VARIANTS:
                case = root / transform / variant
                selected = pd.read_csv(case / "frozen_selections.csv")
                require(len(selected) == 448 and keys(selected, ["method", "target"]) ==
                        {(method, gene) for method in METHODS for gene in genes},
                        f"{transform}/{variant}: frozen selection coverage mismatch")
                require(selected.groupby("target").graph_parents.nunique().eq(1).all(),
                        "Methods do not share target-specific scaffolds")
                for row in selected.itertuples(index=False):
                    parents = str(row.graph_parents).split("|")
                    require(len(parents) == len(set(parents)) == 2 and row.selected_tf not in parents,
                            "Selection does not insert exactly one TF into two-parent scaffold")
                    require(set(parents + [row.selected_tf]).issubset(genes[:16]) and
                            row.target not in parents + [row.selected_tf],
                            "Selection contains a target-only parent or self-edge")
                selection_tables[(transform, variant)] = selected
                selected_index = selected.set_index(["method", "target"]).sort_index()
                labels = pd.read_csv(case / "source_utility_labels.csv")
                expected_label_count = 192 if variant == "motif_only" else 384
                require(len(labels) == expected_label_count,
                        "Unexpected nested source-label count")
                require(not labels.duplicated(["fold", "target", "evidence_source"]).any(),
                        "Duplicate nested source-label row")
                require(labels.equal_edit_budget.eq(True).all() and
                        labels.used_evaluation_reference.eq(False).all() and labels.tf_only_scaffold.eq(True).all(),
                        "Nested source-valuation contract failed")
                dnase = labels.loc[labels.evidence_source == DNASE]
                require(len(dnase) == (0 if variant == "motif_only" else 192),
                        "Missing/extra DNase source labels")
                if len(dnase):
                    require(dnase.source_value.eq(0).all() and dnase.different_edge.eq(0).all() and
                            dnase.source_constant_within_target.eq(True).all() and
                            (dnase.minus_edge == dnase.plus_edge).all(),
                            "Target-constant DNase incorrectly obtained action utility")
                    dnase_rows += len(dnase)
                expression_selection = selected_index.loc["expression_only"].selected_tf.reindex(genes)
                for method in METHODS:
                    matrix = archive[f"{transform}__{variant}__{method}"]
                    require(matrix.shape == (131, 64), "Wrong cell-error dimensions")
                    require(np.isfinite(matrix).all() and (matrix >= 0).all(), "Invalid squared errors")
                    row = summary_index.loc[(transform, variant, method)]
                    target_means = matrix.mean(axis=0)
                    maximum_error = max(maximum_error, close(
                        [matrix.mean(), np.median(target_means), matrix[:, 16:].mean()],
                        [row.macro_nmse_64, row.median_target_nmse_64, row.macro_nmse_48_target_only],
                        "Saved summary versus cell errors"))
                    target_rows = target_index.loc[(transform, variant, method)].reindex(genes)
                    maximum_error = max(maximum_error, close(target_means, target_rows.nmse.to_numpy(),
                                                             "Per-target NMSE versus cell errors"))
                    current_selection = selected_index.loc[method].selected_tf.reindex(genes)
                    require(int(row.changed_targets_vs_expression) == int((current_selection != expression_selection).sum()),
                            "Reported changed-target count is inconsistent")
                    require(current_selection.tolist() == target_rows.selected_tf.tolist(),
                            "Evaluated choices differ from frozen selections")
                    require(selected_index.loc[method].graph_parents.reindex(genes).tolist() == target_rows.graph_parents.tolist(),
                            "Evaluated scaffolds differ from frozen selections")
                for comparator in COMPARATORS:
                    difference = (archive[f"{transform}__{variant}__gcev"] -
                                  archive[f"{transform}__{variant}__{comparator}"]).mean(axis=1)
                    low, high = np.quantile(difference[draws].mean(axis=1), [.025, .975])
                    row = contrast_index.loc[(transform, variant, "gcev_minus_" + comparator)]
                    maximum_error = max(maximum_error, close(
                        [difference.mean(), low, high],
                        [row.delta_macro_nmse, row.cell_bootstrap_ci95_low, row.cell_bootstrap_ci95_high],
                        "Paired contrast and cell-bootstrap interval"))
    for transform in TRANSFORMS:
        fixed_primary = selection_tables[(transform, "primary_mixed")].query("method == 'fixed_sources'").set_index("target").selected_tf.reindex(genes)
        for variant in VARIANTS:
            fixed = selection_tables[(transform, variant)].query("method == 'fixed_sources'").set_index("target").selected_tf.reindex(genes)
            require(fixed.equals(fixed_primary), "DNase offsets changed a fixed-weight within-target choice")
    require(len(source_audit) == 14, "Expected 14 source-by-case audit rows")
    dnase_audit = source_audit.loc[source_audit.evidence_source == DNASE]
    require(len(dnase_audit) == 6 and dnase_audit.different_edge_rows.eq(0).all() and
            dnase_audit.max_abs_source_value.eq(0).all(), "Source-audit DNase cancellation mismatch")

    reference_comparison = None
    if reference is not None:
        reference = Path(reference)
        ref_frozen = json.loads((reference / "freeze_manifest.json").read_text())
        require(frozen["input_sha256"] == ref_frozen["input_sha256"],
                "Reference reproduction did not use identical source bytes")
        require(frozen["protocol_sha256"] == ref_frozen["protocol_sha256"],
                "Reference reproduction did not use identical locked protocol")
        ref_summary = pd.read_csv(reference / "summary.csv")
        ref_contrasts = pd.read_csv(reference / "contrasts.csv")
        reference_error = max(
            compare_table(summary, ref_summary, ["transform", "variant", "method"], "Reference summary"),
            compare_table(contrasts, ref_contrasts, ["transform", "variant", "contrast"], "Reference contrasts"),
        )
        decision_columns = ["method", "target", "graph_parents", "selected_tf"]
        for (transform, variant), selected in selection_tables.items():
            ref_selected = pd.read_csv(reference / transform / variant / "frozen_selections.csv")
            compare_table(selected[decision_columns], ref_selected[decision_columns],
                          ["method", "target"], f"Reference decisions {transform}/{variant}")
        reference_comparison = {"source_hashes_identical": True, "protocol_hash_identical": True,
                                "n_identical_target_decisions": 3584,
                                "maximum_summary_and_contrast_error": reference_error,
                                "same_dataset_reproduction_not_independent_biological_validation": True}
    result = {
        "status": "PASSED", "checked_hash_entries": checked_hashes,
        "summary_rows": 56, "contrast_rows": 24, "per_target_rows": 3584,
        "cell_error_arrays": 56, "target_constant_dnase_rows_zero_utility": dnase_rows,
        "fixed_source_choices_identical_across_variants": True,
        "metric_and_bootstrap_maximum_error": maximum_error,
        "absolute_numeric_tolerance": TOLERANCE, "test_used": False,
        "regulatory_truth_evaluated": False, "reference_comparison": reference_comparison,
    }
    (root / "independent_verification.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.results, args.reference), indent=2))


if __name__ == "__main__":
    main()
