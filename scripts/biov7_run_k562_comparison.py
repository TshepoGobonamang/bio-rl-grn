"""BIO-V7 frozen-cohort, matched one-edit K562 comparison.

The two CLI stages enforce a concrete boundary: ``freeze`` consumes only the
394 fitting cells and writes all decisions; ``evaluate`` verifies that complete
snapshot before computing the 131-cell internal-score endpoint. The 133 test
cells are never transformed, summarized, used to fit, or predicted. Loading and
hashing source archive files is an integrity operation, not test evaluation.

No reference regulatory network or perturbation response enters this driver.
This is a developmental predictive experiment, not biological confirmation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import sparse
import sklearn

from grnbench.core import candidate_pairs
from grnbench.gcev import compose_graph_conditional_rank, learn_graph_conditional_weights
from scripts.biov6_two_evidence_adapter import (
    DNASE_NAME, MOTIF_NAME, load_two_biov6_evidence_sources,
)
from scripts.biov7_constrained_models import (
    evaluate_graphs, final_contexts, nested_valuation,
)


TRANSFORMS = ["log1p_original", "log1p_cp10k"]
VARIANTS = ["primary_mixed", "matched_v1", "matched_v2", "motif_only"]
METHODS = ["expression_only", "fixed_sources", "global_source", "no_graph_state",
           "no_parent_state", "no_action_state", "gcev"]
CONTRASTS = ["no_graph_state", "fixed_sources", "expression_only"]
COUNTS = {"fit": 394, "internal_score": 131, "untouched_test": 133}
SOURCE_H5AD_SHA = "14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6"
MOTIF_SHA = "ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0"
ORIGINAL_MATRIX_SHA = "0e5dd2f2be0e8924ac089c5f0e620ce94950902ecfd65fa596faafa35aa6e897"
PRIOR_TABLE_HASHES = {
    "biov6/primary/k562_tf_target_accessibility_and_motif_context_1008.csv":
        "1d6ab9b17a0ee8d535481451a8688eedb74916675e0c29e247f0bdf06784e5ca",
    "biov6/pipeline_sensitivity/k562_dnase_v1_matched_joint_context_1008.csv":
        "fa8d85d837cab530e5fd8214c22b5f94e2c6b7bffbf32c2eb13e067f88d5c9cc",
    "biov6/pipeline_sensitivity/k562_dnase_v2_matched_joint_context_1008.csv":
        "0a7871079188ed148deed0b2d5aaa182f9ed8bb7ee9b6154fedd94633e78d951",
}
COHORT_HASHES = {
    "k562_tfaware64_original_scale.npz": "5271925dcb486923977383334d39a513e5aa77474de5d16f3c506a7602742d20",
    "selected_16_TFs.csv": "9bff00f39f2f07a1b7a275472c1fe5d551e491879d2d7d86f85f2a21202effec",
    "selected_48_targets.csv": "a3fd75cc760112c5d4664a94aca1d495e66328e174c033a4793e16868ecf7256",
    "selected_gene_universe.csv": "b8dd6122a56cc82b00ed6dfd9754209a61905631b0f122f5394dde7ccb5e16c9",
    "cell_partitions.csv": "c1b562feeb150d299f2c7c331966084ed65471759f114e1963b782f6d68bed73",
}
IMPLEMENTATION_FILES = [
    "scripts/biov7_run_k562_comparison.py", "scripts/biov7_constrained_models.py",
    "scripts/biov6_two_evidence_adapter.py", "src/grnbench/gcev.py",
    "src/grnbench/core.py",
]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def environment() -> dict:
    return {"python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__, "platform": platform.platform()}


def implementation_hashes() -> dict[str, str]:
    root = Path(__file__).resolve().parents[1]
    return {name: sha(root / name) for name in IMPLEMENTATION_FILES}


def load_protocol(path: Path) -> dict:
    protocol = read_json(path)
    expected = {"seed": 20261127, "alpha": 1.0, "max_parents": 2, "boost": .3,
                "n_folds": 3, "shrinkage": .05, "bootstrap_resamples": 10000,
                "transforms": TRANSFORMS, "variants": VARIANTS, "methods": METHODS}
    if not isinstance(protocol.get("protocol"), str) or not protocol["protocol"]:
        raise ValueError("Locked protocol requires an identifier")
    for name, value in expected.items():
        if protocol.get(name) != value:
            raise ValueError(f"Unexpected locked BIO-V7 protocol setting: {name}")
    return protocol


@dataclass
class Inputs:
    matrix: sparse.csr_matrix
    original_matrix: sparse.csr_matrix
    original_indices: list[int]
    partitions: pd.DataFrame
    genes: list[str]
    regulators: list[str]
    target_only: list[str]
    edges: list[tuple[str, str]]
    evidence: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]
    hashes: dict[str, str]

    def frame(self, partition: str, transform: str) -> pd.DataFrame:
        """Slice before any numeric transformation; test access is forbidden."""
        if partition not in ("fit", "internal_score"):
            raise ValueError("BIO-V7 prohibits test expression access")
        if transform not in TRANSFORMS:
            raise ValueError("Unknown predeclared transform")
        rows = np.flatnonzero(self.partitions.partition.to_numpy() == partition)
        selected = self.matrix[rows, :].toarray().astype(float)
        original = self.original_matrix[rows, :][:, self.original_indices].toarray()
        if not np.array_equal(selected, original):
            raise ValueError(f"Selected genes do not match original source on {partition}")
        if not np.isfinite(selected).all() or (selected < 0).any():
            raise ValueError(f"Invalid original-scale values on {partition}")
        if transform == "log1p_cp10k":
            totals = np.asarray(self.original_matrix[rows, :].sum(axis=1)).ravel()
            if not np.isfinite(totals).all() or (totals <= 0).any():
                raise ValueError("Invalid full-transcriptome row totals")
            selected = selected / totals[:, None] * 10000.
        values = np.log1p(selected)
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite transformed expression")
        return pd.DataFrame(values, columns=self.genes,
                            index=self.partitions.cell_id.iloc[rows].tolist())


def load_inputs(cohort: Path, original: Path, biov5: Path, biov6: Path) -> Inputs:
    hashes: dict[str, str] = {}

    def record(prefix: str, root: Path, name: str, expected: str | None = None) -> Path:
        path = root / name
        digest = sha(path)
        if expected is not None and digest != expected:
            raise ValueError(f"Input checksum mismatch: {prefix}/{name}")
        hashes[f"{prefix}/{name}"] = digest
        return path

    cohort_manifest = read_json(record("cohort", cohort, "biov4b_fit_only_manifest.json"))
    if cohort_manifest.get("original_K562_H5AD_sha256") != SOURCE_H5AD_SHA:
        raise ValueError("BIO-V4B source identity mismatch")
    if cohort_manifest.get("sha256") != COHORT_HASHES:
        raise ValueError("BIO-V4B frozen file manifest mismatch")
    for name, digest in COHORT_HASHES.items():
        record("cohort", cohort, name, digest)
    partitions = pd.read_csv(cohort / "cell_partitions.csv", dtype=str)
    universe = pd.read_csv(cohort / "selected_gene_universe.csv", dtype=str)
    genes = universe.gene_symbol.tolist()
    regulators = pd.read_csv(cohort / "selected_16_TFs.csv", dtype=str).gene_symbol.tolist()
    target_only = pd.read_csv(cohort / "selected_48_targets.csv", dtype=str).gene_symbol.tolist()
    if (len(partitions) != 658 or not partitions.cell_id.is_unique or
            partitions.partition.value_counts().to_dict() != COUNTS):
        raise ValueError("Original 394/131/133 partition identities invalid")
    if (len(genes) != 64 or len(set(genes)) != 64 or len(regulators) != 16 or
            len(target_only) != 48 or genes != regulators + target_only or
            universe.role.tolist() != ["regulator"] * 16 + ["target_only"] * 48):
        raise ValueError("Frozen TF/target universe mismatch")

    source_manifest = read_json(record("original", original, "k562_source_audit.json"))
    if (source_manifest.get("source_sha256") != SOURCE_H5AD_SHA or
            source_manifest.get("n_cells") != 658 or source_manifest.get("n_genes") != 60671):
        raise ValueError("Original K562 source provenance mismatch")
    if source_manifest["export_hashes"].get("k562_expression_original_scale.npz") != ORIGINAL_MATRIX_SHA:
        raise ValueError("Original full-transcriptome matrix fingerprint mismatch")
    for name, digest in source_manifest["export_hashes"].items():
        record("original", original, name, digest)
    source_cells = pd.read_csv(original / "k562_cell_ids.csv", dtype=str).cell_id.tolist()
    source_genes = pd.read_csv(original / "k562_gene_ids.csv", dtype=str).gene_id.tolist()
    if source_cells != partitions.cell_id.tolist() or len(source_genes) != 60671 or len(set(source_genes)) != 60671:
        raise ValueError("Source gene/cell identities mismatch")
    lookup = {gene: index for index, gene in enumerate(source_genes)}
    if not set(genes).issubset(lookup):
        raise ValueError("Frozen gene missing from original source")
    source_indices = [lookup[gene] for gene in genes]
    matrix = sparse.load_npz(cohort / "k562_tfaware64_original_scale.npz").tocsr()
    original_matrix = sparse.load_npz(original / "k562_expression_original_scale.npz").tocsr()
    if matrix.shape != (658, 64) or original_matrix.shape != (658, 60671):
        raise ValueError("Frozen expression dimensions mismatch")

    motif_manifest = read_json(record("biov5", biov5, "biov5_motif_prior_manifest.json"))
    if (motif_manifest.get("source_H5AD_sha256") != SOURCE_H5AD_SHA or
            motif_manifest.get("original_biov4b_manifest_sha256") != hashes["cohort/biov4b_fit_only_manifest.json"] or
            motif_manifest["file_sha256"].get("k562_tf_to_promoter_motif_prior_1008.csv") != MOTIF_SHA):
        raise ValueError("BIO-V5 source chain mismatch")
    for name, digest in motif_manifest["file_sha256"].items():
        record("biov5", biov5, name, digest)
    motif_table = pd.read_csv(biov5 / "k562_tf_to_promoter_motif_prior_1008.csv")
    edges = candidate_pairs(genes, regulators)
    if len(edges) != 1008 or motif_table.duplicated(["tf", "target"]).any():
        raise ValueError("Candidate edge set mismatch")
    motif_table = motif_table.set_index(["tf", "target"])
    if set(motif_table.index) != set(edges):
        raise ValueError("Motif candidate universe mismatch")
    motif = motif_table.loc[edges, "prior_score"].to_numpy(float)
    primary_root = biov6 / "primary"
    primary_manifest = read_json(record("biov6/primary", primary_root, "biov6_source_audit.json"))
    for name, digest in primary_manifest["file_sha256"].items():
        record("biov6/primary", primary_root, name, digest)
    primary = load_two_biov6_evidence_sources(primary_root, edges)
    if not np.allclose(primary[MOTIF_NAME][0], motif, rtol=0., atol=1e-12):
        raise ValueError("BIO-V6 motif values differ from frozen BIO-V5")
    evidence = {"primary_mixed": primary}
    sensitivity_root = biov6 / "pipeline_sensitivity"
    sensitivity_manifest = read_json(record("biov6/pipeline_sensitivity", sensitivity_root,
                                            "biov6_secondary_pipeline_sensitivity_manifest.json"))
    if (sensitivity_manifest.get("parent_biov5_prior_sha256") != MOTIF_SHA or
            sensitivity_manifest.get("original_source_experiment") != "ENCSR000EOT" or
            sensitivity_manifest.get("ENCODE_CRISPRi_response_outcomes_accessed") is not False):
        raise ValueError("BIO-V6 sensitivity provenance mismatch")
    for name, digest in sensitivity_manifest["file_sha256"].items():
        record("biov6/pipeline_sensitivity", sensitivity_root, name, digest)
    for version in ("v1", "v2"):
        table = pd.read_csv(sensitivity_root / f"k562_dnase_{version}_matched_joint_context_1008.csv")
        if len(table) != 1008 or table.duplicated(["tf", "target"]).any():
            raise ValueError("Matched pipeline candidate table malformed")
        table = table.set_index(["tf", "target"])
        if set(table.index) != set(edges):
            raise ValueError("Matched pipeline candidate universe differs")
        table = table.loc[edges]
        values = table.dnase_accessibility_mean.to_numpy(float)
        if (not np.isfinite(values).all() or (values < 0).any() or (values > 1).any() or
                not table.accessibility_source_observed.astype(str).str.lower().eq("true").all() or
                table.groupby(level="target").dnase_accessibility_mean.nunique().max() != 1 or
                not np.allclose(table.motif_prior_score.to_numpy(float), motif, rtol=0., atol=1e-12) or
                not np.allclose(table.joint_motif_dnase_product.to_numpy(float), motif * values, rtol=0., atol=1e-12)):
            raise ValueError("Invalid matched-pipeline motif/accessibility evidence")
        evidence[f"matched_{version}"] = {
            MOTIF_NAME: (motif.copy(), np.ones(1008, bool)),
            DNASE_NAME: (values, np.ones(1008, bool)),
        }
    evidence["motif_only"] = {MOTIF_NAME: (motif.copy(), np.ones(1008, bool))}
    for name, digest in PRIOR_TABLE_HASHES.items():
        if hashes.get(name) != digest:
            raise ValueError("Frozen BIO-V6 source version differs: " + name)
    return Inputs(matrix, original_matrix, source_indices, partitions, genes,
                  regulators, target_only, edges, evidence, hashes)


def invariant_action_evidence(evidence: dict, edges: list[tuple[str, str]]) -> dict:
    """Drop exactly constant target offsets for numerically stable argmax only.

    The unmodified scores are still saved. A common offset cannot change the
    target's chosen TF, regardless of a learned source's target-level weight.
    """
    result = {name: (values.copy(), mask.copy())
              for name, (values, mask) in evidence.items()}
    targets = np.asarray([target for _, target in edges])
    for values, mask in result.values():
        for target in np.unique(targets):
            ix = np.flatnonzero(targets == target)
            if np.ptp(values[ix] * mask[ix]) == 0.:
                values[ix] = 0.
    return result


def select_actions(contexts: pd.DataFrame, evidence: dict, weights: pd.DataFrame,
                   methods: list[str], boost: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidate = contexts.loc[contexts.eligible_action].reset_index(drop=True)
    edges = list(zip(candidate.tf, candidate.target))
    action_ev = invariant_action_evidence(evidence, edges)
    base = candidate.expr_rank.to_numpy(float)
    scores, choices = [], []
    for method in methods:
        if method == "expression_only":
            raw_scores = action_scores = base.copy()
        else:
            raw_scores = compose_graph_conditional_rank(base, edges, evidence, weights, method, boost)
            action_scores = compose_graph_conditional_rank(base, edges, action_ev, weights, method, boost)
        frame = candidate[["tf", "target", "graph_parents", "gain", "expr_rank"]].copy()
        frame["method"], frame["score"], frame["action_score"] = method, raw_scores, action_scores
        frame["constant_target_offsets_removed_for_argmax"] = True
        frame["global_score_is_not_source_utility_endpoint"] = True
        scores.append(frame)
        for target, target_frame in frame.groupby("target", sort=False):
            chosen = target_frame.sort_values(["action_score", "tf"], ascending=[False, True], kind="stable").iloc[0]
            choices.append({"method": method, "target": target,
                            "graph_parents": chosen.graph_parents,
                            "selected_tf": chosen.tf, "selected_score": float(chosen.score),
                            "selected_action_score": float(chosen.action_score),
                            "parent_count_before": len(chosen.graph_parents.split("|")),
                            "inserted_edge_count": 1, "total_parent_count": 3})
    return pd.concat(scores, ignore_index=True), pd.DataFrame(choices)


def snapshot_files(output: Path, exclude: set[str] | None = None) -> dict[str, str]:
    exclude = exclude or set()
    return {str(path.relative_to(output)): sha(path)
            for path in sorted(output.rglob("*"))
            if path.is_file() and str(path.relative_to(output)) not in exclude}


def freeze(args: argparse.Namespace) -> dict:
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Freeze output must be empty; refusing to overwrite a snapshot")
    protocol = load_protocol(args.protocol)
    data = load_inputs(args.cohort, args.original, args.biov5, args.biov6)
    args.output.mkdir(parents=True, exist_ok=True)
    initial_code_hashes = implementation_hashes()
    start = time.perf_counter()
    source_audits, timings = [], []
    for transform in protocol["transforms"]:
        fit = data.frame("fit", transform)
        transform_dir = args.output / transform
        transform_dir.mkdir()
        tick = time.perf_counter()
        contexts = final_contexts(fit, data.edges, alpha=protocol["alpha"],
                                  max_parents=protocol["max_parents"], seed=protocol["seed"],
                                  n_folds=protocol["n_folds"])
        if len(contexts) != 1008 or int(contexts.eligible_action.sum()) != 880:
            raise AssertionError("Expected 1,008 prior edges, 128 scaffold edges and 880 edits")
        contexts.to_csv(transform_dir / "candidate_contexts.csv", index=False)
        timings.append({"transform": transform, "variant": "shared", "stage": "final_contexts",
                        "seconds": time.perf_counter() - tick})
        eligible = contexts.eligible_action.to_numpy(bool)
        for variant in protocol["variants"]:
            case = transform_dir / variant
            case.mkdir()
            full_evidence = data.evidence[variant]
            tick = time.perf_counter()
            labels = nested_valuation(fit, data.edges, full_evidence, alpha=protocol["alpha"],
                                      max_parents=protocol["max_parents"], boost=protocol["boost"],
                                      seed=protocol["seed"], n_folds=protocol["n_folds"])
            expected_rows = 3 * 64 * len(full_evidence)
            if len(labels) != expected_rows:
                raise AssertionError("Unexpected source utility label count")
            learned = learn_graph_conditional_weights(labels, shrinkage=protocol["shrinkage"],
                                                      seed=protocol["seed"])
            labels.to_csv(case / "source_utility_labels.csv", index=False)
            learned.per_target_source.to_csv(case / "learned_weights.csv", index=False)
            learned.paired_decisions.to_csv(case / "paired_decisions_with_oof_predictions.csv", index=False)
            write_json(case / "gate_metadata.json", learned.gate_fit_metadata)
            timings.append({"transform": transform, "variant": variant,
                            "stage": "nested_valuation_and_gate", "seconds": time.perf_counter() - tick})
            tick = time.perf_counter()
            candidate_evidence = {name: (values[eligible], mask[eligible])
                                  for name, (values, mask) in full_evidence.items()}
            scores, selections = select_actions(contexts, candidate_evidence,
                                                 learned.per_target_source,
                                                 protocol["methods"], protocol["boost"])
            if len(scores) != 880 * len(METHODS) or len(selections) != 64 * len(METHODS):
                raise AssertionError("Unequal candidate or action budget")
            if not set(selections.selected_tf).issubset(data.regulators):
                raise AssertionError("Non-TF action escaped frozen candidate universe")
            scores.to_csv(case / "candidate_scores.csv", index=False)
            selections.to_csv(case / "frozen_selections.csv", index=False)
            timings.append({"transform": transform, "variant": variant,
                            "stage": "compose_and_select", "seconds": time.perf_counter() - tick})
            for source, group in labels.groupby("evidence_source", sort=True):
                learned_source = learned.per_target_source.loc[
                    learned.per_target_source.evidence_source == source]
                audit = {"transform": transform, "variant": variant, "evidence_source": source,
                         "n_nested_rows": len(group), "different_edge_rows": int(group.different_edge.sum()),
                         "constant_within_target_rows": int(group.source_constant_within_target.sum()),
                         "mean_source_value": float(group.source_value.mean()),
                         "max_abs_source_value": float(group.source_value.abs().max()),
                         "mean_gcev_weight": float(learned_source.gcev_weight.mean()),
                         "positive_gcev_weights": int((learned_source.gcev_weight > 0).sum()),
                         "source_is_target_constant_dnase": source == DNASE_NAME,
                         "positive_dnase_weights_are_gate_extrapolation_not_action_utility": source == DNASE_NAME}
                if source == DNASE_NAME and (audit["different_edge_rows"] != 0 or audit["max_abs_source_value"] != 0.):
                    raise AssertionError("Target-constant DNase generated false action utility")
                source_audits.append(audit)
            print(f"FROZEN {transform}/{variant}: {len(labels)} source rows, 448 matched selections", flush=True)
    pd.DataFrame(source_audits).to_csv(args.output / "source_action_audit.csv", index=False)
    pd.DataFrame(timings).to_csv(args.output / "training_timings.csv", index=False)
    if initial_code_hashes != implementation_hashes():
        raise RuntimeError("Implementation changed during fitting; freeze invalid")
    manifest = {
        "stage": "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE",
        "protocol": protocol["protocol"], "protocol_sha256": sha(args.protocol),
        "created_utc": now(), "environment": environment(),
        "implementation_sha256": initial_code_hashes, "input_sha256": data.hashes,
        "partition_counts": COUNTS, "fit_cells_used": 394, "internal_score_cells_used": 0,
        "test_cells_used": 0, "test_used": False, "perturbation_outcomes_accessed": False,
        "n_genes": 64, "n_TFs": 16, "candidate_prior_edges": 1008,
        "scaffold_edges": 128, "eligible_single_edits": 880,
        "total_graph_edges_per_method": 192, "inserted_edges_per_method": 64,
        "transforms": TRANSFORMS, "variants": VARIANTS, "methods": METHODS,
        "training_seconds_total": time.perf_counter() - start,
        "test_expression_analytics_performed": False,
        "source_original_scale_count_like_not_independently_certified_raw_UMI": True,
        "cp10k_is_predeclared_processing_sensitivity": True,
        "historical_synthetic_method_redefined": False,
        "model_adaptation": "common TF-only scaffolds, same candidate insertions, frozen cell partitions",
        "final_gain_is_selection_statistic_not_unbiased_scaffold_CV": True,
        "output_sha256": snapshot_files(args.output),
    }
    write_json(args.output / "freeze_manifest.json", manifest)
    return manifest


def paired_cell_interval(left: np.ndarray, right: np.ndarray, *, seed: int,
                         n_resamples: int) -> tuple[float, float, float]:
    if left.shape != right.shape or left.ndim != 2 or left.shape[0] != 131:
        raise ValueError("Paired error arrays must align on 131 internal-score cells")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("Nonfinite paired errors")
    difference = (left - right).mean(axis=1)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(difference), size=(n_resamples, len(difference)))
    bootstrap = difference[draws].mean(axis=1)
    low, high = np.quantile(bootstrap, [.025, .975])
    return float(difference.mean()), float(low), float(high)


def evaluate(args: argparse.Namespace) -> dict:
    if (args.output / "evaluation_manifest.json").exists() or (args.output / "summary.csv").exists():
        raise ValueError("Internal-score results already exist; refusing to overwrite")
    protocol = load_protocol(args.protocol)
    frozen = read_json(args.output / "freeze_manifest.json")
    freeze_sha = sha(args.output / "freeze_manifest.json")
    if (frozen.get("stage") != "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE" or
            frozen.get("protocol_sha256") != sha(args.protocol) or
            frozen.get("implementation_sha256") != implementation_hashes() or
            frozen.get("test_used") is not False or frozen.get("internal_score_cells_used") != 0):
        raise ValueError("Frozen protocol/implementation changed before evaluation")
    for name, digest in frozen["output_sha256"].items():
        if sha(args.output / name) != digest:
            raise ValueError("Frozen decisions modified before evaluation: " + name)
    data = load_inputs(args.cohort, args.original, args.biov5, args.biov6)
    if data.hashes != frozen["input_sha256"]:
        raise ValueError("Frozen inputs changed before evaluation")
    start = time.perf_counter()
    target_frames, summary_rows, contrast_rows, errors_archive = [], [], [], {}
    selected_by_case = {}
    target_only_indices = [data.genes.index(gene) for gene in data.target_only]
    for transform in protocol["transforms"]:
        fit, score = data.frame("fit", transform), data.frame("internal_score", transform)
        for variant in protocol["variants"]:
            case = args.output / transform / variant
            selections = pd.read_csv(case / "frozen_selections.csv")
            selected_by_case[(transform, variant)] = selections
            per_target, errors = evaluate_graphs(fit, score, selections, alpha=protocol["alpha"])
            if not per_target.total_parent_count.eq(3).all():
                raise AssertionError("Evaluation action budget mismatch")
            per_target.insert(0, "variant", variant)
            per_target.insert(0, "transform", transform)
            target_frames.append(per_target)
            expression_selection = selections.loc[selections.method == "expression_only"].set_index("target").selected_tf
            for method in protocol["methods"]:
                matrix = errors[method]
                if matrix.shape != (131, 64) or not np.isfinite(matrix).all():
                    raise AssertionError("Malformed cell-error matrix")
                errors_archive[f"{transform}__{variant}__{method}"] = matrix
                method_selection = selections.loc[selections.method == method].set_index("target").selected_tf
                changed = int((method_selection.reindex(data.genes) != expression_selection.reindex(data.genes)).sum())
                summary_rows.append({
                    "transform": transform, "variant": variant, "method": method,
                    "primary_case": transform == "log1p_original" and variant == "primary_mixed",
                    "macro_nmse_64": float(matrix.mean()),
                    "median_target_nmse_64": float(np.median(matrix.mean(axis=0))),
                    "macro_nmse_48_target_only": float(matrix[:, target_only_indices].mean()),
                    "changed_targets_vs_expression": changed,
                    "n_fit_cells": 394, "n_internal_score_cells": 131,
                    "n_targets": 64, "total_graph_edges": 192,
                    "test_used": False,
                })
            for comparator in CONTRASTS:
                delta, low, high = paired_cell_interval(errors["gcev"], errors[comparator],
                                                       seed=protocol["seed"],
                                                       n_resamples=protocol["bootstrap_resamples"])
                contrast_rows.append({
                    "transform": transform, "variant": variant,
                    "contrast": f"gcev_minus_{comparator}",
                    "primary_contrast": transform == "log1p_original" and variant == "primary_mixed" and comparator == "no_graph_state",
                    "delta_macro_nmse": delta, "cell_bootstrap_ci95_low": low,
                    "cell_bootstrap_ci95_high": high, "n_resamples": protocol["bootstrap_resamples"],
                    "n_resampled_cells": 131, "negative_favors_gcev": True,
                    "interval_conditional_on_this_cell_cohort_and_fitted_models": True,
                    "biological_replicate_confidence_interval": False,
                    "multiple_comparisons_adjusted": False,
                })
            print(f"EVALUATED {transform}/{variant}: 7 frozen methods on 131 internal-score cells", flush=True)
    pd.DataFrame(summary_rows).to_csv(args.output / "summary.csv", index=False)
    pd.DataFrame(contrast_rows).to_csv(args.output / "contrasts.csv", index=False)
    pd.concat(target_frames, ignore_index=True).to_csv(args.output / "per_target_results.csv", index=False)
    np.savez_compressed(args.output / "per_cell_errors.npz", **errors_archive)
    write_json(args.output / "per_cell_errors_axes.json", {
        "rows": data.partitions.loc[data.partitions.partition == "internal_score", "cell_id"].tolist(),
        "columns": data.genes, "array_key": "<transform>__<variant>__<method>",
        "value": "squared prediction error / max(fit target variance ddof0, 1e-9)",
    })
    cancellation_checks = []
    for transform in TRANSFORMS:
        primary_fixed = selected_by_case[(transform, "primary_mixed")].query("method == 'fixed_sources'").set_index("target").selected_tf
        for variant in VARIANTS:
            fixed = selected_by_case[(transform, variant)].query("method == 'fixed_sources'").set_index("target").selected_tf
            differences = int((fixed.reindex(data.genes) != primary_fixed.reindex(data.genes)).sum())
            cancellation_checks.append({"transform": transform, "variant": variant,
                                        "fixed_weight_changed_targets_vs_primary": differences})
            if differences:
                raise AssertionError("Constant DNase changed a fixed-weight within-target action")
    audit = {
        "partition_counts": COUNTS, "fit_cells_used": 394, "internal_score_cells_evaluated": 131,
        "test_cells_used": 0, "test_used": False,
        "test_transforms_or_predictions_or_descriptive_statistics_computed": False,
        "perturbation_outcomes_accessed": False, "reference_network_used": False,
        "same_candidate_universe_all_methods": True, "TF_only_parents_and_actions": True,
        "same_two_parent_scaffold_plus_one_edit": True,
        "frozen_before_internal_score": True, "freeze_manifest_sha256": freeze_sha,
        "matched_fixed_weight_cancellation_checks": cancellation_checks,
        "source_action_audit_file": "source_action_audit.csv",
        "uncertainty": "paired cell bootstrap conditional on one fitted source cohort; not independent biological replication",
        "source_weight_interpretation": "predicted action utility, not binding probability; positive DNase weights are extrapolation",
        "study_status": "DEVELOPMENT_INTERNAL_EXPRESSION_GENERALIZATION_NOT_BIOLOGICAL_CONFIRMATION",
        "primary_endpoint": "macro NMSE across 64 targets on original log1p X, primary mixed-caller prior",
        "primary_contrast": "GCEV minus no_graph_state; negative delta favors GCEV",
    }
    write_json(args.output / "audit_summary.json", audit)
    if implementation_hashes() != frozen["implementation_sha256"]:
        raise RuntimeError("Implementation changed during evaluation")
    evaluation = {
        "stage": "BIOV7_INTERNAL_SCORE_EVALUATED_TEST_SEALED",
        "protocol": protocol["protocol"], "protocol_sha256": sha(args.protocol),
        "freeze_manifest_sha256": freeze_sha, "created_utc": now(),
        "environment": environment(), "implementation_sha256": implementation_hashes(),
        "input_sha256": data.hashes, "evaluation_seconds_total": time.perf_counter() - start,
        "n_cases": 8, "n_methods_per_case": 7, "n_internal_score_cells": 131,
        "test_used": False, "perturbation_outcomes_accessed": False,
        "output_sha256": snapshot_files(args.output, {"evaluation_manifest.json"}),
    }
    write_json(args.output / "evaluation_manifest.json", evaluation)
    return evaluation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["freeze", "evaluate"])
    for name in ("cohort", "original", "biov5", "biov6", "protocol", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = freeze(args) if args.stage == "freeze" else evaluate(args)
    print(json.dumps({key: result[key] for key in ("stage", "protocol", "test_used")}, indent=2))


if __name__ == "__main__":
    main()
