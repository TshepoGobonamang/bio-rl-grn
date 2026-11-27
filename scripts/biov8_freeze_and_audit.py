"""BIO-V8: freeze existing action scores and audit source identifiability.

This driver has no expression, model-fitting, evaluation, or assay-download
interface. It reads a narrow allow-list of BIO-V7 fit-only output files plus
three already-qualified prior tables. Whole-directory scanning is deliberately
avoided because evaluation outputs live beside the permitted fit-only files.

Legacy ``score`` is preserved as text for exploratory across-target ranking;
``action_score`` is preserved separately for the original within-target choice.
Neither is a probability. Scaffold edges receive no invented insertion score.
New fixed-form score probes diagnose decisions, never predictive usefulness.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


TRANSFORMS = ["log1p_original", "log1p_cp10k"]
VARIANTS = ["primary_mixed", "matched_v1", "matched_v2", "motif_only"]
METHODS = ["expression_only", "fixed_sources", "global_source", "no_graph_state",
           "no_parent_state", "no_action_state", "gcev"]
REPRESENTATIONS = ["expression_only", "motif", "global_mean_accessibility_motif",
                   "matched_accessibility_motif", "permuted_accessibility_motif",
                   "legacy_additive_motif_accessibility", "accessibility_only"]
PRIOR_KEYS = {
    "primary_mixed": ("primary_biov6_table_sha256",
        "biov6/primary/k562_tf_target_accessibility_and_motif_context_1008.csv"),
    "matched_v1": ("matched_v1_table_sha256",
        "biov6/pipeline_sensitivity/k562_dnase_v1_matched_joint_context_1008.csv"),
    "matched_v2": ("matched_v2_table_sha256",
        "biov6/pipeline_sensitivity/k562_dnase_v2_matched_joint_context_1008.csv"),
}
DEFAULT_V7_PROTOCOL = Path(__file__).resolve().parents[1] / "configs/biov7_k562_comparison_2026_11_27.json"
N_GENES, N_TFS, N_EDGES, N_ACTIONS = 64, 16, 1008, 880


def sha(path: Path) -> str:
    """Hash only the explicit requested file; never discover nearby inputs."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def bool_values(values: pd.Series, name: str) -> np.ndarray:
    text = values.astype(str).str.lower()
    require(text.isin(["true", "false"]).all(), f"Invalid boolean values: {name}")
    return text.eq("true").to_numpy()


def numeric(values: pd.Series, name: str, *, finite: bool = True) -> np.ndarray:
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric values: {name}") from exc
    if finite:
        require(np.isfinite(result).all(), f"Nonfinite values: {name}")
    return result


def edge_set(frame: pd.DataFrame) -> set[tuple[str, str]]:
    return set(zip(frame.tf.astype(str), frame.target.astype(str)))


def load_protocol(path: Path) -> dict:
    protocol = json.loads(Path(path).read_text())
    identifier = protocol.get("protocol", protocol.get("protocol_id"))
    require(isinstance(identifier, str) and bool(identifier),
            "BIO-V8 requires a protocol identifier")
    protocol["protocol"] = identifier
    for field, expected in {"diagnostic_permutation_seed": 20261128, "boost": .3,
                            "minimum_measured_targets_per_tf": 20,
                            "minimum_qualified_tfs": 6}.items():
        require(protocol.get(field) == expected, f"Unexpected BIO-V8 protocol setting: {field}")
    for field in ["biov7_freeze_manifest_sha256", "biov7_protocol_sha256"]:
        value = protocol.get(field, "")
        require(isinstance(value, str) and len(value) == 64 and
                all(c in "0123456789abcdef" for c in value), f"Missing pinned hash: {field}")
    for field in ["no_expression_access", "no_internal_score_access", "untouched_test_remains_sealed",
                  "no_outcome_access", "no_new_prediction_fitting", "no_new_performance_evaluation"]:
        require(protocol.get("constraints", {}).get(field) is True,
                f"Required BIO-V8 access boundary missing: {field}")
    return protocol


def checked_file(path: Path, expected: str, registry: dict, key: str) -> Path:
    require(isinstance(expected, str) and len(expected) == 64,
            f"Missing input checksum: {key}")
    digest = sha(path)
    require(digest == expected, f"Input checksum mismatch: {key}")
    registry[key] = {"path": str(Path(path).resolve()), "sha256": digest}
    return Path(path)


def read_csv(path: Path) -> pd.DataFrame:
    # Parsing everything as text preserves the exact saved legacy score strings.
    # Numeric copies are made only for structural validation and fixed probes.
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def validate_prior(table: pd.DataFrame, expected_edges: set | None = None) -> pd.DataFrame:
    required = {"tf", "target", "tf_gencode_gene_id", "target_gencode_gene_id",
                "motif_prior_score", "dnase_accessibility_mean", "joint_motif_dnase_product",
                "accessibility_source_observed", "direct_tf_binding_verified",
                "independent_regulatory_edge_truth"}
    require(required.issubset(table), "Prior table has missing columns")
    require(len(table) == N_EDGES and not table.duplicated(["tf", "target"]).any(),
            "Prior must contain exactly 1,008 unique edges")
    genes, tfs = set(table.target), set(table.tf)
    universe = {(tf, gene) for tf in tfs for gene in genes if tf != gene}
    require(len(genes) == N_GENES and len(tfs) == N_TFS and tfs.issubset(genes),
            "Prior TF/gene dimensions differ from frozen universe")
    require(edge_set(table) == universe and (expected_edges is None or universe == expected_edges),
            "Prior candidate identities differ from frozen universe")
    require(bool_values(table.accessibility_source_observed, "observed prior").all(),
            "Observed zero accessibility must remain observed")
    require(not bool_values(table.direct_tf_binding_verified, "binding verified").any() and
            not bool_values(table.independent_regulatory_edge_truth, "edge truth").any(),
            "Predictor inputs must not be promoted to independent biological truth")
    result = table.copy()
    for name in ["motif_prior_score", "dnase_accessibility_mean", "joint_motif_dnase_product"]:
        values = numeric(table[name], name)
        require(((values >= 0) & (values <= 1)).all(), f"Prior outside [0,1]: {name}")
        result[name] = values
    require(result.groupby("target").dnase_accessibility_mean.nunique().eq(1).all(),
            "Accessibility must be target constant in this source")
    require(np.allclose(result.joint_motif_dnase_product,
                        result.motif_prior_score * result.dnase_accessibility_mean,
                        rtol=0., atol=1e-12), "Prior product differs from its parents")
    ids = pd.concat([
        table[["tf", "tf_gencode_gene_id"]].rename(columns={"tf": "gene", "tf_gencode_gene_id": "id"}),
        table[["target", "target_gencode_gene_id"]].rename(columns={"target": "gene", "target_gencode_gene_id": "id"})])
    require(ids.groupby("gene").id.nunique().eq(1).all() and
            not ids.id.astype(str).eq("").any(), "Inconsistent stable gene identifiers")
    return result


def validate_contexts(table: pd.DataFrame, prior: pd.DataFrame) -> pd.DataFrame:
    needed = {"tf", "target", "graph_parents", "scaffold_edge", "eligible_action",
              "gain", "expr_rank", "scaffold_selected_on_full_fit", "gain_coefficient_folds"}
    require(needed.issubset(table), "Candidate context columns missing")
    require(len(table) == N_EDGES and not table.duplicated(["tf", "target"]).any() and
            edge_set(table) == edge_set(prior), "Context candidate identities changed")
    out = table.copy()
    scaffold = bool_values(table.scaffold_edge, "scaffold_edge")
    eligible = bool_values(table.eligible_action, "eligible_action")
    require(np.array_equal(scaffold, ~eligible) and int(eligible.sum()) == N_ACTIONS,
            "Context scaffold/action eligibility invalid")
    require(bool_values(table.scaffold_selected_on_full_fit, "scaffold selection").all() and
            np.all(numeric(table.gain_coefficient_folds, "fold count") == 3),
            "Unexpected context construction")
    out["scaffold_edge"], out["eligible_action"] = scaffold, eligible
    tfs = set(prior.tf)
    for target, group in out.groupby("target", sort=True):
        require(group.graph_parents.nunique() == 1, "Target scaffold differs across candidates")
        parents = group.graph_parents.iloc[0].split("|")
        require(len(parents) == 2 and len(set(parents)) == 2 and set(parents).issubset(tfs)
                and target not in parents, "Scaffold must contain two distinct eligible TFs")
        require(set(group.loc[group.scaffold_edge, "tf"]) == set(parents),
                "Scaffold flags disagree with parent identities")
    for name in ["gain", "expr_rank"]:
        require(table.loc[scaffold, name].eq("").all(),
                "Scaffold insertion scores must remain absent")
        vals = numeric(table.loc[eligible, name], name)
        if name == "expr_rank":
            require(((vals >= 0) & (vals <= 1)).all(), "Expression ranks outside [0,1]")
        out[name] = np.nan
        out.loc[eligible, name] = vals
    return out


def winner(tfs: list[str], scores: np.ndarray) -> int:
    require(len(tfs) == len(scores) and len(tfs) >= 2 and np.isfinite(scores).all(),
            "Winner requires at least two finite candidate scores")
    return min(range(len(tfs)), key=lambda i: (-float(scores[i]), tfs[i]))


def beta_boundary(tfs: list[str], expression: np.ndarray, source: np.ndarray) -> tuple:
    """First source-weight crossing of the expression winner; not an outcome.

    A strict larger beta crosses the boundary. At equality, lexical ordering
    decides the winner. ``None`` means no nonnegative weight can change it.
    """
    baseline = winner(tfs, expression)
    crossings = []
    for i, tf in enumerate(tfs):
        if source[i] > source[baseline]:
            boundary = float((expression[baseline] - expression[i]) /
                             (source[i] - source[baseline]))
            crossings.append((boundary, tf, tf < tfs[baseline]))
    return min(crossings, key=lambda item: (item[0], item[1])) if crossings else (None, "", False)


def validate_scores_and_selections(scores: pd.DataFrame, selections: pd.DataFrame,
                                  contexts: pd.DataFrame) -> None:
    required = {"tf", "target", "graph_parents", "gain", "expr_rank", "method", "score",
                "action_score", "constant_target_offsets_removed_for_argmax",
                "global_score_is_not_source_utility_endpoint"}
    require(required.issubset(scores), "Frozen score columns missing")
    require(len(scores) == N_ACTIONS * len(METHODS) and
            not scores.duplicated(["method", "tf", "target"]).any() and
            set(scores.method) == set(METHODS), "Frozen score method/candidate counts differ")
    required_selections = {"method", "target", "graph_parents", "selected_tf", "selected_score",
                           "selected_action_score", "parent_count_before", "inserted_edge_count",
                           "total_parent_count"}
    require(required_selections.issubset(selections), "Frozen selection columns missing")
    require(len(selections) == N_GENES * len(METHODS) and
            not selections.duplicated(["method", "target"]).any() and
            set(selections.method) == set(METHODS), "Frozen selection method/target counts differ")
    require(bool_values(scores.constant_target_offsets_removed_for_argmax, "constant removal").all()
            and bool_values(scores.global_score_is_not_source_utility_endpoint, "score endpoint").all(),
            "Original score semantics changed")
    for col in ["gain", "expr_rank", "score", "action_score"]:
        numeric(scores[col], col)
    for col, expected in [("parent_count_before", 2), ("inserted_edge_count", 1), ("total_parent_count", 3)]:
        require(np.all(numeric(selections[col], col) == expected), "Frozen graph budget changed")
    eligible = contexts.loc[contexts.eligible_action].set_index(["tf", "target"])
    expected_edges = set(eligible.index)
    expected_targets = set(contexts.target)
    for method, group in scores.groupby("method", sort=True):
        require(edge_set(group) == expected_edges, "Methods do not share exact residual candidates")
        ordered = group.set_index(["tf", "target"]).loc[eligible.index]
        require(ordered.graph_parents.tolist() == eligible.graph_parents.tolist(),
                "Frozen score scaffold differs from fit-only context")
        for col in ["gain", "expr_rank"]:
            require(np.array_equal(numeric(ordered[col], col), eligible[col].to_numpy(float)),
                    "Frozen score context values changed")
        chosen = selections.loc[selections.method == method].set_index("target")
        require(set(chosen.index) == expected_targets, "Method lacks a frozen target choice")
        for target, target_group in group.groupby("target", sort=True):
            action = numeric(target_group.action_score, "action_score")
            index = winner(target_group.tf.tolist(), action)
            selected = target_group.iloc[index]
            recorded = chosen.loc[target]
            require(selected.tf == recorded.selected_tf and
                    selected.graph_parents == recorded.graph_parents,
                    "Frozen graph choice does not match exact action-score argmax")
            require(float(recorded.selected_score) == float(selected.score) and
                    float(recorded.selected_action_score) == float(selected.action_score),
                    "Frozen selected score differs from its candidate row")


def diagnose_contexts(contexts: pd.DataFrame, prior: pd.DataFrame, *, seed: int, boost: float
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fixed fit-only probes; no loss, utility, outcome, or model is evaluated."""
    eligible = contexts.loc[contexts.eligible_action].copy()
    merged = eligible.merge(prior[["tf", "target", "motif_prior_score", "dnase_accessibility_mean"]],
                            on=["tf", "target"], validate="one_to_one")
    accessibility = prior.groupby("target", sort=True).dnase_accessibility_mean.first()
    targets = accessibility.index.tolist()
    donor_order = np.random.default_rng(seed).permutation(len(targets))
    donor = {target: targets[index] for target, index in zip(targets, donor_order)}
    mean_accessibility = float(accessibility.mean())  # equal targets, not unequal edge counts
    permutation_rows = [{"target": target, "donor_target": donor[target],
                         "original_accessibility": float(accessibility[target]),
                         "permuted_accessibility": float(accessibility[donor[target]]),
                         "mean_accessibility_equal_target_weight": mean_accessibility,
                         "permutation_seed": seed} for target in targets]
    rows = []
    for target, group in merged.groupby("target", sort=True):
        group = group.sort_values("tf", kind="stable")
        tfs = group.tf.tolist()
        expr = group.expr_rank.to_numpy(float)
        motif = group.motif_prior_score.to_numpy(float)
        a = float(accessibility[target])
        pa = float(accessibility[donor[target]])
        raw = {
            "expression_only": np.zeros(len(group)),
            "motif": motif,
            "global_mean_accessibility_motif": mean_accessibility * motif,
            "matched_accessibility_motif": a * motif,
            "permuted_accessibility_motif": pa * motif,
            "legacy_additive_motif_accessibility": motif + a,
            "accessibility_only": np.full(len(group), a),
        }
        # Preserve exact invariance, including floating-point near ties.
        effective = {name: values.copy() for name, values in raw.items()}
        effective["legacy_additive_motif_accessibility"] = motif.copy()
        effective["accessibility_only"] = np.zeros(len(group))
        for name, values in effective.items():
            if np.ptp(values) == 0.:
                effective[name] = np.zeros(len(group))
        action_scores = {name: expr + boost * values for name, values in effective.items()}
        indices = {name: winner(tfs, scores) for name, scores in action_scores.items()}
        require(indices["legacy_additive_motif_accessibility"] == indices["motif"] and
                indices["accessibility_only"] == indices["expression_only"],
                "Target-constant additive evidence changed an action")
        baseline = indices["expression_only"]
        expression_margin = float(np.sort(expr)[-1] - np.sort(expr)[-2])
        for name in REPRESENTATIONS:
            source = raw[name]
            action_source = effective[name]
            scores = action_scores[name]
            index = indices[name]
            sorted_scores = np.sort(scores)
            boundary, challenger, tie_favors = beta_boundary(tfs, expr, action_source)
            rows.append({
                "target": target, "representation": name,
                "graph_parents": group.graph_parents.iloc[0],
                "candidate_count": len(group), "observed_source_candidates": len(group),
                "source_min": float(source.min()), "source_max": float(source.max()),
                "source_range": float(np.ptp(source)), "source_unique_values": int(np.unique(source).size),
                "effective_action_source_range": float(np.ptp(action_source)),
                "source_structurally_constant": bool(np.ptp(source) == 0.),
                "accessibility": a, "mean_accessibility": mean_accessibility,
                "permutation_donor_target": donor[target], "permuted_accessibility": pa,
                "expression_selected_tf": tfs[baseline], "selected_tf": tfs[index],
                "changed_vs_expression": bool(index != baseline),
                "changed_vs_motif": bool(index != indices["motif"]),
                "changed_vs_global_attenuation": bool(index != indices["global_mean_accessibility_motif"]),
                "expression_margin": expression_margin,
                "action_margin": float(sorted_scores[-1] - sorted_scores[-2]),
                "selected_action_score": float(scores[index]),
                "selected_raw_probe_score": float(expr[index] + boost * source[index]),
                "minimum_beta_boundary": boundary,
                "boundary_challenger_tf": challenger,
                "boundary_tie_favors_challenger": tie_favors,
                "any_nonnegative_source_weight_can_change_expression_action": boundary is not None,
                "beta": boost, "new_model_fitted": False,
                "predictive_usefulness_evaluated": False,
                "independent_third_source": False,
            })
    return pd.DataFrame(rows), pd.DataFrame(permutation_rows)


def eligibility_outputs(contexts: dict[str, pd.DataFrame], prior: pd.DataFrame,
                        minimum_targets: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    tfs = sorted(prior.tf.unique())
    eligibility = prior[["tf", "target", "tf_gencode_gene_id", "target_gencode_gene_id"]].copy()
    eligibility["target_only_gene"] = ~eligibility.target.isin(tfs)
    for transform, context in contexts.items():
        eligible = context[["tf", "target", "eligible_action"]].rename(
            columns={"eligible_action": f"eligible_{transform}"})
        eligibility = eligibility.merge(eligible, on=["tf", "target"], validate="one_to_one")
    eligibility["eligible_both_transforms"] = eligibility[[f"eligible_{t}" for t in TRANSFORMS]].all(axis=1)
    eligibility["all_methods_share_eligibility_within_transform"] = True
    eligibility["scaffold_insertion_scores_fabricated"] = False
    rows = []
    for tf, group in eligibility.groupby("tf", sort=True):
        row = {"tf": tf, "tf_gencode_gene_id": group.tf_gencode_gene_id.iloc[0],
               "full_nonself_candidates": len(group)}
        for transform in TRANSFORMS:
            eligible = group[f"eligible_{transform}"]
            row[f"residual_targets_{transform}"] = int(eligible.sum())
            row[f"residual_target_only_{transform}"] = int((eligible & group.target_only_gene).sum())
            row[f"candidate_count_ge_{minimum_targets}_{transform}"] = int(eligible.sum()) >= minimum_targets
        row["residual_targets_both_transforms"] = int(group.eligible_both_transforms.sum())
        row["residual_target_only_both_transforms"] = int((group.eligible_both_transforms & group.target_only_gene).sum())
        row["external_measured_target_count"] = None
        row["external_assay_qualified"] = False
        rows.append(row)
    return eligibility, pd.DataFrame(rows)


def run(protocol_path: Path, biov7_results: Path, prior_paths: dict[str, Path], output: Path,
        *, biov7_protocol_path: Path = DEFAULT_V7_PROTOCOL) -> dict:
    """Validate the complete allow-list before creating any output directory."""
    output = Path(output)
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())),
            "BIO-V8 output must be empty; existing snapshots are not overwritten")
    protocol = load_protocol(protocol_path)
    registry = {"biov8_protocol": {"path": str(Path(protocol_path).resolve()), "sha256": sha(protocol_path)}}
    v7_path = checked_file(biov7_protocol_path, protocol["biov7_protocol_sha256"], registry, "biov7_protocol")
    v7 = json.loads(v7_path.read_text())
    for key, value in {"transforms": TRANSFORMS, "variants": VARIANTS, "methods": METHODS,
                       "boost": .3, "frozen_genes": N_GENES, "frozen_regulators": N_TFS,
                       "frozen_edges": N_EDGES, "eligible_residual_actions": N_ACTIONS}.items():
        require(v7.get(key) == value, f"Unexpected frozen BIO-V7 config field: {key}")
    manifest_path = checked_file(Path(biov7_results) / "freeze_manifest.json",
                                protocol["biov7_freeze_manifest_sha256"], registry, "biov7/freeze_manifest.json")
    frozen = json.loads(manifest_path.read_text())
    expected_manifest = {
        "stage": "BIOV7_FIT_ONLY_DECISIONS_FROZEN_BEFORE_INTERNAL_SCORE",
        "protocol_sha256": protocol["biov7_protocol_sha256"],
        "candidate_prior_edges": N_EDGES, "eligible_single_edits": N_ACTIONS,
        "n_genes": N_GENES, "n_TFs": N_TFS, "scaffold_edges": 128,
        "inserted_edges_per_method": 64, "total_graph_edges_per_method": 192,
        "fit_cells_used": 394, "internal_score_cells_used": 0, "test_cells_used": 0,
        "test_used": False, "test_expression_analytics_performed": False,
        "perturbation_outcomes_accessed": False,
        "methods": METHODS, "transforms": TRANSFORMS, "variants": VARIANTS,
    }
    for key, value in expected_manifest.items():
        require(frozen.get(key) == value, f"BIO-V7 frozen boundary mismatch: {key}")
    require(set(prior_paths) == set(PRIOR_KEYS), "Exactly three declared prior tables required")
    priors = {}
    for variant, (config_key, input_key) in PRIOR_KEYS.items():
        expected = v7["provenance"][config_key]
        require(frozen["input_sha256"].get(input_key) == expected,
                "BIO-V7 manifest and prior config hash disagree")
        require(protocol.get("input_sha256", {}).get(config_key) == expected,
                "BIO-V8 protocol and BIO-V7 prior hash disagree")
        path = checked_file(prior_paths[variant], expected, registry, input_key)
        prior = validate_prior(read_csv(path), edge_set(priors["primary_mixed"]) if priors else None)
        if priors:
            p0 = priors["primary_mixed"].set_index(["tf", "target"])
            aligned = prior.set_index(["tf", "target"]).loc[p0.index]
            require(np.allclose(aligned.motif_prior_score, p0.motif_prior_score, rtol=0., atol=1e-12),
                    "Motif source changed across processing sensitivity")
            require(aligned[["tf_gencode_gene_id", "target_gencode_gene_id"]].equals(
                p0[["tf_gencode_gene_id", "target_gencode_gene_id"]]), "Gene IDs changed across priors")
        priors[variant] = prior
    require(set(protocol.get("regulators", [])) == set(priors["primary_mixed"].tf),
            "BIO-V8 locked regulators differ from prior identities")
    regulator_ids = priors["primary_mixed"].groupby("tf").tf_gencode_gene_id.first()
    require(protocol.get("regulator_ensembl_ids") ==
            {tf: str(value).split(".")[0] for tf, value in regulator_ids.items()},
            "BIO-V8 exact regulator stable identifiers differ from prior")
    contexts, all_scores, all_selections = {}, [], []
    for transform in TRANSFORMS:
        context_name = f"{transform}/candidate_contexts.csv"
        path = checked_file(Path(biov7_results) / context_name, frozen["output_sha256"].get(context_name),
                            registry, "biov7/" + context_name)
        context = validate_contexts(read_csv(path), priors["primary_mixed"])
        contexts[transform] = context
        for variant in VARIANTS:
            frames = {}
            for name in ["candidate_scores.csv", "frozen_selections.csv"]:
                relative = f"{transform}/{variant}/{name}"
                path = checked_file(Path(biov7_results) / relative, frozen["output_sha256"].get(relative),
                                    registry, "biov7/" + relative)
                frames[name] = read_csv(path)
            scores, selections = frames["candidate_scores.csv"], frames["frozen_selections.csv"]
            validate_scores_and_selections(scores, selections, context)
            for frame in [scores, selections]:
                frame.insert(0, "variant", variant)
                frame.insert(0, "transform", transform)
            scores["external_ranking_score_definition"] = "unaltered_BIOV7_legacy_score_exploratory_across_targets"
            scores["score_is_probability"] = False
            scores["scaffold_edge"] = False
            scores["target_only_gene"] = ~scores.target.isin(set(priors["primary_mixed"].tf))
            scores["primary_external_condition"] = (transform == "log1p_original" and variant == "primary_mixed")
            all_scores.append(scores)
            all_selections.append(selections)
    score_export = pd.concat(all_scores, ignore_index=True)
    selection_export = pd.concat(all_selections, ignore_index=True)
    require(len(score_export) == 49280 and len(selection_export) == 3584,
            "Incomplete frozen score/graph export")
    eligibility, coverage = eligibility_outputs(contexts, priors["primary_mixed"],
                                                protocol["minimum_measured_targets_per_tf"])
    diagnostic_frames, permutation_frames = [], []
    for transform in TRANSFORMS:
        for variant, prior in priors.items():
            diag, permutation = diagnose_contexts(contexts[transform], prior,
                seed=protocol["diagnostic_permutation_seed"], boost=protocol["boost"])
            diag.insert(0, "variant", variant)
            diag.insert(0, "transform", transform)
            diagnostic_frames.append(diag)
            if transform == TRANSFORMS[0]:
                permutation.insert(0, "variant", variant)
                permutation_frames.append(permutation)
    diagnostics = pd.concat(diagnostic_frames, ignore_index=True)
    summary = diagnostics.groupby(["transform", "variant", "representation"], sort=False).agg(
        targets=("target", "size"), structurally_constant_targets=("source_structurally_constant", "sum"),
        changed_vs_expression=("changed_vs_expression", "sum"), changed_vs_motif=("changed_vs_motif", "sum"),
        changed_vs_global_attenuation=("changed_vs_global_attenuation", "sum"),
        median_source_range=("source_range", "median"), median_action_margin=("action_margin", "median")
    ).reset_index()
    summary["structurally_varying_targets"] = summary.targets - summary.structurally_constant_targets
    summary["predictive_usefulness_evaluated"] = False
    summary["biological_validation_performed"] = False
    # Source files and protocol must stay unchanged throughout validation.
    for name, record in registry.items():
        require(sha(Path(record["path"])) == record["sha256"], f"Input changed while reading: {name}")
    output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "frozen_biov7_candidate_scores.csv": score_export,
        "frozen_biov7_graph_selections.csv": selection_export,
        "candidate_eligibility.csv": eligibility,
        "tf_candidate_coverage.csv": coverage,
        "source_action_diagnostics.csv": diagnostics,
        "source_action_summary.csv": summary,
        "accessibility_permutation.csv": pd.concat(permutation_frames, ignore_index=True),
        "frozen_probe_selections.csv": diagnostics[["transform", "variant", "representation", "target",
            "graph_parents", "selected_tf", "selected_action_score", "selected_raw_probe_score",
            "changed_vs_expression", "new_model_fitted", "predictive_usefulness_evaluated"]],
    }
    hashes = {}
    for name, frame in outputs.items():
        path = output / name
        frame.to_csv(path, index=False)
        hashes[name] = sha(path)
    manifest = {
        "protocol": protocol["protocol"], "protocol_sha256": sha(protocol_path),
        "stage": "BIOV8_FROZEN_EXISTING_SCORES_AND_SOURCE_ACTION_QUALIFICATION_NO_OUTCOMES",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha(Path(__file__)), "input_sha256": registry,
        "biov7_freeze_manifest_sha256": sha(manifest_path),
        "frozen_score_rows": len(score_export), "frozen_graph_decisions": len(selection_export),
        "conditions": 56, "residual_candidates_per_condition": N_ACTIONS,
        "source_action_diagnostic_rows": len(diagnostics), "source_action_summary_rows": len(summary),
        "genes": N_GENES, "tfs": N_TFS, "prior_edges": N_EDGES,
        "scaffold_insertion_scores_fabricated": False, "scores_are_probabilities": False,
        "legacy_score_preserved_verbatim": True,
        "external_across_target_ranking_is_exploratory": True,
        "new_models_fitted": 0, "expression_matrices_read": 0,
        "fit_expression_cells_accessed": 0, "internal_score_cells_accessed": 0,
        "untouched_test_cells_accessed": 0, "untouched_test_cells_still_reserved": 133,
        "new_perturbation_outcomes_accessed": False, "new_binding_outcomes_accessed": False,
        "predictive_usefulness_evaluated": False, "biological_validation_performed": False,
        "diagnostic_permutation_seed": protocol["diagnostic_permutation_seed"],
        "permutation_is_biological_replicate": False, "product_is_independent_third_source": False,
        "output_sha256": hashes,
        "interpretation": "Action changes establish decision influence only, not useful or biologically correct changes.",
    }
    write_json(output / "freeze_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--biov7-results", type=Path, required=True)
    parser.add_argument("--primary-prior", type=Path, required=True)
    parser.add_argument("--matched-v1", type=Path, required=True)
    parser.add_argument("--matched-v2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = run(args.protocol, args.biov7_results,
        {"primary_mixed": args.primary_prior, "matched_v1": args.matched_v1, "matched_v2": args.matched_v2},
        args.output)
    print(json.dumps({k: manifest[k] for k in ["stage", "frozen_score_rows", "frozen_graph_decisions",
        "source_action_diagnostic_rows", "new_models_fitted", "expression_matrices_read",
        "untouched_test_cells_accessed", "new_perturbation_outcomes_accessed"]}, indent=2))


if __name__ == "__main__":
    main()
