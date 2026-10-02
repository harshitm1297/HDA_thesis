"""Multiverse, leave-one-patient-out, and pathway robustness for Phase 6."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .pathways import collapse_gene_ranking, enrichment_score, GeneSet
from .phase3 import abundance_inference, detection_inference
from .statistics import benjamini_hochberg


def _scenario_id(representation: str, cohort: str, pairs: int, estimator: str) -> str:
    return f"{representation}__{cohort}__pairs{pairs}__no_imputation__{estimator}"


def build_multiverse(sensitivity: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    scenario_rows, feature_frames, scenario_tables = [], [], {}
    for (representation, cohort, minimum_pairs), frame in sensitivity.groupby(
        ["representation", "cohort", "minimum_complete_pairs"], sort=True
    ):
        for estimator in config["multiverse_estimators"]:
            table = frame.copy()
            p_column = "p_value" if estimator == "moderated_t" else "ordinary_p_value"
            statistic_column = "moderated_t" if estimator == "moderated_t" else "ordinary_t"
            table["multiverse_p_value"] = table[p_column]
            table["multiverse_q_value"] = benjamini_hochberg(table[p_column].to_numpy(dtype=float))
            table["multiverse_statistic"] = table[statistic_column]
            table["core_supported"] = (
                table["eligible"]
                & (table["multiverse_q_value"] <= config["abundance_fdr"])
                & (table["mean_log2_T_minus_N"].abs() >= config["abundance_effect_threshold_log2"])
                & (table["direction_consistency"] >= config["direction_consistency_threshold"])
            )
            eligible = table["eligible"] & table["multiverse_statistic"].notna()
            table["absolute_rank_percentile"] = np.nan
            if eligible.any():
                table.loc[eligible, "absolute_rank_percentile"] = table.loc[eligible, "multiverse_statistic"].abs().rank(
                    method="average", pct=True, ascending=True
                )
            scenario = _scenario_id(str(representation), str(cohort), int(minimum_pairs), estimator)
            table["scenario_id"] = scenario
            table["estimator"] = estimator
            table["missingness_strategy"] = "no_imputation"
            scenario_tables[scenario] = table
            scenario_rows.append({
                "scenario_id": scenario, "representation": representation, "cohort": cohort,
                "minimum_complete_pairs": int(minimum_pairs), "missingness_strategy": "no_imputation",
                "estimator": estimator, "eligible_features": int(eligible.sum()),
                "fdr_supported_features": int((eligible & (table["multiverse_q_value"] <= config["abundance_fdr"])).sum()),
                "core_supported_features": int(table["core_supported"].sum()),
            })
            feature_frames.append(table[[
                "scenario_id", "feature_key", "eligible", "mean_log2_T_minus_N", "multiverse_statistic",
                "multiverse_p_value", "multiverse_q_value", "core_supported", "absolute_rank_percentile",
            ]])
    long = pd.concat(feature_frames, ignore_index=True)
    primary = sensitivity.loc[
        sensitivity["representation"].eq("log2_uncentered")
        & sensitivity["cohort"].eq("all_pairs")
        & sensitivity["minimum_complete_pairs"].eq(30)
    ].set_index("feature_key")["mean_log2_T_minus_N"]
    long["same_primary_sign"] = (
        np.sign(long["mean_log2_T_minus_N"])
        == np.sign(long["feature_key"].map(primary))
    )
    eligible_long = long.loc[long["eligible"]].copy()
    feature_summary = eligible_long.groupby("feature_key").agg(
        multiverse_eligible_runs=("scenario_id", "size"),
        multiverse_sign_agreement=("same_primary_sign", "mean"),
        multiverse_fdr_support_frequency=("multiverse_q_value", lambda value: float((value <= config["abundance_fdr"]).mean())),
        multiverse_core_support_frequency=("core_supported", "mean"),
        multiverse_effect_min=("mean_log2_T_minus_N", "min"),
        multiverse_effect_max=("mean_log2_T_minus_N", "max"),
        multiverse_effect_median=("mean_log2_T_minus_N", "median"),
        multiverse_rank_percentile_sd=("absolute_rank_percentile", "std"),
    ).reset_index()
    feature_summary["multiverse_effect_range"] = feature_summary["multiverse_effect_max"] - feature_summary["multiverse_effect_min"]

    similarity_rows = []
    scenario_names = sorted(scenario_tables)
    for left_index, left_name in enumerate(scenario_names):
        left = scenario_tables[left_name].set_index("feature_key")
        for right_name in scenario_names[left_index + 1:]:
            right = scenario_tables[right_name].set_index("feature_key")
            joined = left[["eligible", "multiverse_statistic", "mean_log2_T_minus_N"]].join(
                right[["eligible", "multiverse_statistic", "mean_log2_T_minus_N"]], lsuffix="_left", rsuffix="_right"
            )
            use = joined["eligible_left"] & joined["eligible_right"]
            similarity_rows.append({
                "scenario_left": left_name, "scenario_right": right_name, "shared_eligible_features": int(use.sum()),
                "statistic_spearman": float(joined.loc[use, "multiverse_statistic_left"].corr(
                    joined.loc[use, "multiverse_statistic_right"], method="spearman"
                )),
                "effect_spearman": float(joined.loc[use, "mean_log2_T_minus_N_left"].corr(
                    joined.loc[use, "mean_log2_T_minus_N_right"], method="spearman"
                )),
            })
    return pd.DataFrame(scenario_rows), feature_summary, pd.DataFrame(similarity_rows), scenario_tables


def leave_one_patient_out(
    delta: np.ndarray,
    states: np.ndarray,
    patient_ids: np.ndarray,
    feature_manifest: pd.DataFrame,
    primary: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    effects, abundance_q, abundance_core, detection_q, detection_core = [], [], [], [], []
    detail_frames = []
    tiered_keys = set(primary.loc[primary["evidence_tier"].ne("not_tiered"), "feature_key"])
    for omitted_index, patient in enumerate(patient_ids):
        keep = np.arange(len(patient_ids)) != omitted_index
        abundance, _ = abundance_inference(delta[keep], feature_manifest, config["lopo_minimum_complete_pairs"], 0.95)
        detection = detection_inference(states[keep], feature_manifest)
        a_core = (
            abundance["eligible"]
            & (abundance["bh_q_value"] <= config["abundance_fdr"])
            & (abundance["mean_log2_T_minus_N"].abs() >= config["abundance_effect_threshold_log2"])
            & (abundance["direction_consistency"] >= config["direction_consistency_threshold"])
        )
        d_core = (
            (detection["bh_q_value"] <= config["detection_fdr"])
            & (detection["detection_fraction_T_minus_N"].abs() >= config["detection_effect_threshold"])
            & (detection["discordant_pairs"] >= config["detection_min_discordant_pairs"])
        )
        effects.append(abundance["mean_log2_T_minus_N"].to_numpy(dtype=float))
        abundance_q.append(abundance["bh_q_value"].to_numpy(dtype=float))
        abundance_core.append(a_core.to_numpy(dtype=bool))
        detection_q.append(detection["bh_q_value"].to_numpy(dtype=float))
        detection_core.append(d_core.to_numpy(dtype=bool))
        detail = pd.DataFrame({
            "omitted_patient": patient, "feature_key": abundance["feature_key"],
            "mean_log2_T_minus_N": abundance["mean_log2_T_minus_N"],
            "abundance_q_value": abundance["bh_q_value"], "abundance_core": a_core,
            "detection_difference": detection["detection_fraction_T_minus_N"],
            "detection_q_value": detection["bh_q_value"], "detection_core": d_core,
        })
        detail_frames.append(detail.loc[detail["feature_key"].isin(tiered_keys)])

    effect_matrix = np.vstack(effects)
    abundance_q_matrix = np.vstack(abundance_q)
    abundance_core_matrix = np.vstack(abundance_core)
    detection_q_matrix = np.vstack(detection_q)
    detection_core_matrix = np.vstack(detection_core)
    primary_effect = primary["mean_log2_T_minus_N"].to_numpy(dtype=float)
    primary_se = primary["moderated_standard_error"].to_numpy(dtype=float)
    change = np.abs(effect_matrix - primary_effect[None, :])
    safe_change = np.where(np.isfinite(change), change, -np.inf)
    worst_index = np.argmax(safe_change, axis=0)
    all_missing_change = ~np.isfinite(change).any(axis=0)
    max_change = np.max(safe_change, axis=0)
    max_change[all_missing_change] = np.nan
    worst_patient = patient_ids[worst_index].astype(object)
    worst_patient[all_missing_change] = None
    valid_effects = np.isfinite(effect_matrix)
    primary_sign = np.sign(primary_effect)
    sign_flips = ((np.sign(effect_matrix) != primary_sign[None, :]) & valid_effects & np.isfinite(primary_effect)[None, :]).sum(axis=0)
    summary = primary[["feature_key", "evidence_tier"]].copy()
    summary["lopo_eligible_runs"] = valid_effects.sum(axis=0)
    effect_min = np.min(np.where(valid_effects, effect_matrix, np.inf), axis=0)
    effect_max = np.max(np.where(valid_effects, effect_matrix, -np.inf), axis=0)
    effect_min[~valid_effects.any(axis=0)] = np.nan
    effect_max[~valid_effects.any(axis=0)] = np.nan
    summary["lopo_effect_min"] = effect_min
    summary["lopo_effect_max"] = effect_max
    summary["lopo_max_absolute_effect_change"] = max_change
    summary["lopo_max_standardized_effect_change"] = np.divide(
        max_change, primary_se, out=np.full_like(max_change, np.nan), where=np.isfinite(primary_se) & (primary_se > 0)
    )
    summary["lopo_worst_patient"] = worst_patient
    summary["lopo_sign_flip_count"] = sign_flips
    summary["lopo_abundance_fdr_fraction"] = np.nanmean(abundance_q_matrix <= config["abundance_fdr"], axis=0)
    summary["lopo_abundance_core_fraction"] = abundance_core_matrix.mean(axis=0)
    summary["lopo_detection_fdr_fraction"] = np.nanmean(detection_q_matrix <= config["detection_fdr"], axis=0)
    summary["lopo_detection_core_fraction"] = detection_core_matrix.mean(axis=0)
    summary["lopo_tier_retention_fraction"] = np.select(
        [summary["evidence_tier"].eq("A_concordant"), summary["evidence_tier"].eq("B_abundance"), summary["evidence_tier"].eq("C_detection_pattern")],
        [np.minimum(summary["lopo_abundance_core_fraction"], summary["lopo_detection_core_fraction"]),
         summary["lopo_abundance_core_fraction"], summary["lopo_detection_core_fraction"]],
        default=np.nan,
    )
    return summary, pd.concat(detail_frames, ignore_index=True)


def pathway_sensitivity(
    scenario_tables: dict[str, pd.DataFrame],
    gene_sets: list[GeneSet],
    primary_ranked: pd.DataFrame,
    minimum: int,
    maximum: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = primary_ranked.set_index("pathway_id")
    gene_set_by_id = {gene_set.stable_id: gene_set for gene_set in gene_sets}
    rows = []
    for scenario, table in sorted(scenario_tables.items()):
        ranking_input = table.copy()
        ranking_input["moderated_t"] = table["multiverse_statistic"]
        ranking, _ = collapse_gene_ranking(ranking_input)
        genes = ranking["gene_symbol"].astype(str).to_numpy()
        weights = ranking["statistic"].to_numpy(dtype=float)
        positions = {gene: index for index, gene in enumerate(genes)}
        for pathway_id, primary_row in primary.iterrows():
            gene_set = gene_set_by_id.get(pathway_id)
            if gene_set is None:
                continue
            overlap = sorted(gene_set.genes.intersection(positions))
            if not minimum <= len(overlap) <= maximum:
                continue
            hits = np.asarray([positions[gene] for gene in overlap], dtype=int)
            score, peak = enrichment_score(weights, hits)
            ordered_hits = np.sort(hits)
            leading_positions = ordered_hits[: peak + 1] if score >= 0 else ordered_hits[peak:]
            leading = set(genes[leading_positions])
            primary_leading = set(str(primary_row["leading_edge_genes"]).split(";"))
            union = leading | primary_leading
            rows.append({
                "scenario_id": scenario, "pathway_id": pathway_id, "pathway_name": primary_row["pathway_name"],
                "enrichment_score": score, "same_primary_direction": bool(np.sign(score) == np.sign(primary_row["normalized_enrichment_score"])),
                "leading_edge_jaccard": len(leading & primary_leading) / len(union) if union else np.nan,
                "leading_edge_gene_count": len(leading),
            })
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["pathway_id", "pathway_name"]).agg(
        pathway_sensitivity_runs=("scenario_id", "size"),
        pathway_direction_agreement=("same_primary_direction", "mean"),
        pathway_leading_edge_jaccard_median=("leading_edge_jaccard", "median"),
        pathway_leading_edge_jaccard_min=("leading_edge_jaccard", "min"),
    ).reset_index()
    return summary, detail
