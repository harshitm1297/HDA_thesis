"""Paired abundance and detection inference for Phase 3."""

from __future__ import annotations

import math
from dataclasses import asdict

import numpy as np
import pandas as pd

from .statistics import (
    benjamini_hochberg, estimate_variance_prior, exact_paired_binary_pvalue,
    student_t_ppf, student_t_two_sided_p,
)


def abundance_inference(delta: np.ndarray, feature_manifest: pd.DataFrame, minimum_pairs: int, confidence_level: float) -> tuple[pd.DataFrame, dict]:
    count = np.isfinite(delta).sum(axis=0)
    effect = np.divide(np.nansum(delta, axis=0), count, out=np.full(delta.shape[1], np.nan), where=count > 0)
    median = np.full(delta.shape[1], np.nan)
    for index in np.flatnonzero(count > 0):
        median[index] = float(np.median(delta[np.isfinite(delta[:, index]), index]))
    centered = delta - effect[None, :]
    variance = np.divide(np.nansum(centered * centered, axis=0), count - 1, out=np.full(delta.shape[1], np.nan), where=count > 1)
    residual_df = count - 1
    eligible = count >= minimum_pairs
    prior = estimate_variance_prior(variance[eligible], residual_df[eligible])
    posterior_variance = (prior.degrees_freedom * prior.scale + residual_df * variance) / (prior.degrees_freedom + residual_df)
    standard_error = np.sqrt(posterior_variance / count)
    moderated_t = effect / standard_error
    ordinary_standard_error = np.sqrt(variance / count)
    ordinary_t = effect / ordinary_standard_error
    total_df = prior.degrees_freedom + residual_df
    p_value = np.full(len(count), np.nan)
    ordinary_p_value = np.full(len(count), np.nan)
    ci_low = np.full(len(count), np.nan)
    ci_high = np.full(len(count), np.nan)
    for index in np.flatnonzero(eligible):
        p_value[index] = student_t_two_sided_p(float(moderated_t[index]), float(total_df[index]))
        ordinary_p_value[index] = student_t_two_sided_p(float(ordinary_t[index]), float(residual_df[index]))
        critical = student_t_ppf(0.5 + confidence_level / 2.0, float(total_df[index]))
        ci_low[index] = effect[index] - critical * standard_error[index]
        ci_high[index] = effect[index] + critical * standard_error[index]
    q_value = benjamini_hochberg(p_value)
    positive = np.nansum(delta > 0, axis=0)
    negative = np.nansum(delta < 0, axis=0)
    nonzero = positive + negative
    consistency = np.divide(np.maximum(positive, negative), nonzero, out=np.full_like(effect, np.nan), where=nonzero > 0)
    result = feature_manifest[["feature_key", "feature_id_original", "duplicate_identifier", "compound_identifier", "all_abundance_missing"]].copy()
    result["eligible"] = eligible
    result["complete_pair_count"] = count
    result["mean_log2_T_minus_N"] = effect
    result["median_log2_T_minus_N"] = median
    result["moderated_standard_error"] = standard_error
    result["moderated_t"] = moderated_t
    result["moderated_df"] = total_df
    result["ordinary_standard_error"] = ordinary_standard_error
    result["ordinary_t"] = ordinary_t
    result["ordinary_p_value"] = ordinary_p_value
    result["ci_low"] = ci_low
    result["ci_high"] = ci_high
    result["p_value"] = p_value
    result["bh_q_value"] = q_value
    result["positive_pair_count"] = positive
    result["negative_pair_count"] = negative
    result["direction_consistency"] = consistency
    return result, asdict(prior)


def detection_inference(states: np.ndarray, feature_manifest: pd.DataFrame) -> pd.DataFrame:
    both_missing = (states == 0).sum(axis=0)
    n_only = (states == 1).sum(axis=0)
    t_only = (states == 2).sum(axis=0)
    both_detected = (states == 3).sum(axis=0)
    discordant = n_only + t_only
    pairs = states.shape[0]
    p_value = np.asarray([exact_paired_binary_pvalue(int(t), int(n)) for t, n in zip(t_only, n_only)])
    result = feature_manifest[["feature_key", "feature_id_original", "duplicate_identifier", "compound_identifier", "all_abundance_missing"]].copy()
    result["both_missing_pairs"] = both_missing
    result["N_only_pairs"] = n_only
    result["T_only_pairs"] = t_only
    result["both_detected_pairs"] = both_detected
    result["discordant_pairs"] = discordant
    result["detection_fraction_T_minus_N"] = (t_only - n_only) / pairs
    result["exact_p_value"] = p_value
    result["bh_q_value"] = benjamini_hochberg(p_value)
    return result


def integrate_evidence(abundance: pd.DataFrame, detection: pd.DataFrame, detection_sensitivity: pd.DataFrame, sensitivity: pd.DataFrame, config: dict) -> pd.DataFrame:
    evidence = abundance.merge(
        detection[["feature_key", "discordant_pairs", "detection_fraction_T_minus_N", "exact_p_value", "bh_q_value"]],
        on="feature_key", how="left", suffixes=("_abundance", "_detection"),
    )
    evidence = evidence.merge(
        detection_sensitivity[["feature_key", "discordant_pairs", "detection_fraction_T_minus_N", "bh_q_value"]].rename(columns={
            "discordant_pairs": "sensitivity_detection_discordant_pairs",
            "detection_fraction_T_minus_N": "sensitivity_detection_fraction_T_minus_N",
            "bh_q_value": "sensitivity_detection_bh_q_value",
        }),
        on="feature_key", how="left",
    )
    eligible_sensitivity = sensitivity.loc[sensitivity["eligible"]].copy()
    primary_sign = np.sign(abundance.set_index("feature_key")["mean_log2_T_minus_N"])
    eligible_sensitivity["same_primary_sign"] = eligible_sensitivity.apply(
        lambda row: np.sign(row["mean_log2_T_minus_N"]) == primary_sign.get(row["feature_key"], 0), axis=1
    )
    stability = eligible_sensitivity.groupby("feature_key").agg(
        sensitivity_runs_eligible=("eligible", "size"),
        sensitivity_sign_agreement=("same_primary_sign", "mean"),
        sensitivity_fdr_fraction=("bh_q_value", lambda values: float((values <= config["abundance_fdr"]).mean())),
    )
    evidence = evidence.merge(stability, on="feature_key", how="left")
    abundance_core = (
        evidence["eligible"]
        & (evidence["bh_q_value_abundance"] <= config["abundance_fdr"])
        & (evidence["mean_log2_T_minus_N"].abs() >= config["abundance_effect_threshold_log2"])
        & (evidence["direction_consistency"] >= config["direction_consistency_threshold"])
    )
    detection_primary_core = (
        (evidence["bh_q_value_detection"] <= config["detection_fdr"])
        & (evidence["detection_fraction_T_minus_N"].abs() >= config["detection_effect_threshold"])
        & (evidence["discordant_pairs"] >= config["detection_min_discordant_pairs"])
    )
    detection_sensitivity_core = (
        (evidence["sensitivity_detection_bh_q_value"] <= config["detection_fdr"])
        & (evidence["sensitivity_detection_fraction_T_minus_N"].abs() >= config["detection_effect_threshold"])
        & (evidence["sensitivity_detection_discordant_pairs"] >= config["detection_min_discordant_pairs"])
    )
    detection_direction_stable = (
        np.sign(evidence["detection_fraction_T_minus_N"])
        == np.sign(evidence["sensitivity_detection_fraction_T_minus_N"])
    )
    detection_core = detection_primary_core & detection_sensitivity_core & detection_direction_stable
    cross_view_direction_concordant = (
        np.sign(evidence["mean_log2_T_minus_N"])
        == np.sign(evidence["detection_fraction_T_minus_N"])
    )
    stable = (
        (evidence["sensitivity_sign_agreement"].fillna(0) >= config["minimum_sensitivity_sign_agreement"])
        & (evidence["sensitivity_fdr_fraction"].fillna(0) >= config["minimum_sensitivity_fdr_fraction"])
    )
    unambiguous = ~evidence["duplicate_identifier"] & ~evidence["compound_identifier"]
    evidence["abundance_core"] = abundance_core
    evidence["detection_primary_core"] = detection_primary_core
    evidence["detection_core"] = detection_core
    evidence["cross_view_direction_concordant"] = cross_view_direction_concordant
    evidence["identifier_unambiguous"] = unambiguous
    evidence["evidence_tier"] = np.select(
        [abundance_core & detection_core & cross_view_direction_concordant & stable & unambiguous, abundance_core & stable & unambiguous, detection_core & unambiguous],
        ["A_concordant", "B_abundance", "C_detection_pattern"],
        default="not_tiered",
    )
    evidence["permitted_claim"] = np.select(
        [evidence["evidence_tier"].eq("A_concordant"), evidence["evidence_tier"].eq("B_abundance"), evidence["evidence_tier"].eq("C_detection_pattern")],
        [
            "Internally concordant abundance-and-detection candidate; external validation required.",
            "Internal paired-abundance candidate; detection evidence is absent or discordant; external validation required.",
            "Stable tissue-associated detection pattern; quantitative abundance is unresolved; external validation required.",
        ],
        default="No priority claim under the frozen Phase 3 rules.",
    )
    return evidence

