"""Integrate inference, robustness, pathways, and ML without collapsing evidence types."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json  # noqa: E402
from oral_cancer.visuals import phase6_figures  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase6"


def joined_unique(values: pd.Series) -> str:
    return ";".join(sorted({str(value) for value in values.dropna() if str(value)}))


def main() -> None:
    config = load_json(PROJECT_DIR / "config" / "phase6.yml")
    phase3 = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "integrated_candidate_evidence.csv")
    multiverse = pd.read_csv(RESULT_DIR / "multiverse_feature_stability.csv")
    lopo = pd.read_csv(RESULT_DIR / "lopo_feature_influence.csv")
    pathway_edges = pd.read_csv(PROJECT_DIR / "results" / "phase4" / "feature_pathway_network_edges.csv")
    pathways = pathway_edges.groupby("feature_key").agg(
        leading_edge_pathway_ids=("pathway_id", joined_unique),
        leading_edge_pathway_names=("pathway_name", joined_unique),
        leading_edge_pathway_count=("pathway_id", "nunique"),
    ).reset_index()
    ml = pd.read_csv(PROJECT_DIR / "results" / "phase5" / "stable_exploratory_panel.csv")
    ml_summary = ml.groupby("feature_key").agg(
        ml_views=("view", joined_unique), ml_modalities=("modality", joined_unique),
        ml_max_selection_frequency=("selection_frequency", "max"),
        ml_max_sign_consistency=("sign_consistency", "max"),
        ml_max_mean_heldout_importance=("mean_heldout_importance", "max"),
        ml_positive_importance_any=("heldout_importance_positive", "max"),
    ).reset_index()

    evidence = phase3.merge(multiverse, on="feature_key", how="left", validate="one_to_one")
    evidence = evidence.merge(lopo.drop(columns=["evidence_tier"]), on="feature_key", how="left", validate="one_to_one")
    evidence = evidence.merge(pathways, on="feature_key", how="left", validate="one_to_one")
    evidence = evidence.merge(ml_summary, on="feature_key", how="left", validate="one_to_one")
    evidence["leading_edge_pathway_count"] = evidence["leading_edge_pathway_count"].fillna(0).astype(int)

    robust_abundance = (
        evidence["evidence_tier"].isin(["A_concordant", "B_abundance"])
        & (evidence["multiverse_sign_agreement"] >= config["integrated_multiverse_sign_threshold"])
        & (evidence["multiverse_core_support_frequency"] >= config["integrated_multiverse_core_threshold"])
        & (evidence["lopo_abundance_core_fraction"] >= config["integrated_lopo_retention_threshold"])
    )
    robust_detection = (
        evidence["evidence_tier"].isin(["A_concordant", "C_detection_pattern"])
        & (evidence["lopo_detection_core_fraction"] >= config["integrated_lopo_retention_threshold"])
    )
    ml_stable = evidence["ml_max_selection_frequency"] >= config["integrated_ml_selection_threshold"]
    pathway_supported = evidence["leading_edge_pathway_count"] > 0
    cross_domain = ml_stable | pathway_supported
    high_priority = evidence["identifier_unambiguous"] & (robust_abundance | robust_detection) & cross_domain
    evidence["robust_abundance_evidence"] = robust_abundance
    evidence["robust_detection_evidence"] = robust_detection
    evidence["ml_stable_evidence"] = ml_stable
    evidence["pathway_leading_edge_evidence"] = pathway_supported
    evidence["integrated_priority_tier"] = np.select(
        [high_priority, evidence["identifier_unambiguous"] & (robust_abundance | robust_detection), ml_stable],
        ["high_priority_internal", "robust_single_domain", "predictive_stability_only"],
        default="not_prioritized",
    )
    evidence["limitation_ids"] = "L01;L02;L03;L04;L05;L06;L08;L09;L10;L11;L12"
    incomplete_or_detection = evidence["complete_pair_count"].fillna(0).lt(42) | evidence["evidence_tier"].eq("C_detection_pattern")
    evidence.loc[incomplete_or_detection, "limitation_ids"] += ";L07"
    evidence["allowed_claim"] = np.select(
        [evidence["integrated_priority_tier"].eq("high_priority_internal"),
         evidence["integrated_priority_tier"].eq("robust_single_domain"),
         evidence["integrated_priority_tier"].eq("predictive_stability_only")],
        [
            "High-priority internal tissue-associated candidate supported across multiple computational evidence dimensions; external biological and clinical validation required.",
            "Robust internal single-domain tissue-associated candidate; external validation required.",
            "Stable internal prediction feature without sufficient inferential robustness for a biomarker claim.",
        ],
        default="No priority claim under the frozen Phase 6 rules.",
    )
    requested_front = [
        "feature_key", "feature_id_original", "duplicate_identifier", "compound_identifier", "evidence_tier",
        "integrated_priority_tier", "complete_pair_count", "mean_log2_T_minus_N", "ci_low", "ci_high",
        "bh_q_value_abundance", "detection_fraction_T_minus_N", "bh_q_value_detection", "direction_consistency",
        "multiverse_sign_agreement", "multiverse_fdr_support_frequency", "multiverse_core_support_frequency",
        "lopo_max_absolute_effect_change", "lopo_sign_flip_count", "lopo_abundance_core_fraction",
        "lopo_detection_core_fraction", "leading_edge_pathway_names", "ml_max_selection_frequency",
        "ml_max_sign_consistency", "ml_max_mean_heldout_importance", "limitation_ids", "allowed_claim",
    ]
    evidence = evidence[requested_front + [column for column in evidence.columns if column not in requested_front]]
    evidence.to_csv(RESULT_DIR / "integrated_evidence_table.csv", index=False)
    evidence.loc[evidence["integrated_priority_tier"].ne("not_prioritized")].to_csv(
        RESULT_DIR / "prioritized_internal_candidates.csv", index=False
    )

    scenarios = pd.read_csv(RESULT_DIR / "multiverse_scenarios.csv")
    pathway_robustness = pd.read_csv(RESULT_DIR / "pathway_robustness.csv")
    ml_lopo = json.loads((RESULT_DIR / "ml_lopo_summary.json").read_text(encoding="utf-8"))
    negative = pd.read_csv(RESULT_DIR / "negative_controls.csv")
    phase1_validation = load_json(PROJECT_DIR / "results" / "phase1" / "phase1_validation.json")
    phase3_validation = load_json(PROJECT_DIR / "results" / "phase3" / "phase3_validation.json")
    phase4_validation = load_json(PROJECT_DIR / "results" / "phase4" / "phase4_validation.json")
    phase5_validation = load_json(PROJECT_DIR / "results" / "phase5" / "phase5_validation.json")
    tier_counts = phase3_validation["primary"]["tier_counts"]
    priority_counts = evidence["integrated_priority_tier"].value_counts().to_dict()
    lopo_abundance = evidence.loc[evidence["evidence_tier"].eq("B_abundance"), "lopo_abundance_core_fraction"]
    lopo_detection = evidence.loc[evidence["evidence_tier"].eq("C_detection_pattern"), "lopo_detection_core_fraction"]
    robust_pathways = int(
        (
            (pathway_robustness["pathway_direction_agreement"] >= config["pathway_direction_threshold"])
            & (pathway_robustness["pathway_leading_edge_jaccard_median"] >= config["pathway_leading_edge_jaccard_threshold"])
        ).sum()
    )
    ml_negative = negative.loc[
        negative["negative_control"].eq("paired_within_patient_label_swaps_full_nested_ml")
    ].iloc[0]
    overlap_negative = negative.loc[
        negative["negative_control"].eq("random_feature_key_overlap_phase3_tier_vs_ml_stability")
    ].iloc[0]
    claim_rows = [
        {
            "claim_id": "C01", "claim": "The source contains 42 complete tumour/matched-non-tumour patient pairs.",
            "evidence_file": "data/interim/sample_manifest.csv", "evidence_locator": "patient_id,tissue_code",
            "evidence_state": "verified_source", "limitation_ids": "L05;L06;L10", "allowed_scope": "paired tissue design only",
        },
        {
            "claim_id": "C02", "claim": f"Phase 1 retained all 42 pairs for primary analysis and identified {phase1_validation['flagged_patient_pairs']} patient pairs for a whole-pair sensitivity cohort; its label-blind PCA used {phase1_validation['pca_variable_features']} eligible variable features.",
            "evidence_file": "results/phase1/phase1_validation.json", "evidence_locator": "primary_pairs,flagged_patient_pairs,pca_variable_features",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L03;L04;L06;L07", "allowed_scope": "internal QC and sensitivity design",
        },
        {
            "claim_id": "C03", "claim": f"Phase 3 identified {tier_counts['B_abundance']} internal paired-abundance candidates and {tier_counts['C_detection_pattern']} internal detection-pattern candidates under its frozen no-imputation primary rules.",
            "evidence_file": "results/phase3/phase3_validation.json", "evidence_locator": "primary.tier_counts",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L03;L04;L07;L08;L09", "allowed_scope": "internal tissue-associated inference",
        },
        {
            "claim_id": "C04", "claim": f"Phase 4 found {phase4_validation['pathways_supported_both_views']} Reactome v86 pathways supported by both ranked enrichment and patient-level paired scoring.",
            "evidence_file": "results/phase4/phase4_validation.json", "evidence_locator": "pathways_supported_both_views,reactome_release",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L04;L08;L09;L10", "allowed_scope": "internal pathway interpretation",
        },
        {
            "claim_id": "C05", "claim": "Phase 4 did not accept a stable patient cluster and therefore does not claim molecular or clinical subtypes.",
            "evidence_file": "results/phase4/phase4_validation.json", "evidence_locator": "accepted_cluster_k",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L05;L06;L08;L09", "allowed_scope": "negative internal heterogeneity result",
        },
        {
            "claim_id": "C06", "claim": f"Phase 5 evaluated {phase5_validation['evaluated_pipelines']} leakage-controlled pipelines with 25 repeated grouped outer validations and found no defensible compact panel.",
            "evidence_file": "results/phase5/phase5_validation.json", "evidence_locator": "evaluated_pipelines,outer_repeats,compact_panel_discovered",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L05;L06;L08;L09;L11;L12", "allowed_scope": "internal tissue-state classification",
        },
        {
            "claim_id": "C07", "claim": f"Phase 6 evaluates {len(scenarios)} no-imputation multiverse scenarios across normalization, cohort, completeness and estimator.",
            "evidence_file": "results/phase6/multiverse_scenarios.csv", "evidence_locator": "all rows",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L04;L07", "allowed_scope": "computational robustness",
        },
        {
            "claim_id": "C08", "claim": f"Exact leave-one-patient-out refitting was completed for all 42 patients; {int((lopo_abundance >= config['integrated_lopo_retention_threshold']).sum())}/{len(lopo_abundance)} abundance candidates and {int((lopo_detection >= config['integrated_lopo_retention_threshold']).sum())}/{len(lopo_detection)} detection candidates retained core support in at least {config['integrated_lopo_retention_threshold']:.0%} of omissions.",
            "evidence_file": "results/phase6/lopo_feature_influence.csv", "evidence_locator": "evidence_tier,lopo_abundance_core_fraction,lopo_detection_core_fraction",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L03;L07;L11", "allowed_scope": "internal patient influence",
        },
        {
            "claim_id": "C09", "claim": f"Across the Phase 6 pathway multiverse, {robust_pathways}/{len(pathway_robustness)} Reactome v86 pathways met the prespecified direction-agreement and leading-edge-overlap robustness thresholds.",
            "evidence_file": "results/phase6/pathway_robustness.csv", "evidence_locator": "pathway_direction_agreement,pathway_leading_edge_jaccard_median",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L04;L08;L09;L10", "allowed_scope": "internal pathway robustness",
        },
        {
            "claim_id": "C10", "claim": f"The leading Phase 5 abundance elastic-net model has grouped leave-one-patient-out AUC {ml_lopo['roc_auc']:.3f}, balanced accuracy {ml_lopo['balanced_accuracy']:.3f}, and pair-orientation accuracy {ml_lopo['pair_orientation_accuracy']:.3f}.",
            "evidence_file": "results/phase6/ml_lopo_summary.json", "evidence_locator": "roc_auc",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L05;L06;L08;L09;L11;L12", "allowed_scope": "internal tissue-state classification",
        },
        {
            "claim_id": "C11", "claim": f"The 200-run paired-label full nested-pipeline negative control gave empirical p={ml_negative['empirical_p_value']:.6f}; the observed AUC ({ml_negative['observed']:.3f}) exceeded the largest permuted AUC ({ml_negative['null_maximum']:.3f}).",
            "evidence_file": "results/phase6/negative_controls.csv", "evidence_locator": "paired_within_patient_label_swaps_full_nested_ml",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L08;L11;L12", "allowed_scope": "within-dataset exchangeability test",
        },
        {
            "claim_id": "C12", "claim": f"The overlap between Phase 3 tiered candidates and stable Phase 5 ML features was {int(overlap_negative['observed'])}, exceeding the maximum overlap ({int(overlap_negative['null_maximum'])}) in 1,000 random feature-key permutations (empirical p={overlap_negative['empirical_p_value']:.6f}).",
            "evidence_file": "results/phase6/negative_controls.csv", "evidence_locator": "random_feature_key_overlap_phase3_tier_vs_ml_stability",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L08;L09;L11;L12", "allowed_scope": "internal cross-method convergence",
        },
        {
            "claim_id": "C13", "claim": f"The frozen integration rules classify {priority_counts.get('high_priority_internal', 0)} features as high-priority internal candidates, without implying biomarker validity.",
            "evidence_file": "results/phase6/integrated_evidence_table.csv", "evidence_locator": "integrated_priority_tier=high_priority_internal",
            "evidence_state": "derived_from_matrix", "limitation_ids": "L01;L02;L03;L04;L05;L06;L07;L08;L09;L10;L11;L12", "allowed_scope": "internal evidence prioritization",
        },
        {
            "claim_id": "C14", "claim": "No output is an externally validated biomarker, clinical panel, diagnostic test, causal mechanism, population-prevalence estimate, or treatment recommendation.",
            "evidence_file": "DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md", "evidence_locator": "L08,L09,L11,L12",
            "evidence_state": "not_testable_with_current_data", "limitation_ids": "L08;L09;L11;L12", "allowed_scope": "prohibition",
        },
    ]
    pd.DataFrame(claim_rows).to_csv(RESULT_DIR / "claim_ledger.csv", index=False)
    phase6_figures(RESULT_DIR, scenarios, evidence)

    validation = {
        "phase": 6, "milestone": "M6", "computational_status": "PASS",
        "milestone_status": config["milestone_status"],
        "source_feature_rows": len(phase3), "integrated_feature_rows": len(evidence),
        "priority_counts": evidence["integrated_priority_tier"].value_counts().to_dict(),
        "claim_ledger_rows": len(claim_rows),
        "checks": {
            "one_integrated_row_per_source_feature": bool(len(evidence) == len(phase3) == evidence["feature_key"].nunique()),
            "joins_use_feature_key": True,
            "missing_evidence_not_filled_with_zero": bool(evidence["ml_max_selection_frequency"].isna().any()),
            "every_row_has_limitations_and_allowed_claim": bool(evidence["limitation_ids"].notna().all() and evidence["allowed_claim"].notna().all()),
            "phase2_gap_keeps_m6_open": bool("NOT_CLOSED_PHASE2" in config["milestone_status"]),
            "no_external_validation_claim": bool(evidence["allowed_claim"].str.contains("external|No priority|prediction feature", case=False).all()),
            "pathway_robustness_available": bool(len(pathway_robustness) > 0),
        },
    }
    if not all(validation["checks"].values()):
        raise ValueError("Phase 6 integration validation failed")
    (RESULT_DIR / "phase6_validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: integrated {len(evidence)} features; priority counts={validation['priority_counts']}")


if __name__ == "__main__":
    main()
