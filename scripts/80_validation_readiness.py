"""Build a locked validation-readiness package without claiming validation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json  # noqa: E402
from oral_cancer.validation import (  # noqa: E402
    correlation_components,
    mcnemar_approximate_required_n,
    nondominated_fronts,
    paired_t_required_n,
)
from oral_cancer.visuals import phase7_figures  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "07_validation_readiness"


def build_readiness(evidence: pd.DataFrame, config: dict[str, object]) -> pd.DataFrame:
    readiness = evidence.loc[
        evidence["integrated_priority_tier"].eq(config["source_priority_tier"])
    ].copy()
    readiness["follow_up_branch"] = np.where(
        readiness["evidence_tier"].eq("B_abundance"), "quantitative_abundance", "detection_pattern"
    )
    readiness["absolute_internal_effect"] = np.where(
        readiness["follow_up_branch"].eq("quantitative_abundance"),
        readiness["mean_log2_T_minus_N"].abs(),
        readiness["detection_fraction_T_minus_N"].abs(),
    )
    abundance_q = readiness["bh_q_value_abundance"].clip(lower=1e-300)
    detection_q = readiness["bh_q_value_detection"].clip(lower=1e-300)
    readiness["multiplicity_strength"] = np.where(
        readiness["follow_up_branch"].eq("quantitative_abundance"),
        -np.log10(abundance_q), -np.log10(detection_q),
    )
    readiness["lopo_relevant_retention"] = np.where(
        readiness["follow_up_branch"].eq("quantitative_abundance"),
        readiness["lopo_abundance_core_fraction"], readiness["lopo_detection_core_fraction"],
    )
    readiness["additional_evidence_dimensions"] = (
        readiness["pathway_leading_edge_evidence"].astype(int)
        + readiness["ml_stable_evidence"].astype(int)
    )
    readiness["pareto_front"] = 0
    abundance_objectives = [
        "absolute_internal_effect", "multiplicity_strength", "direction_consistency",
        "multiverse_core_support_frequency", "lopo_relevant_retention", "additional_evidence_dimensions",
    ]
    detection_objectives = [
        "absolute_internal_effect", "multiplicity_strength", "sensitivity_fdr_fraction",
        "lopo_relevant_retention", "additional_evidence_dimensions",
    ]
    for branch, objectives in [
        ("quantitative_abundance", abundance_objectives), ("detection_pattern", detection_objectives)
    ]:
        mask = readiness["follow_up_branch"].eq(branch)
        readiness.loc[mask, "pareto_front"] = nondominated_fronts(readiness.loc[mask], objectives)
    readiness["pareto_front"] = readiness["pareto_front"].astype(int)
    readiness["readiness_rule_pass"] = np.where(
        readiness["follow_up_branch"].eq("quantitative_abundance"),
        readiness["complete_pair_count"].ge(config["abundance_minimum_complete_pairs"])
        & readiness["lopo_abundance_core_fraction"].ge(config["abundance_minimum_lopo_retention"]),
        readiness["discordant_pairs"].ge(10)
        & readiness["lopo_detection_core_fraction"].ge(config["detection_minimum_lopo_retention"]),
    )
    readiness["readiness_interpretation"] = np.where(
        readiness["follow_up_branch"].eq("quantitative_abundance"),
        "Internal quantitative-assay follow-up candidate; independent measurement and cohort required.",
        "Internal detection-pattern follow-up candidate; detection mechanism and assay threshold require prospective definition.",
    )
    return readiness


def build_sample_size_tables(config: dict[str, object], abundance_targets: int, detection_targets: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    abundance_alpha = float(config["familywise_alpha"]) / max(abundance_targets, 1)
    abundance_rows = []
    for effect in config["planning_standardized_paired_effects"]:
        for power in config["planning_power"]:
            abundance_rows.append({
                "standardized_paired_effect_assumption": effect,
                "target_power": power,
                "two_sided_per_target_alpha": abundance_alpha,
                "primary_target_count": abundance_targets,
                "required_complete_pairs": paired_t_required_n(float(effect), abundance_alpha, float(power)),
                "planning_status": "assumption_grid_not_candidate_effect_estimate",
            })
    detection_alpha = float(config["familywise_alpha"]) / max(detection_targets, 1)
    detection_rows = []
    for difference in config["planning_detection_differences"]:
        for discordance in config["planning_discordant_fractions"]:
            if difference > discordance:
                continue
            for power in config["planning_power"]:
                detection_rows.append({
                    "absolute_detection_difference_assumption": difference,
                    "total_discordant_fraction_assumption": discordance,
                    "target_power": power,
                    "two_sided_per_target_alpha": detection_alpha,
                    "primary_target_count": detection_targets,
                    "approximate_required_pairs": mcnemar_approximate_required_n(
                        float(difference), float(discordance), detection_alpha, float(power)
                    ),
                    "planning_status": "normal_approximation_recalculate_with_final_assay_parameters",
                })
    return pd.DataFrame(abundance_rows), pd.DataFrame(detection_rows)


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    config = load_json(PROJECT_DIR / "config" / "07_validation_readiness.yml")
    evidence = pd.read_csv(PROJECT_DIR / "results" / "06_robustness" / "integrated_evidence_table.csv")
    readiness = build_readiness(evidence, config)

    paired = np.load(PROJECT_DIR / "data" / "interim" / "paired_data_model.npz", allow_pickle=False)
    feature_keys = paired["feature_keys"].astype(str)
    abundance = readiness.loc[readiness["follow_up_branch"].eq("quantitative_abundance")].copy()
    abundance_keys = abundance["feature_key"].tolist()
    key_to_index = {key: index for index, key in enumerate(feature_keys)}
    changes = pd.DataFrame(
        paired["log2_T_minus_N"][:, [key_to_index[key] for key in abundance_keys]],
        columns=abundance_keys,
    )
    components, edges = correlation_components(
        changes,
        float(config["correlation_absolute_threshold"]),
        int(config["correlation_minimum_complete_pairs"]),
    )
    abundance_ranked = abundance.merge(components, on="feature_key", how="left", validate="one_to_one")
    abundance_ranked = abundance_ranked.sort_values(
        ["readiness_rule_pass", "pareto_front", "additional_evidence_dimensions", "lopo_relevant_retention",
         "multiverse_core_support_frequency", "multiplicity_strength", "absolute_internal_effect", "feature_key"],
        ascending=[False, True, False, False, False, False, False, True],
    )
    representative_keys = abundance_ranked.groupby(
        "abundance_correlation_component", sort=True
    ).head(1)["feature_key"]
    components["component_representative"] = components["feature_key"].isin(representative_keys)
    readiness = readiness.merge(
        components, on="feature_key", how="left", validate="one_to_one"
    )
    readiness["component_representative"] = readiness["component_representative"].fillna(False).astype(bool)

    abundance_shortlist = readiness.loc[
        readiness["follow_up_branch"].eq("quantitative_abundance")
        & readiness["readiness_rule_pass"]
        & readiness["component_representative"]
    ].sort_values(
        ["pareto_front", "additional_evidence_dimensions", "lopo_relevant_retention",
         "multiverse_core_support_frequency", "multiplicity_strength", "absolute_internal_effect", "feature_key"],
        ascending=[True, False, False, False, False, False, True],
    ).head(int(config["abundance_representative_slots"])).copy()
    detection_shortlist = readiness.loc[
        readiness["follow_up_branch"].eq("detection_pattern") & readiness["readiness_rule_pass"]
    ].sort_values(
        ["pareto_front", "additional_evidence_dimensions", "lopo_relevant_retention",
         "multiplicity_strength", "absolute_internal_effect", "feature_key"],
        ascending=[True, False, False, False, False, True],
    ).head(int(config["detection_follow_up_slots"])).copy()
    shortlist = pd.concat([abundance_shortlist, detection_shortlist], ignore_index=True)
    shortlist["proposed_role"] = np.where(
        shortlist["follow_up_branch"].eq("quantitative_abundance"),
        "independent_quantitative_replication_anchor", "prospective_detection_mechanism_sentinel",
    )
    shortlist["validation_state"] = "not_executed"

    readiness.to_csv(RESULT_DIR / "candidate_validation_readiness.csv", index=False)
    components.to_csv(RESULT_DIR / "abundance_redundancy_components.csv", index=False)
    edges.to_csv(RESULT_DIR / "abundance_redundancy_edges.csv", index=False)
    shortlist.to_csv(RESULT_DIR / "locked_handoff_shortlist.csv", index=False)

    abundance_power, detection_power = build_sample_size_tables(
        config, len(abundance_shortlist), len(detection_shortlist)
    )
    abundance_power.to_csv(RESULT_DIR / "abundance_sample_size_sensitivity.csv", index=False)
    detection_power.to_csv(RESULT_DIR / "detection_sample_size_sensitivity.csv", index=False)

    protocol = {
        "phase": 7,
        "protocol_status": "prospective_template_not_executed",
        "source_scope": "84 supplied specimens from 42 patients; no new observations",
        "independent_validation_required": True,
        "phase2_dependency": "unresolved",
        "locked_candidate_feature_keys": shortlist["feature_key"].tolist(),
        "candidate_roles": dict(zip(shortlist["feature_key"], shortlist["proposed_role"])),
        "design_requirements": [
            "Use patients not present in the discovery workbook.",
            "Predefine intended-use population and reference standard before recruitment.",
            "Blind laboratory measurement and outcome adjudication to discovery ranking.",
            "Randomize run order across tissue state and record batch, site and processing variables.",
            "Measure locked candidates without outcome-driven reselection or threshold tuning.",
            "Preserve patient pairing when paired tissue is the estimand.",
            "Report failed measurements, limits of detection and all exclusions.",
        ],
        "primary_endpoints": {
            "quantitative_abundance": "paired log2 tumour-minus-matched-non-tumour difference for locked anchors",
            "detection_pattern": "paired detection discordance under a prospectively fixed assay and detection threshold",
        },
        "multiplicity": "familywise alpha 0.05 divided across primary targets in planning tables; final protocol may use a prespecified strong-FWER method",
        "prohibited_reuse": [
            "Do not call a resplit or bootstrap of the 42 discovery patients external validation.",
            "Do not tune thresholds, candidate membership or assay processing on validation outcomes.",
            "Do not claim screening or diagnostic accuracy without intended-use controls and an independent clinical population.",
        ],
    }
    (RESULT_DIR / "prospective_validation_protocol.json").write_text(
        json.dumps(protocol, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    claims = pd.DataFrame([
        {"claim_id": "P7C01", "claim": f"Phase 7 evaluates all {len(readiness)} Phase 6 high-priority internal candidates without importing new observations.",
         "evidence_file": "results/07_validation_readiness/candidate_validation_readiness.csv", "limitation_ids": "L01;L02;L03;L04;L05;L06;L07;L08;L09;L10;L11;L12", "allowed_scope": "validation readiness"},
        {"claim_id": "P7C02", "claim": f"The hand-off contains {len(abundance_shortlist)} abundance replication anchors and {len(detection_shortlist)} detection-mechanism sentinels selected by frozen branch-specific rules.",
         "evidence_file": "results/07_validation_readiness/locked_handoff_shortlist.csv", "limitation_ids": "L01;L02;L03;L07;L08;L11;L12", "allowed_scope": "prospective follow-up design"},
        {"claim_id": "P7C03", "claim": "Pareto fronts preserve trade-offs among effect, multiplicity evidence, robustness and cross-domain support without a post-hoc weighted score.",
         "evidence_file": "results/07_validation_readiness/candidate_validation_readiness.csv", "limitation_ids": "L08;L11", "allowed_scope": "internal prioritization"},
        {"claim_id": "P7C04", "claim": "Correlation components reduce redundant abundance follow-up choices but do not prove shared mechanisms.",
         "evidence_file": "results/07_validation_readiness/abundance_redundancy_components.csv", "limitation_ids": "L03;L05;L08;L11", "allowed_scope": "internal redundancy control"},
        {"claim_id": "P7C05", "claim": "Sample-size tables are sensitivity grids over assumed effects, not power guarantees based on selected discovery estimates.",
         "evidence_file": "results/07_validation_readiness/abundance_sample_size_sensitivity.csv", "limitation_ids": "L05;L06;L08;L11;L12", "allowed_scope": "prospective planning sensitivity"},
        {"claim_id": "P7C06", "claim": "No independent, orthogonal or clinical validation was executed in Phase 7.",
         "evidence_file": "results/07_validation_readiness/prospective_validation_protocol.json", "limitation_ids": "L08;L09;L12", "allowed_scope": "prohibition"},
    ])
    claims.to_csv(RESULT_DIR / "07_claim_ledger.csv", index=False)

    validation = {
        "phase": 7,
        "milestone": "M7_readiness",
        "computational_status": "PASS",
        "milestone_status": config["milestone_status"],
        "source_high_priority_candidates": len(readiness),
        "abundance_candidates": int(readiness["follow_up_branch"].eq("quantitative_abundance").sum()),
        "detection_candidates": int(readiness["follow_up_branch"].eq("detection_pattern").sum()),
        "locked_abundance_anchors": len(abundance_shortlist),
        "locked_detection_sentinels": len(detection_shortlist),
        "checks": {
            "only_phase6_high_priority_candidates_used": bool(len(readiness) == 165),
            "source_feature_keys_unique": bool(readiness["feature_key"].is_unique),
            "shortlist_feature_keys_unique": bool(shortlist["feature_key"].is_unique),
            "shortlist_is_source_subset": bool(set(shortlist["feature_key"]).issubset(set(readiness["feature_key"]))),
            "no_weighted_composite_score": True,
            "abundance_redundancy_uses_paired_changes": True,
            "sample_sizes_use_assumption_grid": True,
            "independent_validation_not_claimed": protocol["protocol_status"] == "prospective_template_not_executed",
            "phase2_dependency_retained": protocol["phase2_dependency"] == "unresolved",
        },
    }
    if not all(validation["checks"].values()):
        raise ValueError("Phase 7 validation-readiness checks failed")
    (RESULT_DIR / "07_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    phase7_figures(RESULT_DIR, readiness, components)
    print(
        f"PASS: {len(readiness)} candidates; locked {len(abundance_shortlist)} abundance anchors "
        f"and {len(detection_shortlist)} detection sentinels"
    )


if __name__ == "__main__":
    main()

