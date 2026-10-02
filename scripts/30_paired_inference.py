"""Thin command-line entry point for the modular Phase 3 analysis."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json, load_project_data, paired_delta  # noqa: E402
from oral_cancer.inference import abundance_inference, detection_inference, integrate_evidence  # noqa: E402
from oral_cancer.visuals import phase3_figures  # noqa: E402


CONFIG_PATH = PROJECT_DIR / "config" / "03_paired_inference.yml"
RESULT_DIR = PROJECT_DIR / "results" / "03_paired_inference"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    config = load_json(CONFIG_PATH)
    data = load_project_data(PROJECT_DIR)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    primary_delta, primary_patients = paired_delta(data, config["primary_representation"], config["primary_cohort"])
    primary, primary_prior = abundance_inference(
        primary_delta, data.feature_manifest, config["primary_min_complete_pairs"], config["confidence_level"]
    )
    primary.to_csv(RESULT_DIR / "paired_abundance_primary.csv", index=False)

    representations = [config["primary_representation"], *config["sensitivity_representations"]]
    cohorts = [config["primary_cohort"], *config["sensitivity_cohorts"]]
    thresholds = [config["primary_min_complete_pairs"], *config["sensitivity_min_complete_pairs"]]
    sensitivity_tables = []
    sensitivity_priors = []
    for representation in representations:
        for cohort in cohorts:
            delta, patients = paired_delta(data, representation, cohort)
            for threshold in thresholds:
                if representation == config["primary_representation"] and cohort == config["primary_cohort"] and threshold == config["primary_min_complete_pairs"]:
                    continue
                table, prior = abundance_inference(delta, data.feature_manifest, threshold, config["confidence_level"])
                table.insert(0, "representation", representation)
                table.insert(1, "cohort", cohort)
                table.insert(2, "minimum_complete_pairs", threshold)
                sensitivity_tables.append(table)
                sensitivity_priors.append({"representation": representation, "cohort": cohort, "minimum_complete_pairs": threshold, **prior})
    sensitivity = pd.concat(sensitivity_tables, ignore_index=True)
    sensitivity.to_csv(RESULT_DIR / "paired_abundance_sensitivity.csv", index=False)

    paired = np.load(PROJECT_DIR / "data" / "interim" / "paired_data_model.npz")
    detection = detection_inference(paired["detection_state"], data.feature_manifest)
    detection.to_csv(RESULT_DIR / "paired_detection_primary.csv", index=False)
    paired_patient_ids = paired["patient_ids"].astype(str)
    unflagged = set(data.cohort_membership.loc[data.cohort_membership["sensitivity_unflagged_pairs_included"], "patient_id"].astype(str))
    sensitivity_detection = detection_inference(paired["detection_state"][np.asarray([patient in unflagged for patient in paired_patient_ids])], data.feature_manifest)
    sensitivity_detection.to_csv(RESULT_DIR / "paired_detection_sensitivity.csv", index=False)

    evidence = integrate_evidence(primary, detection, sensitivity_detection, sensitivity, config)
    evidence.to_csv(RESULT_DIR / "integrated_candidate_evidence.csv", index=False)
    phase3_figures(RESULT_DIR, primary, detection, evidence)

    eligibility = pd.DataFrame({
        "minimum_complete_pairs": list(range(2, len(primary_patients) + 1)),
        "eligible_feature_count": [int((np.isfinite(primary_delta).sum(axis=0) >= threshold).sum()) for threshold in range(2, len(primary_patients) + 1)],
    })
    eligibility.to_csv(RESULT_DIR / "complete_pair_eligibility_curve.csv", index=False)

    summary = {
        "phase": 3,
        "milestone": "M3",
        "status": "PASS",
        "milestone_status": "CONDITIONAL_PENDING_FULL_PHASE2_IMPUTATION_BENCHMARK",
        "config_sha256": sha256(CONFIG_PATH),
        "primary": {
            "representation": config["primary_representation"],
            "cohort": config["primary_cohort"],
            "patient_pairs": int(len(primary_patients)),
            "minimum_complete_pairs": config["primary_min_complete_pairs"],
            "eligible_features": int(primary["eligible"].sum()),
            "abundance_bh_below_0_05": int((primary["bh_q_value"] <= config["abundance_fdr"]).sum()),
            "abundance_core": int(evidence["abundance_core"].sum()),
            "detection_bh_below_0_05": int((detection["bh_q_value"] <= config["detection_fdr"]).sum()),
            "detection_primary_core": int(evidence["detection_primary_core"].sum()),
            "detection_core": int(evidence["detection_core"].sum()),
            "tier_counts": evidence["evidence_tier"].value_counts().to_dict(),
            "variance_prior": primary_prior,
        },
        "sensitivity_variance_priors": sensitivity_priors,
        "checks": {
            "primary_uses_all_42_pairs": len(primary_patients) == 42,
            "primary_has_no_imputation": config["primary_representation"] == "log2_uncentered",
            "primary_effect_sign_is_T_minus_N": config["effect_sign"] == "tumour minus matched non-tumour",
            "abundance_q_values_in_range": bool(primary["bh_q_value"].dropna().between(0, 1).all()),
            "detection_q_values_in_range": bool(detection["bh_q_value"].dropna().between(0, 1).all()),
            "confidence_intervals_ordered": bool((primary.loc[primary["eligible"], "ci_low"] <= primary.loc[primary["eligible"], "ci_high"]).all()),
            "sensitivity_contains_no_primary_duplicate": len(sensitivity) == 11 * len(primary),
            "duplicate_and_compound_identifiers_not_tier_A": bool((~evidence.loc[evidence["evidence_tier"].eq("A_concordant"), ["duplicate_identifier", "compound_identifier"]].any(axis=1)).all()),
        },
    }
    if not all(summary["checks"].values()):
        raise ValueError("Phase 3 validation failed")
    (RESULT_DIR / "03_validation.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: Phase 3 completed with {summary['primary']['eligible_features']} primary abundance features")


if __name__ == "__main__":
    main()

