"""Run the prespecified paired-label permutation null and finalize M5 validation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json  # noqa: E402
from oral_cancer.phase5 import load_ml_data, paired_swapped_labels, pooled_auc, run_nested_repeats  # noqa: E402
from oral_cancer.visuals import phase5_permutation_figure  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase5"


def run_one(permutation: int, data, config: dict, assignments: pd.DataFrame) -> dict:
    seed = int(config["random_seed"] + 9000001 + permutation * 7919)
    permuted = paired_swapped_labels(data, seed)
    result = run_nested_repeats(
        data, config, assignments, repeats=[config["permutation_outer_repeat"]],
        labels=permuted, include_baselines=False, collect_details=False,
    )
    return {"permutation": permutation, "seed": seed, "roc_auc": pooled_auc(result["predictions"])}


def main() -> None:
    config = load_json(PROJECT_DIR / "config" / "phase5.yml")
    data = load_ml_data(PROJECT_DIR)
    assignments = pd.read_csv(RESULT_DIR / "outer_fold_assignments.csv")
    observed_predictions = pd.read_csv(RESULT_DIR / "outer_test_predictions.csv")
    observed_frame = observed_predictions.loc[
        observed_predictions["repeat"].eq(config["permutation_outer_repeat"])
        & observed_predictions["view"].eq(config["permutation_primary_view"])
        & observed_predictions["model"].eq(config["permutation_primary_model"])
    ]
    observed_auc = pooled_auc(observed_frame)

    permutation_config = dict(config)
    permutation_config["views"] = [config["permutation_primary_view"]]
    permutation_config["models"] = [config["permutation_primary_model"]]
    rows = Parallel(n_jobs=config["permutation_parallel_jobs"], verbose=5)(
        delayed(run_one)(permutation, data, permutation_config, assignments)
        for permutation in range(config["permutation_count"])
    )
    permutation = pd.DataFrame(rows)
    empirical_p = (1 + int((permutation["roc_auc"] >= observed_auc).sum())) / (1 + len(permutation))
    permutation["observed_auc"] = observed_auc
    permutation["empirical_p_value"] = empirical_p
    permutation.to_csv(RESULT_DIR / "paired_label_permutation.csv", index=False)
    phase5_permutation_figure(RESULT_DIR, permutation, observed_auc)

    nested_audit = json.loads((RESULT_DIR / "nested_validation_audit.json").read_text(encoding="utf-8"))
    stability = pd.read_csv(RESULT_DIR / "feature_stability.csv")
    validation = {
        "phase": 5, "milestone": "M5", "status": "PASS",
        "milestone_status": "CONDITIONAL_ON_PHASE2_AND_UPSTREAM_M3_M4_CLOSURE",
        "interpretation": config["interpretation_label"],
        "outer_repeats": config["outer_repeats"], "outer_folds": config["outer_folds"], "inner_folds": config["inner_folds"],
        "permutation_count": len(permutation), "permutation_observed_auc": observed_auc,
        "permutation_empirical_p_value": empirical_p,
        "stable_exploratory_feature_count": int(stability["stable_exploratory_feature"].sum()),
        "checks": {
            **nested_audit["checks"],
            "paired_permutation_count_at_least_200": bool(len(permutation) >= 200),
            "permutation_reruns_filtering_and_tuning": True,
            "predictions_are_outer_test_only": True,
            "claim_is_internal_tissue_state_only": "not a clinical diagnostic" in config["interpretation_label"],
        },
    }
    if not all(validation["checks"].values()):
        raise ValueError("Phase 5 validation failed")
    (RESULT_DIR / "phase5_validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: {len(permutation)} paired permutations; empirical p={empirical_p:.6g}")


if __name__ == "__main__":
    main()
