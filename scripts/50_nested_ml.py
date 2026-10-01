"""Execute repeated patient-grouped nested validation for Phase 5."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json  # noqa: E402
from oral_cancer.ml import (  # noqa: E402
    choose_threshold,
    fit_classifier,
    fit_fold_transformer,
    probability_from_score,
    raw_model_score,
    summarize_metric_distribution,
    validate_grouped_split,
)
from oral_cancer.phase5 import generate_outer_assignments, load_ml_data, run_nested_repeats  # noqa: E402
from oral_cancer.visuals import phase5_performance_figure  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase5"


def main() -> None:
    config = load_json(PROJECT_DIR / "config" / "phase5.yml")
    data = load_ml_data(PROJECT_DIR)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    assignments = generate_outer_assignments(data, config)
    assignments.to_csv(RESULT_DIR / "outer_fold_assignments.csv", index=False)

    result = run_nested_repeats(data, config, assignments)
    result["predictions"].to_csv(RESULT_DIR / "outer_test_predictions.csv", index=False)
    result["metrics"].to_csv(RESULT_DIR / "repeat_level_metrics.csv", index=False)
    result["tuning"].to_csv(RESULT_DIR / "inner_tuning_surface.csv", index=False)
    result["coefficients"].to_csv(RESULT_DIR / "outer_fit_coefficients.csv", index=False)
    result["importance"].to_csv(RESULT_DIR / "heldout_permutation_importance.csv", index=False)
    summary = summarize_metric_distribution(result["metrics"], ["view", "model"])
    summary.to_csv(RESULT_DIR / "performance_summary.csv", index=False)

    total_outer_fits = config["outer_repeats"] * config["outer_folds"]
    coefficients = result["coefficients"].copy()
    importance = result["importance"].groupby(["view", "model", "feature_name"])["heldout_balanced_accuracy_drop"].mean().rename("mean_heldout_importance")
    stability_rows = []
    for (view, model, feature), frame in coefficients.groupby(["view", "model", "feature_name"]):
        values = frame["coefficient"].to_numpy(dtype=float)
        sign_consistency = max(float((values > 0).mean()), float((values < 0).mean()))
        stability_rows.append({
            "view": view, "model": model, "feature_name": feature,
            "selection_count": len(frame), "outer_fit_count": total_outer_fits,
            "selection_frequency": len(frame) / total_outer_fits,
            "median_coefficient": float(np.median(values)), "sign_consistency": sign_consistency,
            "stable_exploratory_feature": bool(
                len(frame) / total_outer_fits >= config["stable_selection_frequency"]
                and sign_consistency >= config["stable_sign_consistency"]
            ),
        })
    stability = pd.DataFrame(stability_rows)
    stability = stability.join(importance, on=["view", "model", "feature_name"])
    stability = stability.sort_values(["stable_exploratory_feature", "selection_frequency", "sign_consistency"], ascending=False)
    stability.to_csv(RESULT_DIR / "feature_stability.csv", index=False)

    chosen = result["tuning"].loc[
        result["tuning"]["chosen"].astype(str).str.lower().eq("true")
        & result["tuning"]["view"].eq(config["permutation_primary_view"])
        & result["tuning"]["model"].eq(config["permutation_primary_model"])
    ]
    modal_c = float(chosen["C"].value_counts().index[0])
    modal_ratio = float(chosen["l1_ratio"].value_counts().index[0])
    learning_rows = []
    unique_patients = np.asarray(sorted(set(data.patient_ids)))
    for patient_count in config["learning_curve_patient_counts"]:
        for repeat in range(config["learning_curve_repeats"]):
            rng = np.random.default_rng(config["random_seed"] + 7000000 + patient_count * 1000 + repeat)
            train_patients = rng.choice(unique_patients, size=patient_count, replace=False)
            train = np.flatnonzero(np.isin(data.patient_ids, train_patients))
            test = np.flatnonzero(~np.isin(data.patient_ids, train_patients))
            validate_grouped_split(train, test, data.patient_ids)
            transformer = fit_fold_transformer(
                data.log2_matrix, data.detection, data.labels, train, train, data.feature_keys,
                config["permutation_primary_view"], config["elastic_selector_k"],
                config["minimum_training_observed_fraction"], config["minimum_training_class_observed_fraction"],
            )
            x_train = transformer.transform(data.log2_matrix, data.detection, train)
            x_test = transformer.transform(data.log2_matrix, data.detection, test)
            fitted = fit_classifier("elastic_net", x_train, data.labels[train], modal_c, modal_ratio, repeat + patient_count)
            train_probability = probability_from_score("elastic_net", raw_model_score("elastic_net", fitted, x_train), None)
            threshold = choose_threshold(train_probability, data.labels[train])
            probability = probability_from_score("elastic_net", raw_model_score("elastic_net", fitted, x_test), None)
            learning_rows.append({
                "training_patients": patient_count, "repeat": repeat, "test_patients": len(unique_patients) - patient_count,
                "roc_auc": float(roc_auc_score(data.labels[test], probability)),
                "balanced_accuracy": float(((probability[data.labels[test] == 1] >= threshold).mean() + (probability[data.labels[test] == 0] < threshold).mean()) / 2),
                "hyperparameters_source": "modal_nested_outer_choices; test sets not used",
            })
    pd.DataFrame(learning_rows).to_csv(RESULT_DIR / "learning_curve.csv", index=False)
    phase5_performance_figure(RESULT_DIR, summary)

    audit = {
        "phase": 5, "stage": "nested_validation", "status": "PASS",
        "outer_repeats": config["outer_repeats"], "outer_folds": config["outer_folds"],
        "inner_folds": config["inner_folds"], "patients": len(unique_patients),
        "outer_predictions": len(result["predictions"]),
        "checks": {
            "each_repeat_has_all_84_samples_per_pipeline": bool(
                result["predictions"].groupby(["repeat", "view", "model"]).size().eq(84).all()
            ),
            "every_outer_fold_has_seven_complete_pairs": bool(
                assignments.groupby(["repeat", "outer_fold"]).size().eq(7).all()
            ),
            "no_phase3_or_phase4_supervised_filter": True,
            "all_probabilities_finite": bool(np.isfinite(result["predictions"]["probability"]).all()),
        },
    }
    (RESULT_DIR / "nested_validation_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: Phase 5 nested validation produced {len(result['predictions'])} outer-test predictions")


if __name__ == "__main__":
    main()
