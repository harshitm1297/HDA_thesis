"""Execute Phase 6 ML leave-one-patient-out and negative-control summaries."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json  # noqa: E402
from oral_cancer.ml import fit_classifier, fit_fold_transformer, probability_from_score, raw_model_score, validate_grouped_split  # noqa: E402
from oral_cancer.phase5 import load_ml_data, tune_linear_model  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase6"


def main() -> None:
    config = load_json(PROJECT_DIR / "config" / "phase6.yml")
    phase5_config = load_json(PROJECT_DIR / "config" / "phase5.yml")
    data = load_ml_data(PROJECT_DIR)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    prediction_rows, parameter_rows = [], []
    for patient_number, patient in enumerate(sorted(set(data.patient_ids))):
        test = np.flatnonzero(data.patient_ids == patient)
        train = np.flatnonzero(data.patient_ids != patient)
        validate_grouped_split(train, test, data.patient_ids)
        seed = int(config["random_seed"] + patient_number * 1009)
        parameter, threshold, calibrator, _ = tune_linear_model(
            data, data.labels, train, config["ml_lopo_view"], config["ml_lopo_model"], phase5_config, seed,
        )
        transformer = fit_fold_transformer(
            data.log2_matrix, data.detection, data.labels, train, train, data.feature_keys,
            config["ml_lopo_view"], parameter["selector_k"],
            phase5_config["minimum_training_observed_fraction"], phase5_config["minimum_training_class_observed_fraction"],
        )
        x_train = transformer.transform(data.log2_matrix, data.detection, train)
        x_test = transformer.transform(data.log2_matrix, data.detection, test)
        fitted = fit_classifier(
            config["ml_lopo_model"], x_train, data.labels[train], parameter["C"], parameter["l1_ratio"], seed + 700001,
        )
        probability = probability_from_score(
            config["ml_lopo_model"], raw_model_score(config["ml_lopo_model"], fitted, x_test), calibrator,
        )
        for local, sample_index in enumerate(test):
            prediction_rows.append({
                "left_out_patient": patient, "sample_id": data.sample_ids[sample_index],
                "label": int(data.labels[sample_index]), "probability": float(probability[local]),
                "threshold": threshold, "predicted_label": int(probability[local] >= threshold),
            })
        parameter_rows.append({"left_out_patient": patient, **parameter, "selected_panel_size": int(np.count_nonzero(np.abs(fitted.coef_[0]) > 1e-12))})
        if (patient_number + 1) % 7 == 0:
            print(f"LOPO ML {patient_number + 1}/42", flush=True)
    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(RESULT_DIR / "ml_lopo_predictions.csv", index=False)
    pd.DataFrame(parameter_rows).to_csv(RESULT_DIR / "ml_lopo_parameters.csv", index=False)
    labels = predictions["label"].to_numpy(dtype=int)
    probability = predictions["probability"].to_numpy(dtype=float)
    predicted = predictions["predicted_label"].to_numpy(dtype=int)
    lopo_summary = {
        "view": config["ml_lopo_view"], "model": config["ml_lopo_model"], "patients": 42,
        "roc_auc": float(roc_auc_score(labels, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predicted)),
        "pair_orientation_accuracy": float(np.mean([
            frame.loc[frame["label"].eq(1), "probability"].iloc[0] > frame.loc[frame["label"].eq(0), "probability"].iloc[0]
            for _, frame in predictions.groupby("left_out_patient")
        ])),
    }
    (RESULT_DIR / "ml_lopo_summary.json").write_text(json.dumps(lopo_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    phase5_predictions = pd.read_csv(PROJECT_DIR / "results" / "phase5" / "outer_test_predictions.csv")
    phase5_metrics = pd.read_csv(PROJECT_DIR / "results" / "phase5" / "repeat_level_metrics.csv")
    influence_rows = []
    for (view, model, repeat), frame in phase5_predictions.groupby(["view", "model", "repeat"]):
        original = phase5_metrics.loc[
            phase5_metrics["view"].eq(view) & phase5_metrics["model"].eq(model) & phase5_metrics["repeat"].eq(repeat)
        ].iloc[0]
        for patient, reduced in ((patient, frame.loc[frame["patient_id"].ne(patient)]) for patient in sorted(frame["patient_id"].unique())):
            auc = roc_auc_score(reduced["label"], reduced["probability"])
            ba = balanced_accuracy_score(reduced["label"], reduced["predicted_label"])
            influence_rows.append({
                "view": view, "model": model, "repeat": repeat, "omitted_patient": patient,
                "roc_auc_without_patient": auc, "roc_auc_change": auc - original["roc_auc"],
                "balanced_accuracy_without_patient": ba, "balanced_accuracy_change": ba - original["balanced_accuracy"],
            })
    influence = pd.DataFrame(influence_rows)
    influence.to_csv(RESULT_DIR / "ml_outer_prediction_patient_influence.csv.gz", index=False, compression={"method": "gzip", "compresslevel": 9, "mtime": 0})
    influence.groupby(["view", "model"]).agg(
        maximum_absolute_auc_change=("roc_auc_change", lambda value: float(np.max(np.abs(value)))),
        maximum_absolute_balanced_accuracy_change=("balanced_accuracy_change", lambda value: float(np.max(np.abs(value)))),
        median_absolute_auc_change=("roc_auc_change", lambda value: float(np.median(np.abs(value)))),
    ).reset_index().to_csv(RESULT_DIR / "ml_patient_influence_summary.csv", index=False)

    winner = phase5_metrics.loc[phase5_metrics.groupby("repeat")["roc_auc"].idxmax(), ["repeat", "view", "model", "roc_auc"]]
    robustness = phase5_metrics.groupby(["view", "model"]).agg(
        repeat_count=("repeat", "size"), median_auc=("roc_auc", "median"), minimum_auc=("roc_auc", "min"), maximum_auc=("roc_auc", "max"),
        median_balanced_accuracy=("balanced_accuracy", "median"), minimum_balanced_accuracy=("balanced_accuracy", "min"),
    ).reset_index()
    winner_count = winner.groupby(["view", "model"]).size().rename("repeat_winner_count")
    robustness = robustness.join(winner_count, on=["view", "model"]).fillna({"repeat_winner_count": 0})
    robustness.to_csv(RESULT_DIR / "ml_view_seed_robustness.csv", index=False)

    phase3 = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "integrated_candidate_evidence.csv")
    stable = pd.read_csv(PROJECT_DIR / "results" / "phase5" / "stable_exploratory_panel.csv")
    tiered = set(phase3.loc[phase3["evidence_tier"].ne("not_tiered"), "feature_key"])
    ml_keys = set(stable["feature_key"])
    observed_overlap = len(tiered & ml_keys)
    universe = phase3["feature_key"].to_numpy(dtype=str)
    rng = np.random.default_rng(config["random_seed"] + 9000000)
    null = np.asarray([
        len(tiered & set(rng.choice(universe, size=len(ml_keys), replace=False)))
        for _ in range(config["integration_overlap_permutations"])
    ])
    overlap_p = (1 + int((null >= observed_overlap).sum())) / (1 + len(null))
    phase5_permutation = pd.read_csv(PROJECT_DIR / "results" / "phase5" / "paired_label_permutation.csv")
    controls = pd.DataFrame([
        {
            "negative_control": "paired_within_patient_label_swaps_full_nested_ml",
            "observed": float(phase5_permutation.loc[0, "observed_auc"]),
            "null_median": float(phase5_permutation["roc_auc"].median()),
            "null_maximum": float(phase5_permutation["roc_auc"].max()),
            "permutations": len(phase5_permutation), "empirical_p_value": float(phase5_permutation.loc[0, "empirical_p_value"]),
        },
        {
            "negative_control": "random_feature_key_overlap_phase3_tier_vs_ml_stability",
            "observed": observed_overlap, "null_median": float(np.median(null)), "null_maximum": int(null.max()),
            "permutations": len(null), "empirical_p_value": overlap_p,
        },
    ])
    controls.to_csv(RESULT_DIR / "negative_controls.csv", index=False)
    print(f"PASS: ML LOPO AUC={lopo_summary['roc_auc']:.4f}; stability overlap p={overlap_p:.6g}")


if __name__ == "__main__":
    main()
