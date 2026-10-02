"""Repeated grouped nested validation used by the Phase 5 entry points."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

from .data import load_project_data, log2_representation
from .ml import (
    apply_pca_transform,
    choose_threshold,
    classification_metrics,
    fit_classifier,
    fit_fold_transformer,
    fit_pca_transform,
    fit_platt,
    make_group_folds,
    parameter_grid,
    probability_from_score,
    raw_model_score,
    validate_grouped_split,
)


@dataclass(frozen=True)
class MLData:
    log2_matrix: np.ndarray
    detection: np.ndarray
    labels: np.ndarray
    sample_ids: np.ndarray
    patient_ids: np.ndarray
    feature_keys: np.ndarray


def load_ml_data(project_dir: Path) -> MLData:
    data = load_project_data(project_dir)
    manifest = data.sample_manifest.set_index("sample_id").loc[data.sample_ids]
    labels = manifest["tissue_code"].eq("T").to_numpy(dtype=int)
    patient_ids = manifest["patient_id"].to_numpy(dtype=str)
    return MLData(
        log2_matrix=log2_representation(data, "log2_uncentered"),
        detection=data.detected.astype(bool), labels=labels, sample_ids=data.sample_ids,
        patient_ids=patient_ids, feature_keys=data.feature_keys,
    )


def generate_outer_assignments(data: MLData, config: dict) -> pd.DataFrame:
    rows = []
    unique_patients = np.asarray(sorted(set(data.patient_ids)))
    for repeat in range(config["outer_repeats"]):
        seed = int(config["random_seed"] + repeat * 1009)
        mapping = make_group_folds(unique_patients, config["outer_folds"], seed)
        for patient in unique_patients:
            rows.append({"repeat": repeat, "outer_seed": seed, "patient_id": patient, "outer_fold": mapping[patient]})
    return pd.DataFrame(rows)


def _inner_splits(data: MLData, outer_train: np.ndarray, config: dict, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    train_patients = data.patient_ids[outer_train]
    mapping = make_group_folds(np.unique(train_patients), config["inner_folds"], seed)
    splits = []
    for fold in range(config["inner_folds"]):
        validation_patients = {patient for patient, assigned in mapping.items() if assigned == fold}
        validation = outer_train[np.isin(data.patient_ids[outer_train], list(validation_patients))]
        training = outer_train[~np.isin(data.patient_ids[outer_train], list(validation_patients))]
        validate_grouped_split(training, validation, data.patient_ids)
        splits.append((training, validation))
    return splits


def tune_linear_model(
    data: MLData,
    labels: np.ndarray,
    outer_train: np.ndarray,
    view: str,
    model: str,
    config: dict,
    seed: int,
) -> tuple[dict, float, object | None, pd.DataFrame]:
    params = parameter_grid(model, config)
    splits = _inner_splits(data, outer_train, config, seed)
    maximum_features = max(item["selector_k"] for item in params)
    prepared = []
    for split_number, (training, validation) in enumerate(splits):
        transformer = fit_fold_transformer(
            data.log2_matrix, data.detection, labels, training, outer_train, data.feature_keys,
            view, maximum_features, config["minimum_training_observed_fraction"],
            config["minimum_training_class_observed_fraction"],
        )
        prepared.append((
            training, validation, transformer,
            transformer.transform(data.log2_matrix, data.detection, training),
            transformer.transform(data.log2_matrix, data.detection, validation),
        ))

    surfaces, oof_scores = [], {}
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    for parameter_index, parameter in enumerate(params):
        fold_scores, panel_sizes = [], []
        oof = np.full(len(labels), np.nan)
        for split_number, (training, validation, transformer, x_training_max, x_validation_max) in enumerate(prepared):
            k = min(parameter["selector_k"], x_training_max.shape[1])
            fitted = fit_classifier(
                model, x_training_max[:, :k], labels[training], parameter["C"], parameter["l1_ratio"],
                seed + parameter_index * 101 + split_number,
            )
            score = raw_model_score(model, fitted, x_validation_max[:, :k])
            oof[validation] = score
            predicted = score >= (0.5 if model == "elastic_net" else 0.0)
            fold_scores.append(float(balanced_accuracy_score(labels[validation], predicted)))
            panel_sizes.append(int(np.count_nonzero(np.abs(fitted.coef_[0]) > 1e-12)))
        row = {
            "parameter_index": parameter_index, **parameter,
            "mean_inner_balanced_accuracy": float(np.mean(fold_scores)),
            "mean_selected_panel_size": float(np.mean(panel_sizes)),
        }
        surfaces.append(row)
        oof_scores[parameter_index] = oof

    surface = pd.DataFrame(surfaces)
    surface = surface.sort_values(
        ["mean_inner_balanced_accuracy", "C", "mean_selected_panel_size", "l1_ratio", "selector_k"],
        ascending=[False, True, True, False, True], na_position="last",
    ).reset_index(drop=True)
    best = surface.iloc[0]
    parameter = params[int(best["parameter_index"])]
    for row in surfaces:
        row["chosen"] = bool(row["parameter_index"] == int(best["parameter_index"]))
    raw_oof = oof_scores[int(best["parameter_index"])][outer_train]
    calibrator = fit_platt(raw_oof, labels[outer_train], seed + 900001) if model == "linear_svm" else None
    probability = probability_from_score(model, raw_oof, calibrator)
    threshold = choose_threshold(probability, labels[outer_train])
    return parameter, threshold, calibrator, pd.DataFrame(surfaces)


def tune_pca_baseline(
    data: MLData, labels: np.ndarray, outer_train: np.ndarray, config: dict, seed: int,
) -> tuple[float, float, pd.DataFrame]:
    splits = _inner_splits(data, outer_train, config, seed)
    prepared = []
    for split_number, (training, validation) in enumerate(splits):
        state, x_training = fit_pca_transform(
            data.log2_matrix, labels, training, outer_train, config["pca_candidate_features"],
            config["pca_components"], seed + split_number,
        )
        prepared.append((training, validation, x_training, apply_pca_transform(state, data.log2_matrix, validation)))
    rows, predictions = [], {}
    for parameter_index, c_value in enumerate(config["c_grid"]):
        oof = np.full(len(labels), np.nan)
        fold_scores = []
        for split_number, (training, validation, x_training, x_validation) in enumerate(prepared):
            fitted = LogisticRegression(C=float(c_value), solver="lbfgs", max_iter=3000, random_state=seed + split_number)
            fitted.fit(x_training, labels[training])
            probability = fitted.predict_proba(x_validation)[:, 1]
            oof[validation] = probability
            fold_scores.append(float(balanced_accuracy_score(labels[validation], probability >= 0.5)))
        rows.append({"parameter_index": parameter_index, "C": float(c_value), "mean_inner_balanced_accuracy": float(np.mean(fold_scores))})
        predictions[parameter_index] = oof
    surface = pd.DataFrame(rows).sort_values(["mean_inner_balanced_accuracy", "C"], ascending=[False, True])
    best_index = int(surface.iloc[0]["parameter_index"])
    for row in rows:
        row["chosen"] = bool(row["parameter_index"] == best_index)
    probability = predictions[best_index][outer_train]
    return float(config["c_grid"][best_index]), choose_threshold(probability, labels[outer_train]), pd.DataFrame(rows)


def _heldout_importance(
    model: str, fitted, calibrator, transformer, x_test: np.ndarray, labels: np.ndarray,
    threshold: float, feature_names: np.ndarray, repeats: int, seed: int,
) -> list[dict]:
    raw = raw_model_score(model, fitted, x_test)
    probability = probability_from_score(model, raw, calibrator)
    baseline = balanced_accuracy_score(labels, probability >= threshold)
    coefficients = fitted.coef_[0]
    order = np.argsort(np.abs(coefficients))[::-1][: min(len(coefficients), 25)]
    rng = np.random.default_rng(seed)
    rows = []
    for column in order:
        values = []
        for _ in range(repeats):
            altered = x_test.copy()
            altered[:, column] = altered[rng.permutation(len(altered)), column]
            changed = probability_from_score(model, raw_model_score(model, fitted, altered), calibrator)
            values.append(float(baseline - balanced_accuracy_score(labels, changed >= threshold)))
        rows.append({"feature_name": feature_names[column], "heldout_balanced_accuracy_drop": float(np.mean(values))})
    return rows


def run_nested_repeats(
    data: MLData,
    config: dict,
    assignments: pd.DataFrame,
    repeats: list[int] | None = None,
    labels: np.ndarray | None = None,
    include_baselines: bool = True,
    collect_details: bool = True,
) -> dict[str, pd.DataFrame]:
    labels = data.labels if labels is None else np.asarray(labels, dtype=int)
    repeats = list(range(config["outer_repeats"])) if repeats is None else repeats
    prediction_rows, tuning_frames, coefficient_rows, importance_rows = [], [], [], []

    for repeat in repeats:
        repeat_assignment = assignments.loc[assignments["repeat"] == repeat].set_index("patient_id")["outer_fold"].to_dict()
        for outer_fold in range(config["outer_folds"]):
            test_patients = {patient for patient, fold in repeat_assignment.items() if fold == outer_fold}
            test = np.flatnonzero(np.isin(data.patient_ids, list(test_patients)))
            train = np.flatnonzero(~np.isin(data.patient_ids, list(test_patients)))
            validate_grouped_split(train, test, data.patient_ids)
            seed = int(config["random_seed"] + repeat * 100003 + outer_fold * 1009)

            for view in config["views"]:
                for model in config["models"]:
                    parameter, threshold, calibrator, surface = tune_linear_model(data, labels, train, view, model, config, seed)
                    transformer = fit_fold_transformer(
                        data.log2_matrix, data.detection, labels, train, train, data.feature_keys, view,
                        parameter["selector_k"], config["minimum_training_observed_fraction"],
                        config["minimum_training_class_observed_fraction"],
                    )
                    x_train = transformer.transform(data.log2_matrix, data.detection, train, parameter["selector_k"])
                    x_test = transformer.transform(data.log2_matrix, data.detection, test, parameter["selector_k"])
                    fitted = fit_classifier(model, x_train, labels[train], parameter["C"], parameter["l1_ratio"], seed + 700001)
                    raw = raw_model_score(model, fitted, x_test)
                    probability = probability_from_score(model, raw, calibrator)
                    for local, sample_index in enumerate(test):
                        prediction_rows.append({
                            "repeat": repeat, "outer_fold": outer_fold, "view": view, "model": model,
                            "sample_id": data.sample_ids[sample_index], "patient_id": data.patient_ids[sample_index],
                            "label": int(labels[sample_index]), "probability": float(probability[local]),
                            "threshold": float(threshold), "predicted_label": int(probability[local] >= threshold),
                        })
                    if collect_details:
                        surface.insert(0, "model", model); surface.insert(0, "view", view)
                        surface.insert(0, "outer_fold", outer_fold); surface.insert(0, "repeat", repeat)
                        tuning_frames.append(surface)
                        coefficients = fitted.coef_[0]
                        for name, coefficient in zip(transformer.feature_names[: len(coefficients)], coefficients):
                            if model == "linear_svm" or abs(float(coefficient)) > 1e-12:
                                coefficient_rows.append({
                                    "repeat": repeat, "outer_fold": outer_fold, "view": view, "model": model,
                                    "feature_name": name, "coefficient": float(coefficient),
                                    "selected_panel_size": int(np.count_nonzero(np.abs(coefficients) > 1e-12)),
                                })
                        for row in _heldout_importance(
                            model, fitted, calibrator, transformer, x_test, labels[test], threshold,
                            transformer.feature_names, config["importance_repeats"], seed + 800003,
                        ):
                            importance_rows.append({"repeat": repeat, "outer_fold": outer_fold, "view": view, "model": model, **row})

            if include_baselines:
                # Coverage-only baseline: one training-standardized, prespecified feature.
                coverage = data.detection.mean(axis=1)
                mean, scale = float(coverage[train].mean()), float(coverage[train].std())
                scale = scale if scale > 0 else 1.0
                x_train = ((coverage[train] - mean) / scale).reshape(-1, 1)
                x_test = ((coverage[test] - mean) / scale).reshape(-1, 1)
                coverage_model = LogisticRegression(C=1.0, solver="lbfgs", random_state=seed).fit(x_train, labels[train])
                train_probability = coverage_model.predict_proba(x_train)[:, 1]
                threshold = choose_threshold(train_probability, labels[train])
                probability = coverage_model.predict_proba(x_test)[:, 1]
                for local, sample_index in enumerate(test):
                    prediction_rows.append({
                        "repeat": repeat, "outer_fold": outer_fold, "view": "coverage", "model": "logistic_baseline",
                        "sample_id": data.sample_ids[sample_index], "patient_id": data.patient_ids[sample_index],
                        "label": int(labels[sample_index]), "probability": float(probability[local]),
                        "threshold": float(threshold), "predicted_label": int(probability[local] >= threshold),
                    })

                # Fold-local unsupervised abundance PCA baseline.
                c_value, threshold, surface = tune_pca_baseline(data, labels, train, config, seed + 300007)
                state, pca_train = fit_pca_transform(
                    data.log2_matrix, labels, train, train, config["pca_candidate_features"],
                    config["pca_components"], seed + 300011,
                )
                pca_model = LogisticRegression(C=c_value, solver="lbfgs", max_iter=3000, random_state=seed).fit(pca_train, labels[train])
                probability = pca_model.predict_proba(apply_pca_transform(state, data.log2_matrix, test))[:, 1]
                for local, sample_index in enumerate(test):
                    prediction_rows.append({
                        "repeat": repeat, "outer_fold": outer_fold, "view": "abundance_pca", "model": "logistic_baseline",
                        "sample_id": data.sample_ids[sample_index], "patient_id": data.patient_ids[sample_index],
                        "label": int(labels[sample_index]), "probability": float(probability[local]),
                        "threshold": float(threshold), "predicted_label": int(probability[local] >= threshold),
                    })
                if collect_details:
                    surface.insert(0, "model", "logistic_baseline"); surface.insert(0, "view", "abundance_pca")
                    surface.insert(0, "outer_fold", outer_fold); surface.insert(0, "repeat", repeat)
                    tuning_frames.append(surface)

                for local, sample_index in enumerate(test):
                    prediction_rows.append({
                        "repeat": repeat, "outer_fold": outer_fold, "view": "intercept", "model": "null_baseline",
                        "sample_id": data.sample_ids[sample_index], "patient_id": data.patient_ids[sample_index],
                        "label": int(labels[sample_index]), "probability": 0.5, "threshold": 0.5,
                        "predicted_label": 1,
                    })

    predictions = pd.DataFrame(prediction_rows)
    metric_rows = []
    for (repeat, view, model), frame in predictions.groupby(["repeat", "view", "model"]):
        metrics = classification_metrics(frame["label"].to_numpy(), frame["probability"].to_numpy(), 0.5)
        labels_array = frame["label"].to_numpy(dtype=int)
        predicted_array = frame["predicted_label"].to_numpy(dtype=int)
        sensitivity = float(predicted_array[labels_array == 1].mean())
        specificity = float((predicted_array[labels_array == 0] == 0).mean())
        metrics.update({"sensitivity": sensitivity, "specificity": specificity, "balanced_accuracy": (sensitivity + specificity) / 2})
        metric_rows.append({"repeat": repeat, "view": view, "model": model, **metrics})
    return {
        "predictions": predictions,
        "metrics": pd.DataFrame(metric_rows),
        "tuning": pd.concat(tuning_frames, ignore_index=True) if tuning_frames else pd.DataFrame(),
        "coefficients": pd.DataFrame(coefficient_rows),
        "importance": pd.DataFrame(importance_rows),
    }


def paired_swapped_labels(data: MLData, seed: int) -> np.ndarray:
    labels = data.labels.copy()
    rng = np.random.default_rng(seed)
    for patient in sorted(set(data.patient_ids)):
        indices = np.flatnonzero(data.patient_ids == patient)
        if rng.random() < 0.5:
            labels[indices] = labels[indices[::-1]]
    return labels


def pooled_auc(predictions: pd.DataFrame) -> float:
    return float(roc_auc_score(predictions["label"], predictions["probability"]))
