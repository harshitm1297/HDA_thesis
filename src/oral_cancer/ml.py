"""Leakage-guarded feature processing and linear predictive models for Phase 5."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)
from sklearn.svm import LinearSVC


def make_group_folds(patient_ids: np.ndarray, folds: int, seed: int) -> dict[str, int]:
    unique = np.asarray(sorted(set(patient_ids.astype(str))))
    if len(unique) % folds != 0:
        raise ValueError("Patient count must divide evenly across the prespecified folds")
    shuffled = unique.copy()
    np.random.default_rng(seed).shuffle(shuffled)
    return {patient: int(index % folds) for index, patient in enumerate(shuffled)}


def validate_grouped_split(train: np.ndarray, test: np.ndarray, patient_ids: np.ndarray) -> None:
    train, test = np.asarray(train, dtype=int), np.asarray(test, dtype=int)
    if np.intersect1d(train, test).size:
        raise ValueError("Leakage: sample index occurs in both training and test sets")
    train_patients = set(patient_ids[train].astype(str))
    test_patients = set(patient_ids[test].astype(str))
    if train_patients & test_patients:
        raise ValueError("Leakage: patient occurs in both training and test sets")
    for patient in set(patient_ids.astype(str)):
        indices = np.flatnonzero(patient_ids.astype(str) == patient)
        membership = set(np.isin(indices, test).tolist())
        if len(membership) > 1:
            raise ValueError("Leakage: paired specimens received different fold assignments")


def validate_fit_scope(fit_indices: np.ndarray, allowed_training_indices: np.ndarray) -> None:
    fit_set = set(np.asarray(fit_indices, dtype=int).tolist())
    allowed = set(np.asarray(allowed_training_indices, dtype=int).tolist())
    if not fit_set or not fit_set <= allowed:
        raise ValueError("Leakage: transformer fit includes a non-training sample")


def reject_external_supervised_filter(feature_source: str) -> None:
    if feature_source != "fold_local_training_labels":
        raise ValueError("Leakage: supervised feature filters must be learned inside the training fold")


def reject_outer_test_hyperparameter_selection(outer_test_metric: float | None) -> None:
    if outer_test_metric is not None:
        raise ValueError("Leakage: outer-test performance cannot select hyperparameters")


@dataclass(frozen=True)
class FoldTransformer:
    view: str
    source_indices: np.ndarray
    source_types: np.ndarray
    feature_names: np.ndarray
    medians: np.ndarray
    means: np.ndarray
    scales: np.ndarray
    training_sample_count: int

    def transform(self, log2_matrix: np.ndarray, detection: np.ndarray, indices: np.ndarray, k: int | None = None) -> np.ndarray:
        take = len(self.source_indices) if k is None else min(int(k), len(self.source_indices))
        columns = []
        for source, kind, median in zip(self.source_indices[:take], self.source_types[:take], self.medians[:take]):
            if kind == "A":
                values = log2_matrix[np.asarray(indices, dtype=int), source].astype(float)
                values = np.where(np.isfinite(values), values, median)
            else:
                values = detection[np.asarray(indices, dtype=int), source].astype(float)
            columns.append(values)
        if not columns:
            raise ValueError(f"No eligible features for {self.view}")
        matrix = np.column_stack(columns)
        return (matrix - self.means[:take]) / self.scales[:take]


def fit_fold_transformer(
    log2_matrix: np.ndarray,
    detection: np.ndarray,
    labels: np.ndarray,
    fit_indices: np.ndarray,
    allowed_training_indices: np.ndarray,
    feature_keys: np.ndarray,
    view: str,
    maximum_features: int,
    minimum_observed: float,
    minimum_class_observed: float,
    feature_source: str = "fold_local_training_labels",
) -> FoldTransformer:
    validate_fit_scope(fit_indices, allowed_training_indices)
    reject_external_supervised_filter(feature_source)
    fit_indices = np.asarray(fit_indices, dtype=int)
    y = labels[fit_indices]
    if set(y.tolist()) != {0, 1}:
        raise ValueError("Both tissue labels are required in every training fold")

    records: list[tuple[float, str, int, str, float, float, float]] = []
    if view in {"abundance", "combined"}:
        training = log2_matrix[fit_indices]
        finite = np.isfinite(training)
        eligible = (
            (finite.mean(axis=0) >= minimum_observed)
            & (finite[y == 0].mean(axis=0) >= minimum_class_observed)
            & (finite[y == 1].mean(axis=0) >= minimum_class_observed)
        )
        for source in np.flatnonzero(eligible):
            observed = training[:, source]
            median = float(np.nanmedian(observed))
            filled = np.where(np.isfinite(observed), observed, median)
            mean = float(filled.mean())
            scale = float(filled.std(ddof=0))
            if not np.isfinite(scale) or scale <= 0:
                continue
            standardized = (filled - mean) / scale
            score = abs(float(standardized[y == 1].mean() - standardized[y == 0].mean()))
            name = f"A:{feature_keys[source]}"
            records.append((score, name, int(source), "A", median, mean, scale))

    if view in {"detection", "combined"}:
        training = detection[fit_indices].astype(float)
        prevalence = training.mean(axis=0)
        eligible = (prevalence > 0) & (prevalence < 1)
        for source in np.flatnonzero(eligible):
            values = training[:, source]
            mean = float(values.mean())
            scale = float(values.std(ddof=0))
            if scale <= 0:
                continue
            standardized = (values - mean) / scale
            score = abs(float(standardized[y == 1].mean() - standardized[y == 0].mean()))
            name = f"D:{feature_keys[source]}"
            records.append((score, name, int(source), "D", 0.0, mean, scale))

    records.sort(key=lambda item: (-item[0], item[1]))
    records = records[: int(maximum_features)]
    if not records:
        raise ValueError(f"No fold-local eligible features for {view}")
    return FoldTransformer(
        view=view,
        source_indices=np.asarray([item[2] for item in records], dtype=int),
        source_types=np.asarray([item[3] for item in records], dtype=str),
        feature_names=np.asarray([item[1] for item in records], dtype=str),
        medians=np.asarray([item[4] for item in records], dtype=float),
        means=np.asarray([item[5] for item in records], dtype=float),
        scales=np.asarray([item[6] for item in records], dtype=float),
        training_sample_count=len(fit_indices),
    )


def fit_classifier(model: str, x: np.ndarray, y: np.ndarray, c_value: float, l1_ratio: float | None, seed: int):
    if model == "elastic_net":
        fitted = LogisticRegression(
            penalty="elasticnet", solver="saga", C=float(c_value), l1_ratio=float(l1_ratio),
            max_iter=3000, tol=1e-3, random_state=seed,
        )
    elif model == "linear_svm":
        fitted = LinearSVC(C=float(c_value), dual="auto", max_iter=10000, random_state=seed)
    else:
        raise ValueError(f"Unknown model: {model}")
    return fitted.fit(x, y)


def raw_model_score(model: str, fitted, x: np.ndarray) -> np.ndarray:
    if model == "elastic_net":
        return fitted.predict_proba(x)[:, 1]
    return fitted.decision_function(x)


def fit_platt(scores: np.ndarray, labels: np.ndarray, seed: int) -> LogisticRegression:
    calibrator = LogisticRegression(C=1.0, solver="lbfgs", random_state=seed)
    return calibrator.fit(np.asarray(scores).reshape(-1, 1), labels)


def probability_from_score(model: str, scores: np.ndarray, calibrator: LogisticRegression | None) -> np.ndarray:
    if model == "elastic_net":
        return np.clip(scores, 1e-8, 1 - 1e-8)
    if calibrator is None:
        raise ValueError("Linear SVM probabilities require a training-only calibrator")
    return np.clip(calibrator.predict_proba(np.asarray(scores).reshape(-1, 1))[:, 1], 1e-8, 1 - 1e-8)


def choose_threshold(probability: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(np.clip(probability, 1e-8, 1 - 1e-8))
    candidates = np.unique(np.concatenate(([0.5], unique, (unique[:-1] + unique[1:]) / 2)))
    rows = []
    for threshold in candidates:
        score = balanced_accuracy_score(labels, probability >= threshold)
        rows.append((float(score), -abs(float(threshold) - 0.5), float(threshold)))
    return max(rows)[2]


def classification_metrics(labels: np.ndarray, probability: np.ndarray, threshold: float) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    probability = np.clip(np.asarray(probability, dtype=float), 1e-8, 1 - 1e-8)
    predicted = probability >= threshold
    sensitivity = float(predicted[labels == 1].mean())
    specificity = float((~predicted[labels == 0]).mean())
    result = {
        "roc_auc": float(roc_auc_score(labels, probability)),
        "balanced_accuracy": float((sensitivity + specificity) / 2),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "brier_score": float(brier_score_loss(labels, probability)),
        "log_loss": float(log_loss(labels, probability, labels=[0, 1])),
    }
    logits = np.log(probability / (1 - probability)).reshape(-1, 1)
    try:
        calibration = LogisticRegression(C=1e6, solver="lbfgs").fit(logits, labels)
        result["calibration_intercept"] = float(calibration.intercept_[0])
        result["calibration_slope"] = float(calibration.coef_[0, 0])
    except Exception:
        result["calibration_intercept"] = np.nan
        result["calibration_slope"] = np.nan
    return result


def fit_pca_transform(
    log2_matrix: np.ndarray,
    labels: np.ndarray,
    fit_indices: np.ndarray,
    allowed_training_indices: np.ndarray,
    candidate_features: int,
    components: int,
    seed: int,
) -> tuple[dict, np.ndarray]:
    validate_fit_scope(fit_indices, allowed_training_indices)
    training = log2_matrix[np.asarray(fit_indices, dtype=int)]
    finite = np.isfinite(training)
    eligible = finite.mean(axis=0) >= 0.5
    indices = np.flatnonzero(eligible)
    medians = np.nanmedian(training[:, indices], axis=0)
    filled = np.where(np.isfinite(training[:, indices]), training[:, indices], medians)
    variance = filled.var(axis=0)
    order = np.argsort(variance, kind="stable")[-min(candidate_features, len(indices)):]
    indices, medians, filled = indices[order], medians[order], filled[:, order]
    means, scales = filled.mean(axis=0), filled.std(axis=0)
    nonconstant = scales > 0
    indices, medians, means, scales = indices[nonconstant], medians[nonconstant], means[nonconstant], scales[nonconstant]
    standardized = (filled[:, nonconstant] - means) / scales
    pca = PCA(n_components=min(components, len(fit_indices) - 1, standardized.shape[1]), random_state=seed).fit(standardized)
    state = {"indices": indices, "medians": medians, "means": means, "scales": scales, "pca": pca}
    return state, pca.transform(standardized)


def apply_pca_transform(state: dict, log2_matrix: np.ndarray, indices: np.ndarray) -> np.ndarray:
    values = log2_matrix[np.asarray(indices, dtype=int)][:, state["indices"]]
    filled = np.where(np.isfinite(values), values, state["medians"])
    return state["pca"].transform((filled - state["means"]) / state["scales"])


def parameter_grid(model: str, config: dict) -> list[dict]:
    if model == "elastic_net":
        return [
            {"C": float(c), "l1_ratio": float(ratio), "selector_k": int(config["elastic_selector_k"])}
            for c in config["c_grid"] for ratio in config["l1_ratio_grid"]
        ]
    return [
        {"C": float(c), "l1_ratio": None, "selector_k": int(k)}
        for c in config["c_grid"] for k in config["svm_selector_k_grid"]
    ]


def summarize_metric_distribution(metrics: pd.DataFrame, group_columns: Iterable[str]) -> pd.DataFrame:
    value_columns = [
        "roc_auc", "balanced_accuracy", "sensitivity", "specificity", "brier_score", "log_loss",
        "calibration_intercept", "calibration_slope",
    ]
    rows = []
    for keys, frame in metrics.groupby(list(group_columns), dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        base = dict(zip(group_columns, keys))
        for metric in value_columns:
            values = frame[metric].dropna().to_numpy(dtype=float)
            if len(values):
                rows.append({
                    **base, "metric": metric, "median": float(np.median(values)),
                    "interval_low": float(np.quantile(values, 0.025)),
                    "interval_high": float(np.quantile(values, 0.975)), "repeat_count": len(values),
                })
    return pd.DataFrame(rows)
