"""Patient-change PCA, bootstrap stability, and consensus clustering."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PCAResult:
    scores: np.ndarray
    loadings: np.ndarray
    explained_fraction: np.ndarray
    selected_feature_indices: np.ndarray
    scaled_matrix: np.ndarray


def robust_change_pca(delta: np.ndarray, top_features: int, winsor_limit: float) -> PCAResult:
    complete = np.isfinite(delta).all(axis=0)
    complete_indices = np.flatnonzero(complete)
    values = delta[:, complete_indices]
    medians = np.median(values, axis=0)
    mad = np.median(np.abs(values - medians), axis=0) * 1.4826
    variable = mad > 0
    complete_indices = complete_indices[variable]
    values, mad = values[:, variable], mad[variable]
    selection = np.argsort(mad, kind="stable")[-min(top_features, len(mad)):]
    selected_indices = complete_indices[selection]
    scaled = (values[:, selection] - np.median(values[:, selection], axis=0)) / mad[selection]
    scaled = np.clip(scaled, -winsor_limit, winsor_limit)
    scaled -= scaled.mean(axis=0)
    u, singular, vt = np.linalg.svd(scaled, full_matrices=False)
    explained = singular**2 / np.sum(singular**2)
    return PCAResult(u * singular, vt.T, explained, selected_indices, scaled)


def bootstrap_explained_variance(scaled: np.ndarray, resamples: int, seed: int, components: int = 10) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    values = np.empty((resamples, components))
    for bootstrap in range(resamples):
        sampled = scaled[rng.integers(0, scaled.shape[0], size=scaled.shape[0])]
        sampled = sampled - sampled.mean(axis=0)
        singular = np.linalg.svd(sampled, full_matrices=False, compute_uv=False)
        fraction = singular**2 / np.sum(singular**2)
        values[bootstrap] = fraction[:components]
    rows = []
    for component in range(components):
        rows.append({
            "component": component + 1,
            "observed_fraction": np.nan,
            "bootstrap_median": float(np.median(values[:, component])),
            "bootstrap_ci_low": float(np.quantile(values[:, component], 0.025)),
            "bootstrap_ci_high": float(np.quantile(values[:, component], 0.975)),
        })
    return pd.DataFrame(rows)


def _kmeans_once(matrix: np.ndarray, k: int, rng: np.random.Generator, iterations: int = 100) -> tuple[np.ndarray, float]:
    n = len(matrix)
    centers = [matrix[rng.integers(0, n)]]
    while len(centers) < k:
        distances = np.min(np.column_stack([np.sum((matrix - center) ** 2, axis=1) for center in centers]), axis=1)
        probability = distances / distances.sum() if distances.sum() > 0 else np.full(n, 1 / n)
        centers.append(matrix[rng.choice(n, p=probability)])
    centers_array = np.asarray(centers)
    labels = np.zeros(n, dtype=int)
    for _ in range(iterations):
        distance = np.column_stack([np.sum((matrix - center) ** 2, axis=1) for center in centers_array])
        updated = np.argmin(distance, axis=1)
        if np.array_equal(updated, labels):
            break
        labels = updated
        for cluster in range(k):
            if np.any(labels == cluster):
                centers_array[cluster] = matrix[labels == cluster].mean(axis=0)
            else:
                centers_array[cluster] = matrix[rng.integers(0, n)]
    inertia = float(sum(np.sum((matrix[labels == cluster] - centers_array[cluster]) ** 2) for cluster in range(k)))
    return labels, inertia


def kmeans(matrix: np.ndarray, k: int, seed: int, starts: int = 20) -> np.ndarray:
    rng = np.random.default_rng(seed)
    best_labels, best_inertia = None, np.inf
    for _ in range(starts):
        labels, inertia = _kmeans_once(matrix, k, rng)
        if inertia < best_inertia:
            best_labels, best_inertia = labels.copy(), inertia
    assert best_labels is not None
    return best_labels


def silhouette_score(matrix: np.ndarray, labels: np.ndarray) -> float:
    distances = np.sqrt(np.maximum(0.0, np.sum((matrix[:, None, :] - matrix[None, :, :]) ** 2, axis=2)))
    scores = []
    for index in range(len(matrix)):
        same = labels == labels[index]
        same[index] = False
        a = float(distances[index, same].mean()) if same.any() else 0.0
        b = min(float(distances[index, labels == cluster].mean()) for cluster in np.unique(labels) if cluster != labels[index])
        scores.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(scores))


def consensus_clustering(
    matrix: np.ndarray,
    k_values: list[int],
    resamples: int,
    patient_fraction: float,
    pac_interval: tuple[float, float],
    thresholds: dict,
    seed: int,
) -> tuple[pd.DataFrame, dict[int, np.ndarray], dict[int, np.ndarray]]:
    rng = np.random.default_rng(seed)
    n = len(matrix)
    subset_size = max(2, int(round(patient_fraction * n)))
    summaries, consensus_by_k, labels_by_k = [], {}, {}
    for k in k_values:
        coobserved = np.zeros((n, n), dtype=float)
        coclustered = np.zeros((n, n), dtype=float)
        for iteration in range(resamples):
            selected = np.sort(rng.choice(n, size=subset_size, replace=False))
            labels = kmeans(matrix[selected], k, seed + k * 100000 + iteration, starts=5)
            coobserved[np.ix_(selected, selected)] += 1
            for cluster in range(k):
                members = selected[labels == cluster]
                coclustered[np.ix_(members, members)] += 1
        consensus = np.divide(coclustered, coobserved, out=np.zeros_like(coclustered), where=coobserved > 0)
        np.fill_diagonal(consensus, 1.0)
        base_labels = kmeans(matrix, k, seed + k, starts=50)
        off_diagonal = consensus[np.triu_indices(n, 1)]
        pac = float(((off_diagonal > pac_interval[0]) & (off_diagonal < pac_interval[1])).mean())
        within_mask = (base_labels[:, None] == base_labels[None, :]) & ~np.eye(n, dtype=bool)
        within = float(consensus[within_mask].mean())
        silhouette = silhouette_score(matrix, base_labels)
        minimum_size = int(pd.Series(base_labels).value_counts().min())
        accepted = (
            pac <= thresholds["maximum_pac"]
            and within >= thresholds["minimum_within_consensus"]
            and silhouette >= thresholds["minimum_silhouette"]
            and minimum_size >= thresholds["minimum_cluster_size"]
        )
        summaries.append({
            "k": k,
            "pac": pac,
            "mean_within_cluster_consensus": within,
            "silhouette": silhouette,
            "minimum_cluster_size": minimum_size,
            "accepted": accepted,
        })
        consensus_by_k[k] = consensus
        labels_by_k[k] = base_labels
    return pd.DataFrame(summaries), consensus_by_k, labels_by_k

