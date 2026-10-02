"""Validation-readiness methods that do not masquerade as external validation."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.stats import nct, norm, t


def nondominated_fronts(frame: pd.DataFrame, objectives: list[str]) -> np.ndarray:
    """Assign one-based Pareto fronts, maximizing every named objective."""
    values = frame[objectives].to_numpy(dtype=float)
    values = np.where(np.isfinite(values), values, -np.inf)
    remaining = list(range(len(frame)))
    fronts = np.zeros(len(frame), dtype=int)
    front_number = 1
    while remaining:
        current: list[int] = []
        for candidate in remaining:
            dominated = False
            for other in remaining:
                if other == candidate:
                    continue
                if np.all(values[other] >= values[candidate]) and np.any(values[other] > values[candidate]):
                    dominated = True
                    break
            if not dominated:
                current.append(candidate)
        if not current:
            raise RuntimeError("Pareto sorting did not make progress")
        fronts[current] = front_number
        selected = set(current)
        remaining = [index for index in remaining if index not in selected]
        front_number += 1
    return fronts


def correlation_components(
    paired_changes: pd.DataFrame,
    absolute_threshold: float,
    minimum_complete_pairs: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build deterministic connected components from high absolute Spearman correlations."""
    correlations = paired_changes.corr(method="spearman", min_periods=minimum_complete_pairs)
    keys = list(correlations.columns)
    adjacency = {key: set() for key in keys}
    edge_rows: list[dict[str, object]] = []
    for left_index, left in enumerate(keys):
        for right in keys[left_index + 1:]:
            value = correlations.loc[left, right]
            if np.isfinite(value) and abs(float(value)) >= absolute_threshold:
                adjacency[left].add(right)
                adjacency[right].add(left)
                edge_rows.append({"feature_key_1": left, "feature_key_2": right, "spearman_rho": float(value)})
    component_rows: list[dict[str, object]] = []
    visited: set[str] = set()
    component_number = 0
    for start in sorted(keys):
        if start in visited:
            continue
        component_number += 1
        stack, members = [start], []
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            members.append(node)
            stack.extend(sorted(adjacency[node] - visited, reverse=True))
        component_id = f"AC{component_number:03d}"
        for member in sorted(members):
            component_rows.append({
                "feature_key": member,
                "abundance_correlation_component": component_id,
                "component_size": len(members),
            })
    return pd.DataFrame(component_rows), pd.DataFrame(
        edge_rows, columns=["feature_key_1", "feature_key_2", "spearman_rho"]
    )


def paired_t_required_n(effect_size: float, alpha: float, power: float, maximum_n: int = 5000) -> int:
    """Smallest paired-sample count attaining two-sided noncentral-t power."""
    if effect_size <= 0 or not 0 < alpha < 1 or not 0 < power < 1:
        raise ValueError("effect_size, alpha and power must be in their valid ranges")
    for sample_size in range(3, maximum_n + 1):
        critical = t.ppf(1 - alpha / 2, sample_size - 1)
        noncentrality = effect_size * math.sqrt(sample_size)
        achieved = nct.cdf(-critical, sample_size - 1, noncentrality) + 1 - nct.cdf(
            critical, sample_size - 1, noncentrality
        )
        if achieved >= power:
            return sample_size
    raise ValueError("maximum_n is too small for the requested design")


def mcnemar_approximate_required_n(difference: float, discordant_fraction: float, alpha: float, power: float) -> int:
    """Normal-approximation planning value for a paired binary endpoint."""
    if not 0 < difference <= discordant_fraction <= 1:
        raise ValueError("Require 0 < difference <= discordant_fraction <= 1")
    z_alpha = norm.ppf(1 - alpha / 2)
    z_power = norm.ppf(power)
    numerator = (
        z_alpha * math.sqrt(discordant_fraction)
        + z_power * math.sqrt(max(discordant_fraction - difference**2, 0.0))
    ) ** 2
    return int(math.ceil(numerator / difference**2))

