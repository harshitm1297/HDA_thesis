"""Dependency-light statistical functions used by paired inference."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    maximum_iterations, epsilon, floor = 300, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = floor if abs(d) < floor else d
    d = 1.0 / d
    h = d
    for iteration in range(1, maximum_iterations + 1):
        m2 = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = floor if abs(d) < floor else d
        c = 1.0 + aa / c
        c = floor if abs(c) < floor else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = floor if abs(d) < floor else d
        c = 1.0 + aa / c
        c = floor if abs(c) < floor else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            return h
    raise ArithmeticError("Incomplete beta continued fraction did not converge")


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(a, b, x) / a
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def student_t_cdf(value: float, degrees_freedom: float) -> float:
    if not math.isfinite(value) or degrees_freedom <= 0:
        return math.nan
    x = degrees_freedom / (degrees_freedom + value * value)
    tail = 0.5 * regularized_incomplete_beta(degrees_freedom / 2.0, 0.5, x)
    return 1.0 - tail if value >= 0 else tail


def student_t_two_sided_p(value: float, degrees_freedom: float) -> float:
    if not math.isfinite(value) or degrees_freedom <= 0:
        return math.nan
    x = degrees_freedom / (degrees_freedom + value * value)
    return regularized_incomplete_beta(degrees_freedom / 2.0, 0.5, x)


def student_t_ppf(probability: float, degrees_freedom: float) -> float:
    if probability == 0.5:
        return 0.0
    sign = 1.0 if probability > 0.5 else -1.0
    target = probability if probability > 0.5 else 1.0 - probability
    lower, upper = 0.0, 50.0
    for _ in range(100):
        middle = (lower + upper) / 2.0
        if student_t_cdf(middle, degrees_freedom) < target:
            lower = middle
        else:
            upper = middle
    return sign * (lower + upper) / 2.0


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    indices = np.flatnonzero(valid)
    if not len(indices):
        return adjusted
    order = indices[np.argsort(values[indices], kind="stable")]
    ranked = values[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def digamma(x: float) -> float:
    result = 0.0
    while x < 8.0:
        result -= 1.0 / x
        x += 1.0
    inverse = 1.0 / x
    inverse2 = inverse * inverse
    return result + math.log(x) - 0.5 * inverse - inverse2 * (1 / 12 - inverse2 * (1 / 120 - inverse2 / 252))


def trigamma(x: float) -> float:
    result = 0.0
    while x < 8.0:
        result += 1.0 / (x * x)
        x += 1.0
    inverse = 1.0 / x
    inverse2 = inverse * inverse
    return result + inverse + inverse2 / 2 + inverse2 * inverse / 6 - inverse2 * inverse2 * inverse / 30 + inverse2**3 * inverse / 42


def inverse_trigamma(target: float) -> float:
    lower, upper = 1e-4, 1e6
    for _ in range(100):
        middle = math.sqrt(lower * upper)
        if trigamma(middle) > target:
            lower = middle
        else:
            upper = middle
    return math.sqrt(lower * upper)


@dataclass(frozen=True)
class VariancePrior:
    degrees_freedom: float
    scale: float
    feature_count: int


def estimate_variance_prior(variances: np.ndarray, residual_df: np.ndarray) -> VariancePrior:
    valid = np.isfinite(variances) & (variances > 0) & np.isfinite(residual_df) & (residual_df > 0)
    s2, df = variances[valid], residual_df[valid]
    if len(s2) < 20:
        raise ValueError("At least 20 eligible variances are required for moderation")
    log_s2 = np.log(s2)
    lower, upper = np.quantile(log_s2, [0.05, 0.95])
    winsor = np.clip(log_s2, lower, upper)
    observed_variance = float(np.var(winsor, ddof=1))
    sampling_variance = float(np.mean([trigamma(value / 2.0) for value in df]))
    between_variance = max(observed_variance - sampling_variance, 1e-6)
    prior_df = min(2.0 * inverse_trigamma(between_variance), 1e6)
    adjusted_log_variance = np.mean([
        log_value - digamma(degree / 2.0) + math.log(degree / 2.0)
        for log_value, degree in zip(log_s2, df)
    ])
    log_scale = adjusted_log_variance - math.log(prior_df / 2.0) + digamma(prior_df / 2.0)
    return VariancePrior(prior_df, math.exp(log_scale), int(len(s2)))


def exact_paired_binary_pvalue(t_only: int, n_only: int) -> float:
    discordant = t_only + n_only
    if discordant == 0:
        return math.nan
    tail = min(t_only, n_only)
    probability = sum(math.comb(discordant, index) for index in range(tail + 1)) / (2**discordant)
    return min(1.0, 2.0 * probability)

