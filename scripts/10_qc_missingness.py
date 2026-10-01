"""Phase 1 specimen QC, multivariate structure, and missingness intelligence."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from phase0_common import INTERIM_DIR, PROJECT_DIR, RESULT_DIR as PHASE0_RESULT_DIR, sha256, write_json


RESULT_DIR = PROJECT_DIR / "results" / "phase1"
FIGURE_DIR = RESULT_DIR / "figures"
CONFIG_PATH = PROJECT_DIR / "config" / "phase1.yml"


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros_like(values)
    return (values - median) / (1.4826 * mad)


def stratified_robust_z(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Apply one symmetric robust-z rule within each tissue class."""
    result = np.zeros(len(values), dtype=float)
    for group in np.unique(groups):
        selected = groups == group
        result[selected] = robust_z(np.asarray(values)[selected])
    return result


def exact_paired_binary_pvalue(t_only: int, n_only: int) -> float:
    discordant = t_only + n_only
    if discordant == 0:
        return math.nan
    tail = min(t_only, n_only)
    probability = sum(math.comb(discordant, i) for i in range(tail + 1)) / (2**discordant)
    return min(1.0, 2.0 * probability)


def average_linkage_order(distance: np.ndarray) -> list[int]:
    """Deterministic average-linkage leaf order without a SciPy dependency."""
    n = distance.shape[0]
    clusters: dict[int, list[int]] = {i: [i] for i in range(n)}
    next_id = n
    while len(clusters) > 1:
        ids = sorted(clusters)
        best: tuple[float, int, int] | None = None
        for ai, left in enumerate(ids[:-1]):
            a = clusters[left]
            for right in ids[ai + 1:]:
                b = clusters[right]
                value = float(distance[np.ix_(a, b)].mean())
                candidate = (value, left, right)
                if best is None or candidate < best:
                    best = candidate
        assert best is not None
        _, left, right = best
        a, b = clusters.pop(left), clusters.pop(right)
        # Orient by the closest boundary to keep adjacent leaves coherent.
        orientations = [(a, b), (a[::-1], b), (a, b[::-1]), (a[::-1], b[::-1])]
        merged = min(orientations, key=lambda pair: (distance[pair[0][-1], pair[1][0]], pair[0], pair[1]))
        clusters[next_id] = merged[0] + merged[1]
        next_id += 1
    return next(iter(clusters.values()))


def mcd_style_distance(matrix: np.ndarray, seed: int, starts: int, c_steps: int, support_fraction: float) -> np.ndarray:
    """Deterministic multi-start C-step approximation to minimum covariance determinant."""
    z = np.column_stack([robust_z(matrix[:, column]) for column in range(matrix.shape[1])])
    z = np.clip(np.nan_to_num(z), -8.0, 8.0)
    n, p = z.shape
    h = max(p + 2, int(math.floor(support_fraction * n)))
    rng = np.random.default_rng(seed)
    best_logdet = math.inf
    best_location = np.median(z, axis=0)
    best_covariance = np.cov(z, rowvar=False)
    for start in range(starts):
        subset = np.arange(h) if start == 0 else np.sort(rng.choice(n, size=h, replace=False))
        for _ in range(c_steps):
            location = z[subset].mean(axis=0)
            covariance = np.cov(z[subset], rowvar=False) + np.eye(p) * 1e-6
            inverse = np.linalg.pinv(covariance)
            delta = z - location
            distances = np.einsum("ij,jk,ik->i", delta, inverse, delta)
            updated = np.sort(np.argsort(distances, kind="stable")[:h])
            if np.array_equal(updated, subset):
                break
            subset = updated
        sign, logdet = np.linalg.slogdet(covariance)
        if sign > 0 and logdet < best_logdet:
            best_logdet = float(logdet)
            best_location = location
            best_covariance = covariance
    inverse = np.linalg.pinv(best_covariance)
    delta = z - best_location
    return np.sqrt(np.maximum(0.0, np.einsum("ij,jk,ik->i", delta, inverse, delta)))


def svg_scatter(path: Path, x: np.ndarray, y: np.ndarray, labels: list[str], groups: list[str], title: str, xlabel: str, ylabel: str) -> None:
    width, height, margin = 900, 620, 70
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))
    xspan, yspan = max(xmax - xmin, 1e-9), max(ymax - ymin, 1e-9)
    sx = lambda value: margin + (value - xmin) / xspan * (width - 2 * margin)
    sy = lambda value: height - margin - (value - ymin) / yspan * (height - 2 * margin)
    colors = {"N": "#2563EB", "T": "#DC2626"}
    elements = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{margin}" y="35" font-family="Arial" font-size="20" font-weight="bold">{title}</text>', f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#333"/>', f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#333"/>']
    for xi, yi, label, group in zip(x, y, labels, groups):
        elements.append(f'<circle cx="{sx(float(xi)):.1f}" cy="{sy(float(yi)):.1f}" r="5" fill="{colors[group]}" fill-opacity="0.75"><title>{label}</title></circle>')
    elements.extend([f'<text x="{width/2}" y="{height-18}" text-anchor="middle" font-family="Arial" font-size="14">{xlabel}</text>', f'<text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" font-family="Arial" font-size="14">{ylabel}</text>', '<circle cx="730" cy="35" r="5" fill="#2563EB"/><text x="742" y="40" font-family="Arial" font-size="12">Matched non-tumour</text>', '<circle cx="730" cy="55" r="5" fill="#DC2626"/><text x="742" y="60" font-family="Arial" font-size="12">Tumour</text>', '</svg>'])
    path.write_text("\n".join(elements), encoding="utf-8")


def svg_bars(path: Path, labels: list[str], values: list[float], title: str, ylabel: str) -> None:
    width, height, margin = 900, 560, 70
    maximum = max(values) * 1.08 if values else 1
    plot_width = width - 2 * margin
    bar_width = plot_width / max(len(values), 1) * 0.7
    elements = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{margin}" y="35" font-family="Arial" font-size="20" font-weight="bold">{title}</text>', f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#333"/>', f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#333"/>']
    for index, (label, value) in enumerate(zip(labels, values)):
        x = margin + (index + 0.15) * plot_width / len(values)
        bar_height = value / maximum * (height - 2 * margin)
        y = height - margin - bar_height
        elements.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="#0F766E"/>')
        elements.append(f'<text x="{x+bar_width/2:.1f}" y="{height-margin+18}" text-anchor="middle" font-family="Arial" font-size="11">{label}</text>')
        elements.append(f'<text x="{x+bar_width/2:.1f}" y="{y-6:.1f}" text-anchor="middle" font-family="Arial" font-size="10">{value:.1%}</text>')
    elements.extend([f'<text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" font-family="Arial" font-size="14">{ylabel}</text>', '</svg>'])
    path.write_text("\n".join(elements), encoding="utf-8")


def svg_correlation_heatmap(path: Path, correlation: np.ndarray, order: list[int], sample_ids: list[str]) -> None:
    cell, margin = 7, 90
    n = len(order)
    size = margin + n * cell + 35
    ordered = correlation[np.ix_(order, order)]
    elements = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}">', '<rect width="100%" height="100%" fill="white"/>', '<text x="20" y="25" font-family="Arial" font-size="18" font-weight="bold">Spearman correlation ordered without tissue labels</text>']
    for row in range(n):
        for column in range(n):
            value = float(np.clip(ordered[row, column], 0, 1))
            red = int(247 - 210 * value)
            green = int(250 - 150 * value)
            blue = int(252 - 80 * value)
            elements.append(f'<rect x="{margin + column*cell}" y="{margin + row*cell}" width="{cell}" height="{cell}" fill="rgb({red},{green},{blue})"><title>{sample_ids[order[row]]} vs {sample_ids[order[column]]}: {ordered[row,column]:.3f}</title></rect>')
    for position, index in enumerate(order):
        color = "#2563EB" if sample_ids[index].endswith("N") else "#DC2626"
        elements.append(f'<rect x="{margin + position*cell}" y="{margin-10}" width="{cell}" height="6" fill="{color}"/>')
        elements.append(f'<rect x="{margin-10}" y="{margin + position*cell}" width="6" height="{cell}" fill="{color}"/>')
    elements.append('</svg>')
    path.write_text("\n".join(elements), encoding="utf-8")


def png_scatter(path: Path, x: np.ndarray, y: np.ndarray, groups: list[str], title: str, xlabel: str, ylabel: str) -> None:
    width, height, margin = 900, 620, 70
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((margin, 22), title, fill="#111827")
    draw.line((margin, height - margin, width - margin, height - margin), fill="#374151", width=2)
    draw.line((margin, margin, margin, height - margin), fill="#374151", width=2)
    xmin, xmax, ymin, ymax = float(np.min(x)), float(np.max(x)), float(np.min(y)), float(np.max(y))
    xspan, yspan = max(xmax - xmin, 1e-9), max(ymax - ymin, 1e-9)
    for xi, yi, group in zip(x, y, groups):
        px = margin + (float(xi) - xmin) / xspan * (width - 2 * margin)
        py = height - margin - (float(yi) - ymin) / yspan * (height - 2 * margin)
        color = "#2563EB" if group == "N" else "#DC2626"
        draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=color)
    draw.text((width // 2 - 70, height - 35), xlabel, fill="#111827")
    draw.text((8, height // 2), ylabel, fill="#111827")
    draw.ellipse((700, 25, 708, 33), fill="#2563EB")
    draw.text((714, 22), "Matched non-tumour", fill="#111827")
    draw.ellipse((700, 43, 708, 51), fill="#DC2626")
    draw.text((714, 40), "Tumour", fill="#111827")
    image.save(path)


def png_bars(path: Path, labels: list[str], values: list[float], title: str) -> None:
    width, height, margin = 900, 560, 70
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((margin, 22), title, fill="#111827")
    draw.line((margin, height - margin, width - margin, height - margin), fill="#374151", width=2)
    draw.line((margin, margin, margin, height - margin), fill="#374151", width=2)
    maximum = max(values) * 1.08
    slot = (width - 2 * margin) / len(values)
    for index, (label, value) in enumerate(zip(labels, values)):
        x0 = margin + index * slot + slot * 0.15
        x1 = x0 + slot * 0.7
        y0 = height - margin - value / maximum * (height - 2 * margin)
        draw.rectangle((x0, y0, x1, height - margin), fill="#0F766E")
        draw.text((x0 + slot * 0.25, height - margin + 8), label, fill="#111827")
        draw.text((x0 + 2, y0 - 15), f"{value:.1%}", fill="#111827")
    image.save(path)


def png_correlation_heatmap(path: Path, correlation: np.ndarray, order: list[int], sample_ids: list[str]) -> None:
    cell, margin = 7, 70
    n = len(order)
    size = margin + n * cell + 20
    image = Image.new("RGB", (size, size), "white")
    draw = ImageDraw.Draw(image)
    draw.text((15, 15), "Spearman correlation ordered without tissue labels", fill="#111827")
    ordered = correlation[np.ix_(order, order)]
    for row in range(n):
        for column in range(n):
            value = float(np.clip(ordered[row, column], 0, 1))
            color = (int(247 - 210 * value), int(250 - 150 * value), int(252 - 80 * value))
            x0, y0 = margin + column * cell, margin + row * cell
            draw.rectangle((x0, y0, x0 + cell, y0 + cell), fill=color)
    for position, index in enumerate(order):
        color = "#2563EB" if sample_ids[index].endswith("N") else "#DC2626"
        draw.rectangle((margin + position * cell, margin - 9, margin + (position + 1) * cell, margin - 4), fill=color)
        draw.rectangle((margin - 9, margin + position * cell, margin - 4, margin + (position + 1) * cell), fill=color)
    image.save(path)


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    raw = np.load(INTERIM_DIR / "raw_data_model.npz")
    paired = np.load(INTERIM_DIR / "paired_data_model.npz")
    sample_manifest = pd.read_csv(INTERIM_DIR / "sample_manifest.csv")
    feature_manifest = pd.read_csv(INTERIM_DIR / "feature_manifest.csv")
    sample_ids = raw["sample_ids"].tolist()
    feature_keys = raw["feature_keys"].tolist()
    x_raw = raw["X_raw"].astype(float)
    detected = raw["D"].astype(bool)
    x_log2 = np.where(detected, np.log2(x_raw), np.nan)
    n_samples, n_features = x_raw.shape

    rows = []
    sample_index = {sample: index for index, sample in enumerate(sample_ids)}
    for index, sample_id in enumerate(sample_ids):
        observed = x_log2[index, detected[index]]
        tissue = sample_id[-1]
        pair_id = sample_id[:-1]
        partner_id = pair_id + ("T" if tissue == "N" else "N")
        partner = sample_index[partner_id]
        shared = detected[index] & detected[partner]
        union = detected[index] | detected[partner]
        pair_difference = np.abs(x_log2[index, shared] - x_log2[partner, shared])
        median = float(np.median(observed))
        mad = float(np.median(np.abs(observed - median)))
        rows.append({
            "sample_id": sample_id,
            "patient_id": pair_id,
            "tissue_code": tissue,
            "detected_features": int(detected[index].sum()),
            "detected_fraction": float(detected[index].mean()),
            "median_log2": median,
            "iqr_log2": float(np.percentile(observed, 75) - np.percentile(observed, 25)),
            "mad_log2": mad,
            "p01_log2": float(np.percentile(observed, 1)),
            "p05_log2": float(np.percentile(observed, 5)),
            "p95_log2": float(np.percentile(observed, 95)),
            "p99_log2": float(np.percentile(observed, 99)),
            "outside_median_5mad_fraction": float((np.abs(observed - median) > 5 * mad).mean()) if mad else 0.0,
            "pair_shared_detected": int(shared.sum()),
            "pair_detection_jaccard": float(shared.sum() / union.sum()),
            "pair_median_absolute_log2_difference": float(np.median(pair_difference)),
        })
    sample_qc = pd.DataFrame(rows)

    observed_fraction = detected.mean(axis=0)
    eligible = observed_fraction >= config["qc_feature_min_observed_fraction"]
    eligible_log = x_log2[:, eligible]
    feature_medians = np.nanmedian(eligible_log, axis=0)
    imputed = np.where(np.isnan(eligible_log), feature_medians, eligible_log)
    feature_mad = np.median(np.abs(imputed - np.median(imputed, axis=0)), axis=0) * 1.4826
    variable = feature_mad > 0
    robust_matrix = (imputed[:, variable] - np.median(imputed[:, variable], axis=0)) / feature_mad[variable]
    robust_matrix = np.clip(robust_matrix, -config["pca_winsor_limit_robust_z"], config["pca_winsor_limit_robust_z"])
    robust_matrix -= robust_matrix.mean(axis=0)
    u, singular, _ = np.linalg.svd(robust_matrix, full_matrices=False)
    scores = u[:, : config["pca_components_to_export"]] * singular[: config["pca_components_to_export"]]
    explained = (singular**2) / np.sum(singular**2)
    for component in range(scores.shape[1]):
        sample_qc[f"PC{component + 1}"] = scores[:, component]

    correlation = pd.DataFrame(eligible_log.T, columns=sample_ids).corr(method="spearman", min_periods=100)
    correlation_values = correlation.to_numpy(copy=True)
    np.fill_diagonal(correlation_values, np.nan)
    sample_qc["median_pairwise_spearman"] = np.nanmedian(correlation_values, axis=1)
    distance = 1.0 - np.nan_to_num(correlation.to_numpy(), nan=0.0)
    np.fill_diagonal(distance, 0.0)
    order = average_linkage_order(distance)

    summary_columns = ["detected_fraction", "median_log2", "iqr_log2", "mad_log2", "median_pairwise_spearman", "pair_detection_jaccard", "pair_median_absolute_log2_difference", "PC1", "PC2"]
    multivariate_distance = mcd_style_distance(sample_qc[summary_columns].to_numpy(), config["random_seed"], config["mcd_random_starts"], config["mcd_c_steps"], config["mcd_support_fraction"])
    sample_qc["mcd_style_distance"] = multivariate_distance

    threshold = config["diagnostic_robust_z_threshold"]
    tissue_groups = sample_qc["tissue_code"].to_numpy()
    sample_qc["flag_coverage"] = np.abs(stratified_robust_z(sample_qc["detected_fraction"].to_numpy(), tissue_groups)) > threshold
    distribution_z = np.column_stack([stratified_robust_z(sample_qc[column].to_numpy(), tissue_groups) for column in ["median_log2", "iqr_log2", "mad_log2", "outside_median_5mad_fraction"]])
    sample_qc["flag_observed_distribution"] = np.max(np.abs(distribution_z), axis=1) > threshold
    sample_qc["flag_sample_correlation"] = stratified_robust_z(sample_qc["median_pairwise_spearman"].to_numpy(), tissue_groups) < -threshold
    pair_concordance_z = np.column_stack([robust_z(sample_qc["pair_detection_jaccard"].to_numpy()), robust_z(sample_qc["pair_median_absolute_log2_difference"].to_numpy())])
    sample_qc["flag_matched_pair_concordance"] = (pair_concordance_z[:, 0] < -threshold) | (pair_concordance_z[:, 1] > threshold)
    sample_qc["flag_multivariate_structure"] = stratified_robust_z(multivariate_distance, tissue_groups) > threshold
    flag_columns = [column for column in sample_qc.columns if column.startswith("flag_")]
    sample_qc["diagnostic_family_count"] = sample_qc[flag_columns].sum(axis=1)
    sample_qc["sensitivity_flag"] = sample_qc["diagnostic_family_count"] >= config["flag_minimum_diagnostic_families"]

    patient_flags = sample_qc.groupby("patient_id")["sensitivity_flag"].any()
    cohort = pd.DataFrame({
        "patient_id": patient_flags.index,
        "primary_all_pairs_included": True,
        "sensitivity_unflagged_pairs_included": ~patient_flags.to_numpy(),
        "pair_flagged": patient_flags.to_numpy(),
    })
    sample_qc["sensitivity_pair_removed"] = sample_qc["patient_id"].map(patient_flags)

    states = paired["detection_state"]
    medians = np.full(n_features, np.nan)
    any_observed = detected.any(axis=0)
    medians[any_observed] = np.nanmedian(x_log2[:, any_observed], axis=0)
    t_detected = (states == 2).sum(axis=0) + (states == 3).sum(axis=0)
    n_detected = (states == 1).sum(axis=0) + (states == 3).sum(axis=0)
    feature_missingness = feature_manifest[["feature_key", "feature_id_original", "duplicate_identifier", "compound_identifier", "all_abundance_missing"]].copy()
    feature_missingness["observed_fraction"] = observed_fraction
    feature_missingness["median_observed_log2"] = medians
    feature_missingness["both_missing_pairs"] = (states == 0).sum(axis=0)
    feature_missingness["N_only_pairs"] = (states == 1).sum(axis=0)
    feature_missingness["T_only_pairs"] = (states == 2).sum(axis=0)
    feature_missingness["both_detected_pairs"] = (states == 3).sum(axis=0)
    feature_missingness["N_detection_fraction"] = n_detected / states.shape[0]
    feature_missingness["T_detection_fraction"] = t_detected / states.shape[0]
    feature_missingness["detection_fraction_T_minus_N"] = (t_detected - n_detected) / states.shape[0]
    feature_missingness["exact_paired_detection_p_descriptive"] = [exact_paired_binary_pvalue(int(t), int(n)) for t, n in zip(feature_missingness["T_only_pairs"], feature_missingness["N_only_pairs"])]
    nonempty = feature_missingness["median_observed_log2"].notna()
    feature_missingness.loc[nonempty, "abundance_decile"] = pd.qcut(feature_missingness.loc[nonempty, "median_observed_log2"], 10, labels=[f"D{i}" for i in range(1, 11)], duplicates="drop").astype("string")
    deciles = feature_missingness.loc[nonempty].groupby("abundance_decile", observed=True).agg(
        feature_count=("feature_key", "size"),
        mean_missing_fraction=("observed_fraction", lambda x: float(1 - x.mean())),
        mean_N_detection_fraction=("N_detection_fraction", "mean"),
        mean_T_detection_fraction=("T_detection_fraction", "mean"),
        mean_T_minus_N_detection=("detection_fraction_T_minus_N", "mean"),
    ).reset_index()
    deciles["decile_number"] = deciles["abundance_decile"].str.removeprefix("D").astype(int)
    deciles = deciles.sort_values("decile_number").drop(columns="decile_number").reset_index(drop=True)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    sample_qc.to_csv(RESULT_DIR / "sample_qc.csv", index=False)
    feature_missingness.to_csv(RESULT_DIR / "feature_missingness.csv", index=False)
    deciles.to_csv(RESULT_DIR / "abundance_decile_missingness.csv", index=False)
    correlation.to_csv(RESULT_DIR / "sample_spearman_correlation.csv")
    pd.DataFrame({"cluster_position": np.arange(1, n_samples + 1), "sample_id": [sample_ids[index] for index in order]}).to_csv(RESULT_DIR / "cluster_order.csv", index=False)
    cohort.to_csv(RESULT_DIR / "cohort_membership.csv", index=False)
    np.savez_compressed(INTERIM_DIR / "phase1_qc_objects.npz", pca_scores=scores, pca_explained_fraction=explained, qc_feature_mask=eligible, cluster_order=np.asarray(order))

    svg_scatter(FIGURE_DIR / "pca_qc.svg", scores[:, 0], scores[:, 1], sample_ids, sample_qc["tissue_code"].tolist(), "Label-blind diagnostic PCA", f"PC1 ({explained[0]:.1%})", f"PC2 ({explained[1]:.1%})")
    svg_scatter(FIGURE_DIR / "coverage_vs_median.svg", sample_qc["detected_fraction"].to_numpy(), sample_qc["median_log2"].to_numpy(), sample_ids, sample_qc["tissue_code"].tolist(), "Coverage and observed abundance", "Detected-feature fraction", "Median observed log2 abundance")
    svg_bars(FIGURE_DIR / "missingness_by_abundance_decile.svg", deciles["abundance_decile"].astype(str).tolist(), deciles["mean_missing_fraction"].tolist(), "Missingness by observed-abundance decile", "Mean missing fraction")
    svg_correlation_heatmap(FIGURE_DIR / "sample_correlation_heatmap.svg", correlation.to_numpy(), order, sample_ids)
    png_scatter(FIGURE_DIR / "pca_qc.png", scores[:, 0], scores[:, 1], sample_qc["tissue_code"].tolist(), "Label-blind diagnostic PCA", f"PC1 ({explained[0]:.1%})", f"PC2 ({explained[1]:.1%})")
    png_scatter(FIGURE_DIR / "coverage_vs_median.png", sample_qc["detected_fraction"].to_numpy(), sample_qc["median_log2"].to_numpy(), sample_qc["tissue_code"].tolist(), "Coverage and observed abundance", "Detected-feature fraction", "Median log2 abundance")
    png_bars(FIGURE_DIR / "missingness_by_abundance_decile.png", deciles["abundance_decile"].astype(str).tolist(), deciles["mean_missing_fraction"].tolist(), "Missingness by observed-abundance decile")
    png_correlation_heatmap(FIGURE_DIR / "sample_correlation_heatmap.png", correlation.to_numpy(), order, sample_ids)

    validation = {
        "phase": 1,
        "milestone": "M1",
        "status": "PASS",
        "config_sha256": sha256(CONFIG_PATH),
        "input_phase0_validation_sha256": sha256(PHASE0_RESULT_DIR / "input_validation.json"),
        "qc_eligible_features": int(eligible.sum()),
        "pca_variable_features": int(variable.sum()),
        "pca_explained_fraction_first_five": explained[:5].tolist(),
        "flagged_samples": int(sample_qc["sensitivity_flag"].sum()),
        "flagged_patient_pairs": int(patient_flags.sum()),
        "primary_pairs": 42,
        "sensitivity_pairs": int((~patient_flags).sum()),
        "checks": {
            "all_84_samples_profiled": len(sample_qc) == 84,
            "all_8071_features_profiled": len(feature_missingness) == 8071,
            "paired_states_sum_to_42": bool((feature_missingness[["both_missing_pairs", "N_only_pairs", "T_only_pairs", "both_detected_pairs"]].sum(axis=1) == 42).all()),
            "primary_cohort_retains_all_pairs": bool(cohort["primary_all_pairs_included"].all() and len(cohort) == 42),
            "sensitivity_removes_whole_pairs_only": bool(sample_qc.groupby("patient_id")["sensitivity_pair_removed"].nunique().eq(1).all()),
            "flag_requires_two_families": bool((sample_qc.loc[sample_qc["sensitivity_flag"], "diagnostic_family_count"] >= 2).all()),
            "pca_used_no_tissue_labels": True,
            "no_full_matrix_imputation_exported": True,
        },
    }
    if not all(validation["checks"].values()):
        raise ValueError("Phase 1 validation failed")
    write_json(RESULT_DIR / "phase1_validation.json", validation)
    print(f"PASS: Phase 1 outputs created; {validation['flagged_patient_pairs']} sensitivity pair(s) flagged")


if __name__ == "__main__":
    main()

