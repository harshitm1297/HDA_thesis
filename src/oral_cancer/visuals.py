"""Small dependency-light scientific figures for Phase 3."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw


def _scatter(path: Path, x: np.ndarray, y: np.ndarray, highlighted: np.ndarray, title: str, xlabel: str, ylabel: str) -> None:
    width, height, margin = 920, 640, 75
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y, highlighted = x[finite], y[finite], highlighted[finite]
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))
    xpad, ypad = max((xmax - xmin) * 0.03, 1e-6), max((ymax - ymin) * 0.03, 1e-6)
    xmin, xmax, ymin, ymax = xmin - xpad, xmax + xpad, ymin - ypad, ymax + ypad
    sx = lambda value: margin + (float(value) - xmin) / (xmax - xmin) * (width - 2 * margin)
    sy = lambda value: height - margin - (float(value) - ymin) / (ymax - ymin) * (height - 2 * margin)
    draw.text((margin, 24), title, fill="#111827")
    draw.line((margin, height - margin, width - margin, height - margin), fill="#374151", width=2)
    draw.line((margin, margin, margin, height - margin), fill="#374151", width=2)
    if xmin <= 0 <= xmax:
        draw.line((sx(0), margin, sx(0), height - margin), fill="#D1D5DB", width=1)
    for xi, yi, selected in zip(x, y, highlighted):
        color = "#DC2626" if selected else "#94A3B8"
        radius = 3 if selected else 2
        draw.ellipse((sx(xi) - radius, sy(yi) - radius, sx(xi) + radius, sy(yi) + radius), fill=color)
    draw.text((width // 2 - 60, height - 30), xlabel, fill="#111827")
    draw.text((8, height // 2), ylabel, fill="#111827")
    image.save(path)


def phase3_figures(result_dir: Path, abundance: pd.DataFrame, detection: pd.DataFrame, evidence: pd.DataFrame) -> None:
    figure_dir = result_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    abundance_selected = evidence["abundance_core"].to_numpy(dtype=bool)
    volcano_y = -np.log10(np.maximum(abundance["bh_q_value"].to_numpy(dtype=float), 1e-300))
    _scatter(
        figure_dir / "abundance_volcano.png",
        abundance["mean_log2_T_minus_N"].to_numpy(dtype=float), volcano_y, abundance_selected,
        "Paired abundance with empirical-Bayes moderation", "Mean log2 tumour minus non-tumour", "-log10 BH q-value",
    )
    detection_selected = evidence["detection_core"].to_numpy(dtype=bool)
    detection_y = -np.log10(np.maximum(detection["bh_q_value"].to_numpy(dtype=float), 1e-300))
    _scatter(
        figure_dir / "detection_volcano.png",
        detection["detection_fraction_T_minus_N"].to_numpy(dtype=float), detection_y, detection_selected,
        "Paired detection differences", "Detection fraction tumour minus non-tumour", "-log10 BH q-value",
    )


def phase4_figures(result_dir: Path, patient_scores: pd.DataFrame, pathway_results: pd.DataFrame, consensus: np.ndarray | None) -> None:
    figure_dir = result_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    width, height, margin = 900, 620, 70
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    x, y = patient_scores["PC1"].to_numpy(), patient_scores["PC2"].to_numpy()
    xmin, xmax, ymin, ymax = float(x.min()), float(x.max()), float(y.min()), float(y.max())
    sx = lambda value: margin + (float(value) - xmin) / max(xmax - xmin, 1e-9) * (width - 2 * margin)
    sy = lambda value: height - margin - (float(value) - ymin) / max(ymax - ymin, 1e-9) * (height - 2 * margin)
    draw.text((margin, 24), "Patient tumour-minus-non-tumour change PCA", fill="#111827")
    draw.line((margin, height - margin, width - margin, height - margin), fill="#374151", width=2)
    draw.line((margin, margin, margin, height - margin), fill="#374151", width=2)
    for xi, yi, label in zip(x, y, patient_scores["patient_id"]):
        draw.ellipse((sx(xi) - 4, sy(yi) - 4, sx(xi) + 4, sy(yi) + 4), fill="#0F766E")
        draw.text((sx(xi) + 5, sy(yi) - 5), str(label), fill="#374151")
    draw.text((width // 2 - 20, height - 30), "PC1", fill="#111827")
    draw.text((12, height // 2), "PC2", fill="#111827")
    image.save(figure_dir / "patient_change_pca.png")

    supported = pathway_results.loc[pathway_results["supported_both_views"]].copy()
    display_pool = supported if len(supported) >= 20 else pathway_results.copy()
    display_pool["absolute_display_score"] = display_pool["display_score"].abs()
    top = display_pool.sort_values(
        ["absolute_display_score", "combined_priority", "pathway_name"],
        ascending=[False, False, True],
    ).head(20).sort_values("display_score").copy()
    bar_width, bar_height = 1100, 680
    image = Image.new("RGB", (bar_width, bar_height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((20, 18), "Largest concordant pathway effects (ranked NES)", fill="#111827")
    maximum = max(float(top["display_score"].abs().max()), 1e-9)
    center = 650
    for row, (_, item) in enumerate(top.iterrows()):
        y0 = 48 + row * 30
        label = str(item["pathway_name"])[:72]
        draw.text((20, y0), label, fill="#111827")
        length = int(abs(float(item["display_score"])) / maximum * 380)
        color = "#DC2626" if item["display_score"] > 0 else "#2563EB"
        if item["display_score"] > 0:
            draw.rectangle((center, y0, center + length, y0 + 14), fill=color)
        else:
            draw.rectangle((center - length, y0, center, y0 + 14), fill=color)
    draw.line((center, 42, center, bar_height - 20), fill="#6B7280", width=1)
    image.save(figure_dir / "pathway_evidence.png")

    if consensus is not None:
        order = np.lexsort((patient_scores["PC1"].to_numpy(), patient_scores["exploratory_cluster"].to_numpy()))
        consensus = consensus[np.ix_(order, order)]
        n, cell, heat_margin = len(consensus), 12, 45
        size = heat_margin + n * cell + 20
        image = Image.new("RGB", (size, size), "white")
        draw = ImageDraw.Draw(image)
        draw.text((12, 12), "Consensus matrix for best evaluated k", fill="#111827")
        for row in range(n):
            for column in range(n):
                value = float(np.clip(consensus[row, column], 0, 1))
                color = (int(245 - 210 * value), int(248 - 125 * value), int(250 - 80 * value))
                x0, y0 = heat_margin + column * cell, heat_margin + row * cell
                draw.rectangle((x0, y0, x0 + cell, y0 + cell), fill=color)
        image.save(figure_dir / "consensus_clustering.png")

