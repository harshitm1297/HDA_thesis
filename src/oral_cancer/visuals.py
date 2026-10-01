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

