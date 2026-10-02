"""Data loading and representation construction with stable sample/feature axes."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ProjectData:
    sample_ids: np.ndarray
    feature_keys: np.ndarray
    patient_ids: np.ndarray
    x_raw: np.ndarray
    detected: np.ndarray
    sample_manifest: pd.DataFrame
    feature_manifest: pd.DataFrame
    cohort_membership: pd.DataFrame


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_project_data(project_dir: Path) -> ProjectData:
    raw = np.load(project_dir / "data" / "interim" / "raw_data_model.npz")
    sample_manifest = pd.read_csv(project_dir / "data" / "interim" / "sample_manifest.csv")
    feature_manifest = pd.read_csv(project_dir / "data" / "interim" / "feature_manifest.csv")
    cohort = pd.read_csv(project_dir / "results" / "01_qc_missingness" / "cohort_membership.csv")
    sample_ids = raw["sample_ids"].astype(str)
    return ProjectData(
        sample_ids=sample_ids,
        feature_keys=raw["feature_keys"].astype(str),
        patient_ids=sample_manifest["patient_id"].drop_duplicates().to_numpy(dtype=str),
        x_raw=raw["X_raw"].astype(float),
        detected=raw["D"].astype(bool),
        sample_manifest=sample_manifest,
        feature_manifest=feature_manifest,
        cohort_membership=cohort,
    )


def log2_representation(data: ProjectData, name: str) -> np.ndarray:
    matrix = np.full_like(data.x_raw, np.nan, dtype=float)
    positive = data.detected & (data.x_raw > 0)
    matrix[positive] = np.log2(data.x_raw[positive])
    if name == "log2_uncentered":
        return matrix
    if name == "log2_median_centered":
        medians = np.nanmedian(matrix, axis=1)
        return matrix - medians[:, None]
    raise ValueError(f"Unknown representation: {name}")


def paired_delta(data: ProjectData, representation: str, cohort: str) -> tuple[np.ndarray, np.ndarray]:
    matrix = log2_representation(data, representation)
    sample_index = {sample: index for index, sample in enumerate(data.sample_ids)}
    if cohort == "all_pairs":
        patients = data.patient_ids
    elif cohort == "unflagged_pairs":
        patients = data.cohort_membership.loc[
            data.cohort_membership["sensitivity_unflagged_pairs_included"], "patient_id"
        ].to_numpy(dtype=str)
    else:
        raise ValueError(f"Unknown cohort: {cohort}")
    n_index = np.asarray([sample_index[f"{patient}N"] for patient in patients])
    t_index = np.asarray([sample_index[f"{patient}T"] for patient in patients])
    return matrix[t_index] - matrix[n_index], patients

