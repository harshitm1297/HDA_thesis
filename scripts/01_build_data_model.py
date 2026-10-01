"""Materialise traceable Phase 0 raw, detection, and paired objects."""

from __future__ import annotations

import platform
import sys

import numpy as np
import pandas as pd

from phase0_common import (
    CONFIG_PATH, INTERIM_DIR, PROJECT_DIR, RESULT_DIR, load_json_yaml,
    read_source, sha256, validate_source, write_json,
)


def main() -> None:
    config = load_json_yaml(CONFIG_PATH)
    data, sheets = read_source(config)
    manifest, features, abundance, validation = validate_source(data, sheets, config)

    sample_ids = manifest["sample_id"].tolist()
    feature_keys = features["feature_key"].tolist()
    x_raw = abundance[sample_ids].T.to_numpy(dtype=np.float64)
    detected = np.isfinite(x_raw)
    x_log2 = np.full_like(x_raw, np.nan)
    positive = detected & (x_raw > 0)
    x_log2[positive] = np.log2(x_raw[positive])

    patients = manifest["patient_id"].drop_duplicates().tolist()
    sample_index = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    n_index = np.array([sample_index[f"{patient}N"] for patient in patients])
    t_index = np.array([sample_index[f"{patient}T"] for patient in patients])
    paired_delta = x_log2[t_index] - x_log2[n_index]
    n_detected = detected[n_index]
    t_detected = detected[t_index]
    detection_state = n_detected.astype(np.uint8) + 2 * t_detected.astype(np.uint8)

    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        INTERIM_DIR / "raw_data_model.npz",
        X_raw=x_raw,
        D=detected,
        sample_ids=np.asarray(sample_ids),
        feature_keys=np.asarray(feature_keys),
    )
    np.savez_compressed(
        INTERIM_DIR / "paired_data_model.npz",
        log2_T_minus_N=paired_delta,
        detection_state=detection_state,
        patient_ids=np.asarray(patients),
        feature_keys=np.asarray(feature_keys),
    )

    dictionary = pd.DataFrame(
        [
            ("sample_manifest.csv", "sample_id", "Original sample column", "source"),
            ("sample_manifest.csv", "patient_id", "Patient parsed from anchored sample-ID regex", "derived"),
            ("sample_manifest.csv", "tissue_code", "N=matched non-tumour; T=tumour", "derived"),
            ("feature_manifest.csv", "feature_key", "Stable row key F00001..F08071", "derived"),
            ("feature_manifest.csv", "feature_id_original", "Unmodified PG.Genes label", "source"),
            ("raw_data_model.npz", "X_raw", "84 x 8,071 float64 source values; NaN means missing", "source-shaped"),
            ("raw_data_model.npz", "D", "84 x 8,071 boolean finite-value detection matrix", "derived"),
            ("paired_data_model.npz", "log2_T_minus_N", "42 x 8,071 paired log2 difference; NaN unless both values observed", "derived"),
            ("paired_data_model.npz", "detection_state", "42 x 8,071: 0=neither, 1=N-only, 2=T-only, 3=both", "derived"),
        ],
        columns=["file", "field", "definition", "provenance_type"],
    )
    dictionary.to_csv(INTERIM_DIR / "data_dictionary.csv", index=False)

    complete_counts = np.isfinite(paired_delta).sum(axis=0)
    metadata = {
        "shapes": {
            "X_raw": list(x_raw.shape),
            "D": list(detected.shape),
            "log2_T_minus_N": list(paired_delta.shape),
            "detection_state": list(detection_state.shape),
        },
        "orientation": "rows are samples/patients; columns are stable feature_keys",
        "detection_state_codebook": {"0": "neither", "1": "N-only", "2": "T-only", "3": "both"},
        "log2_rule": config["abundance_transform"],
        "paired_delta_sign": config["tumour_minus_normal_sign"],
        "features_with_at_least_one_complete_pair": int((complete_counts > 0).sum()),
        "features_complete_in_all_42_pairs": int((complete_counts == 42).sum()),
        "source_sha256": sha256(PROJECT_DIR / config["source_workbook"]),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
        },
        "validation_metrics": validation["metrics"],
    }
    write_json(INTERIM_DIR / "data_model_metadata.json", metadata)
    print("PASS: raw and paired data models created with stable axes")


if __name__ == "__main__":
    main()

