"""Shared, deterministic Phase 0 data-contract functions."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_DIR / "config" / "analysis.yml"
CONTRACT_PATH = PROJECT_DIR / "config" / "analysis_contract.yml"
INTERIM_DIR = PROJECT_DIR / "data" / "interim"
RESULT_DIR = PROJECT_DIR / "results" / "00_input_validation"


def load_json_yaml(path: Path) -> dict:
    """Load JSON syntax stored in a .yml file (JSON is valid YAML 1.2)."""
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def excel_column_name(one_based_index: int) -> str:
    letters: list[str] = []
    while one_based_index:
        one_based_index, remainder = divmod(one_based_index - 1, 26)
        letters.append(chr(65 + remainder))
    return "".join(reversed(letters))


def read_source(config: dict) -> tuple[pd.DataFrame, list[str]]:
    workbook = PROJECT_DIR / config["source_workbook"]
    if not workbook.exists():
        raise FileNotFoundError(workbook)
    excel = pd.ExcelFile(workbook, engine="openpyxl")
    data = pd.read_excel(workbook, sheet_name=config["worksheet"], engine="openpyxl")
    return data, excel.sheet_names


def parse_samples(columns: list[str], config: dict) -> pd.DataFrame:
    pattern = re.compile(config["sample_regex"])
    rows = []
    for source_position, sample_id in enumerate(columns, start=2):
        match = pattern.fullmatch(sample_id)
        if match is None:
            raise ValueError(f"Sample ID violates contract: {sample_id!r}")
        patient_number = int(match.group("patient"))
        tissue_code = match.group("tissue")
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": f"P{patient_number}",
                "patient_number": patient_number,
                "pair_id": f"P{patient_number}",
                "tissue_code": tissue_code,
                "tissue_status": "tumour" if tissue_code == "T" else "matched_non_tumour",
                "source_sheet": config["worksheet"],
                "source_excel_column": excel_column_name(source_position),
                "source_column_position": source_position,
            }
        )
    return pd.DataFrame(rows).sort_values(["patient_number", "tissue_code"], kind="stable").reset_index(drop=True)


def coerce_abundance(data: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    source = data.iloc[:, 1:]
    numeric = source.apply(pd.to_numeric, errors="coerce")
    nonblank = source.notna() & source.astype("string").apply(lambda column: column.str.strip().ne(""))
    coercion_failures = int((nonblank & numeric.isna()).to_numpy().sum())
    return numeric, coercion_failures


def build_feature_manifest(data: pd.DataFrame, abundance: pd.DataFrame, config: dict) -> pd.DataFrame:
    identifiers = data.iloc[:, 0].astype("string")
    counts = identifiers.value_counts(dropna=False)
    duplicate_count = identifiers.map(counts).astype("Int64")
    return pd.DataFrame(
        {
            "feature_key": [f"F{i:05d}" for i in range(1, len(data) + 1)],
            "source_sheet": config["worksheet"],
            "source_excel_row": np.arange(2, len(data) + 2),
            "feature_id_original": identifiers,
            "identifier_missing": identifiers.isna() | identifiers.fillna("").str.strip().eq(""),
            "duplicate_identifier_count": duplicate_count,
            "duplicate_identifier": duplicate_count.gt(1),
            "compound_identifier": identifiers.str.contains(";", regex=False, na=False),
            "observed_count": abundance.notna().sum(axis=1),
            "missing_count": abundance.isna().sum(axis=1),
            "missing_fraction": abundance.isna().mean(axis=1),
            "all_abundance_missing": abundance.isna().all(axis=1),
        }
    )


def validate_source(data: pd.DataFrame, sheets: list[str], config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    sample_columns = [str(column) for column in data.columns[1:]]
    manifest = parse_samples(sample_columns, config)
    abundance, coercion_failures = coerce_abundance(data)
    features = build_feature_manifest(data, abundance, config)
    values = abundance.to_numpy(dtype=float)
    finite = np.isfinite(values)
    pairing = manifest.groupby("patient_id")["tissue_code"].apply(lambda x: sorted(x.tolist()) == ["N", "T"])
    checks = {
        "worksheet_exact": sheets == [config["worksheet"]],
        "identifier_column_exact": str(data.columns[0]) == config["identifier_column"],
        "feature_count_exact": len(data) == config["expected_features"],
        "sample_count_exact": len(sample_columns) == config["expected_samples"],
        "sample_ids_unique": len(sample_columns) == len(set(sample_columns)),
        "patient_count_exact": manifest["patient_id"].nunique() == config["expected_patients"],
        "every_patient_has_one_N_and_one_T": bool(pairing.all()),
        "no_numeric_coercion_failures": coercion_failures == 0,
        "no_infinite_values": int(np.isinf(values).sum()) == 0,
        "no_negative_values": int((values < 0).sum()) == 0,
        "no_zero_values": int((values == 0).sum()) == 0,
        "no_missing_identifiers": int(features["identifier_missing"].sum()) == 0,
        "feature_keys_unique": bool(features["feature_key"].is_unique),
    }
    metrics = {
        "feature_rows": len(data),
        "sample_columns": len(sample_columns),
        "patients": int(manifest["patient_id"].nunique()),
        "complete_pairs": int(pairing.sum()),
        "numeric_values": int(finite.sum()),
        "missing_values": int(np.isnan(values).sum()),
        "missing_fraction": float(np.isnan(values).mean()),
        "zeros": int((values == 0).sum()),
        "negative_values": int((values < 0).sum()),
        "infinite_values": int(np.isinf(values).sum()),
        "numeric_coercion_failures": coercion_failures,
        "minimum_positive": float(values[finite & (values > 0)].min()),
        "maximum_observed": float(values[finite].max()),
        "all_missing_feature_rows": int(features["all_abundance_missing"].sum()),
        "duplicate_identifier_rows": int(features["duplicate_identifier"].sum()),
        "distinct_duplicated_identifiers": int(features.loc[features["duplicate_identifier"], "feature_id_original"].nunique()),
        "compound_identifier_rows": int(features["compound_identifier"].sum()),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Phase 0 source validation failed: " + "; ".join(failed))
    return manifest, features, abundance, {"checks": checks, "metrics": metrics}


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

