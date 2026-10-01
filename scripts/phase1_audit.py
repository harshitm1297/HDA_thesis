"""Create reproducible Phase 1 records for the oral-cancer proteomics matrix."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_WORKBOOK = PROJECT_DIR / "OC_Dataset_84Sample_v2_24012024.xlsx"
OUTPUT_DIR = PROJECT_DIR / "phase1_outputs"
SAMPLE_PATTERN = re.compile(r"^P(?P<patient>\d+)(?P<tissue>[NT])$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def excel_column_name(one_based_index: int) -> str:
    value = one_based_index
    letters: list[str] = []
    while value:
        value, remainder = divmod(value - 1, 26)
        letters.append(chr(65 + remainder))
    return "".join(reversed(letters))


def build_manifest(sample_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    unknown_fields = [
        "oral_subsite",
        "histology",
        "stage",
        "grade",
        "node_status",
        "age",
        "sex",
        "tobacco_exposure",
        "smokeless_tobacco_or_betel_quid",
        "alcohol_exposure",
        "hpv_status",
        "treatment_before_collection",
        "tumor_purity",
        "inflammation",
        "non_tumor_distance_and_pathology",
        "collection_details",
        "sample_preparation_batch",
        "instrument_batch",
        "run_order",
        "technical_replicate_id",
        "notes",
    ]

    for source_position, sample_id in enumerate(sample_columns, start=2):
        match = SAMPLE_PATTERN.fullmatch(str(sample_id))
        if match is None:
            raise ValueError(f"Unrecognized sample column: {sample_id}")
        patient_number = int(match.group("patient"))
        tissue_code = match.group("tissue")
        row: dict[str, object] = {
            "sample_id": sample_id,
            "patient_id": f"P{patient_number}",
            "patient_number": patient_number,
            "pair_id": f"P{patient_number}",
            "tissue_code": tissue_code,
            "tissue_status": "tumor" if tissue_code == "T" else "matched_non_tumor",
            "source_sheet": "Sheet1",
            "source_excel_column": excel_column_name(source_position),
            "source_column_position": source_position,
            "included_initially": True,
            "exclusion_reason": "",
        }
        row.update({field: "unknown" for field in unknown_fields})
        rows.append(row)

    manifest = pd.DataFrame(rows)
    return manifest.sort_values(["patient_number", "tissue_code"], kind="stable").reset_index(drop=True)


def build_feature_audit(data: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    identifier_column = str(data.columns[0])
    sample_columns = [str(column) for column in data.columns[1:]]
    abundance = data[sample_columns].apply(pd.to_numeric, errors="coerce")
    identifiers = data[identifier_column].astype("string")
    duplicate_count = identifiers.map(identifiers.value_counts(dropna=False)).astype("Int64")

    patients = sorted(manifest["patient_number"].unique())
    tumor_columns = [f"P{patient}T" for patient in patients]
    non_tumor_columns = [f"P{patient}N" for patient in patients]
    tumor_detected = abundance[tumor_columns].notna().to_numpy()
    non_tumor_detected = abundance[non_tumor_columns].notna().to_numpy()

    contains_group = identifiers.str.contains(";", regex=False, na=False)
    group_member_count = identifiers.fillna("").map(
        lambda value: len([part for part in value.split(";") if part]) if value else 0
    )

    audit = pd.DataFrame(
        {
            "feature_key": [f"feature_{index:05d}" for index in range(1, len(data) + 1)],
            "source_sheet": "Sheet1",
            "source_excel_row": np.arange(2, len(data) + 2),
            "source_identifier_column": identifier_column,
            "feature_id_original": identifiers,
            "identifier_missing": identifiers.isna() | identifiers.fillna("").str.strip().eq(""),
            "duplicate_identifier_count": duplicate_count,
            "duplicate_identifier": duplicate_count.gt(1),
            "semicolon_group": contains_group,
            "group_member_count": group_member_count,
            "non_missing_sample_count": abundance.notna().sum(axis=1),
            "missing_sample_count": abundance.isna().sum(axis=1),
            "missing_fraction": abundance.isna().mean(axis=1),
            "all_abundance_missing": abundance.isna().all(axis=1),
            "tumor_detected_count": tumor_detected.sum(axis=1),
            "non_tumor_detected_count": non_tumor_detected.sum(axis=1),
            "complete_pair_count": (tumor_detected & non_tumor_detected).sum(axis=1),
            "tumor_only_pair_count": (tumor_detected & ~non_tumor_detected).sum(axis=1),
            "non_tumor_only_pair_count": (~tumor_detected & non_tumor_detected).sum(axis=1),
            "both_missing_pair_count": (~tumor_detected & ~non_tumor_detected).sum(axis=1),
            "observed_minimum": abundance.min(axis=1, skipna=True),
            "observed_median": abundance.median(axis=1, skipna=True),
            "observed_maximum": abundance.max(axis=1, skipna=True),
        }
    )

    conditions = [
        audit["all_abundance_missing"],
        audit["identifier_missing"],
        audit["duplicate_identifier"] & audit["semicolon_group"],
        audit["duplicate_identifier"],
        audit["semicolon_group"],
    ]
    labels = [
        "exclude_all_abundance_missing",
        "review_missing_identifier",
        "review_duplicate_and_group",
        "review_duplicate_identifier",
        "review_protein_group",
    ]
    audit["phase1_status"] = np.select(conditions, labels, default="retain_unresolved")
    audit["stable_accession"] = "unresolved"
    audit["identifier_resolution_note"] = "Original PG.Genes row retained; accession evidence unavailable."
    return audit


def build_dictionary() -> pd.DataFrame:
    records = [
        ("sample_manifest.csv", "sample_id", "Source sample column name", "source"),
        ("sample_manifest.csv", "patient_id", "Patient identifier parsed from sample_id", "derived"),
        ("sample_manifest.csv", "patient_number", "Numeric patient identifier used for sorting and validation", "derived"),
        ("sample_manifest.csv", "pair_id", "Identifier linking tumor and matched non-tumor specimens", "derived"),
        ("sample_manifest.csv", "tissue_code", "N for matched non-tumor or T for tumor", "derived"),
        ("sample_manifest.csv", "tissue_status", "Expanded tissue label", "derived"),
        ("sample_manifest.csv", "source_sheet", "Workbook worksheet containing the sample", "source"),
        ("sample_manifest.csv", "source_excel_column", "Original Excel column letter", "derived"),
        ("sample_manifest.csv", "source_column_position", "One-based original Excel column position", "derived"),
        ("sample_manifest.csv", "included_initially", "Initial inclusion flag before QC", "decision"),
        ("sample_manifest.csv", "exclusion_reason", "Reason for later exclusion, blank initially", "decision"),
        ("sample_manifest.csv", "clinical_and_technical_fields", "Fields absent from source and explicitly marked unknown", "placeholder"),
        ("feature_identifier_audit.csv", "feature_key", "Stable row-based identifier assigned in Phase 1", "derived"),
        ("feature_identifier_audit.csv", "source_excel_row", "Original one-based Excel row", "derived"),
        ("feature_identifier_audit.csv", "feature_id_original", "Unmodified PG.Genes value", "source"),
        ("feature_identifier_audit.csv", "identifier_missing", "Original identifier is blank", "derived"),
        ("feature_identifier_audit.csv", "duplicate_identifier_count", "Number of source rows with the same identifier", "derived"),
        ("feature_identifier_audit.csv", "duplicate_identifier", "Identifier occurs in more than one source row", "derived"),
        ("feature_identifier_audit.csv", "semicolon_group", "Identifier contains semicolon-delimited members", "derived"),
        ("feature_identifier_audit.csv", "group_member_count", "Number of non-empty semicolon-delimited members", "derived"),
        ("feature_identifier_audit.csv", "non_missing_sample_count", "Number of numeric sample values", "derived"),
        ("feature_identifier_audit.csv", "missing_sample_count", "Number of missing sample values", "derived"),
        ("feature_identifier_audit.csv", "missing_fraction", "Missing sample values divided by 84", "derived"),
        ("feature_identifier_audit.csv", "all_abundance_missing", "No numeric values in any sample", "derived"),
        ("feature_identifier_audit.csv", "tumor_detected_count", "Tumor specimens with a numeric value", "derived"),
        ("feature_identifier_audit.csv", "non_tumor_detected_count", "Matched non-tumor specimens with a numeric value", "derived"),
        ("feature_identifier_audit.csv", "complete_pair_count", "Pairs measured in both tissues", "derived"),
        ("feature_identifier_audit.csv", "tumor_only_pair_count", "Pairs measured only in tumor", "derived"),
        ("feature_identifier_audit.csv", "non_tumor_only_pair_count", "Pairs measured only in matched non-tumor", "derived"),
        ("feature_identifier_audit.csv", "both_missing_pair_count", "Pairs missing in both tissues", "derived"),
        ("feature_identifier_audit.csv", "observed_minimum", "Minimum observed abundance", "derived"),
        ("feature_identifier_audit.csv", "observed_median", "Median observed abundance", "derived"),
        ("feature_identifier_audit.csv", "observed_maximum", "Maximum observed abundance", "derived"),
        ("feature_identifier_audit.csv", "phase1_status", "Initial identifier and completeness disposition", "decision"),
        ("feature_identifier_audit.csv", "stable_accession", "Stable protein accession; unresolved in Phase 1", "placeholder"),
        ("feature_identifier_audit.csv", "identifier_resolution_note", "Reason the original feature unit is retained", "decision"),
    ]
    return pd.DataFrame(records, columns=["file", "field", "definition", "provenance_type"])


def validate(data: pd.DataFrame, manifest: pd.DataFrame, audit: pd.DataFrame) -> list[str]:
    sample_columns = [str(column) for column in data.columns[1:]]
    checks = {
        "Workbook contains one identifier column and 84 sample columns": len(sample_columns) == 84,
        "All sample names match P<number><N or T>": all(SAMPLE_PATTERN.fullmatch(column) for column in sample_columns),
        "Sample names are unique": len(set(sample_columns)) == len(sample_columns),
        "Manifest contains 84 specimens": len(manifest) == 84,
        "Manifest contains 42 patients": manifest["patient_id"].nunique() == 42,
        "Every patient has exactly one N and one T specimen": manifest.groupby("patient_id")["tissue_code"].apply(lambda values: sorted(values) == ["N", "T"]).all(),
        "Feature audit contains every source feature row": len(audit) == len(data),
        "Feature keys are unique": audit["feature_key"].is_unique,
        "Source Excel rows are unique": audit["source_excel_row"].is_unique,
        "Paired detection categories sum to 42 for every feature": (
            audit[["complete_pair_count", "tumor_only_pair_count", "non_tumor_only_pair_count", "both_missing_pair_count"]].sum(axis=1) == 42
        ).all(),
    }
    failed = [name for name, passed in checks.items() if not bool(passed)]
    if failed:
        raise ValueError("Phase 1 validation failed: " + "; ".join(failed))
    return list(checks)


def write_report(
    data: pd.DataFrame,
    manifest: pd.DataFrame,
    audit: pd.DataFrame,
    checks: list[str],
    output_paths: list[Path],
) -> None:
    abundance = data.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    duplicate_ids = (
        audit.loc[audit["duplicate_identifier"], "feature_id_original"].value_counts().sort_index().to_dict()
    )
    output_hashes = "\n".join(f"- `{path.name}`: `{sha256(path)}`" for path in output_paths)
    check_lines = "\n".join(f"- PASS: {check}" for check in checks)
    report = f"""# Phase 1 validation report

## Source

- Workbook: `{INPUT_WORKBOOK.name}`
- SHA-256: `{sha256(INPUT_WORKBOOK)}`
- Worksheet: `Sheet1`

## Dataset summary

- Feature rows: {len(data):,}
- Sample columns: {abundance.shape[1]:,}
- Independent patients: {manifest['patient_id'].nunique():,}
- Complete tumor/non-tumor pairs: {manifest['pair_id'].nunique():,}
- Missing abundance cells: {int(abundance.isna().sum().sum()):,} of {abundance.size:,} ({abundance.isna().to_numpy().mean():.3%})
- All-missing feature rows: {int(audit['all_abundance_missing'].sum()):,}
- Duplicated identifier rows: {int(audit['duplicate_identifier'].sum()):,}
- Distinct duplicated identifiers: {audit.loc[audit['duplicate_identifier'], 'feature_id_original'].nunique():,}
- Duplicated identifiers and counts: `{duplicate_ids}`
- Semicolon-delimited group rows: {int(audit['semicolon_group'].sum()):,}
- Semicolon-delimited group rows that are also all-missing: {int((audit['semicolon_group'] & audit['all_abundance_missing']).sum()):,}
- Missing identifiers: {int(audit['identifier_missing'].sum()):,}
- Sample detected-feature range: {int(abundance.notna().sum(axis=0).min()):,} to {int(abundance.notna().sum(axis=0).max()):,}

## Validation checks

{check_lines}

## Output hashes

{output_hashes}

## Interpretation

The source contains 42 complete biological pairs. The sample manifest is structurally valid. The feature audit preserves all 8,071 source rows and identifies rows requiring completeness or identifier review. No duplicate or grouped identifier was automatically averaged, split, or mapped to a stable accession.

All 48 semicolon-delimited group rows are among the 114 rows with no abundance measurement. They remain documented but are not candidates for quantitative modelling in the supplied matrix.

Clinical and technical metadata fields in the manifest remain `unknown` because they are not present in the supplied source files. These values must not be interpreted as negative findings.
"""
    (OUTPUT_DIR / "PHASE1_VALIDATION_REPORT.md").write_text(report, encoding="utf-8")


def main() -> None:
    if not INPUT_WORKBOOK.exists():
        raise FileNotFoundError(INPUT_WORKBOOK)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    excel = pd.ExcelFile(INPUT_WORKBOOK)
    if excel.sheet_names != ["Sheet1"]:
        raise ValueError(f"Unexpected worksheets: {excel.sheet_names}")
    data = pd.read_excel(INPUT_WORKBOOK, sheet_name="Sheet1")
    if str(data.columns[0]) != "PG.Genes":
        raise ValueError(f"Unexpected identifier column: {data.columns[0]}")

    sample_columns = [str(column) for column in data.columns[1:]]
    manifest = build_manifest(sample_columns)
    audit = build_feature_audit(data, manifest)
    dictionary = build_dictionary()
    checks = validate(data, manifest, audit)

    manifest_path = OUTPUT_DIR / "sample_manifest.csv"
    audit_path = OUTPUT_DIR / "feature_identifier_audit.csv"
    dictionary_path = OUTPUT_DIR / "phase1_data_dictionary.csv"
    manifest.to_csv(manifest_path, index=False, encoding="utf-8")
    audit.to_csv(audit_path, index=False, encoding="utf-8")
    dictionary.to_csv(dictionary_path, index=False, encoding="utf-8")
    write_report(data, manifest, audit, checks, [manifest_path, audit_path, dictionary_path])

    print(f"Created {manifest_path}")
    print(f"Created {audit_path}")
    print(f"Created {dictionary_path}")
    print(f"Created {OUTPUT_DIR / 'PHASE1_VALIDATION_REPORT.md'}")


if __name__ == "__main__":
    main()
