"""Validate the immutable workbook and emit Phase 0 manifests."""

from __future__ import annotations

from common import (
    CONFIG_PATH, CONTRACT_PATH, INTERIM_DIR, PROJECT_DIR, RESULT_DIR,
    load_json_yaml, read_source, sha256, validate_source, write_json,
)


def main() -> None:
    config = load_json_yaml(CONFIG_PATH)
    contract = load_json_yaml(CONTRACT_PATH)
    data, sheets = read_source(config)
    manifest, features, _, result = validate_source(data, sheets, config)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(INTERIM_DIR / "sample_manifest.csv", index=False)
    features.to_csv(INTERIM_DIR / "feature_manifest.csv", index=False)
    result.update(
        {
            "phase": 0,
            "milestone": "M0",
            "source_workbook": config["source_workbook"],
            "source_sha256": sha256(PROJECT_DIR / config["source_workbook"]),
            "config_sha256": sha256(CONFIG_PATH),
            "contract_sha256": sha256(CONTRACT_PATH),
            "contract_version": contract["version"],
            "validation_status": "PASS",
        }
    )
    write_json(RESULT_DIR / "input_validation.json", result)
    print("PASS: immutable input, sample pairing, identifiers, and numeric domain validated")


if __name__ == "__main__":
    main()

