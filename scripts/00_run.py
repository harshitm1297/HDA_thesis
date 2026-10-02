"""Run the complete Phase 0 pipeline and write a hashed output manifest."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from common import INTERIM_DIR, PROJECT_DIR, RESULT_DIR, sha256, write_json


def run(script: str) -> None:
    subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / script)], check=True)


def main() -> None:
    run("00_validate_input.py")
    run("01_build_data_model.py")
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=PROJECT_DIR, check=True)

    outputs = [
        INTERIM_DIR / "sample_manifest.csv",
        INTERIM_DIR / "feature_manifest.csv",
        INTERIM_DIR / "data_dictionary.csv",
        INTERIM_DIR / "raw_data_model.npz",
        INTERIM_DIR / "paired_data_model.npz",
        INTERIM_DIR / "data_model_metadata.json",
        RESULT_DIR / "input_validation.json",
    ]
    manifest = {
        "phase": 0,
        "milestone": "M0",
        "status": "PASS",
        "outputs": [
            {"path": str(path.relative_to(PROJECT_DIR)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        ],
    }
    write_json(RESULT_DIR / "output_manifest.json", manifest)
    print("PASS: Phase 0 complete; M0 output manifest written")


if __name__ == "__main__":
    main()

