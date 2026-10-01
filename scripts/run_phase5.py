"""Reproduce Phase 5, run tests, and fingerprint deliverables."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = PROJECT_DIR / "results" / "phase5"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-existing", action="store_true", help="Validate and fingerprint completed expensive fits")
    args = parser.parse_args()
    if not args.reuse_existing:
        subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / "50_nested_ml.py")], check=True)
        subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / "51_permutation_ml.py")], check=True)
    subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / "52_finalize_phase5.py")], check=True)
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=PROJECT_DIR, check=True)
    outputs = sorted(path for path in RESULT_DIR.rglob("*") if path.is_file() and path.name != "output_manifest.json")
    payload = {
        "phase": 5, "status": "PASS",
        "outputs": [
            {"path": str(path.relative_to(PROJECT_DIR)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        ],
    }
    (RESULT_DIR / "output_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PASS: Phase 5 output manifest written")


if __name__ == "__main__":
    main()
