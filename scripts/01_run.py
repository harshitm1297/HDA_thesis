"""Run Phase 1, tests, and write reproducibility hashes."""

from __future__ import annotations

import subprocess
import sys

from common import PROJECT_DIR, sha256, write_json


def main() -> None:
    subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / "10_qc_missingness.py")], check=True)
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=PROJECT_DIR, check=True)
    result_dir = PROJECT_DIR / "results" / "01_qc_missingness"
    outputs = sorted(path for path in result_dir.rglob("*") if path.is_file() and path.name != "output_manifest.json")
    write_json(result_dir / "output_manifest.json", {
        "phase": 1,
        "milestone": "M1",
        "status": "PASS",
        "outputs": [{"path": str(path.relative_to(PROJECT_DIR)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs],
    })
    print("PASS: M1 output manifest written")


if __name__ == "__main__":
    main()

