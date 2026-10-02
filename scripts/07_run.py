"""Run, test, reproduce and fingerprint the Phase 7 readiness package."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy
import pandas
import scipy
import sklearn


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = PROJECT_DIR / "results" / "07_validation_readiness"
EXCLUDED = {"output_manifest.json", "reproducibility_gate.json", "environment_snapshot.json"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def core_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(PROJECT_DIR)).replace("\\", "/"): sha256(path)
        for path in sorted(RESULT_DIR.rglob("*"))
        if path.is_file() and path.name not in EXCLUDED and path.suffix.lower() != ".md"
    }


def execute() -> None:
    subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / "80_validation_readiness.py")], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--verify-against-existing", action="store_true")
    args = parser.parse_args()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    exact = None
    if args.verify_against_existing:
        before = core_hashes()
        execute()
        after = core_hashes()
        exact = before == after
    elif not args.reuse_existing:
        execute()
        before = core_hashes()
        execute()
        after = core_hashes()
        exact = before == after
    if exact is not None:
        gate = {
            "phase": 7,
            "same_environment_phase7_exact_reproduction": exact,
            "compared_output_count": len(after),
            "validation_executed": False,
            "scope": "deterministic readiness-package generation only",
        }
        (RESULT_DIR / "reproducibility_gate.json").write_text(
            json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if not exact:
            raise ValueError("Phase 7 deterministic reproduction failed")
    try:
        git_commit = subprocess.check_output(
            ["git", "-c", f"safe.directory={PROJECT_DIR.as_posix()}", "rev-parse", "HEAD"],
            cwd=PROJECT_DIR, text=True,
        ).strip()
    except Exception:
        git_commit = "unavailable"
    environment = {
        "python": sys.version, "platform": platform.platform(), "git_input_commit": git_commit,
        "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__, "command": "scripts/80_validation_readiness.py",
    }
    (RESULT_DIR / "environment_snapshot.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=PROJECT_DIR, check=True)
    outputs = sorted(path for path in RESULT_DIR.rglob("*") if path.is_file() and path.name != "output_manifest.json")
    manifest = {
        "phase": 7, "computational_status": "PASS",
        "outputs": [{
            "path": str(path.relative_to(PROJECT_DIR)).replace("\\", "/"),
            "bytes": path.stat().st_size, "sha256": sha256(path),
        } for path in outputs],
    }
    (RESULT_DIR / "output_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("PASS: Phase 7 readiness package, tests and manifest completed")


if __name__ == "__main__":
    main()
