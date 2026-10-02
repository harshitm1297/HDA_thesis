"""Run Phase 6, verify same-environment determinism, test, and fingerprint outputs."""

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
RESULT_DIR = PROJECT_DIR / "results" / "phase6"
EXCLUDED_FROM_GATE = {"output_manifest.json", "reproducibility_gate.json", "environment_snapshot.json"}


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
        if path.is_file() and path.name not in EXCLUDED_FROM_GATE and path.suffix.lower() not in {".md"}
    }


def execute_phase6() -> None:
    for script in ["60_robustness.py", "61_ml_influence.py", "70_build_evidence_table.py"]:
        subprocess.run([sys.executable, str(PROJECT_DIR / "scripts" / script)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--verify-against-existing", action="store_true")
    args = parser.parse_args()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    if args.verify_against_existing:
        before = core_hashes()
        execute_phase6()
        after = core_hashes()
        exact = before == after
    elif args.reuse_existing:
        exact = None
    else:
        execute_phase6()
        before = core_hashes()
        execute_phase6()
        after = core_hashes()
        exact = before == after

    if exact is not None:
        gate = {
            "phase": 6, "same_environment_phase6_exact_reproduction": exact,
            "compared_output_count": len(after), "numeric_tolerance_needed": False,
            "scope": "Phase 6 scripts rerun in the same project-local environment",
            "full_clean_end_to_end_reproduction": False,
            "full_gate_blockers": ["Phase 2 not implemented", "clean independent environment rerun not executed"],
        }
        if not exact:
            gate["changed_outputs"] = sorted(set(before) ^ set(after) | {key for key in before.keys() & after.keys() if before[key] != after[key]})
            raise ValueError(f"Phase 6 deterministic reproduction failed: {gate['changed_outputs']}")
        (RESULT_DIR / "reproducibility_gate.json").write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    try:
        git_commit = subprocess.check_output([
            "git", "-c", f"safe.directory={PROJECT_DIR.as_posix()}", "rev-parse", "HEAD"
        ], cwd=PROJECT_DIR, text=True).strip()
    except Exception:
        git_commit = "unavailable"
    environment = {
        "python": sys.version, "platform": platform.platform(), "git_input_commit": git_commit,
        "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
        "commands": ["scripts/60_robustness.py", "scripts/61_ml_influence.py", "scripts/70_build_evidence_table.py"],
    }
    (RESULT_DIR / "environment_snapshot.json").write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=PROJECT_DIR, check=True)
    outputs = sorted(path for path in RESULT_DIR.rglob("*") if path.is_file() and path.name != "output_manifest.json")
    manifest = {
        "phase": 6, "computational_status": "PASS",
        "outputs": [{
            "path": str(path.relative_to(PROJECT_DIR)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)
        } for path in outputs],
    }
    (RESULT_DIR / "output_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PASS: Phase 6 tests and output manifest completed")


if __name__ == "__main__":
    main()
