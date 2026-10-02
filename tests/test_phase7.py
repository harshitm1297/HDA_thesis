"""Contracts for the Phase 7 validation-readiness hand-off."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.phase7 import nondominated_fronts, paired_t_required_n  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase7"


class TestPhase7Methods(unittest.TestCase):
    def test_pareto_fronts_respect_dominance(self) -> None:
        frame = pd.DataFrame({"a": [2.0, 1.0, 2.0], "b": [2.0, 1.0, 1.0]})
        np.testing.assert_array_equal(nondominated_fronts(frame, ["a", "b"]), [1, 3, 2])

    def test_sample_size_behaves_monotonically(self) -> None:
        self.assertGreater(paired_t_required_n(0.5, 0.01, 0.9), paired_t_required_n(0.5, 0.01, 0.8))
        self.assertGreater(paired_t_required_n(0.5, 0.01, 0.8), paired_t_required_n(0.8, 0.01, 0.8))


class TestPhase7Outputs(unittest.TestCase):
    def test_candidate_and_shortlist_contracts(self) -> None:
        readiness = pd.read_csv(RESULT_DIR / "candidate_validation_readiness.csv")
        shortlist = pd.read_csv(RESULT_DIR / "locked_handoff_shortlist.csv")
        self.assertEqual(len(readiness), 165)
        self.assertEqual(len(shortlist), 20)
        self.assertTrue(readiness["feature_key"].is_unique)
        self.assertTrue(shortlist["feature_key"].is_unique)
        self.assertEqual(shortlist["validation_state"].unique().tolist(), ["not_executed"])

    def test_protocol_does_not_claim_validation(self) -> None:
        protocol = json.loads((RESULT_DIR / "prospective_validation_protocol.json").read_text(encoding="utf-8"))
        self.assertTrue(protocol["independent_validation_required"])
        self.assertEqual(protocol["protocol_status"], "prospective_template_not_executed")
        self.assertEqual(protocol["phase2_dependency"], "unresolved")

    def test_validation_checks_pass(self) -> None:
        validation = json.loads((RESULT_DIR / "phase7_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(validation["computational_status"], "PASS")
        self.assertIn("VALIDATION_NOT_EXECUTED", validation["milestone_status"])
        self.assertTrue(all(validation["checks"].values()))


if __name__ == "__main__":
    unittest.main()
