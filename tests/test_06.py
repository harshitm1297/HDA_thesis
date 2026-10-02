"""Contracts for Phase 6 robustness, integration, and claim control."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = PROJECT_DIR / "results" / "06_robustness"


class TestPhase6Outputs(unittest.TestCase):
    def test_multiverse_is_frozen_and_complete(self) -> None:
        scenarios = pd.read_csv(RESULT_DIR / "multiverse_scenarios.csv")
        self.assertEqual(len(scenarios), 24)
        self.assertEqual(set(scenarios["missingness_strategy"]), {"no_imputation"})
        self.assertEqual(set(scenarios["estimator"]), {"moderated_t", "ordinary_t"})
        self.assertEqual(scenarios["scenario_id"].nunique(), 24)

    def test_lopo_and_integrated_axes(self) -> None:
        lopo = pd.read_csv(RESULT_DIR / "lopo_feature_influence.csv")
        evidence = pd.read_csv(RESULT_DIR / "integrated_evidence_table.csv")
        self.assertEqual(len(lopo), 8071)
        self.assertEqual(len(evidence), 8071)
        self.assertEqual(evidence["feature_key"].nunique(), 8071)
        self.assertTrue(evidence["limitation_ids"].notna().all())
        self.assertTrue(evidence["allowed_claim"].notna().all())
        self.assertTrue(evidence["ml_max_selection_frequency"].isna().any())

    def test_ml_lopo_preserves_pairs(self) -> None:
        predictions = pd.read_csv(RESULT_DIR / "ml_lopo_predictions.csv")
        self.assertEqual(len(predictions), 84)
        self.assertTrue(predictions.groupby("left_out_patient").size().eq(2).all())
        self.assertTrue(predictions.groupby("left_out_patient")["label"].nunique().eq(2).all())

    def test_phase6_stays_open_for_real_gap(self) -> None:
        validation = json.loads((RESULT_DIR / "06_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(validation["computational_status"], "PASS")
        self.assertIn("NOT_CLOSED_PHASE2", validation["milestone_status"])
        self.assertTrue(all(validation["checks"].values()))

    def test_claim_ledger_is_limitations_aware(self) -> None:
        claims = pd.read_csv(RESULT_DIR / "claim_ledger.csv")
        self.assertGreaterEqual(len(claims), 6)
        self.assertTrue(claims["limitation_ids"].notna().all())
        self.assertTrue(claims["allowed_scope"].notna().all())


if __name__ == "__main__":
    unittest.main()
