"""Unit and integration tests for modular Phase 3 inference."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.statistics import benjamini_hochberg, student_t_cdf, student_t_two_sided_p  # noqa: E402
from oral_cancer.inference import abundance_inference  # noqa: E402


class TestStatistics(unittest.TestCase):
    def test_student_t_reference_values(self) -> None:
        self.assertAlmostEqual(student_t_cdf(0.0, 10), 0.5, places=12)
        self.assertAlmostEqual(student_t_cdf(2.228139, 10), 0.975, places=5)
        self.assertAlmostEqual(student_t_two_sided_p(2.228139, 10), 0.05, places=4)
        self.assertGreater(student_t_two_sided_p(50.0, 40), 0.0)

    def test_bh_monotonicity_and_missing(self) -> None:
        adjusted = benjamini_hochberg(np.asarray([0.01, 0.04, 0.03, math.nan]))
        self.assertTrue(np.isnan(adjusted[3]))
        order = np.argsort(np.asarray([0.01, 0.04, 0.03]))
        self.assertTrue(np.all(np.diff(adjusted[:3][order]) >= -1e-12))
        self.assertTrue(np.all((adjusted[:3] >= 0) & (adjusted[:3] <= 1)))

    def test_synthetic_paired_effect_sign(self) -> None:
        rng = np.random.default_rng(17)
        delta = rng.normal(0, 0.5, size=(35, 25))
        delta[:, 0] += 2.0
        manifest = pd.DataFrame({
            "feature_key": [f"F{i:05d}" for i in range(25)],
            "feature_id_original": [f"G{i}" for i in range(25)],
            "duplicate_identifier": False,
            "compound_identifier": False,
            "all_abundance_missing": False,
        })
        result, _ = abundance_inference(delta, manifest, minimum_pairs=30, confidence_level=0.95)
        self.assertGreater(result.loc[0, "mean_log2_T_minus_N"], 1.5)
        self.assertLess(result.loc[0, "bh_q_value"], 0.001)
        self.assertGreater(result.loc[0, "ci_low"], 0)


class TestPhase3Outputs(unittest.TestCase):
    def test_validation_and_shapes(self) -> None:
        validation = json.loads((PROJECT_DIR / "results" / "03_paired_inference" / "03_validation.json").read_text(encoding="utf-8"))
        abundance = pd.read_csv(PROJECT_DIR / "results" / "03_paired_inference" / "paired_abundance_primary.csv")
        detection = pd.read_csv(PROJECT_DIR / "results" / "03_paired_inference" / "paired_detection_primary.csv")
        evidence = pd.read_csv(PROJECT_DIR / "results" / "03_paired_inference" / "integrated_candidate_evidence.csv")
        self.assertEqual(validation["status"], "PASS")
        self.assertTrue(all(validation["checks"].values()))
        self.assertEqual(len(abundance), 8071)
        self.assertEqual(len(detection), 8071)
        self.assertEqual(len(evidence), 8071)

    def test_primary_sign_and_pair_support(self) -> None:
        abundance = pd.read_csv(PROJECT_DIR / "results" / "03_paired_inference" / "paired_abundance_primary.csv")
        self.assertTrue((abundance.loc[abundance["eligible"], "complete_pair_count"] >= 30).all())
        self.assertTrue((abundance.loc[~abundance["eligible"], "bh_q_value"].isna()).all())


if __name__ == "__main__":
    unittest.main()

