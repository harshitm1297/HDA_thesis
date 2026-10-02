"""Phase 1 QC and missingness contract tests."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "scripts"))
from importlib import import_module  # noqa: E402

phase1 = import_module("10_qc_missingness")


class TestPhase1Methods(unittest.TestCase):
    def test_exact_paired_binary(self) -> None:
        self.assertTrue(math.isnan(phase1.exact_paired_binary_pvalue(0, 0)))
        self.assertEqual(phase1.exact_paired_binary_pvalue(1, 0), 1.0)
        self.assertAlmostEqual(phase1.exact_paired_binary_pvalue(10, 0), 2 / 1024)
        self.assertEqual(phase1.exact_paired_binary_pvalue(5, 5), 1.0)


class TestPhase1Outputs(unittest.TestCase):
    def test_cohorts_preserve_pairs(self) -> None:
        cohort = pd.read_csv(PROJECT_DIR / "results" / "01_qc_missingness" / "cohort_membership.csv")
        qc = pd.read_csv(PROJECT_DIR / "results" / "01_qc_missingness" / "sample_qc.csv")
        self.assertEqual(len(cohort), 42)
        self.assertTrue(cohort["primary_all_pairs_included"].all())
        self.assertTrue(qc.groupby("patient_id")["sensitivity_pair_removed"].nunique().eq(1).all())

    def test_feature_states_and_flags(self) -> None:
        features = pd.read_csv(PROJECT_DIR / "results" / "01_qc_missingness" / "feature_missingness.csv")
        qc = pd.read_csv(PROJECT_DIR / "results" / "01_qc_missingness" / "sample_qc.csv")
        self.assertEqual(len(features), 8071)
        self.assertTrue((features[["both_missing_pairs", "N_only_pairs", "T_only_pairs", "both_detected_pairs"]].sum(axis=1) == 42).all())
        self.assertTrue((qc.loc[qc["sensitivity_flag"], "diagnostic_family_count"] >= 2).all())

    def test_validation_passes(self) -> None:
        validation = json.loads((PROJECT_DIR / "results" / "01_qc_missingness" / "01_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(validation["status"], "PASS")
        self.assertTrue(all(validation["checks"].values()))
        for name in ["pca_qc.svg", "coverage_vs_median.svg", "missingness_by_abundance_decile.svg", "sample_correlation_heatmap.svg"]:
            self.assertGreater((PROJECT_DIR / "results" / "01_qc_missingness" / "figures" / name).stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()

