"""Leakage-failure fixtures and output contracts for Phase 5."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.ml import (  # noqa: E402
    fit_fold_transformer,
    make_group_folds,
    reject_external_supervised_filter,
    reject_outer_test_hyperparameter_selection,
    validate_grouped_split,
)
from oral_cancer.modelling import load_ml_data, paired_swapped_labels  # noqa: E402


class TestLeakageGuards(unittest.TestCase):
    def test_patient_overlap_deliberately_fails(self) -> None:
        patients = np.asarray(["P1", "P1", "P2", "P2"])
        with self.assertRaisesRegex(ValueError, "patient occurs"):
            validate_grouped_split(np.asarray([0, 1, 2]), np.asarray([3]), patients)

    def test_split_pair_deliberately_fails(self) -> None:
        patients = np.asarray(["P1", "P1", "P2", "P2"])
        with self.assertRaisesRegex(ValueError, "patient occurs|paired specimens"):
            validate_grouped_split(np.asarray([0, 2, 3]), np.asarray([1]), patients)

    def test_transformer_fitted_on_test_deliberately_fails(self) -> None:
        log2 = np.asarray([[1.0, 2.0], [1.1, 2.1], [3.0, 4.0], [3.1, 4.1]])
        detection = np.isfinite(log2)
        labels = np.asarray([0, 1, 0, 1])
        with self.assertRaisesRegex(ValueError, "non-training sample"):
            fit_fold_transformer(
                log2, detection, labels, np.asarray([0, 1, 2]), np.asarray([0, 1]),
                np.asarray(["F1", "F2"]), "abundance", 2, 0.5, 0.35,
            )

    def test_full_data_supervised_filter_deliberately_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "inside the training fold"):
            reject_external_supervised_filter("phase3_full_data_significant_features")

    def test_outer_test_tuning_deliberately_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "outer-test"):
            reject_outer_test_hyperparameter_selection(0.99)

    def test_group_fold_size_and_paired_swap(self) -> None:
        data = load_ml_data(PROJECT_DIR)
        mapping = make_group_folds(np.unique(data.patient_ids), 6, 123)
        self.assertTrue(pd.Series(mapping).value_counts().eq(7).all())
        swapped = paired_swapped_labels(data, 99)
        for patient in np.unique(data.patient_ids):
            self.assertEqual(set(swapped[data.patient_ids == patient]), {0, 1})


class TestPhase5Outputs(unittest.TestCase):
    def test_phase5_validation_and_outer_predictions(self) -> None:
        path = PROJECT_DIR / "results" / "05_machine_learning" / "05_validation.json"
        if not path.exists():
            self.skipTest("Phase 5 has not been executed yet")
        validation = json.loads(path.read_text(encoding="utf-8"))
        predictions = pd.read_csv(PROJECT_DIR / "results" / "05_machine_learning" / "outer_test_predictions.csv")
        self.assertEqual(validation["status"], "PASS")
        self.assertTrue(all(validation["checks"].values()))
        self.assertEqual(validation["permutation_count"], 200)
        self.assertTrue(predictions.groupby(["repeat", "view", "model"]).size().eq(84).all())
        self.assertTrue(predictions["probability"].between(0, 1).all())


if __name__ == "__main__":
    unittest.main()
