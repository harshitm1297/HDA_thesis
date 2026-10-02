"""Contract tests for the first implemented project phase."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from common import CONFIG_PATH, load_json_yaml, parse_samples  # noqa: E402


class TestPairingContract(unittest.TestCase):
    def test_sample_regex_is_anchored(self) -> None:
        config = load_json_yaml(CONFIG_PATH)
        with self.assertRaises(ValueError):
            parse_samples(["prefix_P1N"], config)

    def test_manifest_has_complete_pairs(self) -> None:
        manifest = pd.read_csv(PROJECT_DIR / "data" / "interim" / "sample_manifest.csv")
        self.assertEqual(len(manifest), 84)
        self.assertEqual(manifest["patient_id"].nunique(), 42)
        self.assertTrue(manifest.groupby("patient_id")["tissue_code"].apply(lambda x: sorted(x) == ["N", "T"]).all())


class TestDataObjects(unittest.TestCase):
    def test_axes_and_pair_sign(self) -> None:
        raw = np.load(PROJECT_DIR / "data" / "interim" / "raw_data_model.npz")
        paired = np.load(PROJECT_DIR / "data" / "interim" / "paired_data_model.npz")
        self.assertEqual(raw["X_raw"].shape, (84, 8071))
        self.assertEqual(raw["D"].shape, (84, 8071))
        self.assertEqual(paired["log2_T_minus_N"].shape, (42, 8071))
        self.assertEqual(paired["detection_state"].shape, (42, 8071))
        sample_index = {sample: i for i, sample in enumerate(raw["sample_ids"].tolist())}
        feature = next(i for i in range(8071) if np.isfinite(raw["X_raw"][sample_index["P1N"], i]) and np.isfinite(raw["X_raw"][sample_index["P1T"], i]))
        expected = np.log2(raw["X_raw"][sample_index["P1T"], feature]) - np.log2(raw["X_raw"][sample_index["P1N"], feature])
        self.assertAlmostEqual(paired["log2_T_minus_N"][0, feature], expected)

    def test_validation_and_contract_are_frozen(self) -> None:
        validation = json.loads((PROJECT_DIR / "results" / "00_input_validation" / "input_validation.json").read_text(encoding="utf-8"))
        contract = json.loads((PROJECT_DIR / "config" / "analysis_contract.yml").read_text(encoding="utf-8"))
        self.assertEqual(validation["validation_status"], "PASS")
        self.assertTrue(all(validation["checks"].values()))
        self.assertEqual(contract["status"], "frozen-at-M0")
        self.assertEqual(contract["independent_unit"], "patient")


if __name__ == "__main__":
    unittest.main()

