"""Unit and integration tests for Phase 4 pathways and heterogeneity."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.pathways import GeneSet, enrichment_score, pathway_scores, read_gmt  # noqa: E402


class TestPathwayMethods(unittest.TestCase):
    def test_enrichment_score_direction(self) -> None:
        weights = np.asarray([5, 4, 3, 2, 1, -1, -2, -3, -4, -5], dtype=float)
        positive, _ = enrichment_score(weights, np.asarray([0, 1, 2]))
        negative, _ = enrichment_score(weights, np.asarray([7, 8, 9]))
        self.assertGreater(positive, 0)
        self.assertLess(negative, 0)

    def test_reactome_resource_is_nonempty(self) -> None:
        sets = read_gmt(PROJECT_DIR / "resources" / "reactome_v86" / "gmt" / "ReactomePathways.gmt")
        self.assertGreater(len(sets), 2500)
        self.assertTrue(any(item.stable_id.startswith("R-HSA-") for item in sets))

    def test_patient_pathway_score_preserves_paired_zero(self) -> None:
        rng = np.random.default_rng(7)
        delta = rng.normal(1.5, 0.2, size=(42, 5))
        manifest = pd.DataFrame({
            "feature_id_original": [f"G{i}" for i in range(5)],
            "duplicate_identifier": [False] * 5,
            "compound_identifier": [False] * 5,
        })
        gene_set = GeneSet("positive shift", "R-HSA-TEST", frozenset(manifest["feature_id_original"]))
        result, scores = pathway_scores(
            delta, manifest, np.asarray([f"P{i}" for i in range(42)]),
            [gene_set], minimum_genes=5, permutations=1000, seed=13,
        )
        self.assertGreater(result.loc[0, "mean_patient_score"], 0)
        self.assertLess(result.loc[0, "sign_flip_p_value"], 0.01)
        self.assertEqual(len(scores), 42)


class TestPhase4Outputs(unittest.TestCase):
    def test_validation_and_patient_outputs(self) -> None:
        validation = json.loads((PROJECT_DIR / "results" / "04_pathways_heterogeneity" / "04_validation.json").read_text(encoding="utf-8"))
        patients = pd.read_csv(PROJECT_DIR / "results" / "04_pathways_heterogeneity" / "patient_change_scores.csv")
        pathways = pd.read_csv(PROJECT_DIR / "results" / "04_pathways_heterogeneity" / "integrated_pathway_evidence.csv")
        self.assertEqual(validation["status"], "PASS")
        self.assertTrue(all(validation["checks"].values()))
        self.assertEqual(len(patients), 42)
        self.assertEqual(
            len(pathways),
            len(set(pathways["pathway_id"])),
            "Integrated evidence must contain one row per Reactome pathway.",
        )
        self.assertGreaterEqual(len(pathways), validation["ranked_pathways_tested"])
        self.assertGreaterEqual(len(pathways), validation["patient_score_pathways_tested"])
        self.assertTrue(patients["interpretation"].str.contains("not a clinical subtype").all())

    def test_pathway_sizes_and_release(self) -> None:
        ranked = pd.read_csv(PROJECT_DIR / "results" / "04_pathways_heterogeneity" / "reactome_ranked_enrichment.csv")
        self.assertTrue(ranked["measured_gene_count"].between(10, 300).all())
        validation = json.loads((PROJECT_DIR / "results" / "04_pathways_heterogeneity" / "04_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(validation["reactome_release"], 86)


if __name__ == "__main__":
    unittest.main()

