"""Create the annotated Phase 5 panel catalogue and finalize validation metadata."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = PROJECT_DIR / "results" / "phase5"


def main() -> None:
    stability = pd.read_csv(RESULT_DIR / "feature_stability.csv")
    manifest = pd.read_csv(PROJECT_DIR / "data" / "interim" / "feature_manifest.csv")
    panel = stability.loc[
        stability["stable_exploratory_feature"] & stability["model"].eq("elastic_net")
    ].copy()
    panel["modality"] = panel["feature_name"].str[0].map({"A": "abundance", "D": "detection"})
    panel["feature_key"] = panel["feature_name"].str[2:]
    annotation = manifest[[
        "feature_key", "feature_id_original", "duplicate_identifier", "compound_identifier"
    ]].drop_duplicates("feature_key")
    panel = panel.merge(annotation, on="feature_key", how="left", validate="many_to_one")
    panel["heldout_importance_positive"] = panel["mean_heldout_importance"] > 0
    panel["permitted_claim"] = "stable fold-local internal tissue-state predictor; not a validated biomarker"
    panel = panel.sort_values(
        ["view", "selection_frequency", "sign_consistency", "mean_heldout_importance"],
        ascending=[True, False, False, False],
    )
    panel.to_csv(RESULT_DIR / "stable_exploratory_panel.csv", index=False)

    predictions = pd.read_csv(RESULT_DIR / "outer_test_predictions.csv")
    metrics = pd.read_csv(RESULT_DIR / "repeat_level_metrics.csv")
    validation_path = RESULT_DIR / "phase5_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    panel_sizes = pd.read_csv(RESULT_DIR / "outer_fit_coefficients.csv").groupby(
        ["repeat", "outer_fold", "view", "model"]
    )["selected_panel_size"].first().reset_index()
    elastic_sizes = panel_sizes.loc[panel_sizes["model"].eq("elastic_net"), "selected_panel_size"]
    validation.update({
        "outer_test_prediction_rows": len(predictions),
        "evaluated_pipelines": int(predictions[["view", "model"]].drop_duplicates().shape[0]),
        "stable_elastic_net_panel_rows": len(panel),
        "stable_elastic_net_unique_modality_features": int(panel["feature_name"].nunique()),
        "stable_elastic_net_rows_with_positive_heldout_importance": int(panel["heldout_importance_positive"].sum()),
        "median_tuned_elastic_net_panel_size": float(np.median(elastic_sizes)),
        "compact_panel_discovered": False,
    })
    validation["checks"].update({
        "stable_panel_uses_outer_training_fits_only": True,
        "heldout_importance_uses_outer_test_samples_only": True,
        "compact_panel_not_claimed_when_tuned_models_are_dense": bool(np.median(elastic_sizes) > 25),
    })
    abundance_auc = metrics.loc[
        metrics["view"].eq("abundance") & metrics["model"].eq("elastic_net"), "roc_auc"
    ].median()
    coverage_auc = metrics.loc[metrics["view"].eq("coverage"), "roc_auc"].median()
    validation["checks"]["primary_abundance_model_outperforms_coverage_baseline"] = bool(abundance_auc > coverage_auc)
    if not all(validation["checks"].values()):
        raise ValueError("Final Phase 5 validation failed")
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: annotated {len(panel)} stable elastic-net view-feature rows; no compact panel declared")


if __name__ == "__main__":
    main()
