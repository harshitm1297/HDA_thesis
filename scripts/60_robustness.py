"""Execute Phase 6 inference and pathway robustness analyses."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json, load_project_data  # noqa: E402
from oral_cancer.pathways import read_gmt  # noqa: E402
from oral_cancer.robustness import build_multiverse, leave_one_patient_out, pathway_sensitivity  # noqa: E402


RESULT_DIR = PROJECT_DIR / "results" / "phase6"


def main() -> None:
    config = load_json(PROJECT_DIR / "config" / "phase6.yml")
    phase4_config = load_json(PROJECT_DIR / "config" / "phase4.yml")
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    sensitivity = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "paired_abundance_sensitivity.csv")
    primary_abundance = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "paired_abundance_primary.csv")
    primary_abundance.insert(0, "minimum_complete_pairs", 30)
    primary_abundance.insert(0, "cohort", "all_pairs")
    primary_abundance.insert(0, "representation", "log2_uncentered")
    sensitivity = pd.concat([sensitivity, primary_abundance], ignore_index=True)
    expected = (
        len(config["multiverse_representations"]) * len(config["multiverse_cohorts"])
        * len(config["multiverse_minimum_complete_pairs"]) * len(config["multiverse_missingness_strategies"])
        * len(config["multiverse_estimators"])
    )
    scenario_summary, feature_summary, similarity, scenario_tables = build_multiverse(sensitivity, config)
    if len(scenario_summary) != expected:
        raise ValueError(f"Expected {expected} multiverse scenarios, observed {len(scenario_summary)}")
    scenario_summary.to_csv(RESULT_DIR / "multiverse_scenarios.csv", index=False)
    feature_summary.to_csv(RESULT_DIR / "multiverse_feature_stability.csv", index=False)
    similarity.to_csv(RESULT_DIR / "multiverse_scenario_similarity.csv", index=False)
    core_matrix = pd.DataFrame({"feature_key": sensitivity["feature_key"].drop_duplicates().to_numpy()})
    for scenario, table in sorted(scenario_tables.items()):
        core_matrix[scenario] = table.set_index("feature_key").loc[core_matrix["feature_key"], "core_supported"].to_numpy(dtype=bool)
    core_matrix.to_csv(RESULT_DIR / "multiverse_core_support_matrix.csv", index=False)

    paired = np.load(PROJECT_DIR / "data" / "interim" / "paired_data_model.npz")
    data = load_project_data(PROJECT_DIR)
    primary = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "integrated_candidate_evidence.csv")
    lopo_summary, lopo_detail = leave_one_patient_out(
        paired["log2_T_minus_N"].astype(float), paired["detection_state"].astype(np.uint8),
        paired["patient_ids"].astype(str), data.feature_manifest, primary, config,
    )
    lopo_summary.to_csv(RESULT_DIR / "lopo_feature_influence.csv", index=False)
    lopo_detail.to_csv(
        RESULT_DIR / "lopo_tiered_candidate_detail.csv.gz", index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )

    gene_sets = read_gmt(PROJECT_DIR / phase4_config["reactome_gmt"])
    primary_ranked = pd.read_csv(PROJECT_DIR / "results" / "phase4" / "reactome_ranked_enrichment.csv")
    pathway_summary, pathway_detail = pathway_sensitivity(
        scenario_tables, gene_sets, primary_ranked,
        phase4_config["pathway_min_measured_genes"], phase4_config["pathway_max_measured_genes"],
    )
    pathway_summary.to_csv(RESULT_DIR / "pathway_robustness.csv", index=False)
    pathway_detail.to_csv(RESULT_DIR / "pathway_sensitivity_detail.csv.gz", index=False, compression={"method": "gzip", "compresslevel": 9, "mtime": 0})

    audit = {
        "phase": 6, "stage": "inference_pathway_robustness", "status": "PASS",
        "multiverse_scenarios": len(scenario_summary), "missingness_strategy_count": len(config["multiverse_missingness_strategies"]),
        "lopo_patients": len(paired["patient_ids"]), "pathways_sensitivity_tested": len(pathway_summary),
        "checks": {
            "all_42_patients_omitted_once": bool(lopo_detail["omitted_patient"].nunique() == 42),
            "multiverse_is_prespecified_cartesian_grid": bool(len(scenario_summary) == expected),
            "missingness_strategy_not_invented": bool(config["multiverse_missingness_strategies"] == ["no_imputation"]),
            "pathway_release_remains_reactome_v86": bool(phase4_config["reactome_release"] == 86),
        },
    }
    if not all(audit["checks"].values()):
        raise ValueError("Phase 6 robustness audit failed")
    (RESULT_DIR / "robustness_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: {len(scenario_summary)} multiverse scenarios, 42 LOPO runs, {len(pathway_summary)} pathways")


if __name__ == "__main__":
    main()
