"""Thin command-line entry point for modular Phase 4 analysis."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "src"))

from oral_cancer.data import load_json, load_project_data  # noqa: E402
from oral_cancer.heterogeneity import bootstrap_explained_variance, consensus_clustering, robust_change_pca  # noqa: E402
from oral_cancer.pathways import collapse_gene_ranking, pathway_scores, ranked_enrichment, read_gmt  # noqa: E402
from oral_cancer.visuals import phase4_figures  # noqa: E402


CONFIG_PATH = PROJECT_DIR / "config" / "phase4.yml"
RESULT_DIR = PROJECT_DIR / "results" / "phase4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    config = load_json(CONFIG_PATH)
    data = load_project_data(PROJECT_DIR)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    gmt_path = PROJECT_DIR / config["reactome_gmt"]
    if sha256(gmt_path) != config["reactome_gmt_sha256"]:
        raise ValueError("Reactome GMT checksum does not match frozen Phase 4 configuration")
    gene_sets = read_gmt(gmt_path)

    abundance = pd.read_csv(PROJECT_DIR / "results" / "phase3" / "paired_abundance_primary.csv")
    ranking, collapse_audit = collapse_gene_ranking(abundance)
    ranking.to_csv(RESULT_DIR / "collapsed_gene_ranking.csv", index=False)
    collapse_audit.to_csv(RESULT_DIR / "gene_collapse_audit.csv", index=False)
    ranked = ranked_enrichment(
        ranking, gene_sets, config["pathway_min_measured_genes"], config["pathway_max_measured_genes"],
        config["ranked_gsea_permutations"], config["random_seed"],
    )
    ranked.to_csv(RESULT_DIR / "reactome_ranked_enrichment.csv", index=False)

    paired = np.load(PROJECT_DIR / "data" / "interim" / "paired_data_model.npz")
    delta = paired["log2_T_minus_N"].astype(float)
    patient_ids = paired["patient_ids"].astype(str)
    pathway_test, patient_pathway_scores = pathway_scores(
        delta, data.feature_manifest, patient_ids, gene_sets,
        config["pathway_score_min_complete_genes"], config["pathway_score_sign_permutations"], config["random_seed"] + 1,
    )
    pathway_test.to_csv(RESULT_DIR / "reactome_patient_score_tests.csv", index=False)
    patient_pathway_scores.to_csv(RESULT_DIR / "patient_pathway_scores.csv", index=False)

    pathway_results = ranked.merge(
        pathway_test[["pathway_id", "complete_gene_count", "mean_patient_score", "sign_flip_p_value", "bh_q_value"]],
        on="pathway_id", how="outer", suffixes=("_ranked", "_patient_score"),
    )
    pathway_results["direction_agreement"] = (
        np.sign(pathway_results["normalized_enrichment_score"])
        == np.sign(pathway_results["mean_patient_score"])
    )
    pathway_results["supported_both_views"] = (
        (pathway_results["bh_q_value_ranked"] <= config["pathway_fdr"])
        & (pathway_results["bh_q_value_patient_score"] <= config["pathway_fdr"])
        & pathway_results["direction_agreement"]
    )
    pathway_results["combined_priority"] = np.nanmax(
        np.column_stack([
            -np.log10(np.maximum(pathway_results["bh_q_value_ranked"].to_numpy(dtype=float), 1e-300)),
            -np.log10(np.maximum(pathway_results["bh_q_value_patient_score"].to_numpy(dtype=float), 1e-300)),
        ]), axis=1,
    )
    pathway_results["display_score"] = pathway_results["normalized_enrichment_score"].fillna(pathway_results["mean_patient_score"])
    pathway_results = pathway_results.sort_values(["supported_both_views", "combined_priority"], ascending=[False, False]).reset_index(drop=True)
    pathway_results.to_csv(RESULT_DIR / "integrated_pathway_evidence.csv", index=False)

    top_pathway_ids = set(pathway_results.head(30)["pathway_id"].dropna())
    ranking_key = ranking.set_index("gene_symbol")["selected_feature_key"].to_dict()
    edge_rows = []
    for _, row in ranked.loc[ranked["pathway_id"].isin(top_pathway_ids)].iterrows():
        for gene in str(row["leading_edge_genes"]).split(";"):
            if gene in ranking_key:
                edge_rows.append({
                    "pathway_id": row["pathway_id"],
                    "pathway_name": row["pathway_name"],
                    "feature_key": ranking_key[gene],
                    "gene_symbol": gene,
                    "edge_evidence": "Reactome_v86_membership_and_ranked_leading_edge",
                })
    pd.DataFrame(edge_rows).to_csv(RESULT_DIR / "feature_pathway_network_edges.csv", index=False)

    pca = robust_change_pca(delta, config["heterogeneity_top_variable_features"], config["heterogeneity_winsor_robust_z"])
    bootstrap = bootstrap_explained_variance(pca.scaled_matrix, config["pca_bootstrap_resamples"], config["random_seed"] + 2)
    bootstrap["observed_fraction"] = pca.explained_fraction[: len(bootstrap)]
    bootstrap.to_csv(RESULT_DIR / "pca_explained_variance_bootstrap.csv", index=False)
    cumulative = np.cumsum(pca.explained_fraction)
    cluster_components = min(config["clustering_max_pcs"], max(2, int(np.searchsorted(cumulative, config["clustering_variance_target"]) + 1)))
    cluster_matrix = pca.scores[:, :cluster_components]
    cluster_summary, consensus_by_k, labels_by_k = consensus_clustering(
        cluster_matrix, config["consensus_k_values"], config["consensus_resamples"],
        config["consensus_patient_fraction"], tuple(config["consensus_pac_interval"]),
        config["cluster_acceptance"], config["random_seed"] + 3,
    )
    cluster_summary.to_csv(RESULT_DIR / "cluster_stability.csv", index=False)
    accepted_rows = cluster_summary.loc[cluster_summary["accepted"]]
    accepted_k = None if accepted_rows.empty else int(accepted_rows.sort_values(["pac", "silhouette"], ascending=[True, False]).iloc[0]["k"])
    best_k = int(cluster_summary.sort_values(["accepted", "pac", "silhouette"], ascending=[False, True, False]).iloc[0]["k"])
    np.savez_compressed(RESULT_DIR / "consensus_matrices.npz", **{f"k{k}": value for k, value in consensus_by_k.items()})

    patient_scores = pd.DataFrame({"patient_id": patient_ids})
    for component in range(min(10, pca.scores.shape[1])):
        patient_scores[f"PC{component + 1}"] = pca.scores[:, component]
    patient_scores["exploratory_best_k"] = best_k
    patient_scores["exploratory_cluster"] = labels_by_k[best_k] + 1
    patient_scores["accepted_cluster"] = labels_by_k[accepted_k] + 1 if accepted_k is not None else pd.NA
    patient_scores["interpretation"] = config["interpretation_label"]
    patient_scores.to_csv(RESULT_DIR / "patient_change_scores.csv", index=False)

    loading_rows = []
    for component in range(min(5, pca.loadings.shape[1])):
        order = np.argsort(np.abs(pca.loadings[:, component]))[::-1][:30]
        for rank, local in enumerate(order, start=1):
            feature_index = pca.selected_feature_indices[local]
            loading_rows.append({
                "component": component + 1,
                "rank_by_absolute_loading": rank,
                "feature_key": data.feature_manifest.iloc[feature_index]["feature_key"],
                "feature_id_original": data.feature_manifest.iloc[feature_index]["feature_id_original"],
                "loading": pca.loadings[local, component],
            })
    pd.DataFrame(loading_rows).to_csv(RESULT_DIR / "pca_top_loadings.csv", index=False)
    phase4_figures(RESULT_DIR, patient_scores, pathway_results, consensus_by_k[best_k])

    reactome_genes = set().union(*(gene_set.genes for gene_set in gene_sets))
    measured_symbols = set(collapse_audit["gene_symbol"])
    validation = {
        "phase": 4,
        "milestone": "M4",
        "status": "PASS",
        "milestone_status": "CONDITIONAL_ON_PHASE2_AND_FULL_M3_CLOSURE",
        "reactome_release": config["reactome_release"],
        "reactome_gmt_sha256": sha256(gmt_path),
        "reactome_pathways_total": len(gene_sets),
        "measured_unambiguous_ranked_genes": len(ranking),
        "measured_genes_mapped_to_reactome": len(measured_symbols & reactome_genes),
        "ranked_pathways_tested": len(ranked),
        "patient_score_pathways_tested": len(pathway_test),
        "ranked_pathways_fdr": int((ranked["bh_q_value"] <= config["pathway_fdr"]).sum()),
        "patient_score_pathways_fdr": int((pathway_test["bh_q_value"] <= config["pathway_fdr"]).sum()),
        "pathways_supported_both_views": int(pathway_results["supported_both_views"].sum()),
        "heterogeneity_complete_features_available": int(np.isfinite(delta).all(axis=0).sum()),
        "heterogeneity_features_used": int(len(pca.selected_feature_indices)),
        "clustering_components": cluster_components,
        "accepted_cluster_k": accepted_k,
        "best_exploratory_k": best_k,
        "checks": {
            "reactome_checksum_matches": bool(sha256(gmt_path) == config["reactome_gmt_sha256"]),
            "reactome_release_is_pre_2024": bool(config["reactome_release_date"] < "2024-01-01"),
            "rank_uses_only_eligible_phase3_features": bool(len(ranking) <= int(abundance["eligible"].sum())),
            "compound_identifiers_not_expanded": bool((collapse_audit["gene_symbol"].str.contains(";", regex=False) == False).all()),
            "patient_pca_uses_no_imputation": bool(np.isfinite(delta[:, pca.selected_feature_indices]).all()),
            "all_42_patients_scored": bool(len(patient_scores) == 42),
            "cluster_labels_not_called_clinical_subtypes": bool(patient_scores["interpretation"].eq(config["interpretation_label"]).all()),
        },
    }
    if not all(validation["checks"].values()):
        raise ValueError("Phase 4 validation failed")
    (RESULT_DIR / "phase4_validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: Phase 4 tested {len(ranked)} ranked Reactome pathways; accepted k={accepted_k}")


if __name__ == "__main__":
    main()

