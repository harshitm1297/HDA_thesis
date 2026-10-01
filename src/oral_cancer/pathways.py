"""Versioned Reactome parsing, ranked enrichment, and paired pathway scores."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .statistics import benjamini_hochberg


@dataclass(frozen=True)
class GeneSet:
    name: str
    stable_id: str
    genes: frozenset[str]


def read_gmt(path: Path) -> list[GeneSet]:
    sets: list[GeneSet] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip().split("\t")
        if len(fields) >= 3:
            sets.append(GeneSet(fields[0], fields[1], frozenset(fields[2:])))
    return sets


def collapse_gene_ranking(abundance: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tested = abundance.loc[
        abundance["eligible"]
        & ~abundance["compound_identifier"]
        & abundance["feature_id_original"].notna()
    ].copy()
    tested["absolute_statistic"] = tested["moderated_t"].abs()
    tested = tested.sort_values(["feature_id_original", "absolute_statistic", "feature_key"], ascending=[True, False, True])
    audit = tested.groupby("feature_id_original").agg(
        source_row_count=("feature_key", "size"),
        selected_feature_key=("feature_key", "first"),
        selected_moderated_t=("moderated_t", "first"),
    ).reset_index().rename(columns={"feature_id_original": "gene_symbol"})
    audit["collapse_rule"] = np.where(audit["source_row_count"] > 1, "largest_absolute_moderated_t", "single_source_row")
    ranking = audit[["gene_symbol", "selected_feature_key", "selected_moderated_t"]].rename(columns={"selected_moderated_t": "statistic"})
    ranking = ranking.sort_values(["statistic", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    return ranking, audit


def enrichment_score(weights: np.ndarray, hit_indices: np.ndarray) -> tuple[float, int]:
    n = len(weights)
    hits = np.sort(np.asarray(hit_indices, dtype=int))
    size = len(hits)
    if size == 0 or size == n:
        return np.nan, -1
    hit_weights = np.abs(weights[hits])
    total = float(hit_weights.sum())
    hit_weights = np.full(size, 1.0 / size) if total == 0 else hit_weights / total
    cumulative = np.cumsum(hit_weights)
    miss_penalty = 1.0 / (n - size)
    misses_through_hit = hits - np.arange(size)
    after_hit = cumulative - misses_through_hit * miss_penalty
    before_hit = np.concatenate(([0.0], after_hit[:-1])) - np.concatenate(([hits[0]], np.diff(hits) - 1)) * miss_penalty
    maximum_index = int(np.argmax(after_hit))
    minimum_index = int(np.argmin(before_hit))
    maximum, minimum = float(after_hit[maximum_index]), float(before_hit[minimum_index])
    return (maximum, maximum_index) if abs(maximum) >= abs(minimum) else (minimum, minimum_index)


def ranked_enrichment(ranking: pd.DataFrame, gene_sets: list[GeneSet], minimum: int, maximum: int, permutations: int, seed: int) -> pd.DataFrame:
    genes = ranking["gene_symbol"].astype(str).to_numpy()
    weights = ranking["statistic"].to_numpy(dtype=float)
    positions = {gene: index for index, gene in enumerate(genes)}
    rng = np.random.default_rng(seed)
    rows = []
    for gene_set in gene_sets:
        overlap = sorted(gene_set.genes.intersection(positions))
        if not minimum <= len(overlap) <= maximum:
            continue
        hit_indices = np.asarray([positions[gene] for gene in overlap])
        observed, peak_index = enrichment_score(weights, hit_indices)
        null = np.empty(permutations)
        for permutation in range(permutations):
            random_hits = rng.choice(len(weights), size=len(hit_indices), replace=False)
            null[permutation], _ = enrichment_score(weights, random_hits)
        same_sign = null[null * observed > 0]
        denominator = float(np.mean(np.abs(same_sign))) if len(same_sign) else float(np.mean(np.abs(null)))
        nes = observed / denominator if denominator > 0 else np.nan
        empirical_p = (1 + int((np.abs(null) >= abs(observed)).sum())) / (permutations + 1)
        ordered_hits = np.sort(hit_indices)
        if observed >= 0:
            leading_positions = ordered_hits[: peak_index + 1]
        else:
            leading_positions = ordered_hits[peak_index:]
        rows.append({
            "pathway_id": gene_set.stable_id,
            "pathway_name": gene_set.name,
            "measured_gene_count": len(overlap),
            "enrichment_score": observed,
            "normalized_enrichment_score": nes,
            "empirical_p_value": empirical_p,
            "leading_edge_gene_count": len(leading_positions),
            "leading_edge_genes": ";".join(genes[leading_positions]),
            "measured_genes": ";".join(overlap),
        })
    result = pd.DataFrame(rows)
    result["bh_q_value"] = benjamini_hochberg(result["empirical_p_value"].to_numpy())
    return result.sort_values(["bh_q_value", "empirical_p_value", "normalized_enrichment_score"], ascending=[True, True, False]).reset_index(drop=True)


def pathway_scores(
    patient_delta: np.ndarray,
    feature_manifest: pd.DataFrame,
    patient_ids: np.ndarray,
    gene_sets: list[GeneSet],
    minimum_genes: int,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    complete = np.isfinite(patient_delta).all(axis=0)
    symbols = feature_manifest["feature_id_original"].astype("string").to_numpy()
    unambiguous = (~feature_manifest["duplicate_identifier"] & ~feature_manifest["compound_identifier"]).to_numpy()
    eligible_indices = np.flatnonzero(complete & unambiguous & pd.notna(symbols))
    gene_to_index = {str(symbols[index]): index for index in eligible_indices}
    selected_sets = []
    for gene_set in gene_sets:
        genes = sorted(gene_set.genes.intersection(gene_to_index))
        if len(genes) >= minimum_genes:
            selected_sets.append((gene_set, genes))
    # Keep zero anchored at "no tumour-minus-normal change".  Subtracting the
    # observed patient median here would erase the location shift that the
    # paired sign-flip test is intended to detect.  MAD is used only as a
    # robust feature-wise scale; SD and then one are deterministic fallbacks
    # for features whose patient deltas have zero MAD.
    feature_median = np.median(patient_delta[:, eligible_indices], axis=0)
    robust_scale = np.median(
        np.abs(patient_delta[:, eligible_indices] - feature_median), axis=0
    ) * 1.4826
    standard_scale = np.std(patient_delta[:, eligible_indices], axis=0, ddof=1)
    robust_scale = np.where(robust_scale > 0, robust_scale, standard_scale)
    robust_scale = np.where(robust_scale > 0, robust_scale, 1.0)
    scaled = patient_delta[:, eligible_indices] / robust_scale
    local_index = {original: local for local, original in enumerate(eligible_indices)}
    score_matrix = np.column_stack([
        np.nanmean(scaled[:, [local_index[gene_to_index[gene]] for gene in genes]], axis=1)
        for _, genes in selected_sets
    ])
    observed = score_matrix.mean(axis=0)
    rng = np.random.default_rng(seed)
    exceed = np.zeros(len(selected_sets), dtype=int)
    batch_size = 1000
    completed = 0
    while completed < permutations:
        batch = min(batch_size, permutations - completed)
        signs = rng.choice(np.asarray([-1.0, 1.0]), size=(batch, len(patient_ids)))
        permuted = signs @ score_matrix / len(patient_ids)
        exceed += (np.abs(permuted) >= np.abs(observed)[None, :]).sum(axis=0)
        completed += batch
    p_values = (exceed + 1) / (permutations + 1)
    pathway_table = pd.DataFrame({
        "pathway_id": [gene_set.stable_id for gene_set, _ in selected_sets],
        "pathway_name": [gene_set.name for gene_set, _ in selected_sets],
        "complete_gene_count": [len(genes) for _, genes in selected_sets],
        "mean_patient_score": observed,
        "sign_flip_p_value": p_values,
        "genes": [";".join(genes) for _, genes in selected_sets],
    })
    pathway_table["bh_q_value"] = benjamini_hochberg(pathway_table["sign_flip_p_value"].to_numpy())
    score_table = pd.DataFrame(score_matrix, columns=[gene_set.stable_id for gene_set, _ in selected_sets])
    score_table.insert(0, "patient_id", patient_ids)
    return pathway_table.sort_values(["bh_q_value", "sign_flip_p_value"]).reset_index(drop=True), score_table

