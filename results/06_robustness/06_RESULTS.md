# Phase 6 robustness and integrated-evidence results

## Analytical status

Phase 6's computational work is complete. It stress-tests the Phase 3 inference, Phase 4 pathway interpretation and Phase 5 predictive analysis, and then integrates them without treating one evidence type as a substitute for another. All inputs come from the supplied 84-specimen workbook and its derived Phase 0–5 artifacts. No external patient measurements, labels or outcomes were introduced.

The software checks pass, but M6 is deliberately recorded as `IMPLEMENTED_NOT_CLOSED_PHASE2_AND_FULL_CLEAN_ENVIRONMENT_GATES_PENDING`. Phase 2 has not yet supplied the planned preprocessing and missing-data benchmark, so the multiverse can test only the scientifically defensible no-imputation branch. The pipeline also reproduced all 22 compared Phase 6 outputs exactly in the present environment, but it has not yet been rebuilt from the workbook in a newly created independent environment. These are milestone dependencies, not computational failures.

## Prespecified robustness multiverse

The grid contains 24 scenarios: two abundance representations (uncentred log2 and sample-median-centred log2), two cohorts (all 42 pairs and the 35-pair unflagged sensitivity cohort), three minimum complete-pair thresholds (20, 30 and 35), and two estimators (empirical-Bayes moderated and ordinary paired t statistics). Missingness is deliberately represented by one level, `no_imputation`; adding invented imputation variants before Phase 2 would create an appearance of robustness unsupported by the completed workflow.

Across the 8,071 source features, 3,077 have complete sign agreement wherever they are eligible, 1,757 receive FDR support in at least 80% of eligible scenarios, and 368 meet the complete frozen core rule in at least 80% of eligible scenarios. The number of core-supported features changes materially with analytical choice: the strict uncentred, unflagged, 35-pair scenario supports roughly 319 features, whereas the median-centred all-pair 20-pair scenarios support approximately 1,050. This variation is precisely the point of the multiverse: a single preferred analysis would hide sensitivity to representation, cohort and eligibility thresholds. Multiverse analysis is used here as a structured sensitivity analysis, following the principle that defensible analytical decisions should be displayed rather than collapsed into one researcher-selected result ([Steegen et al., 2016](https://doi.org/10.1177/1745691616658637)).

## Patient influence

The primary abundance and detection analyses were refitted 42 times, each time removing both specimens from one patient. This preserves the patient as the independent unit. Among the 614 Phase 3 abundance candidates, 583 retain core abundance support in at least 80% of omissions, and none reverses its estimated abundance direction. Among the 1,827 detection-pattern candidates, 1,804 retain core detection support in at least 80% of omissions. The detailed table retains the maximum effect displacement, sign flips and retention fractions for every feature so candidates influenced by particular patients remain visible rather than being silently discarded.

The leading abundance elastic-net pipeline was also refitted in an exact grouped leave-one-patient-out design. For every held-out pair, all abundance eligibility decisions, medians, scaling, supervised selection and hyperparameter tuning were learned from the other 41 patients with five inner patient-grouped folds. The resulting 84 held-out predictions yield ROC AUC 0.9847, balanced accuracy 0.9643 and correct tumour-versus-non-tumour orientation in 41 of 42 patient pairs. These estimates describe internal tissue-state separation only. Leave-one-out validation does not manufacture an independent population and cannot establish clinical discrimination.

## Pathway and model stability

Reactome release 86 remains frozen throughout the pathway sensitivity analysis. Of 783 tested pathways, 596 meet both prespecified Phase 6 criteria: direction agreement of at least 80% across the 24 scenarios and median leading-edge Jaccard overlap of at least 0.30. In total, 563 pathways retain identical direction in every scenario. The leading-edge overlap criterion matters because a pathway can preserve its name and enrichment direction while being supported by different measured proteins.

The Phase 5 comparison across views, learners, seeds and grouped fold assignments remains stable. Abundance elastic net has median repeated-nested-validation AUC 0.9824, abundance linear SVM 0.9796, the fold-local abundance PCA baseline 0.9694, combined-view models about 0.974, detection models 0.954–0.960, the coverage-only baseline 0.8685 and the intercept null 0.500. Phase 6 does not select a new winner from these results; it verifies that the main interpretation does not depend on one seed or partition. The small performance gap between the high-dimensional abundance model and the PCA baseline, and the dense median elastic-net model size of 97 from Phase 5, still prevent a compact protein-panel claim.

## Negative controls

The already frozen Phase 5 paired-label control reran the entire nested combined elastic-net pipeline for 200 within-patient label swaps. The observed AUC was 0.9768, the null median was 0.4786, the largest null AUC was 0.7319 and the finite-sample empirical p-value was 0.004975. This rejects paired-label exchangeability for that prespecified pipeline inside this dataset; it says nothing about transportability or cancer specificity.

A second Phase 6 control tests whether convergence between inferential and predictive feature sets is greater than expected from set sizes alone. The observed overlap between Phase 3 tiered features and stable Phase 5 ML features is 140. Across 1,000 seeded random feature-key permutations, the null median is 45, the maximum is 64 and the empirical p-value is 0.000999. The result supports non-random internal convergence, but correlated proteins, shared data and the absence of an external cohort mean that it is not independent biological replication.

## Integrated evidence package

The final join uses the immutable `feature_key`, never a possibly duplicated gene label. All 8,071 source rows remain present, and missing evidence remains missing rather than being converted to zero. The frozen rules assign 165 rows to `high_priority_internal`, 1,909 to `robust_single_domain`, 28 to `predictive_stability_only` and 5,969 to `not_prioritized`. The high-priority set comprises 107 abundance candidates and 58 detection-pattern candidates. Quantitative abundance directions among these rows are 75 positive and 89 negative; one detection-only row has no meaningful quantitative direction.

Frequently stable abundance examples include PRELP, SERPINH1, HSPH1, MAOB, KRT4, FRY, COLGALT1, OGN, SLC3A2, LIMA1, GPX3, DCN and FKBP9. These names are examples for follow-up, not a ranked clinical assay. A row reaches high internal priority only when its identifier is unambiguous, its relevant Phase 3 result survives the frozen multiverse or leave-one-patient-out rule, and it also has pathway-leading-edge or stable-ML support. The component columns remain visible so another analyst can reject the composite tier and inspect each evidence dimension separately.

![Phase 6 evidence tiers](figures/integrated_evidence_tiers.png)

## Traceability and interpretation boundary

The 14-row claim ledger connects each substantive conclusion to a saved evidence file, exact locator, evidence state, limitation IDs and allowed scope. The output manifest records hashes, the environment snapshot records Python 3.12.14 and the principal package versions, and the same-environment gate reproduced the compared artifacts exactly. This provenance structure is consistent with FAIR's emphasis on findability and reusability and MIAPE's minimum-information approach to proteomics reporting, while explicitly recording which acquisition metadata are unavailable ([Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18); [Taylor et al., 2007](https://doi.org/10.1038/nbt1329)).

The only supported final conclusion is that the supplied paired tissue proteomics matrix contains internally robust abundance, detection, pathway and predictive signals under the tested no-imputation analyses. Because the study has 42 patients, one laboratory-derived matrix, incomplete technical and clinical metadata, matched non-tumour rather than population controls, and no external or orthogonal validation, no output is an externally validated biomarker, clinical panel, diagnostic test, causal mechanism, prevalence estimate or treatment recommendation. This boundary is also consistent with the distinction between internal model development and clinical validation emphasized by TRIPOD and PROBAST ([Collins et al., 2015](https://doi.org/10.7326/M14-0697); [Wolff et al., 2019](https://doi.org/10.7326/M18-1376)).

## Artifact map

- `multiverse_scenarios.csv` and `multiverse_feature_stability.csv` expose scenario-level counts and feature-level stability.
- `lopo_feature_influence.csv` and `lopo_tiered_candidate_detail.csv.gz` retain patient-deletion influence summaries and full deletion results.
- `pathway_robustness.csv` and `pathway_sensitivity_detail.csv.gz` retain pathway direction and leading-edge stability.
- `ml_lopo_predictions.csv`, `ml_lopo_parameters.csv` and `ml_lopo_summary.json` retain exact grouped leave-one-patient-out predictions and tuning decisions.
- `negative_controls.csv` contains both empirical null experiments.
- `integrated_evidence_table.csv` retains all source features; `prioritized_internal_candidates.csv` is its filtered review view.
- `claim_ledger.csv`, `06_validation.json`, `reproducibility_gate.json`, `environment_snapshot.json` and `output_manifest.json` provide the audit layer.

