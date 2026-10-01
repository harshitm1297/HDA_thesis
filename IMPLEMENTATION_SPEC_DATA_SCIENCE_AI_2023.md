# Low-level implementation specification

## Relationship to the roadmap

This document implements `FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md`. The supplied workbook is the only source of patient-level observations, tissue labels and quantitative measurements. Method choices are restricted to research available by 31 December 2023. The implementation does not import previous analytical results or use an external cohort.

The software design separates four workflows:

1. **Data engineering** creates immutable, validated sample and feature records.
2. **Statistical inference** estimates paired abundance and detection differences.
3. **Systems bioinformatics** organises workbook-derived statistics using gene and pathway knowledge.
4. **Predictive AI** evaluates tissue-state discrimination using unseen patients.

The streams meet only in the final evidence table. A feature's inferential significance is never used to select predictors before a cross-validation split, and machine-learning importance is never treated as a p-value or causal effect.

## Repository and execution architecture

Use Python for workbook ingestion, data validation, exploratory QC, missingness summaries, resampling orchestration, machine learning and reproducibility tests. Use R/Bioconductor for empirical-Bayes differential analysis and ranked gene-set enrichment, where established domain-specific implementations are preferable. Exchange only versioned CSV or Parquet tables with explicit schemas; never pass an undocumented interactive object between languages.

The implementation should converge on this structure:

```text
config/
  analysis.yml
  feature_rules.yml
  ml_grid.yml
data/
  immutable/                 # source pointer and checksum, not a second edited workbook
  interim/                   # validated matrices and manifests
  processed/                 # locked representations
results/
  qc/
  missingness/
  preprocessing/
  inference/
  pathways/
  heterogeneity/
  ml/
  robustness/
reports/
scripts/
  00_validate_input.py
  01_build_data_model.py
  10_qc_missingness.py
  20_preprocessing_benchmark.py
  30_export_paired_matrices.py
  31_limma_paired.R
  32_detection_tests.py
  40_pathway_enrichment.R
  41_heterogeneity.py
  50_nested_ml.py
  51_permutation_ml.py
  60_robustness.py
  70_build_evidence_table.py
tests/
  test_data_contracts.py
  test_pairing.py
  test_no_leakage.py
  test_reproducibility.py
environment/
  python-lock.*
  renv.lock
```

`config/analysis.yml` is the single human-edited control file. It stores the input hash, sheet, identifier column, sample-name expression, random seeds, feature-eligibility candidates, normalisation candidates, outer and inner fold counts, repeat counts and output version. Analytical scripts read this file and must not contain hidden alternative thresholds.

Each table receives `analysis_version`, `input_sha256`, `created_utc` and, where relevant, `random_seed`. Every script stops on a failed assertion rather than continuing with a warning.

## Core data model

Construct these objects once and use them everywhere:

### Sample manifest

One row per specimen with at least:

```text
sample_id, patient_id, tissue, source_column,
source_column_index, pair_complete, detected_features,
missing_fraction, manifest_status
```

`tissue` is a two-level field derived only from the validated sample-name suffix. `patient_id` is the grouping variable for all resampling. Assert that there are exactly 84 unique samples, 42 unique patients, one `N` and one `T` specimen per patient, and no unmatched pair.

### Feature manifest

One row per source row with at least:

```text
feature_key, source_excel_row, supplied_label, label_type,
is_duplicate_label, is_compound_label, all_missing,
n_observed, n_missing, n_complete_pairs,
n_t_only, n_n_only, n_neither
```

`feature_key` is the primary key. `supplied_label` is annotation, not a key. Duplicate labels remain separate features. Compound labels remain one measured source feature and are excluded from any analysis requiring an unambiguous one-gene mapping unless a separate sensitivity rule is declared.

### Quantitative and detection matrices

Create `X_raw`, an 84 × 8,071 floating-point matrix indexed by `sample_id` and `feature_key`, with blanks represented only as `NaN`. Create `D = is_finite(X_raw)` as the binary detection matrix. Do not fill `NaN` with zero.

Create a paired representation for each feature and patient:

```text
delta[p, f] = log2(T[p, f]) - log2(N[p, f])
```

`delta` is missing unless both paired values are observed. Also create the paired detection state `both`, `T_only`, `N_only` or `neither`.

## Phase 0 implementation — validation and scientific contract

### Input validation

`00_validate_input.py` performs these steps in order:

1. Calculate SHA-256 before opening the workbook and compare it with the configured hash.
2. Require exactly one expected worksheet and identifier column.
3. Parse sample names with the anchored expression `^P(?P<patient>[0-9]+)(?P<tissue>[NT])$`.
4. Reject duplicate sample names, unparseable names, unmatched pairs or extra non-sample columns.
5. Convert abundance cells with strict numeric parsing. Store parsing failures with cell coordinates and stop.
6. Count zeros, negative values, infinities, blanks and numeric values. A direct log2 transformation is permitted only if every observed value is positive.
7. Assign a stable `feature_key` based on the source row, not the label.
8. Write manifests and validation JSON before producing any transformed matrix.

Acceptance tests compare the source and exported row/column counts, reconstruct sample order, prove unique keys and verify that the four paired detection states sum to 42 for every feature.

### Analysis contract

Write a machine-readable `analysis_contract.yml` containing:

- primary estimand: mean paired log2 tumour-minus-non-tumour difference;
- secondary estimand: paired difference in detection probability;
- independent unit: patient;
- prediction target: specimen tissue state in held-out patients;
- multiplicity families: abundance and detection tested separately;
- allowed claims and prohibited clinical claims;
- primary and sensitivity cohort definitions;
- frozen evidence cutoff and dataset-only boundary.

M0 passes when all validation tests succeed and the contract is committed before outcome-directed analysis.

## Phase 1 implementation — QC and missingness

### Univariate sample QC

`10_qc_missingness.py` computes, for each specimen:

- number and fraction of detected features;
- median, interquartile range and median absolute deviation of observed `log2` values;
- 1st, 5th, 95th and 99th percentiles;
- fraction of values outside robust median ± 5 MAD after log transformation;
- number of detected features shared with its matched specimen;
- Jaccard similarity of the paired detection vectors.

Do not calculate a pseudocount unless zeros exist. If zeros exist, report them and treat the log rule as a separate decision rather than silently adding a constant.

### Multivariate QC

Build a QC-eligible matrix without imputing the full dataset. Use features observed in a high fraction of specimens solely for correlation and PCA diagnostics. Calculate pairwise-complete Spearman correlation, hierarchical clustering with `1 - correlation` distance, and PCA after a declared temporary diagnostic imputation fitted without tissue labels. The imputed PCA is a visual diagnostic, not an inferential input.

For each specimen, create a compact vector containing coverage, distribution summaries, first robust-PC scores, median sample correlation and matched-pair distance. Robustly scale those variables with median and MAD. Calculate robust Mahalanobis distance using minimum covariance determinant where estimable. An isolation forest may be fit to the compact QC table with a fixed seed and contamination left unspecified or used only for ranking.

### Flagging rule

Create one Boolean flag per diagnostic rather than one unexplained outlier score. A sample enters a sensitivity flag set only if at least two conceptually different diagnostic families agree—for example, extreme coverage plus multivariate distance, or low correlation plus extreme distribution. The threshold for each diagnostic is stored in configuration and applied symmetrically to both tissues.

Do not remove flagged samples from the primary all-pair cohort. Generate named sensitivity cohorts by removing the entire patient pair, never one member of a pair.

### Missingness analysis

For each feature, compute tissue-specific detection proportions and the four paired states. The exact paired detection p-value is:

```text
k = n_T_only
n = n_T_only + n_N_only
p = 2 * min[P(Binomial(n, 0.5) <= k), P(Binomial(n, 0.5) >= k)]
```

Cap the two-sided value at 1.0 and return `NA` when `n = 0`. This calculation is descriptive in Phase 1; final FDR adjustment occurs over the locked detection-testing family in Phase 3.

Assess abundance-dependent missingness by calculating each feature's median observed log abundance, placing features into abundance deciles and comparing their missingness rates and paired detection states. At sample level, examine the relationship between coverage and observed-value distribution summaries. Do not infer a laboratory mechanism from these associations. Produce heatmaps ordered independently of tissue first, with tissue annotation added afterward, to avoid visually forcing separation.

M1 passes when all-pair and flagged-pair cohorts, QC metrics and flag thresholds are frozen before differential testing.

## Phase 2 implementation — preprocessing benchmark

### Candidate representations

Construct these representations without looking at feature-level tumour effects:

```text
R1 = log2 positive observed values, no additional normalisation
R2 = R1 minus each sample's observed-value median
R3 = sensitivity representation selected for a justified alternative normalisation
```

Compare R1 and R2 using sample-distribution spread, relation of median intensity to feature coverage, preservation of within-patient correlation, and the magnitude of the global paired shift. Do not choose the representation that produces the most discoveries or highest non-nested classifier accuracy.

### Feature eligibility

Generate an eligibility curve over candidate minimum complete-pair counts rather than choosing a convenient percentage. For each threshold, record the number of retained features, residual degrees of freedom, median standard error under a paired-difference model and coverage of the intensity distribution. Select the primary threshold using only observation counts and precision diagnostics; lock at least one broader and one stricter sensitivity threshold.

All-missing features are excluded from numerical analysis but retained in the manifest. Features lacking enough complete pairs for abundance analysis remain eligible for detection analysis if they have sufficient discordant pairs.

### Imputation benchmark

Primary inference uses available complete pairs without filling values. Benchmark imputers only for sensitivity analysis and ML requirements:

- training-feature median as the simple baseline;
- KNN using distances calculated from training data only;
- low-rank reconstruction with rank tuned inside training data;
- left-censored random draws only as a sensitivity scenario for abundance-dependent non-detection.

Create artificial missingness masks stratified by observed-abundance decile and feature completeness. Evaluate normalised root-mean-square error, median absolute error, distributional shift, pair-difference distortion and downstream rank stability. A method is not selected solely because it reconstructs randomly masked values well; authentic missingness may follow a different mechanism.

M2 passes when the primary representation, eligibility rule and all permitted sensitivity pipelines are frozen in configuration.

## Phase 3 implementation — paired inference

### Differential abundance

`30_export_paired_matrices.py` exports a feature × patient matrix of paired differences for every locked representation. `31_limma_paired.R` fits an intercept-only model to each difference row:

```r
design <- matrix(1, ncol(delta_matrix), 1)
fit <- limma::lmFit(delta_matrix, design)
fit <- limma::eBayes(fit, robust = TRUE)
result <- limma::topTable(fit, number = Inf, adjust.method = "BH",
                          sort.by = "none", confint = TRUE)
```

The coefficient is the mean paired log2 change. Confirm its sign convention with a synthetic two-patient fixture whose expected effect is known. Export effect, moderated standard error, confidence interval, t-statistic, raw p-value, BH value, usable-pair count and directional consistency. Directional consistency is the larger of the positive- and negative-delta counts divided by non-zero complete deltas; it is descriptive and not another significance test.

Run a conventional one-sample paired-difference t-test and Wilcoxon signed-rank test as transparent sensitivities, with their own adjusted values. Do not combine their p-values.

### Differential detection

`32_detection_tests.py` applies the exact discordant-pair binomial test to features meeting the locked discordance rule, then applies BH correction within that family. Export T-only and N-only counts, detection difference, odds-direction label, raw p-value and adjusted value.

### Candidate tiers

Do not calculate an opaque weighted score. Define rule-based tiers:

- **Tier A:** FDR-supported abundance or detection evidence, minimum pair support, consistent direction and stability across the primary cohort and declared preprocessing sensitivities.
- **Tier B:** strong effect and consistency but insufficient FDR or observation support; explicitly exploratory.
- **Tier C:** pipeline-dependent, identifier-ambiguous or patient-influence-sensitive evidence.

The exact numeric criteria are stored before inspecting candidate identities. Produce a long-form sensitivity table with one row per `feature_key × pipeline × cohort`.

M3 passes when every result is keyed to the source row and can be regenerated with identical signs, counts and adjusted values.

## Phase 4 implementation — systems bioinformatics

### Gene mapping contract

Keep the original feature table intact. For pathway analysis, construct a separate mapping table containing `feature_key`, supplied label, eligibility status, mapping rule and exclusion reason. The primary pathway analysis uses only unambiguous single-gene labels. Compound labels are excluded. Duplicate single-gene labels are resolved by a declared deterministic rule and tested with an alternative rule; they are never treated as multiple independent genes.

### Ranked enrichment

`40_pathway_enrichment.R` uses a signed moderated statistic from the primary abundance analysis. The gene universe is the set of features eligible for the same analysis, not every gene in a pathway database. Run preranked GSEA/`fgsea` with fixed random seed, versioned gene sets, minimum and maximum set sizes, BH adjustment and leading-edge export.

Repeat enrichment using at least one sensitivity statistic, such as effect signed by `-log10(p)` or the primary t-statistic under the stricter completeness representation. A pathway is stable only when its direction and leading interpretation agree across declared representations.

### Patient heterogeneity

`41_heterogeneity.py` uses the patient × feature delta matrix. Filter features using completeness and variability without using external clinical labels. Robustly scale features, fit PCA, and report explained variance with bootstrap intervals. Do not interpret UMAP or t-SNE islands as subtypes; if used, they are display-only and repeated across seeds.

For consensus clustering, evaluate a small prespecified range such as `k = 2..4`, resample patients, and calculate co-clustering probabilities, silhouette width and cluster-wise Jaccard stability. Reject subtype language unless stability is strong; even a stable solution is called a data-derived change pattern because clinical correlates are unavailable.

M4 passes when enrichment inputs are traceable and heterogeneity conclusions survive resampling and preprocessing sensitivity.

## Phase 5 implementation — predictive AI

### Prediction views

Build three separately evaluated model views:

1. **Detection view:** binary `D`; no abundance imputation.
2. **Abundance view:** log abundance with training-only filtering, imputation and scaling.
3. **Combined view:** concatenate controlled abundance and detection features, or combine their out-of-fold model scores using a meta-model trained only inside the inner loop.

Do not compute the combined view from Phase 3 significant features. Supervised feature filtering must be repeated inside each training fold. Fixed published pathway membership may be used, but any centring, scaling, PCA, module construction or feature selection learned from this cohort must also be fold-local.

### Repeated grouped nested cross-validation

`50_nested_ml.py` first generates and saves fold assignments. For each of 25 outer repeats, shuffle the 42 patient IDs using a recorded seed and divide them into six folds of seven patients. Each test fold therefore contains seven tumour and seven non-tumour specimens. Within the remaining 35 patients, use five grouped inner folds of seven patients for hyperparameter tuning.

The outer loop is:

```text
for repeat_seed in outer_seeds:
    partition patient IDs into six outer folds
    for outer_fold:
        define train_patients and test_patients
        fit every filter/transformer/selector only on train_patients
        tune hyperparameters by grouped inner CV
        refit the chosen complete pipeline on all outer-training patients
        predict probabilities/scores for outer-test patients once
        store predictions, chosen parameters and selected features
```

Save fold assignments before fitting. Assert that `train_patients ∩ test_patients = ∅`, both specimens of each patient share a fold, and no preprocessing object reports more fitted samples than the outer-training set.

### Models and tuning

The primary model is elastic-net logistic regression. Start with a prespecified grid:

```text
C: 10^-4, 10^-3, 10^-2, 10^-1, 1, 10, 100
l1_ratio: 0.1, 0.5, 0.9, 1.0
```

The comparator is a linear SVM with the same `C` grid. If probability estimates are required, calibrate within the inner training process; never calibrate on outer-test predictions. Avoid RBF kernels, random forests and boosting in the primary comparison because the patient count is too small for a broad model search. They may appear only as clearly secondary stress tests with the same nested design.

Tune using mean inner balanced accuracy, with a deterministic tie-break favouring the stronger regularisation and smaller selected panel. Record the full inner score surface rather than only the winning parameters.

### Fold-local processing

Within every outer and inner training split:

1. Remove features failing the training-only observation/prevalence rule.
2. Apply the configured log transformation to positive observed values.
3. Fit the imputer on training samples only.
4. Fit centring and scaling on training samples only.
5. Perform supervised filtering, if used, using training patients only.
6. Fit the classifier and choose hyperparameters.
7. Apply the frozen fitted operations to validation or test data.

For the linear SVM, use either embedded regularisation or a training-only paired univariate selector with `k` tuned in the inner loop. A selector fitted once on all 42 patients is prohibited.

### Baselines and metrics

Evaluate an intercept-only probability, a specimen-coverage model and a low-dimensional unsupervised baseline. A complex model must outperform these baselines to justify its complexity.

For each outer repeat, calculate ROC AUC, balanced accuracy, sensitivity, specificity, Brier score and log loss. Choose a probability or decision threshold inside the inner loop; do not optimise it on the outer test fold. Aggregate repeat-level metrics with medians, percentile intervals and full distributions. Calculate calibration intercept and slope from out-of-fold probabilities with an explicit small-sample warning.

### Permutation test

`51_permutation_ml.py` performs a paired label permutation. For every patient, flip the N/T labels together by randomly swapping the two labels with probability 0.5. Rerun the **entire** nested pipeline, including filtering and tuning. Use at least 200 permutations during development and target 1,000 for the final analysis if runtime allows. The empirical p-value is:

```text
p = (1 + number(null_metric >= observed_metric)) / (1 + number_of_permutations)
```

### Stability and interpretability

Across outer fits, record feature-selection frequency, median coefficient, coefficient-sign agreement and the distributions of selected panel sizes. Calculate permutation importance only on each outer-test fold, then aggregate it; never calculate importance on the training data and present it as generalisation evidence.

Define a stable exploratory feature as one exceeding a prespecified selection-frequency threshold and showing consistent coefficient direction. Report correlated feature groups so that elastic-net substitutions are not misread as biological contradictions. SHAP is optional and secondary; it does not replace fold-level performance or stability.

### Leakage unit tests

`test_no_leakage.py` must deliberately fail when:

- a patient appears in both outer train and test sets;
- paired specimens receive different folds;
- a scaler or imputer is fitted before splitting;
- full-data differential results are used as ML filters;
- a test label is available to a preprocessing function;
- pathway or PCA components learned on all samples enter a test fold;
- hyperparameters are selected using outer-test performance.

M5 passes only when these tests pass, the permutation null is reported and all performance values come from outer-test predictions.

## Phase 6 implementation — robustness and evidence integration

### Multiverse grid

`60_robustness.py` enumerates only defensible prespecified alternatives:

```text
cohort × completeness rule × normalisation × missingness strategy × estimator
```

For each feature, calculate effect-sign agreement, rank correlation, FDR-support frequency and maximum leave-one-patient-out effect change. For pathways, calculate direction agreement and leading-edge overlap. For ML, compare feature views, seeds, fold assignments and permitted imputers without selecting the best result after the fact.

### Leave-one-patient-out influence

Rerun the primary paired inference 42 times, omitting one patient each time. Record the range of effects, sign flips, maximum change in standardised effect and loss of candidate tier. For leading pathways and models, perform an analogous influence summary at the patient level.

### Integrated evidence table

`70_build_evidence_table.py` joins only by `feature_key` and produces:

```text
feature_key
supplied_label
identifier_status
n_complete_pairs
mean_log2_change
confidence_interval
abundance_fdr
detection_difference
detection_fdr
direction_consistency
preprocessing_stability
cohort_stability
leave_one_patient_out_stability
pathway_membership
ml_selection_frequency
ml_sign_stability
limitation_ids
allowed_claim
```

Missing evidence remains missing; it is never converted to a zero score. Candidate tiers are applied from the frozen rule file. Inference, pathway and prediction columns remain visible separately so readers can disagree with the final prioritisation.

### Reproducibility gate

Run the pipeline twice in a clean environment with the same configuration. Require exact equality for manifests, counts, fold assignments and discrete decisions. Require configured numerical tolerance for floating-point results. Record source and output hashes, Python and R package locks, random seeds, operating system and commands. Generate the final report from saved tables rather than copying values manually.

M6 passes when a clean run reproduces every final table and figure and every substantive claim links to a workbook-derived result, sensitivity output and limitation ID.

## Initial implementation sprint

The next coding sprint should stop before differential analysis. Its concrete sequence is:

1. Create `config/analysis.yml` and `analysis_contract.yml`.
2. Refactor the existing workbook audit into `00_validate_input.py` and `01_build_data_model.py` without changing validated source facts.
3. Materialise `X_raw`, `D`, the sample manifest, feature manifest and paired-state table.
4. Add unit tests for input hash, dimensions, sample parsing, exact pairing, stable row keys and detection-state sums.
5. Implement Phase 1 sample summaries, pair similarities, missingness tables and plots.
6. Write the QC flag rule in configuration before calculating the final anomaly rankings.
7. Produce a Phase 1 review report containing data contracts, QC findings, flags and locked cohorts.

Only after that report is frozen should Phase 2 preprocessing benchmarking begin.
