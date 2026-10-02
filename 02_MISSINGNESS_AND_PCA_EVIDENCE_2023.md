# Phase 2 missingness and PCA evidence lock through 2023

## Scope

This note freezes the pre-2024 evidence constraints that will govern Phase 2. It is tailored to the supplied protein-group/gene-labelled abundance matrix: 84 specimens from 42 matched tumour/non-tumour patients, 8,071 feature rows, and 35.456% missing abundance cells. No external study supplies replacement values, labels, batches, or biological outcomes. Published work is used only to define plausible methods, diagnostics, and failure modes.

## Missingness model appropriate to this dataset

The project will not assume that all missing entries arise from one mechanism. Label-free proteomics can contain random failures, sample- or feature-associated missingness, and intensity-dependent non-detection. Lazar et al. showed that the appropriate imputation strategy depends on the missingness mechanism and processing stage rather than on a universally best algorithm ([Lazar et al., 2016](https://doi.org/10.1021/acs.jproteome.5b00981)). O'Brien et al. showed that nonignorable missingness can bias contrasts and that single imputation hides uncertainty and inestimable comparisons ([O'Brien et al., 2018](https://doi.org/10.1214/18-AOAS1144)).

The most directly relevant result available by 2023 is the detection-probability work of Li et al. Their label-free proteomics analyses model detection probability as an intensity-dependent logistic function and describe the process as intermediate between simple MAR and hard censoring. They also show that common imputation strategies can reduce power or inflate false discoveries ([Li et al., 2023](https://doi.org/10.1093/bioinformatics/btad200)). This aligns closely with the current dataset, where missingness falls from 78.6% in the lowest observed-abundance decile to 4.4% in the highest.

Accordingly, Phase 2 will retain two complementary data views:

- observed log2 abundance for quantitative paired contrasts;
- finite-value detection for paired detection contrasts.

Blank cells will never be interpreted as measured zeros. Imputation will not be mandatory for the primary paired abundance analysis if the selected model can use available complete pairs. It will be used only where a complete matrix is mathematically required and as an explicit sensitivity analysis.

## Candidate methods that may enter the benchmark

The benchmark will contain restrained representatives of different assumptions rather than many superficially different algorithms.

### No-imputation reference

Available-pair quantitative analysis and the separate detection analysis are the reference. This prevents imputation from becoming an unquestioned source of pseudo-observations. Features with insufficient complete pairs remain eligible for detection analysis rather than being filled solely to force them into an abundance model.

### Local-similarity methods

K-nearest-neighbour or local least-squares imputation may represent locally predictable MAR-like structure. Earlier label-free proteomics evaluations found local-similarity approaches competitive, but no method was uniformly best ([Webb-Robertson et al., 2015](https://pubmed.ncbi.nlm.nih.gov/25855118/)). A large label-free comparison found random-forest imputation strong under its benchmark conditions, while emphasizing that performance depended more on the MNAR proportion than the total missing fraction ([Jin et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33469060/)). These findings justify testing local or nonlinear prediction, not selecting it in advance.

For this `p >> n` dataset, neighbour counts, distance metrics, and feature scaling must be tuned without tissue labels and evaluated for stability. A random-forest method is optional because 84 specimens provide limited information for thousands of feature-wise predictions and could overfit the observed matrix.

### Low-rank and PCA-based methods

Iterative SVD, NIPALS, probabilistic PCA or Bayesian PCA may represent global low-rank MAR-like structure. Comparative work through 2022 evaluated PPCA, NIPALS, SVD and singular-value-thresholding alongside local and censored approaches and found performance to depend on the simulated mechanism and evaluation design ([Wang et al., 2022](https://doi.org/10.1038/s41598-022-05227-0)). `pcaMethods` provides established incomplete-data PCA implementations, including NIPALS, BPCA, PPCA and SVD-impute ([Stacklies et al., 2007](https://doi.org/10.1093/bioinformatics/btm069)).

Low-rank approaches are plausible because patient proteomes contain correlated programmes, but they can oversmooth patient-specific biology, strengthen dominant tissue structure, and manufacture correlations. Component rank must therefore be selected through masked-data cross-validation, never through visual tissue separation.

### Left-censored methods

A left-censored strategy such as QRILC or a restrained low-intensity draw may represent intensity-dependent non-detection. It will not be applied to every blank, because the 2023 detection-probability evidence rejects a simple hard threshold and the current data cannot identify each cell's mechanism. Minimum, half-minimum, zero, or blanket down-shifted Gaussian replacement will not serve as the primary method because these substitutions compress uncertainty and can create artificial group differences.

### Detection-probability modelling

An intensity-dependent detection-probability curve is the most mechanism-aware option supported by 2023 work. The project will estimate an empirical detection-versus-observed-abundance relationship and assess whether a logistic curve is adequate. A full likelihood model will be used only if it can be justified for the supplied feature-level table; Li et al.'s method was developed principally for peptide/precursor-level proteomics, whereas this workbook contains already aggregated protein-group/gene-labelled rows.

## Benchmark design

The project has no complete hidden truth, so evaluation will combine artificial masking of observed entries with stability on the authentic missingness pattern. Pure uniform random masking is insufficient.

Artificial masks will be generated repeatedly at the patient-pair level under at least three regimes:

1. MCAR-style masking sampled uniformly from observed cells.
2. Intensity-dependent masking calibrated to the empirical detection curve.
3. A mixture of random and intensity-dependent masking, stratified by tissue and feature-abundance decile.

Masking will preserve the 42-patient structure and will be repeated over fixed seeds. Evaluation will include log2-scale MAE/RMSE, bias by abundance decile and tissue, error in paired `T−N` differences, variance compression, correlation distortion, detection-pattern preservation where applicable, and stability of downstream ranked effects. No method will be selected solely from entry-wise reconstruction error.

The primary choice will require acceptable performance across multiple masking regimes and must not create implausibly narrow distributions or tissue differences. If methods trade reconstruction accuracy against paired-effect preservation, the paired scientific estimand takes priority.

## PCA design for sparse high-dimensional proteomics

Ordinary PCA requires a complete matrix and can be strongly affected by imputation. Missing-aware PCA methods have therefore been developed for omics data. NIPALS, BPCA and PPCA can estimate components with incomplete data ([Stacklies et al., 2007](https://doi.org/10.1093/bioinformatics/btm069)), and sequential projection-pursuit PCA was designed specifically for non-randomly missing label-free mass-spectrometry data ([Nelson et al., 2018](https://doi.org/10.1186/s12859-018-2338-7)). Proteomics work has also demonstrated missing-value-tolerant NIPALS PCA in incomplete proteomic matrices ([Neumann et al., 2022](https://doi.org/10.1038/s41467-022-31007-x)).

Phase 2 will compare, at minimum:

- complete-feature SVD PCA as a low-missingness reference;
- the Phase 1 high-coverage median-imputed robust-scaled PCA;
- a missing-aware NIPALS or probabilistic PCA;
- an iterative low-rank PCA/SVD reconstruction selected by masked-data cross-validation.

The comparison will vary prespecified observation thresholds such as 70%, 80%, 90% and 100%, while reporting how many features and abundance ranges each threshold retains. Center-only log2 PCA will be compared with robust autoscaling, because autoscaling can give low-variance noisy proteins the same influence as stable high-information proteins. Extreme-value winsorisation will remain a declared diagnostic sensitivity, not an invisible preprocessing step.

High dimensionality is addressed computationally through thin SVD or an equivalent dual/sample-space solution, since 84 samples constrain the maximum centered rank to 83. Statistical stability will be assessed by resampling whole patients, not individual specimens. For each representation, the project will record:

- PC-score and loading stability under patient bootstrap and leave-one-pair-out analysis;
- principal-subspace angles or Procrustes agreement across imputations and filters;
- explained-variance stability;
- matched-pair neighbourhood preservation;
- association of PCs with sample coverage, missing fraction and median observed abundance;
- whether apparent tissue separation persists across no-imputation and missing-aware alternatives.

PCA remains exploratory. A clean tumour/non-tumour plot is not a selection criterion, and no PC will be called biological when it is primarily associated with coverage. Because batch, run order, purity and tissue composition are unavailable, stable tissue separation can be described only as tissue-associated structure, not as a pure tumour programme.

## Prediction-specific constraint

For later machine learning, filtering, normalization, imputation, scaling, PCA and feature selection must be fitted inside each patient-grouped training fold. The outer-test patients cannot influence neighbours, feature medians, component loadings, retained-feature rules or hyperparameters. The Phase 1 PCA is therefore a cohort-level QC visualization and must never be reused as a fitted transformer in Phase 5.

## Locked decision rule

Phase 2 will not choose an imputer or PCA representation because it produces more differentially abundant proteins, clearer visual separation, or better non-nested classifier performance. It will choose the smallest defensible set of primary and sensitivity representations based on missingness mechanism diagnostics, masked-value recovery, preservation of paired contrasts, multivariate stability, and explicit failure analysis. If no imputation method is reliable across plausible mechanisms, the primary abundance analysis will remain available-pair based and imputed results will be reported only as sensitivity evidence.

