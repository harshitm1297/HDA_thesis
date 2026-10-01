# Thesis extension: data-science and AI plan

> Superseded: this thesis-relative plan is retained only as a historical record. The project now starts from first principles under `FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md`; no previous analysis, threshold, exclusion or conclusion is used as a baseline.

## Purpose and evidence cutoff

This document defines the work that is genuinely new relative to the existing MSc thesis, using methodological research available up to and including 31 December 2023. It is the governing analytical roadmap for the data-only continuation of the oral-cancer proteomics project. Later publications and later regulatory developments are not used to justify the statistical or machine-learning methods in this plan.

“Data-only” has a strict meaning here. The only patient-level observations, labels and quantitative measurements used in any calculation will come from `OC_Dataset_84Sample_v2_24012024.xlsx`. Publications available by the cutoff may justify a method and published gene-set definitions may be used for clearly secondary interpretation, but literature cannot supply missing cohort attributes, laboratory metadata, outcomes, labels, validation samples or replacement measurements. No public cohort is part of this roadmap.

The comparison below is grounded in a direct review of `MSc_Project__Sem_4___Updated_.pdf` and in the reproducible structural audit of `OC_Dataset_84Sample_v2_24012024.xlsx`. The audit establishes 8,071 source rows, 84 specimens, 42 complete tumour/non-tumour patient pairs and 35.456% missing abundance cells. These are dataset-derived facts; the publications cited later provide methodological support rather than additional project observations.

The extension is not a rerun of the thesis with additional plots. Its contribution is a reproducible, missingness-aware, multiplicity-controlled and leakage-resistant high-dimensional analysis, followed by interpretable machine learning and stability assessment. The patient, rather than the specimen, remains the independent unit throughout.

## Existing thesis baseline

The thesis already performed substantial exploratory work. It described protein-detection counts, applied a 70% detection cutoff, separated tumour and non-tumour groups, min-max scaled them, selected condition-specific K values for K-nearest-neighbour imputation, log2 transformed the completed matrix, examined correlations and PCA, removed three patient pairs identified through correlation/PCA, applied a nominal paired-test threshold, added fold-change cutoffs, and produced candidate lists and a volcano plot.

These steps provide a baseline but leave several high-dimensional-analysis problems unresolved:

- the report calls the 84 specimen columns 84 patients even though the design contains 42 paired patients;
- `PG.Genes` rows are described interchangeably as genes and proteins, although the source is a protein-group/gene-labelled quantitative matrix;
- KNN imputation is performed separately by known tissue class, allowing the true outcome to influence the completed predictor matrix;
- the missing-at-random assumption is asserted from presumed laboratory uniformity but cannot be verified from the matrix;
- the 70% cutoff and three sample exclusions are not evaluated through systematic sensitivity analysis;
- PCA and correlation are used to exclude biologically unusual observations without a documented technical-failure criterion;
- thousands of protein tests use a nominal p-value threshold without false-discovery-rate control;
- the displayed test formula should be checked against a conventional paired analysis of within-patient differences;
- candidate ranking depends largely on p-value and fold-change cutoffs, without confidence intervals, usable-pair counts, directional consistency or pipeline stability;
- no pathway-level analysis, resampling-based feature stability, formal prediction experiment, calibration analysis or external computational replication is reported;
- the conclusion mentions a Wilcoxon signed-rank test, but the thesis does not provide a corresponding methods-and-results section.

The new work therefore begins from the original workbook, not from the thesis's imputed 3,611-feature table or its 47-protein endpoint.

## Direct comparison: thesis versus extension

| Analytical question | Existing thesis | Data-only extension | Why the change matters |
| --- | --- | --- | --- |
| What is independent? | Describes 84 columns as patients | Treats them as 84 specimens from 42 paired patients | Prevents pseudoreplication and keeps both specimens from a patient together in validation |
| What is a feature? | Uses gene/protein language interchangeably | Preserves each source row as a protein-group/gene-labelled quantitative feature with a stable row key | Avoids silently inventing protein, isoform or gene-level certainty |
| How is missingness handled? | A 70% filter followed by tumour-specific and non-tumour-specific KNN imputation | Models detection patterns, makes unimputed paired abundance analysis primary, and limits imputation to labelled sensitivity analyses | Avoids making the completed matrix depend on the known class and exposes potentially informative non-detection |
| How are unusual samples handled? | Correlation/PCA screening leads to removal of P19, P20 and P23 pairs | Uses several prespecified QC diagnostics; retains all structurally valid pairs in the primary analysis and treats exclusions as sensitivities | Separates technical suspicion from biological extremeness and quantifies influence |
| How is differential abundance tested? | Nominal p-values and fold-change cutoffs after a procedure labelled paired | Explicit patient-blocked model, moderated uncertainty, effect confidence intervals and Benjamini-Hochberg FDR | At 3,611 tests, p < 0.05 alone would yield about 181 false positives on average under a complete null, so multiplicity control is essential |
| How are candidates ranked? | Thresholded p-value/fold-change lists | Combined evidence from effect size, uncertainty, FDR, usable pairs, paired direction, detection pattern and sensitivity stability | Produces an auditable evidence ranking rather than one cutoff-dependent list |
| Is biological organisation assessed? | Mainly protein-by-protein interpretation | Secondary ranked pathway/leading-edge analysis using versioned, pre-2024 research resources | Tests whether weak but coordinated signals agree at a systems level without changing the primary feature-level findings |
| Is AI/ML evaluated? | No formal predictive pipeline | Three feature views, two model families, nested patient-grouped validation, permutation testing and selection stability | Converts AI from a label into a falsifiable internal experiment with leakage controls |
| Is generalisation established? | Candidate claims extend toward diagnosis/prognosis | Only internal discrimination is estimated; diagnosis, prognosis and external validity remain prohibited | Aligns claims with the actual tumour-versus-matched-non-tumour data |

The novelty claimed in this project is **relative to the submitted thesis**, not the invention of a globally new algorithm. The scientific contribution is the integration of paired inference, missingness intelligence, sensitivity analysis and leakage-safe sparse learning into one reproducible evidence framework for this dataset.

## Appreciable data-science and AI contribution

The project will contain a substantial data-science component through reproducible data engineering, quantitative quality control, missingness-pattern analysis, competing preprocessing pipelines, high-dimensional paired inference, sensitivity matrices, ranked pathway learning and patient-level heterogeneity analysis.

The AI/ML component will be a complete prediction experiment rather than a decorative algorithm. It will compare abundance-only, detection-only and combined feature views; use elastic-net logistic regression as the primary sparse model and a linear support-vector machine as a comparator; keep all learned preprocessing and feature selection inside nested patient-grouped resampling; measure selection and coefficient stability; compare the fitted models with simple baselines; and test performance against within-pair label permutations. The output is an internally validated tissue-state model and a stability-ranked panel, not a clinical diagnostic system.

Deep neural networks, autoencoders, transformers and large nonlinear ensembles are excluded. With 42 independent patients and more than 8,000 candidate rows, these methods would add complexity without credible generalisation evidence. The methodological sophistication comes from correct validation, uncertainty and stability—not model size.

The planned technical workload is appreciable and separable into four data-science products: a reproducible data model and validation layer; a multivariate QC and missingness-analysis layer; a high-dimensional inferential layer with multiplicity and sensitivity control; and a leakage-safe predictive-learning layer. The AI component is therefore not being measured by the number of algorithms. It is measured by whether the entire learned pipeline is evaluated on unseen patients, beats simple and permuted-label baselines, and yields features whose selection is stable across resamples.

## Phase 0 - thesis audit and reproducible data foundation

The extension first establishes exactly which thesis operations are being retained, corrected or replaced. The immutable workbook is fingerprinted, sample names are parsed into 42 pairs, every source feature receives a stable row key, missingness and paired-detection counts are calculated, and ambiguous identifiers are preserved instead of silently merged.

This phase is new because the thesis does not provide a complete executable provenance chain or a row-level identifier audit. It also corrects the unit-of-analysis error by treating 42 patients, not 84 specimens, as independent.

Deliverables:

- immutable input checksum and reproducible import script;
- 84-specimen/42-patient manifest;
- 8,071-row feature and identifier audit;
- machine-readable data dictionary and validation report;
- thesis-method comparison table and limitations register.

Milestone M0 is complete when every source sample and row is traceable and all thesis-derived assumptions are separated from workbook-derived facts. The existing Phase 1 audit already satisfies the structural portion of this milestone and becomes the input to Phase 1 below.

## Phase 1 - multivariate QC and missingness intelligence

This phase replaces ad hoc outlier deletion with a multivariate QC system. It calculates per-sample coverage, missing fraction, log-scale distribution summaries, robust correlations, hierarchical structure, PCA coordinates, within-patient distances and robust anomaly scores. An optional isolation forest may be applied only to a compact set of QC summaries, not to the full 8,071-dimensional matrix. It can flag unusual samples for sensitivity analysis but cannot prove technical failure.

Missingness is treated as information rather than an inconvenience to erase. Every feature is represented by tumour and non-tumour detection rates and the four paired states: detected in both, tumour only, non-tumour only and missing in both. Missingness dependence on observed abundance, sample coverage and tissue state is quantified. This directly challenges the thesis assumption that missing values are random.

No patient is removed solely for PCA position or low average correlation. The primary cohort contains all structurally valid pairs; any flagged-patient exclusion is a named sensitivity dataset. This makes it possible to measure whether the thesis's removal of P19, P20 and P23 materially changes the conclusions.

Deliverables:

- sample and feature QC report;
- missingness map and paired-detection report;
- predeclared flag table with concordance across diagnostics;
- primary 42-pair cohort and sensitivity cohorts, including a thesis-comparable 39-pair cohort;
- explicit comparison of thesis and new QC decisions.

Milestone M1 is reached when flags and cohort definitions are frozen before differential results are examined. These outputs determine the preprocessing candidates evaluated in Phase 2.

## Phase 2 - preprocessing benchmark and uncertainty specification

This phase replaces a single irreversible preprocessing chain with a prespecified benchmark. Positive values will be evaluated on a log2 scale. The main comparison is log2 with no additional normalisation versus log2 with sample-median centring. The choice is made using distributional alignment, relation to sample coverage, preservation of pair structure and avoidance of erasing a possible global tissue shift.

Unimputed analysis is primary. Proteomics research before 2024 shows that missingness has multiple mechanisms and that no single imputation method is universally superior; in some classification settings, no imputation performs best. Class-wise KNN completion is therefore not used for primary inference. If KNN, low-value or other imputation is evaluated, it is a labelled sensitivity pipeline; in machine learning it is fitted only from the training fold without using held-out labels. [Webb-Robertson et al., 2015](https://pubmed.ncbi.nlm.nih.gov/25855118/), [Lazar et al., 2016](https://doi.org/10.1021/acs.jproteome.5b00981), and [O'Brien et al., 2018](https://doi.org/10.1214/18-AOAS1144) support this approach.

Feature eligibility is defined by paired information rather than separate class-wise availability. At least 30 complete pairs is the primary abundance threshold and at least 21 is a sensitivity threshold. Features dominated by differential detection are retained for detection analysis rather than forced into an abundance model.

Deliverables:

- preprocessing benchmark with objective comparison criteria;
- locked primary representation and named sensitivity pipelines;
- feature eligibility table with complete-pair and detection information;
- comparison with the thesis's 70% cutoff, min-max scaling and class-wise KNN results.

Milestone M2 is reached when transformations, normalisation, completeness rules and missingness handling are frozen without reference to protein significance. The locked representations are passed to Phase 3.

## Phase 3 - robust paired differential inference

This phase replaces nominal per-protein testing with high-dimensional paired inference. The primary model is `abundance ~ patient + tissue_status` on the locked representation. Empirical-Bayes variance moderation stabilises feature-wise uncertainty when sample size is modest, and Benjamini-Hochberg adjustment controls the false discovery rate. Limma's modelling and empirical-Bayes framework were established well before 2023 and support blocked or paired designs. [Ritchie et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4402510/).

Differential detection is analysed separately with an exact paired test. A robust proteomics model based on the MSqRob family is a sensitivity analysis rather than a mandatory replacement, because the supplied table lacks peptide-level evidence. [Goeminne et al., 2018](https://pubmed.ncbi.nlm.nih.gov/28391044/) provides the pre-2023 proteomics modelling basis.

Each feature receives an effect estimate, confidence interval, raw and adjusted p-values, usable-pair count, detection rates, direction-consistency count, identifier status and stability across preprocessing and sample-set scenarios. A multiverse-style sensitivity matrix shows which candidates survive reasonable choices. The thesis's 47-protein list is re-evaluated rather than accepted or discarded wholesale.

Deliverables:

- complete differential-abundance and differential-detection tables;
- FDR-controlled volcano and MA plots;
- paired patient-level plots for leading candidates;
- thesis-list replication table showing retained, weakened and unsupported candidates;
- robustness-ranked candidate set with explicit limitation IDs.

Milestone M3 is reached when candidates meet effect, uncertainty, multiplicity, observation and stability criteria. The full ranked statistic—not only the significant subset—is passed to Phase 4.

## Phase 4 - systems-level interpretation and data-driven heterogeneity

This phase adds biological organisation absent from the thesis. Ranked gene-set enrichment uses the complete ordered statistic, avoiding an arbitrary significant-protein cutoff. Any Reactome, Gene Ontology or curated-complex definition is treated as a versioned pre-2024 research resource used **after** the workbook-derived statistics have been calculated; it does not add patients, measurements, clinical labels or outcomes, and no primary conclusion may depend on it. The tested-feature universe is explicit, and ambiguous source identifiers are not expanded as though they were separately measured proteins. GSEA provides the ranked-set framework, and the 2022 Reactome resource is within the evidence cutoff. [Subramanian et al., 2005](https://pubmed.ncbi.nlm.nih.gov/16199517/) and [Reactome 2022](https://academic.oup.com/nar/article/50/D1/D687/6426058).

Patient-level paired-change profiles are explored using robust PCA, module scores and resampling-based clustering stability. Because stage, site, HPV, purity and exposure data are unavailable, clusters are called data-derived response patterns rather than clinical or molecular subtypes. Immune, stromal, metabolic or cytoskeletal interpretations remain hypotheses when tissue composition cannot be measured.

Deliverables:

- ranked pathway and leading-edge results;
- protein-pathway network with evidence and identifier status;
- patient paired-change map and cluster-stability assessment;
- interpretation ledger separating observation, literature-supported interpretation and speculation.

Milestone M4 is reached when pathway results are reproducible across reasonable preprocessing scenarios and unstable heterogeneity patterns are excluded from substantive claims. Module scores and robust candidate information are passed to Phase 5, but predictive feature selection is still repeated independently inside resampling.

## Phase 5 - leakage-safe AI and stable panel discovery

This phase is the principal AI contribution. Three prediction views are compared:

1. abundance-only features from sufficiently observed rows;
2. detection-only binary features;
3. a combined view or low-dimensional pathway/module representation.

Elastic-net logistic regression is the primary model because it performs shrinkage and feature selection when predictors greatly outnumber observations and can retain groups of correlated variables. A linear SVM is the comparator. [Zou and Hastie, 2005](https://doi.org/10.1111/j.1467-9868.2005.00503.x).

The resampling unit is patient. Both samples from an outer-test patient remain unseen during filtering, normalisation, imputation, scaling, feature selection and tuning. Repeated six-fold outer cross-validation, with seven patients per outer fold, is proposed; an inner patient-grouped cross-validation selects hyperparameters. Leave-one-patient-out analysis is a sensitivity check. The 2023 `nestedcv` work directly supports nested feature selection for high-dimensional `p >> n` data. [Lewis et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37113250/).

Performance is compared with intercept-only and simple coverage/PCA baselines. Metrics include balanced accuracy, ROC AUC, sensitivity, specificity, Brier score, calibration intercept and slope where estimable, and distributions across repeats rather than a single value. A within-patient label-permutation experiment estimates whether the entire pipeline performs beyond chance without breaking the paired structure.

Feature selection frequency, coefficient-sign stability and panel-size stability are recorded across outer folds and repeats. Stability selection supplies the pre-2023 conceptual basis for treating repeated selection as evidence rather than trusting one fitted list. [Meinshausen and Buhlmann, 2010](https://people.math.ethz.ch/~peterbu/Files/Manuscripts/stabilityselection.pdf). Model explanation is centred on coefficients, selection frequency and fold-wise permutation importance; SHAP is secondary and will not be interpreted causally.

Deliverables:

- fully nested patient-grouped ML pipeline;
- comparison of abundance, detection and combined feature views;
- elastic-net and linear-SVM performance distributions;
- paired-label permutation null distribution;
- feature and coefficient stability report;
- compact exploratory panel with an explicit internal-validation label;
- TRIPOD-style model report and PROBAST risk-of-bias self-assessment, based on the 2015 and 2019 guidance available by the end of 2023: [TRIPOD](https://doi.org/10.1136/bmj.g7594) and [PROBAST](https://pubmed.ncbi.nlm.nih.gov/30596875/).

Milestone M5 is reached when all preprocessing and selection are verified to occur inside resampling, performance uncertainty is reported, and the panel is stable enough to describe as an internal candidate panel. No diagnostic or early-detection claim is permitted. Model outputs and stable candidates pass to Phase 6.

## Phase 6 - internal stress-testing, reproducibility and final synthesis

The final phase uses only the supplied workbook to determine how much of the result is robust and how much is contingent on analytical choices. It reruns the locked analysis across the declared completeness thresholds, normalisation candidates, all-pair and flagged-pair sensitivity cohorts, robust and classical estimators, and reasonable random seeds. Leave-one-patient-out influence analysis identifies candidates or model performance driven by individual patients. Negative controls include paired-label permutation for the ML pipeline and, where applicable, deliberately irrelevant or shuffled feature structures to confirm that the workflow does not manufacture signal.

This phase also performs an executable reproducibility check from the immutable workbook: regenerate all tables and figures in a clean environment, verify row counts and hashes, compare outputs with expected checksums or tolerances, and build a machine-readable claim ledger. External validation is recorded as unavailable by design, not replaced by an external public dataset or by literature-derived numbers.

Deliverables:

- sensitivity grid showing which findings survive each defensible analysis choice;
- leave-one-patient-out influence report;
- negative-control and paired-label-permutation results;
- final thesis-versus-extension comparison;
- integrated candidate evidence table containing statistical, biological and ML stability evidence;
- one-command or one-workflow regeneration of complete code, environment, hashes, tables and figures;
- final limitations and claim table, including the unmet need for genuinely independent validation;
- future experimental-validation specification, clearly separated from completed work.

Milestone M6 is reached when every final claim can be traced to a source row, an analysis version, a sensitivity result and applicable limitation IDs, and when no result depends on an external cohort or on invented project metadata.

## Minimum new contribution required for success

The extension will count as appreciably new only if it produces all of the following:

- a reproducible audit that corrects the 42-patient/84-specimen structure;
- missingness-aware abundance and detection analyses without class-wise primary imputation;
- FDR-controlled paired modelling with effects and confidence intervals;
- a robustness matrix across preprocessing and sample-inclusion scenarios;
- ranked pathway analysis and leading-edge interpretation;
- a fully nested patient-grouped AI pipeline with two model families and three feature views;
- permutation-based chance comparison and fold-wise feature-stability analysis;
- an explicit comparison showing which thesis candidates reproduce under the improved analysis;
- an internal stress-test matrix and leave-one-patient-out influence analysis;
- a final limitations-aware evidence table rather than a single unqualified biomarker list.

This combination provides a substantive data-science and AI extension while remaining credible for a 42-patient, high-dimensional dataset.

## Immediate next milestone

The structural audit is already complete. The next work is Phase 1 multivariate QC and missingness intelligence. Before computing sample rankings, the QC metrics, flag-combination rules, primary all-pair cohort and thesis-comparable sensitivity cohort must be written into a locked Phase 1 analysis specification.

