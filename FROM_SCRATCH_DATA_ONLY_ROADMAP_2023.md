# From-scratch oral-cancer proteomics roadmap

## Governing scope

This is the governing scientific roadmap for analysing `OC_Dataset_84Sample_v2_24012024.xlsx` from first principles. It does not use any previous analysis, candidate list, threshold, sample exclusion or conclusion as a baseline. The empirical evidence is restricted to the supplied workbook. Research published by 31 December 2023 may justify methods and support clearly labelled interpretation, but it may not provide replacement measurements, cohort characteristics, clinical outcomes, laboratory metadata or validation samples.

The executable algorithms, schemas, model grids, resampling design and validation tests are specified in `IMPLEMENTATION_SPEC_DATA_SCIENCE_AI_2023.md`.

The workbook contains 8,071 protein-group/gene-labelled feature rows measured across 84 specimens representing 42 complete tumour/matched-non-tumour patient pairs. There are 240,380 missing abundance cells out of 677,964, or 35.456%. The biological unit is therefore the patient, not the individual specimen.

The source field `PG.Genes` is used as an annotation label attached to a quantitative proteomics feature. It is not gene-expression data and does not establish that every row is one uniquely identified protein or isoform. This project intersects gene-centred bioinformatics and proteomics by using the supplied gene labels to organise protein-group evidence, assess pathways and communicate results. It is not a transcriptomic–proteomic multi-omics study because no transcriptomic measurements are present.

## Scientific objectives

The primary objective is to identify reproducible within-patient differences between tumour and matched non-tumour tissue at the quantitative-feature and detection levels. Reproducibility here means stability across reasonable preprocessing choices, sample influence checks and resampling—not external validation.

The second objective is to determine whether the feature-level changes form coherent biological programmes when interpreted through published pathway knowledge available by 2023. Pathway results remain secondary: the primary statistical evidence must first be calculated entirely from the workbook.

The third objective is to characterise heterogeneity in each patient's tumour-minus-non-tumour change profile without assigning clinical subtypes that the available metadata cannot support.

The fourth objective is to test whether tumour and matched non-tumour specimens can be discriminated internally using a leakage-safe, patient-grouped machine-learning pipeline. This is an exploratory tissue-state classification experiment. It is not a screening, diagnostic, prognostic or treatment-response model.

## Analytical principles

Pairing is preserved in every inferential and validation step. A test, resampling split or sensitivity analysis that treats the 84 specimens as independent is invalid for this project. Both specimens belonging to an outer-test patient must be absent from every operation fitted on training data.

Missingness is both a measurement problem and a possible source of biological signal. Blank cells are not automatically zero and are not assumed to be missing at random. Detection and quantitative abundance will be analysed as complementary views. Proteomics research available before 2024 shows that missingness mechanisms and imputation performance vary by dataset, so no single imputation method can be accepted without diagnostics and sensitivity analysis. [Webb-Robertson et al., 2015](https://doi.org/10.1021/pr501138h), [Lazar et al., 2016](https://doi.org/10.1021/acs.jproteome.5b00981), [O'Brien et al., 2018](https://doi.org/10.1214/18-AOAS1144), and [Kong et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36349819/).

All high-dimensional inference must control multiplicity and report effects and uncertainty. A raw p-value cutoff is not a candidate-discovery system. Benjamini–Hochberg false-discovery-rate control and empirical-Bayes variance moderation supply the core pre-2024 framework. [Benjamini and Hochberg, 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) and [Ritchie et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4402510/).

Prediction is kept separate from explanation. Differential abundance asks which features differ within paired tissues; machine learning asks how well a pipeline discriminates tissue state in patients not used to fit it. A protein can be inferentially strong but predictively redundant, or predictively useful but difficult to interpret. The two evidence streams will be integrated only after both analyses are complete.

## Phase 0 — research contract and immutable data model

Phase 0 converts the workbook into a trustworthy analytical object before any biological result is examined. The source file is fingerprinted; sample names are parsed into patient and tissue identifiers; every source row receives a stable row key; numeric coercion, impossible values and missing cells are checked; duplicate and compound labels are flagged; and a data dictionary defines every derived field. No duplicate identifier is averaged and no compound label is split without evidence that the rows represent equivalent measurements.

The scientific contract fixes the estimands and prohibited claims. The primary abundance estimand is the average within-patient tumour-minus-non-tumour difference on the selected log scale among features meeting a prespecified paired-observation rule. The detection estimand is the paired difference in detection probability. The prediction target is the tissue state of a specimen from a patient excluded from training. Prognosis, recurrence, treatment response, stage-specific effects and population screening are outside scope because their required outcomes or comparison groups are absent.

Deliverables are the immutable-input checksum, 84-specimen/42-patient manifest, 8,071-row feature audit, data dictionary, limitations register, computational environment file and written analysis specification. Milestone M0 is passed only when every sample and feature can be traced to the source and the primary estimands are frozen. The existing structural audit may satisfy this milestone because it was derived directly from the workbook, but no previous analytical decisions are inherited.

## Phase 1 — multivariate quality control and missingness intelligence

Phase 1 describes the observable structure without deleting unusual biology. Per-specimen summaries include detected-feature count, missing fraction, positive-value distribution, median and robust spread on the log scale. Multivariate diagnostics include robust sample correlations, hierarchical clustering, PCA on eligible features, within-patient distances and anomaly scores calculated from a compact QC summary matrix. Isolation forest may be used as an exploratory anomaly detector on those summaries, but it cannot label a specimen as a technical failure.

Missingness is mapped at feature, specimen, tissue and patient-pair levels. Each feature receives counts for detected in both tissues, tumour only, non-tumour only and neither. Paired discordance is tested descriptively and later inferentially. Relationships between detection, observed abundance, specimen coverage and tissue state are examined to determine whether simple random-missingness assumptions are implausible.

No structurally valid patient is removed from the primary cohort solely because of PCA position, clustering or an anomaly score. Prespecified concordance across diagnostics creates a flag, and flagged-patient exclusions define sensitivity cohorts. Deliverables are a QC report, missingness atlas, patient-pair distance report, feature-observation table and frozen flag table. Milestone M1 is passed when the all-pair primary cohort and every sensitivity cohort are locked before differential results are inspected.

## Phase 2 — preprocessing and missing-data benchmark

Phase 2 selects a restrained primary representation using criteria that do not depend on which proteins become significant. Positive measurements are evaluated on a log2 scale. The principal normalisation comparison is log2 without additional centring versus log2 followed by sample-median centring. Diagnostics consider distribution alignment, association with sample coverage, preservation of within-patient contrast and the risk of erasing a real global tissue shift.

The primary abundance analysis should avoid imputation where the chosen model can use available paired observations. Features are admitted according to complete-pair counts, with a stricter threshold as the primary rule and a broader threshold as sensitivity. Features with strong tissue-specific detection but insufficient paired quantitative observations remain eligible for detection analysis rather than being forced into an abundance model.

Imputation is a sensitivity experiment, not a hidden preprocessing default. Candidate strategies may include KNN for locally structured random missingness, low-rank reconstruction for global structure and a left-censored strategy for abundance-dependent non-detection. Artificial masking is stratified by observed abundance and missingness pattern rather than applied only to completely observed features. Evaluation covers reconstruction error, distribution distortion, preservation of pair differences and downstream candidate stability. For predictive modelling, every imputer is fitted independently inside each training fold.

Deliverables are a preprocessing benchmark, locked primary representation, named sensitivity representations, feature-eligibility table and machine-readable parameter file. Milestone M2 is passed only when transformation, normalisation, completeness and missingness rules are fixed without reference to significance or classifier performance.

## Phase 3 — paired differential abundance and detection

Phase 3 estimates tumour-associated change while respecting pairing and high dimensionality. The primary abundance model is a patient-blocked linear model of the form `abundance ~ patient + tissue_status`, with empirical-Bayes variance moderation. An equivalent paired-difference formulation is retained as a transparent check. For every eligible feature the result contains the effect estimate, confidence interval, raw p-value, Benjamini–Hochberg adjusted value, number of usable pairs and direction-consistency count.

Detection differences are analysed separately from observed abundance. For each feature, the discordant pair counts support an exact paired binary test such as McNemar's exact test. Features can therefore be abundance-dominant, detection-dominant, concordant across both views or unresolved. A robust estimator and a missingness-aware or hurdle-style method may be used as sensitivity analyses if they are valid for the available protein-group-level table.

Candidate status is determined by a transparent evidence table rather than a single composite score. Evidence dimensions are effect magnitude, interval precision, FDR, pair support, direction consistency, detection evidence, identifier ambiguity and stability across the Phase 2 representations and Phase 1 cohorts. Deliverables are complete abundance and detection tables, paired plots, FDR-controlled volcano and mean–difference plots, sensitivity matrices and a tiered candidate table. Milestone M3 is passed when candidate tiers can be regenerated from explicit rules and every candidate carries its observation and limitation metadata.

## Phase 4 — gene-centred systems bioinformatics and heterogeneity

Phase 4 moves from individual source rows to coordinated biological patterns without changing the primary evidence. Ranked gene-set enrichment uses the complete signed Phase 3 statistic instead of only a significant subset. Gene-set definitions from Reactome, Gene Ontology or curated complexes must have a recorded release and must predate or be available by the 2023 cutoff. The measured-feature universe is the background, duplicate labels are handled explicitly, and compound protein-group labels are not expanded into multiple independent measurements. The ranked-set approach is supported by [Subramanian et al., 2005](https://pubmed.ncbi.nlm.nih.gov/16199517/), with the versioned pathway resource described by [Reactome 2022](https://academic.oup.com/nar/article/50/D1/D687/6426058).

Patient heterogeneity is evaluated on tumour-minus-non-tumour change profiles, which remove stable between-person baseline differences. Robust PCA and pathway/module scores provide lower-dimensional representations. Consensus clustering is attempted only on this restrained representation and is accepted only if cluster assignments are stable under patient resampling and reasonable preprocessing choices. Without stage, site, HPV, histology, purity or outcome metadata, any stable groups are described as computational response patterns rather than clinical or molecular subtypes.

Deliverables are ranked pathway results, leading-edge feature tables, an evidence-labelled feature–pathway network, patient change maps and cluster-stability diagnostics. Milestone M4 is passed when pathway conclusions survive reasonable sensitivity analyses and every interpretive statement is labelled as direct observation, research-supported interpretation or speculation.

## Phase 5 — leakage-safe AI and stable feature-panel discovery

Phase 5 is a prespecified internal prediction study. Three feature views are compared: quantitative abundance, binary detection and a combined or pathway-score view. Elastic-net logistic regression is the primary learner because it regularises a `p >> n` problem and can retain groups of correlated predictors. A linear support-vector machine is the main comparator. A simple coverage-only model and a low-dimensional PCA or pathway-score model are baselines. Elastic net is supported by [Zou and Hastie, 2005](https://doi.org/10.1111/j.1467-9868.2005.00503.x).

Validation is repeated nested cross-validation grouped by patient. A six-fold outer design gives seven patients—and therefore seven tumour and seven matched non-tumour specimens—in each outer test fold. Inner grouped folds tune hyperparameters. Filtering, normalisation, imputation, scaling, feature selection and pathway-score construction that learns from the cohort all occur inside the outer training data. This design follows the high-dimensional nested-validation principles described by [Lewis et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37113250/).

Performance is reported as distributions and confidence intervals for balanced accuracy, ROC AUC, sensitivity, specificity, Brier score and calibration measures where sample size permits. Repeated within-patient label permutation reruns the complete pipeline and establishes a chance distribution. Leave-one-patient-out validation is a sensitivity analysis, not the only performance estimate.

Interpretability focuses on fold-wise coefficients, selection frequency, coefficient-sign consistency and permutation importance calculated only on held-out data. Stability selection motivates requiring repeated selection rather than trusting one fitted panel. [Meinshausen and Bühlmann, 2010](https://doi.org/10.1111/j.1467-9868.2010.00740.x). SHAP may be included only as a secondary explanation of a validated fold-level model and never as causal evidence.

Deliverables are the nested pipeline, performance and calibration distributions, baseline comparisons, permutation-null report, learning curves, model-card-style documentation and a stability-ranked exploratory panel. Milestone M5 is passed only if an automated leakage audit confirms that no outer-test information entered training and the model is described as internally validated. TRIPOD and PROBAST provide the pre-2024 reporting and bias-assessment framework. [TRIPOD, 2015](https://pubmed.ncbi.nlm.nih.gov/25563062/) and [PROBAST, 2019](https://pubmed.ncbi.nlm.nih.gov/30596875/).

## Phase 6 — integration, stress testing and reproducible synthesis

Phase 6 tests whether the scientific conclusions survive plausible analytical perturbations. A multiverse grid varies the declared completeness rule, normalisation, missing-data strategy, robust versus classical estimator, flagged-patient cohort and random seed. Leave-one-patient-out influence analysis identifies findings driven by one pair. Negative controls include paired-label permutation and shuffled or deliberately irrelevant feature structures where appropriate.

Inference, systems interpretation and predictive evidence are integrated without allowing one to substitute for another. The final candidate table records source-row key, supplied label, quantitative effect, confidence interval, FDR, detection evidence, usable pairs, directional consistency, preprocessing stability, patient-influence stability, pathway membership, ML selection frequency, limitation identifiers and permitted claim. A feature is called a high-priority internal candidate only when multiple evidence dimensions agree; it is not called a biomarker validated for clinical use.

The complete project is regenerated from the immutable workbook in a clean environment. Scripts verify input and output hashes, dimensions, pair structure, random seeds, software versions and expected numerical tolerances. FAIR and MIAPE principles motivate the provenance and reporting structure. [Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18) and [Taylor et al., 2007](https://doi.org/10.1038/nbt1329).

Deliverables are the robustness grid, patient-influence report, negative-control results, integrated evidence table, reproducible code and environment, final figures, claim ledger and limitations-aware scientific report. Milestone M6 is passed when every claim can be traced to workbook rows, an analysis version, sensitivity evidence and applicable limitations.

## Phase hand-offs

| From | Required hand-off | To |
| --- | --- | --- |
| M0 | Validated pair manifest, stable feature keys, estimands and limitations | Phase 1 QC |
| M1 | Locked primary and sensitivity cohorts plus missingness descriptors | Phase 2 preprocessing |
| M2 | Locked representations and feature-eligibility rules | Phase 3 inference |
| M3 | Complete ranked statistics and tiered feature evidence | Phase 4 systems analysis |
| M4 | Stable pathway/module summaries and heterogeneity assessment | Phase 5 predictive modelling |
| M5 | Out-of-fold predictions, null results and feature stability | Phase 6 synthesis |
| M6 | Fully traceable internal evidence package | Future independent validation, outside current scope |

## Minimum success criteria

The project succeeds as a rigorous data-science and AI study only if it delivers all of the following:

- a fully reproducible patient-paired data model and audit;
- an explicit missingness atlas and separate detection analysis;
- a locked preprocessing benchmark that does not select choices using significant results;
- FDR-controlled paired inference with effects, confidence intervals and usable-pair counts;
- sensitivity and single-patient influence evidence for every leading candidate;
- versioned pathway analysis and a stability-tested heterogeneity assessment;
- a patient-grouped nested ML experiment with preprocessing inside folds;
- simple-model and paired-permutation baselines;
- feature-selection and coefficient stability rather than a one-fit panel;
- an integrated evidence and limitations table that prohibits unsupported clinical claims.

## Immediate milestone

The directly derived structural audit already provides the raw material for M0. The next action is to write and lock the Phase 1 QC and missingness analysis specification before calculating sample anomaly rankings or looking for tumour-associated candidates. That specification must define the exact QC metrics, transformation used only for visual diagnostics, flag-concordance rule, primary all-pair cohort, sensitivity cohorts and required figures.
