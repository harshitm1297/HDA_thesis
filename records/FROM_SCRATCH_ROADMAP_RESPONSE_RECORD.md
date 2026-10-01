# From-scratch roadmap response record

## Request recorded

The project was reset so that the existing thesis would not serve as a baseline or influence analytical decisions. A rigorous roadmap was requested at the intersection of gene-centred bioinformatics, quantitative proteomics, data science and AI, while retaining the established data-only scope and research cutoff of 31 December 2023.

## Response recorded

The project now begins with the supplied workbook and first principles. The measured entities are treated as protein-group/gene-labelled quantitative features, not gene-expression measurements. Gene symbols provide an annotation bridge to pathway analysis, but the project is not described as multi-omics because there is no transcriptomic layer.

Four scientific objectives govern the work: estimate reproducible within-patient tumour/non-tumour abundance and detection differences; determine whether changes organise into coherent pathways; characterise patient heterogeneity in paired change profiles; and evaluate tissue-state discrimination through internally validated machine learning. Prognosis, screening, recurrence, response and causal mechanism remain outside scope.

The seven connected phases are an immutable data model and research contract; multivariate QC and missingness intelligence; preprocessing and missing-data benchmarking; paired differential abundance and detection; systems bioinformatics and heterogeneity; leakage-safe AI and stable-panel discovery; and integrated stress testing with reproducible synthesis.

The central AI design compares abundance, detection and combined or pathway-score feature views. Elastic-net logistic regression is the primary sparse learner and a linear support-vector machine is the comparator. Repeated nested cross-validation is grouped by patient, and all learned preprocessing occurs inside training folds. Baselines, paired-label permutations, calibration, held-out importance and selection stability are mandatory. Deep learning is excluded because 42 independent patients cannot support a credible high-parameter model.

The current governing specification is `../FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md`. The older thesis-relative plan is retained only as a superseded historical record and no threshold, exclusion, candidate or conclusion from it is inherited.

## Low-level implementation follow-up

The roadmap has been translated into `../IMPLEMENTATION_SPEC_DATA_SCIENCE_AI_2023.md`. That specification defines the Python/R division of work, repository structure, configuration contracts, sample and feature schemas, raw and detection matrices, paired-difference construction, QC statistics, missingness tests, preprocessing benchmark, limma model, exact detection test, pathway mapping rules, heterogeneity stability checks, nested patient-grouped machine learning, model grids, permutation procedure, leakage unit tests, multiverse analysis, leave-one-patient-out influence analysis and final evidence-table schema.

The AI pipeline uses 25 repeated six-fold outer patient partitions and five-fold inner grouped tuning. Elastic-net logistic regression is primary, linear SVM is the comparator, and detection, abundance and combined views are evaluated separately. All feature filtering, imputation, scaling, selection, calibration and tuning are fold-local. The specification also defines explicit failure tests for patient overlap, split pairs, preprocessing before splitting, full-data feature selection and outer-test-driven hyperparameter choice.

## Immediate milestone

The workbook-derived structural audit supplies the raw material for Phase 0 without importing previous analytical decisions. The next deliverable is a locked Phase 1 QC and missingness specification defining metrics, diagnostic transformation, sample-flag concordance, primary all-pair cohort, sensitivity cohorts and required outputs before anomaly rankings or tumour-associated candidates are inspected.
