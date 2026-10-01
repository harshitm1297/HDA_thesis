# Thesis-extension and 2023 data-only response record

## Request recorded

The project was asked to identify what is genuinely new relative to the existing thesis, ensure that the continuation contains an appreciable amount of data science and AI, redesign phases and deliverables where necessary, use research available only through 2023, and make no analytical dependency on information outside the supplied dataset.

## Response recorded

The new contribution is not another round of correlation, PCA, KNN imputation and thresholded t-tests. It is a reproducible high-dimensional evidence pipeline that starts from the original workbook and explicitly treats the 84 columns as 42 paired patients. Relative to the thesis, it adds a source-row audit, missingness-pattern modelling, prespecified multivariate QC, patient-blocked moderated inference, false-discovery-rate control, effect confidence intervals, differential-detection analysis, sensitivity matrices, data-driven heterogeneity analysis, and a leakage-safe machine-learning experiment.

The thesis remains the baseline comparator. Its 70% filter, class-specific KNN imputation, PCA/correlation exclusions, 39-pair analysis, nominal p-value threshold and 47-candidate endpoint will be reconstructed wherever the available files allow. The extension will then show which findings remain when the analysis retains all structurally valid pairs, avoids outcome-aware primary imputation, uses patient-blocked models and controls multiplicity. At 3,611 simultaneous tests, a nominal 0.05 threshold would produce approximately 181 false positives on average if every null hypothesis were true and the tests were calibrated; this is why the revised project cannot treat raw p-values as sufficient evidence.

The appreciable AI component is a full internal prediction experiment, not the inclusion of a fashionable model. Abundance-only, detection-only and combined or module-level views will be compared. Elastic-net logistic regression is the primary model and a linear support-vector machine is the comparator. All filtering, normalisation, imputation, scaling, selection and tuning will occur inside nested cross-validation, and both specimens from a patient will remain in the same fold. Evaluation will include performance distributions, calibration where estimable, simple baselines, within-patient label permutation, feature-selection frequency and coefficient-sign stability. Deep learning is excluded because 42 independent patients cannot support credible estimation of a highly parameterised network.

Only the supplied workbook contributes patient-level observations, tissue labels or quantitative measurements. Literature published by 31 December 2023 is used to justify methods and to qualify interpretation; it cannot manufacture stage, site, HPV status, batch, instrument settings, outcomes or missing values. Published gene sets may be used only as a secondary, versioned interpretation layer after workbook-derived statistics have been calculated. No external cohort is imported, and external validation remains an unmet future requirement.

The redesigned phases are:

1. Phase 0 audits the thesis and establishes a reproducible, immutable source foundation.
2. Phase 1 performs multivariate QC and builds a paired missingness atlas without automatic outlier removal.
3. Phase 2 benchmarks restrained preprocessing choices and locks primary and sensitivity representations before testing outcomes.
4. Phase 3 performs paired differential-abundance and differential-detection inference with FDR, effects, confidence intervals and stability evidence.
5. Phase 4 adds secondary systems-level interpretation and tests whether patient-level paired-change patterns are stable, without naming clinical subtypes.
6. Phase 5 performs leakage-safe sparse AI modelling, patient-grouped nested validation, permutation controls and stable-panel discovery.
7. Phase 6 performs internal stress tests, leave-one-patient-out influence analysis, clean regeneration and final evidence synthesis using only the supplied workbook.

The governing specification, exact deliverables, milestone hand-offs, citations and thesis-comparison table are maintained in `../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md`. The methodological bibliography restricted to research through 2023 is in `../RESEARCH_EVIDENCE_BASE.md`. The continuing restrictions caused by absent project-group information are in `../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`.

## Decision status

The structural data audit is complete. The immediate analytical milestone is Phase 1 multivariate QC and missingness intelligence. Before result-driven inspection, the QC metrics, sample-flag rule, all-pair primary cohort and thesis-comparable sensitivity cohort must be locked in a written analysis specification.
