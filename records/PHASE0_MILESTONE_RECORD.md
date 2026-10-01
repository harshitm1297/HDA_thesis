# Phase 0 milestone record

## Record date

1 October 2026

## Questions addressed

The Phase 0 discussion established the scientific claim, identified AI and data-science methods that fit the sample size and paired design, reconstructed the best available experimental provenance, and defined how unresolved laboratory information will be handled.

## Decisions

The primary scientific question is paired tumor-versus-matched-non-tumor protein abundance and detection. Predictive modelling is secondary. The project will use the term “matched non-tumor tissue,” not “healthy control,” because distance, pathology, and field effects are unknown.

The supplied matrix will be described as positive, linear-scale, protein-group/gene-labelled abundance values from 42 matched pairs. A Spectronaut-derived, label-free DIA workflow is the leading hypothesis, but acquisition mode, quantity field, upstream normalization, protein-inference settings, identification FDR, and instrument details remain unverified.

Project-specific experimental facts will not be replaced with typical settings from unrelated publications. Instead, uncertainty is recorded with confidence levels, and the downstream analysis will use conservative defaults that remain valid under incomplete provenance.

The primary missingness strategy will preserve missing values, report complete-pair and detection counts, and separate differential abundance from differential detection. Class-wise KNN imputation will not define the primary dataset.

The primary predictive model, if reached, will be elastic-net logistic regression under patient-grouped nested resampling. Deep learning and high-capacity nonlinear ensembles are outside the primary scope.

The European alignment review retained “oral cancer” as the supported cohort term and reserved OCSCC for pathology-confirmed oral-cavity cases. Oral cavity/mobile tongue and oropharynx/base of tongue will remain distinct, and any recovered stage will carry its cTNM or pTNM components, assessment date, and TNM edition. This prevents the 2026 TNM 9 terminology from being imposed retrospectively on an incompletely documented 2024 cohort.

The provenance inventory now includes ethics, consent or another lawful basis, controller roles, permitted secondary use, retention, access, transfer and deposition conditions. Patient codes are treated as pseudonymous pending formal assessment. EHDS is used to shape a future-compatible dataset catalogue, not as an automatic permission to disclose data.

The classifier remains research-only. A later diagnostic intended purpose would require a separate assessment under the IVDR/MDR and, where applicable, the EU AI Act. The current work will preserve versions, intended use, provenance, validation boundaries and limitations without implying clinical or regulatory validation.

## Phase 0 artifacts

- `THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md`: governing redesigned project roadmap and AI/data-science scope
- `PHASE0_RECONSTRUCTED_PROVENANCE.md`: verified facts, evidence-based inferences, unresolved information, source checksums, and operational defaults
- `records/PHASE0_DETAILED_RESPONSE_RECORD.md`: complete user-visible response archive
- `records/PHASE0_DETAILED_RESPONSE_RECORD.md`: detailed user-visible response archive
- `RESEARCH_EVIDENCE_BASE.md`: phase-by-phase supporting literature and data resources
- `EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md`: European clinical vocabulary, data governance, AI/IVDR boundary, and implementation fields current to 1 October 2026

## Completion status

Phase 0 is operationally complete for a transparent secondary analysis. Full experimental reproducibility remains impossible without original raw files, a complete Spectronaut export, an analysis project, or laboratory records. Any recovered material will be incorporated as a versioned provenance update rather than silently replacing this record.

## Transition to Phase 1

Phase 1 will create a reproducible import, sample manifest, feature-level identifier audit, data dictionary, and validation report. The original workbook will remain unchanged.

## Data-only operating decision — 1 October 2026

The project will proceed without expecting further information from the supervisor, laboratory or project group. Phase 0 therefore closes when uncertainty is bounded and propagated, not when full provenance is recovered. Missing raw files, acquisition settings, protein-inference evidence, batch labels, pathology, clinical covariates and governance records remain permanent limitations for this analysis. Their IDs, mitigations and prohibited claims are maintained in `../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`.

The governing sequence is now `../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md`. Phase 0 hands a frozen claim boundary and limitation register to the completed Phase 1 structural audit, which in turn hands stable sample and feature keys to Phase 2 matrix-only QC.

