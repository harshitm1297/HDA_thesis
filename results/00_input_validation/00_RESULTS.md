# Phase 0 implementation results

## Decision and scope

Phase 0 is the first phase in the governing from-scratch roadmap, so it was implemented before multivariate QC. Its purpose is to prevent downstream analysis from silently changing the cohort, pairing, feature identity, numerical scale, or permitted scientific claims. The source workbook is treated as immutable and is identified by SHA-256 `65211B375E366D0CF6078E1016A3D29867EFC5E0EBA10D97F1F0282B25161DA8`.

The implementation uses the patient as the independent unit. Tumour and matched non-tumour specimens are parsed only from the anchored naming rule `^P(?P<patient>[0-9]+)(?P<tissue>[NT])$`. The abundance estimand is frozen as within-patient `log2(T) - log2(N)` for eligible complete pairs, while detection is retained as a separate paired binary process. This separation is necessary because proteomics missingness can mix abundance-dependent non-detection with other mechanisms; treating every blank as zero or applying one untested imputation strategy would manufacture quantitative information. This design is consistent with the missing-value cautions of [Lazar et al. (2016)](https://doi.org/10.1021/acs.jproteome.5b00981) and [Webb-Robertson et al. (2015)](https://doi.org/10.1021/pr501138h).

## What was implemented

The source validator checks the workbook name and hash, exact worksheet and identifier field, row and sample counts, uniqueness and syntax of sample IDs, completeness of all patient pairs, numerical coercion, and the presence of zero, negative, infinite, or missing identifiers. It assigns stable feature keys `F00001` through `F08071` without merging repeated gene labels or splitting semicolon-delimited protein groups.

The data-model builder emits a raw abundance matrix `X_raw` and a finite-value detection matrix `D`, each shaped 84 specimens by 8,071 feature rows. It also emits a 42-patient by 8,071-feature paired log2 tumour-minus-non-tumour matrix and a paired detection-state matrix. Detection states are explicitly encoded as `0=neither`, `1=N-only`, `2=T-only`, and `3=both`. Blank values remain `NaN`; there is no zero filling or imputation. Because all 437,584 observed values are finite and strictly positive, direct log2 transformation of observed values is numerically valid without a pseudocount.

The frozen contract prohibits claims about screening, diagnosis, prognosis, recurrence, treatment response, stage, grade, HPV, oral subsite, tumour purity, exposure effects, causality, or external biomarker validation. Those fields are absent from the supplied data and are not reconstructed from literature. FAIR-style traceability is implemented through input and output hashes, stable row keys, a data dictionary, and regenerable outputs, following the general reproducibility principles of [Wilkinson et al. (2016)](https://doi.org/10.1038/sdata.2016.18). Proteomics reporting discipline is informed by MIAPE ([Taylor et al., 2007](https://doi.org/10.1038/nbt1329)), while explicitly recording which acquisition metadata remain unavailable.

## Validation findings

The workbook contains 8,071 feature rows, 84 uniquely named specimens, and 42 complete `N/T` patient pairs. There are 437,584 numeric observations and 240,380 blanks, giving 35.456% missingness. No observed value is zero, negative, infinite, or non-numeric after coercion. Observed values range from approximately 0.002342 to 7.311 billion, supporting the contract's interpretation of the supplied values as positive linear-scale quantities rather than already log-transformed values.

The feature audit finds 114 rows with no measured abundance in any specimen. Eight rows carry duplicated labels: `CDKN2A`, `CUX1`, `MOCS2`, and `TMPO` each occur twice. Forty-eight rows contain semicolon-delimited compound labels. These findings are flags about identifier ambiguity, not permission to combine or expand rows. Every original row remains independently traceable.

Of the 8,071 feature rows, 7,239 have at least one patient with both tissues observed and can potentially contribute to a paired quantitative analysis; 832 have no complete quantitative pair. Only 1,350 are complete across all 42 pairs. Across all 338,982 patient-feature combinations, 176,223 are detected in both tissues, 18,605 only in matched non-tumour, 66,533 only in tumour, and 77,621 in neither. The large number of tissue-discordant detections is a direct reason to retain the separate detection estimand. It is not yet evidence of biological differential detection because Phase 1 QC and Phase 3 multiplicity-controlled inference have not been run.

## Verification and milestone decision

Four automated tests passed. They verify that the sample regular expression is anchored, the manifest contains exactly 42 complete pairs, every matrix has the specified orientation and dimensions, the paired-difference sign is truly tumour minus matched non-tumour, every validation check is true, and the scientific contract is frozen. The output manifest records the byte size and SHA-256 hash of every generated Phase 0 object.

Milestone M0 is therefore **passed**. This means the data foundation is reproducible and the estimands are fixed; it does not mean that any protein, pathway, outlier, or classifier has been validated. Phase 1 can now consume the sample manifest, feature manifest, raw matrix, detection matrix, and paired objects without reparsing or silently redefining the cohort.

## Immediate Phase 1 hand-off

Phase 1 should use all 42 pairs as the primary cohort and calculate specimen coverage, missing fraction, positive log2 distribution summaries, robust correlations, PCA/clustering on prespecified eligible features, and within-pair distances. Any anomaly method must rank or flag samples for sensitivity analysis only; it must not declare technical failure in the absence of batch and laboratory metadata. Missingness must be tabulated by specimen, feature, tissue, and pair state before any imputation or feature filtering is chosen.

