# Phase 3 paired abundance and detection results

## Analytical status

The Phase 3 implementation is complete for the frozen no-imputation inference branch. Computational validation passes, but the full M3 milestone remains conditional until the wider Phase 2 imputation/PCA benchmark is completed. This distinction prevents the present results from being treated as the final preprocessing-invariant candidate set.

The primary analysis uses all 42 patient pairs, uncentred log2 observed abundances, no imputation, and a minimum of 30 complete pairs per feature. The tumour-minus-matched-non-tumour sign convention is used throughout. A feature is tested only when the source values support the required number of paired differences.

Two separate multiple-testing families are maintained. Quantitative abundance uses a paired-difference model with an inverse-chi-square empirical-Bayes variance prior estimated across eligible features. Detection uses the exact paired binary test based on tumour-only and non-tumour-only detections. Benjamini–Hochberg adjustment is applied separately to each family, consistent with the high-dimensional testing framework of [Benjamini and Hochberg (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). The variance-moderation design follows the empirical-Bayes principle used by limma, while this implementation uses a transparent dependency-light moment estimator rather than claiming byte-for-byte equivalence to limma ([Smyth, 2004](https://doi.org/10.2202/1544-6115.1027); [Ritchie et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4402510/)).

## Quantitative abundance

The primary ≥30-pair rule admits 3,430 of 8,071 source rows. The fitted variance prior has 3.776 prior degrees of freedom and scale 0.718 on the log2-difference variance scale. The output reports ordinary and moderated statistics, moderated confidence intervals, raw and BH-adjusted p-values, pair support, and patient-direction consistency for every source row.

At BH q≤0.05, 2,194 eligible rows show quantitative evidence: 986 have positive tumour-minus-non-tumour effects and 1,208 have negative effects. Applying the prespecified practical filters—absolute mean log2 difference ≥1, direction consistency ≥70%, and ≥30 complete pairs—leaves 669 primary abundance-core rows. Requiring at least 80% sign agreement and BH significance across the eleven sensitivity runs leaves 614 unambiguous stable abundance candidates: 128 positive and 486 negative.

The most statistically precise stable abundance rows include the source labels `PRELP`, `SERPINH1`, `OGN`, `SLC3A2`, `HSPH1`, `ALDH9A1`, `EPHX1`, `PEBP1`, `SOD3`, and `DDAH2`. These names are annotations attached to quantitative source rows, not independently validated protein identities or clinical biomarkers. `PRELP`, for example, has a mean log2 tumour-minus-non-tumour difference of −5.346 with a moderated 95% interval from −6.052 to −4.640; `SERPINH1` has a mean difference of 2.425 with an interval from 2.091 to 2.760. Their interpretation remains internal to this cohort.

## Detection differences

The primary exact paired-detection family contains 2,887 rows with BH q≤0.05. Adding an absolute detection-fraction difference ≥0.20 and at least ten discordant pairs gives 2,549 primary detection-core rows. Requiring the same direction, effect and FDR criteria in the 35-pair unflagged sensitivity cohort leaves 1,831 stable detection patterns. After excluding ambiguous duplicate or compound identifiers, 1,827 rows receive the `C_detection_pattern` label.

Stable detection patterns are strongly asymmetric: 1,695 favour tumour detection and 136 favour matched non-tumour detection. Leading tumour-detection source labels include `SLC38A2`, `NEDD1`, `KPNA7`, `NCAPG`, `NOMO1`, `HMGA2`, `WDR75`, `SLC38A5`, `POP4`, and `SLC39A14`. This broad asymmetry must be interpreted alongside the Phase 1 finding that tumour specimens have higher overall coverage. The result is a reproducible tissue-associated detection pattern, not proof that every missing matched-non-tumour value represents biological absence.

## Integration and sensitivity

The analysis evaluates twelve prespecified combinations of representation, cohort and pair-support threshold. After excluding the primary combination, eleven runs form the abundance sensitivity matrix:

- uncentred versus sample-median-centred log2 abundance;
- all 42 pairs versus the 35-pair unflagged cohort;
- minimum complete-pair thresholds of 20, 30 and 35.

There are 614 stable `B_abundance` rows and 1,827 stable `C_detection_pattern` rows. No row meets the strict `A_concordant` definition. This is not an error: a feature measured quantitatively in at least 30 pairs has little room to accumulate the large number of discordant detections required for strong exact detection evidence. The two evidence types therefore identify largely complementary regions of the data.

The absence of an A-tier result prevents forced agreement between abundance and detection. Detection-only rows are explicitly labelled as patterns with unresolved quantitative abundance. Abundance rows are labelled as internal paired-abundance candidates whose detection evidence may be absent or discordant.

## Modularity

The implementation separates data access, statistical distributions and multiplicity correction, paired inference, evidence integration, visualisation, configuration and command-line orchestration. Statistical backends can therefore be replaced—for example with an R/limma backend—without changing source-row keys, cohort definitions, output schemas or evidence rules. All numerical thresholds reside in `config/phase3.yml` rather than being embedded throughout the analysis code.

The main result tables contain all 8,071 source rows, including ineligible and ambiguous rows. Missing statistical values mean that a row did not satisfy the declared testing rule; they are not silently converted to null results. Duplicate and compound labels cannot enter the highest evidence tier.

## Limitations and next dependency

The current results do not establish clinical biomarkers, causality, diagnosis, prognosis or external validity. Upstream normalization, batch, acquisition order, tumour purity and histopathology remain unknown. The predominance of negative abundance effects and positive tumour-detection effects may reflect genuine tissue biology, tissue composition, measurement coverage, or a mixture.

Before treating M3 as fully closed, Phase 2 must complete the mechanism-aware imputation and PCA benchmark recorded in `PHASE2_MISSINGNESS_AND_PCA_EVIDENCE_2023.md`. The present primary no-imputation branch will remain unchanged; Phase 2 can add sensitivity evidence but cannot be used to select a representation because it increases the discovery count.

