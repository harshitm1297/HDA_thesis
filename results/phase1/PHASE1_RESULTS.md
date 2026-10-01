# Phase 1 QC and missingness results

## Scope and frozen decisions

Phase 1 characterises observable data quality, multivariate structure, pair concordance, and missingness before any differential-abundance result is calculated. The primary cohort remains all 42 patient pairs. A separate sensitivity cohort removes an entire pair when either specimen is flagged by at least two conceptually different diagnostic families. No single unusual PCA position, coverage value, or correlation is sufficient for removal.

The QC feature set was fixed without looking at differential results: a feature had to be observed in at least 80% of the 84 specimens. This admitted 3,759 rows, all of which had non-zero robust variability. For diagnostic PCA only, missing entries were replaced by the feature median without using tissue labels, features were median/MAD scaled, robust z-scores were clipped at ±5, and singular-value decomposition was applied. This temporary matrix is not exported as an inferential or predictive dataset.

Specimen flags use the same absolute robust-z threshold of 3.5 in both tissue classes, calibrated separately within `N` and `T`. This prevents a dataset-wide tissue difference from being mislabeled as a technical outlier. The five families are coverage, observed-value distribution, sample correlation, matched-pair concordance, and multivariate structure. The multivariate score uses a deterministic multi-start C-step approximation to minimum-covariance-determinant distance on the compact QC summaries. Isolation forest was not added because it would supply another correlated ranking without resolving the absent batch and laboratory metadata.

## Specimen coverage and observed-value distributions

Matched non-tumour specimens contain a median of 4,759 detected features, ranging from 2,460 to 6,549. Tumour specimens contain a median of 5,919.5, ranging from 4,463 to 6,370. Within pairs, tumour coverage exceeds matched non-tumour coverage for 37 of 42 patients. The paired coverage difference has a median of 1,023 features and a mean of 1,141, with a range from −1,262 to 3,467.

Median observed log2 abundance moves in the opposite direction: the median across matched non-tumour specimens is 13.084, versus 12.552 across tumour specimens. Coverage and observed median abundance have a pooled Spearman correlation of −0.914; the within-tissue correlations are −0.971 for matched non-tumour and −0.739 for tumour. This is strong evidence that measurement availability depends on observed abundance structure. It is not proof of a particular laboratory missingness mechanism, because acquisition settings, run order, normalization history, and batch labels are unavailable. Missing values therefore remain a separate detection process rather than being silently replaced. This treatment follows the cautions in [Webb-Robertson et al. (2015)](https://doi.org/10.1021/pr501138h), [Lazar et al. (2016)](https://doi.org/10.1021/acs.jproteome.5b00981), and [Kong et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/36349819/).

## Pair concordance and multivariate structure

The median detection-set Jaccard similarity within a patient pair is 0.708, with a range from 0.371 to 0.880. Among features observed in both specimens, the median absolute log2 difference per pair is 0.734, with a range from 0.301 to 1.521. The median matched-pair Spearman correlation is 0.745; individual pairs range from 0.367 to 0.940. These are descriptive concordance measures and do not establish whether a low-correlation pair is technically defective or biologically divergent.

The label-blind diagnostic PCA explains 20.07% of robust-scaled variance on PC1 and 9.89% on PC2; the first five components together explain 46.75%. Tissue centroids differ strongly on both leading components. The clustering order was calculated from pairwise-complete Spearman correlation before tissue annotations were attached. Both results demonstrate large systematic tissue-associated structure, but they cannot distinguish tumour biology from tissue-correlated technical effects because the dataset supplies no batch, run-order, purity, inflammation, or sampling-distance metadata.

## Missingness intelligence

Missingness falls sharply with observed abundance. The lowest observed-abundance decile has mean missingness of 78.64%, whereas the highest decile has mean missingness of 4.45%. Tumour detection exceeds matched non-tumour detection for 5,413 feature rows; the reverse occurs for 986 rows, and 1,672 rows have equal tissue detection fractions. The tumour-minus-non-tumour detection difference is largest on average in the lower-middle abundance deciles rather than at the best-observed end.

Exact paired binary p-values were calculated from tumour-only and non-tumour-only pairs for every feature with discordance. There are 3,317 rows with an unadjusted descriptive value below 0.05. This number is not a discovery result: no feature is called significant in Phase 1, and multiplicity correction is deliberately deferred to the locked detection-testing family in Phase 3. Benjamini–Hochberg false-discovery control will be required before inferential interpretation ([Benjamini and Hochberg, 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)).

## Sensitivity flags and cohorts

Twelve specimens from seven patients meet the two-family flag rule. The flagged pairs are `P1`, `P2`, `P19`, `P21`, `P25`, `P31`, and `P42`. Pair-level concordance and robust multivariate structure jointly flag both specimens for several pairs; `P25T` is instead unusual within the tumour class for coverage, distribution, and sample correlation. `P31N` has the broadest agreement across diagnostic families. These are sensitivity flags, not technical-failure labels.

The frozen cohorts are therefore:

- Primary cohort: all 42 patient pairs and all 84 specimens.
- Sensitivity cohort: 35 complete patient pairs after removing the seven flagged pairs in their entirety.

This conservative construction prevents broken pairing and avoids choosing exclusions after viewing protein-level differential results. Phase 2 and Phase 3 must run the primary cohort first and use the 35-pair cohort only to assess robustness.

## Deliverables and milestone decision

The implementation writes specimen QC, the full feature missingness atlas, abundance-decile summaries, the 84×84 pairwise-complete Spearman matrix, label-blind clustering order, frozen pair-level cohort membership, PCA and clustering objects, four figures, machine-readable validation, and output hashes. Eight automated tests covering Phase 0 and Phase 1 pass.

Milestone M1 is **passed**. The all-pair and sensitivity cohorts, QC metrics, diagnostic feature rule, transformation used only for QC, and flag thresholds are now frozen before differential testing. M1 does not authorize sample deletion, biomarker claims, or imputation. The next phase must benchmark preprocessing and missing-data strategies without optimizing for the number of significant proteins or classifier performance.

