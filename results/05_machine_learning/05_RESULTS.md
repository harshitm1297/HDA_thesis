# Phase 5 leakage-safe AI and stability results

## Analytical status

Phase 5 completes the prespecified internal tumour-versus-matched-non-tumour tissue-state classification experiment. All reported performance comes from outer-test specimens belonging to patients absent from model fitting, preprocessing, feature selection, threshold choice and hyperparameter selection. Computational M5 passes, while its scientific status remains **conditional on Phase 2 and final upstream M3/M4 closure**.

This is not a screening, diagnostic, prognostic or treatment-response model. The dataset contains paired tumour and matched non-tumour tissue from the same 42 patients, not an independent clinical-use population. A model can separate these two tissue preparations while having unknown performance in healthy people, benign disease, premalignancy, other cancers, new laboratories or prospective samples.

## Validation design

Twenty-five independently seeded outer repeats each divide the 42 patients into six folds of seven patients. Every outer test fold therefore contains seven tumour and seven matched non-tumour specimens. The remaining 35 patients enter five grouped inner folds. Both specimens from a patient always share a fold.

Every inner and outer training fit independently learns feature eligibility, abundance medians, scaling parameters, supervised feature ranking and model parameters. Abundance missing values use training-sample medians only. Detection is modelled as a separate binary view without abundance imputation. The combined view concatenates fold-local abundance and detection candidates. No Phase 3 significant list, Phase 4 pathway ranking, or full-data PCA component enters a model filter. This structure follows the need for nested validation when feature selection and tuning occur in high-dimensional, small-sample studies ([Varma and Simon, 2006](https://doi.org/10.1186/1471-2105-7-91); [Lewis et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37113250/)).

Elastic-net logistic regression is the primary learner, using the frozen 28-point grid formed by seven `C` values and four L1 ratios. Elastic net is appropriate for correlated, high-dimensional predictors but does not make the resulting coefficients causal ([Zou and Hastie, 2005](https://doi.org/10.1111/j.1467-9868.2005.00503.x)). Linear SVM is the comparator, with seven `C` values and fold-local candidate counts of 25, 50 or 100. SVM scores receive Platt-style calibration fitted only to inner out-of-fold scores. Decision thresholds are also chosen from inner out-of-fold predictions.

The experiment evaluates nine pipelines: elastic net and linear SVM for abundance, detection and combined views; a specimen-coverage logistic baseline; a five-component abundance PCA logistic baseline whose imputation, scaling, variance filtering and PCA are refitted in every fold; and an intercept-only null baseline. In total, 18,900 outer-test predictions are saved.

## Outer-test performance

The table reports medians and 2.5th–97.5th percentiles across 25 repeated outer partitions. These intervals describe sensitivity to fold assignment. They are not independent-cohort confidence intervals because the same 42 patients occur across repeats.

| View and model | ROC AUC | Balanced accuracy | Sensitivity | Specificity | Brier score |
| --- | ---: | ---: | ---: | ---: | ---: |
| Abundance elastic net | 0.982 (0.977–0.989) | 0.964 (0.964–0.969) | 0.976 (0.976–0.986) | 0.952 (0.952–0.952) | 0.076 (0.059–0.079) |
| Abundance linear SVM | 0.980 (0.972–0.994) | 0.964 (0.940–0.969) | 0.976 (0.952–1.000) | 0.929 (0.929–0.952) | 0.075 (0.064–0.106) |
| Combined elastic net | 0.974 (0.949–0.982) | 0.952 (0.929–0.969) | 0.976 (0.929–1.000) | 0.929 (0.929–0.952) | 0.070 (0.056–0.079) |
| Combined linear SVM | 0.974 (0.962–0.991) | 0.952 (0.919–0.969) | 0.952 (0.933–1.000) | 0.929 (0.905–0.952) | 0.079 (0.066–0.106) |
| Detection elastic net | 0.954 (0.929–0.965) | 0.940 (0.924–0.952) | 0.976 (0.943–0.976) | 0.929 (0.881–0.929) | 0.092 (0.074–0.101) |
| Detection linear SVM | 0.960 (0.932–0.975) | 0.940 (0.917–0.952) | 0.976 (0.929–0.976) | 0.929 (0.905–0.929) | 0.125 (0.096–0.141) |
| Fold-local abundance PCA baseline | 0.969 (0.941–0.992) | 0.940 (0.924–0.964) | 0.952 (0.929–0.976) | 0.929 (0.895–0.952) | 0.046 (0.040–0.082) |
| Coverage-only baseline | 0.868 (0.857–0.874) | 0.786 (0.769–0.810) | 0.833 (0.762–0.843) | 0.762 (0.714–0.786) | 0.144 (0.142–0.148) |

Abundance elastic net has the highest median AUC and balanced accuracy, and it clearly exceeds the coverage-only baseline. The fold-local PCA baseline is nevertheless close in discrimination and has the best Brier score. The high-dimensional abundance model therefore adds modest ranking performance over a low-dimensional unsupervised abundance representation, not a wholly new predictive capability.

Calibration estimates are unstable and show that several complex models produce probabilities that are too conservative in this internal sample. For abundance elastic net, the median calibration intercept is −0.406 and the median slope is 3.53. The PCA baseline has an intercept near zero and slope near one, but all calibration estimates use only 84 repeated outer-test predictions per repeat. These probabilities must not be transported to a clinical setting or interpreted as disease risk.

## Paired-label null

The permutation target was frozen before execution as combined elastic net on outer repeat 0. For each of 200 permutations, the N/T labels were randomly swapped within each patient and the complete six-fold nested pipeline was rerun, including filtering, imputation, scaling, supervised selection and all 28 hyperparameter combinations. The observed AUC is 0.9768. Null AUCs range from 0.3135 to 0.7319, with a median near 0.4785. None reaches the observed value, giving the finite-sample empirical p-value `(1 + 0)/(1 + 200) = 0.00498`.

This result rejects exchangeability of the paired tissue labels for the prespecified pipeline within this dataset. It does not test external transportability, cancer specificity or clinical benefit. The abundance model's slightly higher median AUC was not substituted into this test after results were observed.

## Stability catalogue

Feature stability is calculated from 150 outer training fits per model/view. A stable exploratory row requires selection in at least 60% of outer fits and coefficient-sign agreement of at least 80%. Held-out permutation importance is calculated only on the corresponding outer-test fold.

The tuned elastic nets are not compact: their median selected size is 97 and the modal hyperparameter uses L1 ratio 0.1, which behaves largely as a grouped/ridge-like model. Consequently, Phase 5 does **not** declare a compact panel. It records 222 stable elastic-net view-feature rows representing 148 unique modality-qualified features; only 70 have positive mean held-out permutation importance. Correlated predictors can substitute for one another, and permuting one redundant feature may produce a negligible or negative importance estimate.

Highly recurrent rows include abundance SERPINH1, PRELP and HSPH1; combined abundance PRELP and SERPINH1; and detection NEDD1. Other recurrent rows include GNA11, MAOB, COLGALT1, OGN, SPTAN1, SLC3A2, NOMO1 and SLC38A2. Their coefficient signs indicate tissue-state prediction direction after fold-local standardization. They do not measure causal effects and do not constitute a validated assay panel. Stability selection motivates repeated selection over a single fitted coefficient set, but its formal error guarantees do not automatically apply to this adapted nested workflow ([Meinshausen and Bühlmann, 2010](https://doi.org/10.1111/j.1467-9868.2010.00740.x)).

## Learning behaviour and model-card conclusion

The descriptive learning curve uses the modal nested hyperparameters and repeatedly draws grouped training sets of 14, 21, 28 and 35 patients. Median AUC is already approximately 0.98 at 14 patients and reaches 1.00 in the small seven-patient test complements at 35 training patients. The ceiling behaviour and shrinking test sets make this unsuitable for extrapolating future sample size. It primarily indicates that the tumour/non-tumour contrast is strong within this paired cohort.

The internally validated result satisfies the computational M5 contract: outer-test-only metrics, saved grouped folds, full inner tuning surfaces, paired permutation evidence, held-out importance and explicit leakage-failure tests are present. Reporting follows the separation between model development and clinical validation emphasized by TRIPOD and PROBAST ([Collins et al., 2015](https://doi.org/10.7326/M14-0697); [Wolff et al., 2019](https://doi.org/10.7326/M18-1376)). PROBAST would still identify high concern about applicability and likely risk from the small single-cohort design, absent external validation and incomplete clinical metadata.

The permitted conclusion is narrow: the supplied proteomic abundance and detection data contain a reproducible internal signal that distinguishes paired tumour from matched non-tumour tissue under patient-grouped nested validation. The phase cannot support screening, diagnosis, prognosis, patient-level risk, biological causality or a deployable protein panel.

