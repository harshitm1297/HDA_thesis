# Phase 5 implementation record

## Architecture

`src/oral_cancer/ml.py` contains leakage guards, patient-fold construction, fold-local abundance/detection processing, elastic-net and SVM fitting, inner threshold selection, SVM calibration, PCA processing and metrics. `src/oral_cancer/phase5.py` owns grouped inner/outer execution, complete prediction capture, coefficient and held-out-importance extraction, baselines and paired label swaps. The three numbered scripts separate the expensive nested run, the permutation null and deterministic finalization.

All statistical choices are stored in `config/phase5.yml`. Scikit-learn 1.7.2 is pinned in `environment/requirements.txt`. The input workbook is never edited; Phase 5 consumes the validated Phase 0 matrices.

## Leakage controls

Outer folds are written before model fitting. Split validation rejects sample overlap, patient overlap and split pairs. Transformers require explicit fit indices contained inside the permitted training scope. Supervised selectors reject any feature source other than fold-local training labels. Hyperparameter-selection code has no outer-test metric input, and a guard deliberately rejects such an input. PCA, median imputation, standardization, supervised ranking, SVM calibration and thresholds are refitted from training data.

The tests deliberately exercise failure states for overlapping patients, split pairs, preprocessing that includes a test sample, use of a full-data supervised feature list and outer-test-driven tuning. Output tests require all 84 specimens once per repeat/pipeline, 200 permutations, finite probabilities and the non-clinical interpretation label.

## Reproducibility state

The completed run contains 25 × 6 outer fits for each of six high-dimensional model/view combinations, plus three baselines. It stores 18,900 outer-test predictions, complete inner score surfaces, 200 full paired-label nested reruns, coefficient histories, held-out permutation importance, learning-curve results and SHA-256 output hashes.

M5 computational status is PASS. The milestone remains conditional on Phase 2 and upstream M3/M4 closure. Phase 6 must stress-test these conclusions and must not reinterpret the stability catalogue as a finalized biomarker panel.

