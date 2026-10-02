# Phase 5 model card

## Intended use

The models quantify internal discrimination between tumour and matched non-tumour tissue specimens in the supplied 42-patient proteomics dataset. They are research pipelines for methodological evaluation and candidate stability analysis.

## Prohibited use

Do not use the outputs for screening, diagnosis, prognosis, recurrence, treatment selection, individual risk, clinical decision-making or claims about a general oral-cancer population. No external cohort, healthy-control population, benign-disease controls, prospective sampling or assay transfer study is available.

## Inputs and outputs

Inputs are workbook-derived log2 abundance values and binary detection indicators with patient IDs used only for grouped resampling. Outputs are tumour-class probabilities or SVM-calibrated probabilities, fold-specific decisions, performance distributions, coefficient stability and held-out importance.

## Evaluation

Twenty-five repeated six-fold outer patient partitions and five grouped inner folds are used. The strongest median discrimination is abundance elastic net with AUC 0.982 and balanced accuracy 0.964. A fold-local PCA baseline reaches AUC 0.969 and has better probability error. The predeclared combined elastic-net paired-swap test gives p=0.00498 with 200 complete nested reruns.

## Important limitations

Only 42 independent patients are present. Outer-repeat percentile ranges reuse the same patients and are not external confidence intervals. Sample acquisition, batch, tumour purity, histology, site, HPV, stage and outcomes are unavailable. The tuned elastic nets are dense, so no compact panel is claimed. Probability calibration and transportability remain unestablished.

## Reproduction

Run `python scripts/05_run.py` from the project-local environment. `python scripts/05_run.py --reuse-existing` revalidates completed expensive fits, reruns all tests and refreshes hashes without recomputing the nested and permutation loops.

