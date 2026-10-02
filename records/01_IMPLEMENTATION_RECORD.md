# Phase 1 implementation record

## Implemented workflow

`scripts/10_qc_missingness.py` consumes the immutable Phase 0 objects rather than reparsing sample identities. It calculates per-specimen coverage and positive-log2 distribution summaries, matched-pair detection and abundance concordance, pairwise-complete Spearman correlation, label-blind robust-scaled PCA, deterministic average-linkage ordering, an MCD-style multivariate distance, five independent diagnostic-family flags, pair-safe cohort membership, feature-level paired detection states, exact descriptive paired-detection p-values, and abundance-decile missingness summaries.

`scripts/01_run.py` reruns the analysis, executes the full test suite, and fingerprints the generated outputs. Configuration is frozen in `config/01_qc_missingness.yml`. The workflow needs only the pinned Python, NumPy, pandas, and openpyxl environment already recorded for Phase 0.

The original workbook was not edited. Full-matrix imputation was not performed. The only imputation occurs in a temporary high-coverage matrix used for diagnostic PCA, is fitted without tissue labels, and is not exported for inference. Hierarchical ordering is learned without tissue labels; annotations are applied only after ordering.

## Recorded outcome

M1 passed with 3,759 high-coverage QC features. Twelve specimens from seven patient pairs met the prespecified two-family sensitivity rule. The primary analysis retains all 42 pairs; the sensitivity cohort contains 35 complete pairs. The detailed numerical interpretation and boundaries are recorded in `results/01_qc_missingness/01_RESULTS.md`.

## Reproduction command

```powershell
python scripts/01_run.py
```

