# Phase 3 implementation record

## Architecture

The Phase 3 implementation introduces the reusable `oral_cancer` package under `src/`. `data.py` owns stable-axis loading, log2 representations and pair construction. `statistics.py` owns Student-t probabilities, Benjamini–Hochberg adjustment, exact paired binary probabilities and empirical-Bayes variance-prior estimation. `phase3.py` owns abundance inference, detection inference and evidence integration. `visuals.py` owns static scientific figures. The scripts directory contains only thin orchestration entry points.

This separation keeps data contracts independent from statistical backends. A future limma, MSqRob2 or missingness-aware likelihood backend can produce the same table schema without rewriting cohort or identifier logic. Every numerical decision is stored in `config/phase3.yml`.

## Reproduction

```powershell
python scripts/run_phase3.py
```

The command regenerates all Phase 3 outputs, runs the complete Phase 0–3 test suite, and fingerprints each output. The test suite includes distribution reference values, non-zero extreme-tail probabilities, BH monotonicity, a synthetic paired-effect sign fixture, pair-support rules, output dimensions and validation checks.

## Milestone record

The Phase 3 computation and primary no-imputation branch pass. Full M3 closure is conditional on completing the broader Phase 2 preprocessing benchmark. This condition is machine-readable in `results/phase3/phase3_validation.json`.

