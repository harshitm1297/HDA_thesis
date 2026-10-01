# Oral-cancer paired proteomics project

This repository analyses the supplied 84-specimen proteomics workbook as 42 matched tumour/non-tumour patient pairs. The empirical analysis is data-only: literature available through 2023 may justify methods, but it cannot replace missing measurements or metadata.

Phase 0 is the implemented foundation. It fingerprints and validates the immutable workbook, constructs traceable sample and feature manifests, materialises raw/detection/paired data objects, freezes the scientific contract, and runs automated tests.

## Reproduce Phase 0

From this directory, using Python 3 with the pinned packages in `environment/requirements-phase0.txt`:

```powershell
python scripts/run_phase0.py
python -m unittest discover -s tests -v
```

The source workbook is intentionally excluded from Git. Generated patient-level matrices under `data/interim/` are also excluded; they can be regenerated from the immutable workbook. Machine-readable validation and the human-readable results report are under `results/phase0/`.

## Analysis order

1. Phase 0: research contract and immutable data model — implemented.
2. Phase 1: multivariate QC and missingness intelligence.
3. Phase 2: preprocessing and missing-data benchmark.
4. Phase 3: paired abundance and detection inference.
5. Phase 4: pathway interpretation and heterogeneity.
6. Phase 5: leakage-safe AI and stable-panel discovery.
7. Phase 6: stress testing and reproducible synthesis.

