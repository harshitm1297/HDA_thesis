# Oral-cancer paired proteomics project

This repository analyses the supplied 84-specimen proteomics workbook as 42 matched tumour/non-tumour patient pairs. The empirical analysis is data-only: literature available through 2023 may justify methods, but it cannot replace missing measurements or metadata.

Phases 0 and 1 are implemented and passed. Phase 0 fingerprints and validates the immutable workbook, constructs traceable sample and feature manifests, materialises raw/detection/paired data objects, and freezes the scientific contract. Phase 1 adds specimen QC, label-blind multivariate diagnostics, a feature-level missingness atlas, and frozen primary and sensitivity cohorts. Phase 3's modular no-imputation inference branch and Phase 4's pathway/heterogeneity branch are also implemented; their milestones remain conditional on completion of the Phase 2 preprocessing benchmark.

## Reproduce Phase 0

From this directory, using Python 3 with the pinned packages in `environment/requirements.txt`:

```powershell
python scripts/run_phase0.py
python -m unittest discover -s tests -v
```

The source workbook is intentionally excluded from Git. Generated patient-level matrices under `data/interim/` are also excluded; they can be regenerated from the immutable workbook. Machine-readable validation and the human-readable results report are under `results/phase0/`.

## Analysis order

1. Phase 0: research contract and immutable data model — implemented, M0 passed.
2. Phase 1: multivariate QC and missingness intelligence — implemented, M1 passed.
3. Phase 2: preprocessing and missing-data benchmark.
4. Phase 3: paired abundance and detection inference — modular no-imputation branch implemented; full M3 remains conditional on Phase 2.
5. Phase 4: pathway interpretation and heterogeneity — implemented, computational checks passed; M4 remains conditional on Phase 2 and full M3 closure.
6. Phase 5: leakage-safe AI and stability-ranked feature discovery — implemented, computational checks passed; no compact panel declared and M5 remains conditional on upstream closure.
7. Phase 6: stress testing and reproducible synthesis.

The pre-2024 evidence and locked evaluation rules for Phase 2 imputation and PCA are recorded in `PHASE2_MISSINGNESS_AND_PCA_EVIDENCE_2023.md`.

## Reproduce Phase 4

```powershell
python scripts/run_phase4.py
```

The command verifies the frozen Reactome v86 checksum, regenerates pathway and patient-heterogeneity outputs, runs all Phase 0–4 tests, and writes an output manifest. The detailed interpretation is in `results/phase4/PHASE4_RESULTS.md`; the software and safeguards are recorded in `records/PHASE4_IMPLEMENTATION_RECORD.md`.

## Reproduce Phase 5

```powershell
python scripts/run_phase5.py
```

The full command runs 25 repeated grouped nested validations and 200 complete paired-label permutation pipelines, so it is intentionally computationally expensive. Use `python scripts/run_phase5.py --reuse-existing` only to revalidate and fingerprint already completed fits. Results and model limitations are recorded in `results/phase5/PHASE5_RESULTS.md` and `results/phase5/MODEL_CARD.md`.

