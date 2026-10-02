# Phase 6 implementation record

## Architecture and frozen contract

`config/06_robustness.yml` freezes the multiverse factors, inferential thresholds, integration thresholds, pathway stability criteria, permutation count, leading ML pipeline and milestone status. `src/oral_cancer/robustness.py` contains reusable multiverse, patient-deletion and pathway-stability functions. `scripts/60_robustness.py` owns inference and pathway stress tests, `scripts/61_ml_influence.py` owns grouped predictive influence and negative controls, and `scripts/70_build_evidence_table.py` owns evidence integration, claims and figures. `scripts/06_run.py` provides orchestration, tests, hashing and rerun comparison.

The separation is intentional. Adding a representation, cohort, estimator or downstream integration rule requires a configuration or module-level change rather than edits scattered through one notebook. Heavy detail tables are compressed; review tables, JSON validations and figures remain directly inspectable. The source workbook is never edited.

## Robustness implementation

The inference multiverse is the Cartesian product of two representations, two cohorts, three complete-pair rules, one missingness strategy and two estimators. Every run recalculates feature eligibility and the frozen Phase 3 core rule. Feature summaries retain eligibility frequency, sign agreement, FDR and core support frequencies, effect extrema and rank variability. The singleton missingness dimension is explicit so the record does not falsely suggest that Phase 2 alternatives were tested.

Exact patient-deletion refits rerun both abundance and detection calculations after removing one complete pair. Pathways are reranked for all 24 scenarios against frozen Reactome release 86, with direction agreement and leading-edge Jaccard overlap retained separately. This implements sensitivity analysis in the spirit of multiverse analysis rather than selecting whichever scenario produces the preferred list ([Steegen et al., 2016](https://doi.org/10.1177/1745691616658637)).

The predictive influence routine performs 42 outer leave-one-patient-out fits. For each outer deletion it reconstructs the abundance pipeline on the remaining 41 patients and uses five grouped inner folds for tuning; fold-local eligibility, median filling, scaling, supervised selection and threshold selection never see the deleted pair. `make_group_folds` was generalized to allow unavoidable one-patient differences between fold sizes while preserving patient grouping.

The first negative control reuses the Phase 5 full nested paired-label experiment. The second fixes the observed Phase 3/Phase 5 set sizes and randomizes feature keys 1,000 times. This tests excess overlap, not biological independence. Both controls use the frozen seed.

## Integration and claim controls

The evidence table joins only on `feature_key` and validates one-to-one cardinality. Inferential, pathway and predictive columns remain separate. Missing ML or pathway evidence remains `NaN` rather than being treated as contradictory evidence. The composite tiers are deterministic boolean rules, not a learned score, and each row receives inherited limitation IDs and a scope-limited allowed claim.

The claim ledger contains 14 substantive statements spanning source structure, QC, Phase 3 inference, Phase 4 pathways and heterogeneity, Phase 5 model scope, Phase 6 robustness, both negative controls, final prioritization and prohibited conclusions. Each row names a saved evidence file, field locator, evidence state, limitation IDs and allowed scope. This supports FAIR-style traceability and the minimum-information intent of MIAPE while making absent proteomics metadata explicit ([Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18); [Taylor et al., 2007](https://doi.org/10.1038/nbt1329)).

## Verification state

Thirty Phase 0–6 unit tests pass. Phase 6 validations confirm the 24-scenario Cartesian grid, all 42 patient deletions, frozen pathway release, no invented missingness strategy, unique 8,071-row integration, preservation of missing evidence, limitation-aware claims and the open Phase 2 dependency. A same-environment rerun reproduced 22 compared outputs byte-for-byte. The environment record captures the operating system, Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, SciPy 1.18.1, scikit-learn 1.7.2 and the input Git revision.

M6 is not marked closed. The roadmap requires both a completed Phase 2 preprocessing benchmark and full regeneration from a clean independent environment. The current work satisfies the Phase 6 computational implementation and produces the final data-only evidence package, but retaining those two gates prevents a reproducibility claim broader than the evidence.

