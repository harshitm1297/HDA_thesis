# Phase 4 implementation record

## Frozen inputs and resource provenance

Phase 4 consumes the immutable Phase 0 data model and Phase 3 no-imputation outputs. It never modifies the workbook and never fills a missing protein value. The external annotation resource is Reactome release 86 from 7 September 2023. The downloaded ZIP and extracted GMT are fingerprinted in `resources/reactome_v86/README.md`; the analysis refuses to run if the GMT checksum differs from the value frozen in `config/phase4.yml`.

The Reactome file supplies only gene-set membership. This is an explicit exception to the data-only measurement rule: prior biological knowledge may organize the supplied observations, but it cannot add a participant, abundance, covariate, label or outcome. The release is archived locally so later Reactome changes cannot silently alter results.

## Modular implementation

`src/oral_cancer/pathways.py` owns GMT parsing, deterministic gene-symbol collapse, weighted ranked enrichment, leading-edge extraction, robust patient pathway scoring and paired sign-flip inference. `src/oral_cancer/heterogeneity.py` owns complete-case robust PCA, bootstrap explained-variance summaries, deterministic multi-start k-means, silhouette calculation and consensus clustering. `src/oral_cancer/visuals.py` owns Phase 4 scientific figures. `scripts/40_pathways_heterogeneity.py` is a thin orchestrator, and `scripts/run_phase4.py` regenerates outputs, runs the full test suite and fingerprints deliverables.

Every numerical choice is in `config/phase4.yml`: pathway-size limits, 1,000 ranked permutations, 20,000 sign permutations, q=0.05, complete-case requirements, 500 heterogeneity features, robust clipping, 200 PCA bootstraps, four candidate k values, 200 consensus resamples and all cluster-acceptance thresholds. The seed is fixed at 20231031 and offset deterministically for separate procedures.

## Safeguards added during execution

The patient-level pathway statistic preserves zero as the biological null. Patient deltas are divided by a robust feature scale but are not recentered on their observed median. A synthetic test with uniformly positive paired changes verifies that the output remains positive and significant. This prevents a common analytical error in which preprocessing removes the cohort shift that the subsequent sign-flip test is intended to detect.

Compound identifiers cannot enter the ranked symbol table. Duplicate symbols use a documented largest-absolute-moderated-statistic rule and remain traceable through `gene_collapse_audit.csv`. Patient pathway scoring is stricter: only unambiguous features complete in all 42 pairs enter. PCA likewise uses only all-patient-complete changes and performs no imputation. A machine-readable interpretation label prevents exploratory cluster assignments from being called clinical subtypes.

The test suite covers Phase 0, Phase 1, Phase 3 and Phase 4 contracts. Phase 4-specific tests check enrichment direction, zero-preserving patient scores, Reactome availability and release, pathway-size filters, unique integrated rows, all 42 patient outputs and the non-clinical cluster label. The successful run contains 18 passing tests.

## Result state and hand-off

The phase tests 783 ranked Reactome pathways and 749 complete-case paired-score pathways. There are 271 same-direction pathways passing q≤0.05 in both internal views. Patient-change PCA uses 500 of 1,350 all-patient-complete features. No k=2–5 cluster solution meets all frozen stability criteria, so no stable patient grouping is declared.

The implementation status is `PASS`; the scientific milestone is `CONDITIONAL_ON_PHASE2_AND_FULL_M3_CLOSURE`. Phase 5 can consume the frozen tables as candidate representations, but any supervised selection, scaling, pathway choice or tuning must be learned anew within each training fold. The full-data Phase 4 pathway ranking is descriptive and cannot be inserted unchanged into cross-validation without leakage.

