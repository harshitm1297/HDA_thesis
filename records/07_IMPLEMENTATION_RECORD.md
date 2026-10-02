# Phase 7 implementation record

## Definition and architecture

Phase 7 was added after the original Phase 0–6 roadmap as a validation-readiness extension. `config/07_validation_readiness.yml` freezes candidate-source tier, shortlist capacity, correlation rule, planning alpha, power and assumption grids. `src/oral_cancer/validation.py` implements deterministic nondominated sorting, paired-change correlation components and sample-size calculations. `scripts/80_validation_readiness.py` creates the evidence and protocol package; `scripts/07_run.py` reruns, tests, fingerprints and compares outputs.

The source scope remains unchanged. Phase 7 reads the Phase 6 integrated table and Phase 0 paired-change matrix. It does not query public expression cohorts, add literature-derived patient labels, or fit another classifier to the discovery outcomes.

## Selection safeguards

The implementation has no weighted composite score. Pareto sorting is branch-specific and maximizes explicitly named evidence dimensions. Missing values are treated as unavailable during dominance checks, not silently converted into evidence. Shortlist construction starts only from Phase 6 `high_priority_internal` rows and applies frozen eligibility, Pareto and tie-break rules.

Abundance redundancy uses pairwise-complete Spearman correlations of actual patient-level log2 tumour-minus-non-tumour changes. It does not correlate final summary statistics or gene annotations. Connected components are deterministic under sorted feature keys, and the full component membership and edge list are saved. Detection candidates are not forced into this abundance-correlation structure.

The power tables deliberately use generic effect assumptions rather than the largest selected discovery estimates, reducing winner's-curse optimism. Paired continuous calculations use the noncentral t distribution. Detection calculations are labelled normal approximations and must be replaced after assay-specific discordance assumptions exist.

## Verification

Five Phase 7 tests cover dominance behaviour, sample-size monotonicity, the 165-to-20 candidate contract, the unexecuted-validation protocol and validation checks. They run alongside all earlier tests. The default runner generates Phase 7 twice and requires exact hashes for every non-document core output.

The computational readiness package can pass; independent validation cannot pass without new patients and measurements. This distinction is encoded in every shortlisted row as `validation_state=not_executed`, in the JSON protocol status and in the Phase 7 milestone string.
