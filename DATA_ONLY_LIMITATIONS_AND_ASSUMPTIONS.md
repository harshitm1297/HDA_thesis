# Data-only limitations and assumptions register

## Register policy

This register tracks limitations caused by proceeding without new information from the supervisor, laboratory, clinical team, or project group. A limitation is never removed because a common practice seems likely. Under the current scope, it is closed only by evidence present in the supplied workbook. Research available through 2023 may justify a method or qualify an interpretation, but it cannot fill a project-specific fact. Analytical mitigations reduce bias or test robustness; they do not convert an unknown project fact into a known one.

Phase numbers refer to the from-scratch Phase 0–6 sequence in `FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md`.

Evidence states used throughout the repository are `verified_source`, `derived_from_matrix`, `literature_supported_assumption`, `unresolved`, and `not_testable_with_current_data`.

## Permanent and cross-phase limitations

| ID | Missing information | Consequence | Data-only mitigation | Claims that remain prohibited | Affected phases |
| --- | --- | --- | --- | --- | --- |
| L01 | Raw mass-spectrometry files and complete quantitative export | Identification, interference and instrument-level QC cannot be rechecked | Preserve reduced matrix; use sample-level statistical QC; prioritise stable features; document unverifiable upstream evidence | Reprocessed raw-data validation; instrument-performance certification | 0–3, 6 |
| L02 | Instrument, acquisition settings, Spectronaut version, quantity field and upstream normalisation | Exact measurement generation and scale are uncertain | Use the cautious measurement statement; compare a small set of plausible preprocessing scenarios | Exact DIA method, exact quantity definition, or claim that values are raw intensities | 0–3 |
| L03 | Search database, peptide evidence, accessions, protein inference and identification FDR | Gene labels cannot prove individual proteins or isoforms | Preserve original rows and identifiers; flag ambiguity without silently remapping or splitting them | Isoform-specific claims; reconstructed identification FDR; automatic splitting or averaging | 0, 3–6 |
| L04 | Batch, run order, preparation batch and technical controls | Technical structure cannot be assigned or corrected confidently | Inspect latent structure; compare robust sensitivity scenarios; do not use spreadsheet order as batch | Batch-adjusted causal interpretation; proof that an outlier is a laboratory failure | 0–3 |
| L05 | Exact anatomical site, histology, stage, grade, HPV, purity, inflammation, exposures and treatment | Clinical heterogeneity and confounding cannot be explained or adjusted | Use broad verified labels; restrict analysis to paired tissue status; avoid clinical subgroup claims | Uniform OCSCC label; stage/HPV/exposure effects; adjusted clinical models | 0, 1, 4–6 |
| L06 | Distance and pathology of the non-tumour specimen | Comparison tissue may contain field change, dysplasia or inflammation | Use “matched non-tumour”; interpret contrast as within-patient tissue difference | Healthy-control, cancer-free population, or early-detection claims | 0, 3–6 |
| L07 | Cause of each missing value | Missingness cannot be assigned to low abundance versus technical filtering | Preserve missingness; analyse detection separately; use imputation only as sensitivity | Claim that blank equals zero or biological absence | 1–3, 5 |
| L08 | Independent cohort and orthogonal assays | Generalisability and assay transfer are unproven | Use nested grouped internal validation, permutation tests, sensitivity analysis and influence diagnostics; preserve external validation as future work | Externally validated biomarker, clinical panel or diagnostic test | 3–6 |
| L09 | Clinical outcomes and appropriate population controls | Prognosis, response and screening cannot be studied | Do not create surrogate outcomes; restrict objective to paired discovery and tissue classification | Survival, recurrence, treatment-response, screening or population-risk prediction | 0, 4–6 |
| L10 | Ethics approval, consent scope, controller and sharing conditions in supplied records | Public or cross-border release cannot be authorised from the analytical files | Keep data local and pseudonymous; publish only non-disclosive artifacts unless authority is established | Claim of GDPR/EHDS-compliant release or unrestricted public deposition | 0, 6 |
| L11 | Only 42 independent patients | High-dimensional estimates and prediction are unstable | Paired modelling, moderated variance, strict multiplicity control, sparse models, nested grouped resampling and uncertainty reporting | Clinical-grade model performance or dependable fine-grained subtypes | 2–6 |
| L12 | No healthy, benign inflammatory, premalignant or asymptomatic groups | Specificity for cancer and early disease cannot be measured | Treat classifier as tumour-versus-matched-non-tumour tissue experiment | Screening, early diagnosis or differential diagnosis | 0, 5–6 |

## Operational assumptions

| ID | Working assumption | Basis | Confidence | Required sensitivity or safeguard | Revisit in |
| --- | --- | --- | --- | --- | --- |
| A01 | `P<number>N` and `P<number>T` represent a matched pair from one patient | Consistent naming and supplied report | High | Enforce exactly one N and one T; group all resampling by patient | Every phase |
| A02 | Positive values are not already on a conventional log2 scale | Range from approximately 0.0023 to 7.31 billion and prior report | High | Compare log2 representations; retain original values unchanged | 1–2 |
| A03 | Blank abundance cells are missing or unreported measurements | Workbook structure | High for representation; low for mechanism | Never convert to zero; analyse detection; compare missingness scenarios | 1–3, 5 |
| A04 | Each source row is the safest available quantitative feature unit | `PG.Genes` protein-group terminology and absent accessions/peptides | High | Keep stable row keys; do not split or average; annotate additively | 0, 3–6 |
| A05 | Patient is the independent biological unit | Paired design | Verified | Include patient block in inference and patient grouping in resampling | 1–5 |
| A06 | No additional normalisation and median centring are the two primary plausible candidates | Unknown upstream normalisation and restrained workflow design | Moderate | Compare using prespecified matrix diagnostics before differential testing | 2 |
| A07 | No-imputation abundance analysis is the primary inference | Missingness mechanism unknown and proteomics benchmarks | Moderate-to-high | Separate detection analysis; compare limited imputation/hurdle sensitivities | 2–3 |
| A08 | Literature through 2023 may justify methods and contextualise results but cannot replace project-specific facts | Current data-only scope and provenance principle | High | Keep workbook-derived results separate from literature-supported interpretation; do not import external patient data, labels, outcomes or measurements | All phases |

## Limitation propagation rule

Every phase report must list the limitation IDs that affect it. A later phase cannot silently drop a limitation inherited from an earlier phase. If a candidate protein depends on a duplicated symbol, sparse detection, one preprocessing scenario, or an unstable patient subset, those dependencies remain attached to the candidate in pathway, prediction and final reports.

The candidate evidence table will therefore include at least `feature_key`, original identifier, annotation status, usable pairs, detection pattern, primary effect, adjusted evidence, preprocessing stability, sample-set stability, limitation IDs, and allowed claim. This turns the register into an analytical control rather than a disclaimer appended after results are known.

## Review cadence

The register is reviewed at every milestone. New limitations may be added when analysis exposes a dependency. Existing limitations can be narrowed only when the current supplied files provide evidence. Under the data-only operating position, project-group contact and external datasets are not expected resolution routes.

