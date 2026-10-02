# Phase 1 work log

## 1 October 2026

Phase 1 was opened after completing the initial Phase 0 reconstruction. The agreed objective is to preserve the source workbook and produce reproducible manifests and audits rather than manually cleaning the Excel file.

The initial Phase 1 implementation creates a sample manifest, feature identifier audit, data dictionary, and validation report. Clinical and technical fields that cannot be recovered are carried forward as explicit `unknown` values. This is deliberate: an unknown value is not equivalent to a negative finding or an absent condition.

The source row is retained as the feature unit because the workbook lacks stable protein accessions and peptide evidence. Duplicate symbols and grouped genes are flagged without averaging or splitting. This preserves provenance and prevents unsupported biological precision.

The generated tables have been accepted as the structural hand-off. Current external annotations for unambiguous single-gene rows are optional and may be added in parallel; they do not delay Phase 2 QC.

## Initial implementation result

The import and validation script completed successfully and reproduced identical output hashes on a second run. It confirmed 84 unique specimen columns, 42 patient IDs, and exactly one tumor and one matched non-tumor specimen for every patient.

The feature audit contains all 8,071 source rows. It records 114 all-missing rows, eight rows belonging to four duplicated identifiers, and 48 semicolon-delimited group rows. All 48 grouped rows are also entirely unmeasured in the supplied abundance matrix. They will remain in the audit record but cannot contribute to numerical modelling.

No sample has been excluded. No identifier has been averaged, split, or forcibly mapped. Phase 1 now has a reproducible structural baseline; the remaining work is optional external annotation of unambiguous identifiers and review of the manifest fields before beginning Phase 2 QC.

## European alignment update — 1 October 2026

The Phase 1 structure was checked against current WHO/IARC, EHNS–ESMO–ESTRO, UICC, GDPR, EHDS, and European proteomics-repository practice. No source or generated audit file was overwritten, and the structural counts and hashes remain unchanged.

The next versioned manifest will add exact site, pathology entity and classification edition, cTNM/pTNM components and edition, comparison-tissue pathology, clinically relevant oral-cavity features, exposure and geography fields, and source/verification metadata. It will distinguish `unknown`, `not_applicable`, and `withheld` so absence of knowledge is not confused with inapplicability or privacy controls.

A governance-limitation statement and EHDS-oriented catalogue description are optional documentation extensions. HGNC/UniProt annotation remains additive and release-versioned. Because sharing authority cannot be established from the current files, participant-level repository deposition is outside the executable data-only roadmap; PRIDE/ProteomeXchange remains a future route only if independent authority is later established.

## Documentation update

The complete user-visible response is preserved in `records/01_DETAILED_RESPONSE_RECORD.md`, while the maintained technical protocol is `01_PLAN_AND_DECISIONS.md`. Methodological support for identifier preservation, reproducibility, missingness modelling, paired inference, pathway analysis, machine learning, and validation is maintained in `RESEARCH_EVIDENCE_BASE.md` and cited in the governing roadmap.

European terminology and implementation decisions are maintained in `EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md` and are linked from both detailed Phase 1 records.

## Data-only transition — 1 October 2026

The structural Phase 1 milestone is accepted as complete without waiting for clinical, pathology, batch, acquisition or governance information. Those fields remain explicit unknowns governed by `../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`. Optional HGNC/UniProt annotation may improve current nomenclature but cannot reconstruct original protein groups or block Phase 2.

Phase 1 now hands the immutable workbook hash, 84-row sample manifest, 8,071-row feature audit, data dictionary, paired-detection counts and validation report to Phase 2. Phase 2 must write its QC metrics and flag rules before ranking samples, and must return a locked primary cohort and sensitivity cohorts to Phase 3.

