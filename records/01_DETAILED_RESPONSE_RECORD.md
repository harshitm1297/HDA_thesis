# Phase 1 detailed response record

## Record date

1 October 2026

## User requests

> for phase 0, note down your complete responses in markdown files and store them in the folder repository. and then, focus on phase 1. we will create a record of our conversations on these phases and milestones through markdown files

> also explain what you did in phase 1, is it recorded word to word in your md files? if not, create detailed markdown files for the phases thus far, putting your detailed responses in this chat in them

## Full explanation of Phase 1

Phase 1 created a reproducible structural representation of the oral-cancer dataset. It did not normalize values, impute missing measurements, exclude samples, merge identifiers, or perform differential analysis. Its purpose was to make the source matrix traceable and safe for later QC and modelling.

### Source preservation

The original workbook was treated as immutable. The import script reads `OC_Dataset_84Sample_v2_24012024.xlsx` but does not write to it. Its SHA-256 checksum is recorded as `65211B375E366D0CF6078E1016A3D29867EFC5E0EBA10D97F1F0282B25161DA8`, allowing later checks that the analytical input has not changed.

The temporary Excel lock file beginning with `~$` was excluded from analysis. The script requires the expected worksheet name `Sheet1` and expected identifier header `PG.Genes`; it stops with an error if these structural assumptions no longer hold.

### Reproducible import script

The original single-file audit implementation was superseded by the maintained modular workflow. [01_build_data_model.py](../scripts/01_build_data_model.py) constructs the manifests and immutable data objects, while [10_qc_missingness.py](../scripts/10_qc_missingness.py) performs specimen QC, feature-level missingness analysis, multivariate diagnostics, validation, and output generation.

The process was executed twice and produced identical hashes for every output. This verifies deterministic regeneration from the same source workbook and software environment.

This use of checksums, stable keys, a data dictionary, and scripted regeneration follows the FAIR principles for findable, accessible, interoperable, and reusable scientific data: [Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18).

### Sample-name parsing and pair validation

Every sample column was required to match the pattern `P<number>N` or `P<number>T`. The numeric component defines the patient, `N` means matched non-tumor tissue, and `T` means tumor tissue.

The script verified that all 84 sample names are unique, exactly 42 patient identifiers exist, and every patient has exactly one `N` and one `T` specimen. No malformed, unmatched, or duplicate sample column was found.

The biological sample size is therefore 42 patients, not 84 independent patients. The 84 columns represent correlated specimens arranged into 42 complete pairs.

### Sample manifest

The historical local file `legacy_01_outputs/sample_manifest.csv` contains one row per specimen. It is excluded from version control and superseded by the reproducible current workflow. It records:

- Original sample ID
- Patient and pair ID
- Numeric patient value for stable sorting
- Tissue code and expanded tissue label
- Original worksheet
- Original Excel column letter and position
- Initial inclusion flag
- Exclusion reason field
- Clinical and technical metadata placeholders

The metadata placeholders include oral subsite, histology, stage, grade, node status, age, sex, tobacco, smokeless-tobacco or betel-quid exposure, alcohol, HPV status, treatment before collection, tumor purity, inflammation, non-tumor distance/pathology, collection details, preparation batch, instrument batch, run order, technical replicate, and notes.

These fields contain the literal value `unknown` because the information is absent from the supplied files. An unknown value is not interpreted as “no,” “negative,” or “not applicable.” The explicit placeholders make missing knowledge visible and provide a stable location for future recovered metadata.

### Feature-level provenance

The historical local file `legacy_01_outputs/feature_identifier_audit.csv` contains exactly one row for every one of the 8,071 source feature rows. It is excluded from version control and assigns a stable key such as `feature_00001` while preserving the original Excel row and unmodified `PG.Genes` value.

The stable feature key is row based because the source does not contain stable protein accessions. It ensures that duplicate symbols can remain separate and that future annotation never destroys the connection to the original measurement row.

For each feature, the audit records whether the identifier is missing, how many source rows share it, whether it is duplicated, whether it contains a semicolon-delimited group, and how many group members are represented.

### Completeness and paired-detection calculations

For every feature, the script counts numeric and missing sample values and calculates the missing fraction. It separately records tumor detections, matched-non-tumor detections, and four paired categories:

- Both tumor and non-tumor measured
- Tumor measured and non-tumor missing
- Non-tumor measured and tumor missing
- Both tissues missing

These four counts are required to sum to 42 for every feature. This validation passed for all 8,071 rows.

The counts are descriptive. A blank cell is called missing or unreported, not biologically absent. The paired categories will later support differential-detection and missingness analyses.

For observed values, the audit also records the minimum, median, and maximum. These summaries do not replace the source matrix; they provide quick checks and support later eligibility rules.

### Initial feature disposition

Phase 1 gives each row an initial status without changing the underlying data. Rows with no measurements receive `exclude_all_abundance_missing` because they cannot enter numerical abundance analysis. Duplicate identifiers receive `review_duplicate_identifier`. Semicolon groups and missing identifiers have their own review statuses when applicable. All other rows receive `retain_unresolved`.

This is a structural disposition, not a biological quality score. A retained row has not yet been declared reliable or significant.

The audit found 114 entirely unmeasured rows. It found eight rows representing four duplicated identifiers: `CDKN2A`, `CUX1`, `MOCS2`, and `TMPO`. It also found 48 semicolon-delimited group rows. All 48 grouped rows are among the 114 entirely unmeasured rows, so they remain documented but cannot contribute to quantitative modelling in the supplied matrix.

No feature identifier is blank.

### Identifier policy

No duplicate gene symbols were averaged. Averaging could combine distinct protein groups, accessions, isoforms, or inference outcomes that happen to share a symbol. No semicolon-delimited group was split into independent proteins because the source provides no peptide or accession evidence supporting such a split.

This policy is supported by the [Spectronaut manual](https://biognosys.com/content/uploads/2024/09/Spectronaut-19-manual-v4.pdf), which defines `PG.Genes` at the protein-group level and documents multi-member protein groups. It is also consistent with [HGNC nomenclature guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC7494048/), under which stable HGNC IDs, rather than mutable display symbols alone, provide durable gene identity.

The fields `stable_accession` and `identifier_resolution_note` were added to make the unresolved state explicit. A later annotation step may add current HGNC or UniProt information for unambiguous single-gene rows, but the original row key and identifier will remain authoritative for provenance.

### Data dictionary

The historical local file `legacy_01_outputs/01_data_dictionary.csv` defines each generated manifest and feature-audit field. Each definition is labelled as source, derived, decision, or placeholder provenance. This distinction prevents later users from mistaking an inferred value or project decision for an original laboratory measurement.

### Validation report and hashes

The historical local report `legacy_01_outputs/01_VALIDATION_REPORT.md` records the source hash, dimensions, missingness, identifier issues, validation checks, and hashes of the generated CSV files. It is not part of the maintained GitHub artifact set.

The following checks passed:

- One identifier column and 84 sample columns
- Every sample name matches the required pattern
- Every sample name is unique
- The manifest contains 84 specimens
- The manifest contains 42 patients
- Every patient has exactly one `N` and one `T`
- The feature audit contains all source rows
- Every generated feature key is unique
- Every source Excel row appears once
- Paired detection categories sum to 42 for every feature

The validated dataset summary is:

- 8,071 feature rows
- 84 sample columns
- 42 independent patients
- 42 complete pairs
- 240,380 missing abundance cells out of 677,964, or 35.456%
- 114 all-missing feature rows
- Eight duplicated-identifier rows representing four symbols
- 48 semicolon-delimited group rows, all also all-missing
- Zero missing identifiers
- Sample detected-feature counts ranging from 2,460 to 6,549

### What Phase 1 did not do

Phase 1 did not perform log transformation, normalization, imputation, batch correction, PCA, clustering, sample exclusion, differential abundance, pathway analysis, or machine learning. Those operations depend on later QC and analysis decisions.

It also did not prove stable protein accessions, isoforms, peptide evidence, identification FDR, acquisition mode, or original Spectronaut settings. Those limitations remain governed by the Phase 0 provenance record.

### Phase 1 status and remaining work

The reproducible structural baseline is complete. We now have a validated sample manifest, a row-complete feature audit, a data dictionary, deterministic output hashes, and an executable regeneration script.

Optional Phase 1 follow-up is limited to annotation rather than structural repair. Current HGNC and UniProt annotations may be added to unambiguous single-gene rows, but they do not block Phase 2. Any annotation must be stored in additional fields, include database version and lookup date, and never overwrite `feature_id_original` or imply reconstruction of the original Spectronaut protein inference.

The evidence supporting this and later phase decisions is maintained in [RESEARCH_EVIDENCE_BASE.md](../RESEARCH_EVIDENCE_BASE.md). It includes primary or authoritative sources for proteomics QC, missingness and imputation, paired modelling, pathway analysis, nested validation, prediction-model reporting, and external proteomic repositories.

After that annotation decision, the project can move into Phase 2 technical and sample-level QC.

## Associated records

- [Phase 1 plan and decisions](../01_PLAN_AND_DECISIONS.md)
- [Phase 1 work log](01_WORK_LOG.md)
- Historical Stage 01 validation report: `legacy_01_outputs/01_VALIDATION_REPORT.md` (local and unversioned)
- [Stage 01 QC and missingness implementation](../scripts/10_qc_missingness.py)

## European alignment request — 1 October 2026

> Can you ensure that the research, ideas and terminologies are also closely aligned with what has happened in Europe till now? Also, make the markdown reports so far, more detailed, including what, why and how (without these headlines)

The completed structural audit remains valid and unchanged. European alignment adds a versioned metadata and governance layer rather than rewriting the original CSVs or their hashes.

The current manifest establishes 42 complete patient pairs, but it cannot establish OCSCC, exact oral subsite, HPV-defined disease or a staging system. The next manifest version will add original and mapped site, coding system, pathology wording, WHO/IARC edition, cTNM and pTNM components, stage group, TNM edition and date, depth of invasion, margins, nodal findings, extranodal extension, perineural and lymphovascular invasion, and evidence source. Values remain `unknown` until recovered; they will not be predicted from the proteome. This follows current [WHO/IARC](https://publications.iarc.who.int/Book-And-Report-Series/Who-Classification-Of-Tumours/Head-And-Neck-Tumours-2024), [European multidisciplinary guideline](https://www.annalsofoncology.org/article/S0923-7534%2820%2939949-X/fulltext), and [UICC TNM](https://www.uicc.org/resources/tnm/publications-resources) terminology.

The comparison-tissue record will add exact sampled site, tumour distance, pathology review, dysplasia or another oral potentially malignant disorder, inflammation, margin relationship and contamination or purity. These fields explain why `N` is interpreted as matched non-tumour tissue and not a healthy-population control.

The exposure and transportability layer will retain smoking, alcohol, smokeless tobacco, areca/betel quid, oral submucous fibrosis, HPV method and result, age, sex, geography, centre and socioeconomic variables where justified and lawfully available. It enables later checks of whether a European dataset is genuinely comparable. It does not authorise underpowered subgroup modelling in the 42-patient cohort.

Metadata states will distinguish `unknown`, `not_applicable` and `withheld`, and populated fields will carry a source, verification status and update date. A separate controlled governance inventory will hold ethics approval, consent scope or other basis, controller, permitted secondary use, access level, retention, transfer and deposition restrictions, and re-identification-key arrangements. This separation makes a shareable manifest possible without confusing privacy-driven withholding with absent science. [GDPR](https://eur-lex.europa.eu/eli/reg/2016/679/oj).

An EHDS-oriented dataset catalogue will document purpose, population, geography and dates, molecular data category, specimen design, coding systems, formats, quality and utility limitations, linkage potential, access conditions, holder and version history. It prepares the dataset for European discovery and assessment but does not claim that EHDS supplies permission to release it. [EHDS](https://eur-lex.europa.eu/eli/reg/2025/327/oj) and [2026 dataset-metadata rules](https://eur-lex.europa.eu/eli/reg_impl/2026/2098/oj).

HGNC and UniProt annotation will be an additive, versioned table. It will store release and retrieval date, mapping multiplicity and ambiguity, and will preserve every original row. A gene-to-many-proteins match will not be expanded into several invented quantitative observations. If governance permits eventual deposition, PRIDE/ProteomeXchange is the preferred European proteomics route, ideally with raw spectra and complete search outputs rather than only the reduced Excel matrix. [EMBL-EBI proteomics resources](https://www.ebi.ac.uk/about/teams/proteomics-metabolomics/).

Possible documentation extensions are a `v2` schema containing explicit unknowns, a governance-limitation statement, an EHDS-compatible catalogue description, and versioned identifier annotation. They are not prerequisites for matrix QC and cannot populate unavailable project facts. The original audit outputs remain frozen. The maintained protocol is in [Phase 1 plan and decisions](../01_PLAN_AND_DECISIONS.md), with contextual terminology in [European alignment and terminology](../EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md).

## Data-only roadmap request — 1 October 2026

> With the notion of working only with the data, give me the previous roadmaps and markdown reports again. Ensure every phase roadmap connects to the next phase, include all previous phases and milestones with required modifications, and track all limitations caused by not receiving new project-group information.

The structural Phase 1 milestone is now accepted without waiting for clinical, pathology, batch, acquisition or governance records. The existing manifest, feature audit, data dictionary and validation report hand stable sample and feature identities to Phase 2. Optional external annotation remains additive and cannot reconstruct original protein inference. The full sequence and hand-offs are in the [redesigned thesis-extension plan](../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md); inherited limitations and prohibited claims are in [data-only limitations and assumptions](../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md).
