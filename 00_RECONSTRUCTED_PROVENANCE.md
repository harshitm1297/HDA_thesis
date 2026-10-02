# Phase 0 reconstructed provenance

> Data-only status: no new laboratory, clinical or project-group information is expected. The unresolved fields in this document are permanent limitations of the present analysis and are tracked by ID in `DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`; they are not prerequisites that block matrix-based phases.

## Purpose

This document reconstructs the experimental and computational provenance of the oral-cancer proteomics matrix using the surviving workbook, the 2024 project report, embedded file metadata, public documentation, and current proteomics methods literature. It distinguishes verified facts, evidence-based inferences, unresolved facts, and analysis defaults adopted because the original laboratory record is unavailable.

No literature source can recover a project-specific instrument setting, batch assignment, or search parameter that was never recorded. Such fields remain unknown rather than being replaced by generic defaults.

## Surviving source files

### Quantitative workbook

- File: `OC_Dataset_84Sample_v2_24012024.xlsx`
- SHA-256: `65211B375E366D0CF6078E1016A3D29867EFC5E0EBA10D97F1F0282B25161DA8`
- Embedded creation time: 2024-01-24 20:47:48 UTC
- Embedded modification time: 2024-01-27 18:54:54 UTC
- Embedded creator: `Admin`
- Embedded last modifier: `Anjali Kanojia`
- One worksheet named `Sheet1`

### Project report

- File: `MSc_Project__Sem_4___Updated_.pdf`
- SHA-256: `3D09CF4D6CD4BE44D2B7CF882DEAF8D41EADEAE96F862474243B61CA9D57F612`
- Report date: 22 April 2024

## Reconstructed experimental record

| Field | Best available conclusion | Confidence | Evidence and interpretation |
| --- | --- | ---: | --- |
| Biological study design | Matched tumor and non-tumor tissue from 42 patients with oral cancer | Verified | The matrix contains `P1` through `P42`, each with `T` and `N`; the report describes paired tissue from the same patients. There are 84 specimens, not 84 independent patients. |
| Independent biological unit | Patient | Verified | Each patient contributes two correlated specimens. All inference and resampling must respect the pair. |
| Organism | Human | High | The report describes patients and the identifiers are human gene symbols. |
| Tissue label meaning | Tumor tissue versus matched non-tumor tissue | Verified | Exact distance from the tumor and pathological confirmation of the non-tumor tissue are not recorded. Use “matched non-tumor,” not “healthy control.” |
| Disease subtype | Oral cancer; OSCC is plausible but not demonstrated for every patient | Moderate | The report uses oral cancer terminology. Histology fields are absent. Do not relabel the cohort as uniformly OSCC until pathology is recovered. |
| European-compatible anatomical entity | Not established | Unresolved | The surviving files do not distinguish mobile tongue and other oral-cavity subsites from base of tongue/oropharynx. WHO/IARC 2024 treats these as separate anatomical groupings, so a broad label cannot be translated into a precise entity. |
| Stage and staging edition | Not available | Unresolved | No cTNM, pTNM, stage group, assessment date, or TNM edition survives. TNM 9 took effect in 2026, but it must not be assigned retrospectively to the 2024 cohort without source components. |
| Matrix dimensions | 8,071 feature rows by 84 sample columns plus one identifier column | Verified | Direct workbook audit. |
| Matrix feature level | Protein-group/gene-labelled quantitative features | High | The identifier header is `PG.Genes`, where `PG` is Spectronaut protein-group terminology. Forty-eight rows contain multiple semicolon-delimited genes. |
| Upstream software family | Spectronaut-derived export is highly likely | High inference | `PG.Genes` is a documented Spectronaut protein-group report header. The surviving table has probably been reduced from a richer report, because accessions, protein groups, q-values and quantity-field names are absent. |
| Acquisition strategy | Label-free LC-MS/MS with DIA is the leading hypothesis | Moderate inference | Spectronaut is primarily used for DIA analysis, and the matrix consists of run/sample-level continuous quantities. The workbook does not prove DIA versus a different Spectronaut-supported workflow. |
| Quantification level | Protein-group quantity, probably pivoted to one column per sample | Moderate inference | The table structure is consistent with a wide protein-group quantity export, but the original field could have been `PG.Quantity`, `PG.MS1Quantity`, `PG.MS2Quantity`, or another configured quantity. |
| Numerical scale | Positive linear-scale quantities, not log2 values | High | Values span approximately 0.0023 to 7.31 billion, and the report states that log2 transformation was applied later. |
| Upstream normalization | Unknown | Unresolved | Spectronaut quantity is defined by analysis settings, and the reduced workbook does not preserve those settings. The report’s later median normalization does not establish whether Spectronaut normalization was already active. |
| Missing-value representation | Blank cells interpreted as unreported/unquantified values | Verified for the matrix | The abundance block is 35.456% missing. The reason for each missing value—below detection, identification filtering, interference, failed transfer or other cause—is not recoverable from this table. |
| Protein identification FDR | Unknown and not independently verifiable | Unresolved | The reduced matrix has no protein or precursor q-values. A conventional 1% FDR may have been used upstream, but it must not be asserted as a project fact. |
| Protein inference | Spectronaut/search-engine protein grouping is likely, exact rules unknown | Moderate inference | `PG.Genes` and multi-gene rows are consistent with protein groups. The original `PG.ProteinGroups` and `PG.ProteinAccessions` fields are missing. |
| Sequence database | Unknown | Unresolved | No database name, release date, contaminant database, isoform policy or parsing rule survives. |
| Sample preparation | Unknown | Unresolved | Tissue preservation, lysis, protein extraction, reduction, alkylation, digestion, cleanup and loading quantities are not documented in the available report. |
| Instrument and chromatography | Unknown | Unresolved | Instrument model, LC system, column, gradient, ion source and acquisition windows cannot be reconstructed from abundance values. |
| Spectral-library strategy | Unknown | Unresolved | Library-based DIA, directDIA and hybrid strategies cannot be distinguished from this matrix. |
| Batch and run order | Unknown | Unresolved | Spreadsheet column order must not be treated as acquisition order. No batch correction will be applied without defensible batch labels. |
| Technical replicates and QC injections | None identifiable in the workbook | Low-to-moderate | All 84 columns follow patient/tissue naming. Separate QC or technical-replicate files may once have existed. |
| Clinical covariates | Not available | Verified absence from supplied files | Subsite, histology, stage, grade, node status, age, sex, exposure history, HPV, tumor purity, inflammation and normal-tissue distance/pathology are not included. |
| Ethics, consent and data governance | Not available in the supplied analytical files | Unresolved | Ethics approval, consent scope or other lawful basis, controller/processor roles, retention, secondary-use permission, international transfer conditions, and public-deposition permission must be recovered from institutional records. Their absence from the workbook is not evidence that they do not exist. |
| Public data accession | No matching accession verified | Unresolved | Searches of publications and ProteomeXchange/PRIDE did not identify a definitive public match for this exact 42-patient matrix. |

## Research basis for the reconstruction

Spectronaut documentation defines `PG.Genes` as genes associated with protein groups and states that `PG.Quantity` is the quantitative value defined by the analysis settings. This supports identification of the software family but also explains why the exact quantity and normalization cannot be recovered from the stripped workbook alone.

- Spectronaut 19 manual: https://biognosys.com/content/uploads/2024/09/Spectronaut-19-manual-v4.pdf

MIAPE-MS and HUPO guidance identify instrument configuration, acquisition, identification criteria, protein evidence and data deposition as necessary provenance. These standards justify leaving missing project-specific items unresolved rather than filling them with typical laboratory settings.

- MIAPE-MS reporting guidance: https://doi.org/10.1038/nbt0808-860
- HUPO HPP MS Data Interpretation Guidelines 3.0: https://pmc.ncbi.nlm.nih.gov/articles/PMC6986310/
- ProteomeXchange submission guidance: https://www.proteomexchange.org/docs/guidelines_px.pdf

Recent QC guidance separates system suitability, internal controls, external pooled controls and study samples. Because none of these controls can be identified in the surviving table, sample-level statistical QC can be performed, but instrument performance and preparation reproducibility cannot be retrospectively certified.

- Quantitative proteomics QC framework: https://pmc.ncbi.nlm.nih.gov/articles/PMC11973981/

## European clinical and governance interpretation

The European-compatible disease description is narrower than the colloquial phrase “oral cancer.” WHO/IARC's 2024 fifth-edition classification has a dedicated oral cavity and mobile tongue section, separate from the oropharynx; it lists oral squamous cell carcinoma within the former and HPV-associated and HPV-independent squamous carcinomas within the latter. The missing subsite and histology fields therefore prevent a valid cohort-wide conversion to OCSCC or an HPV-defined category. The recovery method is documentary: obtain pathology reports or a verified clinical table, retain the original wording, add a mapped WHO/IARC entity and mapping confidence, and keep mobile tongue distinct from base of tongue. https://publications.iarc.who.int/Book-And-Report-Series/Who-Classification-Of-Tumours/Head-And-Neck-Tumours-2024

The EHNS–ESMO–ESTRO guideline identifies depth of invasion, resection margins, nodal burden and location, extranodal extension, perineural invasion, and lymphatic invasion as relevant oral-cavity pathological information. These variables should be recovered because they can explain tumour composition and proteomic heterogeneity. They must not be reverse-engineered from protein patterns. https://www.annalsofoncology.org/article/S0923-7534%2820%2939949-X/fulltext

Any recovered stage must include whether it is clinical or pathological, its component T/N/M values, stage group, assessment date, and edition. UICC recommends TNM 9 from 1 January 2026, whereas the project report predates that date. The historical classification remains the primary record; a harmonised TNM 9 field may be added only when sufficient source variables make conversion valid and the conversion rule is documented. https://www.uicc.org/resources/tnm/publications-resources

European data protection cannot be inferred from de-identification of visible names. The sample codes preserve repeated observations for a person and may be linkable through institutional records, so they should be treated as pseudonymous until a formal identifiability assessment proves otherwise. Recovery must therefore include the ethics committee and approval number, consent language or other lawful basis, controller, re-identification-key location, authorised users, retention, transfer and publication conditions. GDPR Articles 6, 9 and 89 provide the legal structure; institutional and applicable national rules determine the project-specific answer. https://eur-lex.europa.eu/eli/reg/2016/679/oj

The EHDS Regulation explicitly includes proteomic data, research cohorts, and biobank-linked data in its secondary-use categories. It applies progressively and does not retroactively cure an unknown lawful basis. Its immediate analytical use is to structure a dataset catalogue and quality description, supported by the September 2026 implementing metadata regulation. https://eur-lex.europa.eu/eli/reg/2025/327/oj and https://eur-lex.europa.eu/eli/reg_impl/2026/2098/oj

## Operational defaults for analysis

### Measurement statement

Until stronger evidence is recovered, describe the input as:

> A supplied matrix of positive, linear-scale, protein-group/gene-labelled abundance values for 42 matched oral-tumor and non-tumor tissue pairs, likely derived from a Spectronaut-based label-free LC-MS/MS workflow. Exact acquisition, quantity definition, upstream normalization, protein-inference settings and identification FDR were unavailable.

Do not shorten this to “raw protein intensities” or “DIA protein abundances” without preserving the uncertainty.

### Identifier policy

Preserve `PG.Genes` as the original feature identifier. Do not split semicolon-delimited groups into independent proteins and do not average duplicate symbols. Add current HGNC/UniProt annotations where a mapping is unambiguous, but keep the original row as the analytical unit. Claims about a specific isoform or protein accession are out of scope unless accessions or peptide evidence are recovered.

### Transformation and normalization policy

The main transformation candidate is log2 because the values are positive and clearly not already on a log2 scale. Upstream normalization is unknown, so the analysis will compare a no-additional-normalization baseline with a restrained sample-median-centering alternative. The selected primary method must improve technical comparability without erasing paired biological changes. Quantile normalization or batch correction will not be applied automatically.

### Missingness policy

Do not assume missing completely at random and do not make class-wise KNN imputation the primary dataset. Primary analyses will preserve missing values, report complete-pair counts and separately model condition-linked detection. A hurdle-model sensitivity analysis can combine abundance and detection evidence without requiring a single imputed ground truth. Proteomics benchmarking has shown that unimputed differential analysis often performs as well as or better than common imputation strategies, while single-value methods perform especially poorly.

- Imputation benchmark: https://pmc.ncbi.nlm.nih.gov/articles/PMC10949645/
- `msqrob2` robust and hurdle-model workflows: https://bioconductor.org/packages/release/bioc/html/msqrob2.html

### Statistical design policy

Patient is the blocking factor. The primary abundance model will estimate the tissue effect after controlling for patient. A moderated paired model and a proteomics-specific robust/hurdle analysis will be compared. All large-scale tests will receive false-discovery correction, and results will include effects, uncertainty, usable-pair counts, detection rates and patient-level direction consistency.

- Limma paired-design guidance: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

### Batch policy

Do not infer batch from spreadsheet column order. Without batch or run-order metadata, no ComBat or other batch-removal method will be used. PCA and other QC can identify suspicious structure, but latent components will not automatically be removed because they may represent the tumor effect or real patient heterogeneity.

### Claim policy

The dataset can support paired biological discovery and cautious internal tissue classification. It cannot independently support claims about early diagnosis, prognosis, treatment response, population screening, saliva/blood performance or clinical utility. Protein-identification confidence cannot be revalidated from gene labels alone.

For European-facing work, “candidate biomarker” means a discovery-stage molecular feature worthy of independent validation; it does not mean an IVDR-validated test. “Internal validation” means resampling within this cohort, not validation in another centre or population. “External validation” requires a genuinely independent dataset and a prespecified analysis matched on intended use, anatomical site, morphology, specimen and comparison group. A future diagnostic assay or medical-purpose classifier would require separate analytical and clinical performance work under the applicable IVDR/MDR and AI Act framework. https://eur-lex.europa.eu/eli/reg/2017/746/oj and https://eur-lex.europa.eu/eli/reg/2024/1689/oj

## Information unavailable under the data-only constraint

The following cannot be established from the supplied analytical files. Current literature can explain why they matter and guide robust alternatives, but it cannot supply their project-specific values:

- Instrument model and acquisition settings
- LC gradient and injection order
- Tissue preservation and sample-preparation protocol
- Spectronaut version and complete analysis settings
- Exact quantity field and normalization configuration
- Spectral-library construction
- Protein/peptide FDR and filtering settings
- Sequence database and release
- Stable protein accessions and peptide evidence
- Batch assignments and technical QC results
- Clinical, pathological and exposure metadata
- Ethics approval and consent identifiers
- Data controller, permitted secondary uses, access level, retention, and transfer conditions
- Exact oral subsite, WHO/IARC pathological entity, TNM components, and TNM edition

These gaps do not prevent a transparent secondary analysis of the supplied matrix. They prevent full experimental reproducibility and permanently limit protein-identification, clinical, batch-adjusted, sharing, screening and diagnostic claims. Their mitigations and downstream propagation are recorded in `DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`.
