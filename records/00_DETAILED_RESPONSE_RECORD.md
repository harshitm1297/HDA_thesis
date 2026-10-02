# Phase 0 detailed response record

## Record date

1 October 2026

## User request

> What AI and data science techniques do you think are in the scope of this project? and for what purposes. post that, zoom in on phase 0. I will need your help to fill in the knowledge gap as the supervisor is not in touch anymore. for phase 0, fill the knowledge gaps with the most appropriate research available

## Full assistant response

The project has meaningful scope for AI and data science, but the strongest contribution will come from rigorous paired statistics, missingness-aware proteomics methods and interpretable machine learning—not from complex deep-learning models. With only 42 independent patients and more than 8,000 features, model complexity must remain tightly controlled.

A research-based Phase 0 reconstruction is stored in [00_RECONSTRUCTED_PROVENANCE.md](../00_RECONSTRUCTED_PROVENANCE.md).

## AI and data-science techniques in scope

| Technique | Purpose | Role in this project |
| --- | --- | --- |
| Reproducible data engineering | Convert the workbook into validated analytical objects | Essential |
| Statistical and visual QC | Detect failed runs, loading differences, missingness and possible batch structure | Essential |
| Missingness modelling | Separate abundance changes from condition-linked detection | Essential |
| Paired statistical modelling | Estimate tumor effects while controlling for patient | Primary analysis |
| Unsupervised learning | Explore sample structure and patient heterogeneity | Exploratory |
| Pathway and network analysis | Convert protein-level results into biological processes | Biological interpretation |
| Stability-based feature selection | Identify proteins repeatedly selected across resamples and preprocessing alternatives | Candidate prioritization |
| Regularized classification | Test whether a compact panel distinguishes the tissue states | Secondary analysis |
| Explainable AI | Describe which proteins drive an internal classifier | Secondary and descriptive |
| External-data integration | Check whether candidates replicate in other oral-cancer cohorts | Validation |
| Deep learning | Train nonlinear high-capacity models | Out of scope with 42 patients |

### Reproducible data engineering

The first data-science component is not modelling; it is converting the workbook into a trustworthy analytical dataset. We need automated checks for dimensions, duplicate sample names, pair completeness, duplicated protein identifiers, multi-gene protein groups, missingness and impossible values. The script should create immutable raw, annotated and analysis-ready objects.

This stage prevents manual spreadsheet operations from silently changing the analytical population. It also creates a complete audit trail from the original workbook to every result.

### Technical and sample-level anomaly detection

PCA, robust sample correlations, hierarchical clustering, detected-protein counts and intensity-distribution diagnostics are appropriate for identifying technically unusual samples. Robust distances or isolation-type anomaly scores may be used as supporting diagnostics, but they should never automatically exclude samples.

An anomaly score answers “does this sample look unlike the others?” It does not answer “is this sample technically invalid?” A sample could be unusual because of genuine tumor biology, inflammation, purity or oral subsite. Exclusion therefore requires corroborating evidence such as extremely poor coverage, a failed run, contamination or a sample-identity problem.

The modern proteomics QC literature recommends distinguishing system suitability controls, internal controls, external pooled controls and biological study samples. None of these controls is identifiable in the supplied matrix, so we can diagnose sample-level irregularities but cannot retrospectively certify LC-MS or preparation reproducibility. See the [quantitative-proteomics QC framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC11973981/).

### Missingness-aware modelling

Missing values are one of the main scientific signals and technical problems in this dataset. We should model at least two outcomes: the measured abundance when a protein is quantified and whether the protein is detected at all.

For abundance, the main analysis should preserve missing values and use sufficiently observed patient pairs. For detection, each pair can be classified as detected in both tissues, neither tissue, tumor only or non-tumor only. Exact paired detection tests or a proteomics hurdle model can then determine whether detection differs systematically by tissue.

This is preferable to treating every missing value as an unknown number that must be imputed. Proteomics benchmarks have found that unimputed differential analysis often performs as well as or better than commonly used imputation strategies. Single-value replacement methods perform particularly poorly, and KNN performance depends on the missingness mechanism. See the [proteomics imputation benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC10949645/).

A strong candidate method for sensitivity analysis is `msqrob2`, whose hurdle workflow combines abundance and detection information without requiring one imputed matrix to be treated as truth. See the [msqrob2 documentation](https://bioconductor.org/packages/release/bioc/html/msqrob2.html).

### Paired differential-abundance modelling

This should remain the statistical backbone. For each protein, the model should estimate the tissue effect after controlling for the patient:

```text
abundance ~ patient + tissue_status
```

On a verified log2 scale, the tissue coefficient represents the average within-patient tumor-minus-non-tumor change. Limma’s empirical-Bayes variance moderation is suitable because it stabilizes variance estimates across thousands of proteins and explicitly supports paired designs. See the [Limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

The results should include effect estimates, confidence intervals, false-discovery-adjusted p-values, usable-pair counts, detection rates and patient-level direction consistency. This is more informative than calling a protein significant solely because an unadjusted p-value crosses 0.05.

### Unsupervised learning and molecular heterogeneity

PCA and hierarchical clustering are clearly in scope. They can reveal whether tumor samples shift in a common direction, whether patient identity dominates variation and whether subsets of patients show different paired responses.

Consensus clustering, non-negative matrix factorization or other subtype-discovery methods may be explored only after strong feature reduction. With 42 patients, any clusters will be hypothesis-generating. Cluster stability should be tested by resampling, and clusters should not be presented as established oral-cancer subtypes.

UMAP and t-SNE may help visual exploration, but they should not be used as evidence of discrete biological groups. Their apparent clusters can change with parameters and initialization.

### Pathway, network and knowledge-based analysis

Ranked pathway enrichment is strongly in scope because biological changes are often distributed across several related proteins rather than concentrated in a few dramatic markers. We should examine Gene Ontology, Reactome, KEGG, protein complexes, subcellular localization and protein-interaction networks.

Network propagation or knowledge-graph-based prioritization can be used to identify biologically connected candidate modules. These methods should prioritize interpretation and external validation; they should not convert a nonsignificant protein into a proven biomarker simply because it connects to known cancer proteins.

### Stability-based feature selection

Feature selection is necessary because thousands of proteins cannot be used directly in a model trained on 42 patients. Selection should be repeated inside resampling rather than performed once on the full dataset.

Appropriate methods include elastic-net regularization, stability selection, repeated univariate filtering inside training folds, selection-frequency and coefficient-sign analysis, and sparse partial least squares as an exploratory comparator.

A useful protein panel should contain features that recur across folds, retain the same direction and remain selected under reasonable preprocessing alternatives. Selection stability is more informative here than a single “top 20” list.

### Regularized tissue classification

Elastic-net logistic regression is the best primary classifier because it handles correlated predictors, performs shrinkage and can produce a compact panel. A linear SVM is a reasonable comparator.

Both samples from a patient must always stay in the same resampling fold. Otherwise, the model may partially recognize a held-out tissue using the same patient’s other tissue in training. Filtering, normalization, imputation, scaling, feature selection and tuning must all occur inside the training fold.

Repeated patient-grouped nested cross-validation is the appropriate framework. Leave-one-patient-out evaluation can be used as a sensitivity analysis. Data leakage is especially damaging in small biomedical datasets, and feature selection before cross-validation can substantially inflate performance. See [nested cross-validation for high-dimensional data](https://pmc.ncbi.nlm.nih.gov/articles/PMC10125905/).

The model remains an exploratory tissue classifier. Forty-two patients are far below what would normally be expected for reliable clinical prediction-model development with several candidate predictors. Contemporary sample-size guidance stresses controlling overfitting and precisely estimating predicted probabilities, not relying on a simple events-per-variable rule. See [clinical prediction-model sample-size guidance](https://www.bmj.com/content/368/bmj.m441).

### Explainable AI

Coefficient direction, selection frequency, partial dependence and SHAP-style explanations may be used to describe a fitted model. For elastic net, standardized coefficients and their resampling stability will usually be clearer than SHAP.

These explanations are not causal. Correlated proteins may substitute for each other across resamples, so a low importance score does not mean a protein is biologically irrelevant, and a high score does not prove that it drives carcinogenesis.

Any model report should follow TRIPOD+AI and be assessed using PROBAST+AI, particularly for participant selection, predictor handling, outcome definition, sample size, overfitting and evaluation. See [TRIPOD+AI](https://www.tripod-statement.org/) and [PROBAST+AI](https://pmc.ncbi.nlm.nih.gov/articles/11931409/).

### Techniques outside the current scope

Deep neural networks, autoencoders, transformers and large ensemble models are not justified as primary methods. They have far more capacity than 42 patients can support and could produce visually impressive but unstable results.

The dataset also cannot currently support prognosis or survival prediction, stage or node-status prediction, treatment-response modelling, causal inference about cancer mechanisms, early-detection claims, saliva or blood biomarker performance, population-level clinical decision support, or confirmed molecular subtype discovery.

## Phase 0: best-available reconstruction

Phase 0 cannot be completed by filling undocumented experimental details with typical values from other studies. Literature can tell us what should have been documented and what analytical defaults are safest, but it cannot reveal which instrument or search settings generated this particular matrix.

The appropriate solution is an evidence hierarchy: verified, high-confidence inference, moderate inference and unresolved.

### Verified conclusions

The biological design comprises 42 independent patients and 84 specimens, with one tumor and one matched non-tumor sample per patient. The thesis sometimes describes these as 84 patients, but the column names and pairing show that these are 84 specimens from 42 patients.

The workbook contains 8,071 feature rows, 35.456% missing abundance cells, 114 entirely unmeasured rows, four duplicated gene labels and 48 multi-gene groups.

The values are positive linear-scale quantities, not log2 measurements. They range from approximately 0.0023 to 7.31 billion, and the report describes log2 transformation as a later analysis step.

The correct specimen terminology is “matched non-tumor tissue.” Exact distance from the tumor, pathological normality and evidence of field change are unavailable, so “healthy control” should not be used.

### High-confidence inference

The column header `PG.Genes` is documented Spectronaut terminology for genes associated with protein groups. This strongly indicates that the workbook originated from a Spectronaut report and was subsequently reduced to one identifier column and 84 sample-quantity columns.

The rows should therefore be treated as protein-group or gene-labelled quantitative features—not as independently verified individual proteins and not as gene-expression measurements. Semicolon-delimited identifiers are likely protein groups that could not be uniquely assigned to one gene.

The official Spectronaut documentation states that the meaning of `PG.Quantity` depends on the configured analysis settings. Consequently, the stripped matrix does not tell us whether the original quantity was total, MS1-based, MS2-based or otherwise configured. See the [Spectronaut manual](https://biognosys.com/content/uploads/2024/09/Spectronaut-19-manual-v4.pdf).

### Probable but unverified conclusion

The leading hypothesis is a label-free LC-MS/MS DIA workflow processed in Spectronaut. This is consistent with the software fingerprint, continuous sample-level quantities, high proteome coverage and extensive missingness.

However, the acquisition mode cannot be promoted to a verified fact. The correct methods wording is:

> The supplied matrix was likely derived from a Spectronaut-based label-free LC-MS/MS workflow, probably using data-independent acquisition. The original acquisition and processing records were unavailable.

We also cannot determine whether Spectronaut normalization had already been applied. The thesis’s later median normalization does not prove that the source values were unnormalized. Applying another normalization automatically could therefore produce double normalization.

### Genuinely unresolved information

The following details cannot be reconstructed from the abundance matrix or general literature:

- Instrument manufacturer and model
- LC system, column and gradient
- Tissue preservation and extraction protocol
- Reduction, alkylation and digestion settings
- Injection quantities and acquisition order
- DIA window design or DDA settings
- Spectronaut version
- Library-based, directDIA or hybrid workflow
- Quantity field used in the export
- Spectronaut normalization settings
- Protein and precursor q-value thresholds
- Protein-inference parameters
- Sequence database and release
- Stable protein accessions and peptide evidence
- Batch assignments and technical replicates
- Internal standards or pooled QC samples
- Clinical, pathology and exposure metadata
- Ethics and consent identifiers

MIAPE-MS and HUPO standards consider this information part of the minimum record needed to reproduce and evaluate a proteomics experiment. See [MIAPE-MS](https://www.nature.com/articles/nbt0808-860) and the [HUPO MS Data Interpretation Guidelines](https://pmc.ncbi.nlm.nih.gov/articles/PMC6986310/).

No public ProteomeXchange or PRIDE accession was verified as a match for this exact 42-patient matrix. An experimental protocol from another oral-cancer study must therefore not be presented as this project’s protocol.

## Research-based operational defaults

The data will be described as a supplied linear-scale protein-group abundance matrix. It will not be called raw intensity, normalized intensity or confirmed DIA output unless stronger evidence appears.

`PG.Genes` will remain the original feature key. Multi-gene groups will remain grouped, and duplicated symbols will remain separate until accessions or peptide evidence allow resolution. Current UniProt/HGNC annotations can be added, but they cannot retroactively establish the original protein inference.

Log2 transformation is justified as an analytical candidate, but normalization will be benchmarked. A no-additional-normalization baseline will be compared with restrained sample-median centering. Quantile normalization will not be automatic, and batch correction will not be attempted without defensible batch labels.

The primary differential analysis will avoid imputation. Missingness will receive its own paired detection analysis, and `msqrob2` hurdle modelling can provide a proteomics-specific sensitivity analysis.

The protein-identification FDR is unverifiable. We can control the downstream statistical FDR across abundance and detection tests, but this is distinct from verifying that every protein group was originally identified at 1% FDR.

No latent component found by PCA will automatically be removed. Without batch metadata, a component may represent biology, technical variation or both.

All biological claims will be limited to paired tissue discovery. Individual protein identities with ambiguous groups, isoform-level claims, diagnostic claims and clinical utility require stronger evidence.

## Immediate Phase 0 completion criteria

Phase 0 is operationally complete when the reconstructed provenance is accepted as version 1 and linked to a sample manifest, data dictionary and identifier audit.

The source workbook and report are fingerprinted with SHA-256 checksums in the reconstructed provenance document. If an old Spectronaut export, `.sne` project, `.tsv` report, raw instrument file, email attachment, laboratory notebook or ethics document surfaces later, the provenance can be updated without changing the existing audit trail.

This is enough to begin a transparent secondary analysis. It is not enough to claim full experimental reproducibility, but it provides an honest and scientifically defensible foundation.

## Associated source records

- [Redesigned thesis-extension and data-science/AI plan](../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md)
- [Phase 0 reconstructed provenance](../00_RECONSTRUCTED_PROVENANCE.md)
- [Phase 0 milestone record](00_MILESTONE_RECORD.md)

## European alignment request — 1 October 2026

> Can you ensure that the research, ideas and terminologies are also closely aligned with what has happened in Europe till now? Also, make the markdown reports so far, more detailed, including what, why and how (without these headlines)

The Phase 0 record was reviewed against European sources current to 1 October 2026. The revision changes the working vocabulary and the information-recovery protocol rather than adding a superficial Europe bibliography.

The cohort continues to be described as oral cancer because the available files do not verify the exact anatomical site and morphology for every person. Under the 2024 WHO/IARC fifth-edition classification, oral cavity/mobile tongue and oropharynx/base of tongue are separate groupings. OCSCC will therefore be used only after verified pathology mapping, and HNSCC will remain a broad comparison category rather than a synonym. The non-tumour specimen remains “matched non-tumour tissue”; it is not a healthy control unless distance, pathology and collection evidence establish that description. [WHO/IARC classification](https://publications.iarc.who.int/Book-And-Report-Series/Who-Classification-Of-Tumours/Head-And-Neck-Tumours-2024) and [EHNS–ESMO–ESTRO guideline](https://www.annalsofoncology.org/article/S0923-7534%2820%2939949-X/fulltext).

The metadata recovery list now includes exact site, original pathology wording, mapped classification and edition, depth of invasion, margins, perineural and lymphovascular invasion, nodal findings, extranodal extension, HPV method, oral potentially malignant disorders, and the pathology of the comparison tissue. Stage will be stored as cTNM and pTNM components, stage group, assessment date, source and edition. This is required because TNM 9 took effect on 1 January 2026 and must not be silently imposed on an older cohort. [UICC TNM resources](https://www.uicc.org/resources/tnm/publications-resources).

European transportability will be evaluated using exact site, morphology, specimen, geography, centre, exposure pattern, stage, treatment and platform. Smoking and alcohol fields are supplemented—not substituted—for smokeless tobacco, areca/betel quid and oral submucous fibrosis. This preserves the likely source-population context while making later European comparisons interpretable. ECIS site categories and Europe's cancer-inequalities work support stratified, population-aware comparison rather than a pooled head-and-neck claim. [ECIS](https://ecis.jrc.ec.europa.eu/database-description) and [Europe's Beating Cancer Plan](https://health.ec.europa.eu/non-communicable-diseases/cancer_en).

The Phase 0 inventory now includes ethics and data governance. Coded patient IDs are treated as pseudonymous, not automatically anonymous. Ethics approval, consent scope or another lawful basis, controller roles, permitted secondary use, retention, access, re-identification-key controls, publication and international transfer must be recovered before wider sharing. This implements the GDPR research framework rather than assuming that removal of names is sufficient. [GDPR](https://eur-lex.europa.eu/eli/reg/2016/679/oj) and [EDPB pseudonymisation guidance](https://www.edpb.europa.eu/system/files/2025-01/edpb_guidelines_202501_pseudonymisation_en.pdf).

The EHDS explicitly includes proteomic data and now has a September 2026 implementing regulation for dataset-description metadata. Because EHDS applies progressively, it is used to design a future-compatible catalogue, not as a present data permit. The catalogue will state population, purpose, geography and period, data categories, coding systems, quality, limitations, access conditions, provenance and responsible holder. [EHDS](https://eur-lex.europa.eu/eli/reg/2025/327/oj) and [Implementing Regulation (EU) 2026/2098](https://eur-lex.europa.eu/eli/reg_impl/2026/2098/oj).

The exploratory classifier remains scientific research software. A later medical intended purpose can move a protein assay or classifier into the IVDR/MDR and can trigger AI Act high-risk requirements where third-party conformity assessment applies. The project will preserve intended use, provenance, versioning, validation separation, limitations and human-oversight assumptions now, while making no claim of clinical performance, CE marking or conformity. [EU AI Act](https://eur-lex.europa.eu/eli/reg/2024/1689/oj), [IVDR](https://eur-lex.europa.eu/eli/reg/2017/746/oj), and [AIB/MDCG guidance](https://health.ec.europa.eu/document/download/b78a17d7-e3cd-4943-851d-e02a2f22bbb4_en?filename=mdcg_2025-6_en.pdf).

The complete operational treatment is in [European alignment and terminology](../EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md), with supporting evidence mapped in [Research evidence base](../RESEARCH_EVIDENCE_BASE.md).

## Data-only roadmap request — 1 October 2026

> With the notion of working only with the data, give me the previous roadmaps and markdown reports again. Ensure every phase roadmap connects to the next phase, include all previous phases and milestones with required modifications, and track all limitations caused by not receiving new project-group information.

The project now treats lack of supervisor and laboratory contact as a fixed design condition. Phase 0 closes by bounding uncertainty rather than recovering full experimental provenance. The governing documents are the [redesigned thesis-extension plan](../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md), [data-only limitations and assumptions](../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md), and [research evidence base through 2023](../RESEARCH_EVIDENCE_BASE.md). Phase 0 hands its claim boundary, evidence-state rules and limitation IDs to Phase 1; every later phase must inherit them.

