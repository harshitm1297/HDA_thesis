# Research evidence base through 2023

This evidence register supports the redesigned oral-cancer proteomics project. Its methodological cutoff is **31 December 2023**: later publications are not used to justify analytical choices. The sources help choose defensible computational methods; they do not supply missing facts about this particular cohort, laboratory workflow, or data export. Those unresolved facts remain explicit in `DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`.

## Reproducibility, provenance, and proteomics reporting

The Minimum Information About a Proteomics Experiment (MIAPE) framework defines the information required to interpret and reproduce proteomics studies. The mass-spectrometry module makes instrument acquisition and processing provenance especially relevant. HUPO-PSI reporting guidance and the FAIR principles further support preserving identifiers, metadata, processing decisions, and machine-readable outputs. In this project, those standards justify the Phase 0 data dictionary, patient-pair checks, immutable raw-data fingerprint, decision log, and explicit separation between observed metadata and inferred assumptions.

- Taylor CF et al. The minimum information about a proteomics experiment (MIAPE). *Nature Biotechnology*. 2007;25:887–893. https://doi.org/10.1038/nbt1329
- Taylor CF et al. Guidelines for reporting the use of mass spectrometry in proteomics. *Nature Biotechnology*. 2008;26:860–861. https://doi.org/10.1038/nbt0808-860
- Deutsch EW et al. Proteomics Standards Initiative: fifteen years of progress and future work. *Journal of Proteome Research*. 2019;18:2634–2635. https://doi.org/10.1021/acs.jproteome.9b00199
- Wilkinson MD et al. The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data*. 2016;3:160018. https://doi.org/10.1038/sdata.2016.18

## Missingness and imputation in mass-spectrometry proteomics

Proteomics missingness can reflect both abundance-dependent non-detection and other technical or stochastic mechanisms. Treating every missing value as missing at random and selecting a single global imputer can distort differential-abundance estimates. The project therefore profiles missingness by tissue, protein, patient pair, and intensity; separates detection indicators from quantitative abundance; benchmarks multiple preprocessing strategies; and treats imputation as an analytical choice whose uncertainty and sensitivity must be reported. It does not assert a laboratory cause that cannot be established from the supplied data.

- Webb-Robertson B-JM et al. Review, evaluation, and discussion of the challenges of missing value imputation for mass spectrometry-based label-free global proteomics. *Journal of Proteome Research*. 2015;14:1993–2001. https://doi.org/10.1021/pr501138h
- Lazar C et al. Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets to compare imputation strategies. *Journal of Proteome Research*. 2016;15:1116–1125. https://doi.org/10.1021/acs.jproteome.5b00981
- O'Brien JJ et al. The effects of nonignorable missing data on label-free mass spectrometry proteomics experiments. *Annals of Applied Statistics*. 2018;12:2075–2095. https://doi.org/10.1214/18-AOAS1144
- Chion M, Carapito C, Bertrand F. Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics. *PLoS Computational Biology*. 2022;18:e1010420. https://doi.org/10.1371/journal.pcbi.1010420

## Paired differential-abundance inference and multiplicity

The tumour and adjacent-tissue samples are paired within patients, so the inferential unit is the patient and the pairing must be preserved. Thousands of protein-wise tests also require false-discovery-rate control. Empirical-Bayes linear modelling is appropriate for small, high-dimensional studies because it stabilises variance estimates; robust proteomics models are useful sensitivity analyses. Results will therefore report effect sizes and confidence intervals alongside Benjamini–Hochberg-adjusted values, rather than treating nominal *p* < 0.05 as sufficient evidence.

- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Ritchie ME et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43:e47. https://doi.org/10.1093/nar/gkv007
- Goeminne LJE et al. Experimental design and data-analysis in label-free quantitative LC/MS proteomics: a tutorial with MSqRob. *Journal of Proteomics*. 2018;171:23–36. https://doi.org/10.1016/j.jprot.2017.04.004

## Pathway and systems-level interpretation

Protein-by-protein findings alone are difficult to interpret and unstable in a small cohort. Ranked-set enrichment avoids requiring an arbitrary significance cutoff, while a versioned pathway resource makes the mapping reproducible. Pathway results remain exploratory because this dataset cannot establish mechanism or causality.

- Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*. 2005;102:15545–15550. https://doi.org/10.1073/pnas.0506580102
- Gillespie M et al. The Reactome pathway knowledgebase 2022. *Nucleic Acids Research*. 2022;50:D687–D692. https://doi.org/10.1093/nar/gkab1028

## Sparse prediction and feature stability

The feature count is much larger than the patient count, so an unconstrained classifier would be unstable. Elastic-net logistic regression supplies shrinkage and sparse selection while accommodating correlated proteins. Stability selection motivates reporting how often a feature is selected across resamples rather than presenting one fitted coefficient list as a discovered biomarker panel. The project uses these methods for a proof-of-concept tissue classifier, not a clinical diagnostic claim.

- Zou H, Hastie T. Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society: Series B*. 2005;67:301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
- Meinshausen N, Bühlmann P. Stability selection. *Journal of the Royal Statistical Society: Series B*. 2010;72:417–473. https://doi.org/10.1111/j.1467-9868.2010.00740.x

## Nested validation, patient grouping, and leakage prevention

Model selection performed on the same resamples used for performance estimation produces optimistic results. In addition, fitting filters, imputers, scalers, feature selectors, or batch corrections before cross-validation leaks information from held-out observations. All learned preprocessing and tuning will therefore occur inside repeated nested cross-validation, with both samples from a patient kept in the same split. A within-patient label-permutation test will assess whether observed performance exceeds a null pipeline run under the same selection process.

- Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*. 2006;7:91. https://doi.org/10.1186/1471-2105-7-91
- Lewis JE et al. nestedcv: an R package for fast implementation of nested cross-validation with embedded feature selection designed for transcriptomics and high-dimensional data. *Bioinformatics Advances*. 2023;3:vbad048. https://pubmed.ncbi.nlm.nih.gov/37113250/
- Davis SE et al. A framework for understanding label leakage in machine learning for health care. *Journal of the American Medical Informatics Association*. Published online 2023;31:274–280. https://doi.org/10.1093/jamia/ocad178
- Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns*. 2023;4:100804. https://doi.org/10.1016/j.patter.2023.100804

## Transparent prediction reporting and risk-of-bias assessment

TRIPOD provides the reporting structure for multivariable prediction studies, and PROBAST exposes common sources of bias in participants, predictors, outcomes, and analysis. Their use here does not make the study a clinical validation. Instead, they constrain the language, require complete reporting of resampling and model construction, and make the absence of an external cohort explicit.

- Collins GS et al. Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis (TRIPOD): the TRIPOD statement. *BMJ*. 2015;350:g7594. https://doi.org/10.1136/bmj.g7594
- Wolff RF et al. PROBAST: a tool to assess the risk of bias and applicability of prediction model studies. *Annals of Internal Medicine*. 2019;170:51–58. https://doi.org/10.7326/M18-1376

## Interpretation rule for this repository

Published research determines the set of defensible methods and sensitivity analyses. It does **not** justify inventing cohort attributes, instrument settings, treatment histories, tumour stages, pathology details, batch labels, or sample-processing facts. Where such information is absent, analyses must either use observable proxies with a stated limitation, test a bounded sensitivity scenario, or omit the claim. Every such decision is linked to the limitation and assumption identifiers in `DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`.
