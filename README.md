# Robust data-only analysis of paired oral-cancer tissue proteomics

This repository contains a reproducible, data-only analysis of an oral-cancer proteomics workbook comprising 84 specimens from 42 matched tumour/non-tumour patient pairs. The project integrates proteomics, bioinformatics, statistical inference, and leakage-controlled machine learning while preserving the patient as the independent biological unit.

The supplied workbook is the only source of participant-level measurements, labels, cohorts, and outcomes. Literature published through 31 December 2023 is used to justify methods and interpret pathways, not to manufacture missing laboratory or clinical metadata. European terminology and translational boundaries are documented separately without changing the analytical evidence cutoff.

## Research objectives

The analysis addresses four connected questions:

1. Which protein-group features show reproducible paired tumour/non-tumour abundance or detection differences?
2. Which Reactome pathways organise those feature-level signals?
3. How stable are the findings across preprocessing choices, patient influence analyses, and grouped resampling?
4. Which candidates are sufficiently robust for prospective validation planning, without being misrepresented as clinically validated biomarkers?

## Study design and scope

| Element | Definition |
|---|---|
| Design | 42 matched tumour/non-tumour patient pairs; 84 specimens |
| Feature unit | One immutable source protein-group row |
| Quantitative estimand | Within-patient log2 abundance difference, tumour minus matched non-tumour |
| Detection estimand | Within-patient detection difference, tumour minus matched non-tumour |
| Prediction target | Tissue state in patients wholly held out from model training |
| External information | Method selection and pathway annotation only |
| Evidence cutoff | 31 December 2023 for analytical methods |
| Intended use | Internal scientific discovery and prospective validation planning |

The repository does **not** establish an externally validated biomarker, a clinical diagnostic or screening test, prognosis, treatment response, causal pathway effects, or CE-mark/IVDR conformity. These boundaries are maintained in [DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md](DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md).

## Workflow status

| Stage | Principal deliverable | Status | Results |
|---:|---|---|---|
| 00 | Immutable data model, provenance, and scientific contract | Complete; M0 passed | `results/00_input_validation/` |
| 01 | Quality control, missingness atlas, and frozen cohorts | Complete; M1 passed | `results/01_qc_missingness/` |
| 02 | Missingness/imputation and PCA benchmark | Prespecified but not implemented | [02_MISSINGNESS_AND_PCA_EVIDENCE_2023.md](02_MISSINGNESS_AND_PCA_EVIDENCE_2023.md) |
| 03 | Paired abundance and detection inference | Computational branch complete; milestone conditional on Stage 02 | `results/03_paired_inference/` |
| 04 | Reactome interpretation and patient heterogeneity | Computational checks passed; milestone conditional | `results/04_pathways_heterogeneity/` |
| 05 | Leakage-controlled machine learning and stability analysis | Computational checks passed; no compact panel declared | `results/05_machine_learning/` |
| 06 | Multiverse, influence analysis, and integrated evidence synthesis | Implemented; independent clean rebuild remains outstanding | `results/06_robustness/` |
| 07 | Validation readiness and locked prospective hand-off | Readiness complete; external validation not executed | `results/07_validation_readiness/` |

Stage 02 is intentionally visible rather than silently skipped. All downstream conclusions inherit this unresolved preprocessing limitation.

## Repository organisation

```text
config/             Frozen analytical decisions and scientific contract
data/interim/       Regenerable patient-level data objects; ignored by Git
environment/        Pinned Python dependencies
overleaf_thesis/    Structured LaTeX thesis project
records/            Decision, implementation, and conversation audit records
resources/          Versioned Reactome v86 pathway membership resource
results/             Machine-readable outputs, summaries, and figures by stage
scripts/             Numbered orchestration and analysis entry points
src/oral_cancer/     Reusable statistical, bioinformatics, and ML modules
tests/               Automated contract and result tests
```

Top-level scientific documents provide the governing roadmap, evidence base, limitations register, implementation specification, and thesis manuscript. Numbered filenames identify the relevant workflow stage without embedding the word “phase” in paths.

## Data availability and governance

The source workbook is intentionally excluded from Git. To reproduce the analysis, place the authorised workbook at the repository root as:

```text
OC_Dataset_84Sample_v2_24012024.xlsx
```

The workflow expects worksheet `Sheet1`, identifier column `PG.Genes`, 8,071 feature rows, and sample identifiers matching `P<number>N` or `P<number>T`. Generated patient-level matrices under `data/interim/` are also excluded because they are reproducible from the source workbook and may contain sensitive specimen-level information.

Do not publish or redistribute patient-level material without confirmation of the relevant ethics, consent, controller, and data-sharing permissions. The current files do not establish GDPR/EHDS compliance or permission for cross-border release.

## Installation

Python 3.12 was used for the recorded environment. Create an isolated environment and install the pinned dependencies:

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r environment/requirements.txt
```

Activate the environment using the command appropriate to the operating system before running the workflow.

## Reproduction

Run the stages in numerical order. Each runner regenerates or validates its outputs and invokes the relevant automated checks.

```bash
python scripts/00_run.py
python scripts/01_run.py
# Stage 02 remains prespecified and unimplemented.
python scripts/03_run.py
python scripts/04_run.py
python scripts/05_run.py
python scripts/06_run.py
python scripts/07_run.py
```

Stage 05 performs repeated grouped nested validation and paired-label permutation analysis and is intentionally computationally expensive. Existing completed fits can be verified without recomputation using:

```bash
python scripts/05_run.py --reuse-existing
python scripts/06_run.py --reuse-existing
```

Run the complete automated test suite independently with:

```bash
python -m unittest discover -s tests -v
```

## Analytical safeguards

- All inferential comparisons respect matched patient pairs.
- Patient grouping is enforced in every predictive split.
- Feature selection, scaling, imputation where applicable, and model tuning occur within training folds.
- Missing abundance is not automatically interpreted as biological absence or zero.
- Quantitative abundance and detection evidence are modelled as separate information domains.
- False-discovery-rate control is applied within prespecified hypothesis families.
- Multiverse, leave-one-patient-out, negative-control, and redundancy analyses qualify candidate prioritisation.
- Results are presented as internal evidence and validation readiness, not clinical validity.

Detailed methods and citations are provided in [THESIS_REPORT.md](THESIS_REPORT.md), with reusable LaTeX sources in [overleaf_thesis/main.tex](overleaf_thesis/main.tex). The thesis project uses A4 KOMA-Script formatting, British English, BibLaTeX/Biber references, numbered chapters, tables, figures, and appendices.

## Reproducibility and provenance

The immutable analysis contract is stored in `config/analysis_contract.yml`. Stage-specific configurations are numerically prefixed under `config/`. Output manifests record artifact sizes and SHA-256 hashes; environment snapshots record the software stack; implementation records document decisions and known limitations.

The Reactome v86 GMT file is retained in extracted form with its provenance and checksum documented in `resources/reactome_v86/README.md`. Redundant archives, rendered PDFs, temporary files, caches, raw workbooks, and regenerable patient-level matrices are excluded from version control.

## Citation

No DOI or archival software release has yet been assigned. Until one is available, cite the repository URL, the exact Git commit used, the access date, and the accompanying thesis title. Any publication should also cite the primary methodological and Reactome sources listed in `overleaf_thesis/references.bib`.

## Licensing

No project-wide reuse licence has been declared. The absence of a licence means the code, text, and derived outputs are not automatically licensed for redistribution or modification. The bundled Reactome annotation resource retains the provenance and CC0 notice documented in its resource README.
