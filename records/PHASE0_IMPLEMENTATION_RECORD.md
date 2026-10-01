# Phase 0 implementation record

## Repository change

The project directory was converted into a local Git repository on branch `main`. Raw workbooks, PDFs, temporary Office files, generated patient-level matrices, Python caches, and the superseded `phase1_outputs` directory are excluded from version control. The generated matrices remain locally available and are reproducible from the source workbook.

## Executable workflow

`scripts/00_validate_input.py` validates the immutable workbook and writes the sample manifest, feature manifest, and machine-readable validation result. `scripts/01_build_data_model.py` creates the raw, log2-paired, and detection objects plus their metadata and data dictionary. `scripts/run_phase0.py` executes both steps, runs the automated tests, and writes an output checksum manifest. Configuration and scientific restrictions are versioned in `config/analysis.yml` and `config/analysis_contract.yml`.

The implementation deliberately does not inherit biological candidates, thresholds, exclusions, or conclusions from the pre-existing thesis. It does reuse directly verifiable facts from the workbook, such as its dimensions and pair structure. The source field `PG.Genes` remains an annotation on a quantitative proteomics row; it is not reinterpreted as a transcriptomic measurement.

## Reproduction command

```powershell
python scripts/run_phase0.py
```

The tested local runtime was Python 3.12.14 with NumPy 2.3.5, pandas 3.0.1, and openpyxl 3.1.5. Exact package pins are recorded in `environment/requirements.txt`.

## Recorded outcome

The run completed successfully, all four tests passed, and M0 passed. The definitive quantitative findings and their interpretation are recorded in `results/phase0/PHASE0_RESULTS.md`; raw machine-readable checks are in `results/phase0/input_validation.json`, and artifact checksums are in `results/phase0/output_manifest.json`.

