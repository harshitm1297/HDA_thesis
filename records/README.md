# Phase and milestone record index

## Governing position

The project proceeds without expecting new information from the former supervisor, laboratory, clinical team, or project group. `../FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md` is the single governing Phase 0–6 roadmap. It starts from the supplied workbook and first principles, without inheriting previous analytical decisions. Its analytical evidence cutoff is 31 December 2023. All patient-level calculations use only the supplied workbook; external cohorts are excluded.

`../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md` is the living register of missing information, mitigations, and prohibited claims. `../RESEARCH_EVIDENCE_BASE.md` gives the supporting methods literature through 2023. `../EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md` is a separate current-context note for European terminology, governance, and translational boundaries; it does not alter the 2023 methodological cutoff.

## Record types

The retained records have distinct purposes:

1. Files ending in `_RECORD.md` preserve detailed user-visible phase discussions and the requests that prompted them. They are conversation archives, not necessarily the current protocol.
2. Milestone and work-log files provide a concise chronological audit of decisions and completed work.
3. Technical phase documents define the maintained protocol, outputs, evidence, and acceptance criteria.
4. Machine-generated validation evidence is stored beside the analytical outputs.

These files do not contain hidden model reasoning or raw tool logs. When an archived response conflicts with the governing roadmap, the governing roadmap and current limitations register control execution.

## Phase 0 records

- `../00_RECONSTRUCTED_PROVENANCE.md`: evidence-based reconstruction of the source experiment and unresolved facts.
- `00_DETAILED_RESPONSE_RECORD.md`: detailed user-visible Phase 0 response archive.
- `00_MILESTONE_RECORD.md`: concise Phase 0 decisions and milestone status.

## Phase 1 records

- `../01_PLAN_AND_DECISIONS.md`: maintained import, pairing, identifier-audit, and provenance protocol.
- `01_DETAILED_RESPONSE_RECORD.md`: detailed user-visible Phase 1 response archive.
- `01_WORK_LOG.md`: chronological Phase 1 work log.
- `../results/01_qc_missingness/01_validation.json`: maintained machine-readable validation evidence. The earlier `legacy_01_outputs/` directory remains local and unversioned.

## Cross-phase controls

- `../FROM_SCRATCH_DATA_ONLY_ROADMAP_2023.md`: current governing roadmap for the from-scratch proteomics, bioinformatics, data-science and AI project.
- `../IMPLEMENTATION_SPEC_DATA_SCIENCE_AI_2023.md`: low-level data contracts, algorithms, resampling design, tests, file outputs and implementation order.
- `../THESIS_EXTENSION_DATA_SCIENCE_AI_PLAN_2023.md`: superseded thesis-relative plan retained as a historical record.
- `../RESEARCH_EVIDENCE_BASE.md`: cited methodological evidence available through 2023.
- `../DATA_ONLY_LIMITATIONS_AND_ASSUMPTIONS.md`: permanent limitations, working assumptions, safeguards, and claim restrictions.
- `../EUROPEAN_ALIGNMENT_AND_TERMINOLOGY.md`: Europe-facing terminology, metadata, governance, and translation context.
- `THESIS_EXTENSION_2023_DATA_ONLY_RESPONSE_RECORD.md`: detailed record of the thesis-comparison, AI/data-science scope, strict data-only boundary, and redesigned Phase 0–6 response.
- `FROM_SCRATCH_ROADMAP_RESPONSE_RECORD.md`: detailed record of the decision to restart from first principles and the resulting proteomics, bioinformatics, data-science and AI design.

## Update convention

Each substantive phase discussion should update its detailed response record or work log. Analytical scripts should regenerate machine-derived reports. New analytical evidence must respect the 2023 cutoff unless the user explicitly changes it; it should be added to `RESEARCH_EVIDENCE_BASE.md` and cited near the relevant decision.

Every future plan, validation report, detailed response record, and work log must name its upstream milestone, downstream hand-off, and inherited limitation IDs. Literature may select a method or define a sensitivity analysis, but it may not manufacture project-specific clinical, pathology, laboratory, batch, or governance facts.
