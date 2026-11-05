# Gap 4 / Gap 3a: probing-based ensemble benchmark and linear-time sampler bias audit

Pipeline of seven stages (A to G) from the Gap 4 recapitulation draft (November 2026).
Results summary: `results/REPORT.md`. Decisions and deviations: `MEMORY.md`. Failed approaches: `ERRORS.md`.

## Layout
```
config.json           all thresholds, fixed before scoring
scripts/00_fetch_data.sh  A  RMDB RDATs, EternaBench clone, MANIFEST.tsv (sha256 + date)
scripts/select_rmdb.py       RMDB entry selection (no Eterna, no >20 MB libraries)
scripts/01_qc.py          B  RDAT + MERGE-RNA .rc -> data/processed/{reactivity.tsv.gz,rnas.tsv}
scripts/02_split.sh       C  cd-hit-est 0.8 clusters, MERGE-RNA fit / held-out split, length strata
scripts/03_ensembles.nf   D  Nextflow: ViennaRNA exact, LinearPartition-V, LinearSampling (beams 100/200/500/0)
scripts/run_tool.py          D worker (one RNA x tool x beam -> results/ensembles/<rna>/<tool>_b<beam>.npz)
scripts/merge_fit.py      D  MERGE-RNA shared-parameter and per-RNA (Rep1 fit, Rep2 score) modes
scripts/04_score.py       E  Spearman, ROC-AUC (bootstrap CI), base-pair F1/MCC, runtime/RAM scaling
scripts/05_bias_audit.py  F  MAD pair/unpaired probability, JSD (<=150 nt), d-Spearman, Kendall tau
scripts/06_report.py      G  leaderboard, error bound, reliable length range -> results/REPORT.md
scripts/07_structure_viz.py G structure figures (QC diagnostic, reactivity on structure, p_unp vs reactivity)
scripts/diskcheck.sh         disk guard (exit 1 at >=50 %, 2 at >=95 % of session allowance)
tests/test_pipeline.py       infinite-beam LinearPartition/LinearSampling must reproduce ViennaRNA
```

## Run
```bash
export MAMBA_ROOT_PREFIX=/opt/mamba PATH=/opt/mamba/envs/gap4/bin:$PATH
python -m pytest -q tests
bash scripts/00_fetch_data.sh && python scripts/01_qc.py && bash scripts/02_split.sh
nextflow run scripts/03_ensembles.nf -with-trace results/nf_trace.tsv
python scripts/merge_fit.py shared && python scripts/merge_fit.py perRNA
python scripts/04_score.py && python scripts/05_bias_audit.py && python scripts/06_report.py && python scripts/07_structure_viz.py
```

## Software (installed in this cloud session, 2026-11-04)
micromamba 2.9.0 (conda-forge; Miniforge installer is on GitHub releases, which were blocked), env `gap4`
(environments/environment.yml): Python 3.11, ViennaRNA 2.7.2, cd-hit, Nextflow 26.10.4, numpy/scipy/pandas/
scikit-learn/matplotlib, GNU time, rdat_kit 1.8.1 (from source). Built from source with g++ -O3 at the commits in
`environments/tool_commits.tsv`: LinearPartition, LinearSampling. MERGE-RNA run from its clone. TVAE-RNA cloned
but not run (see MEMORY.md).

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality
