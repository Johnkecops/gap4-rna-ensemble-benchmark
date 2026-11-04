# Plan: research manuscript from the gap4 codebase

Target journal: *Computational Biology and Chemistry* (Elsevier), research article.
Output: `/home/user/gap4/RNA-structure-Manuscript-Arli.docx`
Source of truth: this repository (scripts, `config.json`, `results/*.tsv`, `results/REPORT.md`, figures,
`MEMORY.md`, `ERRORS.md`) plus the 11 references in the Gap 4 recapitulation draft.

## 0. Constraints found while planning (2026-11-04)

| Item in the prompt | Status in this cloud session | What I will do |
|---|---|---|
| Author guide (sciencedirect.com) | Blocked: host outside allowed network domains; web search did not return the guide text | Option A: you allow `www.sciencedirect.com` and I read the guide. Option B: I use the standard Elsevier research-article structure from memory and list every format rule I could not verify in a cover note |
| /scientific-writing, /plan, /crossref, /goal | Not installed in this session | Plan written by hand (this file); IMRAD drafting done directly; DOI check done against PubMed records; "goal" replaced by a scripted checklist loop (section 6) |
| /humanizer, /remove-ai-marks, /citation-management, /pubmed-database | Installed as `anthropic-skills:humanizer`, `anthropic-skills:remove-ai-watermarks`, `anthropic-skills:citation-management`, `anthropic-skills:pubmed-database`; PubMed MCP works | Used as asked |
| Crossref / doi.org / OpenAlex | Blocked (HTTP 000 / 403) | A DOI counts as valid when PubMed returns the same DOI for the PMID. A reference with no PubMed record cannot be checked here; I will avoid such references, or flag any I keep |
| Microsoft Equation | No Word available | Write equations in LaTeX and convert with pandoc to native Office Math (OMML), which Word opens as editable Equation objects with real Unicode math characters |

## 1. Inputs to read (reverse engineering)
- Methods: `scripts/00_fetch_data.sh` to `06_report.py`, `foldtools.py`, `run_tool.py`, `merge_fit.py`,
  `positive_control_eternabench.py`, `tests/test_pipeline.py`, `config.json`, `environments/*`.
- Results: `leaderboard.tsv`, `scaling.tsv`, `bias_summary.tsv`, `bias_per_rna.tsv`, `reliable_length.tsv`,
  `qc_dropped.tsv`, `positive_control_eternabench.txt`, `test_pipeline_deviations.txt`, `fig_bias.png`, `fig_scaling.png`.
- Decisions and limits: `MEMORY.md`, `ERRORS.md`, the limitations section of `REPORT.md`.

## 2. Manuscript structure (Elsevier research article, numbered sections)
1. Title page: title, author (Arli Aditya Parikesit, i3L University, ORCID 0000-0001-8716-3926), corresponding e-mail
2. Highlights (3 to 5 bullets, each 85 characters or fewer; Elsevier norm, to be confirmed against the guide)
3. Abstract (single paragraph, no citations, 250 words or fewer unless the guide says otherwise) and keywords (up to 6)
4. 1 Introduction: ensembles and chemical probing; existing benchmarks (EternaBench, Ribonanza); the gap; aims (Gap 4, Gap 3a)
5. 2 Materials and methods: 2.1 data sources and QC; 2.2 deduplication and split; 2.3 ensemble generators;
   2.4 MERGE-RNA modes; 2.5 metrics; 2.6 bias audit; 2.7 positive control and consistency test;
   2.8 software and reproducibility
6. 3 Results: 3.1 dataset; 3.2 consistency and positive control; 3.3 reactivity leaderboard; 3.4 runtime and memory;
   3.5 linear-time bias vs length and beam; 3.6 MERGE-RNA
7. 4 Discussion: interpretation, comparison with prior work, limitations (thin coverage above 418 nt, blocked sources,
   unconverged per-RNA MERGE fits, TVAE-RNA not run, temperature mismatch), next steps
8. 5 Conclusions (no citations)
9. CRediT statement, declaration of competing interest, declaration of generative AI use (required by Elsevier),
   data and code availability (the repository is private now; the statement will say "available on request / to be
   made public on acceptance" unless you decide otherwise), funding, acknowledgements
10. References
11. Tables (data set summary, leaderboard, bias summary and reliable length range, software versions) and figures
    (pipeline flowchart redrawn as a figure, `fig_bias.png`, `fig_scaling.png`), each with a caption

## 3. Equations and pseudo-code
- Equations (numbered, OMML, Unicode): Boltzmann probability of a structure; partition function; base-pair and
  unpaired probability; per-profile normalisation; Spearman ρ; ROC-AUC; bootstrap interval; MAD of p(unpaired) and of
  pair probability; Jensen–Shannon divergence; Kendall τ; F1 and MCC (√(PPV·sensitivity)); MERGE-RNA pairing penalty;
  O(N³)/O(N²) cost. Every sentence that explains an equation carries a citation (McCaskill, Ding and Lawrence,
  Mathews, Lin for JSD, the MERGE-RNA paper, statistics sources found in PubMed).
- Pseudo-code in boxed, single-cell, bordered, monospaced tables: Algorithm 1 profile QC and normalisation;
  Algorithm 2 cluster-level split; Algorithm 3 ensemble generation and scoring; Algorithm 4 bias audit and reliable length.

## 4. References
- Start from the 11 references in the recapitulation draft. Add only what the text needs and PubMed can confirm,
  for example Turner 2004 (NNDB), Nextflow (Di Tommaso 2017), CD-HIT (Li and Godzik 2006; Fu 2012),
  LinearPartition (Zhang 2020), RNA Framework, bootstrap and rank statistics in biomedical literature.
- Each reference: fetch PMID metadata with the PubMed tool; confirm title, authors, year, journal and DOI match.
  Record results in `manuscript/reference_validation.tsv` (ref, PMID, DOI, status).
- Style: Elsevier style as required by the guide (numbered or name-year, decided after the guide is read or confirmed by you).
- Rules from the prompt: every paragraph outside the abstract and conclusions has at least one citation; every
  reference in the list is cited in the body; no citations in the abstract or conclusions.

## 5. Build and post-processing
1. Draft the manuscript in Markdown with LaTeX math in `manuscript/` (source kept in Git).
2. Convert with pandoc to `.docx` (OMML equations), then python-docx to apply the format: Times New Roman 12 pt,
   double spacing, A4, 2.5 cm margins, continuous line numbers, page numbers, pseudo-code boxes, table and figure captions.
3. Humanizer pass on body text only (skip references, equations, pseudo-code, tables).
4. Remove-ai-watermarks pass on the final `.docx` (hidden characters, metadata, document properties).
5. Save as `/home/user/gap4/RNA-structure-Manuscript-Arli.docx`; also render a PDF preview for my own visual check.

## 6. Verification loop (stand-in for /goal), repeated until all pass
`manuscript/check_manuscript.py` will assert:
- [ ] all required sections present, in order; highlights count and length; abstract word count; keyword count
- [ ] every body paragraph (not abstract or conclusions) contains at least one citation
- [ ] abstract and conclusions contain no citation
- [ ] every reference is cited, and every citation resolves to a reference
- [ ] every reference has PubMed-confirmed PMID and DOI (or is flagged)
- [ ] equations are OMML objects with Unicode math; each pseudo-code block sits in a bordered box
- [ ] numbers in the text match `results/*.tsv` (scripted cross-check of key values)
- [ ] no zero-width or other hidden watermark characters; clean document properties
- [ ] file saved under the required name

## 7. Honesty constraints carried into the text
- Report limits plainly: 2 distinct RNAs above 418 nt; Ribonanza, RASP and bpRNA-1m unreachable; per-RNA MERGE-RNA
  mode did not converge; TVAE-RNA not run; temperatures differ (37 °C model vs often 24 °C experiments).
- No claim beyond the data; results framed as a pilot benchmark plus a bias audit.

## 8. Estimated effort
About 2 to 3 hours of session time (drafting, about 25 to 35 PubMed lookups, build, humanizer pass, checks).
Disk impact under 50 MB.
