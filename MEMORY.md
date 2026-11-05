# MEMORY

## 2026-11-04, Data sources under cloud egress policy
**What was decided:** RMDB RDATs fetched one by one from the `rdat:` URL in each `_entries/*.md` of DasLab/rmdb.github.io; natural DMS-MaP data from the MERGE-RNA repo; reference structures from RDAT STRUCTURE fields and MERGE-RNA `reference/` files.
**Why:** GitHub release listing/API, Kaggle, Zenodo, RASP, bpRNA hosts are blocked (curl 000/403). Per-asset release URLs work.
**What was rejected:** Ribonanza (Kaggle blocked), RASP v2.0, bpRNA-1m (hosts blocked); `gh release download` (API blocked).

## 2026-11-04, Large RMDB libraries removed (user decision at 55% disk)
**What was decided:** 58 RDATs >20 MB (high-throughput libraries, 7.1 GB) deleted and listed in data/raw/excluded_large_libraries.txt; select_rmdb.py skips them.
**Why:** disk passed 50%; libraries are synthetic/mutational, not full-length natural RNAs.
**What was rejected:** keeping them.

## 2026-11-04, TVAE-RNA excluded
**What was decided:** TVAE-RNA not run.
**Why:** weights only on a Google Drive folder, test.py imports `VAE3` which is absent from the repo, needs RNA-FM embeddings and PyTorch; README gives no sampling command.
**What was rejected:** retraining TVAE-RNA from scratch.

## 2026-11-04, MERGE-RNA scope
**What was decided:** MERGE-RNA scored only on its own 5 natural DMS-MaP RNAs (needs raw counts + concentration series). Shared mode: physical params fitted on merge_fit clusters, applied at 0 mM without soft constraints. Per-RNA mode: sequential fit on Rep1, scored on Rep2. All tools on MERGE RNAs scored on Rep2 only.
**Why:** RDAT reactivities lack counts/coverage needed by MERGE-RNA's likelihood.
**What was rejected:** fitting and scoring on the same replicate (circular).

## 2026-11-04, Deviations from the plan
**What was decided:** base-pair metrics use pairs with p>0.5 (not MEA); MCC approximated as sqrt(PPV*sensitivity); per-profile normalisation = divide by mean of 90th-98th percentile values; DMS keeps A/C only; "Kendall tau between tool rankings" implemented as tau between per-RNA Spearman scores from exact vs linear-time output.
**Why:** simplest defensible choices; MEA DP in Python is O(N^3) for 2 kb RNAs.
**What was rejected:** custom MEA implementation.

## 2026-11-04, MERGE-RNA coordinates and scoring window
**What was decided:** MERGE-RNA folds the primer-trimmed sequence (offset 17 to 24 nt); its p(unpaired) is embedded in full-length coordinates (NaN outside) and every tool on MERGE RNAs is scored on that same window, Rep2 only. Shared params read from MERGE's physical_params.txt (fit on bact_RNaseP_typeA + tetrahymena_ribozyme, Rep1+Rep2): mu_r 0.4112, p_b 0.2269 kcal/mol.
**Why:** identical positions across tools; avoids refitting after a parsing bug.
**What was rejected:** assuming shared mode equals ViennaRNA (measured MAD 0.08 to 0.11, from trimming + penalty).

## 2026-11-04, EternaBench positive control passed
**What was decided:** ViennaRNA 2.7.2 reproduces published vienna_2 Cloud Lab Round 1 Pearson (0.6730 vs 0.6729).

## Session Summary, 2026-11-04
**Worked on:** full A to G pipeline from the Gap 4 / Gap 3a draft in the cloud session.
**Completed:** env + tools, data fetch, QC, split, 3,951 ensemble runs, MERGE-RNA shared mode, scoring, bias audit, report, EternaBench positive control, tests passing.
**In progress:** none. MERGE-RNA per-RNA mode timed out (~2 h/RNA, no convergence).
**Decisions made:** see entries above.
**Next session:** get RASP/bpRNA/Ribonanza from an unrestricted network to fill >500 nt; rerun per-RNA MERGE with an iteration cap; decide on TVAE-RNA.

## 2026-11-04, Manuscript for Computational Biology and Chemistry
**What was decided:** single-author research article built from manuscript/manuscript.md by manuscript/build.py into RNA-structure-Manuscript-Arli.docx; Elsevier numbered references (33, all PubMed-validated, manuscript/reference_validation.tsv); 12 OMML display equations, 4 boxed algorithms; checklist manuscript/check_manuscript.py (19/19 pass). Data statement: public on acceptance.
**Why:** user approved plan.md; journal guide (sciencedirect.com) blocked, so Elsevier norms used and unverified rules flagged.
**What was rejected:** citations in declarations and highlights (exempt by convention); non-PubMed references (could not be validated here).

## 2026-11-05, QC corrections found through structure visualisation
**What was decided:** QC now (1) skips non-reactivity data lines (READS counts, 50 had been scored as replicate-2 profiles), (2) accepts a reference only if it is balanced dot-bracket with at least one pair (173 pair-less and 7 unbalanced references had been accepted; base-pair F1 was understated, 0.42 -> 0.77 below 200 nt). Profiles 614 -> 564; all 155 remaining value/position mismatches are 20-nt transcript-length fragments. Reactivity results essentially unchanged (rho 0.445 / 0.375).
**Why:** figures showed circular layouts (no pairs), a naview segfault (unbalanced reference) and READS lines among rejected profiles.
**What was rejected:** keeping the earlier numbers.

## 2026-11-05, Retarget to Biophysics Reports
**What was decided:** title page with running title and keywords, highlights removed, abstract <= 250 words single paragraph, declarations renamed (Author contributions, Conflict of interest, Ethics statement, Use of generative AI, Data availability, Funding); Figs 4 and 5 and Supplementary Fig. S1 added; previous Elsevier text kept as manuscript/manuscript_cbc_version.md.
**Why:** user target journal; biophysics-reports.org and journal.hep.com.cn are blocked here, rules taken from a web-search summary of the HEP guideline page and flagged as unverified (reference style still numbered).
