::: titlepage
# Benchmarking RNA secondary structure ensemble generators against chemical mapping of natural RNAs, with an error audit of linear-time sampling

Arli Aditya Parikesit^a,\*^

^a^ Department of Biotechnology, School of Life Sciences, Indonesia International Institute for Life Sciences (i3L), Jl. Pulomas Barat Kav. 88, Jakarta Timur 13210, Indonesia

\* Corresponding author. E-mail: arli.parikesit@i3l.ac.id. ORCID: 0000-0001-8716-3926
:::

## Highlights

- Open pipeline scores RNA ensemble generators on 439 RMDB and MERGE-RNA RNAs
- Exact and linear-time ViennaRNA-model tools give the same reactivity agreement
- Beams of 100 and 200 shift unpaired probabilities by up to 0.16 above 1,000 nt
- Beam 500 and the exact sampler stay within 0.014 of exact values up to 2,216 nt
- The pipeline reproduces the published EternaBench ViennaRNA score (0.673)

## Abstract

Chemical mapping reports reactivity averaged over every structure an RNA adopts, so it can test predicted secondary structure ensembles directly. Existing probing-based benchmarks rely on short synthetic constructs or on windows of longer RNAs, and none scores ensemble generators on full-length natural RNAs. We built a reproducible seven-stage pipeline that collects natural-RNA probing profiles from the RNA Mapping Database and from the MERGE-RNA DMS-MaP data, filters and normalises them, clusters sequences at 80% identity, and runs exact ViennaRNA sampling, LinearPartition and LinearSampling at four beam widths, and MERGE-RNA. The final set held 439 RNAs of 40 to 2,216 nucleotides and 614 reactivity profiles. On held-out RNAs shorter than 200 nucleotides, the Spearman correlation between predicted unpaired probability and reactivity was about 0.44 for every ViennaRNA-model tool, and about 0.38 between 200 and 1,000 nucleotides, with overlapping bootstrap intervals. A bias audit against the exact partition function showed that beams of 100 and 200 changed unpaired probabilities by a mean of 0.09 to 0.12 on RNAs of 1,266 and 2,216 nucleotides, whereas a beam of 500 or an infinite beam kept the mean absolute difference below 0.012. Linear-time error did not change which RNAs were predicted well (Kendall τ 0.97 to 0.99). Only two distinct RNAs exceeded 418 nucleotides, so the long-RNA error bound is preliminary. The pipeline reproduced the published EternaBench ViennaRNA score (Pearson 0.6730 against 0.6729), and its code and data are versioned for reuse.

**Keywords:** RNA secondary structure; partition function; stochastic sampling; chemical probing; benchmark; LinearSampling

## 1. Introduction

Many RNAs do not fold into a single structure. At equilibrium they populate several secondary structures, and functional states such as riboswitch conformations can be minor members of that population [@sacco2026; @mccaskill1990]. The partition function introduced by McCaskill gives the probability of every base pair over the full Boltzmann ensemble, and stochastic traceback draws whole structures from the same distribution [@mccaskill1990; @ding2003]. Pair probabilities from the partition function also indicate which predicted pairs are reliable, which makes the ensemble view useful even when a single model is needed [@mathews2004rna; @lu2009].

Chemical mapping measures how accessible each nucleotide is to a reagent such as dimethyl sulfate (DMS) or a SHAPE acylating agent, and the signal is averaged over every molecule in the sample [@siegfried2014; @zubradt2017]. A reactivity profile therefore reflects the whole ensemble and can test predicted ensembles directly [@waymentsteele2022]. Probing data have long been used as soft constraints that improve folding accuracy, from the first chemical modification constraints to SHAPE pseudo-energies [@mathews2004pnas; @deigan2009; @lowweeks2010].

EternaBench compared folding packages with more than 20,000 synthetic constructs and showed that probing data discriminate between packages in ensemble-oriented tasks [@waymentsteele2022]. Its authors also noted that natural-RNA benchmarks tend to be dominated by well-studied tRNAs and rRNAs, and that the synthetic constructs are shorter than many natural RNAs [@waymentsteele2022]. Ribonanza added about two million chemical mapping profiles, but oligonucleotide synthesis limits its constructs to roughly 90 to 130 nucleotides, so longer RNAs enter only as windows [@he2024]. The RNA Mapping Database (RMDB) holds full-length profiles of natural RNAs in a standard format, yet we found no shared benchmark that scores ensemble generators on such RNAs [@cordero2012; @waymentsteele2022].

Long RNAs raise a second problem. The exact partition function costs O(N³) time and O(N²) memory for a sequence of length N, so genome-scale RNAs need approximations [@mccaskill1990; @mathews2004rna]. LinearFold, LinearPartition and LinearSampling replace the bottom-up recursion with a left-to-right scan and beam pruning, which makes the cost linear in N for a fixed beam [@huang2019; @zhang2020; @zhang2023]. Beam pruning discards low-scoring states, so the sampled distribution is an approximation, and how far it departs from the exact ensemble on natural RNAs of different lengths has not been measured against probing data [@zhang2020; @zhang2023]. Newer generators add further variety: MERGE-RNA fits a physics-based model of the probing experiment to DMS data, and TVAE-RNA learns a latent space of structures with a transformer variational autoencoder [@sacco2026; @mei2025].

This study has two aims. The first is an open benchmark that scores ensemble generators by agreement between predicted unpaired probability and measured reactivity on full-length natural RNAs, with base-pair agreement and runtime as secondary measures [@waymentsteele2022; @cordero2012]. The second is an audit of linear-time sampling that compares LinearPartition and LinearSampling with exact ViennaRNA output on the same RNAs and reports the length range over which the approximation stays within a preset tolerance [@lorenz2011; @zhang2020; @zhang2023].

## 2. Materials and methods

### 2.1. Pipeline overview

The workflow has seven stages (A to G), each implemented as one script in a Git repository (Fig. 1). Stages A to D are shared: data acquisition, quality control, deduplication with a train and test split, and ensemble generation. Stage E scores the benchmark and stage F runs the bias audit on the same ensemble files, so no RNA is folded twice [@ditommaso2017; @wilkinson2016]. Stage D runs as a Nextflow workflow with one process per RNA, tool and beam, and every threshold was fixed in a configuration file before any scoring [@ditommaso2017].

![Fig. 1. Workflow of the benchmark (stages A to G). Stages A to D are shared; stage E scores agreement with reactivity and stage F audits linear-time sampling against exact ViennaRNA output on the same ensemble files [@ditommaso2017; @zhang2023].](fig_workflow.png)

### 2.2. Data sources

RMDB entries were read from the metadata of the database's site repository, and each RDAT file was downloaded from the release link given in its entry [@cordero2012]. Entries from the Eterna category and high-throughput design libraries were excluded because they are synthetic, and 58 library files larger than 20 MB were removed for the same reason; 767 RDAT files remained [@cordero2012; @waymentsteele2022]. Natural DMS-MaP data for five structured RNAs (bacterial RNase P type A, the hc16 ligase, the HCV IRES, the *Tetrahymena* ribozyme and the *Vibrio cholerae* glycine riboswitch) came from the MERGE-RNA repository as RNA Framework count files with two replicates and six DMS concentrations each [@sacco2026; @incarnato2018]. The EternaBench repository served as a positive control, and Ribonanza, RASP and bpRNA-1m could not be used because their hosts were unreachable from the computing environment [@waymentsteele2022; @he2024].

### 2.3. Quality control and normalisation

RDAT files were parsed with RDATKit, and only wild-type profiles without ligand annotations were kept, up to two replicates per construct [@cordero2012]. A profile was removed when the source marked it as weak, when the mean reactivity divided by the mean error fell below 1.0, when fewer than 20 positions were measured, or when sequence positions and values could not be aligned unambiguously [@cordero2012; @lowweeks2010]. Five measured positions at each end were masked, and DMS profiles kept only adenosine and cytidine, the bases that DMS modifies at the Watson-Crick face [@zubradt2017]. For the MERGE-RNA data, reactivity was the mutation rate at the highest DMS concentration minus the rate at 0 mM, using positions with at least 1,000 reads in both samples [@sacco2026; @incarnato2018].

Each profile was normalised by the mean of its values between the 90th and 98th percentiles, a variant of the 2%/8% rule used for SHAPE data [@deigan2009; @lowweeks2010]. For profile values $r_1,\dots,r_n$ with percentiles $q_{90}$ and $q_{98}$, the normalised reactivity is

$$\tilde{r}_k = \frac{r_k}{\frac{1}{|T|}\sum_{l \in T} r_l}, \qquad T = \{\, l : q_{90} \le r_l \le q_{98} \,\}. \qquad (1)$$

Dividing by a high but not extreme quantile keeps a few outliers from setting the scale, which is the purpose of the original rule [@lowweeks2010; @deigan2009]. Algorithm 1 summarises the stage.

```{.algorithm}
Algorithm 1. Profile quality control and normalisation (stage B)
Input: RDAT files F, MERGE-RNA counts M; SNR_min = 1.0, e = 5, cov_min = 1000
Output: table of (rna_id, profile, position, nucleotide, probe, reactivity)
 1: for each RDAT file f in F do
 2:     for each construct c in f do
 3:         W <- wild-type profiles of c without ligand, at most 2
 4:         for each profile d in W do
 5:             if flag(d) = weak or SNR(d) < SNR_min then reject d
 6:             if |positions(d)| != |values(d)| then reject d            # ambiguous alignment
 7:             mask first e and last e measured positions
 8:             if probe(d) = DMS then keep only A and C
 9:             s <- mean{ r : q90(d) <= r <= q98(d) };  d <- d / s       # Eq. (1)
10:             append d to table
11: for each MERGE-RNA system m in M and replicate k in {1, 2} do
12:     d <- rate(c_max) - rate(0 mM) at positions with coverage >= cov_min
13:     apply steps 7 to 10 to d
14: return table
```

### 2.4. Deduplication, split and length strata

All sequences were clustered with cd-hit-est at 80% identity (word size 5) [@li2006; @fu2012]. Sets were assigned cluster by cluster: half of the clusters that contained a MERGE-RNA system went to a MERGE-RNA fitting set and the rest to a held-out set, so no fitted RNA or close homologue appears among scored RNAs [@li2006; @sacco2026]. All other clusters were held out. Each RNA was tagged with a length stratum (under 200, 200 to 1,000, 1,000 to 3,000 and over 3,000 nucleotides) for stratified reporting [@waymentsteele2022]. Algorithm 2 gives the procedure.

```{.algorithm}
Algorithm 2. Cluster-level split (stage C)
Input: sequences S, identity threshold t = 0.8, fit fraction f = 0.5, random seed
Output: set label and length stratum for every RNA
 1: C <- cd-hit-est(S, identity = t, word = 5)
 2: K <- clusters in C that contain a MERGE-RNA system;  shuffle K with seed
 3: K_fit <- first round(f * |K|) clusters of K
 4: for each RNA x in S do
 5:     set(x) <- merge_fit if cluster(x) in K_fit else held_out
 6:     stratum(x) <- bin(length(x); 200, 1000, 3000)
 7: return set, stratum
```

### 2.5. Ensemble generation

For an RNA sequence, the Boltzmann probability of a secondary structure $s$ with free energy $\Delta G(s)$ is given by the partition function $Z$ [@mccaskill1990; @mathews2004rna]:

$$P(s) = \frac{e^{-\Delta G(s)/RT}}{Z}, \qquad Z = \sum_{s' \in \Omega} e^{-\Delta G(s')/RT}, \qquad (2)$$

where Ω is the set of all pseudoknot-free structures, R the gas constant and T the absolute temperature [@mccaskill1990; @turner2010]. The probability that nucleotides i and j pair, and the probability that nucleotide i is unpaired, follow from the same ensemble [@mccaskill1990; @mathews2004rna]:

$$p_{ij} = \sum_{s \in \Omega:\,(i,j) \in s} P(s), \qquad p^{\mathrm{u}}_i = 1 - \sum_{j \ne i} p_{ij}. \qquad (3)$$

Stochastic traceback samples structures with probability P(s) once Z and its partial sums are known, and the frequency of a feature among K samples estimates its ensemble probability [@ding2003]. For a sample set $\{s^{(1)},\dots,s^{(K)}\}$ the sampled unpaired probability is

$$\hat{p}^{\mathrm{u}}_i = \frac{1}{K}\sum_{k=1}^{K} \mathbb{1}\!\left[\, i \text{ unpaired in } s^{(k)} \right]. \qquad (4)$$

Exact computation of Z by dynamic programming takes O(N³) time and O(N²) memory, which is the cost the linear-time methods avoid [@mccaskill1990; @mathews2004rna; @zhang2020].

Four generators were run on every RNA. ViennaRNA 2.7.2 computed the exact partition function, pair probabilities and 1,000 stochastic samples with default Turner 2004 parameters at 37 °C, dangles 2, and unique multiloop decomposition switched on as required for sampling [@lorenz2011; @turner2010; @ding2003]. LinearPartition in its ViennaRNA-parameter mode (-V) computed approximate pair probabilities, and LinearSampling drew 1,000 samples, each at beam widths of 100, 200, 500 and infinity (beam 0) [@zhang2020; @zhang2023]. Exact runs (ViennaRNA and the infinite beams) were limited to RNAs of at most 4,000 nucleotides, a limit no RNA in the final set reached [@mathews2004rna; @zhang2020]. Pairs with probability below 10⁻⁴ were not stored, runtime was the wall-clock time of each run, and peak memory was measured with GNU time for the compiled tools and as the resident-memory increase for in-process ViennaRNA [@lorenz2011; @zhang2020]. Before the first full run, a consistency test required both linear tools at infinite beam to reproduce ViennaRNA unpaired probabilities on five random 100-nucleotide sequences, which rules out a mismatch of energy settings [@zhang2020; @zhang2023; @lorenz2011].

### 2.6. MERGE-RNA

MERGE-RNA models the probing experiment explicitly: a reagent concentration c sets a chemical potential, and paired bases receive a free-energy penalty that depends on it [@sacco2026]. With fitted parameters $\mu_r$ and $p_b$, $\beta = 1/k_BT$ and $c$ in mM, the penalty applied to each paired base is

$$\mu_c = \mu_r + k_B T \ln\frac{c + 0.1}{1000}, \qquad m(c) = -k_B T \ln\frac{1 + e^{\beta(\mu_c - p_b)}}{1 + e^{\beta \mu_c}}. \qquad (5)$$

Pair probabilities are then computed under this penalty and, optionally, under per-nucleotide soft constraints λ fitted by a maximum-entropy principle [@sacco2026]. We ran MERGE-RNA in two modes. In the shared-parameter mode, the physical parameters were fitted on both replicates of the systems in the fitting set and applied at 0 mM without soft constraints to the held-out systems [@sacco2026]. In the per-RNA mode, physical parameters and soft constraints were fitted on replicate 1 of each held-out system for scoring on replicate 2, because fitting and scoring the same profile would be circular [@sacco2026; @incarnato2018]. MERGE-RNA folds the primer-trimmed sequence, so its probabilities were placed in full-length coordinates and every tool on these RNAs was scored on the same trimmed window and on replicate 2 only [@sacco2026].

### 2.7. Agreement with reactivity and with reference structures

Scoring used held-out RNAs only. For each profile, the Spearman rank correlation between predicted unpaired probability and normalised reactivity at the measured positions was computed [@schober2018; @waymentsteele2022]:

$$\rho = 1 - \frac{6\sum_{k=1}^{n} d_k^2}{n(n^2 - 1)}, \qquad (6)$$

where $d_k$ is the rank difference at position $k$ and $n$ the number of positions; ties receive mid-ranks [@schober2018]. Nucleotides with reactivity above the 75th percentile of their profile were labelled reactive, and the area under the ROC curve (AUC) measured how well $p^{\mathrm{u}}$ separated reactive from unreactive positions [@hanley1982]. The AUC equals the probability that a randomly chosen reactive position receives a higher score than a randomly chosen unreactive one [@hanley1982]:

$$\mathrm{AUC} = \Pr\!\left( p^{\mathrm{u}}_a > p^{\mathrm{u}}_b \;\middle|\; a \in \text{reactive},\; b \in \text{unreactive} \right). \qquad (7)$$

Scores were averaged over profiles of the same RNA and then over RNAs within each stratum [@waymentsteele2022]. For a stratum with per-RNA scores $x_1,\dots,x_m$, 95% intervals came from $B = 1{,}000$ bootstrap resamples of RNAs [@henderson2005]:

$$\bar{x}^{*}_{b} = \frac{1}{m}\sum_{k=1}^{m} x_{j_{bk}}, \quad j_{bk} \sim \mathcal{U}\{1,\dots,m\}, \qquad \mathrm{CI}_{95} = \left[\, Q_{0.025}(\bar{x}^{*}),\; Q_{0.975}(\bar{x}^{*}) \,\right]. \qquad (8)$$

Resampling whole RNAs rather than positions respects the dependence between neighbouring nucleotides of one molecule [@henderson2005; @waymentsteele2022].

Where the RDAT file or the MERGE-RNA repository supplied a reference dot-bracket structure of matching length, the predicted structure was the set of pairs with $p_{ij} > 0.5$, which is always free of conflicts [@mathews2004rna; @lu2009]. With sensitivity Sen = TP/(TP + FN) and positive predictive value PPV = TP/(TP + FP) over base pairs,

$$F_1 = \frac{2\,\mathrm{PPV}\cdot\mathrm{Sen}}{\mathrm{PPV} + \mathrm{Sen}}, \qquad \mathrm{MCC} \approx \sqrt{\mathrm{PPV}\cdot\mathrm{Sen}}, \qquad (9)$$

where the second expression approximates the Matthews correlation coefficient when true negatives far outnumber the other classes, as they do for base pairs [@matthews1975; @lu2009]. Base-pair scores were treated as secondary to reactivity scores because comparative and annotated structures carry their own errors and exclude pseudoknots [@rivas1999; @lyngso2000].

### 2.8. Bias audit of linear-time sampling

For every RNA with an exact ViennaRNA run, each LinearPartition and LinearSampling output was compared with the exact output [@zhang2020; @zhang2023; @lorenz2011]. The mean absolute difference (MAD) in unpaired probability and in pair probability were

$$\mathrm{MAD}_{\mathrm{u}} = \frac{1}{N}\sum_{i=1}^{N}\left| p^{\mathrm{u,exact}}_i - p^{\mathrm{u,lin}}_i \right|, \qquad \mathrm{MAD}_{\mathrm{p}} = \frac{1}{|U|}\sum_{(i,j) \in U}\left| p^{\mathrm{exact}}_{ij} - p^{\mathrm{lin}}_{ij} \right|, \qquad (10)$$

with $U$ the union of pairs stored by either method [@zhang2020; @mccaskill1990]. For LinearSampling, $p^{\mathrm{lin}}$ denotes sample frequencies from Eq. (4) [@ding2003; @zhang2023]. For RNAs of at most 150 nucleotides, the Jensen-Shannon divergence between the frequency distributions $P$ and $Q$ of distinct sampled structures was also computed, with $M = (P + Q)/2$ [@capra2007]:

$$\mathrm{JSD}(P \,\Vert\, Q) = \tfrac{1}{2} D_{\mathrm{KL}}(P \,\Vert\, M) + \tfrac{1}{2} D_{\mathrm{KL}}(Q \,\Vert\, M), \qquad D_{\mathrm{KL}}(P \,\Vert\, M) = \sum_{s} P(s)\log_2\frac{P(s)}{M(s)}. \qquad (11)$$

At longer lengths almost every sample is unique, so the structure-level estimate carries no information and was not computed [@ding2003; @capra2007].

Two further measures asked whether approximation error changes conclusions. The change in reactivity correlation was $\Delta\rho = \rho_{\mathrm{lin}} - \rho_{\mathrm{exact}}$ per RNA [@schober2018]. Kendall's τ between per-RNA Spearman scores from exact and from linear-time output measured whether the ranking of RNAs by predictability was preserved [@kraemer2006; @rosner2007]:

$$\tau = \frac{n_c - n_d}{\tfrac{1}{2}\,n(n - 1)}, \qquad (12)$$

with $n_c$ concordant and $n_d$ discordant pairs of RNAs [@kraemer2006]. A linear-time setting was called reliable up to length L when every audited RNA no longer than L had $\mathrm{MAD}_{\mathrm{u}}$ at most 0.05, a tolerance fixed before the audit [@zhang2020; @zhang2023]. Runtime and memory scaling were summarised by the slope α of a least-squares fit of log t = α log N + c over RNAs of at least 50 nucleotides [@zhang2020; @huang2019]. Algorithm 3 outlines stages D to F and Algorithm 4 the reliability rule.

```{.algorithm}
Algorithm 3. Ensemble generation, scoring and audit (stages D to F)
Input: RNAs R, tools T = {ViennaRNA, LinearPartition-V, LinearSampling},
       beams B = {100, 200, 500, 0}, K = 1000 samples, L_max = 4000 nt
Output: per-RNA scores, runtime and memory, audit table
 1: for each RNA x in R and each (tool, beam) in T x B do   # Nextflow
 2:     if tool is exact or beam = 0, and length(x) > L_max then skip
 3:     run tool; store p_ij > 1e-4, p_u (Eq. 3 or 4), samples, time, memory
 4: for each held-out RNA x and each stored output o of x do
 5:     for each profile d of x do  rho(d) <- Eq. (6);  AUC(d) <- Eq. (7)
 6:     score(x, o) <- mean over profiles;  if reference exists: F1, MCC <- Eq. (9)
 7: for each stratum and tool: mean score and bootstrap interval (Eq. 8)
 8: for each RNA x with an exact run and each linear output o of x do
 9:     $\mathrm{MAD}_{\mathrm{u}}$, MAD_p <- Eq. (10);  if length(x) <= 150: JSD <- Eq. (11)
10:     delta_rho <- rho(o) - rho(exact)
11: tau(o) <- Kendall tau of per-RNA rho, exact vs o (Eq. 12)
```

```{.algorithm}
Algorithm 4. Reliable length range for a linear-time setting (stage G)
Input: audit rows (length, $\mathrm{MAD}_{\mathrm{u}}$) for one tool and beam, tolerance eps = 0.05
Output: reliable_up_to, first_failure
 1: sort rows by length
 2: F <- rows with $\mathrm{MAD}_{\mathrm{u}}$ > eps
 3: if F is empty then return (max length, none)
 4: first_failure <- min length in F
 5: reliable_up_to <- max length among rows with length < first_failure
 6: return (reliable_up_to, first_failure)
```

### 2.9. Positive control, software and reproducibility

As a positive control, ViennaRNA unpaired probabilities were computed for the EternaBench Cloud Lab Round 1 demonstration set and scored as in the original study: reactivities of all constructs were concatenated, values of zero or below and above the 95th percentile were removed, poly(A) runs of six or more were excluded, and the Pearson correlation was compared with the published value [@waymentsteele2022]. Software ran in a conda environment with Python 3.11, ViennaRNA 2.7.2, CD-HIT, Nextflow 26.10.4, NumPy, SciPy, pandas and scikit-learn; LinearPartition and LinearSampling were compiled from source with g++ -O3, and the commits of all tools are recorded in the repository [@lorenz2011; @li2006; @ditommaso2017; @zhang2020]. Computations used a cloud machine with four CPU cores and 15 GB of memory, with four jobs running in parallel [@ditommaso2017]. TVAE-RNA was installed but not run, because its released weights, a module imported by its test script and the RNA-FM embeddings it needs were not obtainable in this environment [@mei2025].

## 3. Results

### 3.1. Data set

Quality control kept 614 profiles from 439 RNAs (75,568 scored positions) and removed 277 profiles, mostly because positions and values could not be aligned (218 profiles) [@cordero2012; @lowweeks2010]. Most profiles came from SHAPE reagents (1M7, 339 profiles) and DMS (151), with smaller numbers from other probes (Table 1) [@cordero2012; @siegfried2014]. Sequence lengths ranged from 40 to 2,216 nucleotides with a median of 133, and 384 RNAs carried a usable reference structure [@cordero2012; @sacco2026]. Clustering produced 123 clusters; two MERGE-RNA systems (RNase P and the *Tetrahymena* ribozyme) and their cluster members formed the fitting set of 13 RNAs, leaving 426 held-out RNAs [@li2006; @sacco2026]. The held-out set contained 325 RNAs under 200 nucleotides, 98 between 200 and 1,000 and three between 1,000 and 3,000; no RNA exceeded 3,000 nucleotides [@cordero2012; @he2024].

Table 1. Composition of the benchmark after quality control. Probe counts are profiles; length strata count RNAs. Data from RMDB [@cordero2012] and MERGE-RNA [@sacco2026].

| Item | Value |
|---|---|
| RNAs / profiles / scored positions | 439 / 614 / 75,568 |
| Sources | RMDB 434 RNAs; MERGE-RNA DMS-MaP 5 RNAs |
| Profiles by probe | 1M7 339; DMS 151; SHAPE (unspecified) 33; BzCN 32; CMCT 17; 2A3 17; hydroxyl radical 10; other 15 |
| Length (min / median / max) | 40 / 133 / 2,216 nt |
| RNAs with reference structure | 384 |
| Profiles removed in QC | 277 (alignment 218; normalisation 19; weak 18; too few positions 17; SNR < 1.0 5) |
| Clusters (cd-hit-est, 80%) | 123 |
| MERGE-RNA fitting set | 13 RNAs (2 MERGE-RNA systems and homologues) |
| Held-out RNAs by stratum | < 200 nt: 325; 200 to 1,000 nt: 98; 1,000 to 3,000 nt: 3; > 3,000 nt: 0 |

### 3.2. Consistency test and positive control

At infinite beam, LinearPartition reproduced ViennaRNA unpaired probabilities within a maximum absolute difference of 0.0066 on five random 100-nucleotide sequences, and LinearSampling with 2,000 samples stayed within 0.033, inside five binomial standard errors [@zhang2020; @zhang2023; @lorenz2011]. The two tools therefore used the same energy model as ViennaRNA, so differences in the audit reflect beam pruning and sampling noise [@zhang2020; @turner2010]. On the EternaBench Cloud Lab Round 1 set (1,088 constructs, 79,439 nucleotides), ViennaRNA 2.7.2 gave a Pearson correlation of 0.6730 between unpaired probability and reactivity against a published value of 0.6729, which indicates that the scoring steps match the original implementation [@waymentsteele2022; @lorenz2011].

### 3.3. Agreement with reactivity

All ViennaRNA-model tools agreed with reactivity to the same degree (Table 2) [@lorenz2011; @zhang2020; @zhang2023]. Below 200 nucleotides the mean Spearman correlation was 0.444 for exact ViennaRNA and for the linear-time tools at every beam, with bootstrap intervals of about 0.42 to 0.47; the AUC was 0.74 [@schober2018; @hanley1982; @henderson2005]. Between 200 and 1,000 nucleotides, correlation fell to 0.376 for ViennaRNA and to 0.375 to 0.376 for the linear-time tools, with AUC again 0.74 [@schober2018; @hanley1982]. In the 1,000 to 3,000 stratum only three RNA entries representing two distinct RNAs were scored, the correlation was 0.11 to 0.13 and the intervals included zero, so this stratum cannot rank tools [@henderson2005; @waymentsteele2022]. Base-pair F1 of the p > 0.5 structures was 0.42 below 200 nucleotides and 0.35 between 200 and 1,000 for all tools [@lu2009; @matthews1975].

Table 2. Agreement between predicted unpaired probability and reactivity on held-out RNAs (mean over RNAs, 95% bootstrap interval in brackets). Metrics from Eqs. (6) to (9) [@schober2018; @hanley1982; @henderson2005; @lu2009].

| Stratum | Tool | RNAs | Spearman ρ | AUC | bp F1 |
|---|---|---|---|---|---|
| < 200 nt | ViennaRNA (exact) | 299 | 0.444 [0.421, 0.467] | 0.743 | 0.417 |
| < 200 nt | LinearPartition b100 | 299 | 0.444 [0.422, 0.469] | 0.742 | 0.417 |
| < 200 nt | LinearSampling b100 | 299 | 0.444 [0.419, 0.468] | 0.742 | 0.418 |
| < 200 nt | LinearSampling b500 | 299 | 0.444 [0.421, 0.468] | 0.742 | 0.416 |
| 200 to 1,000 nt | ViennaRNA (exact) | 87 | 0.376 [0.325, 0.424] | 0.742 | 0.354 |
| 200 to 1,000 nt | LinearPartition b100 | 87 | 0.376 [0.325, 0.429] | 0.742 | 0.354 |
| 200 to 1,000 nt | LinearSampling b100 | 87 | 0.375 [0.327, 0.423] | 0.742 | 0.353 |
| 200 to 1,000 nt | LinearSampling b500 | 87 | 0.376 [0.323, 0.428] | 0.742 | 0.353 |
| 1,000 to 3,000 nt | ViennaRNA (exact) | 3 | 0.107 [−0.185, 0.291] | 0.590 | n.a. |
| 1,000 to 3,000 nt | LinearSampling b100 | 3 | 0.129 [−0.164, 0.287] | 0.594 | n.a. |

### 3.4. Runtime and memory

On short RNAs the runtime of every tool was dominated by process start-up and stayed below one second up to about 420 nucleotides (Fig. 2) [@ditommaso2017; @zhang2020]. At 2,216 nucleotides, exact ViennaRNA with 1,000 samples took 39 s and 427 MB, LinearSampling with beam 100 took 4.0 s and 67 MB, and LinearPartition at infinite beam took 59 s and 325 MB (Table 3) [@lorenz2011; @zhang2020; @zhang2023]. Fitted log-log slopes were 2.15 for ViennaRNA time and 0.93 for its memory, against 1.45 to 1.55 for LinearSampling time at finite beams [@zhang2023; @mccaskill1990]. These slopes mix start-up overhead with algorithmic cost and come from four concurrent jobs, so they should not be read as asymptotic complexity [@zhang2020; @ditommaso2017].

Table 3. Wall-clock time and peak memory for the two longest RNAs (HOXA9 RNA, 1,266 nt; CMV BN3 RNA, 2,216 nt, mean of two profiles), measured with four jobs in parallel [@lorenz2011; @zhang2020; @zhang2023].

| Tool and beam | Time, 1,266 nt (s) | Memory, 1,266 nt (MB) | Time, 2,216 nt (s) | Memory, 2,216 nt (MB) |
|---|---|---|---|---|
| ViennaRNA exact, 1,000 samples | 8.7 | 141 | 39.0 | 427 |
| LinearPartition b100 | 3.3 | 48 | 6.6 | 97 |
| LinearPartition b500 | 13.8 | 102 | 54.1 | 263 |
| LinearPartition b0 (infinite) | 13.6 | 111 | 59.0 | 325 |
| LinearSampling b100 | 1.9 | 37 | 4.0 | 67 |
| LinearSampling b500 | 9.8 | 94 | 32.0 | 210 |
| LinearSampling b0 (infinite) | 3.2 | 39 | 10.7 | 111 |

![Fig. 2. Runtime (left) and peak memory (right) against sequence length on log-log axes for every tool and beam. Short-sequence times are dominated by process start-up [@zhang2020; @zhang2023].](fig_scaling.png)

### 3.5. Error of linear-time sampling

Below 1,000 nucleotides all beams stayed close to the exact ensemble (Table 4, Fig. 3) [@zhang2020; @zhang2023]. LinearPartition at beam 500 or infinity had a mean $\mathrm{MAD}_{\mathrm{u}}$ of 0.0002 or less, and LinearSampling had a mean $\mathrm{MAD}_{\mathrm{u}}$ of about 0.004 to 0.006 at every beam, which matches the sampling error expected from 1,000 samples [@ding2003; @zhang2023]. Above 1,000 nucleotides beams of 100 and 200 diverged: the mean $\mathrm{MAD}_{\mathrm{u}}$ was 0.124 and 0.086 for LinearPartition and 0.125 and 0.088 for LinearSampling, while beam 500 stayed at 0.008 and 0.011 [@zhang2020; @zhang2023]. With the 0.05 tolerance, beam 100 first failed at 342 nucleotides (a 342-nucleotide HOXA9 RNA construct, $\mathrm{MAD}_{\mathrm{u}}$ 0.055) and beam 200 at 1,266 nucleotides, whereas beam 500 and the infinite beam passed every RNA up to 2,216 nucleotides [@zhang2020; @huang2019].

Table 4. Linear-time error relative to exact ViennaRNA output: mean $\mathrm{MAD}_{\mathrm{u}}$ (Eq. 10) by length stratum, reliable length range (Algorithm 4, tolerance 0.05) and Kendall τ of per-RNA ρ (Eq. 12) [@zhang2020; @zhang2023; @kraemer2006].

| Tool and beam | $\mathrm{MAD}_{\mathrm{u}}$ < 200 nt (n = 328) | $\mathrm{MAD}_{\mathrm{u}}$ 200 to 1,000 nt (n = 108) | $\mathrm{MAD}_{\mathrm{u}}$ 1,000 to 3,000 nt (n = 3) | Reliable up to (nt) | First failure (nt) | Kendall τ |
|---|---|---|---|---|---|---|
| LinearPartition b100 | 0.0001 | 0.0024 | 0.124 | 338 | 342 | 0.993 |
| LinearPartition b200 | 0.0001 | 0.0002 | 0.086 | 418 | 1,266 | 0.994 |
| LinearPartition b500 | 0.0001 | 0.0002 | 0.008 | ≥ 2,216 | none | 0.994 |
| LinearPartition b0 | 0.0001 | 0.0002 | 0.003 | ≥ 2,216 | none | 0.994 |
| LinearSampling b100 | 0.0047 | 0.0062 | 0.125 | 338 | 342 | 0.966 |
| LinearSampling b200 | 0.0048 | 0.0044 | 0.088 | 418 | 1,266 | 0.967 |
| LinearSampling b500 | 0.0047 | 0.0043 | 0.011 | ≥ 2,216 | none | 0.969 |
| LinearSampling b0 | 0.0047 | 0.0045 | 0.007 | ≥ 2,216 | none | 0.973 |

![Fig. 3. Mean absolute difference in unpaired probability between linear-time and exact output against length (left; dashed line, tolerance 0.05) and against beam width (right; 10⁴ marks the infinite beam) [@zhang2020; @zhang2023].](fig_bias.png)

Approximation error did not change which RNAs were predicted well. Kendall's τ between per-RNA correlations from exact and linear-time output was 0.993 to 0.994 for LinearPartition and 0.966 to 0.973 for LinearSampling [@kraemer2006; @rosner2007]. The largest change in an RNA's correlation was 0.05 above 1,000 nucleotides and up to 0.12 for single short RNAs sampled by LinearSampling, where 1,000 samples leave visible noise [@schober2018; @ding2003]. The Jensen-Shannon divergence between sampled structure distributions on RNAs of at most 150 nucleotides was about 0.37 for every beam, including the infinite beam, so it measured the noise floor of two finite samples of 1,000 structures; any beam effect was too small to detect at this sample size [@capra2007; @ding2003].

### 3.6. MERGE-RNA

In the shared-parameter mode, MERGE-RNA fitted $\mu_r = 0.411$ kcal/mol and $p_b = 0.227$ kcal/mol on RNase P and the *Tetrahymena* ribozyme [@sacco2026]. On the three held-out systems, scored on replicate 2 within the trimmed window, its mean Spearman correlation was 0.629 and its AUC 0.804, against 0.660 and 0.837 for exact ViennaRNA on the same positions; MERGE-RNA was higher for hc16 (0.617 against 0.606) and lower for the HCV IRES (0.533 against 0.561) and the glycine riboswitch (0.735 against 0.814) [@sacco2026; @lorenz2011]. Its unpaired probabilities differed from full-length ViennaRNA by a mean of 0.08 to 0.11, mostly because it folds the primer-trimmed sequence [@sacco2026]. The per-RNA mode did not converge within the two-hour limit per RNA set for this study and produced no output [@sacco2026; @incarnato2018].

## 4. Discussion

On natural RNAs up to 1,000 nucleotides, exact ViennaRNA sampling, LinearPartition and LinearSampling agreed with reactivity equally well, and the differences between them were far smaller than the bootstrap intervals [@lorenz2011; @zhang2020; @zhang2023]. This is expected once the consistency test showed that all three use the same energy model, because at these lengths beam pruning removes almost no probability mass [@zhang2020; @huang2019]. The agreement itself was moderate (ρ near 0.44 and AUC near 0.74), in line with the view from EternaBench that thermodynamic parameters leave room for improvement in ensemble-level tasks [@waymentsteele2022; @turner2010]. The decline from the shortest to the 200 to 1,000 stratum is consistent with the known drop in prediction accuracy for longer RNAs, although strata also differ in composition and probe [@mathews2004rna; @lu2009].

The bias audit gives a concrete rule for users of linear-time tools. Beams of 100 and 200 were accurate below about 420 nucleotides but shifted unpaired probabilities by about 0.1 on average for RNAs above 1,000 nucleotides, an error large enough to move individual positions across a reactivity threshold [@zhang2020; @zhang2023]. A beam of 500 kept the mean error near 0.01 up to 2,216 nucleotides, but at that length it was not faster than exact ViennaRNA in our runs, so its advantage only appears at lengths where the exact O(N²) memory becomes limiting [@zhang2020; @mccaskill1990]. Even the inaccurate beams preserved the ranking of RNAs by predictability (τ ≥ 0.97), so conclusions about which RNAs a model handles well are less sensitive to the beam than the probabilities themselves [@kraemer2006; @rosner2007]. For sampling-based estimates, 1,000 samples add an error of about 0.005 in unpaired probability, and structure-level comparisons such as the Jensen-Shannon divergence need far larger samples before they can detect beam effects [@ding2003; @capra2007].

MERGE-RNA in the shared-parameter mode did not outperform ViennaRNA on the three held-out systems, but three RNAs cannot separate the methods, and the soft constraints that give MERGE-RNA most of its fit to data were not evaluated because the per-RNA fits did not converge in time [@sacco2026]. The comparison also has a structural asymmetry: MERGE-RNA folds the primer-trimmed sequence while the other tools fold the full construct, which changes pair probabilities near the ends [@sacco2026; @incarnato2018]. A fair test of probing-fitted models needs more MaP-count data than the five systems available here and a compute budget for the soft-constraint fits [@sacco2026; @zubradt2017].

The main limitation is coverage of long RNAs. Only two distinct RNAs exceeded 418 nucleotides (a 1,266-nucleotide HOXA9 RNA and a 2,216-nucleotide CMV RNA with two profiles), so the reliable-length values for beams of 100 and 200 fall inside a gap of 418 to 1,266 nucleotides that contains no data, and the long-RNA error bound is preliminary [@cordero2012; @zhang2020]. Ribonanza, RASP and bpRNA-1m could not be reached from the computing environment, and Ribonanza constructs are in any case windows of at most about 130 nucleotides [@he2024; @waymentsteele2022]. Other limits apply to every tool equally: the energy model excludes pseudoknots, predictions were made at 37 °C while many RMDB experiments were done near 24 °C, and DMS and SHAPE reactivities also depend on three-dimensional context and buffer [@rivas1999; @lyngso2000; @sacco2026]. Base-pair scores used p > 0.5 structures instead of maximum expected accuracy structures and reference annotations of variable quality, so they are secondary to the reactivity scores [@lu2009; @mathews2004rna]. TVAE-RNA, the only learned generator considered, could not be run, so the benchmark does not yet cover deep-learning ensembles [@mei2025].

The pipeline is designed to remove these limits as data become available. Every stage reads its inputs from versioned files, every threshold is fixed in one configuration file, and new probing archives or generators enter as new inputs to stage A or new processes in stage D [@ditommaso2017; @wilkinson2016]. Adding RNAs between 400 and 30,000 nucleotides, for example from transcriptome-wide probing sets, would turn the provisional error bound into a firm one and test the regime where linear-time methods are needed [@zhang2023; @siegfried2014].

## 5. Conclusions

On 439 natural RNAs with chemical mapping data, exact and linear-time ViennaRNA-model ensemble generators agreed with reactivity to the same moderate degree up to 1,000 nucleotides. Linear-time sampling with beams of 100 or 200 is accurate for RNAs below about 400 nucleotides but distorts unpaired probabilities by about 0.1 above 1,000 nucleotides, whereas a beam of 500 stays within about 0.01 of the exact ensemble up to 2,216 nucleotides. The ranking of RNAs by predictability was stable under every beam. Because only two distinct RNAs exceeded 418 nucleotides, the long-RNA bound should be treated as provisional until longer natural RNAs are added. The pipeline, configuration and results are versioned so the benchmark can be extended with new data and generators.

## CRediT authorship contribution statement

Arli Aditya Parikesit: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Writing – review & editing, Visualization.

## Declaration of competing interest

The author declares no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper. [Author to confirm.]

## Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the author used Claude (Anthropic) to write and run analysis code and to draft the manuscript text. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the publication.

## Data and code availability

Code, configuration, processed data, per-RNA ensemble files and result tables are in the repository github.com/Johnkecops/gap4-rna-ensemble-benchmark, which will be made public on acceptance and is available to reviewers on request. Raw RDAT files are redistributed from RMDB, and the commits of all external tools are recorded in environments/tool_commits.tsv.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors. [Author to confirm.]

## Acknowledgements

The author thanks the contributors of the RNA Mapping Database and the developers of the MERGE-RNA, ViennaRNA, LinearPartition and LinearSampling software for making their data and code openly available.

## References

::: refs
:::
