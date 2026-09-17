# Benchmark: Pérez-Losada et al. (2020), the "Grandma Hypothesis"

A full benchmark-mode run against a published study, alongside the problem-set
and quiz regression tests. Where those check MicroFitGut against a known answer
key, this checks it against a paper, reconstructing the analysis, scoring the
paper's stated claims, and reporting which survive.

| | |
|---|---|
| **Paper** | Pérez-Losada M, et al. Testing the "Grandma Hypothesis": characterizing skin microbiome diversity as a project-based learning approach to genomics. *J Microbiol Biol Educ.* 2020;21(1):21.1.7. [10.1128/jmbe.v21i1.2010](https://doi.org/10.1128/jmbe.v21i1.2010) |
| **Raw data** | NCBI SRA BioProject [PRJNA553551](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA553551) |
| **Analyzed** | `dataset/ps.RDS`, the processed feature table, 619 taxa × 344 samples, 128 participants, with phylogeny |
| **Hypothesis** | Body regions washed less often ("grandma hotspots": behind the ears, between the toes, navel) host different microbial communities from regions washed more often (forearms, calves) |
| **Verdict** | **Partially reproduced** (sweep: 16 configurations, all executed) |

---

## 1. Study card

Extracted from the paper with a verbatim span and location for every stated
field.

**9 fields stated · 3 inferred · 12 absent.**

The twelve absent fields are a finding in their own right. The paper does not
state its prevalence filter, minimum depth, normalization, rarefaction depth,
beta-diversity distance, alpha metric, FDR method, covariates, or random
effects. Each becomes an axis of the sensitivity sweep, because each is a place
where a careful reader could reasonably choose differently.

This is not unusual for the genre, the study is a teaching paper, and students
chose their own tests within MicrobiomeAnalyst. But it means a single
"replication" of this analysis does not exist; a family of defensible ones does.

### The paper's three claims

| id | Claim | Evidence quoted from the paper |
|---|---|---|
| `alpha-higher-in-washed` | More-washed sites have higher alpha diversity | "sites cleaned more frequently (calves and forearms) have significantly (ANOVA F-value > 3.42; p < 0.014) higher alpha-diversity" |
| `beta-differs` | Communities differ between more- and less-washed sites | "differ significantly (PERMANOVA F-value > 6.017; p < 0.001) in beta-diversity" |
| `ten-genera` | Ten genera differ in relative abundance | "ten genera … varied significantly (Log LDA score > 1; p < 0.05)" |

---

## 2. Reconstruction

Quality control at prevalence 0.10 and minimum depth 5,000: **619 → 150 taxa,
344 → 307 samples.**

Intake flagged three things before any test ran:

- Depth varies **124.6-fold** (1,107 to 137,913), so normalization is not optional.
- Sparsity **91.8%**.
- **Repeated measures**: `patient`, up to 5 samples each across 128 participants.

The third is the crux. The paper's analyses treat samples as independent.

---

## 3. Results

### Alpha diversity: the claim holds

Mean Shannon: **less-washed 1.13, more-washed 2.76.**

| Approach | Result |
|---|---|
| Ignoring participant (the paper's approach) | Wilcoxon, p < 1e-16 |
| Simple blocked test | **Refused** |
| Mixed model with `(1 \| patient)` | effect **+1.62** Shannon, p < 1e-16, ICC 0.161 |

MicroFitGut refused the simple blocked test with a specific reason rather than a
generic warning:

> `patient` indicates repeated measures but the design is not a complete block
> (cell counts range 0–3). No simple test is valid: Kruskal-Wallis would treat
> correlated samples as independent and overstate significance, and Friedman
> requires exactly one observation per cell. Fit a mixed model with
> `(1 | patient)` instead.

Fitting that model, **the claim survives comfortably.** 16.1% of Shannon
variance sits between participants, so the correlation is real, but the effect
is large enough that modeling it changes nothing about the conclusion.

**This is the paper's central finding, and it reproduces.**

### Beta diversity: the claim is not licensed

| Test | Result |
|---|---|
| PERMANOVA, unstratified | R² = 0.2013, F = 76.86, p = 0.001 |
| PERMANOVA, permuted within participant | R² = 0.2013, p = 0.001 |
| **betadisper** | **p = 0.001, dispersion is heterogeneous** |

The significance is not the problem: stratifying by participant leaves it
unchanged, so the repeated-measures concern does not overturn this one.

The dispersion check does. PERMANOVA responds to differences in *spread* as well
as differences in *location*. When within-group dispersion differs significantly
between groups, a significant PERMANOVA cannot be read as "the communities
differ in composition", the test cannot separate the two explanations.

The paper reports no dispersion check. MicroFitGut therefore records this claim
as **unadjudicated rather than failed**: the reanalysis does not contradict it,
it establishes that the evidence presented does not support it as stated. A
different design, or a test robust to heteroscedasticity, could still
vindicate it.

### Differential abundance: the count holds, the identities do not

ANCOM-BC2 at genus level with `(1 | patient)`: **11 robust hits of 80 tested**
(16 before the sensitivity filter).

The paper reported ten genera. **Eleven versus ten is agreement on count.**

Agreement on identity is much weaker, **3 of 10 recovered (30%), Jaccard 0.176:**

| | |
|---|---|
| **Recovered (3)** | *Staphylococcus*, *Streptococcus*, *Escherichia-Shigella* |
| **Not recovered (7)** | *Lactobacillus*, *Bacillus*, *Micrococcus*, *Pseudomonas*, *Lawsonella*, *Acinetobacter*, *Enhydrobacter* |
| **Never tested (1 of those 7)** | *Bacillus*, removed by the prevalence filter, so a filtering difference rather than a disagreement about significance |
| **New in reanalysis (7)** | *Corynebacterium*, *Veillonella*, *Actinomyces*, *Haemophilus*, *Ferruginibacter*, *Reyranella*, SM1A02 |

Two causes are clearly in play and the sweep separates them. The paper used
LEfSe with a Log LDA threshold; the reanalysis used ANCOM-BC2 with a participant
random effect. Per the discrepancy rubric, `statistics.da_method` predicts
moderate recovery with agreement on direction, and `statistics.repeated_measures`
predicts that the independence-assuming analysis reports more. Both fit.

---

## 4. Verdict

```
=== Benchmark verdict: PARTIALLY REPRODUCED ===

Reproduced (2):
  + alpha-higher-in-washed : +1.62 Shannon, p < 1e-16, LMM with (1 | patient)
  + ten-genera             : paper approx 10, reanalysis 11 (tolerance 3)

Not reproduced (0)

Not licensed (1):
  ! beta-differs           : significant at p = 0.001, but betadisper
                             p = 0.001, so PERMANOVA does not license a
                             composition claim

Outside the denominator (0)

Scorable: 3 of 3.  Reproduced: 2 of 3.
```

**This verdict block was rewritten when the scoring scheme was revised.** Two
changes:

1. **The DA recovery item is no longer a claim.** It previously appeared as the
   study's only failure, reading "only 30% of published taxa recovered". That is
   a concordance measure against a hard-coded threshold of 0.7, not something the
   paper asserted, and under the current scheme a concordance measure cannot
   decide an outcome. The finding itself is unchanged and is reported in section
   3 and in the sweep, where it belongs.
2. **`beta-differs` moved from unadjudicated to `not_licensed`.** It had been
   pooled with claims that could not be tested at all. The test ran and returned
   a result; what is missing is the licence to read that result as a composition
   claim. That is a finding about the paper, and it now sits inside the
   denominator rather than beside it.

The `ten-genera` verdict rests on a **tolerance of 3 against a claim of 10**, a
window of 7 to 13, which was an analyst choice and is recorded here because it
decides the outcome. The sweep in section 4b shows the same claim holding in only
3 of 16 configurations.

### What this means

The headline of the paper, that less-washed skin regions carry less diverse
microbial communities, **holds**, and holds under a model that accounts for
participants contributing several samples each.

The beta-diversity claim is **not contradicted but not supported as stated**. The
test used cannot distinguish a shift in community composition from a difference
in variability, and the paper does not report the check that would separate them.

The biomarker list is **the fragile part**. The number of differing genera
reproduces; which genera they are largely does not. Seven of the ten named do
not survive an analysis that models the repeated sampling, and seven genera the
paper does not mention appear instead.

None of this says the paper is wrong. It is a teaching paper reporting what
students found with the tools they were given, and it says so. What the benchmark
establishes is which of its conclusions are robust to defensible analytical
choices and which depend on them, and that is a different, more useful question.

---

## 4b. Sensitivity sweep: 16 configurations

The paper leaves 15 analytical parameters unstated or only inferable. The sweep
varies them one at a time from the reconstruction baseline, so each axis's effect
is isolated. **All 16 configurations executed**, including weighted and unweighted
UniFrac, this dataset ships a phylogeny with 619 tips for 619 taxa, asserted as
genuine via `mfg_mark_tree_real()` and recorded in the run log.

| Run | Recovery | Significant genera |
|---|---|---|
| baseline (ANCOM-BC2, prevalence 0.10) | 0.40 | 19 |
| `da_method=aldex2` | 0.60 | 40 |
| `da_method=deseq2` | 0.60 | 33 |
| `da_method=kruskal` | 0.60 | 60 |
| `random_effects=subject` | **0.20** | **9** |
| `prevalence_filter=0` | **0.80** | 43 |
| `prevalence_filter=0.05` | 0.70 | 33 |
| `prevalence_filter=0.2` | 0.30 | 8 |
| `fdr_method=holm` | 0.30 | 10 |
| `fdr_method=BY` | 0.30 | 15 |
| `min_depth`, `rarefaction_depth`, `beta_distance` | 0.40 | 18–19 |

Across defensible choices, recovery of the paper's genera ranges from **20% to
80%**, and the count of significant genera from **8 to 60**.

### Attribution

| Metric | Dominant axis | Spread | Second | Ratio | Confidence |
|---|---|---|---|---|---|
| Significant genera | `da_method` | 41 | `prevalence_filter` (35) | 1.17 | **weak, not distinguishable** |
| Recovery rate | `prevalence_filter` | 0.50 | `random_effects` (0.20) | 2.50 | **clear** |

The tool names a cause for one metric and refuses for the other in the same run.
For the count of significant genera, the method and the filter move it almost
equally, so no single cause is named. For recovery of the paper's named genera,
the prevalence filter dominates by a factor of 2.5 and is named.

### What the sweep changes about the interpretation

**The prevalence filter, not the statistics, drives most of the disagreement.**
At the baseline filter of 0.10, three of the paper's ten genera are recovered. At
no filter, **eight of ten are**. The genera that "failed to reproduce" were
largely never tested, they are rare enough to be removed before any test ran.

This is the `filtering.prevalence` signature in the discrepancy rubric exactly:
*missed taxa turn out never to have been tested*. It is a different finding from
a statistical disagreement, and a much less damaging one for the paper. The
earlier section of this report, written from the baseline alone, understated
this, which is precisely why the sweep exists.

**Modelling repeated measures costs the most of any single correction.** Adding
`(1 | patient)` halves the significant genera, 19 → 9, and halves recovery,
0.40 → 0.20. It is the only axis that moves results in the conservative
direction as strongly as the filter moves them in the permissive direction.

### Configuration fragility

The paper's count claim, approximately ten genera, **holds in 3 of 16
defensible configurations (19%)**, and only under choices more conservative than
the baseline: a participant random effect, Holm correction, or a strict
prevalence filter.

That is the honest summary of this claim: not refuted, but reproducible only in
a minority of reasonable analyses, and the specific ten genera named are mostly
a product of where the prevalence threshold was set.

---

## 5. Defects this run exposed

Benchmark mode had never been run against a published study before this. Doing so
found eight defects, all fixed. Every one biased results toward reporting
**non-reproduction**, which is the direction that flatters the tool.

| # | Defect | Effect |
|---|---|---|
| 1 | Sweep grid scheduled UniFrac on datasets with no tree | Runs that could only fail, displacing valid ones |
| 2 | Harmonisation applied synonyms before truncating to genus | Real matches missed; recovery understated |
| 3 | Recovery and Jaccard used different matchers | One output reported 75% recovery and "0 shared" |
| 4 | The verdict could not see the paper's claims | Taxon overlap silently stood in for whether a claim held |
| 5 | Sweep baseline paired a normalization with a method that refuses it | **All 15 configurations refused; DA benchmarking never worked** |
| 6 | Empty results crashed instead of scoring zero | The strongest form of non-reproduction was unrecordable |
| 7 | ANCOM-BC2 emits no `significant` column, unlike ALDEx2 and DESeq2 | Code treating results interchangeably silently reads zero |
| 8 | `Escherichia/Shigella` ≠ `Escherichia-Shigella` | Counted as both a miss and a new finding; recovery here 20% → 30% |

Defects 1–6 are fixed in the shipped code. Defect 7 is a documented design
difference, not a bug, but it is a trap and is noted here. Defect 8 is fixed:
compound genus labels now normalize their separator.

---

## 6. Reproducing this

```r
source("skills/microfitgut/scripts/00-packages.R")
for (f in list.files("skills/microfitgut/scripts", "[.]R$", full.names = TRUE)) source(f)

ps <- readRDS("dataset/ps.RDS")
ps <- filter_prevalence(prune_low_depth(ps, 5000), 0.10)
ps <- mfg_set_normalization(ps, "raw")

# alpha, on rarefied counts
psr <- rarefy_to_depth(ps, min(phyloseq::sample_sums(ps)))
a   <- compute_alpha(psr, measures = c("Shannon", "Observed"))
df  <- data.frame(a, mfg_meta(psr)[rownames(a), c("wash", "patient")])
fit_lmm(Shannon ~ wash, data = df, random = "(1 | patient)")

# beta, on relative abundance
d <- compute_distance(tss_transform(ps), "bray")
run_permanova(tss_transform(ps), ~ wash, dist_obj = d, strata = "patient")
check_dispersion(d, mfg_meta(ps)$wash)

# differential abundance, on raw counts at genus level
psg <- mfg_set_normalization(phyloseq::tax_glom(ps, "Genus", NArm = FALSE), "raw")
da_ancombc2(psg, fix_formula = "wash", group = "wash",
            rand_formula = "(1 | patient)", prv_cut = 0.10)
```

Seeds are fixed at 42 throughout. Every number above is written to
`run_log.csv` in the run directory.
