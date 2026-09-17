# Benchmark: Kostic et al. (2015), infant gut microbiome and type 1 diabetes

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **S3** of the validation matrix, the
repeated-measures guard on shotgun data, and the study that shows what a curated
reprocessing database costs.

| | |
|---|---|
| **Paper** | Kostic AD, Gevers D, Siljander H, et al. The dynamics of the human infant gut microbiome in development and in progression toward type 1 diabetes. *Cell Host Microbe.* 2015;17(2):260-273. [10.1016/j.chom.2015.01.001](https://doi.org/10.1016/j.chom.2015.01.001) |
| **Raw data** | DIABIMMUNE, via SRA; reprocessed in curatedMetagenomicData |
| **Analyzed** | `curatedMetagenomicData`, `2021-03-31.KosticAD_2015.relative_abundance`, `counts = TRUE`, 330 species × 120 samples |
| **Design** | Longitudinal infant cohort, repeated stool sampling through infancy. 19 subjects with up to 10 samples each, of whom 4 are T1D progressors |
| **Hypothesis** | The infant gut microbiome develops dynamically through infancy, and progression to type 1 diabetes is marked by a drop in diversity between seroconversion and diagnosis |
| **Verdict** | **Incomplete.** 1 of 1 scorable claim reproduced, 3 outside the denominator. No sweep |

---

## 1. Study card

**9 fields stated · 2 inferred · 12 absent.**

The stated fields are unusually specific about the cohort and the sampling
schedule, because the paper's argument is about timing. The absent fields are the
downstream analytical choices.

One field matters more than the rest and is neither stated nor absent but
*changed*: the cohort. See section 2.

### The paper's four claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `diversity-develops` | Diversity increases through infancy | the paper's framing claim about dynamic development |
| `alpha-drop-progressors` | Diversity drops in progressors between seroconversion and diagnosis | a 25% drop in alpha diversity in that window |
| `strain-stable-within` | Strains are stable within an individual | stated in the results |
| `pathways-constant` | Pathway abundance is comparatively constant | stated in the results |

---

## 2. Reconstruction

330 species × 120 samples. The repeated-measures guard fired correctly:

> Repeated measures detected: `subject_id` has up to 10 samples per level across
> 19 levels.

This is the shape test from defect #13 working in the intended direction. 19
levels holding at most 10 samples each is a subject identifier, so it is
accepted, where A3's 4-level `host` column was rejected.

Group sizes: control 89 samples across 15 subjects, T1D 31 samples across 4
subjects. Sparsity 82.9%.

### The cohort is a subset, and it decides what can be tested

The paper studied **33 infants**. curatedMetagenomicData holds **19**, of whom
**4** are T1D progressors.

That is not a defect in either source. cMD reprocessed what was deposited. But it
determines what can be tested, because the paper's headline claim is a
*within-window* statement: diversity drops in progressors specifically between
seroconversion and diagnosis. Partitioning the four progressors' samples by that
window gives:

| Window | Samples |
|---|---|
| before seroconversion | 16 |
| **seroconversion to diagnosis** | **3** |
| after diagnosis | 12 |

**Three samples.** Median Shannon in that window is 2.957 against 2.378 before
seroconversion, numerically the opposite direction, Wilcoxon p = 0.359.

Reporting that as a failure to reproduce would be wrong. Three samples from at
most three infants cannot test a claim built on 33. The correct answer is that
the claim is **not adjudicable from this subset**, and the scoring machinery
recorded it as such rather than scoring it.

---

## 3. Results

### 3.1 Diversity increases through infancy

Mixed model, `Shannon ~ study_condition + infant_age + (1 | subject_id)`:

| Term | Estimate | p |
|---|---|---|
| `infant_age` | **+0.00055 per day** | **6.5e-04** |
| `study_condition` (T1D) | +0.092 | 0.641 |
| **ICC** | **0.360** | |

Diversity rises with age across infancy, about +0.20 Shannon per year, which is
the paper's framing claim about dynamic development. It reproduces.

**36% of Shannon variance sits between subjects.** With up to 10 samples per
infant, ignoring that structure would be a substantial error. Here it does not
change the conclusion, because the naive Wilcoxon on condition at p = 0.905 and
the mixed model at p = 0.641 agree there is no overall group difference. The
agreement is a result, not an assumption, and it could only be established by
fitting both.

### 3.2 No overall diversity difference between groups

Median Shannon: control 2.601, T1D 2.510. The paper makes no such claim, so this
is an observation rather than a claim test. It is consistent with the paper's
position that the difference is confined to a specific time window.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Reproduced (1):
  + diversity-develops : diversity increases with infant age, p = 6.5e-04

Not reproduced (0)
Not licensed  (0)

Outside the denominator (3):
  ? alpha-drop-progressors : contrast_mismatch - claim is about the
                             seroconversion window; the evidence describes
                             infant_age, and the window holds 3 samples
  ? strain-stable-within   : schema_gap        - no claim type expresses a
                             strain-level stability claim
  ? pathways-constant      : not_attempted     - HUMAnN pathway data not
                             downloaded

Scorable: 1 of 4.  Reproduced: 1 of 1.
```

**This is the lowest adjudicability in the series: 1 claim of 4.** That is the
study's most important result, and it is a finding about data availability rather
than about the paper.

### The contrast guard did its job

`alpha-drop-progressors` was scored `contrast_mismatch`, not `not_reproduced`,
with the reason stated:

> claim is about 'seroconversion_window'; evidence describes 'infant_age'

The analysis that *was* run measured diversity against age across the whole
cohort. Scoring a window-specific claim against it would have produced a
confident answer to a question the analysis never asked. This is the constraint
added after S1 working exactly as intended, on a case where the wrong answer
would have looked entirely reasonable: the numbers exist, they point the other
way, and quoting them would have read as a contradiction of the paper.

### What this means

One claim of four could be tested, and it reproduced. Three could not, for three
different reasons, and pooling them would hide the distinction. The window claim
fails on cohort coverage, the strain claim on what the instrument can express,
and the pathway claim on work not done in this run. Only the last is under the
analyst's control.

The broader point for the series: **any benchmark drawing on a curated
reprocessing database inherits its cohort.** The data is whatever was deposited
and reprocessed, not what the paper analyzed, and a claim about a subgroup is the
first thing that becomes untestable.

### Limitations

1. **19 of 33 subjects.** The curated subset cannot test the paper's central
   window claim.
2. **Reconstructed counts.** Relative abundance times read depth, so no true
   sampling variance. Depth values run 14 to 86 million.
3. **Strain-level claims not tested.** Would need the `marker_presence` product.
4. **Pathway claims not tested.** Would need `pathway_abundance`.
5. **No sweep.** Single reconstruction.
6. **Only 4 T1D subjects.** Any group contrast here is underpowered regardless of
   method, which is why section 3.2 is reported as an observation and not scored.

---

## 5. Defects this run exposed

**None.** S3 is the first study in the series to run end to end without exposing
a bug, which is itself informative after four consecutive studies that each found
two or three.

The fixes from S1 and S4 carried. The MetaPhlAn path, the `counts = TRUE`
reconstruction, the repeated-measures detector and the contrast guard all behaved
correctly on data none of them had seen.

**The mixed-model path works on shotgun data.** This had never been exercised;
every prior mixed model in the project ran on 16S counts. `fit_lmm` with
`(1 | subject_id)` fitted cleanly on reconstructed MetaPhlAn counts and returned
a sensible ICC.

---

## 6. Reproducing this

```r
se <- curatedMetagenomicData("2021-03-31.KosticAD_2015.relative_abundance",
                             dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]
rd  <- as.data.frame(SummarizedExperiment::rowData(se))
tax <- as.matrix(rd[, c("superkingdom","phylum","class","order","family","genus","species")])
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

ps <- phyloseq(otu_table(SummarizedExperiment::assay(se), taxa_are_rows = TRUE),
               tax_table(tax),
               sample_data(as.data.frame(SummarizedExperiment::colData(se))))
ps <- mfg_set_normalization(ps, "raw")
validate_inputs(ps, group_var = "study_condition", subject_var = "subject_id")

a  <- compute_alpha(tss_transform(ps), measures = "Shannon")
df <- data.frame(Shannon = a$Shannon, mfg_meta(ps)[rownames(a), ])
fit_lmm(Shannon ~ study_condition + infant_age, data = df,
        random = "(1 | subject_id)")
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
