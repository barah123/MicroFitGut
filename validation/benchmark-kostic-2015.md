# Benchmark: Kostic et al. (2015), infant gut microbiome and type 1 diabetes

Study **S3** of the validation matrix — the repeated-measures guard on shotgun
data, and the study that shows what a curated reprocessing database costs.

| | |
|---|---|
| **Paper** | Kostic AD, Gevers D, Siljander H, et al. The dynamics of the human infant gut microbiome in development and in progression toward type 1 diabetes. *Cell Host Microbe.* 2015;17(2):260–273. |
| **Data** | curatedMetagenomicData `2021-03-31.KosticAD_2015.relative_abundance`, `counts = TRUE` |
| **Design** | Longitudinal infant cohort, DIABIMMUNE; repeated stool sampling through infancy |
| **Verdict** | **Incomplete** — 1 claim held, 0 failed, 1 not evaluable, 2 unadjudicated |

---

## 1. Reconstruction

330 species × 120 samples, **19 subjects with up to 10 samples each**. The
repeated-measures guard fired correctly:

> Repeated measures detected: `subject_id` has up to 10 samples per level across
> 19 levels.

This is the shape test from defect #13 working in the intended direction: 19
levels holding at most 10 samples each is a subject identifier, so it is
accepted, where A3's 4-level `host` column was rejected.

Group sizes: control 89 samples (15 subjects), T1D 31 samples (4 subjects).
Sparsity 82.9%.

---

## 2. The cohort is a subset, and it matters

The paper studied **33 infants**. curatedMetagenomicData holds **19**, of whom
**4** are T1D progressors.

That is not a defect in either source — cMD reprocessed what was deposited — but
it determines what can be tested. The paper's headline claim is a *within-window*
statement: diversity drops in progressors specifically between seroconversion and
diagnosis. Partitioning the 4 progressors' samples by that window gives:

| Window | Samples |
|---|---|
| before seroconversion | 16 |
| **seroconversion → diagnosis** | **3** |
| after diagnosis | 12 |

**Three samples.** Median Shannon in that window is 2.957 against 2.378 before
seroconversion — numerically the opposite direction, Wilcoxon p = 0.359.

Reporting that as a failure to reproduce would be wrong. Three samples from at
most three infants cannot test a claim built on 33. The correct answer is that
the claim is **not evaluable from this subset**, and the scoring machinery
recorded it as such rather than scoring it.

---

## 3. Results

### 3.1 Diversity increases through infancy — held

Mixed model, `Shannon ~ study_condition + infant_age + (1 | subject_id)`:

| Term | Estimate | p |
|---|---|---|
| `infant_age` | **+0.00055 per day** | **6.5e-04** |
| `study_condition` (T1D) | +0.092 | 0.641 |
| **ICC** | **0.360** | — |

Diversity rises with age across infancy — about +0.20 Shannon per year — which
is the paper's framing claim about dynamic development. It holds.

**36% of Shannon variance sits between subjects.** With up to 10 samples per
infant, ignoring that structure would be a substantial error. Here it does not
change the conclusion, because the naive Wilcoxon on condition (p = 0.905) and
the mixed model (p = 0.641) agree there is no overall group difference — but the
agreement is a result, not an assumption.

### 3.2 No overall diversity difference between groups

Median Shannon: control 2.601, T1D 2.510. The paper makes no such claim, so this
is an observation rather than a claim test — and it is consistent with the
paper's position that the difference is confined to a specific time window.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Claims that held:
  + diversity-develops : diversity increases with infant age, p = 6.5e-04

Claims NOT adjudicated:
  ? alpha-drop-progressors : NOT EVALUABLE — claim is about the
                             seroconversion window; the evidence supplied
                             describes infant_age
  ? strain-stable-within   : strain-level analysis not performed
  ? pathways-constant      : requires HUMAnN pathway data, not downloaded
```

### The contrast guard did its job

`alpha-drop-progressors` was scored `not_evaluable`, not `failed`, with the
reason stated:

> claim is about 'seroconversion_window'; evidence describes 'infant_age'

The analysis that *was* run measured diversity against age across the whole
cohort. Scoring a window-specific claim against it would have produced a
confident answer to a question the analysis never asked. This is the constraint
added after S1 working exactly as intended, on a case where the wrong answer
would have looked entirely reasonable.

---

## 5. What this study adds to the series

**No new defects.** S3 is the first study in the series to run end to end without
exposing a bug — which is itself informative after four consecutive studies that
each found two or three.

The fixes from S1 and S4 carried: the MetaPhlAn path, the `counts = TRUE`
reconstruction, the repeated-measures detector and the contrast guard all behaved
correctly on data none of them had seen.

**The mixed-model path works on shotgun data.** This had never been exercised —
every prior mixed model in the project ran on 16S counts. `fit_lmm` with
`(1 | subject_id)` fitted cleanly on reconstructed MetaPhlAn counts and returned
a sensible ICC.

---

## 6. Limitations

1. **19 of 33 subjects.** The curated subset cannot test the paper's central
   window claim. Any benchmark drawing on the cMD arm inherits this: the cohort
   is whatever was deposited and reprocessed, not what the paper analysed.
   (`METHODS-benchmarking.md` §3.6 item 5.)
2. **Reconstructed counts.** Relative abundance × read depth, so no true sampling
   variance. Depth values run 14–86 million.
3. **Strain-level claims not tested.** Would need the `marker_presence` product.
4. **Pathway claims not tested.** Would need `pathway_abundance`; S5 covers the
   HUMAnN path.
5. **No sweep.** Single reconstruction.
6. **Only 4 T1D subjects.** Any group contrast here is underpowered regardless of
   method.

---

## 7. Reproducing this

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
