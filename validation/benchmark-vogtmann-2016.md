# Benchmark: Vogtmann et al. (2016), colorectal cancer and shotgun metagenomics

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **S1** of the validation matrix, the
first shotgun benchmark, and the first study in the series to reproduce in full.

| | |
|---|---|
| **Paper** | Vogtmann E, Hua X, Zeller G, et al. Colorectal Cancer and the Human Gut Microbiome: Reproducibility with Whole-Genome Shotgun Sequencing. *PLOS ONE.* 2016;11(5):e0155362. [10.1371/journal.pone.0155362](https://doi.org/10.1371/journal.pone.0155362) |
| **Raw data** | ENA [PRJEB12449](https://www.ebi.ac.uk/ena/data/view/PRJEB12449) |
| **Analyzed** | `curatedMetagenomicData`, `2021-03-31.VogtmannE_2016.relative_abundance`, `counts = TRUE`, 540 species × 110 samples |
| **Design** | Case-control: 52 CRC cases, 52 controls, Washington DC 1985–87, frequency matched on sex and BMI |
| **Hypothesis** | Associations between the gut microbiome and colorectal cancer found by 16S are recovered when the same samples are sequenced by whole-genome shotgun |
| **Verdict** | **Reproduced.** 8 of 8 claims reproduced, 0 not reproduced, 0 unscorable. No sweep |

---

## 1. Study card

**12 fields stated · 0 inferred · 12 absent.**

Unusually complete on stated fields. This is a methodologically careful paper
with a detailed statistics section. The twelve absent fields are all downstream
choices the paper had no occasion to state because it did not perform those
analyses: prevalence filter, normalization for the comparison, FDR method,
beta-diversity distance.

The paper is itself a reproducibility study, which makes it a sharper test of the
benchmark machinery than a conventional paper. Most of its claims are
**negative**, of the form "we did not reproduce X", and a negative claim can fail
in both directions: the reanalysis can find an effect where the paper found none,
or it can fail to test the quantity at all.

### The paper's eight claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `alpha-shannon-ns` | Shannon does not differ by case status | Kruskal-Wallis p = 0.252 |
| `alpha-richness-ns` | Richness does not differ by case status | Kruskal-Wallis p = 0.157 |
| `alpha-evenness-ns` | Evenness does not differ by case status | Kruskal-Wallis p = 0.438 |
| `fusobacterium-cases` | *Fusobacterium* is more often present in cases | 75.0% of cases vs 48.1% of controls, p = 0.006 |
| `porphyromonas-cases` | *Porphyromonas* is more often present in cases | 61.5% vs 40.4%, p = 0.032 |
| `atopobium-null` | *Atopobium* does not differ by case status | 53.8% vs 44.2%, p = 0.328 |
| `clostridia-null` | Clostridia relative abundance does not differ | 33.9% case vs 39.0% control, p = 0.092 |
| `bacteroidia-dominant` | Bacteroidia is the most abundant class | 53.2% case, 50.9% control |

Values are the figures the paper reports, not verbatim prose. Two phrases are
quoted directly below where the wording carries the claim.

### One upstream difference, recorded rather than swept

| | Paper | Reanalysis |
|---|---|---|
| Profiler | **MOCAT / mOTU** | **MetaPhlAn**, via curatedMetagenomicData |
| Alpha diversity input | mOTU relative abundance, downsampled to 2,000 inserts | MetaPhlAn relative abundance |
| Counts | true read counts | reconstructed, relative abundance × read depth |

This is a genuine upstream divergence and the single most important
qualification on everything below. The reanalysis does **not** reproduce the
paper's pipeline. It reanalyzes the same samples through a different taxonomic
profiler. In the discrepancy rubric this is closest to
`taxonomy.database_version`, whose predicted signature is *detection rates differ
while conclusions hold*, and that is what was observed.

---

## 2. Reconstruction

540 species × 110 samples on intake. Six samples carry no `study_condition` and
were **excluded explicitly and logged** before analysis, leaving **52 CRC and 52
control**, matching the paper exactly.

Intake flagged three things before any test ran:

- Depth varies **34.6-fold**, 4.1M to 141.9M reconstructed counts.
- Sparsity **78.7%**.
- **No repeated measures**: one sample per subject across 110 subjects. This
  confirms S1's placement in the cross-sectional cell of the validation matrix,
  and it is the reason no mixed model appears anywhere below.

A character column takes an alphabetical reference level unless set deliberately,
so `study_condition` was set explicitly rather than left to default.

---

## 3. Results

### 3.1 Alpha diversity: all three null results reproduce

| Index | Paper (Kruskal-Wallis) | Reanalysis | Conclusion |
|---|---|---|---|
| Shannon | p = 0.252 | p = 0.899 | both null |
| Richness | p = 0.157 | **p = 0.143** | both null |
| Evenness | p = 0.438 | p = 0.082 | both null |

Richness agrees almost exactly. Evenness moves closer to the threshold but does
not cross it. The paper's observation that "in general the controls had slightly
higher alpha diversity" also holds for evenness, 0.320 control against 0.280
case, though not for Shannon.

### 3.2 Presence and absence: every association reproduces, at half the detection

| Taxon | Paper: case % / control % / p | Reanalysis: case % / control % / p |
|---|---|---|
| Fusobacteria (phylum) | 76.9 / 48.1 / **0.003** | 40.4 / 17.3 / **0.017** |
| *Fusobacterium* (genus) | 75.0 / 48.1 / **0.006** | 40.4 / 17.3 / **0.017** |
| *Porphyromonas* | 61.5 / 40.4 / **0.032** | 28.8 / 5.8 / **0.004** |
| *Atopobium* | 53.8 / 44.2 / 0.328 (ns) | 7.7 / 3.8 / 0.674 (ns) |

Every conclusion holds: three significant associations in the same direction, one
null result still null.

**But detection rates are roughly halved.** Fusobacteria is detected in 76.9% of
cases by mOTU and 40.4% by MetaPhlAn. *Atopobium* falls from 53.8% to 7.7%.
MetaPhlAn requires marker-gene evidence before calling a species present, which
is a stricter criterion than mOTU's. The prevalence estimates are therefore not
interchangeable between profilers even though the case-control comparisons agree.

This matters beyond this study. **A paper reporting a prevalence figure is
reporting a property of its profiler as much as of its cohort.**

### 3.3 Relative abundance: both null results reproduce, closely

| Taxon | Paper (Wilcoxon) | Reanalysis |
|---|---|---|
| Clostridia (class) | 33.9% case / 39.0% control, p = 0.092 | 36.6% / 41.7%, **p = 0.079** |
| Bacteroidia, highest abundance class | 53.2% case / 50.9% control | **50.7%, top class** |

Clostridia reproduces to within 0.013 on the p-value and preserves the direction
the paper noted, that it "tended to be lower in cases". Bacteroidia is confirmed
as the most abundant class.

---

## 4. Verdict

```
=== Benchmark verdict: REPRODUCED ===

Reproduced (8):
  + alpha-shannon-ns     : paper: none; reanalysis: none (p = 0.899)
  + alpha-richness-ns    : paper: none; reanalysis: none (p = 0.143)
  + alpha-evenness-ns    : paper: none; reanalysis: none (p = 0.082)
  + fusobacterium-cases  : higher in CRC, p = 0.017
  + porphyromonas-cases  : higher in CRC, p = 0.004
  + atopobium-null       : tested, not significant (p = 0.674)
  + clostridia-null      : tested, not significant (p = 0.079)
  + bacteroidia-dominant : reanalysis top class Bacteroidia (50.7%)

Not reproduced (0)
Not licensed  (0)
Outside the denominator (0)

Scorable: 8 of 8.  Reproduced: 8 of 8.
```

**8 of 8 claims reproduced, through a different taxonomic profiler.**

Under the significance-first rule the five directional and null claims were
rechecked against their reanalysis p-values, and all five hold: the three alpha
claims and both null taxon claims are non-significant as the paper asserted, and
both positive taxon claims are significant at p = 0.017 and p = 0.004.

### What this means

For the paper: its conclusions are robust to the choice of profiler. That is a
stronger result than the paper itself could claim, since it compared 16S against
only one shotgun pipeline.

For the instrument: **a clean study produces a clean verdict.** This is the
negative-control function the validation matrix was designed around.
Pérez-Losada returned a claim the evidence could not license; Vogtmann returns
everything holding. The tool is not manufacturing disagreement.

### Limitations

1. **Different profiler.** MetaPhlAn, not the paper's MOCAT/mOTU. Agreement is
   evidence that conclusions survive a profiler change, not that the paper's
   pipeline was reproduced.
2. **Reconstructed counts.** `counts = TRUE` multiplies relative abundances by
   read depth. The values are integers but carry no true sampling variance, so
   count-based methods assuming a multinomial process are on weaker ground here
   than with 16S data.
3. **Gene, module and pathway claims not tested.** The paper's functional
   findings require the HUMAnN products, which were not downloaded.
4. **No sweep was run.** This is a single reconstruction, not a sensitivity
   analysis. Given that Pérez-Losada's sweep materially changed its
   interpretation, a sweep here would strengthen the reproduced verdict rather
   than merely confirm it.
5. **Reference phylogeny.** curatedMetagenomicData ships a tree of 10,430 tips
   for 540 taxa, and 8 rows were dropped for having no tree match. That tree is a
   reference topology, not one estimated from these samples. MicroFitGut's tree
   guard asks "was this estimated from the sequence data?", which is the right
   question for 16S and unanswerable for shotgun species calls. No phylogenetic
   metric was computed here, so nothing depends on it, but the guard's semantics
   need revisiting before shotgun UniFrac is attempted.

---

## 5. Defects this run exposed

Two, both fixed.

### #9: unassigned samples were counted as a group

`validate_inputs()` computed group balance with `table(g, useNA = "ifany")`, so
the six samples with no `study_condition` were treated as a third level:

> Group sizes are unbalanced **8.7:1**, which is 52 divided by 6.

The design is **52 against 52, perfectly balanced**. Two harms followed. A
fabricated imbalance warning would push an analyst toward an unnecessary
correction, and, worse, the message that mattered was never issued at all.
Samples with no group assignment are dropped silently by almost every test,
changing *n* and the multiple-testing denominator without any record.

Balance and minimum-size checks now run over observed levels only, and missing
assignments raise their own warning naming the count and the consequence.

### #10: no way to encode "tested and not significant"

The claim schema had `da_direction` for "this taxon differs" and `da_count` for
"N taxa differ", but no type for **"this taxon was tested and did not differ"**.
That is the commonest claim in a reproducibility study and two of this paper's
eight. They could only be recorded as unscored, which correctly blocked a clean
verdict but did so for the wrong reason: the evidence to settle them was present.

Added `da_null`, which distinguishes three states that matter:

| Reanalysis result | Status |
|---|---|
| taxon tested, not significant | `reproduced` |
| taxon tested, significant | `not_reproduced` |
| taxon absent from the result table | `not_tested`: absent is not the same as tested and null |

**Disclosure.** The verdict changed as a result of this addition. The run log
holds two verdict runs forty seconds apart, `held=6` then `held=8`, either side
of the change. The addition encodes a kind of claim the schema could not
previously express, which is a defensible reason, but the sequence was score,
inspect the result, change the scorer, rescore. It is recorded here rather than
left in the logs.

---

## 6. Reproducing this

```r
library(curatedMetagenomicData)
se <- curatedMetagenomicData("2021-03-31.VogtmannE_2016.relative_abundance",
                             dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]

# Taxonomy comes from rowData, NOT from parsing rownames: "short" rownames are
# species names, not lineages.
rd  <- as.data.frame(SummarizedExperiment::rowData(se))
tax <- as.matrix(rd[, c("superkingdom","phylum","class","order","family","genus","species")])
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

ps <- phyloseq(otu_table(SummarizedExperiment::assay(se), taxa_are_rows = TRUE),
               tax_table(tax),
               sample_data(as.data.frame(SummarizedExperiment::colData(se))))
ps <- prune_samples(!is.na(mfg_meta(ps)$study_condition), ps)   # drops 6
ps <- mfg_set_normalization(ps, "raw")
validate_inputs(ps, group_var = "study_condition", subject_var = "subject_id")
```

Seeds fixed at 42. Every value above is written to `run_log.csv` in the run
directory.
