# Benchmark: Vogtmann et al. (2016), colorectal cancer and shotgun metagenomics

The first **shotgun** benchmark, and the first study in the series to reproduce
in full. Study S1 of the validation matrix.

| | |
|---|---|
| **Paper** | Vogtmann E, Hua X, Zeller G, et al. Colorectal Cancer and the Human Gut Microbiome: Reproducibility with Whole-Genome Shotgun Sequencing. *PLOS ONE.* 2016;11(5):e0155362. |
| **Raw data** | ENA [PRJEB12449](https://www.ebi.ac.uk/ena/data/view/PRJEB12449) |
| **Analysed** | `curatedMetagenomicData` — `2021-03-31.VogtmannE_2016.relative_abundance`, `counts = TRUE` |
| **Design** | Case-control: 52 CRC cases, 52 controls, Washington DC 1985–87, frequency matched on sex and BMI |
| **Verdict** | **Reproduced** — 8 of 8 claims held |

---

## 1. Why this study is an unusual benchmark

Vogtmann et al. is itself a reproducibility study. It asks whether associations
found by 16S in a cohort are recovered when the same samples are sequenced by
whole-genome shotgun. Most of its claims are therefore **negative** — "we did not
reproduce X" — which makes it a sharper test of the benchmark machinery than a
conventional paper, because a negative claim can fail in both directions.

It is also the first study run through MicroFitGut's shotgun path. **No shotgun
data had ever been processed by the tool before this run.**

---

## 2. Study card

**12 fields stated, 0 inferred, 12 absent.**

Unusually complete for stated fields — this is a methodologically careful paper
with a detailed statistics section. The twelve absent fields are downstream
choices (prevalence filter, normalization for the comparison, FDR method,
beta-diversity distance) that the paper had no occasion to state because it did
not perform those analyses.

### Upstream difference, recorded not swept

| | Paper | Reanalysis |
|---|---|---|
| Profiler | **MOCAT / mOTU** | **MetaPhlAn** (via curatedMetagenomicData) |
| Alpha diversity input | mOTU relative abundance, downsampled to 2,000 inserts | MetaPhlAn relative abundance |
| Counts | true read counts | reconstructed (relative abundance × read depth) |

This is a genuine upstream divergence and the single most important qualification
on everything below. The reanalysis does **not** reproduce the paper's pipeline;
it reanalyses the same samples through a different taxonomic profiler. Per the
discrepancy rubric this is closest to `taxonomy.database_version`, whose predicted
signature is *detection rates differ while conclusions hold* — which is exactly
what was observed.

---

## 3. Reconstruction

540 species × 110 samples on intake. Six samples carry no `study_condition` and
were **excluded explicitly and logged** before analysis, leaving **52 CRC and 52
control** — matching the paper exactly.

Intake flagged:

- Depth varies **34.6-fold** (4.1M to 141.9M reconstructed counts).
- Sparsity **78.7%**.
- **No repeated measures** — one sample per subject, 110 subjects. Confirms S1's
  placement in the cross-sectional cell of the validation matrix.
- Character columns take an alphabetical reference level unless set deliberately.

---

## 4. Results

### 4.1 Alpha diversity — all three null results reproduce

The paper reports no case-control difference on any of three indices.

| Index | Paper (Kruskal-Wallis) | Reanalysis | Conclusion |
|---|---|---|---|
| Shannon | p = 0.252 | p = 0.899 | both null |
| Richness | p = 0.157 | **p = 0.143** | both null |
| Evenness | p = 0.438 | p = 0.082 | both null |

Richness agrees almost exactly. Evenness moves closer to the threshold but does
not cross it. The paper's observation that "in general the controls had slightly
higher alpha diversity" also holds for evenness (0.320 control vs 0.280 CRC),
though not for Shannon.

### 4.2 Presence/absence — every association reproduces, at half the detection

| Taxon | Paper: case % / control % / p | Reanalysis: case % / control % / p |
|---|---|---|
| Fusobacteria (phylum) | 76.9 / 48.1 / **0.003** | 40.4 / 17.3 / **0.017** |
| *Fusobacterium* (genus) | 75.0 / 48.1 / **0.006** | 40.4 / 17.3 / **0.017** |
| *Porphyromonas* | 61.5 / 40.4 / **0.032** | 28.8 / 5.8 / **0.004** |
| *Atopobium* | 53.8 / 44.2 / 0.328 (ns) | 7.7 / 3.8 / 0.674 (ns) |

Every conclusion holds — three significant associations in the same direction,
one null result still null.

**But detection rates are roughly halved.** Fusobacteria is detected in 76.9% of
cases by mOTU and 40.4% by MetaPhlAn; *Atopobium* falls from 53.8% to 7.7%.
MetaPhlAn requires marker-gene evidence before calling a species present, which
is a stricter criterion than mOTU's. The prevalence estimates are therefore not
interchangeable between profilers even though the case-control comparisons agree.

This matters beyond this study: **a paper reporting a prevalence figure is
reporting a property of its profiler as much as of its cohort.**

### 4.3 Relative abundance — both null results reproduce, closely

| Taxon | Paper (Wilcoxon) | Reanalysis |
|---|---|---|
| Clostridia (class) | 33.9% case / 39.0% control, p = 0.092 | 36.6% / 41.7%, **p = 0.079** |
| Bacteroidia — highest abundance class | 53.2% case / 50.9% control | **50.7%, top class** |

Clostridia reproduces to within 0.013 on the p-value and preserves the direction
the paper noted ("tended to be lower in cases"). Bacteroidia is confirmed as the
most abundant class.

---

## 5. Verdict

```
=== Benchmark verdict: REPRODUCED ===

Claims that held:
  + alpha-shannon-ns     : paper: none; reanalysis: none
  + alpha-richness-ns    : paper: none; reanalysis: none
  + alpha-evenness-ns    : paper: none; reanalysis: none
  + fusobacterium-cases  : paper says higher in CRC; reanalysis -> CRC
  + porphyromonas-cases  : paper says higher in CRC; reanalysis -> CRC
  + bacteroidia-dominant : reanalysis top class Bacteroidia (50.7%)
  + atopobium-null       : paper tested, not significant; reanalysis not significant
  + clostridia-null      : paper tested, not significant; reanalysis not significant
```

**8 of 8 claims held**, through a different taxonomic profiler.

### What this establishes

For the paper: its conclusions are robust to the choice of profiler, which is a
stronger result than the paper itself could claim, since it only compared 16S
against one shotgun pipeline.

For the instrument: **a clean study produces a clean verdict.** This is the
negative-control function the validation matrix was designed around. Pérez-Losada
returned *partially reproduced* with a claim the evidence could not license;
Vogtmann returns *reproduced* with everything holding. The tool is not
manufacturing disagreement.

---

## 6. Defects this run exposed

Two, both fixed.

### #9 — Unassigned samples were counted as a group

`validate_inputs()` computed group balance with `table(g, useNA = "ifany")`, so
the six samples with no `study_condition` were treated as a third level:

> Group sizes are unbalanced **8.7:1** — from 52 ÷ 6.

The design is **52 vs 52, perfectly balanced**. Two harms: a fabricated imbalance
warning that would push an analyst toward an unnecessary correction, and — worse —
the message that mattered was never issued at all. Samples with no group
assignment are dropped silently by almost every test, changing *n* and the
multiple-testing denominator without any record.

Balance and minimum-size checks now run over observed levels only, and missing
assignments raise their own warning naming the count and the consequence.

### #10 — No way to encode "tested and not significant"

The claim schema had `da_direction` ("this taxon differs") and `da_count`
("N taxa differ") but no type for **"this taxon was tested and did not differ"** —
the commonest claim in a reproducibility study and two of this paper's eight.
They could only be recorded as `unscored`, which correctly blocked a clean
verdict but did so for the wrong reason: the evidence to settle them was present.

Added `da_null`, which distinguishes three states that matter:

| Reanalysis result | Status |
|---|---|
| taxon tested, not significant | **held** |
| taxon tested, significant | **failed** |
| taxon absent from the result table | **not_evaluable** — absent is not the same as tested-and-null |

With `da_null` in place, both claims scored `held` and the verdict moved from
`incomplete` to `reproduced` — the correct answer, reached by encoding the
evidence rather than by relaxing a rule.

---

## 7. Limitations

1. **Different profiler.** MetaPhlAn, not the paper's MOCAT/mOTU. Agreement is
   therefore evidence that conclusions survive a profiler change, not that the
   paper's pipeline was reproduced.
2. **Reconstructed counts.** `counts = TRUE` multiplies relative abundances by
   read depth. The values are integers but carry no true sampling variance, so
   count-based methods assuming a multinomial process are on weaker ground than
   with 16S data.
3. **Gene, module and pathway claims not tested.** The paper's functional
   findings (4 of 10 genes, 3 of 9 modules, 7 of 17 pathways reproduced from
   population F) require the HUMAnN products, which were not downloaded. Study
   **S5** covers the functional path.
4. **No sweep was run.** This is a single reconstruction, not a sensitivity
   analysis. Given that Pérez-Losada's sweep materially changed its
   interpretation, a sweep here would strengthen the reproduced verdict.
5. **Reference phylogeny.** curatedMetagenomicData ships a tree of 10,430 tips
   for 540 taxa, and 8 rows were dropped for having no tree match. That tree is a
   reference topology, not one estimated from these samples — MicroFitGut's tree
   guard asks "was this estimated from the sequence data?", which is the right
   question for 16S and unanswerable for shotgun species calls. No phylogenetic
   metric was computed here, so nothing depends on it, but the guard's semantics
   need revisiting before shotgun UniFrac is attempted.

---

## 8. Reproducing this

```r
library(curatedMetagenomicData)
se <- curatedMetagenomicData("2021-03-31.VogtmannE_2016.relative_abundance",
                             dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]

# taxonomy comes from rowData, NOT from parsing rownames — "short" rownames are
# species names, not lineages
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
