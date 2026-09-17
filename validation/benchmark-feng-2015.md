# Benchmark: Feng et al. (2015), the colorectal adenoma-carcinoma sequence

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **S2** of the validation matrix, the
multi-group test, and the strongest reproduction in the series.

| | |
|---|---|
| **Paper** | Feng Q, Liang S, Jia H, et al. Gut microbiome development along the colorectal adenoma-carcinoma sequence. *Nat Commun.* 2015;6:6528. [10.1038/ncomms7528](https://doi.org/10.1038/ncomms7528) |
| **Raw data** | ENA [ERP008729](https://www.ebi.ac.uk/ena/browser/view/ERP008729) |
| **Analyzed** | `curatedMetagenomicData`, `2021-03-31.FengQ_2015.relative_abundance`, `counts = TRUE`, 606 species × 154 samples |
| **Design** | Cross-sectional, three groups: control 61, advanced adenoma 47, carcinoma 46. One sample per subject |
| **Hypothesis** | The gut microbiome changes progressively along the adenoma-carcinoma sequence, with richness rising and specific oral anaerobes enriching toward carcinoma |
| **Verdict** | **Partially reproduced.** 7 of 8 claims reproduced, 1 not reproduced, 0 unscorable. No sweep |

---

## 1. Study card

**10 fields stated · 2 inferred · 11 absent.**

S2 exists in the matrix to exercise the **multi-group path**: an omnibus test
followed by post-hoc comparisons, with the post-hoc gated on the omnibus.
Two-group studies never reach that code.

The absent fields are the usual downstream choices, plus one that matters more
here than elsewhere: the paper does not state a multiple-comparison method for
its per-species tests, so the reanalysis had to choose one.

### The paper's eight claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `richness-crc-above-control` | Richness rises along the sequence | "gene richness ... is higher in advanced adenoma than in control, and higher in carcinoma than in advanced adenoma", genus count p = 3.23e-07 |
| `fusobacterium-nucleatum-crc` | *F. nucleatum* enriched in carcinoma | named as carcinoma-enriched, Fig. 3 |
| `gemella-morbillorum-crc` | *G. morbillorum* enriched in carcinoma | named as carcinoma-enriched |
| `parvimonas-micra-crc` | *P. micra* enriched in carcinoma | named as carcinoma-enriched |
| `peptostreptococcus-stomatis-crc` | *P. stomatis* enriched in carcinoma | named as carcinoma-enriched |
| `alistipes-putredinis-crc` | *A. putredinis* enriched in carcinoma | named as carcinoma-enriched |
| `escherichia-coli-crc` | *E. coli* enriched in carcinoma | named as carcinoma-enriched |
| `bilophila-wadsworthia-crc` | *B. wadsworthia* enriched in carcinoma | named as carcinoma-enriched |

### One upstream difference, recorded rather than swept

The paper built a **3.5-million-gene catalog** and used MLG, metagenomic linkage
group, analysis. This reanalysis used **MetaPhlAn species calls**. The routes are
very different, which makes agreement more informative than it would be between
two runs of the same pipeline, and makes any disagreement harder to attribute.

---

## 2. Reconstruction

606 species × 154 samples on intake. Intake flagged:

- Depth 32 to 80 million reconstructed counts.
- **No repeated measures**: one sample per subject across 154 subjects, so no
  mixed model is needed anywhere below.
- Three groups, so the group factor was ordered deliberately as
  control, adenoma, CRC rather than left to alphabetical default. Alphabetical
  order would have placed CRC first and reversed every reported contrast.

---

## 3. Results

### 3.1 Richness rises along the sequence

Median observed species:

| control | adenoma | carcinoma |
|---|---|---|
| **103.0** | **111.0** | **126.0** |

Kruskal-Wallis omnibus **p = 2e-07**, against the paper's p = 3.23e-07 for genus
count. Those are close, for a different profiler working at a different rank.

Post-hoc (Dunn, Benjamini-Hochberg), run only because the omnibus passed:

| Comparison | adjusted p |
|---|---|
| control vs carcinoma | **2.2e-07** |
| adenoma vs carcinoma | **8.0e-05** |
| control vs adenoma | 0.272 |

**The carcinoma comparisons hold decisively. The adenoma step does not.** The
paper's claim is an ordered three-step progression, and the first step, adenoma
above control, is not significant here. The ordering of the medians is correct,
103 < 111 < 126, but the control-to-adenoma difference does not survive
multiple-comparison correction.

The claim as encoded is carcinoma above control, and that reproduces. The weaker
adenoma step is recorded as a qualification rather than a separate failure,
because it was not encoded as a separate claim. Encoding it separately would have
been the better choice and is noted as a study-card lesson.

Shannon diversity shows **no difference at all**, median 3.0 in every group,
omnibus p = 0.565. Richness and evenness diverge completely here: the carcinoma
gut carries more species, not a more even community. The paper's own framing is
about richness, and this reanalysis supports that specificity.

### 3.2 Six of seven named species reproduce

| Species | control | adenoma | carcinoma | Kruskal-Wallis |
|---|---|---|---|---|
| *Fusobacterium nucleatum* | 0.000% | 0.000% | 0.012% | **1.3e-07** |
| *Gemella morbillorum* | 0.001% | 0.001% | 0.020% | **1.1e-05** |
| *Parvimonas micra* | 0.001% | 0.001% | 0.034% | **3.8e-05** |
| *Alistipes putredinis* | 0.869% | 0.708% | 1.644% | **0.002** |
| *Escherichia coli* | 2.438% | 2.719% | 6.648% | **0.005** |
| *Peptostreptococcus stomatis* | 0.000% | 0.002% | 0.032% | **0.005** |
| *Bilophila wadsworthia* | 0.014% | 0.024% | 0.041% | 0.153 |

**All seven are highest in carcinoma. Six of seven are significant.**

The three oral anaerobes the paper singles out, *F. nucleatum*, *G. morbillorum*
and *P. micra*, joined by *P. stomatis*, are essentially absent from control and
adenoma guts and present in carcinoma. These are 10-fold to 30-fold enrichments
off a near-zero baseline, which is why they reach small p-values on modest
abundances.

*Bilophila wadsworthia* is the one that does not reach significance, at p = 0.153,
although its direction is correct and consistent across all three groups.

---

## 4. Verdict

```
=== Benchmark verdict: PARTIALLY REPRODUCED ===

Reproduced (7):
  + richness-crc-above-control      : higher in CRC, omnibus p = 2e-07
  + fusobacterium-nucleatum-crc     : higher in CRC, p = 1.3e-07
  + gemella-morbillorum-crc         : higher in CRC, p = 1.1e-05
  + parvimonas-micra-crc            : higher in CRC, p = 3.8e-05
  + alistipes-putredinis-crc        : higher in CRC, p = 0.002
  + escherichia-coli-crc            : higher in CRC, p = 0.005
  + peptostreptococcus-stomatis-crc : higher in CRC, p = 0.005

Not reproduced (1):
  - bilophila-wadsworthia-crc       : direction agrees, p = 0.153, not significant

Not licensed  (0)
Outside the denominator (0)

Scorable: 8 of 8.  Reproduced: 7 of 8.
```

**This verdict changed on rescoring.** The report originally read 8 of 8, because
*B. wadsworthia* was scored on the direction of the effect alone. Under the
significance-first rule a claim of difference requires the difference to be
detectable, and p = 0.153 does not meet it. The report's own body already said so:
"six of seven are significant".

### What this means

Seven of the paper's eight claims reproduce through a completely different
analytical route. The paper built a gene catalog and used MLG analysis; this used
MetaPhlAn species calls. That two such different routes agree on six named
species and on the richness gradient is a stronger result than either alone.

The one claim that does not reproduce is a failure to detect, not a
contradiction. *B. wadsworthia* is higher in carcinoma in this reanalysis too,
by roughly threefold off a small base. What is absent is evidence that the
difference is distinguishable from noise.

### Limitations

1. **Different profiler.** The paper used MLGs from its own gene catalog; this
   used MetaPhlAn. Agreement is evidence the conclusions survive a profiler
   change, not that the pipeline was reproduced.
2. **Reconstructed counts.** Depth values 32 to 80 million, derived from relative
   abundance times read depth.
3. **The adenoma step is weaker here than in the paper.** Control-to-adenoma
   richness is not significant after correction. Whether that reflects the
   profiler, the subset, or the original result is not established, and it was
   not encoded as its own claim, which it should have been.
4. **Functional claims not tested.** The paper's amino-acid transport and
   synthesis module findings need HUMAnN data.
5. **No sweep.** Single reconstruction. The *B. wadsworthia* result in particular
   would benefit from one, since a p-value of 0.153 is exactly the kind of value
   that moves under a different prevalence filter.

---

## 5. Defects this run exposed

**None in the shipped code.** Second consecutive study to run clean.

The multi-group path worked as designed: Kruskal-Wallis omnibus, then Dunn
post-hoc with Benjamini-Hochberg correction, gated on the omnibus passing. No
post-hoc was computed for Shannon, where the omnibus did not pass, which is the
gating behaving correctly rather than an absence of output.

### One harness error of mine, recorded

My first attempt reported all seven species as "NOT IN PROFILE". I had converted
the species names to underscores before matching, but curatedMetagenomicData's
short rownames use spaces with a `species:` prefix. Every species was present.

This is the third naming-convention error I made across the series, after parsing
`rownames = "short"` as if the values were full lineages in S1, and the
`Escherichia/Shigella` against `Escherichia-Shigella` separator in A0. The
pattern is worth stating plainly: **taxon-name matching fails silently and looks
like a biological result.** An empty match set reads as "this organism was not
found", which is a finding, when it actually means "my string did not match",
which is a bug. Any pipeline doing this at scale needs a positive control: a name
known to be present, asserted before the matching loop runs.

---

## 6. Reproducing this

```r
se <- curatedMetagenomicData("2021-03-31.FengQ_2015.relative_abundance",
                             dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]
rd  <- as.data.frame(SummarizedExperiment::rowData(se))
tax <- as.matrix(rd[, c("superkingdom","phylum","class","order","family","genus","species")])
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

ps <- phyloseq(otu_table(SummarizedExperiment::assay(se), taxa_are_rows = TRUE),
               tax_table(tax),
               sample_data(as.data.frame(SummarizedExperiment::colData(se))))
ps <- mfg_set_normalization(ps, "raw")
# Order the factor deliberately: alphabetical would put CRC first and reverse
# every contrast reported above.
ps <- mfg_add_meta(ps, "grp",
        factor(mfg_meta(ps)$study_condition, levels = c("control","adenoma","CRC")))

a  <- compute_alpha(tss_transform(ps), measures = c("Observed","Shannon"))
df <- data.frame(a, grp = mfg_meta(ps)[rownames(a), "grp"])
run_alpha_test(df, "Observed", "grp")     # omnibus
posthoc_dunn(df, "Observed", "grp")       # gated on the omnibus

# Species names carry a "species:" prefix and use spaces, not underscores.
grep("Fusobacterium nucleatum", taxa_names(ps), value = TRUE)
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
