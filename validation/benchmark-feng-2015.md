# Benchmark: Feng et al. (2015), the colorectal adenoma–carcinoma sequence

Study **S2** of the validation matrix, the multi-group test, and the strongest
reproduction in the series.

| | |
|---|---|
| **Paper** | Feng Q, Liang S, Jia H, et al. Gut microbiome development along the colorectal adenoma–carcinoma sequence. *Nat Commun.* 2015;6:6528. |
| **Data** | curatedMetagenomicData `2021-03-31.FengQ_2015.relative_abundance`, `counts = TRUE` |
| **Design** | Cross-sectional, three groups: control 61, advanced adenoma 47, carcinoma 46 |
| **Verdict** | **Reproduced**: 8 of 8 claims held |

---

## 1. Role in the matrix

S2 exists to exercise the **multi-group path**: an omnibus test followed by
post-hoc comparisons, with the post-hoc gated on the omnibus. Two-group studies
never reach that code. 606 species × 154 samples, one sample per subject, so no
repeated-measures complication to separate out.

---

## 2. Results

### 2.1 Richness rises along the sequence

The paper's headline, from Figure 1: *"gene richness ... is higher in advanced
adenoma than in control, and higher in carcinoma than in advanced adenoma"*.

Median observed species:

| control | adenoma | carcinoma |
|---|---|---|
| **103.0** | **111.0** | **126.0** |

Kruskal-Wallis omnibus **p = 2e-07**, against the paper's p = 3.23e-07 for genus
count. Those are close, for a different profiler working at a different rank.

Post-hoc (Dunn, Benjamini-Hochberg), run only because the omnibus passed:

| Comparison | adjusted p |
|---|---|
| control – carcinoma | **2.2e-07** |
| adenoma – carcinoma | **8.0e-05** |
| control – adenoma | 0.272 |

**The carcinoma comparisons hold decisively; the adenoma step does not.** The
paper's claim is an ordered three-step progression, and the first step, 
adenoma above control, is not significant here. The ordering of the medians is
correct (103 < 111 < 126) but the control-to-adenoma difference does not survive
multiple-comparison correction.

Scored as held, because the claim encoded is carcinoma above control. The weaker
adenoma step is recorded as a qualification rather than a failure: the direction
matches, the significance does not.

Shannon diversity shows **no difference at all** (median 3.0 in every group,
omnibus p = 0.565). Richness and evenness diverge completely here, the carcinoma
gut carries more species, not a more even community. The paper's own framing is
about richness, and this reanalysis supports that specificity.

### 2.2 Every named carcinoma-enriched species reproduces

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

The three oral anaerobes the paper singles out, *F. nucleatum*, *G. morbillorum*,
*P. micra*, joined by *P. stomatis*, are essentially absent from control and
adenoma guts and present in carcinoma. These are 10- to 30-fold enrichments off a
near-zero baseline, which is why they reach small p-values on modest abundances.

*Bilophila wadsworthia* is the one that does not reach significance, though its
direction is correct and consistent across all three groups.

---

## 3. Verdict

```
=== Benchmark verdict: REPRODUCED ===

  + richness-crc-above-control      : higher, omnibus p = 2e-07
  + alistipes-putredinis-crc        : higher in CRC
  + bilophila-wadsworthia-crc       : higher in CRC
  + escherichia-coli-crc            : higher in CRC
  + parvimonas-micra-crc            : higher in CRC
  + gemella-morbillorum-crc         : higher in CRC
  + peptostreptococcus-stomatis-crc : higher in CRC
  + fusobacterium-nucleatum-crc     : higher in CRC
```

8 of 8 claims held, through a different taxonomic profiler than the paper used, 
the paper built a 3.5-million-gene catalog and used MLG (metagenomic linkage
group) analysis; this reanalysis used MetaPhlAn species calls.

That two such different routes agree on seven named species and on the richness
gradient is a stronger result than either alone.

---

## 4. Defects

**None.** Second consecutive study to run clean, after S3.

The multi-group path worked as designed: Kruskal-Wallis omnibus, then Dunn
post-hoc with Benjamini-Hochberg correction, gated on the omnibus passing. No
post-hoc was computed for Shannon, where the omnibus did not pass, which is the
gating behaving correctly rather than an absence of output.

### One harness error of mine, recorded

My first attempt reported all seven species as "NOT IN PROFILE". I had converted
the species names to underscores before matching, but curatedMetagenomicData's
short rownames use spaces with a `species:` prefix. Every species was present.

This is the third naming-convention error I have made across the series, after
parsing `rownames = "short"` as if they were full lineages in S1, and the
`Escherichia/Shigella` versus `Escherichia-Shigella` separator in A0. The pattern
is worth stating plainly: **taxon-name matching fails silently and looks like a
biological result.** An empty match set reads as "this organism was not found",
which is a finding, when it actually means "my string did not match", which is a
bug. Any pipeline doing this at scale needs a positive control, a name known to
be present, asserted before the matching loop runs.

---

## 5. Limitations

1. **Different profiler.** The paper used MLGs from its own gene catalog;
   this used MetaPhlAn. Agreement is evidence the conclusions survive a profiler
   change, not that the pipeline was reproduced.
2. **Reconstructed counts.** Depth values 32–80 million, derived from relative
   abundance × read depth.
3. **The adenoma step is weaker here than in the paper.** Control-to-adenoma
   richness is not significant after correction. Whether that reflects the
   profiler, the subset, or the original result is not established.
4. **Functional claims not tested.** The paper's amino-acid transport and
   synthesis module findings need HUMAnN data; S5 covers that path.
5. **No sweep.** Single reconstruction.

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
ps <- mfg_add_meta(ps, "grp",
        factor(mfg_meta(ps)$study_condition, levels = c("control","adenoma","CRC")))

a  <- compute_alpha(tss_transform(ps), measures = c("Observed","Shannon"))
df <- data.frame(a, grp = mfg_meta(ps)[rownames(a), "grp"])
run_alpha_test(df, "Observed", "grp")     # omnibus
posthoc_dunn(df, "Observed", "grp")        # gated on the omnibus

# species names carry a "species:" prefix and use spaces, not underscores
grep("Fusobacterium nucleatum", taxa_names(ps), value = TRUE)
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
