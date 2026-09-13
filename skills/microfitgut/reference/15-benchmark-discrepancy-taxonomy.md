# Benchmark: discrepancy taxonomy

Known divergence causes paired with **what each looks like in the metrics**. This
is how the agent goes from observing a difference to naming one, instead of
reasoning it out fresh each time and producing something plausible.

The machine-readable form is `MFG_DISCREPANCY_RUBRIC` in
`scripts/14-benchmark-sweep.R` (15 entries); `discrepancy_rubric(search = "...")`
retrieves them.

## How to use it

1. Run the sweep. `attribute_divergence(results)` ranks axes by the **spread**
   each induces in each concordance metric. This is computed, not asserted.
2. Read the ranking against the table below. The dominant axis points at a stage;
   the signature confirms or refutes it.
3. Check the **confidence** column of the attribution. A ratio below 1.3 between
   the top two axes means the sweep does not distinguish them, and **no single
   cause should be named**.
4. Record it with `attribution_record()`, which requires an
   `unexplained_residual`.

## The rubric

### Sequencing / upstream

| id | Cause | Signature |
|---|---|---|
| `upstream.sample_selection` | Different runs or samples included | Sample counts differ before any filtering; the shared-sample intersection is much smaller than either table. Profile concordance can be high on shared samples while group results differ. |
| `upstream.pretrimmed_reads` | Reads already trimmed differently | Depths differ by a roughly constant factor; richness shifts in one direction across all samples. |

### Denoising

| id | Cause | Signature |
|---|---|---|
| `denoise.otu_vs_asv` | 97% OTUs versus exact ASVs | Feature counts differ several-fold. Richness and Chao1 differ strongly, Shannon much less, and **beta-diversity structure is largely preserved** (high Mantel rho) because the abundant community is the same. |

### Taxonomy

| id | Cause | Signature |
|---|---|---|
| `taxonomy.database_version` | Different reference release | **The signature cause of "disappearing" taxa.** Many published names fail to match, and the unmatched list is dominated by renamed clades. Recovery rate is low while effect-size correlation on the matched subset stays high. Confirm by checking whether the unmatched names appear in `MFG_TAXON_SYNONYMS`. |
| `taxonomy.confidence_threshold` | Different classifier confidence | The share of reads unassigned at genus/species differs markedly; agglomerating loses a different fraction on each side. Check `reads_unresolved_pct`. |

### Normalization

| id | Cause | Signature |
|---|---|---|
| `normalization.rarefaction_depth` | Different depth, or rarefaction versus none | Richness and Chao1 differ with a **systematic Bland-Altman bias whose slope scales with depth**, while Shannon and beta barely move. Sample counts may differ if one depth dropped samples. |
| `normalization.method` | Compositional versus count-based | DA lists diverge substantially while alpha and beta agree. Sign agreement stays high; magnitudes are on different scales. |

### Filtering

| id | Cause | Signature |
|---|---|---|
| `filtering.prevalence` | Different prevalence or abundance filter | The multiple-testing denominator differs, so borderline taxa move in and out of significance. **Missed taxa turn out never to have been tested** — check `missed_because_not_tested`. |

### Statistics

| id | Cause | Signature |
|---|---|---|
| `statistics.da_method` | Different DA method | **Moderate recovery with high effect-size rank correlation.** The methods agree about direction and size and disagree about where to draw the line. LEfSe and plain Wilcoxon call more taxa than ANCOM-BC2 on identical counts. |
| `statistics.covariates` | Covariates in one analysis only | PERMANOVA R² for the main term differs substantially while the distance matrices agree (high Mantel rho). Some DA hits vanish when a covariate absorbs their variance. |
| `statistics.repeated_measures` | Repeated measures treated as independent | The independent-samples analysis reports **many more** significant results. The p-value histogram is anti-conservative. Adding a subject random effect removes most of the difference. |
| `statistics.fdr_method` | Different correction | Raw p-values agree closely; adjusted values differ. Holm is far more conservative than BH, so hit counts can differ severalfold with identical underlying statistics. |

### Metadata

| id | Cause | Signature |
|---|---|---|
| `metadata.recoded_groups` | Groups recoded, merged, or reference level changed | **Effect signs flip wholesale while magnitudes are preserved.** `direction_flip` fires while correlation stays high. |
| `metadata.outlier_exclusion` | Unreported outlier removal | Sample counts differ by a handful with no stated rule. Removing the same few samples reconciles the results. |

### Stochasticity

| id | Cause | Signature |
|---|---|---|
| `stochastic.unseeded` | Unseeded subsampling or too few permutations | Results differ slightly and **irreproducibly between runs of the same configuration**. p-values near the threshold flip; effect sizes barely move. |

**Check this one first when the divergence is small.** Re-run the identical
baseline twice. If the two baselines differ, the divergence is noise and no axis
will explain it — and a sweep run against a noisy baseline will attribute that
noise to whichever axis happens to spread most.

## The unexplained residual

`attribution_record()` **refuses** a non-reproducing verdict without an
`unexplained_residual`. That field is the escape valve that stops a plausible
cause being invented when none of the tested axes explain the gap.

> "No tested factor accounts for this" is a valid and valuable answer.

Prompt for it explicitly, and be specific about what remains unaccounted for:

> "Jaccard stays below 0.62 in every configuration, so the reanalysis calls taxa
> the paper did not report in all cases. No tested axis accounts for that excess."

## Coverage limits

An axis with no run cannot be ruled out. `build_sweep_grid()` enumerates untested
levels in `$untested` and names axes that got no run at all. Both must appear in
the benchmark report — this is the same "no silent caps" rule as elsewhere.

## Every cause needs three things

`attribution_record()` enforces the shape:

```r
list(rubric_id = "statistics.da_method",
     axis = "da_method",
     evidence = "recovery 0.33 (aldex2) to 1.00 (deseq2); spread 0.67, ratio 2.4 over the next axis",
     explanation = "...")
```

- **`rubric_id`** — which known cause this is, so the claim is anchored.
- **`axis`** — which sweep axis demonstrated it.
- **`evidence`** — the run ids and the metric movement. Not "the DA method
  differed" but the numbers that show it.

A cause without a run id behind it is a story. Under the traceability rule in
`10`, it does not get reported.
