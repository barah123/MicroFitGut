# QC and filtering

The governing rule: **no sample and no taxon leaves the dataset without a QC log
entry saying which rule dropped it and why.** `run_qc()` produces that record;
`qc_exclusions()` turns it into the table the report must contain.

## Order matters

```
1. unwanted taxa   (they are not data)
2. prevalence      (before depth, see below)
3. singletons
4. low-depth samples
5. empty taxa      (cleanup after sample removal)
```

Removing taxa changes every sample's depth. Pruning samples on a depth computed
*before* taxon removal drops the wrong samples. `run_qc()` enforces this order.

## 1. Unwanted taxa

| Target | Why it is not data |
|---|---|
| `Family == "mitochondria"` | Host mitochondrial 16S amplifies with bacterial primers. It is host DNA, not a community member. |
| `Order == "Chloroplast"` | Plant/algal plastid 16S. Same reason. Common in diet, soil and plant-associated samples. |
| `Kingdom` unassigned / NA / Eukaryota | Chimeras and off-target amplification. Cannot be interpreted, so they only add multiple-testing burden. |

Retaining these inflates richness and spends FDR budget on sequences no
conclusion can be drawn about. In the bundled skin dataset this removes 17 taxa
of 619 (7 mitochondrial families, 10 chloroplast orders).

For shotgun profiles this step is usually a no-op — profilers do not emit
organellar sequences — but host reads can survive as `Eukaryota`.

## 2. Prevalence filtering

**The filter that most changes a differential abundance result.**

```r
filter_prevalence(ps, min_prevalence = 0.10)   # fraction of samples
filter_prevalence(ps, min_prevalence = 5)      # absolute sample count
```

A taxon present in 3 of 300 samples cannot support a group-level inference, and
it costs a full multiple-testing penalty. Removing it *increases* power on
everything else by shrinking the FDR denominator.

**Choosing the threshold.** Plot it first:

```r
plot_prevalence(ps, rank = "Phylum", threshold = 0.10)
```

Look at the shape of the prevalence-versus-abundance cloud. Most microbiome
datasets show a dense low-prevalence, low-abundance cluster that is clearly
separable from the rest; the cut goes between them.

**The threshold is fixed before any test is run.** Choosing it after seeing the
p-values is the pitfall in `11` — it is a garden-of-forking-paths problem and
the resulting FDR does not mean what it says.

Common values: 0.10 is a reasonable default. 0.05 is permissive. 0.20 is
aggressive and appropriate when n is large. ANCOM-BC2 applies its own `prv_cut`
internally, so filtering twice compounds.

On the skin dataset, 0.10 removes 459 of 602 taxa — and retains **89.3% of the
reads**. That asymmetry is the whole point: the removed taxa were nearly all of
the multiple-testing burden and almost none of the data.

## 3. Singletons

```r
remove_singletons(ps, min_total_reads = 1)   # drops taxa with <= 1 read study-wide
```

A taxon observed once across the entire study cannot be distinguished from a
sequencing error. Some pipelines use a stricter cut (`> 2`) or a per-sample rule.

Skipped automatically on relative-abundance input: there are no read counts to
threshold, and applying it to proportions would remove taxa by an arbitrary
numeric coincidence.

## 4. Low-depth samples

```r
prune_low_depth(ps, min_depth = 5000,
                justification = "rarefaction curves plateau above ~5000 reads in this study")
```

`min_depth` has **no default on purpose**. There is no universal threshold, and
supplying one silently would make an analytical decision invisible. The function
errors without it and warns without a justification.

A shallow sample does not hold a small amount of the community — it holds a
*biased* sample of it, weighted toward abundant taxa. Its observed richness is an
artefact of sequencing effort.

**How to choose it.** Two pieces of evidence:

1. `choose_rarefaction_depth(ps)` — the sample-retention versus read-retention
   trade-off at candidate depths.
2. `mfg_rarefaction_curve(ps)` then `plot_rarefaction(rc, depth = d)` — does the
   curve flatten at your candidate depth? The function computes the slope over
   the final fifth of each curve against the initial fifth and reports how many
   samples plateau.

Values in the course material: 5,000 (Demo 8, Quiz 8) and 20,000 (Demo 6) for
different studies. Both are defended by their curves, not by convention.

## 5. Abundance filtering — optional

```r
filter_abundance(ps, min_total_reads = 10, min_mean_relative = 1e-4)
```

Less commonly needed than prevalence filtering, and more prone to removing real
low-abundance biology. Prefer prevalence.

## What must be logged

`qc_exclusions(log)` returns one row per dropped sample with the rule that
dropped it. `10-reporting-standards.md` requires this table verbatim in the
report — not a count, the actual list. A reader who cannot see which samples left
cannot judge whether the exclusions were selective.

`qc_table(log)` gives the step-by-step: taxa before/after, samples before/after,
percentage of reads retained, and the stated rationale for each step.

## What QC does not fix

- **Batch effects.** Filtering does not remove them. If samples were sequenced in
  runs, the run is a covariate and belongs in the model (`07`), not in the filter.
- **Contamination.** A negative-control-based approach (decontam and similar) is a
  different procedure from prevalence filtering and is not implemented here. If
  the study has blanks, say so and handle them explicitly.
- **Host reads in shotgun data.** Removed upstream, not here.
