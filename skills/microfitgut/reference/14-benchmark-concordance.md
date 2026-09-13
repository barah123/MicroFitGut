# Benchmark: concordance metrics

## Do not reduce agreement to one number

A single "concordance score" hides the thing you want to know: **where** the
agreement breaks.

Alpha diversity can correlate at r = 0.98 per sample and still flip the
group-level significance. A DA list can recover 40% of the named taxa while the
effect sizes rank-correlate at 0.9. Those are different findings and they license
different conclusions.

So each domain gets its own metrics, and the verdict is driven by **conclusion
flips**, not by correlation magnitudes — because the question a benchmark answers
is "does the paper's claim hold", not "how similar are the numbers".

## Harmonise before comparing

`12-benchmark-harmonize.R`, and it is not optional. Skipping it **manufactures**
disagreement: a taxon that "disappeared" is most often the same organism renamed
(`11`, item 9). Counting those as misses inflates the divergence, and then the
attribution step invents a cause for it.

```r
tm <- match_taxa(published_taxa, reanalysis_taxa, published_tax, reanalysis_tax)
```

Three matching rules, tried in order, each recorded per pair:

| Rule | Confidence |
|---|---|
| exact match on harmonised names | high |
| genus-level truncation (paper gives species, reanalysis resolves to genus) | medium |
| via the taxonomy table (reanalysis uses ASV ids) | high at Species/Genus, medium above |

"Matched after genus truncation" is a weaker claim than "matched exactly", and the
report should say which it relied on.

**Unresolvable labels are reported separately.** A published label of
`"unclassified"` cannot be matched to anything, and that is a finding about the
paper's reporting — not a replication failure. Keep the two categories apart.

`align_abundance_tables()` also reports **abundance coverage**: what share of each
table's abundance the matched taxa account for. A concordance computed on the 12
taxa two tables happen to share is not a concordance between the tables, and the
function warns below 70%.

## Alpha diversity — `concordance_alpha()`

Three separate questions:

1. **Do the per-sample values agree?** Pearson and Spearman.
2. **Is there a systematic offset?** Bland-Altman bias and limits of agreement.
   *A high correlation with a constant offset means the values are not
   interchangeable even though they track perfectly* — the signature of a
   different rarefaction depth or normalization. The proportional-bias slope says
   whether the difference grows with magnitude.
3. **Does the group-level conclusion survive?** Significance flip and direction
   flip.

The third is what determines whether the paper's claim reproduced, and a
correlation coefficient cannot answer it.

A **direction flip** with high correlation almost always means a recoded group or
a changed reference level (`11`, item 13) — check that before attributing it to a
pipeline difference.

## Beta diversity — `concordance_beta()`

| Metric | What it compares |
|---|---|
| **Mantel rho** | the two distance matrices directly, independent of ordination |
| **Procrustes M²** | the two ordination configurations after optimal rotation, scaling and reflection |
| **protest p** | whether the configurations are more similar than chance |
| **delta R²** | the PERMANOVA effect size |
| significance flip | the PERMANOVA conclusion |

Procrustes is the right comparison for ordinations because an ordination has **no
intrinsic axis orientation** — comparing coordinates directly would report a
difference where there is none. Lower M² is better.

The informative split: if Mantel rho is high and the PERMANOVA flips, the
distances agree and the *test* diverged (different permutation scheme, different
covariates, stratification). If Mantel rho is low, the divergence is upstream of
the ordination.

## Differential abundance — `concordance_da()`

Two numbers that are routinely conflated, and both are reported:

- **Recovery rate** — what fraction of the paper's named taxa the reanalysis found.
- **Jaccard** — set overlap over the union of both lists.

They differ whenever the lists are different sizes. A reanalysis finding 200 taxa
including all 20 of the paper's has **perfect recovery and a Jaccard of 0.10**.
Reporting only one misrepresents the result.

Also computed:

- **Effect-size Spearman** over taxa present in both result tables *regardless of
  significance*. This is the more robust comparison, because significance is a
  thresholded quantity and near-threshold taxa flip easily.
- **Sign agreement.** A taxon enriched in the opposite group is a contradiction,
  not a smaller effect, and matters more than magnitude.
- **`missed_because_not_tested`** — a taxon the reanalysis filtered out before
  testing is a *filtering* difference, not a disagreement about significance.
  Separating these is what makes the attribution correct.

**The diagnostic pattern to recognise:** high effect-size correlation with low
recovery means the methods agree about direction and magnitude and disagree about
where to draw the line. That is a threshold or power difference, not a different
biological answer.

## Abundance profiles — `concordance_profile()`

Only possible when the paper deposited its table — which is the case worth
**starting** a benchmark with, because it isolates the comparison engine from the
execution layer.

- **Bray-Curtis per sample** — differences in the abundant community.
- **Aitchison per sample** — the compositionally coherent counterpart, sensitive
  to the rare tail.

Report both: they disagree informatively. A difference confined to abundant taxa
shows in Bray-Curtis; one confined to rare taxa shows in Aitchison.

- **Per-taxon Spearman** identifies which taxa the two tables disagree most about,
  which often points straight at the cause.

If mean Bray-Curtis is below ~0.1, the profiles are near-identical and any
downstream disagreement comes from the statistical layer. Above ~0.3, the
reanalysis is not reconstructing the published table, and **attributing downstream
divergence to the statistical method would be wrong** — the cause is upstream.

## The verdict — `concordance_verdict()`

Categorical, not a score:

| Verdict | Condition |
|---|---|
| `reproduced` | no conclusion flipped |
| `partially_reproduced` | some claims held, others did not |
| `diverged` | the central claims did not hold |

It also emits **qualifications** that are not flips but change the reading — for
example an R² that changed by more than 50% while staying significant. The
conclusion reproduced; the effect size did not. Both belong in the report.
