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

## Scoring claims — `score_claims()`

The concordance metrics above describe *how close* two analyses are. They do not
say whether the paper's claim survived, and they must not be allowed to decide
that. A recovery rate of 0.68 against a threshold of 0.7 is not a contradicted
finding; it is a number near a line someone drew.

So scoring runs on **claims**, one at a time, and every claim ends in exactly one
status.

### Four scorable statuses

These, and only these, form the denominator of any reproduction rate.

| Status | Meaning |
|---|---|
| `reproduced` | the reanalysis supports the claim as stated |
| `not_reproduced` | the reanalysis contradicts it |
| `not_licensed` | the test ran, but its assumptions do not license the claim |
| `not_tested` | the quantity was filtered out before testing |

`not_licensed` is the one most easily lost. A significant PERMANOVA under
heterogeneous dispersion is a real result about *something*, and it is not a
composition claim. Pooling that with a missing metadata column throws away a
finding about the paper.

### Seven unscorable reasons

These are reported with their own counts and are **never pooled**, because they
mean different things about different objects.

| Reason | What it is about |
|---|---|
| `schema_gap` | the instrument: no claim type expresses this |
| `out_of_scope` | the paper: not a microbiome measurement |
| `data_absent` | data sharing: the needed data was never deposited |
| `contrast_mismatch` | the deposit: the paper's contrast cannot be built from it |
| `not_attempted` | this run: the analysis was not done |
| `pending_adjudication` | this run: awaiting a decision |

A study whose claims are mostly `data_absent` and one whose claims are mostly
`not_reproduced` are opposite findings. A single "unadjudicated" bucket makes
them look alike.

### The directional rule

A directional claim is scored **significance first, then direction**. A
non-significant result in the claimed direction is `not_reproduced`, because the
paper asserted a difference and the reanalysis found none. Scoring on sign alone
scores noise: an effect at p = 0.955 points somewhere, and that direction means
nothing.

Both sides are symmetric. A taxon absent from the reanalysis is `not_tested`
whether the claim was positive or null, and `da_count` comparisons require an
explicit `tolerance` whenever the comparator is approximate, since an unstated
tolerance is a threshold chosen after seeing the number.

## The endpoint — `concordance_verdict()`

```r
verdict <- concordance_verdict(claims = scored, primary_claim_id = "<id>")
```

`claims` is required. The domain concordance objects are accepted and **ignored**
with a warning; they are descriptive, and a benchmark endpoint that moves when a
recovery threshold moves is measuring the threshold.

The outcome is **one pre-designated claim**, named before the reanalysis runs.
The verdict object also carries:

- `proportion_reproduced` over `n_scorable`, a secondary summary
- `status_counts`, every status reported separately
- `category`, one of `reproduced`, `partially_reproduced`, `diverged`

The category is flagged `category_is_descriptive_only = TRUE` and is **not
ordinal**. Three claims held of three is not twice as good as three of six when
the other three were never deposited.

It also emits **qualifications** that are not flips but change the reading, for
example an R² that changed by more than half while staying significant. The
conclusion reproduced; the effect size did not. Both belong in the report.

An unhandled claim type is an error, not a skipped claim. A scorer that fails
open turns a coverage gap into a silent pass.
