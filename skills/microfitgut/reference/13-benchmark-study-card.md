# Benchmark: the study card

Benchmark mode reanalyses data whose published result is known, then explains any
difference. The study card is the first step and the one that most determines
whether the rest is meaningful.

## What it is

A structured record of **what the paper says, what it does not say, and where
each value came from.** Built with `study_card()` and `study_field()` in
`scripts/14-benchmark-sweep.R`.

```r
card <- study_card(
  study_id = "Author-Year-shortname",
  group_variable = study_field("region", "stated",
                               "samples from three body regions", "Methods 2.1"),
  normalization  = study_field("tss", "inferred",
                               "relative abundances were compared", "Methods 2.4"),
  da_method      = study_field(confidence = "absent"),
  reported_da_taxa = study_field(c("Staphylococcus", "Corynebacterium"), "stated",
                                 "two genera differed by region", "Results 3.2"))
```

## The three confidence levels do the real work

| confidence | Meaning | Consequence |
|---|---|---|
| `stated` | The paper says it. | **Fixed.** Not swept. |
| `inferred` | Deduced from something else the paper says. | **Swept**, with the inference recorded. |
| `absent` | The paper does not say. | **Swept** across all schema options. |

**"Not stated" is a first-class value and a finding in its own right.** A paper
that never says which SILVA release it used generates several replication
candidates instead of one, and that unstated parameter belongs in the benchmark
report as a reporting finding, separate from any replication failure.

Anything the schema knows about and the card omits defaults to `absent` — so
forgetting a field makes it a sweep axis rather than silently fixing it at a
default.

## Evidence and location

`evidence` should be a **verbatim span of at most ~25 words**, and `location` a
pointer such as `"Methods 2.3"` or `"Supp. Table S1"`, so the extraction itself
can be checked. A `stated` field without evidence cannot be verified and
`study_field()` warns.

This matters because study-card extraction is where a benchmark system fails
quietly: a confidently wrong extraction produces a confidently wrong attribution
downstream, and nothing in the pipeline catches it. The `verifier` subagent checks
the card against the paper before the sweep runs.

## Human approval before anything runs

The card is a **checkpoint**. Present it to the user and get confirmation before
executing the sweep. The extraction is the one step that cannot be validated
against an artifact — everything downstream can.

## The fields

### Identification (not swept)
`accessions`, `n_samples`, `design`, `group_variable`, `subject_variable`

`subject_variable` is load-bearing: if the paper has repeated measures and did
not model them, that is a likely divergence cause (`15`,
`statistics.repeated_measures`).

### Upstream (recorded, not swept downstream)
`amplicon_region`, `denoiser`, `reference_db`

Not swept because MicroFitGut starts from a feature table — it does not re-run
denoising. But recorded, because they are the commonest cause of taxon-name
mismatches, and `reference_db` in particular predicts the
`taxonomy.database_version` signature.

### Downstream (every one a sweep axis when unstated)
`prevalence_filter`, `min_depth`, `normalization`, `rarefaction_depth`,
`taxonomic_level`, `alpha_metric`, `beta_distance`, `da_method`, `fdr_method`,
`alpha_threshold`, `covariates`, `random_effects`

### The ground truth
`reported_alpha`, `reported_permanova`, `reported_da_taxa`, `reported_claims`

**This is what the benchmark scores against**, so extract it precisely:

- `reported_alpha` — the means or medians per group, the test, the p-value.
- `reported_permanova` — R² **and** p. R² alone is not enough; the delta in R² is
  a headline concordance metric (`14`).
- `reported_da_taxa` — the named taxa with their effect sizes and direction.
  Names exactly as the paper writes them; `harmonize_taxon()` handles the
  renaming, and it can only do so if the original label is preserved.
- `reported_claims` — the sentences the paper actually concludes with. A benchmark
  scores whether the *claim* reproduced, not only whether the numbers matched.

## From card to sweep

```r
sg <- build_sweep_grid(card, ps, max_runs = 20)
```

Baseline is the best reconstruction of the paper. Each further run changes exactly
one axis, so its effect is isolated.

Axes are ordered by expected impact — `da_method`, `normalization`,
`random_effects` first — and levels are scheduled **round-robin**, one per axis
per pass. That way a run cap thins every axis evenly instead of eliminating the
most important one because it came later in the schema.

Every axis gets at least one run. If `max_runs` is too small for that, the
function warns and names the axes that got none, because an axis with no run
cannot be attributed or ruled out.

Untested levels are enumerated explicitly in `sg$untested` and must appear in the
benchmark report. A silent cap reads as "we covered everything" when it did not.

## Why one-at-a-time rather than factorial

Full factorial over five axes is 200+ runs. One-at-a-time from the baseline is
12–15 and enough for **marginal** attribution, which is what the question needs.

Add a small factorial only for axes known to interact. In practice
denoiser × DA method is the one that genuinely does — and since denoiser is not a
downstream axis here, interactions are rarely worth the run budget.
