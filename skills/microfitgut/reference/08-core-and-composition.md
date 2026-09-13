# Core microbiome and composition

## "Core microbiome" means nothing without two thresholds

A core set is defined by:

- **detection** — how abundant before a taxon counts as present, on the
  relative-abundance scale. `0.2` means 20% of a sample.
- **prevalence** — in what fraction of samples it must clear detection.

Both must be quoted with any core result. `core_taxa()` returns them in the
object and `threshold_sentence` writes the sentence for the report.

Demo 6 uses `detection = 0.2, prevalence = 0.9`, which is strict. Much of the
literature uses detection near zero with prevalence 0.5, which finds a far larger
core. Both are defensible; neither is a default.

## Show the surface, not one number

```r
core_threshold_scan(ps, rank = "Genus")
```

On `ps9.RDS` at genus level:

| detection | prevalence | n core |
|---|---|---|
| 0 | 0.5 | 11 |
| 0 | 0.9 | 2 |
| 1e-04 | 0.5 | 9 |
| 1e-03 | 0.5 | 5 |
| 1e-03 | 0.9 | 2 |
| 1e-02 | 0.9 | 1 |

Quoting "the core comprised 11 genera" or "the core comprised 2 genera" are both
true and they are the same data. The scan is the honest presentation, and quoting
a single cut without it is a choice the reader cannot see.

## What share of the community is the core?

A core of 5 taxa holding 78% of every sample is a different finding from 5 taxa
holding 5%. `core_taxa()` reports `mean_core_share`. On `ps9` at
detection = 0.001, prevalence = 0.5: 5 of 29 genera, accounting for a mean
**77.7%** of each sample — *Streptococcus*, *Escherichia-Shigella*,
*Staphylococcus*, *Corynebacterium*, *Lawsonella*.

## Agglomeration loses reads, and the loss must be stated

```r
ag <- agglomerate_and_summarise(ps, rank = "Genus")
```

`tax_glom` drops taxa with NA at the target rank. "We analysed at genus level"
often means "we discarded the 30% of reads unclassified at genus". The function
reports how many taxa were unresolved and what share of reads went with them, and
warns above 5%.

Rank availability varies by source: amplicon taxonomies often stop at Genus,
shotgun profiles usually resolve to Species. `safe_rank()` returns NULL for an
absent rank rather than pushing an unresolvable column name into phyloseq.

**Agglomerating before differential abundance is a legitimate power move** — 143
genera is a smaller multiple-testing burden than 619 ASVs — but the choice is made
before testing, not after seeing which level gives more hits.

## Abundance versus prevalence

`taxa_summary()` reports both, because they answer different questions and are
routinely confused:

- high mean, low prevalence — abundant in a few samples, absent elsewhere
- low mean, high prevalence — present everywhere at low level
- `mean_when_present` separates "rare because absent" from "rare because scarce"

## Stacked bars

`collapse_to_top(ps, rank, top_n)` collapses everything beyond the top N into
`Other` and records what share that is. A stacked bar showing the top 10 of 140
taxa is not showing the community, and the caption must say so —
`plot_stacked_bar()` puts the collapsed count and mean `Other` share there
automatically.

`Other` belongs last in the stack, not sorted among the named taxa.

Always on relative abundance (`03`): a stacked bar of raw counts shows sequencing
depth.

## Dominance

```r
dominance_summary(ps, group_var = "region")
```

Reports the share held by the single most abundant taxon, by the top five, and
how many taxa it takes to reach 50% of a sample. This makes "dominated"
quantitative.

On `ps9`: mean top-1 share 56.3%, mean top-5 share 89.9%, median **1** taxon to
reach 50%. That is a strongly dominated community, and it explains why Simpson and
Shannon diverge from richness on this data — and why the `Other` category in a
top-8 bar plot is nearly empty.

## Functional profiles

PICRUSt2 and HUMAnN tables load through `build_phyloseq_functional()` with
pathways occupying the taxonomy slot, giving a `Pathways` rank (Demo 6).

Composition, differential abundance and ordination all work unchanged. Diversity
metrics need more care: "pathway richness" is bounded by the reference pathway
database, not by biology, and a plateau in a rarefaction curve of pathways means
the database is exhausted rather than the community sampled. Report functional
composition and differential pathway abundance; be cautious with functional alpha
diversity.
