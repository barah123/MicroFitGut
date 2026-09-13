# Beta diversity

Between-sample community structure. Distances, ordination, and the three
permutation tests.

## The central trap

**A significant PERMANOVA does not mean the groups differ in composition.**

`adonis2` tests whether group membership explains variation in the distance
matrix. Heterogeneous within-group dispersion produces a significant result with
identical group centroids — the groups differ in how *variable* they are, not in
where they sit.

So every PERMANOVA gets a `betadisper` check, and `run_permanova()` computes it
automatically and qualifies the conclusion. The demos do not emphasise this.
Reviewers ask for it.

On the bundled skin data this fires: PERMANOVA on `region` gives R² = 0.30,
p = 0.005, and betadisper gives p = 0.005 with mean distance-to-centroid of 0.272
(BE), 0.325 (BT), 0.468 (FA). The significant PERMANOVA **cannot** be reported as
a composition shift. `permanova_sentence()` writes the qualified version.

## Distance choice

| Metric | Uses | Sensitive to | When |
|---|---|---|---|
| Bray-Curtis | abundance | abundant taxa | Default. Interpretable, well understood. |
| Jaccard | presence/absence | rare taxa | When occurrence matters more than abundance. |
| JSD | abundance (information-theoretic) | whole distribution | Alternative to Bray-Curtis; bounded, symmetric. |
| Unweighted UniFrac | presence + phylogeny | rare lineages | When lineage identity matters and rare taxa are trustworthy. |
| Weighted UniFrac | abundance + phylogeny | abundant lineages | The phylogenetic counterpart of Bray-Curtis. |
| Aitchison | CLR values | compositional structure | The compositionally coherent choice. Requires CLR. |
| Euclidean | CLR or VST values | variance structure | Only on variance-stabilised data. |

**Unweighted versus weighted UniFrac is a real analytical choice, not a detail.**
Unweighted is driven by rare lineages, so it is the more sensitive test and the
one most affected by depth and filtering. Weighted tracks the abundant community.
They routinely disagree, and that disagreement is informative: a signal in
unweighted but not weighted means the difference lies in which rare lineages are
present, not in community structure.

Both require a **real tree** — `compute_distance()` refuses otherwise (`01`).

## Ordination

```r
d   <- compute_distance(ps_tss, "bray")
ord <- run_ordination(ps_tss, "PCoA", dist_obj = d)
```

| Method | Constrained | Reports | Note |
|---|---|---|---|
| PCoA | no | variance explained per axis | Metric MDS. The default. |
| NMDS | no | **stress** | Rank-based; no variance explained. |
| PCA | no | variance explained | Needs CLR or VST. |
| CCA | **yes** | constrained fraction of inertia | Unimodal species response. |
| RDA | **yes** | constrained fraction | Linear species response. |

### The constrained-versus-unconstrained decision

Unconstrained (PCoA, NMDS, PCA) asks: what is the dominant structure in this
data? It is blind to the design, so separation by group in an unconstrained plot
is real evidence.

Constrained (CCA, RDA) asks: how much variation do *these specific variables*
explain? It is told the answer, so separation in a constrained plot is not
evidence of anything by itself — the constrained fraction and its permutation
test are the evidence.

Use unconstrained by default. Reach for constrained when you have several
candidate explanatory variables and want to partition variance between them.
`run_ordination(method = "CCA")` requires a formula and errors without one.

### Stress — the only thing that says whether an NMDS can be trusted

| Stress | Reading |
|---|---|
| < 0.05 | excellent |
| < 0.10 | good |
| < 0.20 | usable, interpret with care |
| ≥ 0.20 | **the plot is misleading — do not read it as a map of distance** |

At stress ≥ 0.20 the two-dimensional solution cannot represent the distances, and
the picture will suggest structure that is an artefact of the projection.
`run_ordination()` warns. Raise `k` or use PCoA.

An NMDS plot without its stress value cannot be assessed, which is why
`plot_ordination_std()` puts it in the subtitle.

### Variance explained

A PCoA plot of axes explaining 17% looks identical to one explaining 80%. The
percentages belong on the axis labels — `plot_ordination_std()` does this — and
the reader needs them to judge whether the visible separation is the dominant
structure or a slice of it.

## The three tests

### PERMANOVA — `run_permanova()`

The primary test. Partitions variation in the distance matrix.

```r
pm <- run_permanova(ps_tss, ~ region + wash, dist_obj = d,
                    permutations = 999, by = "margin", strata = "patient")
```

- **`by = "margin"` is the default** because it gives each predictor its own
  p-value independent of the order terms were entered. Sequential (`by = "terms"`)
  attributes shared variance to whichever term came first, which makes the result
  depend on how you typed the formula. Demo 7 notes marginal as preferable.
- **R² is the effect size** and must be reported. A significant term explaining
  4% of variation is a different finding from one explaining 30%.
- **`strata` blocks permutations** within a level. This is how repeated measures
  are handled: permuting freely across subjects treats within-subject samples as
  exchangeable when they are not. `run_permanova()` warns when repeated measures
  are detected and no strata is given.

### betadisper — `check_dispersion()`

Not optional. Reduces each sample to its distance from its group centroid and
tests whether those distances differ by group.

A significant betadisper does not invalidate the PERMANOVA — it changes what the
PERMANOVA licenses you to say:

- **Homogeneous dispersion:** "community composition differed by group" is
  defensible; the effect is a difference in centroids.
- **Heterogeneous dispersion:** only "the groups differ in their distance-matrix
  structure" is defensible. Report which group is more variable and give the mean
  distances to centroid.

`plot_dispersion()` is the visual form and belongs alongside the ordination.

### ANOSIM — `run_anosim()`

R is bounded −1 to 1. Positive means between-group distances exceed within-group
distances, which is the expected direction for real structure.

**A negative R means samples within a group are more dissimilar than samples from
different groups** — the opposite of grouping having biological meaning (Demo 7).

ANOSIM is *more* sensitive to dispersion differences than PERMANOVA, so it is
reported as a supporting test, never as the primary one. Rough reading of R:
< 0.25 very weak, < 0.5 small to moderate, < 0.75 moderate to strong.

### MRPP — `run_mrpp()`

Compares observed within-group dissimilarity (delta) against its permutation
expectation. Lower observed delta than expected means within-group samples are
more similar than chance — consistent with real grouping. A within the range of
0.1–0.3 is a moderate effect.

## Clustering

```r
hc <- run_hclust(d, method = "ward.D2", k = 3)
```

`ward.D2` is the default: it is the correct Ward implementation for squared
distances and gives the compact clusters the literature conventionally shows.
Linkage changes the tree materially — single linkage chains, complete linkage
forces compactness, average sits between — so the choice is recorded.

**Cophenetic correlation** says how faithfully the tree represents the original
distances. Below about 0.7 the dendrogram is a poor summary and should not be
presented as one. On the skin data it is 0.698, which is borderline and worth
stating.

```r
cluster_vs_group(hc$clusters, meta$region)
```

**Adjusted Rand index** turns "the groups separate" into a number: 0 is chance
agreement, 1 is identical. On the skin data, ARI = 0.50 — real structure, not a
clean partition.

k-means (`run_kmeans()`) is the non-hierarchical alternative and needs `k` chosen
in advance; it reports the between-cluster share of variance.

## Reporting

`permanova_sentence(pm)` writes one defensible sentence per term and
**deliberately refuses to write "community composition differed" when the
dispersion check does not license it.** Use it rather than composing the sentence
by hand — that is where the qualification gets dropped.

Every beta-diversity result must state: the distance metric, the normalization it
was computed on, the number of permutations, whether permutations were stratified,
R², and the betadisper result.
