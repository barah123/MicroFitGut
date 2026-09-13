# Normalization

**Read this before any comparative statistic.** This is the choice that most
often silently changes conclusions, because a wrong pairing produces a clean
result that answers a different question.

## The problem

Sequencing depth is an artefact of the run, not biology. In the bundled skin
dataset it varies **125-fold** across samples (1,107 to 137,913 reads). Any
statistic that responds to depth — richness, Bray-Curtis distance, raw count
comparisons — is partly measuring library size.

Worse, microbiome data is **compositional**: a table of proportions is
constrained to sum to one, so an increase in one taxon forces the others down
arithmetically. Correlations between proportions are spurious by construction,
and a "decrease" in one taxon can be nothing but another taxon's increase.

Different methods address different halves of this, and a method that is right
for one downstream analysis is wrong for another.

## The decision table

`mfg_check_normalization(ps, analysis)` enforces this. It **refuses** invalid
pairings rather than letting them run. The machine-readable form is
`MFG_VALID_NORMALIZATION` in `scripts/03-normalize.R`; if the two ever disagree,
this file is the authority.

| Downstream analysis | Valid normalization | Why nothing else |
|---|---|---|
| Richness: Observed, Chao1, ACE, Fisher | `rarefied` | These are estimated from taxa seen exactly once or twice. On unequal depth they measure sequencing effort, not richness. |
| Faith's PD | `rarefied` | Same — PD is a richness-family metric on the tree. |
| Shannon, Simpson, InvSimpson | `rarefied`, `tss` | Defined on proportions and far less depth-sensitive than richness. Rarefying is still the safer default. |
| Bray-Curtis, JSD | `tss`, `rarefied`, `vst`, `clr` | Raw counts encode library size as community difference. |
| Jaccard | `tss`, `rarefied`, `presence` | Presence-based; depth changes what counts as present. |
| UniFrac (weighted/unweighted) | `rarefied`, `tss` | Plus a **real tree** — see `01`. |
| Aitchison | `clr` only | Aitchison distance *is* Euclidean distance on CLR values, by definition. |
| PCA, Euclidean distance | `clr`, `vst` | Euclidean geometry needs variance independent of the mean. |
| PCoA | `tss`, `rarefied`, `clr`, `vst` | Whatever the distance accepts. |
| NMDS | `tss`, `rarefied` | Rank-based, so monotone transforms are wasted; keep it interpretable. |
| CCA/RDA | `tss`, `raw` | vegan's `cca` applies its own chi-square standardisation. |
| **DESeq2** | `raw` only | It estimates its own size factors and models the mean–variance relationship. Pre-normalized input breaks the dispersion shrinkage its power depends on. |
| **ANCOM-BC2** | `raw` only | Its whole purpose is estimating and correcting the sampling-fraction bias from counts. |
| **ALDEx2** | `raw` only | Monte Carlo Dirichlet-multinomial sampling needs integer counts to sample from. |
| Wilcoxon / Kruskal-Wallis per taxon | `clr`, `tss`, `log10`, `rarefied` | Rank tests tolerate any monotone transform; CLR is preferred because it removes the compositional constraint. |
| ZIP/ZINB/ZHP/ZHNB, count GLMM | `raw` + log-depth offset | Count models. The offset handles depth; pre-normalizing double-corrects. |
| Gaussian LM/LMM | `clr`, `log10`, `tss` | A Gaussian model on raw counts assumes constant variance, which counts violate. |
| Stacked bar plots, core microbiome | `tss` | A stacked bar of raw counts shows depth, not composition. |
| Heatmaps | `clr`, `vst`, `log10` | A heatmap of proportions is dominated by the few abundant taxa and shows nothing else. |

**The most important row is DESeq2/ANCOM-BC2/ALDEx2 taking raw counts only.** The
instinct to "normalize first, then test" is exactly wrong for these three, and it
is the single commonest error in this stage.

## The methods

### TSS — total sum scaling (relative abundance)

```r
ps_rel <- tss_transform(ps)
```

Divide each sample by its total. Simplest depth correction, right for composition
plots and most distances. Does **not** address compositionality: proportions still
sum to one.

### Rarefaction

```r
choose_rarefaction_depth(ps)                  # see the trade-off first
ps_r <- rarefy_to_depth(ps, depth = 5000, rngseed = 42,
                        justification = "curves plateau above 5000")
```

Subsample every sample to the same depth without replacement. The only defensible
input for richness estimators.

**Choosing the depth is a real trade-off and there is no correct answer.** On
`ps9.RDS`:

| depth | samples kept | reads retained |
|---|---|---|
| 518 (minimum) | 209 / 209 (100%) | **1.8%** |
| 1,429 | 198 (94.7%) | 4.8% |
| 11,269 | 157 (75.1%) | 30.0% |
| 24,707 | 105 (50.2%) | 44.0% |

Rarefying to the minimum keeps every sample and throws away 98% of the data.
Rarefying high keeps the data and drops a quarter of the samples. Pick one,
record it, and show the curves.

**Two things are mandatory.** The seed — rarefaction is random, so an unseeded run
is not reproducible and two runs will disagree on richness. And `replace = FALSE`,
which is what rarefaction means; phyloseq's default `TRUE` is faster but is not
the classical procedure.

### Checking the depth plateaus

```r
rc <- mfg_rarefaction_curve(ps_rarefied, step = 100)
plot_rarefaction(rc, depth = 5000)
```

**Call this on the rarefied object, not the raw one.** On an unrarefied object
every curve runs out to that sample's own total depth, so each one necessarily
ends flat and "has it plateaued" is trivially true. The function warns when the
curve end-depths vary, which is the signature of that mistake.

Plateau is judged by the **terminal slope in taxa per 1000 reads** — how many new
taxa another 1000 reads would reveal — with a default threshold of 1.

It is deliberately *not* judged by the terminal slope relative to the initial
slope. The start of any rarefaction curve is steep, because the first reads are
nearly all new taxa, so a ratio against it declares a plateau almost regardless
of saturation. On the bundled `ps9.RDS` rarefied to 518 reads, the ratio
criterion called 37 of 40 curves plateaued while the terminal slope is a median
**3.0 taxa per 1000 reads** — still clearly climbing, which is what the course
material also concludes. The absolute criterion separates that case from
`Demo9.RDS` at 5000 reads (median 0.57, genuinely flat) and at 50 reads
(median 20.0, steeply climbing).

Observed/Chao1 coverage does *not* work as a substitute here: at very shallow
depth there are no doubletons, Chao1 collapses onto Observed, and coverage reads
1.000 on the most undersampled data.

**The rarefaction dispute.** McMurdie & Holmes (2014) argued rarefying is
statistically inadmissible for differential abundance because it discards data
and adds variance. That critique is correct *for differential abundance*, which
is why every DA method in this library takes raw counts. It does not extend to
richness estimation, where equal effort is the only way the estimator means
anything. Both positions are defensible in their own domain; state which you used
and for what. See `11`.

### CLR — centred log-ratio

```r
ps_clr <- clr_transform(ps)
```

Each value divided by the sample's geometric mean, then logged. This moves
compositional data into Euclidean space, where ordinary statistics apply.
`microbiome::transform("clr")` adds a pseudocount internally for zeros.

**Interpretation caveat that must reach the report:** a CLR value is a statement
about a taxon *relative to the rest of that sample*, not an absolute abundance. A
positive CLR coefficient means the taxon rose relative to the sample's geometric
mean — which can happen because it increased or because everything else fell.

The pseudocount is not innocent at high sparsity. With 90%+ zeros, most of the
table is pseudocount, and CLR values for rare taxa are dominated by that choice.
ANCOM-BC2's sensitivity analysis exists for exactly this reason and reports which
findings survive it (`diff_robust_*`).

### log10

```r
ps_log <- log10_transform(ps)
```

`log10(1 + x)` of relative abundance. Cheaper than CLR and easier to explain, but
does not remove the compositional constraint. Fine for heatmaps and visual
scaling; weaker than CLR for inference.

### DESeq2 VST

```r
ps_vst <- vst_transform(ps, blind = TRUE, fitType = "local")
```

Removes the dependence of variance on the mean, which is what makes PCA and
Euclidean distance behave on count data.

`blind = TRUE` keeps the transformation independent of the design — which is what
you want when the values feed an unsupervised method, otherwise the ordination
has partly been told the answer.

Negative values in the output are expected and valid.

### DESeq2 median-of-ratios counts

```r
ps_n <- deseq2_normalize(ps)
```

Useful for plots and for handing normalized counts to something that expects
them. **Do not then run DESeq2 on this object** — it needs raw counts.

The geometric mean is computed over positive counts only. DESeq2's default needs
at least one taxon present in every sample, which microbiome data rarely has, so
the default gives zero everywhere and estimation fails. This is the standard
phyloseq workaround (Demo 10, ps10):

```r
gm_mean <- function(x, na.rm = TRUE) exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
```

### Presence/absence

```r
ps_pa <- presence_transform(ps, detection = 0)
```

For Jaccard and for questions about occurrence rather than abundance.

## Not implemented, and why

**CSS (cumulative sum scaling, metagenomeSeq).** Reasonable and widely used. Not
included because it is not in the source material and adding a normalization the
reference files cannot defend would widen the library without deepening it. If a
paper used CSS, record that in the study card and note it as an untested axis.

**GMPR, TMM, upper-quartile.** Same reasoning.

## Tracking

Every function here stamps the object, and `mfg_normalization(ps)` reads it back.
An object arriving from outside the pipeline reports `*_untracked`, and the
pairing guard **refuses to proceed** on an untracked object — because a report
cannot state what was done if nothing recorded it.

To assert the state of an object you loaded yourself:

```r
ps <- mfg_set_normalization(ps, "raw")
```

`mfg_normalization_sentence(ps)` writes the methods-section sentence, including
the rarefaction depth, seed and sample loss.

**Normalization is chosen per analysis, not once per run.** A typical run holds
three states of the same data simultaneously: rarefied for alpha, TSS for
ordination, raw for differential abundance. That is correct, not wasteful.
