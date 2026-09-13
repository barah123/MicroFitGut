# Alpha diversity

Within-sample diversity. The stage where an analysis most often looks rigorous
and is not.

## Metrics and what each one answers

| Metric | Question | Needs counts | Depth-sensitive |
|---|---|---|---|
| Observed | How many taxa were seen? | no | **very** |
| Chao1 | How many taxa exist, including unseen ones? | **yes** | **very** |
| ACE | Same, weighting abundant and rare differently | **yes** | **very** |
| Fisher's alpha | Parametric richness, assuming a log-series distribution | **yes** | high |
| Shannon | Richness and evenness together | no | low |
| Simpson | Probability two random draws are the same taxon (lower = more diverse) | no | low |
| InvSimpson | Inverse of Simpson (higher = more diverse) | no | low |
| Faith's PD | Total branch length spanning the sample's taxa | no | high |
| Pielou, Camargo, Evar, Bulla | Evenness alone | no | low |

Chao1, ACE and Fisher are estimated from **singletons and doubletons** — taxa
seen exactly once or twice. That is why they need integer counts and why they
require rarefaction: on unequal depth they measure sequencing effort.
`compute_alpha()` drops them automatically on relative-abundance input and says
so.

**Report richness and evenness, not one of them.** A community can be rich and
dominated, or poor and even. Shannon confounds the two; reporting Observed
alongside Pielou separates them.

### Faith's PD

```r
pd <- faith_pd(ps, include_root = TRUE)
```

Total phylogenetic branch length. No non-phylogenetic index substitutes for it,
which is why `picante` is a hard requirement rather than optional.

Two gates, both refusals rather than warnings:

- **The tree must be marked real.** PD on a random topology is a number with no
  meaning (`01`). Confirm with `mfg_mark_tree_real(ps, TRUE)`.
- **`include_root = TRUE` needs a rooted tree.** Root with `phangorn::midpoint()`.

## The workflow

```r
ps_r <- rarefy_to_depth(ps, depth = 5000, rngseed = 42, justification = "...")
ad   <- alpha_table(ps_r, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
```

`alpha_table()` returns one data frame holding every index **and** every metadata
column. Tests and plots both take this object, so a figure and a test can never
disagree about the numbers.

## The assumption sequence

```r
a <- check_assumptions(ad, "Shannon", "region")
```

1. **Shapiro-Wilk** for normality — pooled *and* per group. Per-group matters
   more for a group comparison; a pooled test can reject because the groups
   differ in mean, not because either is non-normal.
2. **Bartlett** for equal variances — more powerful, but assumes normality itself.
3. **Fligner-Killeen** for equal variances — robust to non-normality.

When Bartlett and Fligner disagree, **Fligner decides**, because it does not
assume the thing that is in question.

The output names the branch: `parametric`, `parametric_welch` (normal, unequal
variance), or `nonparametric`.

### What passing does not license

Failing to reject normality is not proving it. With n < 5 per group, Shapiro-Wilk
has almost no power to reject anything, so "normal" means "we could not tell".
`check_assumptions()` sets a `caveat` in that case and recommends the
non-parametric branch. This is a defensible default, not a proof.

## Test selection

`choose_alpha_test()` decides and returns the reasoning, which goes in the report.

| Situation | Test |
|---|---|
| 2 groups, normal, equal variance | t-test |
| 2 groups, normal, unequal variance | Welch's t-test |
| 2 groups, otherwise | Wilcoxon rank-sum |
| >2 groups, normal, equal variance | ANOVA |
| >2 groups, otherwise | Kruskal-Wallis |
| Complete unreplicated block design | Friedman |
| **Repeated measures, not a complete block** | **no simple test is valid → mixed model (`07`)** |

That last row is the one that matters, and it is a refusal rather than a
fallback. `run_alpha_test()` returns `test = "REFUSED"` with an explanation.

Friedman requires exactly one observation per group × block cell. Kruskal-Wallis
on repeated measures treats correlated samples as independent draws and
overstates significance. When the design is repeated-measures but unbalanced —
which is the common case — neither applies, and the answer is
`fit_lmm(Shannon ~ region + (1 | patient))`.

On the bundled skin data, `patient` has up to 3 samples per level across the
`ps9` subset with cell counts of 0–2, so the test is correctly refused.

## Effect sizes

`10-reporting-standards.md` requires an effect size, not just a p-value, so one
is always computed:

| Test | Effect size |
|---|---|
| t-test, Welch | Cohen's d, plus the CI for the difference |
| Wilcoxon | rank-biserial r |
| ANOVA | eta-squared |
| Kruskal-Wallis | eta-squared from H: `(H - k + 1) / (N - k)` |
| Friedman | Kendall's W |

## Post-hoc

**Only after a significant omnibus test.** Running pairwise comparisons after a
non-significant omnibus inflates the false-positive rate; `posthoc_dunn()`
refuses unless you pass `force = TRUE`, and then the report must say so.

| Function | Method | Correction |
|---|---|---|
| `posthoc_dunn()` | Dunn's test (FSA) | BH by default |
| `posthoc_nemenyi()` | Nemenyi (DescTools) | Tukey distribution |
| `posthoc_tukey()` | Tukey HSD, after ANOVA | built in |

Dunn's is the conventional follow-up to Kruskal-Wallis and is the default here.

## Testing several indices at once

```r
r <- alpha_test_all(ad, "region")
```

Returns one row per index with the test chosen, statistic, p, effect size, the
assumption p-values, and the reasoning.

**It applies BH across the indices.** Six indices against one variable is six
tests. They are strongly correlated, so BH across them is conservative in the
right direction rather than exact — but reporting six raw p-values as if each
were a separate discovery is not defensible. The correction is logged with that
caveat attached.

## What a null result means

"No significant difference in alpha diversity" is weak evidence of no difference
unless the study had power to detect one. Before writing it, check
`12-power-and-sample-size.md`. With n = 10 per group, a Wilcoxon test detects
only large effects, and a null result is close to uninformative.
