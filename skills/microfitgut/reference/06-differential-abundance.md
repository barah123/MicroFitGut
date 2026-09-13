# Differential abundance

Which taxa differ between groups. The stage with the widest spread between
methods on identical data, and therefore the one where the method is chosen
before looking at results and stated plainly afterwards.

## Design → method

`choose_da_method(ps, group_var)` applies this and returns the reasoning.

| Design feature | Method | Why |
|---|---|---|
| **Repeated measures** | ANCOM-BC2 with `rand_formula` | The only one of the four that models within-subject correlation. The others assume independence and overstate significance. Alternative: per-taxon negative-binomial GLMM (`07`). |
| **Covariates to adjust for** | ANCOM-BC2 or DESeq2 | Both accept a multi-term design. Plain Wilcoxon cannot adjust for anything. |
| **Smallest group < 10** | ALDEx2 | Monte Carlo sampling propagates the uncertainty small counts carry rather than assuming it away, and it reports dispersion so a borderline call can be judged. DESeq2 has more power here but its shrinkage is doing heavy lifting. |
| **> 85% zeros** | ANCOM-BC2 with `struc_zero = TRUE` | It distinguishes structural zeros (absent) from sampling zeros (present, unobserved). The other three conflate them. |
| Balanced, adequate n, independent | ANCOM-BC2 | Bias correction addresses the compositional problem that makes raw counts misleading. |
| Quick sanity check alongside a proper method | Wilcoxon/Kruskal-Wallis on CLR | Assumption-light, no model of compositionality. Never on its own. |

**All four require raw counts.** See `03` — this is the commonest error in this
stage.

## The four methods

### ALDEx2 — `da_aldex2()`

Dirichlet-multinomial Monte Carlo sampling, then CLR per instance, then a rank or
t test on the transformed values. Because it models the counting uncertainty
explicitly, it behaves best at small n and it reports dispersion alongside effect
size.

**Which p-value column.** ALDEx2 returns several and they mean different things:

| Column | Test | Assumes normality |
|---|---|---|
| `we.ep` / `we.eBH` | Welch's t | **yes** |
| `wi.ep` / `wi.eBH` | Wilcoxon rank-sum | no |
| `kw.ep` / `kw.eBH` | Kruskal-Wallis (>2 groups) | no |
| `glm.ep` / `glm.eBH` | GLM (>2 groups) | — |

Use `wi.eBH` for two groups and `glm.eBH` for more than two. `da_aldex2()` picks
these automatically and records which it used.

**`mc.samples`.** The authors recommend at least 128, and **1,000 for a rigorous
effect size**. The default here is 128, matching the course material; raise it
before reporting effect sizes as findings. The function warns when it has not been
raised.

`aldex.plot(type = "MW")` shows difference against dispersion — the taxa where
between-group difference exceeds within-group variability are the trustworthy
calls.

### ANCOM-BC2 — `da_ancombc2()`

Log-linear model estimating each taxon's true log absolute abundance while
correcting the sampling-fraction bias that makes relative abundances misleading.
Wald test on the log fold change.

Parameters that must be reported:

| Parameter | Default here | Effect |
|---|---|---|
| `prv_cut` | 0.10 | Taxa below this prevalence are dropped before testing, changing the multiple-testing denominator. |
| `lib_cut` | 1000 | Samples below this depth are dropped. |
| `struc_zero` | TRUE | A taxon absent from an entire group is declared a structural zero and handled separately rather than tested as a fold change against zero. |
| `neg_lb` | FALSE | Whether to use the lower bound in declaring structural zeros — more aggressive. |
| `p_adj_method` | holm | Holm is much more conservative than BH. Demo 9 uses holm; BH is commoner in the literature. State which. |
| `rand_formula` | NULL | `"(1 | patient)"` for repeated measures. |

**Reading the output.** `res` carries one column set per model term. The intercept
columns are a baseline, not a comparison — a "significant" intercept only says the
taxon's abundance differs from zero in the reference group. `da_ancombc2()`
excludes intercept columns from the finding counts.

`diff_<term>` is the primary call. `diff_robust_<term>` is the subset that also
survives the pseudo-count sensitivity analysis, and it is the more defensible
number to headline. On the test run: 8 significant, 3 robust. Report both.

`out$res_global` is the whole-factor test for >2 groups; `res_pair` and `res_dunn`
are the follow-ups, and require `pairwise = TRUE` / `dunnet = TRUE`.

`zero_ind` lists the structural zeros. A structural zero is a strong claim — the
taxon is absent from that group, not merely unobserved — and at low sequencing
depth it is often undersampling instead. See `11`.

### DESeq2 — `da_deseq2()`

Negative binomial GLM with dispersion shrunk toward a fitted trend. The shrinkage
is where its power at small n comes from, and also where its assumptions bite.

```r
d <- da_deseq2(ps, ~ region, test = "Wald", fitType = "local", alpha = 0.1)
```

- **`fitType = "local"`** is the default because microbiome dispersion–mean
  relationships are often not well described by the parametric fit. Demo 10 uses
  local explicitly.
- **`alpha = 0.1` is DESeq2's own default and it is load-bearing.** `results()`
  optimises independent filtering for the alpha you pass, so passing 0.1 and then
  calling 0.05 significant is not the same as passing 0.05. Pass the threshold you
  will actually use.
- **Two-level factors** use the Wald test. **Multi-level factors** need
  `test = "LRT", reduced = ~ 1` for the overall effect, then Wald contrasts per
  pair (Demo 10).
- **Shrunken log fold changes** are the ones to report. Unshrunken LFCs for
  low-count taxa are wildly inflated and make volcano plots misleading.
  `shrink_lfc = TRUE` does this via apeglm, falling back to normal.

**Size factors, and when DESeq2's default is degenerate.** `da_deseq2()` always
estimates size factors from the positive-count geometric mean over *all* taxa —
the documented phyloseq workaround (Demo 10, ps10 Q1) — rather than letting
DESeq2 use its default.

That is not just defensive. DESeq2's default computes its reference geometric
mean only over rows with **no zeros at all**, so on a sparse table a handful of
taxa can end up carrying every sample's size factor. On the bundled `ps10.RDS`,
exactly **1 of 51 taxa** is present in every sample: the default therefore
derives all 209 size factors from that single taxon, giving a range of
0.002–18.75 with a median of 1.77 — where a median-of-ratios estimate should
centre near 1. The all-taxa estimate gives a median of 0.976, and the two
vectors correlate at **−0.10**. They are not small variations on each other; they
are different normalizations, and they produce different hit lists (30 vs 34
taxa at q < 0.1, with the direction split moving from 27/3 to 17/17).

`da_deseq2()` logs how many taxa were complete and says so in its message when
the default would have been degenerate, because that is the fact that explains a
divergence against any analysis which used the default.

**Independent filtering must be reported.** DESeq2 sets `padj` to NA for
low-mean-count taxa. That is not a missing value — it is a deliberate exclusion
that increases power on the rest. On the test run, 8 of 49 taxa were filtered, so
the FDR denominator is 41, not 49. "31 of 41 significant" and "31 of 49
significant" are different claims. `da_deseq2()` reports both and carries the
caveat.

**The diagnostics are not decoration.** `da_deseq2_diagnostics()` produces three
plots and quantifies what the eye is meant to check:

- `plotDispEsts` — black points should scatter around the red trend. If they do
  not, the negative binomial assumption is in trouble and the p-values are not
  trustworthy.
- `plotMA` — mean abundance against log fold change. Blue points are significant.
- **p-value histogram** — should be roughly flat with a spike near zero. The
  function computes the fraction of p-values above 0.5 (expected ~50% under the
  null) and reads it: too *few* large p-values means anti-conservative, which
  points at unmodelled structure such as batch or repeated measures inflating
  significance. On the skin data this fires, correctly, because `patient` is a
  repeated measure that `~ region` does not model.

### Wilcoxon / Kruskal-Wallis per taxon — `da_kruskal()`

The loop from Demo 7 and Demo 9, made explicit about its limits. No model of
compositionality, no covariate adjustment, no dispersion estimate, assumes
independent samples.

Run it on CLR values, apply BH, and report it **alongside** a compositional
method, never instead of one.

## Multiple testing

Every method corrects, but they differ in default and in denominator:

| Method | Default correction | Denominator |
|---|---|---|
| ALDEx2 | BH (built into the `eBH` columns) | taxa passing the zero-row filter |
| ANCOM-BC2 | holm (configurable) | taxa passing `prv_cut` and `lib_cut` |
| DESeq2 | BH | taxa surviving independent filtering |
| `da_kruskal` | BH | taxa passing `min_prevalence` |

**The denominator is as important as the method.** Filtering 600 taxa down to 143
before testing is not cheating — it increases power legitimately — but the report
must state how many taxa were tested, not how many were in the raw table.

## Comparing methods — `da_compare()`

The honest presentation. On the test run:

```
kruskal      29 significant, 1 unique
aldex2       19 significant, 0 unique
deseq2       31 significant, 6 unique
consensus (all three): 16 taxa
union: 35 taxa

kruskal vs aldex2   Jaccard 0.655  Spearman(effect) 0.983
kruskal vs deseq2   Jaccard 0.714  Spearman(effect) 0.815
aldex2  vs deseq2   Jaccard 0.471  Spearman(effect) 0.850
```

**Report the consensus set as the primary finding and the method-specific sets as
method-dependent.** A taxon found by one method only is a hypothesis; a taxon
found by all is a result.

The pattern above — high effect-size correlation with moderate hit-list overlap —
is the signature of methods agreeing about direction and magnitude while
disagreeing about where to draw the significance line. That is a threshold and
power difference, not a different biological answer, and saying so is more useful
than picking the method with the most hits.

## Effect sizes

Never report a hit list without them. `da_significant()` normalises all four
methods to `taxon / effect / q_value`:

| Method | Effect | Scale |
|---|---|---|
| ALDEx2 | median CLR difference, or `effect` | CLR units |
| ANCOM-BC2 | `lfc_<term>` | log fold change, bias-corrected |
| DESeq2 | `log2FoldChange` (shrunken) | log2 fold change |
| `da_kruskal` | rank-biserial r or eta-squared | standardised |

These are **not on the same scale**, so do not average them across methods.
Compare ranks, not values.

## Direction

A log fold change is relative to the factor's reference level, which is
alphabetical unless set (`01`). Before interpreting any sign, confirm which level
is the reference — `DESeq2::resultsNames(dds)` and the ANCOM-BC2 column names
both show it. A flipped reference level flips every reported direction and does
not error.
